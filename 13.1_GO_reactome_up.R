###################################################
## Libraries
###################################################
library(DESeq2)
library(qvalue)
library(annotables)
library(dplyr)
library(tidyr)
library(tidyverse)
library(clusterProfiler)
library(ReactomePA)
organism <- "org.Hs.eg.db"
library(organism, character.only = TRUE)
library(ggplot2)
library(data.table)
library(reshape)
library(stringr)
library(ggrepel)
library(scales)
library(readr)   # needed for write_rds
library(grid)    # for unit() in ggplot

####################################################
rm(list = ls())

setwd("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/GO/updown/up/")
outFolder <- getwd()

dataDir <- "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/deseqres/"

platePrefix <- "ALL.0.1.13."
plateSufix  <- "-RNA-CTRL.noCombat_DESeq"

timestamp()

####################################################
## Variables to loop over
variables      <- c("PSS_all_mean", "ISEL_Mean", "cytocomp")
variable_names <- c("Psychological Stress", "Social Support", "Cytokines")

####################################################
##################### Outer loop #####################
for (i in seq_along(variables)) {

  tryCatch({

    myvar    <- variables[i]
    var_name <- variable_names[i]

    cat("\n######################################################################\n")
    cat("## Processing:", var_name, "\n")

    ##### Directories #####
    enrichGODir   <- file.path(outFolder, "enrichGO")
    enrichGOdata  <- file.path(enrichGODir, "stats", "celltype")
    enrichGOrds   <- file.path(enrichGODir, "data",  "celltype")
    enrichGOplot  <- file.path(enrichGODir, "plots", "celltype")
    DEGdirGO      <- file.path(enrichGODir, "DEGtable", "celltype")

    enrichPathDir   <- file.path(outFolder, "enrichPath")
    enrichPathdata  <- file.path(enrichPathDir, "stats", "celltype")
    enrichPathrds   <- file.path(enrichPathDir, "data",  "celltype")
    enrichPathplot  <- file.path(enrichPathDir, "plots", "celltype")
    DEGdirPath      <- file.path(enrichPathDir, "DEGtable", "celltype")

    dir.create(enrichGODir,     recursive = TRUE, showWarnings = FALSE)
    dir.create(enrichGOdata,    recursive = TRUE, showWarnings = FALSE)
    dir.create(enrichGOrds,     recursive = TRUE, showWarnings = FALSE)
    dir.create(enrichGOplot,    recursive = TRUE, showWarnings = FALSE)
    dir.create(DEGdirGO,        recursive = TRUE, showWarnings = FALSE)
    dir.create(enrichPathDir,   recursive = TRUE, showWarnings = FALSE)
    dir.create(enrichPathdata,  recursive = TRUE, showWarnings = FALSE)
    dir.create(enrichPathrds,   recursive = TRUE, showWarnings = FALSE)
    dir.create(enrichPathplot,  recursive = TRUE, showWarnings = FALSE)
    dir.create(DEGdirPath,      recursive = TRUE, showWarnings = FALSE)

    ##### Initialize #####
    dfvar <- list()
    dfsig <- list()

    clusters  <- 0:4
    celltypes <- c("R0 T CD4+", "R1 T CD8+", "R2 NK", "R3 Monocyte", "R4 B")

    ####################################################
    ## Load DESeq results
    ####################################################
    for (cl in clusters) {

      fname <- paste0(
        dataDir,
        platePrefix, "C", cl,
        ".deseqres_", myvar, plateSufix, ".txt"
      )

      res_data <- read.table(fname, sep = "\t", header = TRUE,
                             quote = '"', comment.char = "")

      dfvar[[as.character(cl)]] <- res_data
      dfsig[[as.character(cl)]] <- res_data %>% filter(padj < 0.1 & logFC > 0)

      cat("Cluster", cl, ":", nrow(dfvar[[as.character(cl)]]),
          "genes,", nrow(dfsig[[as.character(cl)]]), "UP DEGs\n")
    }

    ####################################################
    ## Filter out clusters with <50 DEGs
    ####################################################
    deg_counts <- data.frame()
    for (cl_name in names(dfvar)) {

      n_genes   <- nrow(dfvar[[cl_name]])
      deg_count <- sum(dfvar[[cl_name]]$padj < 0.1 &
                       dfvar[[cl_name]]$logFC > 0, na.rm = TRUE)

      deg_counts <- rbind(deg_counts, data.frame(
        Cluster = cl_name,
        nGenes  = n_genes,
        DEGs    = deg_count,
        Var     = myvar
      ))

      if (deg_count < 50) {
        dfvar[[cl_name]] <- NULL
        dfsig[[cl_name]] <- NULL
      }
    }

    write.table(deg_counts,
      file = file.path(DEGdirGO, paste0(myvar, "_number_of_DEGs_per_clusterv2.txt")),
      sep = "\t", quote = FALSE, row.names = FALSE
    )

    print(deg_counts)
    print("Remaining clusters:")
    print(names(dfvar))

    ####################################################
    ## Enrichment initialization
    ####################################################
    enrichlist    <- list()
    enrich_filt   <- list()
    pathwaylist   <- list()
    pathway_filt  <- list()

    ####################################################
    ## Loop over valid clusters
    ####################################################
    for (cl_name in names(dfvar)) {

      cluster  <- cl_name
      c_index  <- as.numeric(cluster)
      clustname <- celltypes[c_index + 1]

      tryCatch({

        # add symbol column
        dfvar[[cluster]]$symbol <- dfvar[[cluster]]$identifier
        dfsig[[cluster]]$symbol <- dfsig[[cluster]]$identifier

        # annotation
        anno <- grch38 %>%
          dplyr::select(symbol, entrez) %>%
          distinct()

        # corrected filter syntax
        annoVar <- anno %>% filter(symbol %in% dfvar[[cluster]]$symbol)

        dfvar[[cluster]] <- inner_join(dfvar[[cluster]], annoVar, by = "symbol") %>%
                            filter(!is.na(entrez))

        dfsig[[cluster]] <- inner_join(dfsig[[cluster]], annoVar, by = "symbol") %>%
                            filter(!is.na(entrez))

        if (nrow(dfsig[[cluster]]) == 0) {
          cat("No UP DEGs with Entrez IDs for cluster", cluster, "\n")
          next
        }

        ### Prepare gene lists
        genes    <- as.character(dfsig[[cluster]]$entrez)
        geneUniv <- as.character(dfvar[[cluster]]$entrez)

        ### enrichGO
        ego <- enrichGO(
          gene     = genes,
          universe = geneUniv,
          OrgDb    = org.Hs.eg.db,
          ont      = "BP"
        )

        enrichlist[[paste0("enrichGO_", cluster)]] <- ego

        ### enrichPathway
        epath <- enrichPathway(
          gene     = genes,
          universe = geneUniv,
          organism = "human"
        )

        pathwaylist[[paste0("pathway_", cluster)]] <- epath

        ### Process enrichGO
        if (!is.null(ego) && nrow(ego@result) > 0) {

          df_GO <- ego@result %>% filter(qvalue <= 0.1)

          if (nrow(df_GO) > 0) {

            df_GO$GeneRatio <- sapply(df_GO$GeneRatio, function(x){
              n <- unlist(strsplit(x,"/"))
              as.numeric(n[1])/as.numeric(n[2])
            })

            df_GO$var_cell <- paste0(myvar, "_", cluster)
            df_GO$cell     <- paste0("cluster_", cluster)
            df_GO$celltype <- clustname

            enrich_filt[[paste0("df_enrichGO_", cluster)]] <- df_GO
          }

        } else {
          cat("No enrichGO terms for cluster", cluster, "\n")
        }

        ### Process pathway
        if (!is.null(epath) && nrow(epath@result) > 0) {

          df_path <- epath@result %>% filter(qvalue <= 0.1)

          if (nrow(df_path) > 0) {

            df_path$GeneRatio <- sapply(df_path$GeneRatio, function(x){
              n <- unlist(strsplit(x,"/"))
              as.numeric(n[1])/as.numeric(n[2])
            })

            df_path$var_cell <- paste0(myvar, "_", cluster)
            df_path$cell     <- paste0("cluster_", cluster)
            df_path$celltype <- clustname

            pathway_filt[[paste0("df_pathway_", cluster)]] <- df_path
          }

        } else {
          cat("No enriched Reactome terms for cluster", cluster, "\n")
        }

        # save per-cluster RDS
        write_rds(
          ego,
          file.path(enrichGOrds, paste0("C_", cluster, "_", myvar, "_enrichGO_results.rds"))
        )
        write_rds(
          epath,
          file.path(enrichPathrds, paste0("C_", cluster, "_", myvar, "_enrichPathway_results.rds"))
        )

      },
      error = function(e) {
        print(paste("Error inside cluster", cluster, ":", e$message))
      },
      finally = {
        print(paste("Completed cluster", cluster))
      })

    } ## END cluster loop

    ####################################################
    ## Save all-cluster lists (correct)
    ####################################################
    write_rds(enrichlist,
      file.path(enrichGOrds, paste0(myvar, "_enrichGO_results_all_clust.rds")))

    write_rds(pathwaylist,
      file.path(enrichPathrds, paste0(myvar, "_enrichPathway_results_all_clust.rds")))

    ####################################################
    ## Combine GO results
    ####################################################
    if (length(enrich_filt) > 0) {

      df_GO_all <- bind_rows(enrich_filt)

      write.table(df_GO_all,
        file = file.path(enrichGOdata, paste0(myvar, "_CTRL_allclust_enrichGO_results_qvalue_0.1_DEG50.txt")),
        sep = "\t", quote = FALSE, row.names = FALSE
      )

      df_GO_all <- df_GO_all[order(df_GO_all$p.adjust), ]

      ## Top GO terms per cluster
      top_GO_list <- list()
      for (cl_id in unique(df_GO_all$cell)) {
        tmp <- df_GO_all %>%
          filter(cell == cl_id, p.adjust < 0.1) %>%
          head(10)
        if (nrow(tmp) > 0) top_GO_list[[cl_id]] <- tmp
      }

      if (length(top_GO_list) > 0) {

        top_GO <- bind_rows(top_GO_list) %>%
          drop_na(p.adjust)

        top_GO$celltype <- factor(top_GO$celltype, levels = celltypes)

        ## Dotplot
        p1 <- ggplot(top_GO, aes(x = celltype, y = str_wrap(Description, 60))) +
          geom_point(aes(size = GeneRatio, color = p.adjust)) +
          scale_color_gradient(low = "red", high = "blue") +
          ggtitle(var_name) +
          ylab(NULL) + xlab(NULL) +
          theme_bw() +
          theme(
            axis.text.x = element_text(angle = 45, hjust = 1, size = 50),
            axis.text.y = element_text(size = 50),
            text        = element_text(size = 50),
            legend.text = element_text(size = 40),
            legend.key.height = unit(2, "lines"),
            legend.title = element_text(size = 48)
          )

        png(file.path(enrichGOplot,
              paste0(myvar, "_CTRL_top10_enrichGO_DotPlot_loop_DEG50.png")),
            width = 3600, height = 3500, res = 120)
        print(p1)
        dev.off()
      }

    } else {
      cat("No GO terms to combine for", myvar, "\n")
    }

    ####################################################
    ## Combine Reactome
    ####################################################
    if (length(pathway_filt) > 0) {

      df_path_all <- bind_rows(pathway_filt)

      write.table(df_path_all,
        file = file.path(enrichPathdata, paste0(myvar, "_CTRL_allclust_enrichPath_results_qvalue_0.1_DEG50.txt")),
        sep = "\t", quote = FALSE, row.names = FALSE
      )

      df_path_all <- df_path_all[order(df_path_all$p.adjust), ]

      ## Top Pathway terms
      top_path_list <- list()
      for (cl_id in unique(df_path_all$cell)) {
        tmp <- df_path_all %>%
          filter(cell == cl_id, p.adjust < 0.1) %>%
          head(10)
        if (nrow(tmp) > 0) top_path_list[[cl_id]] <- tmp
      }

      if (length(top_path_list) > 0) {

        top_path <- bind_rows(top_path_list) %>% drop_na(p.adjust)
        top_path$celltype <- factor(top_path$celltype, levels = celltypes)

        ## Dotplot
        p1 <- ggplot(top_path, aes(x = celltype, y = str_wrap(Description, 60))) +
          geom_point(aes(size = GeneRatio, color = p.adjust)) +
          scale_color_gradient(low = "red", high = "blue") +
          ggtitle(var_name) +
          ylab(NULL) + xlab(NULL) +
          theme_bw() +
          theme(
            axis.text.x = element_text(angle = 45, hjust = 1, size = 50),
            axis.text.y = element_text(size = 50),
            text        = element_text(size = 50),
            legend.text = element_text(size = 40),
            legend.key.height = unit(2, "lines"),
            legend.title = element_text(size = 48)
          )

        png(file.path(enrichPathplot,
              paste0(myvar, "_CTRL_top10_enrichPath_DotPlot_loop_DEG50.png")),
            width = 3600, height = 3500, res = 120)
        print(p1)
        dev.off()
      }

    } else {
      cat("No pathway terms to combine for", myvar, "\n")
    }

    cat("Finished processing variable:", myvar, "\n")

  },
  error = function(e) {
    print(paste("END ERROR:", variables[i], ":", e$message))
  },
  finally = {
    print(paste("END FINALLY:", variables[i]))
  })

} # end outer loop

cat("\n\nFinished processing ALL variables.\n")

