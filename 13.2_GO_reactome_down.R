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
library(readr)    # for write_rds / read_rds
library(grid)     # for unit() in theme

####################################################
rm(list = ls())

setwd("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/GO/updown/down/")
outFolder <- "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/GO/updown/down/"

dataDir <- "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/deseqres/"

platePrefix <- "ALL.0.1.13."
plateSufix  <- "-RNA-CTRL.noCombat_DESeq"

timestamp()

## Variables to loop over
variables      <- c("PSS_all_mean", "ISEL_Mean")
variable_names <- c("Psychological Stress", "Social Support")

timestamp()

##################### Outer loop to process all variables #####################
for (i in seq_along(variables)) {
  tryCatch({

    myvar    <- variables[i]
    var_name <- variable_names[i]

    cat("######################################################################\n")
    cat("## Processing: ", var_name, "\n")

    outdir <- getwd()

    ## ----- Create directories -----
    enrichGODir   <- file.path(outdir, "enrichGO")
    enrichGOdata  <- file.path(enrichGODir, "stats", "celltype")
    enrichGOrds   <- file.path(enrichGODir, "data",  "celltype")
    enrichGOplot  <- file.path(enrichGODir, "plots", "celltype")
    DEGdirGO      <- file.path(enrichGODir, "DEGtable", "celltype")

    enrichPathDir   <- file.path(outdir, "enrichPath")
    enrichPathdata  <- file.path(enrichPathDir, "stats", "celltype")
    enrichPathrds   <- file.path(enrichPathDir, "data",  "celltype")
    enrichPathplot  <- file.path(enrichPathDir, "plots", "celltype")
    DEGdirPath      <- file.path(enrichPathDir, "DEGtable", "celltype")

    dir.create(enrichGODir,     showWarnings = FALSE, recursive = TRUE)
    dir.create(enrichGOdata,    showWarnings = FALSE, recursive = TRUE)
    dir.create(enrichGOrds,     showWarnings = FALSE, recursive = TRUE)
    dir.create(enrichGOplot,    showWarnings = FALSE, recursive = TRUE)
    dir.create(DEGdirGO,        showWarnings = FALSE, recursive = TRUE)
    dir.create(enrichPathDir,   showWarnings = FALSE, recursive = TRUE)
    dir.create(enrichPathdata,  showWarnings = FALSE, recursive = TRUE)
    dir.create(enrichPathrds,   showWarnings = FALSE, recursive = TRUE)
    dir.create(enrichPathplot,  showWarnings = FALSE, recursive = TRUE)
    dir.create(DEGdirPath,      showWarnings = FALSE, recursive = TRUE)

    ## Initialize list to hold the data frames for the current variable
    dfvar <- list()
    dfsig <- list()

    clusters  <- 0:4
    celltypes <- c("R0 T CD4+", "R1 T CD8+", "R2 NK", "R3 Monocyte", "R4 B")

    ################ Inner Loop: load all gene expression data for each cell type ################
    for (cl in clusters) {

      fname <- paste0(
        dataDir,
        platePrefix, "C", cl,
        ".deseqres_", myvar, plateSufix, ".txt"
      )

      res_data <- read.table(fname, sep = "\t", header = TRUE,
                             quote = '"', comment.char = "")

      dfvar[[as.character(cl)]] <- res_data
      dfsig[[as.character(cl)]] <- res_data %>% filter(padj < 0.1 & logFC < 0)

      cat("## total genes present in cluster", cl, "for", var_name, ":", nrow(dfvar[[as.character(cl)]]), "\n")
      cat("## numb DEGs present in cluster",   cl, "for", var_name, ":", nrow(dfsig[[as.character(cl)]]), "\n")
    }

    ############ How many DEGs? Remove df from the list that are under 50 DEGs ##################
    deg_counts <- data.frame(
      Cluster = character(),
      nGenes  = numeric(),
      DEGs    = numeric(),
      Var     = character(),
      stringsAsFactors = FALSE
    )

    for (cl_name in names(dfvar)) {
      n_genes   <- nrow(dfvar[[cl_name]])
      deg_count <- sum(dfvar[[cl_name]]$padj < 0.1 & dfvar[[cl_name]]$logFC < 0, na.rm = TRUE)

      deg_counts <- rbind(
        deg_counts,
        data.frame(Cluster = cl_name,
                   nGenes = n_genes,
                   DEGs   = deg_count,
                   Var    = myvar)
      )

      if (deg_count < 50) {
        dfvar[[cl_name]] <- NULL
        dfsig[[cl_name]] <- NULL
      }
    }

    output_file <- file.path(DEGdirGO, paste0(myvar, "_number_of_DEGs_per_clusterv2.txt"))
    write.table(deg_counts, file = output_file, sep = "\t",
                row.names = FALSE, col.names = TRUE, quote = FALSE)

    print(deg_counts)
    print(names(dfvar))
    print(names(dfsig))

    ## Initialize lists to hold the enrichment results
    enrichlist    <- list()
    enrich_filt   <- list()
    pathwaylist   <- list()
    pathway_filt  <- list()

    ################ Inner Loop: Perform enrichGO and Reactome for each cluster ###################
    for (cl_name in names(dfvar)) {

      cluster  <- cl_name
      c_index  <- as.numeric(cluster)
      clustname <- celltypes[c_index + 1]

      tryCatch({

        ## Add symbol column
        dfvar[[cluster]]$symbol <- dfvar[[cluster]]$identifier
        dfsig[[cluster]]$symbol <- dfsig[[cluster]]$identifier

        ## Annotation
        anno <- grch38 %>%
          dplyr::select(symbol, entrez) %>%
          distinct()

        annoVar <- anno %>%
          filter(symbol %in% dfvar[[cluster]]$symbol)

        dfvar[[cluster]] <- inner_join(dfvar[[cluster]], annoVar, by = "symbol")
        dfsig[[cluster]] <- inner_join(dfsig[[cluster]], annoVar, by = "symbol")

        ## Drop rows without Entrez IDs
        dfvar[[cluster]] <- dfvar[[cluster]] %>% filter(!is.na(entrez))
        dfsig[[cluster]] <- dfsig[[cluster]] %>% filter(!is.na(entrez))

        if (nrow(dfsig[[cluster]]) == 0) {
          cat("## No DEGs with Entrez IDs for cluster", cluster, "for", var_name, "\n")
          next
        }

        cat("## annotations assigned for cluster", cluster, "for", var_name, "\n")

        ## Prepare gene lists
        genes    <- as.character(dfsig[[cluster]]$entrez)
        geneList <- -log10(dfsig[[cluster]]$pvalue)
        names(geneList) <- dfsig[[cluster]]$entrez
        geneList <- sort(geneList, decreasing = TRUE)

        geneUniv <- as.character(dfvar[[cluster]]$entrez)

        cat("## genes and geneUniv assigned for cluster", cluster, "for", var_name, "\n")

        ## ------ enrichGO ------
        message(".................................")
        message("enrichGO")
        ego <- enrichGO(
          gene      = genes,
          universe  = geneUniv,
          OrgDb     = org.Hs.eg.db,
          ont       = "BP"
        )
        enrichlist[[paste0("enrichGO_", cluster)]] <- ego

        cat("## enrichGO run completed for cluster", cluster, "for", var_name, "\n")

        ## ------ Reactome Pathway ------
        epath <- enrichPathway(
          gene     = genes,
          universe = geneUniv,
          organism = "human"
        )
        pathwaylist[[paste0("pathway_", cluster)]] <- epath

        cat("## enrichPathway run completed for cluster", cluster, "for", var_name, "\n")

        ## ------ Process enrichGO result table ------
        if (is.null(ego) || nrow(ego@result) == 0) {
          cat("No enriched GO terms found for cluster", cluster, "\n")
        } else {
          df_GO_enrich <- ego@result %>% filter(qvalue <= 0.1)

          if (nrow(df_GO_enrich) == 0) {
            cat("No enriched GO terms after qvalue filter for cluster", cluster, "\n")
          } else {
            df_GO_enrich$GeneRatio <- sapply(df_GO_enrich$GeneRatio, function(x) {
              numden <- unlist(strsplit(x, "/"))
              as.numeric(numden[1]) / as.numeric(numden[2])
            })

            cat("## enrichGO filtered for cluster", cluster, "for", var_name, "\n")

            df_GO_enrich$var_cell <- paste0(myvar, "_", cluster)
            df_GO_enrich$cell     <- paste0("cluster_", cluster)
            df_GO_enrich$celltype <- clustname

            enrich_filt[[paste0("df_enrichGO_", cluster)]] <- df_GO_enrich

            cat("Completed enrichGO loop for ", myvar, " ", cluster, "\n")
          }
        }

        ## ------ Process Reactome result table ------
        if (is.null(epath) || nrow(epath@result) == 0) {
          cat("No enriched pathway terms found for cluster", cluster, "\n")
        } else {
          df_GO_pathway <- epath@result %>% filter(qvalue <= 0.1)

          if (nrow(df_GO_pathway) == 0) {
            cat("No enriched pathway terms after qvalue filter for cluster", cluster, "\n")
          } else {
            df_GO_pathway$GeneRatio <- sapply(df_GO_pathway$GeneRatio, function(x) {
              numden <- unlist(strsplit(x, "/"))
              as.numeric(numden[1]) / as.numeric(numden[2])
            })

            cat("## enrichPathway filtered for cluster", cluster, "for", var_name, "\n")

            df_GO_pathway$var_cell <- paste0(myvar, "_", cluster)
            df_GO_pathway$cell     <- paste0("cluster_", cluster)
            df_GO_pathway$celltype <- clustname

            pathway_filt[[paste0("df_pathway_", cluster)]] <- df_GO_pathway

            cat("Completed enrichPathway loop for ", myvar, " ", cluster, "\n")
          }
        }

        ## ------ Save per-cluster RDS results ------
        opfn_go_clust <- file.path(
          enrichGOrds,
          paste0("C_", cluster, "_", myvar, "_enrichGO_results.rds")
        )
        write_rds(ego, opfn_go_clust)

        opfn_path_clust <- file.path(
          enrichPathrds,
          paste0("C_", cluster, "_", myvar, "_enrichPathway_results.rds")
        )
        write_rds(epath, opfn_path_clust)

      },
      error = function(e) {
        print(paste("first mess: Error at iteration for variable", myvar,
                    "in cell type", cluster, ":", e$message))
      },
      finally = {
        print(paste("first mess: Completed iteration for variable", myvar,
                    "in cell type", cluster))
      })

    } ## end for(cl_name in names(dfvar))

    ## ------ Save all-cluster RDS (lists of enrichResult objects) ------
    cat("saving rds objects ", myvar, "\n")

    opfn_go_all <- file.path(enrichGOrds, paste0(myvar, "_enrichGO_results_all_clust.rds"))
    write_rds(enrichlist, opfn_go_all)

    opfn_path_all <- file.path(enrichPathrds, paste0(myvar, "_enrichPathway_results_all_clust.rds"))
    write_rds(pathwaylist, opfn_path_all)

    ## ------ Combine filtered enrichGO results ------
    if (length(enrich_filt) > 0) {

      df_enrich_GO_res <- bind_rows(enrich_filt)
      write.table(df_enrich_GO_res,
                  file = file.path(enrichGOdata, paste0(myvar, "_CTRL_allclust_enrichGO_results_qvalue_0.1_DEG50.txt")),
                  quote = FALSE, sep = "\t", row.names = FALSE)

      df_enrich_GO_res <- df_enrich_GO_res[order(df_enrich_GO_res$p.adjust), ]

      ## Top terms per cluster
      topenrichGO_list <- list()
      for (cl_id in unique(df_enrich_GO_res$cell)) {
        top_cluster <- df_enrich_GO_res %>%
          filter(p.adjust < 0.1, cell == cl_id) %>%
          head(10)
        if (nrow(top_cluster) > 0) {
          topenrichGO_list[[cl_id]] <- top_cluster
        }
      }

      if (length(topenrichGO_list) > 0) {
        topenrichGO_GO <- bind_rows(topenrichGO_list) %>%
          drop_na(p.adjust)

        topenrichGO_GO$celltype <- factor(topenrichGO_GO$celltype, levels = celltypes)

        ## Dotplot for GO
        p1 <- ggplot(topenrichGO_GO,
                     aes(x = celltype, y = str_wrap(Description, width = 60))) +
          geom_point(aes(size = GeneRatio, color = p.adjust)) +
          scale_color_gradient(low = "red", high = "blue", space = "Lab") +
          labs(size = "GeneRatio", color = "p.adjust") +
          ylab(NULL) + xlab(NULL) +
          ggtitle(var_name) +
          theme_bw() +
          theme(
            axis.text.x    = element_text(angle = 45, hjust = 1, size = 50),
            axis.text.y    = element_text(size = 50),
            text           = element_text(size = 50),
            legend.text    = element_text(size = 40),
            legend.key.height = grid::unit(2, "lines"),
            legend.title   = element_text(size = 48, margin = margin(b = 30))
          )

        figfn <- file.path(enrichGOplot,
                           paste0(myvar, "_CTRL_top10_enrichGO_DotPlot_loop_DEG50.png"))
        png(figfn, width = 3600, height = 3500, res = 120)
        print(p1)
        dev.off()
      } else {
        cat("No GO terms passed p.adjust < 0.1 for plotting for", myvar, "\n")
      }

    } else {
      cat("No GO enrichments stored in enrich_filt for", myvar, "\n")
    }

    ## ------ Combine filtered pathway results ------
    if (length(pathway_filt) > 0) {

      df_pathway_GO_res <- bind_rows(pathway_filt)
      write.table(df_pathway_GO_res,
                  file = file.path(enrichPathdata, paste0(myvar, "_CTRL_allclust_enrichPath_results_qvalue_0.1_DEG50.txt")),
                  quote = FALSE, sep = "\t", row.names = FALSE)

      df_pathway_GO_res <- df_pathway_GO_res[order(df_pathway_GO_res$p.adjust), ]

      ## Top terms per cluster
      toppathway_list <- list()
      for (cl_id in unique(df_pathway_GO_res$cell)) {
        top_cluster <- df_pathway_GO_res %>%
          filter(p.adjust < 0.1, cell == cl_id) %>%
          head(10)
        if (nrow(top_cluster) > 0) {
          toppathway_list[[cl_id]] <- top_cluster
        }
      }

      if (length(toppathway_list) > 0) {
        toppathway_GO <- bind_rows(toppathway_list) %>%
          drop_na(p.adjust)

        toppathway_GO$celltype <- factor(toppathway_GO$celltype, levels = celltypes)

        ## Dotplot for pathways
        p1 <- ggplot(toppathway_GO,
                     aes(x = celltype, y = str_wrap(Description, width = 60))) +
          geom_point(aes(size = GeneRatio, color = p.adjust)) +
          scale_color_gradient(low = "red", high = "blue", space = "Lab") +
          labs(size = "GeneRatio", color = "p.adjust") +
          ylab(NULL) + xlab(NULL) +
          ggtitle(var_name) +
          theme_bw() +
          theme(
            axis.text.x    = element_text(angle = 45, hjust = 1, size = 50),
            axis.text.y    = element_text(size = 50),
            text           = element_text(size = 50),
            legend.text    = element_text(size = 40),
            legend.key.height = grid::unit(2, "lines"),
            legend.title   = element_text(size = 48, margin = margin(b = 30))
          )

        figfn <- file.path(enrichPathplot,
                           paste0(myvar, "_CTRL_top10_enrichPath_DotPlot_loop_DEG50.png"))
        png(figfn, width = 3600, height = 3500, res = 120)
        print(p1)
        dev.off()
      } else {
        cat("No pathway terms passed p.adjust < 0.1 for plotting for", myvar, "\n")
      }

    } else {
      cat("No pathway enrichments stored in pathway_filt for", myvar, "\n")
    }

    cat("Finished processing for ", myvar, "\n")

  },
  error = function(e) {
    print(paste("end mess: Error at outer loop for variable", variables[i], ":", e$message))
  },
  finally = {
    print(paste("end mess: Completed outer loop for variable", variables[i]))
  })
}

cat("Finished processing all variables.\n")

