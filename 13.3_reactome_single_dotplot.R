##
library(tidyverse)
library(purrr)
library(qvalue)
##
library(clusterProfiler)
library(org.Hs.eg.db)
###
library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)
library(ggExtra)
library(RColorBrewer)
library(gtable)
library(ggsignif)
library(pheatmap)
library(corrplot)
library(viridis)
theme_set(theme_grey())

rm(list=ls())

#setwd("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/adjusted/0.2.11/cell20filt/SES_PCs_sex_age_and_treats_generem/figures/GO/updown/down/")
#dataDir <- "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/adjusted/0.2.11/cell20filt/SES_PCs_sex_age_and_treats_generem/deseqres/"
#timestamp()

setwd("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/GO/updown/down/")
dataDir <- "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/deseqres/"

platePrefix <- "ALL.0.1.13." 
plateSufix <- "-RNA-CTRL.noCombat_DESeq" 

timestamp()


#clusters <- c("0", "1", "2", "3", "4", "5")#, "6")
#celltypes <- c("CD4+ T cells 1", "CD4+ T cells 2", "CD8+ T cells", "NK cells", "Monocytes", "B cells")#, "Dentritic Cells")
#celltypes <- c("R0 CD4+ T cell", "R1 CD4+ T cell", "R2 CD8+ T cell", "R3 NK cell", "R4 Monocyte", "R5 B cell")#, "R6 Dentritic Cell")
#celltypes <- c("R0 T CD4+", "R1 T CD4+", "R2 T CD8+", "R3 NK", "R4 Monocyte", "R5 B")#, "R6 Dendritic cell")


clusters <- c("0", "1", "2", "3", "4")#, "5")#, "6")
celltypes <- c("R0 T CD4+", "R1 T CD8+", "R2 NK", "R3 Monocyte", "R4 B")#, "R6 DC")

################################# load down regulated enriched
setwd("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/GO/updown/down/")

variables <- c("PSS_all_mean", "ISEL_Mean")#, "PSS_all_mean", "cytocomp")
variable_names <- c("Psychological Stress","Social Support")#,"Psychological Stress", "Cytokines")


#variables <- c("cytocomp", "PSS_all_mean", "isel")#, "PSS_all_mean", "cytocomp")
#variable_names <- c("Cytokines", "Psychological Stress","Social Support")#,"Psychological Stress", "Cytokines")

#clusters <- c("0", "1", "2", "3", "4", "5")#, "6")
#celltypes <- c("R0 CD4+ T cell", "R1 CD4+ T cell", "R2 CD8+ T cell", "R3 NK cell", "R4 Monocyte", "R5 B cell")#, "R6 Dentritic Cell")
#celltypes <- c("R0 T CD4+", "R1 T CD4+", "R2 T CD8+", "R3 NK", "R4 Monocyte", "R5 B")#, "R6 Dendritic cell")


clusters <- c("0", "1", "2", "3", "4")#, "5")#, "6")
celltypes <- c("R0 T CD4+", "R1 T CD8+", "R2 NK", "R3 Monocyte", "R4 B")#, "R6 DC")

cgdown_list <- list()
cgdown_var <- list()
cgdown_toplist <- list()
cgdown_topvar <- list()

timestamp()
##################### Outer loop to process all variables #####################
for(i in 1:length(variables)){
  #tryCatch({
  myvar <- variables[i]
  var_name <- variable_names[i]
  
  cat("######################################################################\n")
  cat("## Processing: ", var_name, "\n")
               outdir <- paste(getwd(), sep='')

              enrichGODir <- paste(outdir,'/enrichGO', sep='')
              system(paste("mkdir -p",enrichGODir))

              enrichGOdata <- paste(enrichGODir,'/stats/celltype', sep='')
              system(paste("mkdir -p",enrichGOdata))

              enrichGOrds <- paste(enrichGODir,'/data/celltype', sep='')
              system(paste("mkdir -p",enrichGOrds))

              enrichGOplot <-  paste(enrichGODir,'/plots/celltype', sep='')
              system(paste("mkdir -p",enrichGOplot))

              DEGdirGO <- paste(enrichGODir,'/DEGtable/celltype', sep='')
              system(paste("mkdir -p",DEGdirGO))
 
              enrichPathDir <- paste(outdir,'/enrichPath', sep='')
              system(paste("mkdir -p",enrichPathDir))

              enrichPathdata <- paste(enrichPathDir,'/stats/celltype', sep='')
              system(paste("mkdir -p",enrichPathdata))

              enrichPathrds <- paste(enrichPathDir,'/data/celltype', sep='')
              system(paste("mkdir -p",enrichPathrds))

              enrichPathplot <- paste(enrichPathDir,'/plots/celltype', sep='')
              system(paste("mkdir -p",enrichPathplot))

              DEGdirPath<- paste(enrichPathDir,'/DEGtable/celltype', sep='')
              system(paste("mkdir -p",DEGdirPath))


  ################ Inner Loop: Perform erichGO for each cluster ###################
   for(j in 1:length(clusters)){
     cluster <- clusters[j]
     celltype <- celltypes[j]
   
  	#tryCatch({
	 fn <- paste0(enrichPathrds, "/", "C_", cluster, "_", myvar, "_enrichPathway_results.rds")

#GO/updown/down/enrichPath/data/celltype/C_4_ISEL_Mean_enrichGO_results.rds

    # Check if the RDS file exists; if not, skip to the next cluster
    if (!file.exists(fn)) {
      cat("## File not found for cluster: ", cluster, " and variable: ", myvar, "\n")
      next  # Skip to the next iteration (next cluster)
    }
    
    # If the file exists, load it and add it to the cgdown_list
    cgdown_list[[paste0("enrichGO_", cluster, "_", myvar)]] <- readRDS(fn) 
    cgdown_list[[paste0("enrichGO_", cluster, "_", myvar)]]@result$cluster <- paste0("cluster_", cluster)
    cgdown_list[[paste0("enrichGO_", cluster, "_", myvar)]]@result$var <- paste(myvar)
    cgdown_list[[paste0("enrichGO_", cluster, "_", myvar)]]@result$Variable <- paste(var_name)
    cgdown_list[[paste0("enrichGO_", cluster, "_", myvar)]]@result$treat <- "CTRL"
    cgdown_list[[paste0("enrichGO_", cluster, "_", myvar)]]@result$CellType <- paste(celltype)
    #cgdown_list[[paste0("enrichGO_", cluster, "_", myvar)]]@result$updown <- paste0("Down - ", var_name)
    cgdown_list[[paste0("enrichGO_", cluster, "_", myvar)]]@result$updown <- paste0(var_name, " - Down")
    cgdown_list[[paste0("enrichGO_", cluster, "_", myvar)]]@result$varupdown <- paste0(myvar, ".Down")

    topdown <- cgdown_list[[paste0("enrichGO_", cluster, "_", myvar)]]@result[cgdown_list[[paste0("enrichGO_", cluster, "_", myvar)]]@result$p.adjust < 0.1, ]
    topdown <- head(topdown, 10)
    cgdown_toplist[[paste0("top10_", cluster, "_", myvar)]] <- topdown

  }


    cat("## Loading of cluster: ", cluster, " finished\n")
  }

  cgdown_var[[paste0("enrichGO_", myvar)]] <- do.call(rbind, lapply(cgdown_list, function(x) x@result))
  cgdown_topvar[[paste0("top10_", myvar)]] <- do.call(rbind, cgdown_toplist)
  cgdown_topvar[[paste0("top10_", myvar)]] <- cgdown_topvar[[paste0("top10_", myvar)]] %>% drop_na(qvalue)

#}



  cgdown <- do.call(rbind, cgdown_var)
  cgdowntop <- do.call(rbind, cgdown_topvar)




################################# load up regulated enriched
setwd("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/GO/updown/up/")

variables <- c("PSS_all_mean", "ISEL_Mean")#, "PSS_all_mean", "cytocomp")
variable_names <- c("Psychological Stress","Social Support")#,"Psychological Stress", "Cytokines")


#variables <- c("cytocomp", "PSS_all_mean", "isel")#, "PSS_all_mean", "cytocomp")
#variable_names <- c("Cytokines", "Psychological Stress","Social Support")#,"Psychological Stress", "Cytokines")


#clusters <- c("0", "1", "2", "3", "4", "5")#, "6")
#celltypes <- c("R0 CD4+ T cell", "R1 CD4+ T cell", "R2 CD8+ T cell", "R3 NK cell", "R4 Monocyte", "R5 B cell")#, "R6 Dentritic Cell")
#celltypes <- c("R0 T CD4+", "R1 T CD4+", "R2 T CD8+", "R3 NK", "R4 Monocyte", "R5 B")#, "R6 Dendritic cell")


clusters <- c("0", "1", "2", "3", "4")#, "5")#, "6")
celltypes <- c("R0 T CD4+", "R1 T CD8+", "R2 NK", "R3 Monocyte", "R4 B")#, "R6 DC")

cgup_list <- list()
cgup_var <- list()
cgup_toplist <- list()
cgup_topvar <- list()
timestamp()
##################### Outer loop to process all variables #####################
for(i in 1:length(variables)){
  #tryCatch({
  myvar <- variables[i]
  var_name <- variable_names[i]
  
  cat("######################################################################\n")
  cat("## Processing: ", var_name, "\n")
               outdir <- paste(getwd(), sep='')

              enrichGODir <- paste(outdir,'/enrichGO', sep='')
              system(paste("mkdir -p",enrichGODir))

              enrichGOdata <- paste(enrichGODir,'/stats/celltype', sep='')
              system(paste("mkdir -p",enrichGOdata))

              enrichGOrds <- paste(enrichGODir,'/data/celltype', sep='')
              system(paste("mkdir -p",enrichGOrds))

              enrichGOplot <-  paste(enrichGODir,'/plots/celltype', sep='')
              system(paste("mkdir -p",enrichGOplot))

              DEGdirGO <- paste(enrichGODir,'/DEGtable/celltype', sep='')
              system(paste("mkdir -p",DEGdirGO))
 
              enrichPathDir <- paste(outdir,'/enrichPath', sep='')
              system(paste("mkdir -p",enrichPathDir))

              enrichPathdata <- paste(enrichPathDir,'/stats/celltype', sep='')
              system(paste("mkdir -p",enrichPathdata))

              enrichPathrds <- paste(enrichPathDir,'/data/celltype', sep='')
              system(paste("mkdir -p",enrichPathrds))

              enrichPathplot <- paste(enrichPathDir,'/plots/celltype', sep='')
              system(paste("mkdir -p",enrichPathplot))

              DEGdirPath<- paste(enrichPathDir,'/DEGtable/celltype', sep='')
              system(paste("mkdir -p",DEGdirPath))


  ################ Inner Loop: Perform erichGO for each cluster ###################
   for(j in 1:length(clusters)){
     cluster <- clusters[j]
     celltype <- celltypes[j]
   
   #tryCatch({
    fn <- paste0(enrichPathrds, "/", "C_", cluster, "_", myvar, "_enrichPathway_results.rds")
   
    # Check if the RDS file exists; if not, skip to the next cluster
    if (!file.exists(fn)) {
      cat("## File not found for cluster: ", cluster, " and variable: ", myvar, "\n")
      next  # Skip to the next iteration (next cluster)
    }
    
    # If the file exists, load it and add it to the cgup_list
    cgup_list[[paste0("enrichGO_", cluster, "_", myvar)]] <- readRDS(fn) 
    cgup_list[[paste0("enrichGO_", cluster, "_", myvar)]]@result$cluster <- paste0("cluster_", cluster)
    cgup_list[[paste0("enrichGO_", cluster, "_", myvar)]]@result$var <- paste(myvar)
    cgup_list[[paste0("enrichGO_", cluster, "_", myvar)]]@result$Variable <- paste(var_name)
    cgup_list[[paste0("enrichGO_", cluster, "_", myvar)]]@result$treat <- "CTRL"
    cgup_list[[paste0("enrichGO_", cluster, "_", myvar)]]@result$CellType <- paste(celltype)
    #cgup_list[[paste0("enrichGO_", cluster, "_", myvar)]]@result$updown <- paste0("Up - ", var_name)
    cgup_list[[paste0("enrichGO_", cluster, "_", myvar)]]@result$updown <- paste0(var_name, " - Up")
    cgup_list[[paste0("enrichGO_", cluster, "_", myvar)]]@result$varupdown <- paste0(myvar, ".Up")


    topup <- cgup_list[[paste0("enrichGO_", cluster, "_", myvar)]]@result[cgup_list[[paste0("enrichGO_", cluster, "_", myvar)]]@result$p.adjust < 0.1, ]
    topup <- head(topup, 10)
    cgup_toplist[[paste0("top10_", cluster, "_", myvar)]] <- topup

    cat("## Loading of cluster: ", cluster, " finished\n")
  }

  cgup_var[[paste0("enrichGO_", myvar)]] <- do.call(rbind, lapply(cgup_list, function(x) x@result))
  cgup_topvar[[paste0("top10_", myvar)]] <- do.call(rbind, cgup_toplist)
  cgup_topvar[[paste0("top10_", myvar)]] <- cgup_topvar[[paste0("top10_", myvar)]] %>% drop_na(qvalue)

}

  cgup <- do.call(rbind, cgup_var)
  cguptop <- do.call(rbind, cgup_topvar)

####

############################### combined up and down data
#setwd("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/adjusted/0.2.11/SES_PCs_sex_age_and_treats_generem/figures/GO/updown/plots/")
#setwd("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/adjusted/0.2.11/cell20filt/SES_PCs_sex_age_and_treats_generem/figures/GO/updown/plots/")
setwd("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/GO/updown/plots/")

cg <- rbind(cgdown, cgup)
cgtop <- rbind(cgdowntop, cguptop)

table(cg$varupdown)

cgupsig <-  cgup %>% filter(p.adjust < 0.1)
cgdownsig <-  cgdown %>% filter(p.adjust < 0.1)

table(cgupsig$Variable)
table(cgdownsig$Variable)


cgsig <- cg %>% filter(p.adjust < 0.1) #603 total pathways for 2vars and 6celltypes
cgsigord <- cgsig[order(cgsig$p.adjust), ]


write.table(cgsigord, paste0(getwd(), "/", "enriched_pathways_FDR10_allcelltypes_2var_ISEL_updated.txt"), quote = F, row.names = FALSE, sep = "\t")

cgsig$Variable <- as.character(as.vector(cgsig$Variable))
cgsig$CellType <- as.character(as.vector(cgsig$CellType))

cgsig <- as.data.frame(cgsig)

summary_table <- as.data.frame(table(cgsig$Variable, cgsig$CellType))
colnames(summary_table) <- c("Variable", "CellType", "Count")
library(tidyr)
summary_table <- as.data.frame(pivot_wider(summary_table, names_from = CellType, values_from = Count, values_fill = 0))
write.table(summary_table, paste0(getwd(), "/", "enriched_pathways_numb_per_celltype_pervariable_FDR10_ISEL_updated.txt"), quote = F, row.names = FALSE, sep = "\t")

cgsig %>% filter(Description=="Interferon alpha/beta signaling")


table(cgtop$Variable, cgtop$CellType)
table(cgtop$varupdown, cgtop$CellType)
table(cgtop$updown, cgtop$CellType)

variable_names <- c("Psychological Stress", "Social Support")#, "Cytokines")
varsupdown = c("Psychological Stress - Up" , "Social Support - Down")#,  , "Cytokines - Up")

clusters <- c("0", "1", "2", "3", "4")#, "5")#, "6")
celltypes <- c("R0 T CD4+", "R1 T CD8+", "R2 NK", "R3 Monocyte", "R4 B")#, "R6 DC")

cg$updown <- factor(cg$updown, levels = c("Psychological Stress - Up","Psychological Stress - Down", "Social Support - Up", "Social Support - Down"))#,   "Cytokines - Up"))
cg$varupdown <- factor(cg$varupdown, levels = c("PSS_all_mean.Up", "PSS_all_mean.Down", "ISEL_Mean.Up", "ISEL_Mean.Down"))#,   "cytocomp.Up"))
cg$CellType <- factor(cg$CellType, levels = rev(celltypes))

#dd <- meta%>%group_by(orig.ident)%>%
tbl <- cgtop%>%group_by(Description)%>%
             summarise(freq=n()) %>% as.data.frame() %>% arrange(desc(freq))
#pathway_ls <- tbl$Description[1:10]

cgtop <- cgtop %>% arrange(p.adjust)
pathway_ls <- unique(cgtop$Description)[1:10]

#write.table(tbl, paste0(getwd(), "/top10_2var/new/", "top10_pathways_frequencty_allcellypes_ISEL_updated.txt"), quote = F, row.names = FALSE, sep = "\t")
write.table(tbl, paste0(getwd(), "/top10_pathways_frequencty_allcellypes.txt"), quote = F, row.names = FALSE, sep = "\t")

############# make a table reporting the number of genes in each variable, celltype pathway
#library(writexl)
library(openxlsx)
library(dplyr)
library(tidyr)
library(stringr)

cg <- cg %>%
  mutate(ngenes = str_count(geneID, "/") + 1)

# Create a list to store pathway-specific data
pathway_data <- list()

# Clean Excel sheet names
clean_sheet_name <- function(name) {
  name <- str_replace_all(name, "[\\[\\]:*?/\\\\]", "_")  # replace invalid characters
  name <- str_trunc(name, 31, ellipsis = "")              # truncate to 31 characters
  return(name)
}

# Generate pivot tables per pathway
for (pathway in pathway_ls) {
  df_filtered <- cg %>% filter(Description == pathway)
  
  df_filtered <- df_filtered %>%
    mutate(CellType = factor(CellType, levels = celltypes),
           Variable = factor(Variable, levels = variable_names))
  
  df_pivot <- df_filtered %>%
    group_by(CellType, Variable) %>%
    summarise(Total_Genes = sum(Count, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = Variable, values_from = Total_Genes, values_fill = 0)
  
  sheet_name <- clean_sheet_name(pathway)
  pathway_data[[sheet_name]] <- df_pivot
}

# Write to Excel using openxlsx
wb <- createWorkbook()

for (sheet in names(pathway_data)) {
  addWorksheet(wb, sheet)
  writeData(wb, sheet, pathway_data[[sheet]])
}

# Save the workbook
fname <- file.path(getwd(), "top10_pathways_ngenes_summary_fromall_paths.xlsx")
saveWorkbook(wb, fname, overwrite = TRUE)

######################################################

cg$updown <- factor(cg$updown, levels = c("Psychological Stress - Up","Psychological Stress - Down", "Social Support - Up", "Social Support - Down"))#,   "Cytokines - Up"))
cg$varupdown <- factor(cg$varupdown, levels = c("PSS_all_mean.Up", "PSS_all_mean.Down", "ISEL_Mean.Up", "ISEL_Mean.Down"))#,   "cytocomp.Up"))
cg$CellType <- factor(cg$CellType, levels = rev(celltypes))

#cgbackup <- cg

cg <- cg %>% filter(p.adjust < 0.1)
#######################

###
odds.fun <-  function(df){
###    
   res <- map_dfr(1:nrow(df), function(i){
      Diff <- as.numeric(df[i, c("Diff.in", "Diff.not")])
      Bg <- as.numeric(df[i, c("Bg.in", "Bg.not")])
      dat <- data.frame(Diff=Diff, Bg=Bg)
      rownames(dat) <- c("in.category", "not.category")
      fish <- fisher.test(dat)
      res0 <- data.frame(odds=as.numeric(fish$estimate),
                         CI.low=fish$conf.int[1],
                         CI.high=fish$conf.int[2])
      res0
   })
###
  df$odds <- res$odds
  df$CI.low <- res$CI.low
  df$CI.high <- res$CI.high
  df  
}


variable_names <- c("Psychological Stress","Social Support")#,"Psychological Stress", "Cytokines")
var_color <- c("#E34234", "#3A84D8")#, "#FFB6C1")#, "#E34234", "#ADD8E6",) 

group_colors <- c("Psychological Stress" = "#E34234", "Social Support" = "#3A84D8")#, "Shared" = "#2B8C8C")  # you can adjust the Shared color as needed

###
ExampleGOplot <- function(cg){
   cg <- cg%>%drop_na(odds) 
   fig0 <- ggplot(cg, aes(x=updown, y=CellType))+
      #geom_point(aes(size=odds, colour=p2))+
      geom_point(aes(size=odds, color=Variable))+   
      #scale_x_discrete(labels=c("1.LPS"="LPS.Up", "2.LPS"="LPS.Down",
      #   "1.LPS+DEX"="LPS+DEX.Up", "2.LPS+DEX"="LPS+DEX.Down",
      #   "1.PHA"="PHA.Up", "2.PHA"="PHA.Down",
      #   "1.PHA+DEX"="PHA+DEX.Up", "2.PHA+DEX"="PHA+DEX.Down"))+
      #scale_colour_gradient(name="p.adjust",                           
      #   low="blue", high="red", trans="reverse", na.value=NA,
      #   guide=guide_colourbar(order=1), n.breaks=8)+    #"#ffa500"
    scale_color_manual(values = group_colors )+
    #scale_color_manual(values = c("black", "black"))+
         #low="blue", high="red", trans="reverse", na.value=NA,
         #guide=guide_colourbar(order=1), n.breaks=4)+    #"#ffa500"

      scale_size_binned("odds ratio",
         guide=guide_bins(show.limits=TRUE, axis=TRUE,
            axis.show=arrow(length=unit(1.5,"mm"), ends="both"), order=2),
            n.breaks=4)+
      theme_bw()+
      #scale_x_discrete(position = position_dodge2(preserve = "single", padding = 0.2)) +  # Expansion for labels!
      theme(axis.title=element_blank(),
            axis.text.x=element_text(angle=45, size=15, hjust=1, vjust=1),
            #axis.text.x = element_text(angle = -45, hjust = 1, vjust = 0.5, margin = margin(0, 5, 0, 0, "pt")), # Adjust margin values as needed
            axis.text.y=element_text(size=16),
            legend.background=element_blank(),
            legend.title=element_text(size=12),
            legend.text=element_text(size=12),
            legend.key.size=grid::unit(0.5, "lines")
            ) 

   fig0
}


### enrichment for DEGs 
#cgtop <- read_rds("./10_RNA.Variance_output/tmp9/GSE.ClusterProfiler/3_phiNew_enrichGO.rds")
cg <- cg %>%as.data.frame()%>%
   mutate(Diff.in=as.numeric(gsub("/.*","",GeneRatio)),
          Diff.total=as.numeric(gsub(".*/","",GeneRatio)),
          Diff.not=Diff.total-Diff.in,
          Bg.in=as.numeric(gsub("/.*","", BgRatio)),
          Bg.total=as.numeric(gsub(".*/","", BgRatio)),
          Bg.not=Bg.total-Bg.in)
##



#pathway_ls <- tbl$Description[1:10]

 i <- 1

  for (p in pathway_ls){

  ### pathway   
   cg2 <- cg%>%filter(Description==p)
   cg2 <- odds.fun(cg2)
   #cg2 <- cg2%>%full_join(tmp, by=c("Cluster"="rn"))
   cg2$p2 <- cg2$p.adjust
   cg2$p2[cg2$p2>0.1] <- NA


   fig2 <- ExampleGOplot(cg2)+
      #ggtitle(paste(p))+
      ggtitle(str_wrap(p, width = 50)) +
      theme(plot.title=element_text(hjust=0.5, size=17),
            legend.key.size=grid::unit(0.6, "lines"))

    sanitized_p <- gsub("[/ ]", "-", p)

    prefix <- sprintf("%02d", i)

      figfn <- paste0(getwd(), "/v1_top_padj/", prefix, "_figure1_OnePath_", sanitized_p, "_2var", ".png")
      #figfn <- paste0(getwd(), "/v2_top_freq/", prefix, "_figure1_OnePath_", sanitized_p, "_2var", ".png")

  png(figfn, width = 1300, height = 1100, res = 240)
  print(fig2)
  dev.off()

    i <- i + 1
}


######################################## facet wrap 4 dot plots in one

library(ggplot2)
library(dplyr)
library(purrr)
library(stringr)
library(tidyr)
library(forcats)

pathway_ls <- unique(cgsigord[1:30,]$Description)

#subpathway_ls <- c("Interferon Signaling", "Interferon alpha/beta signaling", "Interferon gamma signaling", "ISG15 antiviral mechanism")
subpathway_ls <- pathway_ls[1:4]
subpathway_ls <- pathway_ls[5:8]
subpathway_ls <- pathway_ls[9:12]


# Function to calculate odds ratios and CIs
odds.fun <- function(df) {
  res <- map_dfr(1:nrow(df), function(i) {
    Diff <- as.numeric(df[i, c("Diff.in", "Diff.not")])
    Bg <- as.numeric(df[i, c("Bg.in", "Bg.not")])
    dat <- data.frame(Diff = Diff, Bg = Bg)
    rownames(dat) <- c("in.category", "not.category")
    fish <- fisher.test(dat)
    data.frame(odds = as.numeric(fish$estimate),
               CI.low = fish$conf.int[1],
               CI.high = fish$conf.int[2])
  })
  df$odds <- res$odds
  df$CI.low <- res$CI.low
  df$CI.high <- res$CI.high
  df
}

# Colors and desired facet order
group_colors <- c("Psychological Stress" = "#E34234", "Social Support" = "#3A84D8")
subpathway_ls_ordered <- subpathway_ls


# Prepare data
cg <- cg %>%
  as.data.frame() %>%
  mutate(Diff.in = as.numeric(gsub("/.*", "", GeneRatio)),
         Diff.total = as.numeric(gsub(".*/", "", GeneRatio)),
         Diff.not = Diff.total - Diff.in,
         Bg.in = as.numeric(gsub("/.*", "", BgRatio)),
         Bg.total = as.numeric(gsub(".*/", "", BgRatio)),
         Bg.not = Bg.total - Bg.in)

cg_combined <- map_dfr(subpathway_ls_ordered, function(p) {
  cg2 <- cg %>% filter(Description == p)
  cg2 <- odds.fun(cg2)
  cg2$p2 <- cg2$p.adjust
  cg2$p2[cg2$p2 > 0.1] <- NA
  cg2$Description <- p
  cg2
}) %>%
   mutate(
    Description = factor(Description, levels = subpathway_ls_ordered),
    Description_wrapped = str_wrap(Description, width = 20),
    Description_wrapped = factor(Description_wrapped, 
                                 levels = str_wrap(subpathway_ls_ordered, width = 20))
  )


# Plot
fig_facet <- ggplot(cg_combined, aes(x = updown, y = CellType)) +
  geom_point(aes(size = odds, color = Variable)) +
  scale_color_manual(values = group_colors) +
  scale_size_binned("odds ratio",
    guide = guide_bins(show.limits = TRUE, axis = TRUE,
                       axis.show = arrow(length = unit(1.5, "mm"), ends = "both"), order = 2),
    n.breaks = 4) +
  facet_wrap(~ Description_wrapped, ncol = 2) +
  #facet_wrap(~ Description_wrapped, ncol = 3) +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, size = 12, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 12),
    strip.text = element_text(size = 12, face = "bold"),
    strip.background = element_blank(),  # <- removes gray box
    plot.title = element_text(hjust = 0.5, size = 14),
    legend.background = element_blank(),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.key.size = grid::unit(0.5, "lines"),
    panel.spacing = unit(1, "lines")
  )

# Save plot
#dir.create("v4", showWarnings = FALSE)
png("v2_facet/03_Facet_Subpathways_DotPlot9-12.png", width = 2100, height = 1800, res = 300)
print(fig_facet)
dev.off()

