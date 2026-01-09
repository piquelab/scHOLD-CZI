# sc-atac-cziHOLD data analysis with CTRL and LPS (exclude DEX)

Use the R version `/wsu/el7/groups/piquelab/R/4.3.2/bin/R`\
Start with **cellranger-atac output** `/rs/rs_grp_schold/CZI/ATAC/counts_cellranger_atac/` \
and **demuxlet output** `/rs/rs_grp_schold/CZI/ATAC/counts_cellranger_atac/fastdemux/fdout/` \
The working path `/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_correct_cellRanger_2025_11_14/` excluding DEX and using up-to-date cellranger-atac results 

Current covariate files in the path `/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/`
- `cp /rs/rs_grp_schold/covariates/dbgap/HOLD_covariates_n165_dbgapIDs_updated_WHR_05_28_2025.txt HOLD-CZI_covariates_Ali_updated_20250528.txt`
- `cp /rs/rs_grp_schold/CZI/RNA/analysis/genotypePCnokin/ALL.eigenvec_pc.txt HOLD-CZI_geno_pc_nokin.txt`
- `cp /rs/rs_grp_schold/covariates/dbgap/HOLD_library_metadata_prelim_n165_dbgapIDs_batch_09_17_2024.txt HOLD-CZI_dbgapIDs_batch_20240917.txt` 

## Folder organization-script for analyzing scATAC-seq data
The folder `1.2_ArchR_process` use the ArchR to analyze sc-atac data from all libraries from CTRL and LPS (exclude LPS-DEX), including qc, cluster, cell type annotation and recall peaks group by cluster,
- `0_demux.R` combine fastdemux results and qc for droplets (SNG, individuals in covariates lists (165) and correct individual_batch), resulting in 798,363 cells. Output-`./0_demux.outs/` 
- `1.2_ArchR_processing.R` create arrow files, merge ArchR project and cluster analysis. ArchR object saved in `./ArchR_all_output/`
- `1_summary_results.R` summarize qc and umap. Output-`./1_results_all/Batch_correct_Harmony/`
- `2.0_generate_scRNA_seurat.R` prepare scRNA-seq data as input data of integration procedure. 
- `2_Integrate_scRNA.R` and `2_Integrate_submit.sh` integrate scRNA with different resolution. Output-`./2_Integrate_output/harmony_default/`
- `2_summary_annot.R` summarize the cell type annotation. 
- `3_reCallPeaks_cluster.R` reCall peaks by the defined cluster
- `3.0_getPeak.R` and `3.0_submit.sh` After reCall peaks, get peak reads matrix for each library. Output-`./3_reCallPeaks_output/option_nFeature15K/`
- `3.1_getPeak.R`, `3.1_submit.sh` and `3.1_summary.R` get peaks reads for each cluster (Cell type)   
- `4b_motif_analysis.R` obtain motif match matrix in the peaks and motif activity. Output-`./4b_motif_output/`
- `4.2_motif_annot.R` obtain motif-peak-gene for each motif and motif-gene across motifs for each cluster (cell type)
- `888_pub_figures.R` and `888_pub_supp.figures.R` plots for manuscript or presentation
  
The folder `2_Differential_analysis` script for differential accessibility and motif activity analysis
- `0_generate_pseudobulk_perLib.R` and `0_submit_perLib.sh` generate pseudobulk ATAC-seq data for Cell type-individual-treatment 
- `0_generate_pseudobulk.R` combine the pseudobulk data across all the libraries and summarize
- `1_differentialPeaks.R` and `1_submit_DESeq.sh` run DESeq2 to perform differential accessibility analysis  
- `1_summary_diff.R` summarize the results of DAR
- `2.0_generate_pseudobulk.R` calculate average motif activity for cell type-individual-treatment
- `2.1_differentialMotifs.R` and `2.1_submit_diff.R` perform differential motif activity analysis associated with psychosocial variables 
- `2.1_summary_diff.R` and `2.1_summary_diff2_LPS.R` summarize the results of DAM associated psychosocial variales in control and LPS treatment
- `2.2_diffMotifs_treat.R` and `2.2_submit_diff_treat` perform differential motif activity between control and treatment (LPS-DAMs)
- `2.2_summary_diff_treat.R` summarize the results of LPS-DAMs
- `3_regulatory2.R` perform enrichment analysis if DEGs are enriched in DAMs
- `3.0_enrich_motif.R` and `3.0_submit_enrich.sh` perform individual motif enrichment analysis
- `888_pub_figure.R` and `888_pub_supp.figures.R` differential plots for manuscript or presentation
 

## Data files used in communication with others, plots or following analysis
Results of processing scATAC-seq 
- `./1.2_ArchR_process/1_results_all/Batch_correct_Harmony/1.0_cluster.rds` cluster results at 0.12 resolution and umap infor
- `./1.2_ArchR_process/2_Integrate_output/scHOLD_RNA_all.outs/sc_rna.hmn13.umap.rds` umap DimReduc objects and `sc_rna.meta.rds` including meta data and umap infor for scRNA-seq datasets
- `./1.2_ArchR_process/3_reCallPeaks_output/Peak_matrix_cluster/1_celltype_peak.txt.gz` 1% cell type peaks
- `./1.2_ArchR_process/4b_motif_output/` contains files including, 
  - `Motif_jaspar2022.activity.mat.rds` for motif activity
  - `Motif_jaspar2022.match.mat.rds` for motif-peak occupancy matrix 
  - `Motif_jaspar2022.motifinfor.ArchR.txt` for motif detailed information
- `./1.2_ArchR_process/4b_motif_output/motif_anno/` contains motif-peak-gene annotation (100 kb to gene start)
  - `1_peak.annotation.rds` for peak-gene annotation
  - `1_motif.infor.txt` contains 692 motif information. 
  - `2_comb_cluster_*_motif.annot.txt.gz` for each cluster include motifs-to-genes 
  - The subfolder `motif_list` contain motif-peak-gene annotation for each motif

Results of differential analysis
- `./2_Differential_analysis/0_pseudobulk/option_nFeature15K_cluster_res0.12/2_YtX_comb.th20.clean.rds` pseudobulk ATAC-seq counts data for cell type-individual-treatment
- `./2_Differential_analysis/1_DiffPeak.outs/summary_option_nFeature15K_cluster_res0.12_psycho_CTRL/2_psycho_plotData.comb.txt.gz` results of DARs associated with psychosocial variable
- `2_Differential_analysis/2_motif.outs/option_nFeature15K_cluster_res0.12/1.3_YtX_ave.th20.clean.rds` average motif activity for cell type-individual-treatment
- `./2_Differential_analysis/2_motif.outs/summary_option_nFeature15K_cluster_res0.12_psycho_CTRL/2_psycho_plotData.comb.txt.gz` results of DAMs associated with psychosocial variable
- `./2_Differential_analysis/2_motif.outs/summary_option_nFeature15K_cluster_res0.12_psycho_treat/2_th20_plotData.comb.txt.gz` results of LPS-DAMs
- `./2_Differential_analysis/3_regulatory.outs/enrich_motif/2_comb_enrich.motif.txt.gz` individual motif enrichment results
- `./2_Differential_analysis/3_regulatory.outs/list_cluster_var_nDEGs_nTF_new.txt` corresponding cluster information for DEG and DAMs  



