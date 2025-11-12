# single-cell RNA-seq for HOLD samples - reprocessing

The reprocessing of scRNA-seq samples with the correct CellRanger folder, including the ~40k reads per library `/rs/rs\_grp\_schold/CZI/RNA/counts\_cellranger\_2024-04-19`. The partial incorrect CellRanger output move to `rs/rs\_grp\_schold/CZI/RNA/torm`

* GitHub repository cloned in `/rs/rs\_grp\_schold/CZI/RNA/CZI-2024-paper`
* Working directory `/rs/rs\_grp\_schold/CZI/RNA/analysis/`
* The old analysis directory renamed to `/rs/rs\_grp\_schold/CZI/RNA/analysis\_bk/`
* CellRanger folder `/rs/rs\_grp\_schold/CZI/RNA/counts\_cellranger\_2024-04-19`
* meta-data file with dbgap.ID batch assignments `/rs/rs\_grp\_schold/covariates/dbgap/HOLD\_library\_metadata\_prelim\_n165\_dbgapIDs\_batch\_09\_17\_2024.txt`


Analysis steps and notes:

* `1\_demux\_alternative.R` repeated and outputs compared to analysis\_bk. No difference in the newly generated files
* `1b\_demux\_anal3` ran to remove mismatches from fastdemux outputs
* `2a\_mergeCellRangerAndDemuxlet\_AR.R`
* `2a` script splitted into `2a1_mergeCellRangerAndDemuxlet_AR.R` and `2a2_mergeCellRangerAndDemuxlet_AR.R`
* `2b_seuratPCA.R` was ran with modification to `ScaleData` to include only the variab; features for scaling: `ScaleData(sc, features = VariableFeatures(object = sc))` ## RPR AR: we think we only need to do it for var features
* `2c_seuratumap.R` includes the loop for dim 13 and 50
* `2d_seuratclustering.R` this script loops through `FindClusters` for multiple resolutions. It was ran twice separately, once for dim13 and once for dim50
* `2e_seuratpreharmony_umap.R` skipped

