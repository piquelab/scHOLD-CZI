# single-cell RNA-seq for HOLD samples - reprocessing

The reprocessing of scRNA-seq samples with the correct CellRanger folder, including the ~40k reads per library /rs/rs\_grp\_schold/CZI/RNA/counts\_cellranger\_2024-04-19. The partial incorrect CellRanger output move to rs/rs\_grp\_schold/CZI/RNA/torm



* GitHub repository cloned in /rs/rs\_grp\_schold/CZI/RNA/CZI-2024-paper
* Working directory /rs/rs\_grp\_schold/CZI/RNA/analysis/
* The old analysis directory renamed to /rs/rs\_grp\_schold/CZI/RNA/analysis\_bk/
* CellRanger folder /rs/rs\_grp\_schold/CZI/RNA/counts\_cellranger\_2024-04-19
* meta-data file with dbgap.ID batch assignments /rs/rs\_grp\_schold/covariates/dbgap/HOLD\_library\_metadata\_prelim\_n165\_dbgapIDs\_batch\_09\_17\_2024.txt





Analysis steps and notes:

* 1\_demux\_alternative.R repeated and outputs compared to analysis\_bk. No difference in the newly generated files
* 1b\_demux\_anal3 ran to remove mismatches from fastdemux outputs
* 2a\_mergeCellRangerAndDemuxlet\_AR.R
