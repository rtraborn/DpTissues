setwd("/home/rraborn/scratch/DpTissues_analysis/tsr_combined")

library(TSRchitect)
library(TSRexploreR)
library(GenomicRanges)
library(DESeq2)
library(ChIPseeker)
library(ggseqlogo)
library(ggplot2)

load("/home/rraborn/scratch/DpTissues_analysis/tsr_antennae/PdSTRIPE_complete.RData")

#creating the annotation and assembly files
#update both paths below
Dp.annot <- "/data/LynchLabCME/Daphnia/DaphniaTissues/GoSTRIPES_sing_DpTissues/STRIPES/DpGENOME/PA42.4.0.gff"
Dp.assembly <- "/data/LynchLabCME/Daphnia/DaphniaTissues/GoSTRIPES_sing_DpTissues/STRIPES/DpGENOME/PA42.4.1.fasta"

# Generate sample sheet
sample_sheet <- data.frame(
  sample_name=stringr::str_glue("DpTissues_{seq_len(6)}"),
  file_1=NA, file_2=NA,
  condition=c(rep("antennae", 3),rep("gut",3))
)


#writing the tss files to the workspace
tss.1 <- PdSTRIPE@tssCountData[[1]]
tss.2 <- PdSTRIPE@tssCountData[[2]]
tss.3 <- PdSTRIPE@tssCountData[[3]]

colnames(tss.1) <- c("seq","TSS", "strand", "score")
colnames(tss.2) <- c("seq","TSS", "strand", "score")
colnames(tss.3) <- c("seq","TSS", "strand", "score")

#making granges files from tss data frames
tss.1.gr <- makeGRangesFromDataFrame(tss.1,
                                     keep.extra.columns = TRUE,
                                     seqnames.field="seq",
                                     start.field="TSS",
                                     end.field="TSS",
                                     strand.field = "strand"
)

tss.2.gr <- makeGRangesFromDataFrame(tss.2,
                                     keep.extra.columns = TRUE,
                                     seqnames.field="seq",
                                     start.field="TSS",
                                     end.field="TSS",
                                     strand.field = "strand"
)

tss.3.gr <- makeGRangesFromDataFrame(tss.3,
                                     keep.extra.columns = TRUE,
                                     seqnames.field="seq",
                                     start.field="TSS",
                                     end.field="TSS",
                                     strand.field = "strand"
)

### Now, loading the gut data

load("/home/rraborn/scratch/DpTissues_analysis/tsr_gut/PdSTRIPE_complete.RData")

#writing the tss files to the workspace
tss.4 <- PdSTRIPE@tssCountData[[1]]
tss.5 <- PdSTRIPE@tssCountData[[2]]
tss.6 <- PdSTRIPE@tssCountData[[3]]

colnames(tss.4) <- c("seq","TSS", "strand", "score")
colnames(tss.5) <- c("seq","TSS", "strand", "score")
colnames(tss.6) <- c("seq","TSS", "strand", "score")

#making granges files from tss data frames
tss.4.gr <- makeGRangesFromDataFrame(tss.4,
                                     keep.extra.columns = TRUE,
                                     seqnames.field="seq",
                                     start.field="TSS",
                                     end.field="TSS",
                                     strand.field = "strand"
)

tss.5.gr <- makeGRangesFromDataFrame(tss.5,
                                     keep.extra.columns = TRUE,
                                     seqnames.field="seq",
                                     start.field="TSS",
                                     end.field="TSS",
                                     strand.field = "strand"
)

tss.6.gr <- makeGRangesFromDataFrame(tss.6,
                                     keep.extra.columns = TRUE,
                                     seqnames.field="seq",
                                     start.field="TSS",
                                     end.field="TSS",
                                     strand.field = "strand"
)

Dp.tss <- list(tss.1.gr, tss.2.gr, tss.3.gr, tss.4.gr, tss.5.gr, tss.6.gr)
names(Dp.tss) <- c("antennae1", "antennae2", "antennae4", "gut1", "gut2", "gut3")

#Creating the TSR explorer object
exp <- tsr_explorer(TSSs=Dp.tss, 
                    genome_annotation=Dp.annot, genome_assembly=Dp.assembly,
                    sample_sheet = sample_sheet
)

#Initial TSS processing

exp <- format_counts(exp, data_type="tss")

#Normalize TSSs
exp <- normalize_counts(exp, data_type = "tss", method = "DESeq2")

## TSS annotation
exp <- annotate_features(exp, data_type = "tss", feature_type="transcript")

##### insert threshold exploration


## Correlation plot: replicates
plot_correlation(
  exp, data_type="tss",
  use_normalized=TRUE, font_size=12,
  heatmap_colors=viridis::viridis(100)
)

ggsave(file="antennae_gut_correlation_matrix_tss.png")
