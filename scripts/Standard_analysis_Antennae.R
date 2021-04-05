setwd("/home/rraborn/scratch/DpTissues_analysis/tsr_antennae")


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
  sample_name=stringr::str_glue("DpTissues_{seq_len(3)}"),
  file_1=NA, file_2=NA,
  condition=rep("antennae", 3)
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

Dp.tss <- list(tss.1.gr, tss.2.gr, tss.3.gr)
names(Dp.tss) <- c("antennae1", "antennae2", "antennae4")

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

##### Sequence analysis
## creating a truncated object for sequence analysis
## some intervals are too close to the edges of short scaffolds, so this was my workaround

plot_threshold_exploration(exp, samples="antennae1", point_size=1) + scale_color_viridis_c()
ggsave(file="threshold_antennae1.png")

plot_threshold_exploration(exp, samples="antennae2", point_size=1) + scale_color_viridis_c()
ggsave(file="threshold_antennae2.png")

plot_threshold_exploration(exp, samples="antennae4", point_size=1) + scale_color_viridis_c()
ggsave(file="threshold_antennae4.png")

## Correlation plot: replicates
plot_correlation(
  exp, data_type="tss",
  use_normalized=TRUE, font_size=12,
  heatmap_colors=viridis::viridis(100)
)

ggsave(file="antennae_correlation_matrix_tss.png")

### Genomic distribution analysis:
#### To update with new UTR-less annotation

plot_genomic_distribution(exp, data_type="tss", samples=c("antennae1", "antennae2","antennae4")) +
  scale_fill_viridis_d(direction=-1, name="Annotation")

ggsave(file="genomic_distribution_tss_antennae.png")

### promoter fraction plot
#### To update with new UTR-less annotation
plot_detected_features(exp, data_type="tss", samples=c("antennae1", "antennae2","antennae4")) +
  scale_fill_viridis_d(direction=-1)

ggsave(file="promoter_fraction_tss_antennae.png")

### Density plot
plot_density(exp, data_type="tss", samples="antennae1")
ggsave(file="TSS_density_CDS_antennae1.png")

plot_density(exp, data_type="tss", samples="antennae2")
ggsave(file="TSS_density_CDS_antennae2.png")

plot_density(exp, data_type="tss", samples="antennae4")
ggsave(file="TSS_density_CDS_antennae4.png")

### TSS pileup heatmap
## need to update annotation/colour
plot_heatmap(
  exp, data_type="tss", samples="antennae1",
  upstream=250, downstream=250,
  use_normalized=TRUE,
  rasterize=TRUE, raster_dpi=150
)

ggsave(file="TSS_pileup_antennae1.png")

plot_heatmap(
  exp, data_type="tss", samples="antennae2",
  upstream=250, downstream=250,
  use_normalized=TRUE,
  rasterize=TRUE, raster_dpi=150
)

ggsave(file="TSS_pileup_antennae2.png")

plot_heatmap(
  exp, data_type="tss", samples="antennae4",
  upstream=250, downstream=250,
  use_normalized=TRUE,
  rasterize=TRUE, raster_dpi=150
)

ggsave(file="TSS_pileup_antennae4.png")

### Sequence logo analysis- all three replicates
plot_sequence_logo(exp, samples="antennae1")
ggsave(file="antennae1_seq_logo.png")

plot_sequence_logo(exp, samples="antennae2")
ggsave(file="antennae2_seq_logo.png")

plot_sequence_logo(exp, samples="antennae4")
ggsave(file="antennae4_seq_logo.png")

### Dinucleotide frequency- all three replicates

plot_dinucleotide_frequencies(exp, samples="antennae1") +
  scale_fill_viridis_c()
ggsave(file="dinucleotide_frequency_antennae1.png")

plot_dinucleotide_frequencies(exp, samples="antennae2") +
  scale_fill_viridis_c()
ggsave(file="dinucleotide_frequency_antennae2.png")

plot_dinucleotide_frequencies(exp, samples="antennae4") +
  scale_fill_viridis_c()
ggsave(file="dinucleotide_frequency_antennae4.png")


### TSS Sequence colour map

plot_sequence_colormap(exp, samples="antennae1", rasterize=TRUE)
ggsave(file="sequence_colormap_antennae1.png")

plot_sequence_colormap(exp, samples="antennae2", rasterize=TRUE)
ggsave(file="sequence_colormap_antennae2.png")

plot_sequence_colormap(exp, samples="antennae4", rasterize=TRUE)
ggsave(file="sequence_colormap_antennae4.png")

### identify TSRs using clustering
exp <- tss_clustering(exp, threshold=3, n_samples=3, max_distance = 25)

# Associate TSSs with TSRs
exp <- associate_with_tsr(exp)

# Annotate TSRs
exp <- annotate_features(exp, data_type="tsr", upstream=250, downstream=100,
                         feature_type="transcript")

# Mark dominant TSS per TSR
exp <- mark_dominant(exp, data_type="tss")


# Calculate TSR metrics
exp <- tsr_metrics(exp)

###### Add tsr analysis from Standard Analysis documentation

## Correlation plot: replicates
plot_correlation(
  exp, data_type="tsr",
  use_normalized=TRUE, font_size=12,
  heatmap_colors=viridis::viridis(100)
)

ggsave(file="antennae_correlation_matrix_tsr.png")

### Genomic distribution analysis:
#### To update with new UTR-less annotation

plot_genomic_distribution(exp, data_type="tsr", samples=c("antennae1", "antennae2","antennae4")) +
  scale_fill_viridis_d(direction=-1, name="Annotation")

ggsave(file="genomic_distribution_tsr_antennae.png")

### promoter fraction plot
#### To update with new UTR-less annotation
plot_detected_features(exp, data_type="tsr", samples=c("antennae1", "antennae2","antennae4")) +
  scale_fill_viridis_d(direction=-1)

ggsave(file="promoter_fraction_tsr_antennae.png")

### Density plot
plot_density(exp, data_type="tsr", samples="antennae1")
ggsave(file="tsr_density_CDS_antennae1.png")

plot_density(exp, data_type="tsr", samples="antennae2")
ggsave(file="tsr_density_CDS_antennae2.png")

plot_density(exp, data_type="tsr", samples="antennae4")
ggsave(file="tsr_density_CDS_antennae4.png")

### tsr pileup heatmap
## need to update annotation/colour
plot_heatmap(
  exp, data_type="tsr", samples="antennae1",
  upstream=250, downstream=250,
  use_normalized=TRUE,
  rasterize=TRUE, raster_dpi=150
)

ggsave(file="tsr_pileup_antennae1.png")

plot_heatmap(
  exp, data_type="tsr", samples="antennae2",
  upstream=250, downstream=250,
  use_normalized=TRUE,
  rasterize=TRUE, raster_dpi=150
)

ggsave(file="tsr_pileup_antennae2.png")

plot_heatmap(
  exp, data_type="tsr", samples="antennae4",
  upstream=250, downstream=250,
  use_normalized=TRUE,
  rasterize=TRUE, raster_dpi=150
)

ggsave(file="tsr_pileup_antennae4.png")

### plot selected tsr metrics
#### TODO: to re-generate using custom script
plot_tsr_metric(exp, tsr_metrics=c("score", "width"), log2_transform=TRUE, samples="antennae1")
ggsave(file="plot_tsr_metrics_antennae1.png")

# Dinucleotide motifs by TSR shape

plot_sequence_logo(exp, dominant=TRUE, samples="antennae1",
                   data_conditions=conditionals(data_grouping=shape_class)
)
ggsave(file="dinucl_motif_plot_shape_antennae1.png")



plot_sequence_logo(exp, dominant=TRUE, samples="antennae2",
                   data_conditions=conditionals(data_grouping=shape_class)
)

ggsave(file="dinucl_motif_plot_shape_antennae2.png")

plot_sequence_logo(exp, dominant=TRUE, samples="antennae4",
                   data_conditions=conditionals(data_grouping=shape_class)
)

ggsave(file="dinucl_motif_plot_shape_antennae4.png")


### Plot selected gene track
#### TODO: Taylor check with Bob about rendering more clearly
gene_tracks(
  exp, feature_name="dp_gene18479",
  samples=c(TSS="antennae1", TSR="antennae1")
)

ggsave(file="gene_track_antennae_gene18479.png")
