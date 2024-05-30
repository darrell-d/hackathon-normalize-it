#!/home/klkoser/miniforge3/envs/cyclone_env/bin/Rscript
# Script to perform batch correction and make checkpoint_1 for input to cyclone for the liver metastases samples from CRC patients. This was altered to run local as it had few enough samples that this is possible.


## By Kelvin Koser

## 5/30/24

############################## Library Packages ################################

suppressMessages({
  library(cyCombine)
  library(data.table)
  library(cowplot)
  library(tidyverse) 
  library(configr)
  library(optparse)
})


# INPUT_DIR=/service/data/input
# OUTPUT_DIR=/service/data/output

###################### Define & parse command-line options #####################
# Define command-line options
option_list <- list(
  make_option(c("-f", "--fcs_file_directory"), type="character", default=NULL,
              help="Path to directory of the FCS files that will be used by cytof.R.", metavar="character"),
  make_option(c("-o", "--output_directory"), type="character", default=NULL,
              help="Path to directory for checkpoint_1.RData output.", metavar="character"),
  make_option(c("-i", "--file_metadata"), type="character", default=NULL,
              help="Path to file_metadata", metavar="character")
)

# Parse command-line options
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
fcs_dir <- Sys.getenv("INPUT_DIR")
output_dir <- Sys.getenv("OUTPUT_DIR")
input_dir <- Sys.getenv("INPUT_DIR")
################### END: Define & parse command-line options ###################

############################ Load metadata files ###############################

#generate variables based on command arguments
file_metadata_path <- paste0(input_dir, "file_metadata.csv")

# Read in file_metadata
file_metadata <- read_csv(file_metadata_path)

# Make sample_id column for downsampling
file_metadata$sample_id <- str_remove(file_metadata$file_name, "\\.fcs$")

########################## End Load metadata files #############################

####################### Read in FCS Files to FlowFrame #########################

# Prepare expression data (read in fcs files into a flow set, extracting expression data)
trans_exp_uncorrected <- prepare_data(data_dir = fcs_dir,
                      metadata = file_metadata,
                      filename_col = "file_name",
                      batch_ids = "pool_id",
                      derand = TRUE,
                      down_sample = FALSE,
                      transform = TRUE,
                      cofactor = 5,
                      seed = 473)

# Save uncorrected flowframe
#save(trans_exp_uncorrected, file = paste0(output_dir, "/trans_exp_uncorrected_rdata.RData"))

print(paste0("Transformed flowframe created at ", Sys.time()))


####################### Pre-Check for Batch Effects ############################

print(paste0("Performing express check for batch effects at ", Sys.time()))

# Use quick function to check if there are batch effects. Output generated in batch_correction_qc folder.
trans_exp_uncorrected %>%
  dplyr::select(!id) %>%
  detect_batch_effect_express(out_dir = paste0(output_dir, "/batch_correction_qc"),
                               seed = 12234)

print(paste0("Express check finished at ", Sys.time(), ". Output in batch_correction_qc folder."))

print(paste0("Performing complete check for batch effects at ", Sys.time()))

# Complete check for batch effects
batch_effect_full <- trans_exp_uncorrected %>%
   detect_batch_effect(batch_col = 'batch',
                       out_dir = paste0(output_dir, '/full_batch_effect_check'), 
                       xdim = 8,
                       downsample = NULL,
                       ydim = 8,
                       seed = 382)

print(paste0("Complete check finished at ", Sys.time(), ". Output in full_batch_effect_check folder."))

##################### End Pre-Check for Batch Effects ##########################

########################## Begin Batch Correction ##############################

print(paste0("Performing batch correction at ", Sys.time()))

# Perform batch correction 
corrected <- trans_exp_uncorrected %>%
  batch_correct(covar = NULL,
                xdim = 8,
                ydim = 8,
                seed = 473,
                norm_method = "scale")

#save
save(corrected, file = paste0(output_dir, "/batch_corrected_raw.RData"))

print(paste0("Batch correction finished at ", Sys.time()))

########################## End Batch Correction ################################

################## Begin QC Evaluation of Batch Correction #####################

print(paste0("Assessing batch correction at ", Sys.time()))

## Make metrics to assess batch correction
labels <- corrected %>%
  cyCombine::create_som(rlen = 10,
                        xdim = 8,
                        ydim = 8)

# Add labels
corrected <- corrected %>%
  dplyr::mutate(som = labels)

# Set column for evaluation of EMD (per-cluster)
celltype_col <- "som"

# Transfer labels to uncorrected data
trans_exp_uncorrected <- corrected %>%
  dplyr::select(id, all_of(celltype_col)) %>%
  dplyr::left_join(trans_exp_uncorrected, by = "id")

# Evaluation using EMD
 emd_val <- trans_exp_uncorrected %>%
   cyCombine::evaluate_emd(corrected,
                           binSize = 0.1,
                           cell_col = celltype_col)

#Make plots
emd_plots <- cowplot::plot_grid(emd_val$violin, emd_val$scatterplot)

 pdf(file = paste0(output_dir, "/full_batch_effect_check/emd_comparison.pdf"), height = 10, width = 10)
 emd_plots
 dev.off()


# Evaluation using MAD
 mad_val <- trans_exp_uncorrected %>%
   cyCombine::evaluate_mad(corrected,
                           filter_limit = NULL,
                           cell_col = "som")


#save mad scores
 print(mad_val$score)
 write_csv(mad_val$mad, file = paste0(output_dir, "/full_batch_effect_check/mad_values.csv"))

#plot distributions per batch comparing original to corrected
 density_plots <- plot_density(trans_exp_uncorrected, corrected, ncol = 4)
 pdf(file = paste0(output_dir, "/full_batch_effect_check/density_distributions.pdf"), height = 20, width = 15)
 density_plots
 dev.off()

# UMAP plots comparisons
 inds <- split(1:length(trans_exp_uncorrected$batch), trans_exp_uncorrected$batch)
 set.seed(6157)
 sample <- unlist(lapply(inds, sample, 10000)) #downsample 7000 per batch
 
 trans_uncorrected_sub <- trans_exp_uncorrected[sample,]
 corrected_sub <- corrected[sample,]
 
 uncorrected_umap <- plot_dimred(trans_uncorrected_sub, name = 'Uncorrected', type = 'umap', return_coord = FALSE)
 uncorrected_umap_coords <- plot_dimred(trans_uncorrected_sub, name = 'Uncorrected', type = 'umap', return_coord = TRUE)$dimred %>%
   as.data.frame()
 corrected_umap <- plot_dimred(corrected_sub, name = 'Corrected 8x8', type = 'umap')
 corrected_umap_coords <- plot_dimred(corrected_sub, name = 'Corrected', type = 'umap', return_coord = TRUE)$dimred %>%
   as.data.frame()
 
 
 umap_plots <- cowplot::plot_grid(uncorrected_umap, corrected_umap)
 pdf(file = paste0(output_dir, "/full_batch_effect_check/umap_comparisons_batch.pdf"), height = 8, width = 10)
 umap_plots
 dev.off()
 
# Add uncorrected umap coords to trans_uncorrected and Plot umap, coloring by "som"
 uncorrected_som_clusters <- trans_uncorrected_sub %>%
   mutate(som = factor(som, levels = seq(1:225))) %>%
   cbind(uncorrected_umap_coords) %>%
   ggplot(aes(x = UMAP1, y = UMAP2, color = som)) +
   geom_point(alpha = .3) +
   theme(legend.position = "none") +
   ggtitle("Uncorrected SOM Clusters (8x8 dim)")
 
 
 
# Add corrected umap coords to corrected and Plot umap, coloring by "som"
 corrected_som_clusters <- corrected_sub %>%
   mutate(som = factor(som, levels = seq(1:225))) %>%
   cbind(corrected_umap_coords) %>%
   ggplot(aes(x = UMAP1, y = UMAP2, color = som)) +
   geom_point(alpha = .3) +
   theme(legend.position = "none") +
   ggtitle("Corrected SOM Clusters (8x8 dim)")
 
 
 # Save a plot of UMAP colored by SOM clusters
 pdf(file = paste0(output_dir, "/full_batch_effect_check/umap_comparisons_som.pdf"), height = 10, width = 15)
 cowplot::plot_grid(uncorrected_som_clusters, corrected_som_clusters)
 dev.off()

print(paste0("QC of batch correction finished at ", Sys.time()))

################### End QC Evaluation of Batch Correction ######################
