options(scipen = 999)
library(ggplot2)
library(readr)
library(dplyr)
library(optparse)

source("~/software/malaria-hub/utils/helpers.R")

option_list = list(
  make_option(c("--workdir"), type = "character", default = NULL,
              help = "Specify main directory",
              metavar = "character"),
  make_option(c("--country_label"), type = "character", default = "country",
              help = "Specify country label from metadata",
              metavar = "character"),
  make_option(c("-m", "--metadata_file"), type = "character",
              help = "Absolute path to metadata file",
              metavar = "character"),            
  make_option(c("--region_label"), type = "character", default = "region",
              help = "Specify region label from metadata",
              metavar = "character"),
  make_option(c("--specific_country"), type = "character", default = NULL,
              help = "Use if only wanting to plot hmmIBD output for a specific country",
              metavar = "character"),
  make_option(c("--specific_region"), type = "character", default = NULL,
              help = "Corresponding region of country in --specific_country",
              metavar = "character"),            
  make_option("--gene_product", type = "character",
              help = "Gene product file",
              metavar = "character"),
  make_option(c("-r", "--ref_index"), type = "character",
              help = "File name for reference index",
              metavar = "character"),
  make_option("--suffix", type = "character", default = format(Sys.time(), "%d_%m_%Y"),
            help = "Suffix for output files",
            metavar = "character"),
  make_option("--remove_chr", type = "character", default = "PvP01_API_v1,PvP01_MIT_v1",
              help = "Chromosomes to remove ex. Pf3D7_API_v3,Pf_M76611",
              metavar = "character"),
  make_option("--regex_chr", type = "character", default = "(.*?)_(.+)_(.*)",
              help = "Regex pattern for chromosome detection. Default matches Pf3D7_01_v3",
              metavar = "character"),
  make_option("--regex_groupid", type = "numeric", default = 3,
              help = "Regex pattern group",
              metavar = "numeric"),
  make_option("--output_name_prefix", type = "character", default = "output",
              help = "Prefix for .pdf output",
              metavar = "character")
);

# Parse the command line arguments
opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

# Workdir
workdir <- opt$workdir
# Country label
country_label <- opt$country_label
# Region label 
region_label <- opt$region_label
# Metadata
metadata_file <- opt$metadata_file
# Specific country to plot
specific_country <- opt$specific_country
# Region of specific country to 
specific_region <- opt$specific_region
# Annotation file
gene_product_file <- opt$gene_product
# Ref index
ref_index <- opt$ref_index
# Suffix
suffix <- opt$suffix
# Remove chromosomes
rm_chr <- opt$remove_chr
# Pattern for chromosome detection
pattern <- opt$regex_chr
# Pattern group
groupid <- opt$regex_groupid
# Output .pdf prefix name
output_prefix <- opt$output_name_prefix

# Loading IBD and IBD fractions (Combined Country Output)
combined_ibd_r <- read_tsv(file.path(workdir, sprintf("%s_hmmIBD_ibd_results_combined.tsv", suffix)), col_types = cols())
fraction_ibd_r <- read_tsv(file.path(workdir, sprintf("%s_hmmIBD_fraction_results_combined.tsv", suffix)), col_types = cols())

# Loading Metadata
metadata <- read_tsv(metadata_file, col_types = cols()) %>%
  select(all_of(country_label), all_of(region_label))
# Reference Index
fai <- read.table(ref_index, stringsAsFactors = FALSE) %>%
  rename(chr = V1, end_chr = V2) %>%
  select(chr, end_chr) %>%
  mutate(start_chr = 1) %>%
  filter(!chr %in% rm_chr)

fai$chr <- as.numeric(stringr::str_match(fai$chr, pattern)[, groupid])
chrom_ends <- c(0, fai$end_chr[-max(NROW(fai$end_chr))])
transpose_chr <- data.frame(chr = fai$chr, tr_chr = chrom_ends) %>%
                     mutate(ind = seq(1, nrow(.)))

# Check if --specific_country and --specific_region are provided
if (!is.null(opt$specific_country) && !is.null(opt$specific_region)) {
  # Specific country and region are provided
  specific_country_ibd_r <- subset(combined_ibd_r, category %in% opt$specific_country)
  specific_country_ibd_r <- specific_country_ibd_r %>%
    left_join(metadata, by = c("country" = country_label, "region" = region_label))
  
  # Subsetting Metadata at Country-Level if --specific_country provided
  specific_country_ibd_r <- specific_country_ibd_r %>% left_join(metadata, by = c("category" = country_label))
  
  # Use specific_country_ibd_r for subsequent operations
  ibd_frac_tr <- specific_country_ibd_r %>% group_by(chr) %>%
    mutate(trans = get_chrom_transposition(transpose_chr, chr)) %>%
    mutate(pos_bp_ed = as.numeric(as.numeric(start) + trans),
           Fraction = as.numeric(fraction)) %>%
    ungroup()
} else {
  # No specific country and region provided, use combined output
    #Combine results with region
  combined_ibd_r <- combined_ibd_r %>% left_join(metadata, by = c("category" = country_label))

  ibd_frac_tr <- combined_ibd_r %>% group_by(chr) %>%
    mutate(trans = get_chrom_transposition(transpose_chr, chr)) %>%
    mutate(pos_bp_ed = as.numeric(as.numeric(start) + trans),
           Fraction = as.numeric(fraction)) %>%
    ungroup()
}

# Establish order in plot
# Arrange according to order 
ibd_frac_tr_gg <- ibd_frac_tr %>% mutate(region = factor(!!sym(region_label))) %>% arrange(region)
ibd_frac_tr_gg <- ibd_frac_tr_gg %>% mutate(category = factor(category, levels = unique(category))) %>% as.data.frame()

# Plotting IBD pairwise fraction in 10kb windows
fraction_plot <- ggplot(data = ibd_frac_tr_gg) +
    geom_line(aes(x = pos_bp_ed, y = fraction, color = region)) +
    scale_y_continuous(limits = c(0, 0.2), breaks = c(0, 0.2), labels = c("0.0", "0.2")) +
    facet_grid(category ~ ., space = "free_x") +
    labs(x = "Chromsome", y = "IBD Fraction") +
    guides(color = guide_legend(title = "Region:", nrow = 2, byrow = FALSE)) +
    scale_color_manual(values = c(
      "East_Africa" = "#E87D72",
      "South_America" = "#B79F00",
      "South_Asia" = "#00BA38",
      "South_East_Asia" = "#56BCC2",
      "Southern_SEA" = "#6F9AF8"
    )) + 
    theme_classic() +
    theme(axis.line.x = element_line(color = "black"),
          axis.line.y = element_line(color = "black"),
          axis.text.x = element_text(size = 8, color = "black", angle = 0, vjust = -0.5),
          axis.text.y = element_text(size = 8, color = "black"),
          axis.title.x = element_text(size = 8, color = "black"),
          axis.title.y = element_text(size = 8, color = "black"),
          plot.title = element_text(size = 15, color = "black", hjust = 0.5),
          strip.placement = "outside",
          strip.text.y = element_text(angle = 0, face = "bold", size = 9),
          strip.background = element_blank(),
          legend.position = "bottom") +
    geom_vline(data = ibd_frac_tr_gg, aes(xintercept = trans), color = "black",
               alpha = 0.5, linetype = "longdash", linewidth = 0.1) +
    scale_x_continuous(breaks = unique(ibd_frac_tr_gg$trans),
                       labels = unique(ibd_frac_tr_gg$chr))
if (!is.null(opt$output_name_prefix)) {
  if (!is.null(opt$specific_country) && !is.null(opt$specific_region)) {
    ggsave(paste0(opt$output_name_prefix, "_pairwise_fractions_10kb_windows.pdf"), plot=fraction_plot, dpi=600, width=12, height=7)
  } else {
    ggsave(paste0(opt$output_name_prefix, "_pairwise_fractions_10kb_windows.pdf"), plot=fraction_plot, dpi=600)
  }
}
