# Load necessary packages
library(tidyverse)
library(circlize)
library(RColorBrewer)
library(gtools)

# Define data directory and file names
data_dir <- "/projectnb/mullenl/alqassar/synteny/circos/"
ref_karyotype <- "karyotype_mel_cinx" 
qry_karyotype <- "karyotype_wcm_pseudochr" 
synteny_data <- "wcm_MelCinx_satsuma_summary.chained.out"

setwd(data_dir)

# Read and process the reference karyotype
Ref_pal <- read.delim(paste0(data_dir, ref_karyotype), sep = "\t", header = TRUE) 
colnames(Ref_pal) <- c("Chr1", "Start", "End", "fill", "species", "size", "color1")

# Extract chromosome numbers and clean data
Ref_pal$Chr <- sapply(strsplit(Ref_pal$Chr1, ' '), function(x) x[3])
Ref_pal <- Ref_pal %>%
  select(Chr, Start, End) %>%
  relocate("Chr", .before = "Start")

# trying to sort dataframe by descending chromosome number
custom_sort <- function(x) {
  x_factor <- factor(x, levels = mixedsort(unique(x), decreasing = TRUE))
  return(order(x_factor))
}
sorted_indices <- custom_sort(Ref_pal$Chr)
Ref_pal <- Ref_pal[sorted_indices, ]

Ref_pal$color <- NULL

# Define custom color palette
x <- c("0 0 0", "173 216 230", "0 191 255", "30 144 255", "0 0 255", "0 0 139", "72 61 139", 
       "123 104 238", "138 43 226", "128 0 128", "218 112 214", "255 0 255", 
       "255 20 147", "176 48 96", "220 20 60", "240 128 128", "255 69 0", 
       "255 165 0", "244 164 96", "240 230 140", "128 128 0", "139 69 19", 
       "255 255 0", "154 205 50", "124 252 0", "144 238 144", "143 188 143", 
       "34 139 34", "0 255 127", "0 255 255", "0 139 139")
color_palette <- sapply(strsplit(x, " "), function(x) rgb(x[1], x[2], x[3], maxColorValue = 255))

Ref_pal <- Ref_pal %>%
  mutate(color = color_palette)

# Read and process the query karyotype
Qry_pal <- read.delim(paste0(data_dir, qry_karyotype), sep = "\t", header = TRUE) 
Qry_pal <- Qry_pal %>%
  select(Chr, Start, End) %>%
  mutate(color = "#FFFFFF")

Qry_pal$Chr <- mixedsort(Qry_pal$Chr)
Qry_pal$Chr <- paste0("Tineola_chr", Qry_pal$Chr)

# Combine reference and query karyotypes
Ref_pal$Chr <- paste0("Melitaea_chr", Ref_pal$Chr)
Pal1 <- rbind(Ref_pal, Qry_pal)

# Read and process the synteny data
synteny_alignment <- read.delim(paste0(data_dir, synteny_data), sep = "\t", header = FALSE)
colnames(synteny_alignment) <- c("ref_chrom_temp", "Start_ref", "End_ref", "qry_chrom_temp", 
                                 "Start_qry", "End_qry", "V6", "V7")

synteny_alignment <- synteny_alignment %>%
  select(-V6, -V7)

# Extract chromosome numbers
synteny_alignment <- synteny_alignment %>%
  mutate(Ref = sapply(strsplit(ref_chrom_temp, '_'), function(x) x[7]),
         Qry = sapply(strsplit(qry_chrom_temp, '_'), function(x) x[3])) %>%
  select(-ref_chrom_temp, -qry_chrom_temp) %>%
  relocate("Ref", .before = "Start_ref") %>%
  relocate("Qry", .before = "Start_qry") %>%
  filter(!grepl('mitochondrion', Ref), !grepl('scaffold', Ref)) %>%
  drop_na() 

#split into two dataframes so that one contains matched chromosome data and the other contains unmatched
synteny_matches <- synteny_alignment %>%
  filter(Ref == Qry)

synteny_unmatched <- synteny_alignment %>%
  filter(Ref != Qry)

#now make into one dataframe
synteny_alignment_reordered <- bind_rows(synteny_matches, synteny_unmatched)

# Add prefixes to chromosome numbers
synteny_alignment_final <- synteny_alignment_reordered %>%
  mutate(Ref = paste0("Melitaea_chr", Ref),
         Qry = paste0("Tineola_chr", Qry))

# Circos visualization
png(filename = "Mel_cinxia_wcm_pseudochr_circos_vis_Sep16.png",
    width = 2000, height = 2000, units = "px")

circos.genomicInitialize(Pal1, plotType = NULL)

circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  chr = CELL_META$sector.index
  xlim = CELL_META$xlim
  circos.rect(xlim[1], 0, xlim[2], 0.5, col = Pal1$color[Pal1$Chr == chr])
  circos.text(mean(xlim), 1, chr, cex = 2, facing = "clockwise",
              niceFacing = TRUE, adj = c(0.4, 0.5))
}, bg.border = NA)

for (i in seq_len(nrow(synteny_alignment_final))) {
  reference <- synteny_alignment_final$Ref[i]
  query <- synteny_alignment_final$Qry[i]
  RS <- synteny_alignment_final$Start_ref[i]
  RE <- synteny_alignment_final$End_ref[i]
  QS <- synteny_alignment_final$Start_qry[i]
  QE <- synteny_alignment_final$End_qry[i]
  
  circos.link(reference, c(RS, RE),
              query, c(QS, QE), col = Pal1$color[Pal1$Chr == reference])
}

circos.clear()
dev.off()

