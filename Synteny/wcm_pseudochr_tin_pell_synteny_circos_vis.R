# Script written by Jasmine Alqassar 2024 for circular visulization of synteny alignment data from Satsuma 
# Inspiration taken from: 
# https://github.com/rstewa03/Pieris_macdunnoughii_genome/blob/8eaafd3a8e4153a6015b5e589d3f39b42b595164/SI/Pmac07_Synteny.R/Pmac07_Synteny.R
# https://github.com/annaorteu/Hypolimnas_genome_Wchr_evolution/blob/44b5d865eba05ae2ac8fbd6da54663e1100acc2f/scripts/Satsuma_comparisons_multiple_species_against_each_other_SuppFigure.R

# load necessary packages  
library(tidyverse)
library(circlize)
library(RColorBrewer)
library(gtools)

data_dir <- "/projectnb/mullenl/alqassar/synteny/circos/"
ref_karyotype <- "karyotype_tin_pell" # ref = your comparison genome 
qry_karyotype <- "karyotype_wcm_pseudochr" #qry = your genome 
synteny_data <- "satsuma_summary.chained.out"
setwd(data_dir)

Ref_pal <- read.delim(paste0(data_dir,ref_karyotype), sep = "\t", header = TRUE) 
colnames(Ref_pal) <- c("Chr1", "name", "Start", "End", "fill", "species", "size", "color1")

# extracting just the chromosome numbers and deleting the former columns with full chromosome names
for(i in 1:nrow(Ref_pal)){
  Ref_pal$Chr[i] <-  strsplit(Ref_pal$Chr1[i], ' ') [[1]][3]}


Ref_pal = Ref_pal %>%
  select(Chr,Start,End) %>%
  relocate("Chr", .before = "Start")
  

Ref_pal$Chr <- mixedsort(Ref_pal$Chr, decreasing=TRUE) #to order the chromosomes numerically #added decreasing to make synteny plot that doesn't have crossing lines

Ref_pal$color <- NULL #creating an empty column to fill later


# adding colors to the reference ideogram 

#first create custom color palette and convert to hexadecimal code

x <- c("173 216 230 ", "0 191 255 ", "30 144 255", "0 0 255", "0 0 139", "72 61 139", "123 104 238", "138 43 226", "128 0 128", "218 112 214", "255 0 255", "255 20 147", "176 48 96", "220 20 60", "240 128 128", "255 69 0", "255 165 0", "244 164 96", "240 230 140", "128 128 0", "139 69 19", "255 255 0", "154 205 50", "124 252 0", "144 238 144", "143 188 143", "34 139 34", "0 255 127", "0 255 255", "0 139 139")
color_palette <- sapply(strsplit(x, " "), function(x)
  rgb(x[1], x[2], x[3], maxColorValue=255))

Ref_pal = Ref_pal %>%
  mutate(color = color_palette) #Paired is the name of a colorblindness friendly palette 
#to display colorblind friendly palettes to choose from: display.brewer.all(colorblindFriendly = TRUE) 



# reading in the query genome's karyotype, and adding just a white color to the ideogram for the query genome
Qry_pal <- read.delim(paste0(data_dir,qry_karyotype), sep = "\t", header = TRUE) 

Qry_pal = Qry_pal %>%
  select(Chr,Start,End) 
Qry_pal$color <- NULL #creating an empty column to fill later

Qry_pal = Qry_pal %>%
  mutate(color = "#FFFFFF")

Qry_pal$Chr <- mixedsort(Qry_pal$Chr) #to order the chromosomes numerically

# before visualizing in an ideogram you need to add different prefixes to the chromosomes
Ref_pal[ ,1] = paste0("Tinea_chr", Ref_pal[, 1])
Qry_pal[ ,1] = paste0("Tineola_chr", Qry_pal[, 1])

Pal1 <- union(Ref_pal, Qry_pal) #combine these files 

synteny_alignment <- read.delim(paste0(data_dir, synteny_data), sep = "\t", header = FALSE)
colnames(synteny_alignment) <- c("ref_chrom_temp", "Start_ref", "End_ref", "qry_chrom_temp", 
                                 "Start_qry", "End_qry", "V6", "V7")

synteny_alignment = synteny_alignment %>%
  select(-V6, -V7) 

#extracting just chromosome numbers
for(i in 1:nrow(synteny_alignment)){
  synteny_alignment$Ref[i] <-  strsplit(synteny_alignment$ref_chrom_temp[i], '_')[[1]][7]
  synteny_alignment$Qry[i] <-  strsplit(synteny_alignment$qry_chrom_temp[i], '_')[[1]][3]}

target <- c("matched", "unmatched")

synteny_alignment2 = synteny_alignment %>%
  select(-ref_chrom_temp, -qry_chrom_temp) %>%
  relocate("Ref", .before = "Start_ref") %>%
  relocate("Qry", .before = "Start_qry") %>%
  filter(!grepl('mitochondrion', Ref)) %>%
  filter(!grepl('scaffold', Ref)) %>%
  drop_na() %>%
  mutate(alignment_status=factor(case_when(Ref != Qry ~ 'unmatched', Ref == Qry ~ 'matched'))) %>%
  arrange(factor(alignment_status, levels=target))
    

#synteny_alignment_reorder <- synteny_alignment %>% filter(synteny_alignment$Ref != synteny_alignment$Qry) #hopefully this removes the synteny points with matches to diff chromosomes

#synteny_alignment_matched <- synteny_alignment %>% filter(synteny_alignment$Ref == synteny_alignment$Qry) 

#synteny_alignment_final <- union(synteny_alignment_matched, synteny_alignment_reorder) #now combine the two synteny dataframes, basically appending the synteny points with diff chr matches to the end i.e. drawing them last

#synteny_alignment_final <- rbind(synteny_alignment_matched, synteny_alignment_reorder) #now combine the two synteny dataframes, basically appending the synteny points with diff chr matches to the end i.e. drawing them last


#adding prefixes to chromosome numbers 
synteny_alignment2[ ,1] = paste0("Tinea_chr", synteny_alignment2[, 1])
synteny_alignment2[ ,4] = paste0("Tineola_chr", synteny_alignment2[, 4])

synteny_alignment_final = synteny_alignment2 %>%
  select(-alignment_status)

# circos visualization

#svg(file="wcm_pseudochr_tinpell_circos.svg") #making an empty svg file to be filled 

#circos.initializeWithIdeogram(Pal1,plotType = NULL)

png(filename = "pseudochr_tin_pell_reordered_unique_colors_June3_2024.png",
    width = 2000, height = 2000, units = "px")

circos.genomicInitialize(Pal1,plotType = NULL)
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  chr = CELL_META$sector.index
  xlim = CELL_META$xlim
  ylim = CELL_META$ylim
  circos.rect(xlim[1], 0, xlim[2], 0.5, col = Pal1$color[Pal1$Chr == chr])
  circos.text(mean(xlim), 1, chr, cex = 2, facing = "clockwise",
              niceFacing = TRUE,  adj = c(0.4, 0.5)) }
  , bg.border = NA)


#text(-0.9, -0.8, "Tinea pellionella\ngenome") #you could end labels to the circos plot for where the respective genomes are
#text(0.9, 0.8, "Tineola bisselliella\ngenome")


for ( i in c(1:dim(synteny_alignment_final)[1])){
  reference <- synteny_alignment_final$Ref[i]
  query <- synteny_alignment_final$Qry[i]
  RS <- synteny_alignment_final$Start_ref[i]
  RE <- synteny_alignment_final$End_ref[i]
  QS <- synteny_alignment_final$Start_qry[i]
  QE <- synteny_alignment_final$End_qry[i]
  
  circos.link(reference, c(RS,RE),
              query, c(QS, QE), col = Pal1$color[Pal1$Chr == reference]
  )}


circos.clear() #only really necessary if you want to make another Circos plot right away, but good to do
dev.off() #need this so that the SVG file or any image file is properly written


# then on your computer convert SVG file 
# rsvg-convert -d 1000 -p 1000 file.svg > file.png 
# could also try to make a high quality 1000 x 1000 pixel png from Circos directly maybe




