# Script written by Sara Michele Schaal December 16, 2022
# Script used and edited by Jasmine Alqassar 2023

### INSTALL PACKAGES & LOAD FUNCTIONS
packages_needed <- c("ggplot2", "plotly", "ggpubr", "tidyverse","rsvg", "RIdeogram")

for(i in 1:length(packages_needed)){
  if(!(packages_needed[i] %in% installed.packages())){install.packages(packages_needed[i])}
  library(packages_needed[i], character.only = TRUE)
}

g_legend <- function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

######################################
## DIRECTORIES AND FILE NAMES

DATADIR <- "/projectnb/mullenl/alqassar/synteny/RIdeogram/"
SYNTENYDATA <- "wcm_MelCinx_satsuma_summary.chained.out"
KARYOTYPEDATA <- "karyotype_wcm_MelCinx.txt"

######################################
df_synteny <- read.delim(paste0(DATADIR, SYNTENYDATA), sep = "\t", header = FALSE)
df_karyotype <- read.delim(paste0(DATADIR, KARYOTYPEDATA), sep = "\t", header = TRUE)

colnames(df_karyotype) <- c("Chr", "Start", "End", "fill", "species", "size", "color1")
# change both of the sp. names to be black in karyotype file 
df_karyotype$color <- "000000"

df_karyotype = df_karyotype %>%
select(-color1)
# now want to change the fill of the two species chromosomes 
df_karyotype$fill[df_karyotype$species=='Melitaea_cinxia'] <- 'cccccc'
df_karyotype$fill[df_karyotype$species=='Tineola_bisselliella '] <- 'cccccc'

#change size of sp name 
df_karyotype$size <- '10'

#rename the karyotype species names because they won't fit on the file

df_karyotype <- mutate_if(df_karyotype, 
is.character, 
str_replace_all, 
pattern = "Tineola_bisselliella", 
replacement = "T. bisselliella")

df_karyotype <- mutate_if(df_karyotype, 
is.character, 
str_replace_all, 
pattern = "Melitaea_cinxia", 
replacement = "M. cinxia")

#reorder the karyotype so biselliella is first 
df_karyotype = df_karyotype %>%
  arrange(desc(species))

#remove Z chr from karyotype
df_karyotype = df_karyotype %>%
  mutate(Chr = na_if(Chr, "Z")) %>%
  mutate(Chr = na_if(Chr, "MT"))

df_karyotype <- df_karyotype[complete.cases(df_karyotype),]

#rename the columns 
colnames(df_synteny) <- c("species2_chrom", "Start_2", "End_2", "species1_chrom", 
                          "Start_1", "End_1", "V6", "V7")
#remove columns 6 and 7 and re-order the columns so bombyx is last
df_synteny = df_synteny %>%
  select(-V6, -V7) %>%
  relocate("species2_chrom", .after = "End_1") %>%
  relocate("Start_2", .after = "species2_chrom") %>%
  relocate("End_2", .after = "Start_2")


df_synteny$species1_chrom_num <- NULL #creating an empty column to fill later
df_synteny$species2_chrom_num <- NULL


for(i in 1:nrow(df_synteny)){
  df_synteny$Species_1[i] <-  strsplit(df_synteny$species1_chrom[i], '_')[[1]][3]
  df_synteny$Species_2[i] <-  strsplit(df_synteny$species2_chrom[i], '_')[[1]][7]}


#need to remove unmapped MelCinx scaffold from dataframe 
df_synteny_chroms_only = df_synteny %>% 
  filter(!grepl('HG992240.1', Species_2))


df_synteny_ideogram = df_synteny_chroms_only %>%
  select(Species_1, Start_1, End_1, Species_2, Start_2, End_2) %>% 
  drop_na()


df_synteny_ideogram = df_synteny_ideogram %>%
  add_column(fill = 'cccccc') #here added the new fill column 
#adding colors based on the condition that the chromosome numbers match will be grey and those that don't will have gradient fill
df_synteny_ideogram = df_synteny_ideogram %>%
  mutate(fill=case_when(df_synteny_ideogram$Species_1 == df_synteny_ideogram$Species_2 ~ 'cccccc'
                        ,TRUE ~ '000000')) %>%
  arrange(desc(fill)) #this put the black lines drawn first 

# Now we are going to try to seperate the unique chr:chr matches and provide each different chr a color
#use the search tool to see what is 000000 fill from previous step
df_synteny_ideogram = df_synteny_ideogram %>%
  mutate(fill=case_when(Species_1 == 1 & fill == '000000' ~ 'FF0000', Species_1 == 2 & fill == '000000' ~ 'FF8C00'
                        ,Species_1 == 3 & fill == '000000' ~ 'FFFF00', Species_1 == 4 & fill == '000000' ~ '00FF00'
                        ,Species_1 == 5 & fill == '000000' ~ '0000FF', Species_1 == 6 & fill == '000000' ~ 'E6E6FA'
                        ,Species_1 == 7 & fill == '000000' ~ 'FF69B4', Species_1 == 9 & fill == '000000' ~ 'F08080'
                        ,Species_1 == 10 & fill == '000000'~ 'FFA500', Species_1 == 11 & fill == '000000'~ '2E8B57'
                        ,Species_1 == 12 & fill == '000000'~ '473C8B', Species_1 == 13 & fill == '000000'~ '63B8FF'
                        ,Species_1 == 14 & fill == '000000'~ 'CD8500', Species_1 == 15 & fill == '000000'~ '8B2500'
                        ,Species_1 == 16 & fill == '000000'~ 'FF83FA' , Species_1 == 17 & fill == '000000'~ '98FB98'
                        ,Species_1 == 19 & fill == '000000'~ 'CDAF95', Species_1 == 20 & fill == '000000' ~ '8B2252'
                        ,Species_1 == 21 & fill == '000000' ~ '9ACD32',Species_1 == 22 & fill == '000000' ~ '63B8FF'  
                        ,Species_1 == 23 & fill == '000000'~ '7FFFD4', Species_1 == 24 & fill == '000000' ~'EE7621'
                        ,Species_1 == 25 & fill == '000000' ~ 'FF7256', Species_1 == 26 & fill == '000000' ~ '8B008B'
                        ,Species_1 == 27 & fill == '000000' ~ 'EE9A00', Species_1 == 28 & fill == '000000' ~ 'CD919E'
                        ,Species_1 == 29 & fill == '000000' ~ '551A8B'
                        ,TRUE ~ 'cccccc'))

# Now match the chromosomes that the lines are coming from with the same color
df_karyotype = df_karyotype %>%
  mutate(fill=case_when(Chr == 1 & species == 'T. bisselliella ' ~ 'FF0000', Chr == 2 & species == 'T. bisselliella ' ~ 'FF8C00'
                        ,Chr == 3 & species == 'T. bisselliella ' ~ 'FFFF00', Chr == 4 & species=='T. bisselliella ' ~ '00FF00'
                        ,Chr == 5 & species == 'T. bisselliella ' ~ '0000FF', Chr == 6 & species == 'T. bisselliella ' ~ 'E6E6FA'
                        ,Chr == 7 & species == 'T. bisselliella ' ~ 'FF69B4', Chr == 9 & species == 'T. bisselliella ' ~ 'F08080'
                        ,Chr == 10 & species == 'T. bisselliella '~ 'FFA500', Chr == 11 & species == 'T. bisselliella '~ '2E8B57'
                        ,Chr == 12 & species == 'T. bisselliella '~ '473C8B', Chr == 13 & species == 'T. bisselliella '~ '63B8FF'
                        ,Chr == 14 & species == 'T. bisselliella '~ 'CD8500', Chr == 15 & species == 'T. bisselliella '~ '8B2500'
                        ,Chr == 16 & species == 'T. bisselliella '~ 'FF83FA' , Chr == 17 & species == 'T. bisselliella '~ '98FB98'
                        ,Chr == 19 & species == 'T. bisselliella '~ 'CDAF95', Chr == 20 & species == 'T. bisselliella ' ~ '8B2252'
                        ,Chr == 21 & species == 'T. bisselliella ' ~ '9ACD32', Chr == 22 & species == 'T. bisselliella ' ~ '63B8FF'  
                        ,Chr == 23 & species == 'T. bisselliella '~ '7FFFD4', Chr == 24 & species == 'T. bisselliella ' ~'EE7621'
                        ,Chr == 25 & species == 'T. bisselliella ' ~ 'FF7256', Chr == 26 & species == 'T. bisselliella ' ~ '8B008B'
                        ,Chr == 27 & species == 'T. bisselliella ' ~ 'EE9A00', Chr == 28 & species == 'T. bisselliella ' ~ 'CD919E'
                        ,Chr == 29 & species == 'T. bisselliella ' ~ '551A8B'
                        ,TRUE ~ 'cccccc'))


#now any other non-specified chromosomes will be filled with grey
# df_synteny_chroms_only = df_synteny_chroms_only %>%
# replace(is.na(.), "cccccc")
# df_synteny_chroms_only$fill <- "cccccc" # do this if you want all grey synteny connections 

#Species_1  Start_1    End_1 Species_2 Start_2   End_2   fill
# if we wanted to segregate chromosomes: df_synteny_order <- df_synteny_chromsOnly[c(8,2,3,9,5,6,7)]
# colnames(df_synteny_chroms_only) <- c("Species_1", "Start_1", "End_1", "Species_2", "Start_2",   "End_2",   "fill")
df_synteny_ideogram$Species_1 <- as.numeric(df_synteny_ideogram$Species_1)
df_synteny_ideogram$Species_2 <- as.numeric(df_synteny_ideogram$Species_2)
df_synteny_ideogram <- df_synteny_ideogram[complete.cases(df_synteny_ideogram),]

df_synteny_ideogram = df_synteny_ideogram %>%
  arrange(desc(fill)) #hopefully this draws the colored lines last so they are not hidden behind the grey 

ideogram(karyotype = df_karyotype, synteny = df_synteny_ideogram, output = "wcm_MelCinx_Oct6_color.svg")


# On your computer use rsvg-convert file.svg > file.png
# it seems pdfs are better quality so: rsvg-convert -f pdf file.svg > file.pdf

#if you use the gradient option need to convert to pdf
# rsvg-convert -f pdf file.svg > file.pdf


