# Gathering helpful script chunks for R-shiny figures

### !!!Notes!!! ################################################################

# 1. I have some custom functions in here that work for me, but may not for shiny
# googling alternatives may be wise, but I wouldn't spend too much time on it
# We can focus on them more when something actually breaks

# 2. To keep your rshinydraft.R uncluttered, it might help to put custom functions
# in a separate "source" file. I usually name mine like "FiReMAGE_functions.R"
# then, at the top of your script you can call it in after the packages like this:
# source("rshiny_functions.R")
# source() will read all the functions into the environment so you can use them
# like normal in your plot codes. It also helps if you use the same functions in
# multiple R scripts, so you don't have to copy functions over everytime. An example
# from the rshinydraft.R would be to put our bigPalette object over there. Then
# you would only have to make the subPalette in your draft for x number of traits.
# This is totally optional, just an idea if the script starts to become too big.

# 3. I'll list the packages my code chunks use, but if they aren't needed when 
# you convert the code to rshiny make sure to remove them. We don't want to add 
# unnecessary packages to the app.

################################################################################

# For these I set my working dir to where my github clone was

setwd("C:/Users/lehar/Documents/GitHub/FiReMAGE_Rshiny/FiReMAGE_results/")

### Dot plots from FiReMAGE ####################################################

## Q: How many candidates are returned per species/trait/ortholog group size?

## Default: Looking at all traits at once is messy, see the graphs in FM_references.
# I made them clearer by making a plot per ortholog group size (3/5, 4/5, 5/5).
# For summaries and posters, 3 figures is too much, so I will subset to make a
# comprehensive graph for all species/trait/ortholog group size (see Summary_subset.png).
# I'm not sure what the default should be. Should we have an option for users to
# choose their format with a description for best practices in each? Should we let
# them choose the x and y facets? Right now I have them as x: species, y: ortholog
# group size. Once you have the code converted play around with it to figure out
# the least confusing option for users.

## Options: subset by trait, species, and ortholog group size. Ability to select
# graph layout?

## Required files: AllSummaries.csv

## Packages:

library(data.table)
library(plyr)
library(dplyr)
library(ggplot2)
library(scales) # used for pretty_breaks function in plot, can remove if not needed for plotly

# AllSummaries has gene/loci counts for each species/trait combo in the actual data
# and random permutation data

AllSummaries <- fread("./AllSummaries.csv",
                      header = T,
                      stringsAsFactors = F)

# summarize 1000 permutations with mean/quantiles for each species/trait/present
# type 3 rounds the quantile to the nearest whole values

graphingDF <-
  ddply(
    AllSummaries,
    .(org, trait, dataset, present),
    summarize,
    GenesMean = mean(GeneCount),
    LociMean = mean(LociCount),
    Gene95 = quantile(GeneCount, .95, type = 3),
    Loci95 = quantile(LociCount, .95, type = 3),
    Gene05 = quantile(GeneCount, .05, type = 3),
    Loci05 = quantile(LociCount, .05, type = 3)
  )

# Actual dataset only has 1 observation, so quantiles don't apply

graphingDF$Gene95[graphingDF$dataset == "Actual"] <- NA
graphingDF$Loci95[graphingDF$dataset == "Actual"] <- NA
graphingDF$Gene05[graphingDF$dataset == "Actual"] <- NA
graphingDF$Loci05[graphingDF$dataset == "Actual"] <- NA

# makes it clear to user that present is species present (ie. 3/5 species present)

graphingDF$present <-
  as.factor(gsub("(.*)", "\\1/5", graphingDF$present))

# Summary subset plot version, I'm subsetting to only traits in the 5/5

trait5s <-
  as.character(unique(graphingDF$trait[graphingDF$present == 5 &
                                         graphingDF$dataset == "Actual"]))

graphingDF_5s <- graphingDF[graphingDF$trait %in% trait5s,]

ggplot(graphingDF_5s,
              aes(
                x = trait,
                y = GenesMean,
                group = dataset,
                color = dataset
              )) +
  geom_errorbar(aes(ymax = Gene95,
                    ymin = Gene05,
                    width = 0.75),
                color =  "black") +
  geom_point(aes(fill = dataset,
                 size = dataset),
             color = "black",
             pch = 21) +
  scale_size_manual(values = c(4, 3.5)) +
  labs(y = "orthologous candidate genes") +
  theme_bw(base_size = 14) +
  #theme(panel.spacing = unit(1, "lines")) +
  scale_fill_manual(values = c("#7fcdbb", "#2c7fb8")) +
  scale_y_continuous(breaks = pretty_breaks()) +
  facet_grid(present ~ org,
             scales = "free",
             labeller = label_bquote(col = italic(.(org))))

# All traits by different species present version (Species_combo pngs in FM_references)
# p will print a figure for each level of species present
# in app this subsetting would be taken care of by user's subset choice, right?

for(p in unique(graphingDF$present)) {
  
  # for the for loop I need to use print, but shouldn't be needed in app since we
  # have lapply and other looping options already, like for the boxplot
  
  print(
    ggplot(
      subset(graphingDF, present == p),
      aes(
        x = GenesMean,
        y = trait,
        group = dataset,
        color = dataset
      )
    ) +
      geom_errorbarh(aes(
        xmax = Gene95,
        xmin = Gene05,
        height = 0.75
      )) +
      geom_point(aes(color = dataset, size = dataset)) +
      scale_size_manual(values = c(4, 3)) +
      labs(
        x = "Overlapped Genes Returned",
        title = paste("Orthologs in groups with", p, "species representation")
      ) +
      theme_bw(base_size = 18) +
      theme(plot.title = element_text(hjust = 0.5)) +
      scale_color_manual(values = c("#7fcdbb", "#2c7fb8")) +
      facet_grid(~ org, scales = "free")
  )
  
  # warning messages for removing missing values is normal ggplot output
}

### FDR violin plots ###########################################################

## Q: What's the false discovery rate distribution by species/trait/species present?

## Default: display by species present on x-axis, with labels for avg trait value.
# My current code covers the default, will need to be edited for options.

## Options: subset species and trait in figure. Option to break plots out into a 
# graph for each species so they can see trait FDRs on a species by species level.

## Require files: full candidate list

## Packages:

library(data.table)
library(doParallel)
library(reader)
library(ggplot2)
library(ggrepel)

# my function for combining gene lists across species and traits

concatenate_lists<-function(path_name, pattern=NULL, select=NULL){
  combined<-foreach(i=list.files(path_name, full.names = T, pattern=pattern),
                    .packages = packages.loaded(),
                    .combine = rbind) %do% {
                      
                      l<-data.table::fread(i, sep = ",", header = T, select=select)
                      return(l)
                    }
  return(combined)
}

# combining full candidate lists

candidateList <- concatenate_lists("./candidate_lists/")

# makes it clear to user that present is species present (ie. 3/5 species present)

candidateList$present<-gsub("^([0-9])$","\\1/5",candidateList$present)

# order the levels so 5/5 is on the left

candidateList$present<-factor(candidateList$present, 
                             levels=c("5/5","4/5","3/5"))

# creat trait labels for plot, by trait and species present.
# I wonder if labels will be movable in plotly since it's interactive?

trait_labels <-
  aggregate(candidateList$pFDR,
            list(candidateList$trait,
                 candidateList$present),
            FUN = mean)
colnames(trait_labels) <- c("trait", "present", "pFDR")

# plot code

ggplot(candidateList, aes(x = present,
                          y = pFDR,
                          fill = present)) +
  geom_violin(color = NA) +
  #facet_grid( ~ present) +
  geom_text_repel(
    data = trait_labels,
    aes(label = trait),
    #nudge_x = -0.1,
    size = 4.5,
    nudge_y = 0.08,
    direction = "both",
    min.segment.length = 2,
    force = 3,
    max.overlaps = 40
  ) +
  scale_fill_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb"),
                    guide = "none") +
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1)) +
  ylab("FDR: mean(Permutations)/Actual") + xlab("species in orthogroup")

### Upset figure ###############################################################

## Q: What is the species composition of the ortholog groups with hits in FM?

## Default: display all species and traits selected

## Options: subset display by trait. Possibility to remove 1 or 2 species?
# groups need 3+ species by default, so if we allow to subset species will have 
# to have a message that let's the user know so they understand why the figure
# doesn't work with only 1 species selected.

## Required files: full candidate list

## Packages:

library(data.table)
library(doParallel)
library(reader)
library(ComplexUpset)
library(ggplot2)

# my function for combining gene lists across species and traits

concatenate_lists<-function(path_name, pattern=NULL, select=NULL){
  combined<-foreach(i=list.files(path_name, full.names = T, pattern=pattern),
                    .packages = packages.loaded(),
                    .combine = rbind) %do% {
                      
                      l<-data.table::fread(i, sep = ",", header = T, select=select)
                      return(l)
                    }
  return(combined)
}

# combining full candidate lists

candidateList <- concatenate_lists("./candidate_lists/")

# creating species presence/absence table for each orthogroup and trait

casted_groups <-
  dcast(
    candidateList,
    Orthogroup + trait ~ org,
    fill = 0,
    value.var = "loci",
    fun.aggregate = function(x) {
      length(x) > 0
    }
  )

# plot code, only send true/false columns to upset function
# text and point size may need to be edited for rshiny plot

upset(
  casted_groups[, -c(1, 2)],
  as.character(colnames(casted_groups)[-c(1, 2)]),
  name = "Orthogroup species composition",
  sort_intersections_by = c("degree", "cardinality"),
  sort_sets = FALSE,
  stripes = c("white", "white"),
  matrix = (
    intersection_matrix(geom = geom_point(size = 4),
                        segment = geom_segment(linetype = NA))
  ),
  base_annotations = list(
    'Intersection size' = intersection_size(
      colour = "black",
      text_colors = c("black", "white"),
      text = list(size = 6)
    )
  ),
  themes = upset_default_themes(text = element_text(size = 25))
)
