### overall questions

# How often is any gene found in multiple traits?

MG <- reshape(
  as.data.frame(
    table(candidateList$`Gene Name`,
          candidateList$trait),
    stringsAsFactors = F
  ),
  idvar = "Var1",
  timevar = "Var2",
  direction = "wide"
)

colnames(MG) <- gsub("Freq.","", colnames(MG))
MG_sums <- data.frame(genes=MG[,1],sums =rowSums(MG[,-1]))

MG_quants <- quantile(MG_sums$sums)

MG_perc <- lapply(c(2:max(MG_sums$sums)), 
                  FUN = function(x){(length(which(MG_sums$sums >= x))/nrow(MG_sums))*100})
names(MG_perc) <- c(2:max(MG_sums$sums))
MG_perc

# How often is a group of trait found more often than others?

# example user selections
traits <- c("Fe","Zn","Mn")
species <- c("Athaliana_TAIR10","Zmays_RefGen_V4","Gmax_Wm82.a2.v1","Sbicolor_v3.1.1","Osativa_v7.0")

MG_subsets <- MG[,c(1,which(colnames(MG) %in% traits))]

print(paste0((MG_perc[[which(names(MG_perc)==paste(length(traits)))]]/100)*nrow(MG_sums),
             " (", MG_perc[[which(names(MG_perc)==paste(length(traits)))]], "%) genes in ", 
             length(traits), " traits or more"))

complete_matches <- MG_subsets[which(rowSums(MG_subsets[,-1])==3),]

print(paste0(nrow(complete_matches),
             " (", (nrow(complete_matches)/nrow(MG_sums))*100, "%) genes found for ",
      paste(traits, collapse = ", ")))

# How often is any orthogroup found in multiple traits?

casted_groups <- data.table::dcast(candidateList, Orthogroup + trait ~ org, fill = 0, value.var = "loci", fun.aggregate = function(x) {
  length(x) > 0
})

MOG <- reshape(
  data.frame(
    table(casted_groups$Orthogroup,
          casted_groups$trait),
    stringsAsFactors = F
  ),
  idvar = "Var1",
  timevar = "Var2",
  direction = "wide"
)

colnames(MOG) <- gsub("Freq.","",colnames(MOG))

MOG_sums <- data.frame(genes=MOG[,1],sums =rowSums(MOG[,-1]))

MOG_quants <- quantile(MOG_sums$sums)

MOG_perc <- lapply(c(2:max(MOG_sums$sums)), 
                  FUN = function(x){(length(which(MOG_sums$sums >= x))/nrow(MOG_sums))*100})
names(MOG_perc) <- c(2:max(MOG_sums$sums))
MOG_perc

# How often is a group of trait found more often than others? (Orthogroup version)

# example user selections
traits <- c("Fe","Zn","Mn")
species <- c("Athaliana_TAIR10","Zmays_RefGen_V4","Gmax_Wm82.a2.v1","Sbicolor_v3.1.1","Osativa_v7.0")

MOG_subsets <- MOG[,c(1,which(colnames(MOG) %in% traits))]

print(paste0((MOG_perc[[which(names(MOG_perc)==paste(length(traits)))]]/100)*nrow(MOG_sums),
             " (", MOG_perc[[which(names(MOG_perc)==paste(length(traits)))]], "%) orthogroups in ", 
             length(traits), " traits or more"))

complete_OGmatches <- MOG_subsets[which(rowSums(MOG_subsets[,-1])==3),]

print(paste0(nrow(complete_OGmatches),
             " (", (nrow(complete_OGmatches)/nrow(MOG_sums))*100, "%) orthogroups found for ",
             paste(traits, collapse = ", ")))

# option 1: All OG in selected traits (default)

OG_output <- casted_groups[casted_groups$Orthogroup %in%
                             MOG_subsets$Var1[(rowSums(MOG_subsets[,-1])>1)]
                           & casted_groups$trait %in% traits,]

# Add a new column for the count of traits

trait_counts <- as.data.frame(table(OG_output$Orthogroup))
colnames(trait_counts) <- c("Orthogroup", "TraitCount")

OG_output <- merge(OG_output, trait_counts, by = "Orthogroup")

# option 2: 

# 

head(table(table$Orthogroup, table$trait))






