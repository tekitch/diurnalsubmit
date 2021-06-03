# needs to extract the eic and profile, write an mzrt file, do it per protein. 
library(RPostgreSQL)
library(dplyr)
library(tidyr)
library(stringr)
library(seqinr)
se <- function(x) sd(x)/sqrt(length(x))

setwd("d:/diurnalfive")

databaseName = "diurnalfive"
daydatanames = c("eod5r1atdto", "eod5r2atdto", "eod5r3atdto")
nightdatanames = c("eon5r1atdto", "eon5r2atdto", "eon5r3atdto")
proteinInfonames = c("eodidatdproteins")

#should run a peptide prophet with all of the pepxmls, labelled and unlabelled
fastaFile = "at.fasta"
anno = read.delim("anno.txt", sep = "\t")
anno$IDENTIFIER = toupper(anno$IDENTIFIER)
anno$IDENTIFIER = gsub("\\..*","",anno$IDENTIFIER)

drv <- dbDriver("PostgreSQL")

con <- dbConnect(drv, dbname = databaseName,
                 host = "localhost", port = 5432,
                 user = "postgres", password = "postgres")
sql1 = "select * from "
sql2 = ' WHERE "nnlsDeviance" < \'.0009\';'

daydata = list()
nightdata = list()
proteinInfo = list()

for(i in 1:length(daydatanames)){
  daydata[[i]] = dbGetQuery(con, paste(sql1, daydatanames[i], sql2, sep = ""))
  nightdata[[i]] = dbGetQuery(con, paste(sql1, nightdatanames[i], sql2, sep = ""))
}

names(daydata) = daydatanames
names(nightdata) = nightdatanames
proteinInfo = dbGetQuery(con, paste(sql1, proteinInfonames, sep = ""))
dbDisconnect(con)


#########
#GET PROTEIN PROPHET PEPTIDE ASSIGNMENTS
#########

protGroupVec = c()
protGroupCount = c()

for(i in 1:length(daydata)){
  protGroupVec = c()
  protGroupCount = c()
  for(j in 1:nrow(daydata[[i]])){
    pepseq = daydata[[i]][j,]$peptide
    proteinIndexVec = grep(pepseq, proteinInfo$unipepseq)
    groupFks = c()
    groupProbs = c()
    for(k in 1:length(proteinIndexVec)){
      groupFks[k] = proteinInfo[proteinIndexVec[k],]$group_fk
    }
    groupFks = unique(groupFks)
    protGroupVec[j] = paste0(groupFks, sep = "", collapse = ";")
    protGroupCount[j] = length(groupFks)
  }
  daydata[[i]]$protGroups = protGroupVec
  daydata[[i]]$protGroupCount = protGroupCount
}

for(i in 1:length(nightdata)){
  protGroupVec = c()
  protGroupCount = c()
  for(j in 1:nrow(nightdata[[i]])){
    pepseq = nightdata[[i]][j,]$peptide
    proteinIndexVec = grep(pepseq, proteinInfo$unipepseq)
    groupFks = c()
    groupProbs = c()
    for(k in 1:length(proteinIndexVec)){
      groupFks[k] = proteinInfo[proteinIndexVec[k],]$group_fk
    }
    groupFks = unique(groupFks)
    protGroupVec[j] = paste0(groupFks, sep = "", collapse = ";")
    protGroupCount[j] = length(groupFks)
  }
  nightdata[[i]]$protGroups = protGroupVec
  nightdata[[i]]$protGroupCount = protGroupCount
}

proteinVec = unique(proteinInfo$group_fk)

#########
#FILTERS#
#########
#daydata = subset(daydata, daydata$protGroupCount == 1)
#nightdata = subset(nightdata, nightdata$protGroupCount == 1) #need to make this more sophisticated by sorting peptides to their most likely group by probability
for(i in 1:length(daydata)){
daydata[[i]] = subset(daydata[[i]], daydata[[i]]$unlabToSynthZero > 0.97)
nightdata[[i]] = subset(nightdata[[i]], nightdata[[i]]$unlabToSynthZero > 0.97)
daydata[[i]] = subset(daydata[[i]], daydata[[i]]$LPFsum < 1.2) #add up all of the nnls pops
nightdata[[i]] = subset(nightdata[[i]], nightdata[[i]]$LPFsum < 1.2)
daydata[[i]] = subset(daydata[[i]], daydata[[i]]$fintensityapex > 50000)
nightdata[[i]] = subset(nightdata[[i]], nightdata[[i]]$fintensityapex > 50000)
daydata[[i]] = subset(daydata[[i]], daydata[[i]]$nnlsDeviance < .0009)
nightdata[[i]] = subset(nightdata[[i]], nightdata[[i]]$nnlsDeviance < .0009)
daydata[[i]] = subset(daydata[[i]], daydata[[i]]$isomassd == 0)
nightdata[[i]] = subset(nightdata[[i]], nightdata[[i]]$isomassd == 0)
daydata[[i]] = subset(daydata[[i]], lastpop < .03)
nightdata[[i]] = subset(nightdata[[i]], lastpop < .03)
daydata[[i]] <- daydata[[i]][order(daydata[[i]]$peptide, -abs(daydata[[i]]$fintensitysum) ), ] #sort by id and reverse of abs(value)
daydata[[i]] <- daydata[[i]][ !duplicated(daydata[[i]]$peptide), ]# take the first row within each id
nightdata[[i]] <- nightdata[[i]][order(nightdata[[i]]$peptide, -abs(nightdata[[i]]$fintensitysum) ), ]
nightdata[[i]] <- nightdata[[i]][ !duplicated(nightdata[[i]]$peptide), ]
daydata[[i]] = daydata[[i]] %>% group_by(protGroups) %>% top_n(n = 3, wt = fintensitysum)
nightdata[[i]] = nightdata[[i]] %>% group_by(protGroups) %>% top_n(n = 3, wt = fintensitysum)
}

daydatamelt = bind_rows(daydata, .id = "column_label")
nightdatamelt = bind_rows(nightdata, .id = "column_label")
results = as.data.frame(matrix(ncol = 15, nrow = length(proteinVec)))
colnames(results) = c("protein", "seqsource", "seqsource2", "daylpf", "daysd", "dayrsd", "dayse", "nightlpf", "nightsd", "nightrsd", "nightse", "daypeptides", "nightpeptides", "unidaypeptides", "uninightpeptides")

for(j in 1:length(proteinVec)){
  protein = as.numeric(proteinVec[j])
  daysubset = subset(daydatamelt, daydatamelt$protGroups == protein)
  nightsubset = subset(nightdatamelt, nightdatamelt$protGroups == protein)
  seqsource = daysubset[1,]$seqsource
  seqsource2 = nightsubset[1,]$seqsource
  dayMeanLPF = mean(daysubset$LPF)
  daySD = sd(daysubset$LPF)
  dayRSD = daySD/dayMeanLPF
  dayse = se(daysubset$LPF)
  nightMeanLPF = mean(nightsubset$LPF)
  nightSD = sd(nightsubset$LPF)
  nightRSD = nightSD/nightMeanLPF
  nightse = se(nightsubset$LPF)
  daypeptides = length(daysubset$peptide)
  nightpeptides = length(nightsubset$peptide)
  unidaypeptides = length(unique(daysubset$peptide))
  uninightpeptides = length(unique(nightsubset$peptide))
  
  results[j,]$protein = protein
  results[j,]$seqsource = seqsource
  results[j,]$seqsource2 = seqsource2
  results[j,]$daylpf = round(dayMeanLPF, digits = 3)
  results[j,]$daysd = round(daySD, digits = 3)
  results[j,]$dayrsd = round(dayRSD, digits = 3)
  results[j,]$dayse = round(dayse, digits = 3)
  results[j,]$nightlpf = round(nightMeanLPF, digits = 3)
  results[j,]$nightsd = round(nightSD, digits = 3)
  results[j,]$nightrsd = round(nightRSD, digits = 3)
  results[j,]$nightse = round(nightse, digits = 3)
  results[j,]$daypeptides = daypeptides
  results[j,]$nightpeptides = nightpeptides
  results[j,]$unidaypeptides = unidaypeptides
  results[j,]$uninightpeptides = uninightpeptides
}

#formatting the results table - removing total NA rows, annotation 
for(i in 1:length(results$seqsource)){
  if(is.na(results$seqsource[i])){results$seqsource[i] = results$seqsource2[i]}
}
results = results[!is.na(results$seqsource),]
results$seqsource2 = NULL
results$AGI = gsub("\\..*","",results$seqsource)
anno$DESCRIPTION = NULL
anno$TYPE = NULL
results = merge(results, anno, by.x = "AGI", by.y = "IDENTIFIER")
results$BINCODE = as.character(results$BINCODE)
onebin = c()
for(i in 1:length(results$BINCODE)){
  onebin[i] = paste0(str_split(results[i,]$BINCODE, "^*[.]", n=3, simplify = TRUE)[1])
}
results$onebin = onebin
twobins = c()
for(i in 1:length(results$BINCODE)){
  twobins[i] = paste0(str_split(results[i,]$BINCODE, "^*[.]", n=3, simplify = TRUE)[1], ".", str_split(results[i,]$BINCODE, "^*[.]", n=3, simplify = TRUE)[2])
}
results$twobins = twobins
twobinnames = c()
results$NAME = as.character(results$NAME)
for(i in 1:length(results$NAME)){
  twobinnames[i] = paste0(str_split(results[i,]$NAME, "^*[.]", n=3, simplify = TRUE)[1], ".", str_split(results[i,]$NAME, "^*[.]", n=3, simplify = TRUE)[2])
}
results$twobinnames = twobinnames

# these resutls are in either day and/or night 
paperresults = subset(results, (results$daysd < 0.06 | results$dayrsd < 0.25) & results$unidaypeptides >2 & results$daypeptides > 4 | (results$nightsd < 0.06 | results$nightrsd < 0.25) & results$uninightpeptides > 2 & results$nightpeptides > 4)

for(i in 1:nrow(paperresults)){
  if(!((paperresults[i,5] < 0.06 | paperresults[i,6] < 0.25) & paperresults[i,14] > 2 & paperresults[i,12] > 4)){
    paperresults[i,c(4,5,6,7,12,14)] = NA
  }
  if(!((paperresults[i,9] < 0.06 | paperresults[i,10] < 0.25) & paperresults[i,15] > 2 & paperresults[i,13] > 4)){
    paperresults[i,c(8,9,10,11,13,15)] = NA
  }
}
paperresults = paperresults[!duplicated(paperresults$AGI),]


# these results occur in both the day and the night with sensible sd or rsd and at least 3 different peptide sequences measured at least 5 times
coreset = subset(results, (results$daysd < 0.06 | results$dayrsd < 0.25) & results$unidaypeptides >2 & results$daypeptides > 4 & (results$nightsd < 0.06 | results$nightrsd < 0.25) & results$uninightpeptides > 2 & results$nightpeptides > 4)
coreset = coreset[!duplicated(coreset$AGI),]

write.csv(results, file = "results.csv")
write.csv(paperresults, file = "partialfulllist.csv")
write.csv(coreset, file = "partialcorelist.csv")
saveRDS(daydatamelt, file = "daydatamelt.rds")
saveRDS(nightdatamelt, file = "nightdatamelt.rds")
saveRDS(proteinInfo, file = "proteininfo.rds")
