# needs to extract the eic and profile, write an mzrt file, do it per protein. 
library(RPostgreSQL)
library(dplyr)
library(tidyr)
library(stringr)
library(seqinr)
se <- function(x) sd(x)/sqrt(length(x))

setwd("d:/diurnalfive")

databaseName = "tst"
daydatanames = c("eod0r1to", "eod0r2to", "eod0r3to")
nightdatanames = c("eon0r1to", "eon0r2to", "eon0r3to")
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

dayquantdata = list()
nightquantdata = list()
quantproteinInfo = list()

for(i in 1:length(daydatanames)){
  dayquantdata[[i]] = dbGetQuery(con, paste(sql1, daydatanames[i], sql2, sep = ""))
  nightquantdata[[i]] = dbGetQuery(con, paste(sql1, nightdatanames[i], sql2, sep = ""))
}

names(dayquantdata) = daydatanames
names(nightquantdata) = nightdatanames
quantproteinInfo = dbGetQuery(con, paste(sql1, proteinInfonames, sep = ""))
dbDisconnect(con)


#########
#GET PROTEIN PROPHET PEPTIDE ASSIGNMENTS
#########

protGroupVec = c()
protGroupCount = c()

for(i in 1:length(dayquantdata)){
  protGroupVec = c()
  protGroupCount = c()
  for(j in 1:nrow(dayquantdata[[i]])){
    pepseq = dayquantdata[[i]][j,]$peptide
    proteinIndexVec = grep(pepseq, quantproteinInfo$unipepseq)
    groupFks = c()
    groupProbs = c()
    for(k in 1:length(proteinIndexVec)){
      groupFks[k] = quantproteinInfo[proteinIndexVec[k],]$group_fk
    }
    groupFks = unique(groupFks)
    protGroupVec[j] = paste0(groupFks, sep = "", collapse = ";")
    protGroupCount[j] = length(groupFks)
  }
  dayquantdata[[i]]$protGroups = protGroupVec
  dayquantdata[[i]]$protGroupCount = protGroupCount
}

for(i in 1:length(nightquantdata)){
  protGroupVec = c()
  protGroupCount = c()
  for(j in 1:nrow(nightquantdata[[i]])){
    pepseq = nightquantdata[[i]][j,]$peptide
    proteinIndexVec = grep(pepseq, quantproteinInfo$unipepseq)
    groupFks = c()
    groupProbs = c()
    for(k in 1:length(proteinIndexVec)){
      groupFks[k] = quantproteinInfo[proteinIndexVec[k],]$group_fk
    }
    groupFks = unique(groupFks)
    protGroupVec[j] = paste0(groupFks, sep = "", collapse = ";")
    protGroupCount[j] = length(groupFks)
  }
  nightquantdata[[i]]$protGroups = protGroupVec
  nightquantdata[[i]]$protGroupCount = protGroupCount
}

proteinVec = unique(quantproteinInfo$group_fk)

#########
#FILTERS#
#########
#daydata = subset(daydata, daydata$protGroupCount == 1)
#nightdata = subset(nightdata, nightdata$protGroupCount == 1) #need to make this more sophisticated by sorting peptides to their most likely group by probability
dayquantnorm = c()
nightquantnorm = c()
for(i in 1:length(dayquantdata)){
  #daydata[[i]] = subset(daydata[[i]], daydata[[i]]$unlabToSynthZero > 0.97)
  #nightdata[[i]] = subset(nightdata[[i]], nightdata[[i]]$unlabToSynthZero > 0.97)
  dayquantdata[[i]] = subset(dayquantdata[[i]], dayquantdata[[i]]$LPFsum < 1.2) #add up all of the nnls pops
  nightquantdata[[i]] = subset(nightquantdata[[i]], nightquantdata[[i]]$LPFsum < 1.2)
  dayquantdata[[i]] = subset(dayquantdata[[i]], dayquantdata[[i]]$fintensityapex > 50000)
  nightquantdata[[i]] = subset(nightquantdata[[i]], nightquantdata[[i]]$fintensityapex > 50000)
  dayquantdata[[i]] = subset(dayquantdata[[i]], dayquantdata[[i]]$nnlsDeviance < .0009)
  nightquantdata[[i]] = subset(nightquantdata[[i]], nightquantdata[[i]]$nnlsDeviance < .0009)
  dayquantdata[[i]] = subset(dayquantdata[[i]], dayquantdata[[i]]$isomassd == 0)
  nightquantdata[[i]] = subset(nightquantdata[[i]], nightquantdata[[i]]$isomassd == 0)
  #daydata[[i]] = subset(daydata[[i]], lastpop < .03)
  #nightdata[[i]] = subset(nightdata[[i]], lastpop < .03)
  dayquantdata[[i]] <- dayquantdata[[i]][order(dayquantdata[[i]]$peptide, -abs(dayquantdata[[i]]$fintensitysum) ), ] #sort by id and reverse of abs(value)
  dayquantdata[[i]] <- dayquantdata[[i]][ !duplicated(dayquantdata[[i]]$peptide), ]# take the first row within each id
  nightquantdata[[i]] <- nightquantdata[[i]][order(nightquantdata[[i]]$peptide, -abs(nightquantdata[[i]]$fintensitysum) ), ]
  nightquantdata[[i]] <- nightquantdata[[i]][ !duplicated(nightquantdata[[i]]$peptide), ]
  dayquantdata[[i]] = dayquantdata[[i]] %>% group_by(protGroups) %>% top_n(n = 3, wt = fintensitysum)
  nightquantdata[[i]] = nightquantdata[[i]] %>% group_by(protGroups) %>% top_n(n = 3, wt = fintensitysum)
  dayquantnorm[i] = .5/((sum(dayquantdata[[i]][,'LPF'], na.rm = TRUE)/length(which(!is.na(dayquantdata[[i]][,'LPF'])))))
  nightquantnorm[i] = .5/((sum(nightquantdata[[i]][,'LPF'], na.rm = TRUE)/length(which(!is.na(nightquantdata[[i]][,'LPF'])))))
  dayquantdata[[i]][,'nLPF'] = dayquantdata[[i]][,'LPF'] * dayquantnorm[i]
  nightquantdata[[i]][,'nLPF'] = nightquantdata[[i]][,'LPF'] * nightquantnorm[i]
}



daydatamelt = bind_rows(dayquantdata, .id = "column_label")
nightdatamelt = bind_rows(nightquantdata, .id = "column_label")
results = as.data.frame(matrix(ncol = 15, nrow = length(proteinVec)))
colnames(results) = c("protein", "seqsource", "seqsource2", "daylpf", "daysd", "dayrsd", "dayse", "nightlpf", "nightsd", "nightrsd", "nightse", "daypeptides", "nightpeptides", "unidaypeptides", "uninightpeptides")

for(j in 1:length(proteinVec)){
  protein = as.numeric(proteinVec[j])
  daysubset = subset(daydatamelt, daydatamelt$protGroups == protein)
  nightsubset = subset(nightdatamelt, nightdatamelt$protGroups == protein)
  seqsource = daysubset[1,]$seqsource
  seqsource2 = nightsubset[1,]$seqsource
  dayMeanLPF = mean(daysubset$nLPF)
  daySD = sd(daysubset$nLPF)
  dayRSD = daySD/dayMeanLPF
  dayse = se(daysubset$nLPF)
  nightMeanLPF = mean(nightsubset$nLPF)
  nightSD = sd(nightsubset$nLPF)
  nightRSD = nightSD/nightMeanLPF
  nightse = se(nightsubset$nLPF)
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

protgroups = results$protein
testres = matrix(data = NA, ncol = 5, nrow = length(protgroups))
colnames(testres) = c("protein", "p.value", "stderr", "daymean", "nightmean")
for(i in 1:length(protgroups)){
  tres = NULL
  daygroup = subset(daydatamelt, daydatamelt$protGroups == protgroups[i])
  nightgroup = subset(nightdatamelt, nightdatamelt$protGroups == protgroups[i])
  tres = tryCatch(t.test(daygroup$nLPF, nightgroup$nLPF, na.rm=TRUE), error=function(cond){tres = NA})
  testres[i,1] = protgroups[i]
 if(!is.na(tres[1])){
   testres[i,2] = tres$p.value
   testres[i,3] = tres$stderr
   testres[i,4] = tres$estimate[1]
   testres[i,5] = tres$estimate[2]
   }
}
hist(testres[,'p.value'], breaks = 100)

quantresults = merge(results, testres, by.x = "protein", by.y = "protein")

saveRDS(daydatamelt, file = 'dayquantpeptides.rds')
saveRDS(nightdatamelt, file = 'nightquantpeptides.rds')
write.csv(quantresults, file = 'quantresults.csv')
