library(RPostgreSQL)
library(plyr)
library(xcms)
library(nnls)
library(seqinr)
library(stringr)
library(data.table)

databaseName = "spikechk"
maxLabel = .98
numberOfDistributionsForNnls = 20
fastaFile = "at.fasta"

fasta = read.fasta(fastaFile, seqtype = "AA", as.string = TRUE, seqonly = TRUE)
fasta = unlist(fasta)

source('C:/Scripts/peptide-property-functions.R')

mzMLFileNames = list.files(pattern = "*.mzML")                               #label, reps, fractions
pepXMLFileNames = list.files(pattern = "*-interact.ipro.pep.xml$")           #label, reps
protXMLFileNames = list.files(pattern = "*.prot.xml")                        #label, reps
mzMLNames = sapply(strsplit(mzMLFileNames,"-"), `[`, 1)
mzMLNames = unique(mzMLNames)
pepXMLNames = sapply(strsplit(pepXMLFileNames,"-"), `[`, 1)
databasePepNames = paste(pepXMLNames, "prophetpeptides", sep="")
databaseProtNames = paste(pepXMLNames, "proteins", sep ="")

#for (i in 1:length(mzMLFileNames)){
#  string  = "java -jar c:/Scripts/Dinosaur.jar"
#  system(paste(string, mzMLFileNames[i]))
#}

featureFileNames = list.files(pattern = "*.features.tsv")                    #label, reps, fractions
featureNames = sapply(strsplit(featureFileNames,"\\."), `[`, 1)
fracTemp = sapply(strsplit(featureNames,"-"), `[`, 2)
fracTemp = unique(fracTemp)

featuresMList = list()
tempList = list()
tempString = ".features.tsv"
for (i in 1:length(mzMLNames)){
  for (j in 1:length(fracTemp)){
    tempFeatureTable = read.delim(paste(mzMLNames[i], "-", fracTemp[j], tempString, sep = ""), sep = "\t", header = TRUE)
    tempFeatureTable$rtStart = tempFeatureTable$rtStart * 60
    tempFeatureTable$rtApex = tempFeatureTable$rtApex * 60
    tempFeatureTable$rtEnd = tempFeatureTable$rtEnd * 60
    tempFeatureTable$fwhm = tempFeatureTable$fwhm * 60
    tempFeatureTable$length = tempFeatureTable$rtEnd - tempFeatureTable$rtStart
    tempFeatureTable = subset(tempFeatureTable, tempFeatureTable$intensityApex > 30000 & tempFeatureTable$length > 2 & tempFeatureTable$nIsotopes > 2)      #filtering out features less than 6 sec long and under 5000 intensity, removes 2/3 of features
    tempFeatureTable$id = row.names(tempFeatureTable)
    tempList[[j]] = tempFeatureTable
  }
  names(tempList) = fracTemp
  featuresMList[i] = list(tempList)
}

rm(tempList)
names(featuresMList) = mzMLNames

detectUnlabelledFeatureFileNames = grepl("^...ID.*", featureFileNames)
unlabelledFeatureFileNames = featureFileNames[detectUnlabelledFeatureFileNames]

drv <- dbDriver("PostgreSQL")
con <- dbConnect(drv, dbname = databaseName,
                 host = "localhost", port = 5432,
                 user = "", password = "")
sql1 = "select * from "
sql2 = " where spectrumid like "

pepXMLList = list()
pepXMLMList = list()
tempNames = NULL
for (i in 1:length(databasePepNames)){
  for (z in 1:length(fracTemp)){
    tempPepSql = dbGetQuery(con, paste(sql1, databasePepNames[i], sql2, "'%-", fracTemp[z], "%'", sep = ""))
    tempPepSql$obsmz = (tempPepSql$precneutmass + tempPepSql$charge)/tempPepSql$charge
    pepXMLList[[z]] = tempPepSql
    print(paste(sql1, databasePepNames[i], sql2, "'%-", fracTemp[z], "%'", sep = ""))
  }
  for(h in 1:length(fracTemp)){tempNames[h] = paste(pepXMLNames[i], fracTemp[h], sep = "-")}
  names(pepXMLList) = fracTemp
  pepXMLMList[i] = list(pepXMLList)
}
rm(pepXMLList)
names(pepXMLMList) = pepXMLNames
dbDisconnect(con)
mzT = 100

fracFeatides = list()
tempList = list()
for (i in 1:length(pepXMLMList)){
  for (j in 1:length(pepXMLMList[[i]])){
    tempFeatureDF = featuresMList[[i]][[j]]
    tempPepDF = pepXMLMList[[i]][[j]]
    mzVec = vector(mode = "numeric", length = nrow(tempPepDF))
    mostabundantmzVec = vector(mode ="numeric", length = nrow(tempPepDF))
    chargeVec = vector(mode ="numeric", length = nrow(tempPepDF))
    rtStartVec = vector(mode ="numeric", length = nrow(tempPepDF))
    rtApexVec = vector(mode ="numeric", length = nrow(tempPepDF))
    rtEndVec = vector(mode ="numeric", length = nrow(tempPepDF))
    fwhmVec = vector(mode ="numeric", length = nrow(tempPepDF))
    nIsotopesVec = vector(mode ="numeric", length = nrow(tempPepDF))
    nScansVec = vector(mode ="numeric", length = nrow(tempPepDF))
    averagineCorrVec = vector(mode ="numeric", length = nrow(tempPepDF))
    massVec = vector(mode ="numeric", length = nrow(tempPepDF))
    massCalibVec = vector(mode ="numeric", length = nrow(tempPepDF))
    intensityApexVec = vector(mode ="numeric", length = nrow(tempPepDF))
    intensitySumVec = vector(mode ="numeric", length = nrow(tempPepDF))
    lengthVec = vector(mode ="numeric", length = nrow(tempPepDF))
    featureidVec = vector(mode ="numeric", length = nrow(tempPepDF))
    rtLoVec = vector(mode ="numeric", length = nrow(tempPepDF))
    rtHiVec = vector(mode ="numeric", length = nrow(tempPepDF))
    mzLoVec = vector(mode ="numeric", length = nrow(tempPepDF))
    mzHiVec = vector(mode ="numeric", length = nrow(tempPepDF))
    IDrtLoVec = vector(mode ="numeric", length = nrow(tempPepDF))
    IDrtHiVec = vector(mode ="numeric", length = nrow(tempPepDF))
    IDmzLoVec = vector(mode ="numeric", length = nrow(tempPepDF))
    IDmzHiVec = vector(mode ="numeric", length = nrow(tempPepDF))
    for (k in 1:nrow(tempPepDF)){
      rt = tempPepDF[k,]$retentiontime
      mz = tempPepDF[k,]$obsmz
      mzTol = mz * mzT * 1e-6
      tempFeatureDF$rtStart = tempFeatureDF$rtStart - 5  #allowing the ID to be just before the start of the feature - ID retention time tends to be very early in the feature
      tempMatchPepFeature = tempFeatureDF[tempFeatureDF$mz > mz - mzTol & tempFeatureDF$mz < mz + mzTol & rt > tempFeatureDF$rtStart & rt < tempFeatureDF$rtEnd,]
      tempFeatureDF$rtStart = tempFeatureDF$rtStart + 5
      IDmzLoVec[k] = mz
      IDmzHiVec[k] = mz + 0.2
      IDrtLoVec[k] = rt
      IDrtHiVec[k] = rt + 0.1
      if(nrow(tempMatchPepFeature)>1){
        sorttemp = arrange(tempMatchPepFeature, tempMatchPepFeature$nIsotopes, decreasing = TRUE)
        tempMatchPepFeature = sorttemp[1,]}
      if(nrow(tempMatchPepFeature)>0){
        mzVec[k] = tempMatchPepFeature$mz
        mostabundantmzVec[k] = tempMatchPepFeature$mostAbundantMz
        chargeVec[k] = tempMatchPepFeature$charge
        rtStartVec[k] = tempMatchPepFeature$rtStart
        rtApexVec[k] = tempMatchPepFeature$rtApex
        rtEndVec[k] = tempMatchPepFeature$rtEnd
        fwhmVec[k] = tempMatchPepFeature$fwhm
        nIsotopesVec[k] = tempMatchPepFeature$nIsotopes
        nScansVec[k] = tempMatchPepFeature$nScans
        averagineCorrVec[k] = tempMatchPepFeature$averagineCorr
        massVec[k] = tempMatchPepFeature$mass
        massCalibVec[k] = tempMatchPepFeature$massCalib
        intensityApexVec[k] = tempMatchPepFeature$intensityApex
        intensitySumVec[k] = tempMatchPepFeature$intensitySum
        lengthVec[k] = tempMatchPepFeature$length
        featureidVec[k] = tempMatchPepFeature$id
        rtLoVec[k] = tempMatchPepFeature$rtStart
        rtHiVec[k] = tempMatchPepFeature$rtEnd
        mzLoVec[k] = tempMatchPepFeature$mz
        mzHiVec[k] = tempMatchPepFeature$mz + (1/tempMatchPepFeature$charge) * tempMatchPepFeature$nIsotopes
      }
    }
    tempPepDF$fmz = mzVec
    tempPepDF$fmstabundnt = mostabundantmzVec
    tempPepDF$fcharge = chargeVec
    tempPepDF$frtstart = rtStartVec
    tempPepDF$frtapex = rtApexVec
    tempPepDF$frtend = rtEndVec
    tempPepDF$ffwhm = fwhmVec
    tempPepDF$fnisotopes = nIsotopesVec
    tempPepDF$fscans = nScansVec
    tempPepDF$faveraginecorr = averagineCorrVec
    tempPepDF$fmassvec = massVec
    tempPepDF$fmasscalibvec = massCalibVec
    tempPepDF$fintensityapex = intensityApexVec
    tempPepDF$fintensitysum = intensitySumVec
    tempPepDF$flength = lengthVec
    tempPepDF$rtLo = rtLoVec
    tempPepDF$rtHi = rtHiVec
    tempPepDF$mzLo = mzLoVec
    tempPepDF$mzHi = mzHiVec
    tempList[[j]] = tempPepDF
    tempMzRtCsv = tempPepDF[,53:56]
    tempMzRtCsv = unique(tempMzRtCsv)
    tempMzRtCsv$rtLo = tempMzRtCsv$rtLo/60
    tempMzRtCsv$rtHi = tempMzRtCsv$rtHi/60
    mzrtName = paste0(mzMLNames[i], fracTemp[j], ".mzrt.csv")
    write.csv(tempMzRtCsv, mzrtName)
    tempIDMzRtCsv = tempPepDF[,53:56]
    tempIDMzRtCsv$mzLo = IDmzLoVec
    tempIDMzRtCsv$mzHi = IDmzHiVec
    tempIDMzRtCsv$rtLo = IDrtLoVec/60
    tempIDMzRtCsv$rtHi = IDrtHiVec/60
    mzrtIDName = paste0(mzMLNames[i], fracTemp[j], "ID.mzrt.csv")
    tempIDMzRtCsv = unique(tempIDMzRtCsv)
    write.csv(tempIDMzRtCsv, mzrtIDName)
  }
  names(tempList) = fracTemp
  fracFeatides[[i]] = tempList
}
names(fracFeatides) = pepXMLNames

temp = list()
dataToExtract = list()
for (i in 1:length(fracFeatides)){
  for (j in 1:length(fracFeatides[[i]])){
    temp[[j]] = subset(fracFeatides[[i]][[j]], (fracFeatides[[i]][[j]]$rtLo != 0))
    mzrt = data.frame(temp[[j]]$frtstart/60, temp[[j]]$frtend/60, temp[[j]]$mzLo, temp[[j]]$mzHi)
    colnames(mzrt) = c("rtLo", "rtHi", "mzLo", "mzHi")
    mzrtName = paste0(mzMLNames[i], fracTemp[j], ".mzrt.csv")
    #write.csv(mzrt, mzrtName)
  }
  dataToExtract[[i]] = temp
}


con <- dbConnect(drv, dbname = databaseName, host = "localhost", port = 5432, user = "postgres", password = "postgres")
sqlNames = c()

par(mfrow=c(1,1))
for (i in 2: length(dataToExtract)){
  print(paste0("working on sample ", i, " of ", length(dataToExtract)))
  for(j in 1: length(dataToExtract[[i]])){
    print(paste0("working on fraction ", j, " of ", length(dataToExtract[[i]])))
    nnlsResult = c()
    nnlsResiduals = c()
    nnlsDeviance = c()
    LPF = c()
    LPFsum = c()
    lastpop = c()
    occurences = c()
    unlabToSynthZero = c()
    base64Plots = c()
    tempmzMLName = paste0(mzMLNames[i], "-", fracTemp[j], ".mzML")
    currentmzML = xcmsRaw(tempmzMLName)

    for(k in 1:nrow(dataToExtract[[i]][[j]])){
      iterMzVec = c(dataToExtract[[i]][[j]][k,]$fmz)
      iterPepSeq = dataToExtract[[i]][[j]][k,]$peptide
      iterCharge = dataToExtract[[i]][[j]][k,]$fcharge
      iterFWHM = dataToExtract[[i]][[j]][k,]$ffwhm
      isoSettings = labelledAtomicAbundances(15, "N", maxLabel)
      tf = isotopicDistribution(iterPepSeq, list(N = isoSettings))
      numEICs = tail(which(0.05 < tf), n=1)
      for (l in 1:numEICs){
        if(iterMzVec[1]+(1/iterCharge*l)<1500){iterMzVec = append(iterMzVec, iterMzVec[1]+(1/iterCharge*l))} #added an if statement so max EIC <1500 mz, else errors thrown
      }
      mzMatrix = as.matrix(iterMzVec-.01) #this needs to be checked - it's quite a wide window
      mzMatrix = cbind(mzMatrix, iterMzVec + .01)
      rtstarttemp = dataToExtract[[i]][[j]][k,]$frtstart
      if(rtstarttemp < 0){rtstarttemp = 1} #stops the EIC being extracted before the start of the run - results in error
      rtMatrix = matrix(rep(c(rtstarttemp, dataToExtract[[i]][[j]][k,]$frtend), nrow(mzMatrix)), ncol = 2, byrow = TRUE)
      currentEIC = getEIC(currentmzML, mzrange = mzMatrix, rtrange = rtMatrix)
      
      
      rtidx = which.max(currentEIC@eic$xcmsRaw[[1]][,2])#this works out where the mono is at max and takes the half max window around it
      rtatMax = currentEIC@eic$xcmsRaw[[1]][rtidx,1]
      rtatMaxLo = rtatMax - iterFWHM/2
      rtatMaxHi = rtatMax + iterFWHM/2

      nnlsVec = c()

      for (x in 1:length(currentEIC@eic$xcmsRaw)){
        Loidx = which(abs(currentEIC@eic$xcmsRaw[[x]][,1]-rtatMaxLo)==min(abs(currentEIC@eic$xcmsRaw[[x]][,1]-rtatMaxLo)))
        Hiidx = which(abs(currentEIC@eic$xcmsRaw[[x]][,1]-rtatMaxHi)==min(abs(currentEIC@eic$xcmsRaw[[x]][,1]-rtatMaxHi)))
        nnlsVec[x] = mean(currentEIC@eic$xcmsRaw[[x]][,2][Loidx:Hiidx])
      }
      
      if(sum(nnlsVec != 0)){tempnormnnlsVec = nnlsVec / sum(nnlsVec)}else{tempnormnnlsVec = nnlsVec}
      tempnormunnnlsVec = nnlsVec / sum(nnlsVec)
      numberOfRatios = length(currentEIC@eic$xcmsRaw)
      nnlsMatrix = matrix(NA, ncol = numberOfDistributionsForNnls +1, nrow = numberOfRatios)
      peptideNatural =  isotopicDistribution(iterPepSeq)
      nnlsMatrix[,1] = peptideNatural[1:numberOfRatios]
      nnlsMatrix[is.na(nnlsMatrix)] = 0 #replaces na with 0 - sometimes the peptidenatural script doesn't return enough values
      nnlsXLabs = c()
      nnlsPopFactor = maxLabel / (ncol(nnlsMatrix)-1)
      for (u in 1:ncol(nnlsMatrix)){
        h=u-1
        enrichmentThisIteration = h*nnlsPopFactor
        vtypor = labelledAtomicAbundances(15, "N", enrichmentThisIteration)
        tf = isotopicDistribution(iterPepSeq, list(N = vtypor))
        nnlsMatrix[,u] = tf[1:nrow(nnlsMatrix)]
      }
      nnlsMatrix[is.na(nnlsMatrix)] = 0
      nnlsout = nnls(nnlsMatrix, tempnormnnlsVec)
      nnlsX = nnlsout$x
      nnlsRes = nnlsout$residuals
      nnlsDev = nnlsout$deviance
      nnlsResult[k] = paste0(nnlsX, sep = "", collapse = ";")
      nnlsResiduals[k] = paste0(nnlsRes, sep = "", collapse = ";")
      nnlsDeviance[k] = nnlsDev
      LPF[k] = sum(nnlsX[2:length(nnlsX)])
      LPFsum[k] = sum(nnlsout$x)
      lastpop[k] = nnlsX[length(nnlsX)]
      occurencestemp = sum(str_count(fasta, iterPepSeq))
      occurences[k] = occurencestemp

    }
    dataToExtract[[i]][[j]]$nnlsResult = nnlsResult
    dataToExtract[[i]][[j]]$nnlsresiduals = nnlsResiduals
    dataToExtract[[i]][[j]]$nnlsDeviance = nnlsDeviance
    dataToExtract[[i]][[j]]$LPF = LPF
    dataToExtract[[i]][[j]]$LPFsum = LPFsum
    dataToExtract[[i]][[j]]$lastpop = lastpop
    dataToExtract[[i]][[j]]$occurences = occurences 

  }
  tempdata = rbindlist(dataToExtract[[i]])
  sqlName = paste0(mzMLNames[i], "TO") #should change this to get the names from the labFeatides df - this breaks if there are extra samples in the analysis ie quant samples
  sqlName = tolower(sqlName)
  dbWriteTable(con, sqlName, tempdata)
  dbGetQuery(con, paste0("ALTER TABLE ", sqlName, " ADD COLUMN pk SERIAL PRIMARY KEY;"))
  sqlNames[i] = sqlName
  dataToExtract[[i]] = "Data in postgres, removed from here to save memory which was limiting total dataset size possible"
  gc()
}
dbDisconnect(con)
print("work complete")
