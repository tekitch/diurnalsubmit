#requires 3 reps, naming format is 'd' = deuterium, '5' = 5% max label, 24 = hours after label addition, 1 = rep, -F01 = fraction ie d5241-F01.mzML
#3 reps because i lazily subsetted the labelled and unlabelled by doing a 1:3 and 4:6
#need to modify to work with greater than 9 % deutrium

medianWithoutNA<-function(x) {
  median(x[which(!is.na(x))])
}
  
library(sqldf)
library(stringr)
library(RPostgreSQL)
library(plyr)
library(xcms)
library(MASS)
library(reshape2)
library(profvis)
library(fitdistrplus)
library(nnls)
library(data.table)
library(base64enc)
library(png)
library(seqinr)

###################################################
databaseName = "diurnalfive"
###################################################
mzT = 100
###################################################
maxLabel = 0.05
###################################################
numberofSecondsAroundMaxMono = 4
###################################################
numberOfDistributionsForNnls = 10
###################################################
numEICs = 12
###################################################
fastaFile = "at.fasta"
###################################################

source('~/ms/scripts/rdeuterium/november2017/peptide-property-functions.R')

fasta = read.fasta(fastaFile, seqtype = "AA", as.string = TRUE, seqonly = TRUE)
fasta = unlist(fasta)

mzMLFileNames = list.files(pattern = "*.mzML")                               #label, reps, fractions
pepXMLFileNames = list.files(pattern = "*-interact.ipro.pep.xml$")           #label, reps
protXMLFileNames = list.files(pattern = "*.prot.xml")                        #label, reps

mzMLNames = sapply(strsplit(mzMLFileNames,"-"), `[`, 1)

mzMLNames = unique(mzMLNames)
pepXMLNames = sapply(strsplit(pepXMLFileNames,"-"), `[`, 1)
databasePepNames = paste(pepXMLNames, "prophetpeptides", sep="")
databaseProtNames = paste(pepXMLNames, "proteins", sep ="")
labelledSampleNames = subset(pepXMLNames, grepl("^.*R.*", pepXMLNames))                                   #"^[a-z][^0].*", pepXMLNames))


for (i in 1:length(mzMLFileNames)){
  string  = "java -jar /home/anon/ms/tools/dinosaur/Dinosaur.jar"
  system(paste(string, mzMLFileNames[i]))
}

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
                 user = "postgres", password = "postgres")
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

##############################
#Match features to peptide IDs
##############################

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

#This line should be adjusted to chop off the front of the gradient before peptide elution starts in the case of equilibration before samples start. In this case there was 1400 seconds before the start of the gradient

for(i in 1:length(fracFeatides)){
  for(j in 1:length(fracFeatides[[i]])){
    fracFeatides[[i]][[j]] = subset(fracFeatides[[i]][[j]], fracFeatides[[i]][[j]]$retentiontime>1400)
  }
}

#split fracFeatides into lablelled and unlabelled samples

unlabFeatides = fracFeatides[c(4, 8)]   #only using one unlabelled sample added in twice to the list to correspond with th etwo labelled samples. 
labFeatides = fracFeatides[c(1,2,3,5,6,7)]
#templabelelmzML = grepl("^.[1-9].*", mzMLNames)
templabelelmzML = grepl("^.*R.*", mzMLNames)
labelledmzMLnames = mzMLNames[templabelelmzML]
tempunlabelelmzML = grepl("^...ID.*", mzMLNames)
unlabelledmzMLnames = mzMLNames[tempunlabelelmzML] #only one unlabelled mzml changed for diurnalyoung test data

#############################################################
#Important lines for configuring unlabelled to labelled samples
##############################################################
corresondingUnlabelled = append(rep(unlabelledmzMLnames[1], 3), rep(unlabelledmzMLnames[2], 3))
unlabelledmzMLnames = corresondingUnlabelled

#############
#RT alignment
#############

#Generate Models
temprtModels = list()
rtModels = list()

#this is a hack to get the numbers of unlabelled samples to match the number of labelled samples so that the lists can be iterated over together. There's probably a better way using the names or something
#reps = length(labFeatides)/length(unlabFeatides)
correspondingunlabFeatides = append(rep(unlabFeatides[1], 3), rep(unlabFeatides[2], 3))
#correspondingunlabFeatides2 = rep(unlabFeatides[2], reps)
#correspondingunlabFeatides = append(correspondingunlabFeatides1, correspondingunlabFeatides2)

for(i in 1:length(labFeatides)){
  for (j in 1:length(labFeatides[[i]])){
    labUniFea = labFeatides[[i]][[j]][!duplicated(labFeatides[[i]][[j]]$peptide), ]
    unlabUniFea = correspondingunlabFeatides[[i]][[j]][!duplicated(correspondingunlabFeatides[[i]][[j]]$peptide), ]
    conjoined = merge(labUniFea, unlabUniFea, by.x = "peptide", by.y = "peptide")
    conjoined = arrange(conjoined, conjoined$retentiontime.y)
    loessModel = loess(retentiontime.x~retentiontime.y, data = conjoined, span = .1, degree = 1, family="symmetric")
    plot(conjoined$retentiontime.x~conjoined$retentiontime.y, main = j, xlim = c(0,2500), ylim = c(0,2500))
    hat1 <- predict(loessModel)
    lines(conjoined$retentiontime.y[order(conjoined$retentiontime.y)], hat1[order(conjoined$retentiontime.y)], col="red", lwd=2)
    temprtModels[j] = list(loessModel)
    print(j)
  }
  rtModels[i] = list(temprtModels)
}


#add corrected RTs for labelled samples to the unlabelled dataframe 
print('add corrected RTs to the unlabelled dataframe')
for (i in 1:length(correspondingunlabFeatides)){
  for (j in 1:length(correspondingunlabFeatides[[i]])){
    corfRTStartVec = c()
    corfRTEndVec = c()
    for (k in 1:nrow(correspondingunlabFeatides[[i]][[j]])){
      fRTStart = correspondingunlabFeatides[[i]][[j]][k,]$frtstart
      fRTEnd = correspondingunlabFeatides[[i]][[j]][k,]$frtend
      corfRTStartVec[k] = predict(rtModels[[i]][[j]], newdata = fRTStart)
      corfRTEndVec[k] = predict(rtModels[[i]][[j]], newdata = fRTEnd)
    }
    print(j)
    correspondingunlabFeatides[[i]][[j]]$corfrtstart = corfRTStartVec
    correspondingunlabFeatides[[i]][[j]]$corfrtend = corfRTEndVec
  }
}
unlabFeatides = correspondingunlabFeatides

########################################
#get isotope profiles from labelled data
########################################

#subset dataframes so only matched peptides are EIC extracted
temp = list()
dataToExtract = list()
for (i in 1:length(unlabFeatides)){
  for (j in 1:length(unlabFeatides[[i]])){
    temp[[j]] = subset(unlabFeatides[[i]][[j]], (!is.na(unlabFeatides[[i]][[j]]$corfrtstart)) & (!is.na(unlabFeatides[[i]][[j]]$corfrtend)))
    labmzrt = data.frame(temp[[j]]$corfrtstart/60, temp[[j]]$corfrtend/60, temp[[j]]$mzLo, temp[[j]]$mzHi)
    colnames(labmzrt) = c("rtLo", "rtHi", "mzLo", "mzHi")
    unlabmzrt = data.frame(temp[[j]]$frtstart/60, temp[[j]]$frtend/60, temp[[j]]$mzLo, temp[[j]]$mzHi)
    colnames(unlabmzrt) = c("rtLo", "rtHi", "mzLo", "mzHi")
    unlabmzrtName = paste0(mzMLNames[i], fracTemp[j], "unlab.mzrt.csv")
    labmzrtName = paste0(mzMLNames[i], fracTemp[j], "lab.mzrt.csv")
    write.csv(unlabmzrt, unlabmzrtName)
    write.csv(labmzrt, labmzrtName)
  }
  dataToExtract[[i]] = temp
}

#subset the data to only look at intense peaks to minimise processing time
rtDeltaDataToExtract = list()
for (i in 1: length(dataToExtract)){
  tempj = list()
  for(j in 1: length(dataToExtract[[i]])){
    tempj[[j]] = subset(dataToExtract[[i]][[j]], dataToExtract[[i]][[j]]$fintensityapex > 5000000)
  }
  rtDeltaDataToExtract[[i]] = tempj
}


#######################################################
# work out where peak should occur in heavier envelopes
#######################################################
graphColours = c("red", "green1", "green2", "green3", "green4", "lightseagreen", "mediumseagreen", "seagreen", "seagreen1",  "seagreen2", "seagreen3", "turquoise", "turquoise1", "turquoise2", "turquoise3", "turquoise4", "steelblue", "steelblue1", "steelblue2", "steelblue3", "steelblue4", "lightslategrey", "lightsteelblue", "gray74", "gray48")
par(mfrow=c(1,2))
rtDeltaSummary = list()
rtDeltaList = list()
for (i in 1: 2){#length(rtDeltaDataToExtract)){
  for(j in 1: 2){#length(rtDeltaDataToExtract[[i]])){
    tempLabmzMLName = paste0(labelledmzMLnames[i], "-", fracTemp[j], ".mzML")
    currentLabmzML = xcmsRaw(tempLabmzMLName)
    secondsAfterMaxMono = list()
    for(k in 1:nrow(rtDeltaDataToExtract[[i]][[j]])){   #could do random sampling here
      print(k)
      iterMzVec = c(rtDeltaDataToExtract[[i]][[j]][k,]$fmz)
      iterPepSeq = rtDeltaDataToExtract[[i]][[j]][k,]$peptide
      iterCharge = rtDeltaDataToExtract[[i]][[j]][k,]$fcharge
      isoSettings = labelledAtomicAbundances(2, "H", maxLabel)
      tf = isotopicDistribution(iterPepSeq, list(H = isoSettings))
      #numEICs = max(which(tf > .05))+1
      for (l in 1:numEICs){
        iterMzVec = append(iterMzVec, iterMzVec[1]+(1/iterCharge*l))        
      }
      mzMatrix = as.matrix(iterMzVec-.01)
      mzMatrix = cbind(mzMatrix, iterMzVec + .01)
      #next line is a hack to prevent the RT correction model from collapsing a feature to very short amount of time and thus extracting only a couple of datapoints
      if(rtDeltaDataToExtract[[i]][[j]][k,]$corfrtend - rtDeltaDataToExtract[[i]][[j]][k,]$corfrtstart < 6 ){
        labrtMatrix = matrix(rep(c(rtDeltaDataToExtract[[i]][[j]][k,]$corfrtstart, rtDeltaDataToExtract[[i]][[j]][k,]$corfrtstart + 15), nrow(mzMatrix)), ncol = 2, byrow = TRUE)
      }else{
        labrtMatrix = matrix(rep(c(rtDeltaDataToExtract[[i]][[j]][k,]$corfrtstart, rtDeltaDataToExtract[[i]][[j]][k,]$corfrtend), nrow(mzMatrix)), ncol = 2, byrow = TRUE)
      }
      unlabrtMatrix = matrix(rep(c(rtDeltaDataToExtract[[i]][[j]][k,]$frtstart, rtDeltaDataToExtract[[i]][[j]][k,]$frtend), nrow(mzMatrix)), ncol = 2, byrow = TRUE)
      currentLabEIC = getEIC(currentLabmzML, mzrange = mzMatrix, rtrange = labrtMatrix)
      tempVec = c()
      loessPlotList = list()
      for (m in 1:length(currentLabEIC@eic$xcmsRaw)){
        #lines(currentLabEIC@eic$xcmsRaw[[m]][,1], currentLabEIC@eic$xcmsRaw[[m]][,2], col = graphColours[m])
        x = currentLabEIC@eic$xcmsRaw[[m]][,1]
        if(sum(currentLabEIC@eic$xcmsRaw[[m]][,2])!=0){y <- currentLabEIC@eic$xcmsRaw[[m]][,2]/max(currentLabEIC@eic$xcmsRaw[[m]][,2])}else{y <- currentLabEIC@eic$xcmsRaw[[m]][,2]}
        lo <- loess(y~x, span = .43) # fails when only two data points in EIC
        xl <- seq(min(x),max(x), (max(x) - min(x))/1000)
        if(is.nan(lo$enp)){}else{loessIntVec = predict(lo,xl)
        loessPlotList[m] = list(loessIntVec)
        intPeak = which.max(loessIntVec)
        EICPeaktime = xl[intPeak]
        tempVec[m] = EICPeaktime
        }
      }
      zerodElutionPeaksPerPeptide = tempVec - tempVec[1]
      idxVec = 1:length(tempVec)

      secondsAfterMaxMono[k] = list(zerodElutionPeaksPerPeptide)
    }#for each peptide
    rtDeltaDF = ldply(secondsAfterMaxMono, rbind) #collapses list of variable length vectors into a single DF
    rtDeltaResult = sapply(rtDeltaDF, function(x){median(x)})
    rtDeltaList[j] = list(rtDeltaResult)
  }#for each fraction
  rtDeltaSummary[[i]] = list(rtDeltaList)
}

loongrtDelta = list()
for(i in 1:length(rtDeltaSummary)){
  for(j in 1:length(rtDeltaSummary[[i]])){
    loongrtDelta[i] = list(as.data.frame(ldply(rtDeltaSummary[[i]][[j]], rbind)))
  }
}
rtDeltaOneDF = ldply(loongrtDelta)
deltaRTs = sapply(rtDeltaOneDF, function(x) medianWithoutNA(x))
vc = c(1,2,3,4,5)
newdata = data.frame(deltaRTs[1:5], vc)
names(newdata) = c("deltart", "vc")
plot(newdata$vc, newdata$deltart, ylim = c(-10,0), xlim = c(1,10))
fit_lm = lm(deltart~I(vc*vc^3), data = newdata)
lines(c(1:10), predict.lm(fit_lm, data.frame(vc=c(1:10))), type = "l")
secondsPer2H = predict.lm(fit_lm, data.frame(vc=c(1:60)))
secondsPer2H[1] = 0
#secondsPer2H = c("0 -0.0418796117250993", "-0.088471215178537", "-0.213910147553177", "-0.478407096388847", "-0.959375802808181", "-1.75143306151662", "-2.96639872080242", "-4.73329568253664", "-7.19834990217314", "-10.5249903887486", "-14.8938492048825", "-20.5027614667771", "-27.5667653442176", "-36.3181020605717", "-47.0062158927903", "-59.8977541714069", "-75.2765672805378", "-93.4437086578821", "-114.717434794722", "-139.433205235922", "-167.943682579929", "-200.618732478775", "-237.845423638071", "-280.028027817015", "-327.588019828385", "-380.964077538543", "-440.612081867432", "-507.005116788581", "-580.633469329099", "-662.004629569679", "-751.643290644597", "-850.091348741711", "-957.907903102462", "-1075.66925602188", "-1203.96891284856", "-1343.4175819847", "-1494.64317488607", "-1658.29080606202", "-1835.0227930755", "-2025.51865654302", "-2230.4751201347", "-2450.6061105742", "-2686.64275763882", "-2939.33339415939", "-3209.44355602035", "-3497.75598215971", "-3805.07061456909", "-4132.20459829367", "-4479.9922814322", "-4849.28521513704", "-5240.95215361412", "-5655.87905412295", "-6094.96907697664", "-6559.14258554186", "-7049.33714623887", "-7566.50752854153", "-8111.62570497725", "-8685.68085112706", "-9289.67934562554")

############################################
#Add the labelled IDs to the unlabelled ones
#making a corrected time to match the unlabs
############################################


#correspdataToExtracta = rep(dataToExtract[1], reps)
#correspdataToExtractb = rep(dataToExtract[1], reps)
#correspdataToExtract = append(correspdataToExtracta, correspdataToExtractb)

for(i in 1:length(labFeatides)){
  for(j in 1:length(labFeatides[[i]])){
    labFeatides[[i]][[j]]$corfrtstart = labFeatides[[i]][[j]]$frtstart
    labFeatides[[i]][[j]]$corfrtend = labFeatides[[i]][[j]]$frtend
    labFeatides[[i]][[j]] = subset(labFeatides[[i]][[j]], labFeatides[[i]][[j]]$corfrtstart > 5)
    dataToExtract[[i]][[j]] = rbind(dataToExtract[[i]][[j]], labFeatides[[i]][[j]])
  }
}
#dataToExtract = correspdataToExtract 

###########################################
#cleanup before mainloop
###########################################
#rm(featuresMList, correspdataToExtract, correspdataToExtracta, correspdataToExtractb, correspondingunlabFeatides, correspondingunlabFeatides1, correspondingunlabFeatides2, fracFeatides, labFeatides, rtModels, temp, unlabFeatides)
gc()
###########################################
#The Main Loop
###########################################
con <- dbConnect(drv, dbname = databaseName, host = "localhost", port = 5432, user = "postgres", password = "postgres")
sqlNames = c()

par(mfrow=c(1,1))
for (i in 1: length(dataToExtract)){
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
    tempLabmzMLName = paste0(labelledmzMLnames[i], "-", fracTemp[j], ".mzML")
    currentLabmzML = xcmsRaw(tempLabmzMLName)
    tempUnlabmzMLName = paste0(unlabelledmzMLnames[i], "-", fracTemp[j], ".mzML")
    currentUnLabmzML = xcmsRaw(tempUnlabmzMLName)
    #for each?
    for(k in 1:nrow(dataToExtract[[i]][[j]])){
      iterMzVec = c(dataToExtract[[i]][[j]][k,]$fmz)
      iterPepSeq = dataToExtract[[i]][[j]][k,]$peptide
      iterCharge = dataToExtract[[i]][[j]][k,]$fcharge
      iterFWHM = dataToExtract[[i]][[j]][k,]$ffwhm
      isoSettings = labelledAtomicAbundances(2, "H", maxLabel)
      tf = isotopicDistribution(iterPepSeq, list(H = isoSettings))
      for (l in 1:numEICs){
        iterMzVec = append(iterMzVec, iterMzVec[1]+(1/iterCharge*l))        
      }
      mzMatrix = as.matrix(iterMzVec-.01)
      mzMatrix = cbind(mzMatrix, iterMzVec + .01)
      if(dataToExtract[[i]][[j]][k,]$corfrtend - dataToExtract[[i]][[j]][k,]$corfrtstart < 6 ){
        labrtMatrix = matrix(rep(c(dataToExtract[[i]][[j]][k,]$corfrtstart, dataToExtract[[i]][[j]][k,]$corfrtstart + 15), nrow(mzMatrix)), ncol = 2, byrow = TRUE)
      }else{
        labrtMatrix = matrix(rep(c(dataToExtract[[i]][[j]][k,]$corfrtstart, dataToExtract[[i]][[j]][k,]$corfrtend), nrow(mzMatrix)), ncol = 2, byrow = TRUE)
      } #these lines extract a minimum of 15 sec if the feature found is less than 6 seconds long
      #had a problem where the corrected feature rt was 0 so with the -2H correction it's extracting EICs from negative rt values making a final EIC which had no values in it - couldn't find the minimum
      secondsPer2HTemp = secondsPer2H[1:nrow(labrtMatrix)]
      labrtMatrix = labrtMatrix + secondsPer2HTemp
      labrtMatrix[,2] = abs(labrtMatrix[,2]) #this is the fix - keep the end value of the EIC positive; it's all nonsense anyway and will be filtered out. 
      unlabrtMatrix = matrix(rep(c(dataToExtract[[i]][[j]][k,]$frtstart, dataToExtract[[i]][[j]][k,]$frtend), nrow(mzMatrix)), ncol = 2, byrow = TRUE)
      currentLabEIC = getEIC(currentLabmzML, mzrange = mzMatrix, rtrange = labrtMatrix)
      currentUnLabEIC = getEIC(currentUnLabmzML, mzrange = mzMatrix, rtrange = unlabrtMatrix)
      labrtidx = which.max(currentLabEIC@eic$xcmsRaw[[1]][,2])
      labrtatMax = currentLabEIC@eic$xcmsRaw[[1]][labrtidx,1]
      labrtatMaxLo = labrtatMax - iterFWHM/2
      labrtatMaxHi = labrtatMax + iterFWHM/2
      
      unlabrtidx = which.max(currentUnLabEIC@eic$xcmsRaw[[1]][,2])
      unlabrtatMax = currentUnLabEIC@eic$xcmsRaw[[1]][unlabrtidx,1]
      unlabrtatMaxLo = unlabrtatMax - iterFWHM/2
      unlabrtatMaxHi = unlabrtatMax + iterFWHM/2
      
      labnnlsVec = c()
      unlabnnlsVec = c()
      for (x in 1:length(currentLabEIC@eic$xcmsRaw)){
        labLoidx = which(abs(currentLabEIC@eic$xcmsRaw[[x]][,1]-labrtatMaxLo)==min(abs(currentLabEIC@eic$xcmsRaw[[x]][,1]-labrtatMaxLo)))
        labHiidx = which(abs(currentLabEIC@eic$xcmsRaw[[x]][,1]-labrtatMaxHi)==min(abs(currentLabEIC@eic$xcmsRaw[[x]][,1]-labrtatMaxHi)))
        labnnlsVec[x] = mean(currentLabEIC@eic$xcmsRaw[[x]][,2][labLoidx:labHiidx])
        
        unlabLoidx = which(abs(currentUnLabEIC@eic$xcmsRaw[[x]][,1]-unlabrtatMaxLo)==min(abs(currentUnLabEIC@eic$xcmsRaw[[x]][,1]-unlabrtatMaxLo)))
        unlabHiidx = which(abs(currentUnLabEIC@eic$xcmsRaw[[x]][,1]-unlabrtatMaxHi)==min(abs(currentUnLabEIC@eic$xcmsRaw[[x]][,1]-unlabrtatMaxHi)))
        unlabnnlsVec[x] = mean(currentUnLabEIC@eic$xcmsRaw[[x]][,2][unlabLoidx:unlabHiidx])
      }
      if(sum(labnnlsVec != 0)){tempnormlabnnlsVec = labnnlsVec / sum(labnnlsVec)}else{tempnormlabnnlsVec = labnnlsVec}
      tempnormunlabnnlsVec = unlabnnlsVec / sum(unlabnnlsVec)
      #plot(tempnormlabnnlsVec, type = "h", col = "red", main = paste(dataToExtract[[i]][[j]][k,]$peptide, "sec around max mono"))
      #lines(tempnormunlabnnlsVec, type = "l")
      numberOfRatios = length(currentLabEIC@eic$xcmsRaw)
      nnlsMatrix = matrix(NA, ncol = numberOfDistributionsForNnls, nrow = numberOfRatios)
      peptideNatural =  isotopicDistribution(iterPepSeq)
      nnlsMatrix[,1] = peptideNatural[1:numberOfRatios]
      nnlsMatrix[is.na(nnlsMatrix)] = 0 #replaces na with 0 - sometimes the peptidenatural script doesn't return enough values
      nnlsXLabs = c()
      nnlsPopFactor = maxLabel / ncol(nnlsMatrix)
      for (u in 1:ncol(nnlsMatrix)){
        h=u-1
        enrichmentThisIteration = h*nnlsPopFactor
        #print(enrichmentThisIteration)
        vtypor = labelledAtomicAbundances(2, "H", enrichmentThisIteration)
        tf = isotopicDistribution(iterPepSeq, list(H = vtypor))
        nnlsMatrix[,u] = tf[1:nrow(nnlsMatrix)]
      }
      nnlsMatrix[is.na(nnlsMatrix)] = 0
      
      nnlsout = nnls(nnlsMatrix, tempnormlabnnlsVec)
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
      unlabToSynthZero[k] = cor(tempnormunlabnnlsVec, nnlsMatrix[,1])
      
      png("temp.png")
      par(fig = c(0,1,0,1))
      plot(tempnormlabnnlsVec, type = 'l', col = 'red', ylim = c(0,.8), main = paste0("DEV ", round(nnlsout$deviance,4), " LPF ", round(LPF[k], 3), " CHG ", iterCharge, "\n", iterPepSeq, " OCC ", occurencestemp))
      lines(tempnormunlabnnlsVec, type = 'l', col = 'blue')
      for(v in 1:ncol(nnlsMatrix)){
        lines(nnlsMatrix[,v], col = "grey")
      }
      
      #plot(nnlsout$x, type = 'h')
      par(fig = c(0.3,1, 0.3, 1), new = T)
      plot(currentLabEIC@eic$xcmsRaw[[1]], type = 'l', col = 'red', xaxt = 'n', ann = FALSE)
      for(e in 2:length(currentLabEIC@eic$xcmsRaw)){
        lines(currentLabEIC@eic$xcmsRaw[[e]], col = graphColours[e])
      }
      abline(v = labrtatMaxLo)
      abline(v = labrtatMaxHi)
      par(fig = c(0.01,0.42, 0.4, .99), new = T)
      plot(nnlsout$x, type = 'h', xaxt = 'n', ann = FALSE, yaxt = 'n')
      dev.off()
      encPlot = base64encode("temp.png")
      base64Plots[k] = encPlot
    }
    dataToExtract[[i]][[j]]$nnlsresult = nnlsResult
    dataToExtract[[i]][[j]]$nnlsresiduals = nnlsResiduals
    dataToExtract[[i]][[j]]$nnlsdeviance = nnlsDeviance
    dataToExtract[[i]][[j]]$lpf = LPF
    dataToExtract[[i]][[j]]$lpfsum = LPFsum
    dataToExtract[[i]][[j]]$lastpop = lastpop
    dataToExtract[[i]][[j]]$occurences = occurences 
    dataToExtract[[i]][[j]]$encplot = base64Plots
    dataToExtract[[i]][[j]]$unlabtosynthzero = unlabToSynthZero
  }
  tempdata = rbindlist(dataToExtract[[i]])
  sqlName = paste0(labelledSampleNames[i], "TO") #should change this to get the names from the labFeatides df - this breaks if there are extra samples in the analysis ie quant samples
  sqlName = tolower(sqlName)
  dbWriteTable(con, sqlName, tempdata)
  dbGetQuery(con, paste0("ALTER TABLE ", sqlName, " ADD COLUMN pk SERIAL PRIMARY KEY;"))
  sqlNames[i] = sqlName
  dataToExtract[[i]] = "Data in postgres, removed from here to save memory which was limiting total dataset size possible"
  gc()
}
print("work complete")