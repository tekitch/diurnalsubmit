setwd("d:/diurnalfive")

partialfulllist = read.csv("results.csv")
quantresults = read.csv("quantresults.csv")
datamerge = merge(partialfulllist, quantresults, by.x = 'protein', by.y = 'protein', suffixes = c('pr','qt'))

fcpdayvec = c()
fcpdayvec1 = c()
fcpnightvec = c()
rawFCPdayvec = c()
rawFCPnightvec = c()
for(i in 1:nrow(datamerge)){
  fcpdayvec[i] = ((1-datamerge[i,'daylpfqt'])/datamerge[i,'daylpfqt'])/((1-datamerge[i,'nightlpfqt'])/datamerge[i,'nightlpfqt'])
  fcpnightvec[i] = ((1-datamerge[i,'nightlpfqt'])/datamerge[i,'nightlpfqt'])/((1-datamerge[i,'daylpfqt'])/datamerge[i,'daylpfqt'])
  rawFCPdayvec[i] = ((1-datamerge[i,'daylpfqt'])/datamerge[i,'daylpfqt'])/((1-datamerge[i,'nightlpfqt'])/datamerge[i,'nightlpfqt'])
  rawFCPnightvec[i] = ((1-datamerge[i,'nightlpfqt'])/datamerge[i,'nightlpfqt'])/((1-datamerge[i,'daylpfqt'])/datamerge[i,'daylpfqt'])
  try(if(datamerge[i, 'p.value'] > 0.05){
   fcpdayvec[i] = 1
   fcpnightvec[i] = 1
  })
}

t = 1
Kddayvec = c()
Kdnightvec = c()
Ksdayvec = c()
Ksnightvec = c()
for(i in 1:nrow(datamerge)){
  Kddayvec[i] = tryCatch( -((log(fcpdayvec[i]*(1-datamerge[i,'daylpfpr']))/t)), finally = NA)
  Kdnightvec[i] = tryCatch( -((log(fcpnightvec[i]*(1-datamerge[i,'nightlpfpr']))/t)), finally = NA)
  Ksdayvec[i] = tryCatch( ((fcpdayvec[i]-exp(-Kddayvec[i]*t))/(1-exp(-Kddayvec[i]*t)))*Kddayvec[i], finally = NA)
  Ksnightvec[i] = tryCatch( ((fcpnightvec[i]-exp(-Kdnightvec[i]*t))/(1-exp(-Kdnightvec[i]*t)))*Kdnightvec[i], finally = NA)
}

datamerge$rawdayFCP = rawFCPdayvec
datamerge$rawnightFCP = rawFCPnightvec
datamerge$dayFCP = fcpdayvec
datamerge$nightFCP = fcpnightvec
datamerge$dayKd = Kddayvec
datamerge$nightKd = Kdnightvec
datamerge$dayKs = Ksdayvec
datamerge$nightKs = Ksnightvec

coreset = datamerge[complete.cases(datamerge),]
coreset = subset(datamerge, (datamerge$daysdpr < 0.06 | datamerge$dayrsdpr < 0.25) & datamerge$unidaypeptidespr >2 & datamerge$daypeptidespr > 4 & (datamerge$nightsdpr < 0.06 | datamerge$nightrsdpr < 0.25) & datamerge$uninightpeptidespr > 2 & datamerge$nightpeptidespr > 4)
coreset = coreset[!duplicated(coreset$AGI),]
coreset = coreset[!is.na(coreset[,'p.value']),]

write.csv(datamerge, file = "unflitKdKs.csv")
write.csv(coreset, file = "coreset.csv")
