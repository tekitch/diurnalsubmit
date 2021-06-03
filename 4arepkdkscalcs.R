setwd("d:/diurnalfive")
represults = readRDS("represults.rds")
quantresults = read.csv("quantresults.csv")
datamerge = list()
for(i in 1:length(represults)){
  datamerge[[i]] = merge(represults[[i]], quantresults, by.x = 'protein', by.y = 'protein', suffixes = c('pr','qt'))
}

for(j in 1:length(datamerge)){
  fcpdayvec = c()
  fcpnightvec = c()
  rawFCPdayvec = c()
  rawFCPnightvec = c()
  for(i in 1:nrow(datamerge[[j]])){
     fcpdayvec[i] = ((1-datamerge[[j]][i,'daylpfqt'])/datamerge[[j]][i,'daylpfqt'])/((1-datamerge[[j]][i,'nightlpfqt'])/datamerge[[j]][i,'daylpfqt'])
     fcpnightvec[i] = ((1-datamerge[[j]][i,'nightlpfqt'])/datamerge[[j]][i,'nightlpfqt'])/((1-datamerge[[j]][i,'daylpfqt'])/datamerge[[j]][i,'daylpfqt'])
     rawFCPdayvec[i] = ((1-datamerge[[j]][i,'daylpfqt'])/datamerge[[j]][i,'daylpfqt'])/((1-datamerge[[j]][i,'nightlpfqt'])/datamerge[[j]][i,'daylpfqt'])
     rawFCPnightvec[i] = ((1-datamerge[[j]][i,'nightlpfqt'])/datamerge[[j]][i,'nightlpfqt'])/((1-datamerge[[j]][i,'daylpfqt'])/datamerge[[j]][i,'daylpfqt'])
     try(if(datamerge[[j]][i, 'p.value'] > 0.05){
       fcpdayvec[i] = 1
       fcpnightvec[i] = 1
     })
  }
  t = 1
  Kddayvec = c()
  Kdnightvec = c()
  Ksdayvec = c()
  Ksnightvec = c()
  for(i in 1:nrow(datamerge[[j]])){
    Kddayvec[i] = tryCatch( -((log(fcpdayvec[i]*(1-datamerge[[j]][i,'daylpfpr']))/t)), finally = NA)
    Kdnightvec[i] = tryCatch( -((log(fcpnightvec[i]*(1-datamerge[[j]][i,'nightlpfpr']))/t)), finally = NA)
    Ksdayvec[i] = tryCatch( ((fcpdayvec[i]-exp(-Kddayvec[i]*t))/(1-exp(-Kddayvec[i]*t)))*Kddayvec[i], finally = NA)
    Ksnightvec[i] = tryCatch( ((fcpnightvec[i]-exp(-Kdnightvec[i]*t))/(1-exp(-Kdnightvec[i]*t)))*Kdnightvec[i], finally = NA)
  }
  
  datamerge[[j]]$rawdayFCP = rawFCPdayvec
  datamerge[[j]]$rawnightFCP = rawFCPnightvec
  datamerge[[j]]$dayFCP = fcpdayvec
  datamerge[[j]]$nightFCP = fcpnightvec
  datamerge[[j]]$dayKd = Kddayvec
  datamerge[[j]]$nightKd = Kdnightvec
  datamerge[[j]]$dayKs = Ksdayvec
  datamerge[[j]]$nightKs = Ksnightvec
}
coreset = list()

for(i in 1:length(datamerge)){
  coreset[[i]] = datamerge[[i]][complete.cases(datamerge[[i]]),]
  coreset[[i]] = coreset[[i]][!duplicated(coreset[[i]]$AGI),]
  coreset[[i]] = coreset[[i]][!is.na(coreset[[i]][,'p.value']),]
}
corereps = coreset

saveRDS(corereps, file = "kdksreps.rds")