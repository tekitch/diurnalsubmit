setwd("d:/diurnalfive")
library(stringr)
library(bayestestR)

unfiltkds = read.csv("unflitKdKs.csv")
paxdb = read.table("paxdb.txt", header = TRUE)
paxdb$internal_id = NULL
paxdb$string_external_id = as.character(paxdb$string_external_id)

for(i in 1:nrow(paxdb)){
  paxdb$string_external_id[i] = strsplit(paxdb$string_external_id[i], "\\.")[[1]][2]
}

step1 = merge(unfiltkds, paxdb, by.x = 'AGI', by.y = 'string_external_id')
names(step1)[names(step1) == 'abundance'] = 'paxabundance'
protinfo = read.delim('protdata.txt', sep = '|', header = FALSE)
protinfo$V1 = sapply(protinfo$V1, str_match, pattern = "AT.G\\d*")

pk = c()
pd = c()
for(i in 1:nrow(protinfo)){
  pk[i] = as.numeric(str_split(protinfo$V4[i], pattern = '=')[[1]][2])
  pd[i] = str_split(protinfo$V2[i], pattern = ':')[[1]][2]
}
protinfo$V4 = pk
protinfo$V2 = pd
names(protinfo) = c("AGI", "Symbol", "Genename", "protlen")
step2 = merge(step1, protinfo, by.x = 'AGI', by.y = 'AGI')
print(colnames(step2))
step3 = step2[,c("AGI","protein","seqsourcepr","Symbol","Genename","BINCODE","NAME","onebin","twobins","twobinnames","protlen","daylpfpr","daysdpr","dayrsdpr","daysepr","nightlpfpr","nightsdpr","nightrsdpr","nightsepr","daypeptidespr","nightpeptidespr","unidaypeptidespr","uninightpeptidespr","daylpfqt","daysdqt","dayrsdqt","dayseqt","nightlpfqt","nightsdqt","nightrsdqt","nightseqt","daypeptidesqt","nightpeptidesqt","unidaypeptidesqt","uninightpeptidesqt","p.value","stderr","daymean","nightmean", "rawdayFCP", "rawnightFCP", "dayFCP","nightFCP","dayKd","nightKd","dayKs","nightKs","paxabundance")]

averageAAweightinProtein = 110 #g/mol
average_weight_of_rosette = 0.14 #g @ 21 days
protein_content_of_rosette = average_weight_of_rosette * .02
aa_content_in_moles = protein_content_of_rosette/averageAAweightinProtein

step3$paxpropor = step3$paxabundance / 1000000 #taking the ppm out of the number making it just a proportion of the total number of proteins are here(should sum to 1 if all present)
propor_total_amino = sum(step3$paxpropor)
aarepresented = aa_content_in_moles * propor_total_amino #how many amino acids are in our measured protein pool
step3$len_by_paxpropor = step3$paxpropor * step3$protlen #proportion of total protein by length
sum_of_paxproportion_by_prot_length = sum(step3$len_by_paxpropor)
moles_aa_per_pax_by_length = aarepresented / sum_of_paxproportion_by_prot_length
step3$moles_aa_in_pool = moles_aa_per_pax_by_length *step3$len_by_paxpropor
step3$micromoles_aa_in_pool = step3$moles_aa_in_pool * 1000000

step3$daydegcost = step3$micromoles_aa_in_pool * step3$dayKd * 1.25 # abs number of amino acids * proportion of total degraded in 12h light * 1.25 (number of ATP per AA degraded)
step3$nightdegcost = step3$micromoles_aa_in_pool * step3$nightKd * 1.25
step3$daysynthcost = step3$micromoles_aa_in_pool * step3$dayKs * 5.25
step3$nightsynthcost = step3$micromoles_aa_in_pool * step3$nightKs *5.25
step3$daypepbondsformed = step3$micromoles_aa_in_pool * step3$dayKs
step3$nightpepbondsformed = step3$micromoles_aa_in_pool * step3$nightKs

coreset = step3[complete.cases(step3),]
coreset = subset(coreset, (coreset$daysdpr < 0.06 | coreset$dayrsdpr < 0.25) & coreset$unidaypeptidespr >2 & coreset$daypeptidespr > 4 & (coreset$nightsdpr < 0.06 | coreset$nightrsdpr < 0.25) & coreset$uninightpeptidespr > 2 & coreset$nightpeptidespr > 4)
coreset = coreset[!duplicated(coreset$AGI),]
coreset = coreset[!is.na(coreset[,'p.value']),]
coreset$daycost = coreset$daydegcost + coreset$daysynthcost
coreset$nightcost = coreset$nightdegcost + coreset$nightsynthcost
coreset$onebin = sapply(strsplit(as.character(coreset$NAME),"\\."), '[', 1)
coreset$onebin = as.factor(coreset$onebin)
#partydata = readRDS('3rdpartydata.rds')
#coreset = merge(coreset, partydata, by.x = 'AGI', by.y = 'AGI', all.x = TRUE)
coreset$kSDN = coreset$dayKs/coreset$nightKs
# ksdn = c()
# costdn = c()
# for(i in 1:nrow(coreset)){
#   if(coreset$dayKs[i] > coreset$nightKs[i]){ksdn[i] = (coreset$dayKs[i]/coreset$nightKs[i])-1}else{ksdn[i] = (-coreset$nightKs[i]/coreset$dayKs[i])+1}
#   if(coreset$daycost[i] > coreset$nightcost[i]){costdn[i] = (coreset$daycost[i]/coreset$nightcost[i])-1}else{costdn[i] = (-coreset$nightcost[i]/coreset$daycost[i])+1}
# }
# coreset$kSDN = ksdn
# coreset$costDN = costdn
coreset$costDN = coreset$daycost/coreset$nightcost
suba = read.csv(file = "Suba4-2019-8-21_12-9.csv")
suba = suba[,c('locus', 'location_consensus')]
coreset = merge(coreset, suba, by.x = 'seqsourcepr', by.y = 'locus')
print(colnames(coreset))
coreset = coreset[,c("seqsourcepr", "AGI", "protein", "Symbol", "Genename", "BINCODE", "NAME", "onebin", "twobins", "twobinnames", "location_consensus", "protlen", "daylpfpr", "daysdpr", "dayrsdpr", "daysepr", "nightlpfpr", "nightsdpr", "nightrsdpr", "nightsepr", "daypeptidespr", "nightpeptidespr", "unidaypeptidespr", "uninightpeptidespr", "daylpfqt", "daysdqt", "dayrsdqt",  "dayseqt", "nightlpfqt", "nightsdqt", "nightrsdqt", "nightseqt", "daypeptidesqt", "nightpeptidesqt", "unidaypeptidesqt", "uninightpeptidesqt", "p.value", "stderr", "daymean", "nightmean", "rawdayFCP", "rawnightFCP",  "dayFCP", "nightFCP", "dayKd", "nightKd", "dayKs", "nightKs", "paxpropor", "moles_aa_in_pool", "daypepbondsformed", "nightpepbondsformed", "daydegcost", "nightdegcost", "daysynthcost", "nightsynthcost", "daycost", "nightcost", "costDN", "kSDN")]
coreset$location_consensus = as.character(coreset$location_consensus)
for(i in 1:nrow(coreset)){
  if(grepl('^ATC.*', coreset[i,'seqsourcepr'])){coreset[i,'translated'] = 'plastid translated'}
  if(grepl('^AT[1-9].*', coreset[i,'seqsourcepr'])){coreset[i,'translated'] = 'cytosol translated'}
}

for(i in 1:nrow(coreset)){
  if(grepl('.*plastid.*', coreset[i,'location_consensus']) & grepl('^AT[1-9].*', coreset[i,'seqsourcepr'])){coreset[i,'location_consensus'] = 'nuclear encoded plastid'}
  if(grepl('.*cytosol.*', coreset[i,'location_consensus'])){coreset[i,'location_consensus'] = 'cytosol'}
  if(grepl('.*mitochondiron.*', coreset[i,'location_consensus'])){coreset[i,'location_consensus'] = 'mitochondrion'}
  if(grepl('^ATC.*', coreset[i,'seqsourcepr'])){coreset[i,'location_consensus'] = 'plastid encoded'}
  coreset[i,'Symbol'] = strsplit(coreset[i,'Symbol'],',')[[1]][1]
  coreset[i,'Symbol'] = strsplit(coreset[i,'Symbol'],'/')[[1]][1]
  if(grepl('^no *.', coreset[i,'Symbol'])){coreset[i,'Symbol'] = ' '}
}

coreset$location_consensus = as.factor(coreset$location_consensus)
#write.csv(coresetv2, file = 'coreset.csv')
# 
# tl = read.csv("tl.csv")
# tl[,1] = as.character(tl[,1])
# tl = tl[,c(1,14:17)]
# 
# tlaucday = c()
# tlaucnight = c()
# tlratio = c()
# #tldaymax = c()
# for(i in 1: nrow(tl)){
#   tlaucday[i] = area_under_curve(c(1,2,3), tl[i,c(2,3,4)])
#   tlaucnight[i] = area_under_curve(c(1,2,3), tl[i,c(4,5,2)])
#   tlratio[i] = tlaucday[i] / tlaucnight[i]
#   # if(tlaucday[i] > tlaucnight[i]){
#   #   tlratio[i] = (tlaucday[i]/tlaucnight[i])-1
#   #   tldaymax[i] = TRUE
#   # }else{
#   #   tlratio[i] = (-tlaucnight[i]/tlaucday[i])+1
#   #   tldaymax[i] = FALSE
#   # }
# }
# tx = read.csv("tx.csv")
# tx[,1] = as.character(tx[,1])
# tx = tx[,c(1,14:17)]
# txaucday = c()
# txaucnight = c()
# txratio = c()
# txdaymax = c()
# for(i in 1: nrow(tx)){
#   txaucday[i] = area_under_curve(c(1,2,3), tx[i,c(2,3,4)])
#   txaucnight[i] = area_under_curve(c(1,2,3), tx[i,c(4,5,2)])
#   txratio[i] = txaucday[i]/txaucnight[i]
#   # if(txaucday[i] > txaucnight[i]){
#   #   txratio[i] = (txaucday[i]/txaucnight[i])-1
#   #   txdaymax[i] = TRUE
#   # }else{
#   #   txratio[i] = (-txaucnight[i]/txaucday[i])+1
#   #   txdaymax[i] = FALSE
#   # }
# }
# tldid = data.frame(matrix(ncol = 4, nrow = nrow(tl)))
# tldid[,1] = as.character(tl[,1])
# tldid[,2] = tlaucday
# tldid[,3] = tlaucnight
# tldid[,4] = tlratio
# 
# txdid = data.frame(matrix(ncol = 4, nrow = nrow(tx)))
# txdid[,1] = as.character(tx[,1])
# txdid[,2] = txaucday
# txdid[,3] = txaucnight
# txdid[,4] = txratio
# 
# 
# colnames(tldid) = c("AGI", "tlaucday", "tlaucnight", "tlratio")
# missratlcore = merge(coreset, tldid, by = "AGI", all.x = TRUE)
# 
# colnames(txdid) = c("AGI", "txaucday", "txaucnight", "txratio")
# missratltxcore = merge(missratlcore, txdid, by = "AGI", all.x = TRUE)

tx = read.csv("tx.csv")
tl = read.csv("tl.csv")
tl[,1] = as.character(tl[,1])
tl = tl[,c(1,14:17)]
tx[,1] = as.character(tx[,1])
tx = tx[,c(1,14:17)]
colnames(tl)[1] = "AGI"
colnames(tx)[1] = "AGI"
tltx = merge(tl, tx, by.x='AGI', by.y='AGI')
tltx[10:13] = tltx[2:5] * tltx[6:9]
tlxtx = tltx[,c(1,10,11,12,13)] 
colnames(tlxtx) = c("agi", "6am", "12pm", "6pm", "12am")
library(bayestestR)
tlxtxaucday = c()
tlxtxaucnight = c()
tlxtxratio = c()
for(i in 1: nrow(tlxtx)){
  tlxtxaucday[i] = area_under_curve(c(1,2,3), tlxtx[i,c(2,3,4)])
  tlxtxaucnight[i] = area_under_curve(c(1,2,3), tlxtx[i,c(4,5,2)])
  tlxtxratio[i] = tlxtxaucday[i]/tlxtxaucnight[i]
}
tlxtxout = tlxtx
colnames(tlxtxout) = c("agi", "tltxday", "tltxnight", "tltxratio")
tlxtxout = tlxtxout[1:4]
tlxtxout$agi = tlxtx$agi
tlxtxout$tltxday = tlxtxaucday
tlxtxout$tltxnight = tlxtxaucnight
tlxtxout$tltxratio = tlxtxratio

colnames(tlxtxout) = c("AGI", "missraday", "missranight", "missraratio")
missratlxtxcore = merge(coreset, tlxtxout, by = "AGI", all.x = TRUE)

blasing = read.csv("blasing1212rna.csv")
colnames(blasing) = c("AGI", "00h", "04h", "08h", "12h", "16h", "20h")
blasingaucday = c()
blasingaucnight = c()
blasingratio = c()
for(i in 1: nrow(blasing)){
  blasingaucday[i] = area_under_curve(c(1,2,3,4), blasing[i,c(2,3,4,5)])
  blasingaucnight[i] = area_under_curve(c(1,2,3,4), blasing[i,c(5,6,7,2)])
  blasingratio[i] = blasingaucday[i]/blasingaucnight[i]
  # if(blasingaucday[i] > blasingaucnight[i]){
  #   blasingratio[i] = (blasingaucday[i]/blasingaucnight[i])-1
  #   blasingdaymax[i] = TRUE
  # }else{
  #   blasingratio[i] = (-blasingaucnight[i]/blasingaucday[i])+1
  #   blasingdaymax[i] = FALSE
  # }
}
blasingauc = as.data.frame(cbind(as.numeric(blasingaucday), as.numeric(blasingaucnight), as.numeric(blasingratio)))
blasingauc[,c(2,3,4)]= blasingauc
blasingauc[,1] = blasing$AGI
colnames(blasingauc) = c("AGI", "blasingaucday", "blasingaucnight", "blasingratio") 
blasingcore = merge(missratlxtxcore, blasingauc, by = "AGI", all.x = TRUE)

seaton = read.csv("seaton.csv")
seatonratio = seaton$mean_18h/seaton$mean_6h
# for(i in 1: nrow(seaton)){
#   if(seaton$mean_18h[i] > seaton$mean_6h[i]){
#     seatonratio[i] = (seaton$mean_18h[i]/seaton$mean_6h[i])-1
#   }else{
#     seatonratio[i] = (-seaton$mean_6h[i]/seaton$mean_18h[i])+1
#   }
# }
seatonout = cbind.data.frame(seaton$ï..locus_id, seatonratio)
seatonout$`seaton$ï..locus_id` = as.character(seatonout$`seaton$ï..locus_id`)
colnames(seatonout) = c("AGI", "seatonratio")
seatoncore = merge(blasingcore, seatonout, by = "AGI", all.x = TRUE)

piques = read.csv("piques.csv")
piques$piquesratio = piques$Piques_lightmol.h.1.g.1FW/piques$Piques_darkmol.h.1.g.1FW
piques = piques[,c(1,4)]
colnames(piques) = c("AGI", "piquesratio")
piquescore = merge(seatoncore, piques, by = "AGI", all.x = TRUE)

saveRDS(piquescore, file = "coresetv3.rds")
