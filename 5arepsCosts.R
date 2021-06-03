setwd("d:/diurnalfive")
library(stringr)

unfiltkds = readRDS(file = "kdksreps.rds")
paxdb = read.table("paxdb.txt", header = TRUE)

paxdb$internal_id = NULL
paxdb$string_external_id = as.character(paxdb$string_external_id)

for(i in 1:nrow(paxdb)){
  paxdb$string_external_id[i] = strsplit(paxdb$string_external_id[i], "\\.")[[1]][2]
}

step1 = list()
for(i in 1:length(unfiltkds)){
  step1[[i]] = merge(unfiltkds[[i]], paxdb, by.x = 'AGI', by.y = 'string_external_id')
  names(step1[[i]])[names(step1[[i]]) == 'abundance'] = 'paxabundance'
}

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

step2 = list()
step3 = list()
coreset = list()
for(j in 1:length(step1)){
  #steps below determie the total number of amino acids in a rosette by weighing the whole thing and saying ~2% is protein
  averageAAweightinProtein = 110 #g/mol
  average_weight_of_rosette = 0.14 #g @ 21 days
  protein_content_of_rosette = average_weight_of_rosette * .02
  aa_content_in_moles = protein_content_of_rosette/averageAAweightinProtein
  
  #steps below determine what proportion of these amino acids are in each protein group pool
  step2[[j]] = merge(step1[[j]], protinfo, by.x = 'AGI', by.y = 'AGI')
  step3[[j]] = step2[[j]][,c("AGI","protein","seqsourcepr","Symbol","Genename","BINCODE","NAME","onebin","twobins","twobinnames","protlen","daylpfpr","daysdpr","dayrsdpr","daysepr","nightlpfpr","nightsdpr","nightrsdpr","nightsepr","daypeptidespr","nightpeptidespr","unidaypeptidespr","uninightpeptidespr","daylpfqt","daysdqt","dayrsdqt","dayseqt","nightlpfqt","nightsdqt","nightrsdqt","nightseqt","daypeptidesqt","nightpeptidesqt","unidaypeptidesqt","uninightpeptidesqt","p.value","stderr","daymean","nightmean", "rawdayFCP", "rawnightFCP", "dayFCP","nightFCP","dayKd","nightKd","dayKs","nightKs","paxabundance")]
  step3[[j]]$paxpropor = step3[[j]]$paxabundance / 1000000 #tkaing the ppm out of the number making it just a proportion of the total number of proteins are here(should sum to 1 if all present)
  #then find how many amino acids from the gross pool measuremetn are represented by this numbner of proteins by adding them together
  propor_total_amino = sum(step3[[j]]$paxpropor)
  aarepresented = aa_content_in_moles * propor_total_amino #how many amino acids are in our measured protein pool
  step3[[j]]$len_by_paxpropor = step3[[j]]$paxpropor * step3[[j]]$protlen #proportion of total protein by length
  sum_of_paxproportion_by_prot_length = sum(step3[[j]]$len_by_paxpropor)
  moles_aa_per_pax_by_length = aarepresented / sum_of_paxproportion_by_prot_length
  step3[[j]]$moles_aa_in_pool = moles_aa_per_pax_by_length *step3[[j]]$len_by_paxpropor
  step3[[j]]$micromoles_aa_in_pool = step3[[j]]$moles_aa_in_pool * 1000000
  
  step3[[j]]$daydegcost = step3[[j]]$micromoles_aa_in_pool * step3[[j]]$dayKd * 1.25 # abs number of amino acids * proportion of total degraded in 12h light * 1.25 (number of ATP per AA degraded)
  step3[[j]]$nightdegcost = step3[[j]]$micromoles_aa_in_pool * step3[[j]]$nightKd * 1.25
  step3[[j]]$daysynthcost = step3[[j]]$micromoles_aa_in_pool * step3[[j]]$dayKs * 5.25
  step3[[j]]$nightsynthcost = step3[[j]]$micromoles_aa_in_pool * step3[[j]]$nightKs *5.25
  step3[[j]]$daypepbondsformed = step3[[j]]$micromoles_aa_in_pool * step3[[j]]$dayKs
  step3[[j]]$nightpepbondsformed = step3[[j]]$micromoles_aa_in_pool * step3[[j]]$nightKs

  coreset[[j]] = step3[[j]][complete.cases(step3[[j]]),]
  coreset[[j]] = coreset[[j]][!is.na(coreset[[j]][,'p.value']),]
  coreset[[j]]$daycost = coreset[[j]]$daydegcost + coreset[[j]]$daysynthcost
  coreset[[j]]$nightcost = coreset[[j]]$nightdegcost + coreset[[j]]$nightsynthcost
  coreset[[j]]$onebin = sapply(strsplit(as.character(coreset[[j]]$NAME),"\\."), '[', 1)
  coreset[[j]]$onebin = as.factor(coreset[[j]]$onebin)
  partydata = readRDS('3rdpartydata.rds')
  coreset[[j]] = merge(coreset[[j]], partydata, by.x = 'AGI', by.y = 'AGI', all.x = TRUE)
  coreset[[j]]$kSDN = coreset[[j]]$dayKs/coreset[[j]]$nightKs
  coreset[[j]]$costDN = coreset[[j]]$daycost/coreset[[j]]$nightcost
  suba = read.csv(file = "Suba4-2019-8-21_12-9.csv")
  suba = suba[,c('locus', 'location_consensus')]
  coreset[[j]] = merge(coreset[[j]], suba, by.x = 'seqsourcepr', by.y = 'locus')
  coreset[[j]] = coreset[[j]][,c("seqsourcepr", "AGI", "protein", "Symbol", "Genename", "BINCODE", "NAME", "onebin", "twobins", "twobinnames", "location_consensus", "protlen", "daylpfpr", "daysdpr", "dayrsdpr", "daysepr", "nightlpfpr", "nightsdpr", "nightrsdpr", "nightsepr", "daypeptidespr", "nightpeptidespr", "unidaypeptidespr", "uninightpeptidespr", "daylpfqt", "daysdqt", "dayrsdqt",  "dayseqt", "nightlpfqt", "nightsdqt", "nightrsdqt", "nightseqt", "daypeptidesqt", "nightpeptidesqt", "unidaypeptidesqt", "uninightpeptidesqt", "p.value", "stderr", "daymean", "nightmean", "rawdayFCP", "rawnightFCP",  "dayFCP", "nightFCP", "dayKd", "nightKd", "dayKs", "nightKs", "paxpropor", "moles_aa_in_pool", "daypepbondsformed", "nightpepbondsformed", "daydegcost", "nightdegcost", "daysynthcost", "nightsynthcost", "daycost", "nightcost", "costDN", "missraTADN", "missraTSDN", "TAD", "TAN", "blasingDN", "seatonLS", "kSDN")]
  coreset[[j]]$location_consensus = as.character(coreset[[j]]$location_consensus)
  for(i in 1:nrow(coreset[[j]])){
    if(grepl('^ATC.*', coreset[[j]][i,'seqsourcepr'])){coreset[[j]][i,'translated'] = 'plastid translated'}
    if(grepl('^AT[1-9].*', coreset[[j]][i,'seqsourcepr'])){coreset[[j]][i,'translated'] = 'cytosol translated'}
  }
  
  for(i in 1:nrow(coreset[[j]])){
    if(grepl('.*plastid.*', coreset[[j]][i,'location_consensus']) & grepl('^AT[1-9].*', coreset[[j]][i,'seqsourcepr'])){coreset[[j]][i,'location_consensus'] = 'nuclear encoded plastid'}
    if(grepl('.*cytosol.*', coreset[[j]][i,'location_consensus'])){coreset[[j]][i,'location_consensus'] = 'cytosol'}
    if(grepl('.*mitochondiron.*', coreset[[j]][i,'location_consensus'])){coreset[[j]][i,'location_consensus'] = 'mitochondrion'}
    if(grepl('^ATC.*', coreset[[j]][i,'seqsourcepr'])){coreset[[j]][i,'location_consensus'] = 'plastid encoded'}
    coreset[[j]][i,'Symbol'] = strsplit(coreset[[j]][i,'Symbol'],',')[[1]][1]
    coreset[[j]][i,'Symbol'] = strsplit(coreset[[j]][i,'Symbol'],'/')[[1]][1]
    if(grepl('^no *.', coreset[[j]][i,'Symbol'])){coreset[[j]][i,'Symbol'] = ' '}
  }
  coreset[[j]]$location_consensus = as.factor(coreset[[j]]$location_consensus)
}

saveRDS(coreset, file = "repscoreset.rds")

