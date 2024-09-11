#### Script to Run GSI on Inseason Samples from BBSRI ####
library(rubias)
library(coda)
library(doFuture)
library(doRNG)

CrntYr <- format(Sys.Date(),"%Y") #Get the current year
dir.create(file.path(paste("Output/",CrntYr,"/",CrntSW,sep="")))


#functions
statweek = function(dates, format="%d-%m-%Y", ...) {
  # convert to Date
  dates = as.Date(dates, format=format, ...) 
  # get correction for the first week of the year (0 if 1-Jan not a Sunday)
  firstweek = 1 - as.numeric(format(as.Date(cut(dates, "year")), "%U")) 
  output = as.numeric(format(dates, "%U")) + firstweek
  return(output)
}



### Process AKRO Files ####
AKRO_Files <- list.files(file.path(paste("./MainData/AKRO/",CrntYr,sep="")),
                         include.dirs = FALSE) #list AKRO files sent by Glenn

#Find the file to process
FileWeeks <- AKRO_Files %>%
  str_split(.,"_")%>%
  lapply(.,"[[",1)%>%
  unlist() %>% 
  as.Date(.,format="%Y%m%d")%>%
  format(.,"%d-%m-%Y")%>%
  statweek(.)

File2Read <- AKRO_Files[which.min(abs(CrntSW-FileWeeks)) ]

AKRO_Obs <- readxl::read_xls(file.path(paste("./MainData/AKRO/",CrntYr,"/",File2Read,sep=""))) %>%
  mutate(DATE = format(DELIVERY_END_DATE,"%d-%m-%Y"))%>%
  mutate(SW = statweek(DATE)) %>%
  filter(SW <= CrntSW)

AKRO_ObsSum <- AKRO_Obs %>%
  group_by(SW,NMFS_AREA) %>%
  summarize(nChum = sum(CHUM_RETENTION_COUNT,na.rm=T),
            nDeliveries = n(),
            nUniqVess = length(unique(DELIVERING_VESSEL_PERMIT)))%>%
  mutate(nChumPerDeliv = nChum/nDeliveries)%>%
  filter(nUniqVess >=3)

#### Process BBSRI Files ####

BBSRI_Dlvry <- readxl::read_xlsx(file.path(paste("./MainData/BBSRI/",CrntYr,"/SW",CrntSW,"/","Catch_Sampling_Master_SW",CrntSW,".xlsx",sep="")),sheet = 'DELIVERYINFO') %>%
  mutate(DATE = as.Date(as.character(DATE),"%Y-%m-%d"))%>%
  mutate(DATE2 = format(DATE,"%d-%m-%Y"))%>%
  mutate(SW = statweek(DATE2))%>%
  rename('PROCESSOR_PERMIT' = 'PROCESSOR',
         'OFFLOAD_NUMBER' = 'OFFLOAD_NO')

#Bind to obs data
AKRO_BBSRI_Dlvry <- AKRO_Obs %>%
  full_join(.,BBSRI_Dlvry, by = c('CRUISE','PROCESSOR_PERMIT',"OFFLOAD_NUMBER","SW"))%>% #don't merge by NMFS area, use NMFS_AREA.x from AKRO data in subsequent code
  filter(SW == CrntSW)

if(any(is.na(AKRO_BBSRI_Dlvry$CRUISE) & (AKRO_BBSRI_Dlvry %>% filter(is.na(CRUISE)) %>% select(CHUM_LANDED) %>% sum(.,na.rm=T) %>% {.>0}))==TRUE) {
  "WARNING: Positive chum catches for cruises not in AKRO data"
}else{ AKRO_BBSRI_Dlvry <- AKRO_BBSRI_Dlvry %>% filter(!is.na(CRUISE))
}


#### Process Genotype Files ####
GenoFiles <- list.files(file.path(paste("./MainData/BBSRI/",CrntYr,"/SW",CrntSW,"/",sep="")))%>%
  grep(pattern='Chip|Run',.,value=T)

GenoMetaDat <- readxl::read_xlsx(file.path(paste("./MainData/BBSRI/",CrntYr,"/SW",CrntSW,"/","Catch_Sampling_Master_SW",CrntSW,".xlsx",sep="")),
                                 sheet = "BIOSAMPLING") 

Genos_Fldm<-lapply(1:length(GenoFiles),function(x) read_delim(file.path(paste("./MainData/BBSRI/",CrntYr,"/SW",CrntSW,"/",GenoFiles[x],sep="")),delim=",",skip = 16,col_names = F) %>%
                     rename('Chamber'= 1 ,
                            'Locus' = 2,
                            'Allele_x' = 3,
                            'Allele_y' = 4,
                            'Sample_Name' = 5,
                            'Sample_Type' = 6,
                            'CallInfo_Auto' = 7,
                            'Call_Confidence' = 8,
                            'CallInfo_Final' = 9,
                            'CallInfo_Converted' = 10,
                            'Intensity_x' = 11,
                            'Intensity_y' = 12) %>%
                     mutate(Chip = str_split(string = GenoFiles[x],pattern="_|Chip") %>% lapply(.,"[[",4)%>% str_trim()%>% unlist()) %>%
                     filter(!grepl(pattern="NTC",Sample_Name)))%>%
  bind_rows(.)%>%
  select(1:12,Chip)

Genos_Fldm_Genotype <- Genos_Fldm %>%
  mutate(CallInfo_Converted = recode(CallInfo_Converted,
                                     'Invalid' = "NA:NA",
                                     'No Call' = "NA:NA",
                                     'NTC' = "NTC:NTC"))%>%
  select(Sample_Name,Chip, Locus,CallInfo_Converted) %>%
  reshape2::dcast(Sample_Name + Chip ~ Locus , value.var= c("CallInfo_Converted"))

rubiasCols <- data.frame(sample_type = "mixture",
                         repunit = NA,
                         collection = paste("SW",CrntSW,sep=""),
                         indiv = Genos_Fldm_Genotype$Sample_Name)

Genos_Fldm_Alleles <- lapply(3:ncol(Genos_Fldm_Genotype),function(x)
  str_split(Genos_Fldm_Genotype[,x],pattern=":") %>%
    {cbind(lapply(.,"[[",1) %>% unlist(),lapply(.,"[[",2) %>% unlist())}
)%>%
  do.call(cbind,.)%>%
  {cbind(rubiasCols,.)}

colnames(Genos_Fldm_Alleles)[5:(ncol(Genos_Fldm_Alleles))]<-paste(rep(colnames(Genos_Fldm_Genotype)[-c(1:2)],each=2),c("1","2"),sep="_")

Genos_Fldm_Alleles[Genos_Fldm_Alleles=="NA"]<-NA


#check for samples genotyped multiple times
if(any(duplicated(Genos_Fldm_Alleles$indiv %>% {.[.!="NTC"]}))){
  paste("WARNING: Duplicated samples in genotype file (",paste(Genos_Fldm_Alleles$indiv[duplicated(Genos_Fldm_Alleles$indiv,fromLast = T)] %>% {.[.!="NTC"]},collapse=","),").",sep="")
}   

nDups<- sum(duplicated(Genos_Fldm_Alleles$indiv %>% {.[.!="NTC"]}))

Baseline <- read.csv("./MainData/ChumBaseline_ABL84_382pops_42071inds_KotzSnd.csv")

Mixtures <- Genos_Fldm_Alleles

Mixtures[,(2:ncol(Mixtures))]<-lapply(((2:ncol(Mixtures))), function(x) as.character(Mixtures[,x]))

Baseline[,2:ncol(Baseline)]<-lapply(2:ncol(Baseline), function(x) as.character(Baseline[,x]))

Mixtures[Mixtures==0]<-NA

# Confirm that the loci names match  
BaseLoci <- colnames(Baseline)[-(1:4)] %>%
  gsub(pattern="\\.1$|_1$|_2$",replacement="",.) %>%
  gsub(pattern="\\.",replacement="_",.) %>%
  tolower()

colnames(Baseline)[-(1:4)] <- paste(BaseLoci,c("",".1"),sep="")

MixLoci <- colnames(Mixtures)[-(1:4)] %>%
  gsub(pattern="_1$|_1.1$|_2$",replacement="",.) %>%
  gsub(pattern="-",replacement="_",.) %>%
  tolower()

colnames(Mixtures)[-(1:4)] <- paste(MixLoci,c("",".1"),sep="")

Loci2RM <- MixLoci[!MixLoci %in% BaseLoci][c(T,F)]

Mixtures <- Mixtures %>%
  select(-(contains(Loci2RM)))

BaseLoci <- colnames(Baseline)[-(1:4)]
MixLoci <- colnames(Mixtures)[-(1:4)]

Test<-cbind(BaseLoci,MixLoci)

if(any(BaseLoci != MixLoci)){ 
  LociOrder <- lapply(BaseLoci[c(T,F)],function(BL) which(MixLoci %in% BL)) %>%
    unlist() +4
} else {
  LociOrder <- 5:(84*2+4)
}

Mixtures<-Mixtures[,c(1:4,LociOrder)]

if(any(colnames(Baseline)!= colnames(Mixtures))){
  "WARNING!!!! rubias baseline and mixture files have different locus order."
}

#remove NTC sample
Mixtures <- Mixtures %>% 
  filter(indiv != "NTC")

#Samples with low genotyping success?
nloci<- (ncol(Mixtures)-4) / 2

BadInds <- Mixtures %>%
  mutate(CntNA = rowSums(is.na(select(., starts_with("oke")))),
         GTrate = round(1-(CntNA / (nloci*2)),2))%>%
  filter(GTrate<0.8)

if(nrow(BadInds)>0){
  Mixtures<- Mixtures %>%
    filter(!(indiv%in%BadInds$indiv))
}


CloseMatches <- suppressWarnings( rubias::close_matching_samples(D = Mixtures %>% 
                                                                   as_tibble() %>% 
                                                                   mutate(sample_type = "reference", 
                                                                          collection = "CrntSW",
                                                                          repunit = "S"),
                                                                 gen_start_col = 5,
                                                                 min_frac_non_miss = 0.8,
                                                                 min_frac_matching = 0.9) )


#starting proportions for the chains
### Starting proportions of each of the chains is 95% of 1 reporting group

Baseline_Reps<-Baseline %>%
  group_by(repunit) %>%
  select(repunit,collection) %>%
  unique()

pi_init_List<-lapply(1:length(unique(Baseline_Reps$repunit)),function(x)
  Baseline_Reps %>%
    mutate(RepPi = ifelse(repunit==unique(Baseline_Reps$repunit)[x],1,2))%>%
    mutate(pi_init = ifelse(RepPi==1,
                            0.95/(Baseline_Reps %>%
                                    filter(repunit==unique(Baseline_Reps$repunit)[x]) %>%
                                    nrow()),
                            0.05/(Baseline_Reps %>%
                                    filter(repunit!=unique(Baseline_Reps$repunit)[x]) %>%
                                    nrow())) ) %>%
    ungroup() %>%
    select(collection,pi_init))

### Prior on stock proportions GCg
# A 1/k prior on baseline populations puts a larger probability on 
# reporting groups with many populations in the baseline. So we use a 
# 1/GCg prior to have a more uniform prior over the reporting groups 

RepColls<-unique(Baseline[,c(2,3)])

G<-length(unique(RepColls$repunit))
PriorRep<-RepColls%>%
  group_by(repunit)%>%
  summarise(nPops=length(unique(collection)))%>%
  mutate(GCg=nPops*G)%>%
  mutate(pi=1/GCg)

prior_GCg<-RepColls%>%
  mutate(pi_param=as.numeric(plyr::mapvalues(x=repunit,
                                             from=PriorRep$repunit,
                                             to=PriorRep$pi,)))%>%
  select(collection,pi_param)

AnalyNames <- Mixtures$collection %>% unlist() %>% unique()


#### Conditional Model #### 

dir.create(file.path(paste("Output/",CrntYr,"/",CrntSW,"/Rubias_NonBootstrapModel",sep="")))

MCMCreps <- 100000
BurnIn <- MCMCreps/2
chains <- 7

RubiasRes<-list()
GR_diag<-list()
GR_Mat_List<-list()

set.seed(11)  
Chn2Samp<-sample(1:chains,1)  

for (mx in 1:length(AnalyNames)){
  
  MixTemp<-Mixtures %>%
    filter(collection == AnalyNames[mx]) #subset on the mix
  
  
  #Check for duplicates
  if(any(duplicated(MixTemp[,-c(1:4)]))){
    cat("Warning there are duplicated genotypes, one duplicate removed!!!")
    DupRows<-which(duplicated(MixTemp[,-c(1:4)]))
    MixTemp<-MixTemp[-DupRows,]
  }
  
  #check to see if any loci are completely ungenotyped
  UnGenoLoci<- MixTemp %>%
    reshape2::melt(.,id.vars=c("sample_type","repunit","collection","indiv")) %>%
    mutate(variable = gsub(pattern="\\.1$",replacement="",variable)) %>%
    group_by(variable) %>%
    summarize(NAct = sum(is.na(value)),
              NAprop = sum(is.na(value))/n())
  Loci2Drop<-UnGenoLoci %>% 
    filter(NAprop == 1) %>%
    select(variable) %>%
    unlist() %>%
    as.vector()
  
  MixTemp <- MixTemp %>%
    select(!starts_with(Loci2Drop))
  
  BaseTemp <- Baseline%>%
    select(!starts_with(Loci2Drop))
  
  registerDoFuture()
  plan(multisession) #multisession on windows
  
  X <- 1:chains
  
  mix_est_list <- foreach(x = X) %dorng% { #dorng vs dopar
    mix_est <- infer_mixture(reference = BaseTemp, 
                             mixture = MixTemp, 
                             gen_start_col = 5,
                             reps = MCMCreps, 
                             burn_in = BurnIn,
                             pi_init = pi_init_List[[x]],                          
                             pi_prior = prior_GCg)#,
    #method="PB")
    mix_est
    
  }
  
  MCMC_chains<-lapply(1:chains,function(x)
    mix_est_list[[x]]$mix_prop_traces %>%
      filter(sweep > BurnIn) %>%
      group_by(sweep, repunit) %>%
      summarise(repprop = sum(pi))%>%
      reshape2::dcast(.,formula = sweep ~ repunit , value.var = "repprop") %>%
      select(-sweep)%>%
      coda::mcmc(.)
  )%>%
    coda::mcmc.list()
  
  GR_diag <- coda::gelman.diag(MCMC_chains,autoburnin = F, multivariate = F)
  
  GR_Mat_List[[mx]]<-GR_diag$psrf %>%
    as.data.frame() %>%
    mutate(repunit = row.names(.))%>%
    mutate(Analysis=AnalyNames[mx])
  
  RG_mix_ests <- mix_est_list[[Chn2Samp]]$mixing_proportions %>%
    group_by(mixture_collection, repunit) %>%
    summarise(repprop = sum(pi)) 
  
  trace_subset <- mix_est_list[[Chn2Samp]]$mix_prop_traces %>%
    filter(sweep > BurnIn) %>%
    group_by(sweep, repunit) %>%
    summarise(repprop = sum(pi)) 
  
  CI_rub <- trace_subset %>%
    group_by(repunit) %>%
    summarise(loCI = quantile(repprop, probs = 0.025),
              hiCI = quantile(repprop, probs = 0.975))
  
  RubiasRes[[mx]]<-cbind(RG_mix_ests[,2:3],CI_rub[,2:3])%>%
    mutate(Analysis=AnalyNames[mx])%>%
    left_join(.,GR_Mat_List[[mx]],by=c("Analysis","repunit"))
  
  saveRDS(mix_est_list[[Chn2Samp]],file=file.path(paste("./Output/",CrntYr,"/",CrntSW,"/","Rubias_NonBootstrapModel/",AnalyNames[mx],".rds",sep="")))
  
  write.csv(RubiasRes[[mx]],
            file=file.path(paste("./Output/",CrntYr,"/",CrntSW,"/","Rubias_NonBootstrapModel/",AnalyNames[mx],"_Results.csv",sep="")),
            quote=F,
            row.names=F)
  
  plan(sequential)
  
}#over mx mixtures

RubiasResMat<-do.call(rbind,RubiasRes)



#### Bootstrap Model #### 
dir.create(file.path(paste("Output/",CrntYr,"/",CrntSW,"/Rubias_BootstrapModel",sep="")))
#New Code to run BS models  
MCMCreps <- 100000
BURNin <- MCMCreps/2
BootStrapIts<-100  
registerDoFuture()
plan(multisession)


X <- 1:length(AnalyNames)   

ModelRes_List <- foreach(x = X) %dorng% {
  
  MixTemp <- Mixtures%>%
    
    filter(collection==!!AnalyNames[[x]])
  
  if(any(duplicated(MixTemp[,-c(1:4)]))){
    
    DupRows<-which(duplicated(MixTemp[,-c(1:4)]))
    MixTemp<-MixTemp[-DupRows,]
  }
  
  mix_est <- infer_mixture(reference = Baseline, 
                           mixture = MixTemp, 
                           gen_start_col = 5,
                           reps = MCMCreps, 
                           burn_in = BURNin,
                           pi_init = pi_init_List[[1]],
                           pi_prior = prior_GCg,
                           method = "PB",
                           pb_iter = BootStrapIts)
  
  RG_mix_ests <- mix_est$mixing_proportions %>%
    group_by(repunit) %>%
    summarise(repprop = sum(pi)) %>% 
    left_join(mix_est$bootstrapped_proportions) 
  
  BS_diff<- mix_est$mixing_proportions %>%
    group_by(repunit) %>%
    summarise(repprop = sum(pi)) %>% 
    left_join(mix_est$bootstrapped_proportions)%>%
    group_by(repunit)%>%
    summarize(BSdiff = bs_corrected_repunit_ppn - repprop)
  
  trace_subset <- mix_est$mix_prop_traces %>%
    filter(sweep > BURNin) %>%
    group_by(sweep, repunit) %>%
    summarise(repprop = sum(pi))%>%
    mutate(DiffBS = as.numeric(plyr::mapvalues(from=BS_diff$repunit,
                                               to=BS_diff$BSdiff,
                                               x=repunit)))%>%
    mutate(BS_repprop = repprop + DiffBS) 
  
  
  trace_subset[which(trace_subset$BS_repprop < 0),5]<-0
  trace_subset[which(trace_subset$BS_repprop > 1),5]<-1
  
  BS_CI_rub <- trace_subset %>%
    group_by(repunit) %>%
    summarise(loCI_BS = quantile(BS_repprop, probs = 0.025),
              hiCI_BS = quantile(BS_repprop, probs = 0.975),
              sd_BS = sd(BS_repprop),
              median_BS = median(BS_repprop))
  
  CI_rub <- trace_subset %>%
    group_by(repunit) %>%
    summarise(loCI = quantile(repprop, probs = 0.025),
              hiCI = quantile(repprop, probs = 0.975),
              sd = sd(repprop),
              median = median(repprop))
  
  BSMod<-left_join(RG_mix_ests,CI_rub,by="repunit")%>%
    relocate(mixture_collection,.before=repunit)%>%
    relocate(starts_with("bs_corrected_repunit"),.after=median)%>%
    left_join(.,BS_CI_rub,by="repunit")
  
  P0 <- trace_subset %>%
    group_by(repunit)%>%
    mutate(LT_1in1Mill = ifelse(BS_repprop < 0.0000005,1,0)) %>%
    summarise(P0 = sum(LT_1in1Mill)/(MCMCreps-BURNin))
  
  BSMod<-right_join(BSMod,P0,by="repunit")
  
  write.csv(x=BSMod,
            file=file.path(paste("./Output/",CrntYr,"/",CrntSW,"/","Rubias_BootstrapModel/",AnalyNames[[x]],"BSmodRes.csv",sep="")),
            row.names = F) 
  
  saveRDS(mix_est,file=file.path(paste("./Output/",CrntYr,"/",CrntSW,"/","Rubias_BootstrapModel/",AnalyNames[[x]],".rds",sep="")))
  
  BSMod
}


