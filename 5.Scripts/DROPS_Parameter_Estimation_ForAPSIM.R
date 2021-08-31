
## FUNCTIONS : 
source('ThermalTime.R')

## Genetic groups ----
Groups <- read.table(paste0(PATHS[3],'/Admixture6Group.csv'), 
                     sep=";",                                           
                     header=TRUE,dec=".",
                     row.names=NULL,
                     as.is=FALSE, 
                     fill=TRUE)
GenGroups <- subset(Groups,select=c('Hybrid_ID','group'))

## ZA16 Phenology ----
Phenology <- read.table(paste0(PATHS[3],'/ZA16_Phenology_2016-06-28.csv'), 
                        sep=";",                                           
                        header=TRUE,dec=".",
                        row.names=NULL,
                        as.is=FALSE, 
                        fill=TRUE)

Pheno <- subset(Phenology,select=c('Hybrid_ID','Scenario','doy','visi_l','ligul_l'))

Meteo <- read.table(paste0(PATHS[3],'/Meteo_ZA16.csv'), 
                    sep=";",                                           
                    header=TRUE,dec=".",
                    row.names=NULL,
                    as.is=FALSE, 
                    fill=TRUE)

Met <- subset(Meteo,date>=138,select=c('date','heure','tc.moy'))
Met$datehour <- (Met$date+(round(Met$heure/100)/24)+((Met$heure)-round(Met$heure/100)*100)/(60*24))-138
Tht <- ThermalTime('Apsim',Met$date,Met$tc.moy)
Tht$doy <- unique(Met$date)
Pheno_Tht <- merge(Pheno,Tht,by=c('doy'))
PhenoTht <- subset(Pheno_Tht,visi_l!=0 & visi_l<14 & Scenario =='WW' )

OutPhenoPanel <- data.frame(Phyllo = numeric(),Ntip_emerg = numeric(),
                            Ligules = numeric(),Tht_1st_Lig = numeric(),
                            Hybrid_ID=character())

## Calculate Phyllo/Ligules : ----
for (geno in unique(PhenoTht$Hybrid_ID))
{
  DataToTreat <- subset(PhenoTht,Hybrid_ID==geno)
  MyModel_Vis <- lm(data=DataToTreat,formula='visi_l~cumtht')
  MyModel_Lig <- lm(data=DataToTreat,formula='ligul_l~cumtht')
  
  ToBind <-  data.frame(Phyllo = MyModel_Vis$coefficients[2],Ntip_emerg = MyModel_Vis$coefficients[1],
                        Ligules = MyModel_Lig$coefficients[2],
                        Tht_1st_Lig = (1-MyModel_Lig$coefficients[1])/MyModel_Lig$coefficients[2] + 70 ,
                        Hybrid_ID=geno)
  OutPhenoPanel <- rbind(OutPhenoPanel,ToBind)
}


## Correspondance Hybrid_ID vs Code_ID : ----
Hybrids <- read.table(paste0(PATHS[3],'/Hybrids.csv'), 
                      sep=";",                                           
                      header=TRUE,dec=".",
                      row.names=NULL,
                      as.is=FALSE, 
                      fill=TRUE)

## Nfinal: ----
Nfinal <- read.table(paste0(PATHS[3],'/Nfinal.csv'), 
                     sep=";",                                           
                     header=TRUE,dec=".",
                     row.names=NULL,
                     as.is=FALSE, 
                     fill=TRUE)
Nfinal$Nf <- round(Nfinal$NbFF) 
colnames(Nfinal) <- c('Code_ID',"NbFF","Nf")
NfMax <- merge(Nfinal,Hybrids,by=c('Code_ID'))


## Sensitivity to PSI : ----
Sensitivity_PSI <- read.table(paste0(PATHS[3],'/Sensitivity.csv'), 
                              sep=",",                                           
                              header=TRUE,dec=".",
                              row.names=NULL,
                              as.is=FALSE, 
                              fill=TRUE)
SensiPSI <- subset(Sensitivity_PSI,select=c('Genotype','slope'))
colnames(SensiPSI) <- c('Hybrid_ID','c_norm')


## A max : ---- 
a_max <- read.table(paste0(PATHS[3],'/a.csv'), 
                    sep=",",                                           
                    header=TRUE,dec=".",
                    row.names=NULL,
                    as.is=FALSE, 
                    fill=TRUE)

a <- subset(a_max,select=c('Genotype','a'))
a$a_per <- a$a / a$a[which(a$Genotype=='B73_H')] 
a$a_dyn <- a$a_per*4.25 # a B73_H on phenodyn ==4.25 !
colnames(a) <- c('Hybrid_ID','a','aper','adyn')
a_ToBind <- subset(a,select=c('Hybrid_ID','adyn'))

# Response to Tnight Mimi :
# GNij = GN15i + B (Tnightij - 15)
# So the response is : B (Tnightij - 15) :
# The response is also normalized for a RI0 for which there is "no" efect on GN of RADi :
# GNmax we want would be : grain.nb15.seq + [(300-ri.zero.seq)*resp.ri.sequ]
M.Tnight <- read.table(paste0(PATHS[3],'/RespTnight.csv'), 
                       sep=",",                                           
                       header=TRUE,dec=".",
                       row.names=NULL,
                       as.is=FALSE, 
                       fill=TRUE)

M.Tnight$GN15cor <-  M.Tnight$grain.nb15.seq + ((300-M.Tnight$ri.zero.seq)*M.Tnight$resp.ri.sequ)
RespTnight <- select(M.Tnight,Hybrid_ID,resp.tnight.seq,GN15cor,resp.ri.sequ)
colnames(RespTnight) <- c('Hybrid_ID','resp.tnight','GN15cor','resp.Ri')

## Merge to get the last set of parameters :
Params <- merge(NfMax,OutPhenoPanel,all.x=T)[c(1,4:8)]
Params_v2 <- merge(Params,SensiPSI,all.x=T)
#Params_v3 <- merge(Params_v2,GenGroups)
Params_v4 <- merge(Params_v2,a_ToBind,all.x=T,all.y=F)
Params_v4$c <- Params_v4$c_norm*Params_v4$adyn*10
Params_v5 <- merge(Params_v4,GrainParams_toBind,all.x=T,all.y=F)
Params_v6 <- merge(Params_v5,RespTnight,by='Hybrid_ID',all.x=T)

## Maintenant on rajoute SensiRAD : 
SensiRAD <- read.table(paste0(PATHS[3],'/SensiRAD.csv'), 
                       sep=",",                                           
                       header=TRUE,dec=".",
                       row.names=NULL,
                       as.is=FALSE, 
                       fill=TRUE)
colnames(SensiRAD) <- c('Hybrid_ID','Leaf','SensiRAD','Oo')
SensiRAD <- select(SensiRAD,Hybrid_ID,SensiRAD)
Params_v7 <- merge(Params_v6,SensiRAD,by='Hybrid_ID',all.x=T)


## Final set of parameters : ----
Final_Params_Panel <- subset(Params_v7, select=c('Hybrid_ID','Nf','Phyllo',
                                                 'Ntip_emerg','Ligules','Tht_1st_Lig',
                                                 'c','adyn','potKernelWt','GNmax','GrainNo','resp.tnight','GN15cor','resp.Ri','SensiRAD'))

colnames(Final_Params_Panel) <- c('Hybrid_ID','Nfinal',
                                  'Phyllochron','Leaf_tip_emerg',
                                  'First_LiguloChrone','tt_firstligule_sinceapp',
                                  'c','a',
                                  'potKernelWt','GNmax',
                                  'GrainNoMax','Res.Tnight',
                                  'resp.Ri','GN15cor','SensiRAD')

write.csv(Final_Params_Panel,file ='D:/Work/Sims - DROPS/4.Out/Parameters_DROPS_Panel_ForAPSIM_WithNA.csv')

# Calculate b sensi to VPD of LER : 
Final_Params_Panel$b <- (-0.0758*Final_Params_Panel$c)-0.3672

Params_WithoutNAs <- Final_Params_Panel[rowSums(is.na(Final_Params_Panel))<1, ]
Params_WithoutNAs$b <- (-0.0758*Params_WithoutNAs$c)-0.3672

# Final CSV :
write.csv(Params_WithoutNAs,file ='D:/Work/Sims - DROPS/4.Out/Parameters_DROPS_Panel_ForAPSIM.csv')