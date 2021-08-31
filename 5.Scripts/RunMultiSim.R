
###POUR LES SCENARIOS WELL WATERED ET RAINFED
source('D:/Work/Sims  - DROPS GxE/5.Scripts/MultiAPSIMfunctions.r')
library(XML)
library(data.table)


CommandLine = args[1]
# If you send it from the command line : 
#args <- commandArgs(TRUE)

if (CommandLine=='True')
{
  Clim<-args[2] #'ClimF1'        ##A remplir
  Sce<- args[3] #'Rainfed'  ## A remplir
  gcm <- as.numeric(args[4]) #1 ## A remplir
} else 
{
  Clim<-'ClimA'        ## A remplir
  Sce<-'Well_Watered'  ## A remplir
  gcm <-'NoG1' ## A remplir
}

Directory.Apsim<-"D:/Work/Sims  - DROPS GxE/1.APSIM/Model"
Directory.param<-"D:/Work/Sims  - DROPS GxE/3.Data"
Directory.Out<-"D:/Work/Sims  - DROPS GxE/4.Out"
PATH.APSIM <- "D:/APSIM/Apsim76/Model/"
  
comb<-paste(Clim,Sce,sep='_')
Directory.Outputs<-paste(Directory.Out,Clim,Sce,sep='/')

#Valeur param pour changement dur?e du cycle
setwd(Directory.param)
param<-read.table('valeurs_param_cycle.csv', sep=',',header=TRUE,dec=".",row.names=NULL, 
                  as.is=FALSE)

# Optimums : 
load('D:/Work/Sims  - DROPS GxE/3.Data/ClimA_Optimum.RData')
Optimums <- subset(ClimA_Optimum,Scenario=='Well_Watered')

# Valeur params pour l'azote : 
Nrates <- read.table('azote.opt climA.csv', sep=',',header=TRUE,dec=".",row.names=NULL, 
                     as.is=FALSE)
Nrates$NO3N <- Nrates$Nopt *(16*3+28)/(28)

## Load and treat parameters for simulation :
Parameters <- read.table('D:/Work/Sims  - DROPS GxE/3.Data/Parameters_DROPS_Panel_ForAPSIM.csv', 
                         sep=",", header=TRUE,dec=".",row.names=NULL,as.is=FALSE, fill=TRUE)

Parameters$respTnight <- -Parameters$Res.Tnight
MyParams <- subset(Parameters,select=colnames(Parameters)[c(2,7:11,17,18)])

ParamsB73 <- subset(MyParams,Hybrid_ID=='B73_H')
Scale_a <- 3.5/ParamsB73$a
Scale_c <- 3.5/ParamsB73$c
MyParams$a <- MyParams$a*Scale_a
MyParams$c <- MyParams$c*Scale_c
MyParams$b <- (-0.0758*MyParams$c)-0.3672


#S?lection des fichiers APSIM :
Apsim_Files<-list.files(path=Directory.Apsim, pattern=comb)

GCM <- gcm  
Simulation_Name<- paste(GCM,Clim,Sce,sep='_')


Stopped_At = 1  # Put the genotype on which you finished simulating (Nf)! 
WorkOnErrors <- F

# 1/ Loop over genotypes : ----
Start <- Sys.time()

for (i in Stopped_At:length(unique(MyParams$Hybrid_ID)))
{
  
  Genotypes <- unique(MyParams$Hybrid_ID)
  geno <- as.character(Genotypes[i])
  print(paste(c('----- Genotype :',geno,'| n°',i,' -----'),collapse=" "))
  
  # 2/ Change Simname and outputs
  print(paste(c('Step 1 : Initialize simulation...'),collapse=" "))
  ParamsToUse <- subset(MyParams,Hybrid_ID==as.character(geno))
  SimulationName<-paste(ParamsToUse[,'Hybrid_ID'])
  NewApsimName<-paste(c(SimulationName,"apsim"),collapse=".")
  ApsimFile<-paste(Directory.Apsim,Simulation_Name,sep='/')
  ApsimFile<-paste(ApsimFile,'apsim',sep='.')
  
  #On charge le fichier .apsim avec la simulation
  SimFileTree<-xmlTreeParse(ApsimFile,getDTD = F , useInternalNodes = T)
  
  # CHeck if we have values for the genotype :
  if (nrow(ParamsToUse)>0)
  {
    # Loop over genotypic parameters and change them :
    for (var in colnames(ParamsToUse[2:ncol(ParamsToUse)]))
    {
      VarNode<-getNodeSet(SimFileTree, paste(c("//",var),collapse="") )
      NewValue<- ParamsToUse[,var]
      for (node in 1:length(VarNode))
      {
        ChangeXmlValue(VarNode[node],NewValue) 
      }
    }
    
    #On va chercher les noeuds des paramètres initiaux utiles dans les calculs de la nouvelle valeur des paramètres
    Nf_init_Node<-getNodeSet(SimFileTree, paste(c("//","Nfinal"),collapse="") )
    Nfemerg_init_Node<-getNodeSet(SimFileTree, paste(c("//","Leaf_tip_emerg"),collapse="") )
    Phyllo_Node<-getNodeSet(SimFileTree, "//Phyllochron" )
    ttei_init_Node<-getNodeSet(SimFileTree, paste(c("//","tt_endjuv_to_init"),collapse="") )
    
    # List of SimNames :
    #SimulationNameNode<-getNodeSet(SimFileTree, "//simulation")
    #xattrs <- xpathSApply(SimFileTree, "//*/simulation/@name")
    #Sims <- data.frame("Title"=xattrs)
    
    # Put the optimum genotype for each simulation :
    for (site in 1:nrow(Optimums))
    {
      Nf_init_Value<-as.numeric(xmlValue(Nf_init_Node[[site]]))
      Nfemerg_init_Value<-as.numeric(xmlValue(Nfemerg_init_Node[[site]]))
      Phyllo_Value<-as.numeric(xmlValue(Phyllo_Node[[site]]))
      ttei_init_Value<-as.numeric(xmlValue(ttei_init_Node[[site]])) #On récupère les valeurs initiales
      
      ttef_init<-(Nf_init_Value-Nfemerg_init_Value)/Phyllo_Value #On calcule le temps thermique entre émergence et flag initial
      Nf_new<-trunc(((ttef_init+param[Optimums$Genotype[site]-9,]$Value)*Phyllo_Value)+Nfemerg_init_Value) #On rajoute un tt en fonction de la durée du cycle voulue puis on calcule le NF final correspondant
      ChangeXmlValue(Nf_init_Node[site],Nf_new) #On change la valeur du NFfinal dans le xml
      
      ttei_new<-ttei_init_Value*((ttef_init+param[Optimums$Genotype[site]-9,]$Value)/ttef_init)#Idem pour le temps thermique entre endjuv et init
      ChangeXmlValue(ttei_init_Node[site],ttei_new)
    } 
    
    
    # On traite l'azote : 
    fert_amount_sow_Node<-getNodeSet(SimFileTree, paste(c("//","fert_amount_sow"),collapse="") )
    
    for (node in 1:length(fert_amount_sow_Node))
    {
      ChangeXmlValue(fert_amount_sow_Node[node],round(Nrates$Nopt[node],digits=0))
    }  
    
    print(paste('==> Done !'),collapse=" ")
    
    # Save the new apsimfile in Model Directory :
    print(paste(c('Step 2 : Saving simulation file...'),collapse=" "))
    saveXML(xmlRoot(SimFileTree), 
            file= paste0(PATH.APSIM,NewApsimName), 
            prefix= '<?xml version="1.0" encoding="utf-8"?>',indent= T)
    print(paste('==> Done !'),collapse=" ")
    
    # Simulate :
    print(paste(c('Step 3 : Simulation...'),collapse=" "))
    setwd(PATH.APSIM)
    
    command<-paste("apsim",NewApsimName)
    
    system(command)
    
    
    # Get/Save outputs  : 
    print(paste(c('Step 4 : Loading outputs...'),collapse=" "))
    
    Output_Sim <- loadApsim(PATH.APSIM, 
                            loadAll = TRUE, 
                            ext = ".out", 
                            returnFrame = TRUE, 
                            n = 0, fill = TRUE, 
                            addConstants = TRUE)
    
    HarvestDataOutputs <- Get_Harvest_Data(Output_Sim,list(Output_Sim$Title,Output_Sim$year))
    rm(Output_Sim)
    
    HarvestDataOutputs$Genotype<-geno
    
    if (i ==1)
    {Global_Output <- HarvestDataOutputs
    }else{
      Global_Output <- rbind(Global_Output,HarvestDataOutputs)
    }
    rm(HarvestDataOutputs)
  }
  
  print(paste('==> Done !'),collapse=" ")  
  
}

Stop <- Sys.time()
TimeDiff <- difftime(Stop, Start, units="secs")
MyTimeDiff <- format(.POSIXct(TimeDiff,tz="GMT"), "%H:%M:%S")
print(paste(c(''),collapse=" "))
print(paste(c('----- Simulation loops finished -----'),collapse=" "))
print(paste(c(''),collapse=" "))
print(paste(c('Time of execution is '),MyTimeDiff,collapse=" "))
