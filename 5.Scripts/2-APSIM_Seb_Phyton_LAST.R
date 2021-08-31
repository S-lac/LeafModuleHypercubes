
# Initial setup  : ----
rm(list = ls())
 
# Sources
source('D:/SAP-LEPSE/APSIM/5. Scripts/MultiAPSIMfunctions.r')

# R - librairies : 
library(dplyr)
library(reshape2)
library(data.table)
library(reticulate)
 

# Python - librairies : 
ApsimPython <- import_from_path('pyapsim', 
                                'D:/SAP-LEPSE/APSIM/Apsim76/Python', 
                                convert = FALSE)

######################################################

# Paths to be used as 'global' variables : ----

RootDirectory <- 'D:/SAP-LEPSE/APSIM/Apsim76/Model'



ParameterFile <- "D:/SAP-LEPSE/APSIM/Datos_DROPS_AMAIZING/Parameters_APSIM_May_2020.csv" 

###############
##- ANALYSIS

MyParamsToUse <- read.table(ParameterFile, sep = ",", header = T)

######################################################


# Get list of sims in rep : ----
setwd(RootDirectory)
Simfiles <- list.files(pattern = "\\.sim$") #??? Else put it yourself :

# Inputs of the new module in python : 
whichSim  <- Simfiles
whereSim  <-  'D:/SAP-LEPSE/APSIM/Apsim76/Model'
list_of_params_names  <-  colnames(MyParamsToUse)
whereTo <- 'D:/SAP-LEPSE/APSIM/Apsim76/Model'

# Test variables: ----
#geno <- 'B104_H'
#var <- colnames(MyParams)[2]
Stopped_At <- 1 #Aca va el numero de genotipo que iba corriendo al trabarse

# 1/ Loop over genotypes : ----
Start <- Sys.time()
for (i in Stopped_At:length(unique(MyParams$Hybrid_ID)))
{
  #i <- 1
  
  
  ######## Inputs of the new python module to run APSIM : 
  list_of_params_values  <- as.character(MyParamsToUse[i,])
  Genotype <- MyParams$Hybrid_ID[i]
  
  
  
  print(paste(c('----- Genotype :',as.character(Genotype),'| n°',i,' -----'),collapse = " "))
  
  # 2/ Change Simname and outputs
  print(paste(c('Step 1 : Write sims...'),collapse = " "))
  
  # Check if we have values for the genotype :
  if (length(list_of_params_values) > 0)
  {
    
    ApsimPython$read_and_write_sims(whichSim,whereSim,list_of_params_names,list_of_params_values,whereTo)
    print(paste('==> Done !'),collapse = " ")
    
    # Simulate :
    print(paste(c('Step 2 : Simulation...'),collapse=" "))
    setwd(whereTo)
    ApsimPython$run_apsim(whereTo,whichSim)
    
    
    # Get/Save outputs  : 
    print(paste(c('Step 3 : Loading outputs...'),collapse=" "))
    
    Output_Sim <- loadApsim(whereTo, 
                            loadAll = TRUE, 
                            ext = ".out", 
                            returnFrame = TRUE, 
                            n = 0, fill = TRUE, 
                            addConstants = TRUE)
    
    Output_Sim$Genotype <- Genotype
    
    if (exists('Global_Output')==FALSE)
    {Global_Output <- Output_Sim
    }else{
      Global_Output <- rbind(Global_Output,Output_Sim)
    }
    rm(Output_Sim)
  }
  
  print(paste('==> Done !'),collapse = " ")  
  
}

Stop <- Sys.time()
TimeDiff <- difftime(Stop, Start, units="secs")
MyTimeDiff <- format(.POSIXct(TimeDiff,tz="GMT"), "%H:%M:%S")
print(paste(c(''),collapse=" "))
print(paste(c('----- Simulation loops finished -----'),collapse=" "))
print(paste(c(''),collapse=" "))
print(paste(c('Time of execution is '),MyTimeDiff,collapse=" "))


# 2/ Save and treat outputs : ----
#save(file = 'D:/Documents/SAP-LEPSE/APSIM/3. Out/AMAIZING_ALL/Global_Output_Parameters_APSIM10_390_geno_AMAIZING_DROPS_TS_Satolas_Souprosse_DEN_Nieder_profden_20.02.2018.Rdata',Global_Output)
#save(file = 'D:/Documents/SAP-LEPSE/APSIM/3. Out/WW_AMAIZING/Global_Output_Parameters_APSIM10_245_geno_WD_ParamPheno_AMAIZING_15.02.2018.Rdata',Global_Output)
#save(file = 'D:/SAP-LEPSE/APSIM/3. Out/Global_Output_NEW_Resptnight_AMAIZING_DROPS_Ita_x_6_v16_27.01.2019.Rdata',Global_Output)
save(file = 'D:/SAP-LEPSE/APSIM/3. Out/Global_Output_NEW_Resptnight_AMAIZING_DROPS_HighCV_Camp_v8_22.05.2020.Rdata',Global_Output)
