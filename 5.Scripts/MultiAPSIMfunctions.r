######################################################################
######################################################################
######################################################################

## Libraries : 
library(beepr)
library(XML)

######################################################################
######################################################################
######################################################################

### FUNCTION : 

## Changing a value in an xml Parsed Tree :
# INPUT : 
# OUTPUT :
# EXAMPLE :
ChangeXmlValue <- function(XNode,NewValue)
{
  sapply(XNode,function(G){
    text = paste(NewValue)
    xmlValue(G) = text
  })
}

###################################

# Alarm when fnished multisim :
# INPUT : 
# OUTPUT :
# EXAMPLE :
AlarmFunction <- function(TimeInSec,TimeStep,Sound)
{
  Stop = TimeInSec/TimeStep
  for(i in 1:Stop){beep(sound=Sound)
                Sys.sleep(TimeStep)
               }
}

###################################

## Getting all data for harvest dates only :
# INPUT : 
# OUTPUT :
# EXAMPLE :
Get_Harvest_Data <- function(Output_Sim,By)
{
  HarvestDateSim <-   
    do.call("rbind",lapply(split(Output_Sim,By,drop=T),
                         function(x)
                         {
                           Temp<-subset(x, day==max(x$day))
                           return(Temp)
                         }
                       ) ## End lapply
                      ) ## End do.call rbind  
  return(HarvestDateSim)
}

# Param: dir - Directory to look for files, not recursive.
loadApsim <- function (dir, loadAll = TRUE, ext = ".out", returnFrame = TRUE, 
                       n = 0, fill = FALSE, addConstants = TRUE) 
{
  if (loadAll) {
    wd <- getwd()
    setwd(dir)
    files <- list.files(dir, paste(ext, "$", sep = ""))
  }
  else {
    files <- dir
  }
  if (length(files == 0)) 
    return
  allData <- list(NULL)
  fileCount <- 0
  for (f in files) {
    print(f)
    con <- file(f, open = "r")
    count <- 0
    size <- 1
    constants <- NULL
    namesFound <- FALSE
    unitsFound <- FALSE
    while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 
           0) {
      if (grepl("factors = ", oneLine)) {
        if (addConstants) {
          oneLine <- stringr::str_replace(oneLine, "factors = ", 
                                          "")
          split <- unlist(strsplit(oneLine, ";", fixed = "TRUE"), 
                          use.names = FALSE)
          for (s in split) {
            constants[length(constants) + 1] <- strsplit(s, 
                                                         "=", fixed = "TRUE")
          }
        }
      }
      else if (grepl("=", oneLine)) {
        if (addConstants) {
          constants[length(constants) + 1] <- strsplit(oneLine, 
                                                       "=", fixed = "TRUE")
        }
      }
      else {
        if (!namesFound) {
          colNames <- unlist(strsplit(stringr::str_trim(oneLine), 
                                      " ", fixed = TRUE), use.names = FALSE)
          colNames <- subset(colNames, colNames != "")
          namesFound <- TRUE
        }
        else if (!unitsFound) {
          units <- unlist(strsplit(stringr::str_trim(oneLine), 
                                   " ", fixed = TRUE), use.names = FALSE)
          units <- subset(units, units != "")
          if (length(units) != length(colNames)) 
            stop(paste("Error reading", f, "number of columns does match number of headings."))
          unitsFound <- TRUE
        }
        else {
          break
        }
      }
      count <- count + 1
    }
    close(con)
    data <- read.table(f, skip = count, header = FALSE, col.names = colNames, 
                       na.strings = "?", stringsAsFactors = FALSE)
    for (c in constants) {
      data[[ncol(data) + 1]] <- c[2]
      colNames[length(colNames) + 1] <- c[1]
    }
    colNames <- stringr::str_trim(colNames)
    names(data) <- colNames
    data$fileName <- f
    allData[[length(allData) + 1]] <- data
    fileCount <- fileCount + 1
    if (fileCount == n) 
      break
  }
  allData <- allData[!sapply(allData, is.null)]
  if (returnFrame) {
    allData <- as.data.frame(data.table::rbindlist(allData, 
                                                   fill = fill))
  }
  else {
    allData <- data.table::rbindlist(allData, fill = fill)
  }
  suppressWarnings(for (i in 1:ncol(allData)) {
    if (!any(is.na(as.numeric(allData[[i]])))) 
      allData[[i]] <- as.numeric(allData[[i]])
  })
  setwd(wd)
  return(allData)
}

###################################

# Adapt phenology to Nfinal (per apsim node):
# INPUT : 
# OUTPUT :
# EXAMPLE :
ChangePhenologyWithNfinal <- function(SimFileTree,param)
{
  
  #On va chercher les noeuds des paramètres initiaux utiles dans les calculs de la nouvelle valeur des paramètres
  Nf_init_Node<-getNodeSet(SimFileTree, paste(c("//","Nfinal"),collapse="") )
  Nfemerg_init_Node<-getNodeSet(SimFileTree, paste(c("//","Leaf_tip_emerg"),collapse="") )
  Phyllo_Node<-getNodeSet(SimFileTree, "//Phyllochron" )
  ttei_init_Node<-getNodeSet(SimFileTree, paste(c("//","tt_endjuv_to_init"),collapse="") )
  
  for (i in 1:length(Nf_init_Node)){
    Nf_init_Value<-as.numeric(xmlValue(Nf_init_Node[[i]]))
    Nfemerg_init_Value<-as.numeric(xmlValue(Nfemerg_init_Node[[i]]))
    Phyllo_Value<-as.numeric(xmlValue(Phyllo_Node[[i]]))
    ttei_init_Value<-as.numeric(xmlValue(ttei_init_Node[[i]])) #On récupère les valeurs initiales
    
    ttef_init<-(Nf_init_Value-Nfemerg_init_Value)/Phyllo_Value #On calcule le temps thermique entre émergence et flag initial
    Nf_new<-trunc(((ttef_init+param[simul,'Value'])*Phyllo_Value)+Nfemerg_init_Value) #On rajoute un tt en fonction de la durée du cycle voulue puis on calcule le NF final correspondant
    ChangeXmlValue(Nf_init_Node[i],Nf_new) #On change la valeur du NFfinal dans le xml
    
    ttei_new<-ttei_init_Value*((ttef_init+param[simul,'Value'])/ttef_init)#Idem pour le temps thermique entre endjuv et init
    ChangeXmlValue(ttei_init_Node[i],ttei_new)
  }
}

###################################

# Change Nrates (per apsim node):
# INPUT : 
# OUTPUT :
# EXAMPLE :
ChangeFertilisation <- function(SimFileTree,Nrates)
{
  fert_amount_sow_Node<-getNodeSet(SimFileTree, paste(c("//","fert_amount_sow"),collapse="") )
  
  for (i in 1:length(fert_amount_sow_Node))
  {
    ChangeXmlValue(fert_amount_sow_Node[i],round(Nrates[i],digits=0))
  }  
}

###################################

# Change irrigation (per apsim node) [carefull start and end irrigation !]
# INPUT : 
# OUTPUT :
# EXAMPLE :
ChangeIrrigation <- function(SimFileTree,start_irr,end_irr)
{
  # Debut Irr :
  StartIrrNode<-getNodeSet(SimFileTree, paste(c("//","start"),collapse="") )
  StartIrrValue<-list()
  for (j in 1:length(StartIrrNode)){
    StartIrrValue[j]<-as.character(start_irr[j,simul])
    ChangeXmlValue(StartIrrNode[j],StartIrrValue[j])
  }
  # Fin Irr :
  EndIrrNode<-getNodeSet(SimFileTree, paste(c("//","end"),collapse=""))
  EndIrrValue<-list()
  for (k in 1:length(EndIrrNode)){
    EndIrrValue[k]<-as.character(end_irr[k,simul])
    ChangeXmlValue(EndIrrNode[k],EndIrrValue[k])
  }
}


######## Use if problem with APSIM :

Fatality <- function(filename)
{
  #Get all tasks running :
  Tasks <- system('tasklist',show.output.on.console=T,intern=T)
  # Kill jobrunner :
  while (sum(grepl("JobRunner.exe", Tasks))>0)
  {
    print('...FATALITY ! Killing JobRunner...')
    system('taskkill /IM JobRunner.exe /F',
           show.output.on.console=F,intern=T)
    Tasks <- system('tasklist',show.output.on.console=T,intern=T)
  }
  # Kill all Apsim.exe :
  while (sum(grepl("Apsim.exe", Tasks))>0)
  {
    print('...FATALITY ! Killing APSIM...')
    system('taskkill /IM Apsim.exe /F',
           show.output.on.console=F,intern=T)
    Tasks <- system('tasklist',show.output.on.console=T,intern=T)
  }
  print('...NO MORE FATALITY !')
} 

