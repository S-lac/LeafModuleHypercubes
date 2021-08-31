
## Changing a value in an xms Parsed Tree :
ChangeXmlValue <- function(XNode,NewValue)
{
  sapply(XNode,function(G){
    text = paste(NewValue)
    xmlValue(G) = text
  })
}

## Getting all data for harvest dates only :
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
# Param: ext - file extension filter e.g. "out" for APSIM output files.
LoadAll <- function(dir, ext)
{
  wd <- dir #set working directory
  setwd(wd)
  
  files <- list.files(wd, paste(ext, "$", sep="")) # create a list of files
  
  DF <- NULL
  for (f in files) {
    fileName <- strsplit(f, ".", fixed = TRUE)
    print(fileName[[1]][1]) #print file information to console
    err <- tryCatch(file.title <- names(read.table(f,skip=2,header=T, nrows=1)), error = function(e) print(e)) #read the titles
    
    if(!inherits(err, "error")){
      file <- read.table(f, skip=4, ,header=F, col.names=file.title, na.strings = "?", stringsAsFactors=FALSE) # read the data
      file$FileName <- fileName[[1]][1]
      
      DF <- rbind(DF, file) # bind the file to the main data frame
    }
    else
      print(paste("Could not read file: ", fileName[[1]][1]))
  }
  names(DF)   <- c(file.title, "File")
  #DF$Date <- as.Date(DF$Date, format="%d/%m/%Y")
  return (DF)
}

# Similar to LoadAll; can be used when column names are different
LoadDifferent <- function(dir, ext)
{
  wd <- dir #set working directory
  setwd(wd)
  
  files <- list.files(wd, paste(ext, "$", sep="")) # create a list of files
  if (length(files == 0))
    return
  
  # store each output file in a seperate DF; put all the frames in a list.
  data <- list()
  
  for (f in files) {
    fileName <- strsplit(f, ".", fixed = TRUE)
    print(fileName[[1]][1]) #print file information to console
    err <- tryCatch(file.title <- names(read.table(f,skip=2,header=T, nrows=1)), error = function(e) print(e)) #read the titles
    
    if(!inherits(err, "error")){
      file <- read.table(f, skip=4, ,header=F, col.names=file.title, na.strings = "?", stringsAsFactors=FALSE) # read the data
      file$FileName <- fileName[[1]][1]
      
      names(file)   <- c(file.title, "File")
      if("dd.mmm.yyyy" %in% colnames(file)){
        names(file)[names(file) == "dd.mmm.yyyy"] <- "Date"
      }
      file$Date <- as.Date(file$Date, format="%d/%b/%Y")
      
      data[[fileName[[1]][1]]] <- file # add the DF to the main list for this directory
    }
    else
      print(paste("Could not read file: ", fileName[[1]][1]))
  }
  
  return (data)
}

# My loading function :
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

Fatality <- function(filename)
{
  #Get all tasks running :
  Tasks <- system('tasklist',show.output.on.console=T,intern=T)
  # Kill jobrunner :
  while (sum(grepl("JobRunner.exe", Tasks))>0)
  {
    print('Killing JobRunner...')
    system('taskkill /IM JobRunner.exe /F',
           show.output.on.console=F,intern=T)
    Tasks <- system('tasklist',show.output.on.console=T,intern=T)
  }
  # Kill all Apsim.exe :
  while (sum(grepl("Apsim.exe", Tasks))>0)
  {
    print('Killing APSIM...')
    system('taskkill /IM Apsim.exe /F',
           show.output.on.console=F,intern=T)
    Tasks <- system('tasklist',show.output.on.console=T,intern=T)
  }
  print('NO FATALITY !')
} 
