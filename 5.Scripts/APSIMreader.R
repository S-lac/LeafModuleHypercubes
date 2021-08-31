############################
## THIS is the Function/app : 
## 

## This app is build to analyse adequation between data and simulation for the APSIM model
## The input is actually a .out file from apsim outputs 
## (the "loadApsim" function loads all .out files in )


##------------------------------------------------------------------------------

## USE THE LIBRARY : 
####################

library(ggplot2)
library(RColorBrewer)
library(shiny)
library(shinyBS)
library(XML)
library(DT)
library(rhandsontable)
library(reshape2)


# LoadApsim does everything without thinking ;: 
## The function used for loadApsim : 
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

ChangeXmlValue <- function(XNode,NewValue)
{
  sapply(XNode,function(G){
    text = paste(NewValue)
    xmlValue(G) = text
  })
}

##------------------------------------------------------------------------------

# THIS IS THE WEBPAGE STYLE ! 
##---------------------------

ui <- fluidPage(
  
  bsAlert("alert"),
  
  navbarPage("The Shiny APSIM",
             
  ################################
  ######### SIMULATION ###########
  ################################
  
    tabPanel("Home",icon=icon('home'),
             img(src='http://www.apsim.info/Validation/apsim_logo.png', height = 100, width = 200),
             h4("Welcome to the Shiny APSIM"),
             p("This application has been developed by Sebastien Lacube for the crop model APSIM."),
             p('Its purpose is to visualize APSIM files, Sim files and Output files,to construct plots and save data'),
             p('directly with a Shiny interface compatible with Rstudio.'),
             p('Have fun !')
             ),
  
    tabPanel("Simulation",icon=icon("terminal"),
             
      
      sidebarPanel(
        

        
        actionButton("ChooseSimFolder","  Choose folder",
                     icon=icon('folder-open'),
                     width = 200),
        uiOutput("choose_ApsimFile"),

        
        actionButton("LoadApsimFile","  Load Apsim file",
                     icon=icon('upload '),
                     style="color: #fff; background-color: #337ab7; border-color: #2e6da4",
                     width = 200),
        
        uiOutput("choose_sim"),

        # Define the checkbox to facet with sims : 
        checkboxInput(inputId='ApplytoAllSim',
                          label='Apply changes to all sims ?',
                          value=FALSE),

        textInput(inputId='APSIMFileName', 
                  'Enter a filename :', 
                  value = "NewAPSIMFile",
                  width = 200,
                  placeholder = NULL),
        
        actionButton("CreateNewApsimFile","  Create New APSIM file",
                     icon=icon('download'),
                     width = 200),
        
        actionButton("RunAPSIM","  Run APSIM",
                     icon=icon('play'),
                     width = 200)

        ), # End sidebarPanel
        

      mainPanel(
        
        tabsetPanel("TabSet_Params_Values",
                    
                    tabPanel("Constants",
                             icon=icon('tree'),
                             rHandsontableOutput('MaizeConstants')
                    ),
                    
                    tabPanel("Genotype",
                             icon=icon('leaf'),
                             rHandsontableOutput('MaizeGenotype')
                    ),
                    
                    tabPanel("Irrigation",
                             icon=icon('tint'),
                             rHandsontableOutput('MyIrrigation')
                    ),
                    
                    tabPanel("Surface OM",
                             icon=icon('tint'),
                             rHandsontableOutput('MySurfaceOM')
                    )   
                             
                    
        ) # End tabset panel
        
      )  # End main panel 
      
    ),

  ################################
  #########   OUTPUTS  ###########
  ################################
  
    tabPanel("Outputs",icon=icon("line-chart"),
      
      sidebarPanel(
        
        actionButton("ChooseOutputFolder","Choose folder",
                     icon=icon('folder-open'),
                     width = 200
                     ),
        
        actionButton("LoadOutputFiles","Load Output files",
                     icon=icon('upload '),
                     style="color: #fff; background-color: #337ab7; border-color: #2e6da4",
                     width = 200
                     ),
        

        tabsetPanel("TabSetOutputs",
          
          tabPanel("Sims",uiOutput("choose_simulation")
                   ),
          
          tabPanel("Vars",
                   uiOutput("choose_variable_x"),
                   uiOutput("choose_variable_y"),
                   uiOutput("Facet_Sims"),
                   p("EITHER..."),
                   uiOutput("OnlyThisYear"),
                   uiOutput("choose_year"),
                   p("OR..."),
                   uiOutput("Facet_Years")
                   ),
          tabPanel("Params",
                   uiOutput("choose_ymax"),
                   uiOutput("choose_xmax"),
                   uiOutput("choose_size"),
                   uiOutput("choose_geom_type")
                   )
          
          ) # End TabSet Panel
          
        ), # End sidebar panel
    
      mainPanel(
        
        tabsetPanel("TabSetPlotsDataSave",
          
          tabPanel("Plot",icon=icon("bar-chart-o"), plotOutput(outputId = "Apsim_Out",width = "100%", height = "500px")), 
          
          tabPanel("Data",icon=icon("database"), rHandsontableOutput('MyDataTable')), 
          
          tabPanel("Save",icon=icon("floppy-o"), 
                     textInput(inputId='PathSaveFile', 
                               'Enter a path to save output :', 
                               value = "D:/", 
                               width = '50%', 
                               placeholder = NULL),
                     textInput(inputId='FileName', 
                               'Enter a filename (no extension):', 
                               value = "NewOutputFile", 
                               width = '50%', 
                               placeholder = NULL),
                     actionButton("SaveFile","Save csv file",
                                  icon=icon('floppy-o '),
                                  style="color: #fff; background-color: #337ab7; border-color: #2e6da4",
                                  width = 200
                   )
                   
                   )
          
          ) # End tabset panel
        
        )  # End main panel 
    
    ) # End 1st NavBarPanel
  
  ) # End NavBarPage
  
) # End fluidpage
##------------------------------------------------------------------------------

# THIS IS THE SERVER RUNNING THE CODE TO CHANGE INPUT INTO OUTPUTS ! 
##-------------------------------------------------------------------

server <- function(input,output,session) 
{             
  
  #############################################################################
  ############################### SIMULATOR ###################################
  #############################################################################
  
  pathSims <- NULL
  
  # Define Folder to read APSIM file from :
  observe({
    
    if(input$ChooseSimFolder > 0){

      pathSims <<- choose.dir(default = 'D:/Work/Widths/1.APSIM',
                     caption = "Select folder from which to import '.APSIM' :")
      pathSims <<- gsub("\\\\", "/", pathSims)
      setwd(pathSims)
      files <<- list.files(pattern = "\\.apsim$")  

    }
    if(!is.null(pathSims)){

      # Define the listbox to choose the sim file : 
      output$choose_ApsimFile <- renderUI(
        {
          selectInput(inputId="choose_ApsimFile", 
                        label = "Which APSIM file to load :", 
                        files,
                        selected = 1, 
                        multiple = FALSE, 
                        selectize = FALSE, 
                        width = 200, 
                        size = 1) 
          
        }) # End Render choose_variable
        

    }
  })
  
  # When there is an APSIM file read it and show constants : 
  observe({
    
    if(input$LoadApsimFile > 0){
      
      SimFile <<- input$choose_ApsimFile
      SimFileTree<<-xmlTreeParse(paste(pathSims, SimFile, sep="/"), useInternalNodes = T)
      SimulationNameNode<-getNodeSet(SimFileTree, "//simulation")
      xattrs <- xpathSApply(SimFileTree, "//*/simulation/@name")
      Df.Simnames <<- data.frame("simulation"=xattrs)

      output$choose_sim <- renderUI(
        {
          selectInput(inputId="choose_sim", 
                      label = "Which SIM to show :", 
                      as.character(Df.Simnames$simulation),
                      selected = NULL, 
                      multiple = FALSE, 
                      selectize = FALSE, 
                      width = 200, 
                      size = 1) 
          
        }) # End Render choose_simulation
    }
    
  
  })
  
  
  #Df.Constants <- reactiveValues()
  #Df.Genotype <- reactiveValues()
  #Df.Constants$data <- data.frame(Parameter = character(0), Value = character(0))
  #Df.Genotype$data <- data.frame(Parameter = character(0), Value = character(0))
  
  observe({
    
      if (!is.null(input$choose_sim))
        
        {
        
          ID_SimName <- as.numeric(which(Df.Simnames==input$choose_sim))
          TopNode <- xmlRoot(SimFileTree)
          
          # Find tag with simulation : 
          TagNames <- xmlElementsByTagName(TopNode, 'simulation', recursive = FALSE)
          TagNames <- TagNames[[ID_SimName]]
          TagNames <- xmlElementsByTagName(TagNames, 'area', recursive = FALSE)
          TagNames_Irri <- xmlElementsByTagName(TagNames[[1]], 'irrigation', recursive = TRUE)
          TagNames_Out <- xmlElementsByTagName(TagNames[[1]], 'outputfile', recursive = FALSE)
          TagNames_VarOut <- xmlElementsByTagName(TagNames_Out[[1]], 'variables', recursive = FALSE)
          #TagNames_SurfaceOM <- xmlElementsByTagName(TagNames_Out[[1]], 'surfaceom', recursive = FALSE)
          
          TagNames <- xmlElementsByTagName(TagNames[[1]], 'maize', recursive = FALSE)
          TagNames_Csts <- xmlElementsByTagName(TagNames[[1]], 'maizeConstants', recursive = FALSE)
          TagNames_Geno <- xmlElementsByTagName(TagNames[[1]], 'maizeGenotype', recursive = FALSE)

          ## Retrieve maize Constants and genotypic parameters : 
          Constants <<- melt(getChildrenStrings(TagNames_Csts[[1]],asVector=F,addNames=T))
          Genotype <<- melt(getChildrenStrings(TagNames_Geno[[1]],asVector=F,addNames=T))
          Irrigation <<- melt(getChildrenStrings(TagNames_Irri[[1]],asVector=F,addNames=T))
          OutputFile <<- melt(getChildrenStrings(TagNames_VarOut[[1]],asVector=F,addNames=T))
          #SurfaceOM <<- melt(getChildrenStrings(TagNames_SurfaceOM[[1]],asVector=F,addNames=T))
          
          if (length(colnames(Constants))>0 && length(colnames(Genotype))>0 &&
              length(colnames(Irrigation))>0 && length(colnames(OutputFile))>0 
              #&& length(colnames(SurfaceOM)) 
              )
            
          {          
            
            colnames(Constants) <<- c("Value","Parameter")
            colnames(Genotype) <<- c("Value","Parameter")
            colnames(Irrigation) <<- c("Value","Parameter")
            colnames(OutputFile) <<- c("Variable","Status")
            #colnames(SurfaceOM) <<- c("Variable","Status")
            
            Df.Constants <<- reactiveValues(data=Constants[,c(2,1)])
            Df.Genotype <<- reactiveValues(data=Genotype[,c(2,1)])
            Df.Irrigation <<- reactiveValues(data=Irrigation[,c(2,1)])
            Df.OutputFile <<- reactiveValues(data=OutputFile[,c(2,1)])
            #Df.SurfaceOM <<- reactiveValues(data=OutputFile[,c(2,1)])
            
            output$MaizeConstants <- renderRHandsontable({
                            
              rhandsontable(Df.Constants$data,height=600,width=600) %>%
                hot_col(col='Parameter',halign='htCenter',strict=FALSE,copyable=T,readOnly=F,source=NULL,type = "autocomplete") %>%
                hot_cols(col=colnames(Df.Constants$data),colWidths=250,manualColumnMove=F,manualColumnResize=T) %>%
                hot_context_menu(allowRowEdit=F,allowColEdit=F) %>%
                hot_col(col='Value',strict=FALSE,copyable=T,readOnly=F,source=NULL,type = "autocomplete") %>%
                hot_table(highlightRow = TRUE)
                          
            })
                          
            output$MaizeGenotype <- renderRHandsontable({
                            
              rhandsontable(Df.Genotype$data,height=600,width=600) %>%
                hot_col(col='Parameter',halign='htCenter',strict=FALSE,copyable=T,readOnly=F,source=NULL,type = "autocomplete") %>%
                hot_cols(col=colnames(Df.Genotype$data),colWidths=250,manualColumnMove=F,manualColumnResize=T) %>%
                hot_context_menu(allowRowEdit=F,allowColEdit=F) %>%
                hot_col(col='Value',strict=FALSE,copyable=T,readOnly=F,source=NULL,type = "autocomplete") %>%
                hot_table(highlightRow = TRUE)
              
            })
            
            output$MyIrrigation <- renderRHandsontable({
              
              rhandsontable(Df.Irrigation$data,height=600,width=600) %>%
                hot_col(col='Parameter',halign='htCenter',strict=FALSE,copyable=T,readOnly=F,source=NULL,type = "autocomplete") %>%
                hot_cols(col=colnames(Df.Irrigation$data),colWidths=250,manualColumnMove=F,manualColumnResize=T) %>%
                hot_context_menu(allowRowEdit=F,allowColEdit=F) %>%
                hot_col(col='Value',strict=FALSE,copyable=T,readOnly=F,source=NULL,type = "autocomplete") %>%
                hot_table(highlightRow = TRUE)
              
            })

            output$MyOutputsVars <- renderRHandsontable({
              
              rhandsontable(Df.OutputFile$data,height=600,width=600) %>%
                hot_col(col='Parameter',halign='htCenter',strict=FALSE,copyable=T,readOnly=F,source=NULL,type = "autocomplete") %>%
                hot_cols(col=colnames(Df.OutputFile$data),colWidths=250,manualColumnMove=F,manualColumnResize=T) %>%
                hot_context_menu(allowRowEdit=F,allowColEdit=F) %>%
                hot_col(col='Value',strict=FALSE,copyable=T,readOnly=F,source=NULL,type = "autocomplete") %>%
                hot_table(highlightRow = TRUE)
              
            })
            
            #output$MySurfaceOM <- renderRHandsontable({
              
            #  rhandsontable(Df.OutputFile$data,height=600,width=600) %>%
            #    hot_col(col='Parameter',halign='htCenter',strict=FALSE,copyable=T,readOnly=F,source=NULL,type = "autocomplete") %>%
            #    hot_cols(col=colnames(Df.OutputFile$data),colWidths=250,manualColumnMove=F,manualColumnResize=T) %>%
            #    hot_context_menu(allowRowEdit=F,allowColEdit=F) %>%
            #    hot_col(col='Value',strict=FALSE,copyable=T,readOnly=F,source=NULL,type = "autocomplete") %>%
            #    hot_table(highlightRow = TRUE)
              
            #})
            
          }
          else {
            
            createAlert(session, "alert", "LoadingAlert", title = "Oops...error loading maize parameters...",
                        content = paste0("No maize constants/genotype found for ",input$choose_sim,'. Check if nodes are linked together in APSIM file !'), append = FALSE)  
          }
                
        }
  })          
   

  observe({
    
    if ( !is.null(input$MaizeConstants) 
         && !is.null(input$MaizeGenotype) 
         && !is.null(input$MyIrrigation)
         && !is.null(input$MyOutputsVars)
         && !is.null(input$MySurfaceOM)
         )
      
    {
      Df.Constants$data <<- hot_to_r(input$MaizeConstants,c('Parameter','Value'))
      Df.Genotype$data <<- hot_to_r(input$MaizeGenotype,c('Parameter','Value'))
      Df.Irrigation$data <<- hot_to_r(input$MyIrrigation,c('Parameter','Value'))
      Df.Irrigation$data <<- hot_to_r(input$MyOutputsVars,c('Variable','Status'))
      #Df.SurfaceOM$data <<- hot_to_r(input$MySurfaceOM,c('Variable','Status'))
      
    }
    
  })
  
  
  observe({
    
    if ( !is.null(input$MaizeConstants) 
         && !is.null(input$MaizeGenotype) 
         && !is.null(input$MyIrrigation)
         #&& !is.null(input$MyOutputsVars)
         )
      
    {
      
    }
    
  })
  
  
            
  
  
  #############################################################################
  ############################# OUTPUT READER #################################
  #############################################################################  
  
  pathOutputs <- NULL  
  Sim_Names <- NULL

  # Define folder to read outputs files from :
  observe({
    if(input$ChooseOutputFolder > 0){
      pathOutputs <<- choose.dir(default = 'D:/Work/Widths/1.APSIM',
                                 caption = "Select folder from which to import '.OUT' :")
      pathOutputs <<- gsub("\\\\", "/", pathOutputs)
    }
  })
  
  observe({
    if(input$LoadOutputFiles > 0){
      
      setwd(pathOutputs)
      files <- list.files(pattern = "\\.out$")
      
      if (length(files)!=0)
        {
        
        try(
          # Read APSIM outputs in the directory : 
          Apsim_Outputs <<- loadApsim(pathOutputs,
                                      loadAll= TRUE,
                                      ext='.out',
                                      returnFrame=TRUE,
                                      n=0,
                                      fill=TRUE,
                                      addConstants = TRUE)
          )
        
        if (!is.null(Apsim_Outputs))
        
        
          {
          Sim_Names <<- unique(Apsim_Outputs$Title)
          
        }
        
        else {
          createAlert(session, "alert", "LoadingAlert", title = "Oops...error loading output files...",
                    content = "There is an error with a '.out' file in the selected folder !", append = FALSE)
        }
         
        }
      else {
        
        createAlert(session, "alert", "LoadingAlert", title = "Oops...error loading output files...",
                    content = "No '.out' files from APSIM in the selected folder !", append = FALSE)
        }
      }
    if (!is.null(Sim_Names))
    {
      ## Read STUFF : 
      ## ------------
      
      # Define the variable to choose depending on the dataframe read : 
      df.Apsim_Outputs <- Apsim_Outputs[, !names(Apsim_Outputs) 
                                        %in% c("Title", 
                                               "fileName",
                                               "ApsimVersion",
                                               "Date",
                                               "stagename")]    
      
      # Define the checkboxes to choose the simulations : 
      output$choose_simulation <- renderUI(
        {
          checkboxGroupInput(inputId='choose_simulation',
                             label='Simulations to show on graphs',
                             choices=Sim_Names,
                             selected=NULL,
                             inline=FALSE)
        }) # End Render choose_simulation
    

      
      # Define the listbox to choose the variable  y: 
      output$choose_variable_y <- renderUI(
        {
          selectInput(inputId="choose_variable_y", 
                      label = "Choose the y variable :", 
                      colnames(df.Apsim_Outputs),
                      selected = "yield", 
                      multiple = FALSE, 
                      selectize = FALSE, 
                      width = '80%', 
                      size = 1) 
        }) # End Render choose_variable
  
      # Define the listbox to choose the variable x: 
      output$choose_variable_x <- renderUI(
        {
          selectInput(inputId="choose_variable_x", 
                      label = "Choose the x variable :", 
                      colnames(df.Apsim_Outputs),
                      selected = "TT", 
                      multiple = FALSE, 
                      selectize = FALSE, 
                      width = '80%', 
                      size = 1) 
        }) # End Render choose_variable
      
      # Change the maximum value of the plot y : 
      output$choose_ymax <- renderUI(
        {
          
          Var_name_y<- input$choose_variable_y
                
          if (!is.null(Var_name_y))
            {
            minimum <- floor(min(df.Apsim_Outputs[,names(df.Apsim_Outputs) %in% Var_name_y],na.rm=TRUE ))
            maximum <- ceiling(1.3 *max(df.Apsim_Outputs[,names(df.Apsim_Outputs) %in% Var_name_y],na.rm=TRUE ))
              
            # Get a slider that will help you choose the zoom on y ! 
            sliderInput(inputId = "choose_ymax",
                        label = "Change y axis :",
                        min   = minimum, 
                        max   = maximum,
                        value = c(minimum,maximum),
                        step = (floor(maximum-minimum))/100
                        )
            }
        } # End Render Yaxis
      )
      # Same goes for X axis ! ,
      output$choose_xmax <- renderUI(
        {
          
          Var_name_x<- input$choose_variable_x
          
          if (!is.null(Var_name_x))
          {
            minimum <- floor(min(df.Apsim_Outputs[,names(df.Apsim_Outputs) %in% Var_name_x],na.rm=TRUE ))
            maximum <- ceiling(1.3 *max(df.Apsim_Outputs[,names(df.Apsim_Outputs) %in% Var_name_x],na.rm=TRUE ))
            
            # Get a slider that will help you choose the zoom on x ! 
            sliderInput(inputId = "choose_xmax",
                        label = "Change x axis :",
                        min   = minimum, 
                        max   = maximum,
                        value = c(minimum,maximum),
                        step = 1
            )
          }
        } # End Render Xaxis
      ) 
      
      # Change the geom type size : 
      output$choose_size <- renderUI(
        {
          selectInput(inputId="choose_size", 
                      label = "Choose the geom size :", 
                      seq(0.5,6,0.5),
                      selected = "2.5", 
                      multiple = FALSE, 
                      selectize = FALSE, 
                      width = '80%', 
                      size = 1) 
        }
      )
      
      # Change the line/points size : 
      output$choose_geom_type <- renderUI(
        {
          selectInput(inputId="choose_geom_type", 
                      label = "Choose the geom type :", 
                      c('Lines','Points'),
                      selected = "Points", 
                      multiple = FALSE, 
                      selectize = FALSE, 
                      width = '80%', 
                      size = 1) 
        }
      )
      
      # Define the checkbox to facet with sims : 
      output$OnlyThisYear <- renderUI(
        {
          checkboxInput(inputId='OnlyThisYear',
                        label='Choose a specific year :',
                        value=FALSE)
        }) # End Render Facet_Sims
      
      
      # Choose a year if one sim is used : 
      output$choose_year <- renderUI(
        {
          ToCheck <-Apsim_Outputs[Apsim_Outputs$Title %in% input$choose_simulation,]
          selectInput(inputId="choose_year", 
                      label = "which one ?", 
                      unique(ToCheck$year),
                      selected = "", 
                      multiple = FALSE, 
                      selectize = FALSE, 
                      width = '80%', 
                      size = 1) 
        }
      )    
      # Define the checkbox to facet with sims : 
      output$Facet_Sims <- renderUI(
        {
          checkboxInput(inputId='Facet_Sims',
                             label='Plot by simulation',
                             value=FALSE)
        }) # End Render Facet_Sims
      
      # Define the checkbox to facet with sims : 
      output$Facet_Years <- renderUI(
        {
          checkboxInput(inputId='Facet_Years',
                        label='Plot by year',
                        value=FALSE)
        }) # End Render Facet_Sims
      
      # Begin RenderDataTable
      output$MyDataTable <- renderRHandsontable({
        
        df <- Apsim_Outputs[Apsim_Outputs$Title %in% input$choose_simulation,]
        df <- subset(df,select=c(input$choose_variable_x,input$choose_variable_y,'Title' ))
        rhandsontable(df
                      ,options=list(),
                      height=400,
                      width=600) %>%
        hot_cols(col=colnames(df),
                 colWidths=130,
                 manualColumnMove=F,
                 manualColumnResize=T)
          
        })
      # End render DataTable
      
      # Now we begin the code to generate the plot !
      output$Apsim_Out <- renderPlot(
        {
          
    #       Sims <- c(" No_Lag"," Current_Version")
    #       Var_name<-'LAI'
          Sims <<- input$choose_simulation
          Year <<- input$choose_year
          OnlyThisYear <- input$OnlyThisYear 
          x_Var <<- input$choose_variable_x 
          y_Var <<- input$choose_variable_y
          
          xlims <<- input$choose_xmax
          ylims <<- input$choose_ymax
          
          Geom_Type <<- input$choose_geom_type
          Size <<- input$choose_size
          Facet_Param <<- ceiling(sqrt(length(Sims)))
  
          if (!is.null(x_Var) && !is.null(y_Var) && !is.null(xlims) && !is.null(ylims) && length(Sims)>=1)
            {
              
              df.ToSave <<- Apsim_Outputs[Apsim_Outputs$Title %in% Sims,]
              df.ToSave <<- subset(df.ToSave,select=c(input$choose_variable_x,input$choose_variable_y,'Title' ))

            
              
              df.ToPlot <<- Apsim_Outputs[Apsim_Outputs$Title %in% Sims,]
              
              if(OnlyThisYear)
              {df.ToPlot <<- df.ToPlot[df.ToPlot$year %in% Year,]}
              
              p <- ggplot()
                
                          if (input$Facet_Sims && input$Facet_Years)
                          {
                          p <- p + facet_wrap(~ Title+year,ncol=Facet_Param,nrow=Facet_Param)
                          
                                  if (Geom_Type=="Lines") 
                                  {
                                    p <- p + geom_line(data=df.ToPlot,
                                             aes(x=get(x_Var),y=get(y_Var)),
                                             size=as.numeric(Size),
                                             col="blue")
                                  }
                                  else {
                                    p <- p + geom_point(data=df.ToPlot,
                                             aes(x=get(x_Var),y=get(y_Var)),
                                             size=as.numeric(Size),
                                             col="blue")
                                  }
                                  
                          }
                          else if (input$Facet_Sims) {
                            p <- p + facet_wrap(~ Title,ncol=Facet_Param,nrow=Facet_Param)
                            
                            if (Geom_Type=="Lines") 
                            {
                              p <- p + geom_line(data=df.ToPlot,
                                                 aes(x=get(x_Var),y=get(y_Var)),
                                                 size=as.numeric(Size),
                                                 col="blue")
                            }
                            else {
                              p <- p + geom_point(data=df.ToPlot,
                                                  aes(x=get(x_Var),y=get(y_Var)),
                                                  size=as.numeric(Size),
                                                  col="blue")
                            }
                          } 
                          else if (input$Facet_Years) {
                            p <- p + facet_wrap(~ year,ncol=Facet_Param,nrow=Facet_Param)
                            
                            if (Geom_Type=="Lines") 
                            {
                              p <- p + geom_line(data=df.ToPlot,
                                                 aes(x=get(x_Var),y=get(y_Var)),
                                                 size=as.numeric(Size),
                                                 col="blue")
                            }
                            else {
                              p <- p + geom_point(data=df.ToPlot,
                                                  aes(x=get(x_Var),y=get(y_Var)),
                                                  size=as.numeric(Size),
                                                  col="blue")
                            }
                          }   
                          else {
                                  if (Geom_Type=="Lines") 
                                  {
                                    p <- p + geom_line(data=df.ToPlot,
                                             aes(x=get(x_Var),y=get(y_Var),col=Title),
                                             size=as.numeric(Size))
                                  }
                                  else {
                                    p <- p + geom_point(data=df.ToPlot,
                                             aes(x=get(x_Var),y=get(y_Var),col=Title),
                                             size=as.numeric(Size))                         
                                  }
  
                          }
                          p <- p + labs(x=x_Var, y=y_Var)+
                                   scale_x_continuous(limits = xlims)+
                                   scale_y_continuous(limits = ylims)
              
              print(p)
          }
        }
      ) # End RenderPlot
    }
  })
  
  # Saving the file at the location you want :
  observe({
    if (input$SaveFile > 0 && !is.null(df.ToSave))
    {
      
      pathSaveCsv <<- input$PathSaveFile
      pathSaveCsv <<- gsub("\\\\", "/", pathSaveCsv)
      
      setwd(pathSaveCsv)
      write.table(df.ToSave,paste0(input$FileName,'.csv'),sep=";",row.names=F)      
      
    }
  })
  
    
}


## ----------------------------------------------------------------------



# THIS IS ACTUALLY LAUCHING YOUR APP ! 
##------------------------------------

shinyApp(ui=ui,server=server)



