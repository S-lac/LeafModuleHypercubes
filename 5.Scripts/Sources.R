
## Code to load sources for multi usage !

# Libraries ----
library(akima)
library(lattice)
library(latticeExtra)
library(RColorBrewer)
library(ggplot2)
library(grid)
library(gtable)
library(reshape2)
library(hydroGOF)
library(XML)
library(APSIM)
library(data.table)
library(R.utils)

## Try the command and give NULL if after time limit : 
TimeOutHandling <- function()
{
  #Get all tasks running :
  Tasks <- system('tasklist',show.output.on.console=F,intern=T)
  # Kill all Apsim.exe :
  while (sum(grepl("Apsim.exe", Tasks))>0)
  {
    print('Killing APSIM...')
    system('taskkill /IM Apsim.exe /F',show.output.on.console=F,intern=T)       
    Tasks <- system('tasklist',show.output.on.console=F,intern=T)
  }
  #                       # Kill jobrunner :
  #                       while (sum(grepl("JobRunner.exe", Tasks))>0)
  #                       {  
  #                         print('Killing JobRunner...')
  #                         system('taskkill /IM JobRunner.exe /F',show.output.on.console=F,intern=T)
  #                         Tasks <- system('tasklist',show.output.on.console=F,intern=T)
  #                       }
}


## Getting all data for harvest dates only :
Get_Harvest_Data <- function(Output_Sim,By)
{
  HarvestDateSim <-   
    do.call("rbind",lapply(split(Output_Sim,By,drop=T),
                           function(x)
                           {
                             Temp<-subset(x, DaysAfterSowing==max(x$DaysAfterSowing))
                             return(Temp[1,])
                           }
    ) ## End lapply
    ) ## End do.call rbind  
  return(HarvestDateSim)
}

GetDataframeToSave <- function(RData,By,Colnames,Scenario) 
{
  
  load(RData)
  HarvestData <- Get_Harvest_Data(Global_Output,By)
  rm(Global_Output)
  HarvestData_Sub <- subset(HarvestData,select=Colnames)
  rm(HarvestData)
  HarvestData_Sub$Scenario <- Scenario
  return(HarvestData_Sub)
}

# Theme base ----
My_Theme <-   
  theme(axis.title.x = element_text(size=20, vjust=-0.35),
        axis.title.y = element_text(size=20, vjust=1.5,hjust = 0.5),
        axis.text.x=element_text(size=20, vjust=0.5),
        axis.text.y=element_text(size=20, vjust=0.5),
        panel.background = element_rect(fill = "white"),
        panel.grid.major= element_line(colour="white"),
        panel.grid.minor= element_line(colour="white"),
        plot.title = element_text(size=20, face="bold", vjust=2) ,
        panel.border = element_blank(), axis.line = element_line())

#Theme Carte ----
theme_map <- function(base_size = 8, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          panel.margin = unit(0, "lines"),
          plot.background = element_blank(),
          legend.position = "bottom",
          legend.direction="horizontal",
          legend.title=element_blank(),
          legend.text=element_text(size=5),
          legend.text.align=0.5)
}

stat_smooth_func <- function(mapping = NULL, data = NULL,
                             geom = "smooth", position = "identity",
                             ...,
                             method = "auto",
                             formula = y ~ x,
                             se = TRUE,
                             n = 80,
                             span = 0.75,
                             fullrange = FALSE,
                             level = 0.95,
                             method.args = list(),
                             na.rm = FALSE,
                             show.legend = NA,
                             inherit.aes = TRUE,
                             xpos = NULL,
                             ypos = NULL) {
  layer(
    data = data,
    mapping = mapping,
    stat = StatSmoothFunc,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      method = method,
      formula = formula,
      se = se,
      n = n,
      fullrange = fullrange,
      level = level,
      na.rm = na.rm,
      method.args = method.args,
      span = span,
      xpos = xpos,
      ypos = ypos,
      ...
    )
  )
}


StatSmoothFunc <- ggproto("StatSmooth", Stat,
                          
                          setup_params = function(data, params) {
                            # Figure out what type of smoothing to do: loess for small datasets,
                            # gam with a cubic regression basis for large data
                            # This is based on the size of the _largest_ group.
                            if (identical(params$method, "auto")) {
                              max_group <- max(table(data$group))
                              
                              if (max_group < 1000) {
                                params$method <- "loess"
                              } else {
                                params$method <- "gam"
                                params$formula <- y ~ s(x, bs = "cs")
                              }
                            }
                            if (identical(params$method, "gam")) {
                              params$method <- mgcv::gam
                            }
                            
                            params
                          },
                          
                          compute_group = function(data, scales, method = "auto", formula = y~x,
                                                   se = TRUE, n = 80, span = 0.75, fullrange = FALSE,
                                                   xseq = NULL, level = 0.95, method.args = list(),
                                                   na.rm = FALSE, xpos=NULL, ypos=NULL) {
                            if (length(unique(data$x)) < 2) {
                              # Not enough data to perform fit
                              return(data.frame())
                            }
                            
                            if (is.null(data$weight)) data$weight <- 1
                            
                            if (is.null(xseq)) {
                              if (is.integer(data$x)) {
                                if (fullrange) {
                                  xseq <- scales$x$dimension()
                                } else {
                                  xseq <- sort(unique(data$x))
                                }
                              } else {
                                if (fullrange) {
                                  range <- scales$x$dimension()
                                } else {
                                  range <- range(data$x, na.rm = TRUE)
                                }
                                xseq <- seq(range[1], range[2], length.out = n)
                              }
                            }
                            # Special case span because it's the most commonly used model argument
                            if (identical(method, "loess")) {
                              method.args$span <- span
                            }
                            
                            if (is.character(method)) method <- match.fun(method)
                            
                            base.args <- list(quote(formula), data = quote(data), weights = quote(weight))
                            model <- do.call(method, c(base.args, method.args))
                            
                            m = model
                            eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                                             list(a = format(coef(m)[1], digits = 3), 
                                                  b = format(coef(m)[2], digits = 3), 
                                                  r2 = format(summary(m)$r.squared, digits = 3)))
                            func_string = as.character(as.expression(eq))
                            
                            if(is.null(xpos)) xpos = min(data$x)*0.9
                            if(is.null(ypos)) ypos = max(data$y)*0.9
                            data.frame(x=xpos, y=ypos, label=func_string)
                            
                          },
                          
                          required_aes = c("x", "y")
)
