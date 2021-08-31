ThermalTime <- function(
  
##############################################################       

## IN :


  Method,  
  
# Method : (1) 'Linear' or (2) 'PT' or (3) 'Apsim'

  DayData,TempData,  

# DayData : DAYS (w/e format - doesn't really matter)
#           JUST ONE VALUE PER DAY !
# TempData : TEMPERATURE (° Celsius)
#           SAME LENGTH AS DayData ! 
  
  Params  )
  
# Params : data.frame with function parameters to apply 
#          see each thermal time function below for parameters

################################################################ 

## OUT :


# data.frame with daily values of : 
  
  # $cumtht = cumulated thermal time from timestep before start (0) to end
  
  # Can be:
  # (1)calculated linear thermal time based Tbase or 
  # (2)calculated thermal time as Parent & Tardieu ,New Phytologist, 2012 or 
  # (3)calculated thermal time based APSIM crop model (see def above)  )
  
  # $tht = thermal time value of each timestep requested
  
  # $MeanTemp = Mean temperature for each day 
  

##############################################################       

## THERMAL TIME CALCULATION FUNCTIONS

#  -----------------------------------------------
## 1.Linear :
#  ----------

# Params$Tbase

# thermaltime=(T - Tbase ) / TimeScale

#  -----------------------------------------------
## 2.Tardieu-Parent function (days at 20 degree) :
#  -----------------------------------------------

# Params$DH ,Params$alpha ,Params$Topt ,Params$R

# thermaltime=(T*(1+exp(-DH/(R*293)) (alpha * (1 - (293 /Topt)))) / (293 * exp (-DH / (293 * R))) * exp((-DH./(T.*R))) ) / (1+(exp((-DH./(ToUse.*R))) ^((1-(T./Topt)).*alpha)

# example for a thermal time tardieu-parent:
# Topt= 33.4 °C,  
# Tbase= 20 °C, 
# alpha = 30.8 °C, 
# DH = 73.9 kj/mol.

#  ------------------------------------------------
## 3.APSIM function : APSIM/Model ...from Maize.xml
#  ------------------------------------------------

# No need of Params, any data named Params will do ! 

# <x_temp units="oC" description="cardinal temps.">0  18.0 26.0 34.0 44.0</x_temp>
# <y_tt units="oC" description="thermal time">0  10.0 18.0 26.0  0.0</y_tt>
# ttParams.read(  scienceAPI, "x_temp",   "y_tt");

# Parameters :
# X=temp° | Y =tt 
# ==> for X = 0 : 18, Y =(10/18)*x
# ==> for X = 18: 34, Y =((34-18) / (26-10))*x+(-8)
# ==> for X = 34:44, Y = (44-34) / (0-26) * x +(440/26)

# Apsim uses a three hour degree days simulation by using a
# simulated 3 hour air temp between Tmax and Tmin.



############################################################## 
############################ CODE ############################
############################################################## 

{
  
  #-------------------------  
  #     INITILISATION
  #-------------------------
  
  # Default parameters :
  if(missing(Params))
  {
    Params<-data.frame(DH=73.9,Topt=33.4,Tbase=20,alpha=30.8)
  }

  # Initialise the data : 
  ThT <- vector()
  cumTht <- vector()
  Tmean <- vector()
  ThermaltimeDiscrete <- vector()

  #------------------------- 
  #   DAILY CALCULATION
  #------------------------- 
  
  # Loop over days : 
  ID <- 1  
  for (day in unique(DayData))
  {

    # Data of today :
    Treat <- TempData[DayData==day]

    
    # Scaling from input data (Quarter/Hour/Day) :
    Scaling <- length(Treat)  
    
    # Calculate each timestep : 
    if (Method=='Linear'){    
      
      # Calculate linear thermal time : 
      ThermaltimeDiscrete<-(Treat-Params$Tbase)/ Scaling
      
      # Get to zeor if temp is lower than tbase (no negative) :
      ThermaltimeDiscrete[ThermaltimeDiscrete<0]<- 0    
      
    } else if (Method=='PT'){
      
      DH <- Params$DH*1000
      alpha <- Params$alpha
      Topt <- Params$Topt+273 # Change in Kelvin
      R <- 8.314 
      
      TempDatainK <- Treat+273
      
      A <- ( 1 + exp (-DH / (R*293)) ^ (alpha * (1 - (293 / Topt)))) / (293 * exp (-DH / (293 * R)))
      
      exponentielle<-exp((-DH/(TempDatainK*R)))
      Numerateur<-TempDatainK*A*exponentielle
      puissance<-(1-(TempDatainK/Topt))*alpha
      Denominateur<-1+(exponentielle^puissance)
      
      ThermaltimeDiscrete<- (Numerateur/Denominateur)/ Scaling
    
    } else if (Method=='Apsim'){
      
      Nt <- 1
      
      for (Temp in Treat)
      {
        
        if (Temp <= 0)
        
        {
          
          ThermaltimeDiscrete[Nt]<-0
          
        } else if (Temp > 0 && Temp <= 18){
        
          a<-(10/18)
          b<-0
          ThermaltimeDiscrete[Nt]<-(a*Temp+b)/ Scaling
        
        } else if (Temp > 18 && Temp <= 34){
          
          a<-((34-18) / (26-10))
          b<-(-8)
          ThermaltimeDiscrete[Nt]<-(a*Temp+b)/ Scaling
        
        } else if (Temp > 34 && Temp <= 44){
        
          a<-((44-34) / (0-26))
          b<-(440/26)
          ThermaltimeDiscrete[Nt]<-(a*Temp+b) / Scaling
        
        } else {
          
          ThermaltimeDiscrete[Nt]<-0
          
        }
        
        Nt <- Nt+1
        
      }
      
    }
    
    # Put it in the vector : 
    if (ID == 1) {
      
      cumTht[ID] <- sum(ThermaltimeDiscrete,na.rm=T)
      
    } else {
      
      cumTht[ID] <- cumTht[ID-1]+sum(ThermaltimeDiscrete,na.rm=T)
      
    }
    
    ThT[ID]<- sum(ThermaltimeDiscrete,na.rm=T)
    Tmean[ID] <- mean(Treat,na.rm=T)  
    
    # Switch day : 
    ID <- ID+1
  }

  #------------------------- 
  #         OUTPUTS
  #-------------------------   
  
  ThermalTime <- data.frame(cumtht = cumTht,
                            tht = ThT,
                            MeanTemp = Tmean)
  
  return(ThermalTime)
  
}

