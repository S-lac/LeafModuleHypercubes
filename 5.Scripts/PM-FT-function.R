##########                                                                       ##########
##########          ET0 PENMAN-MONTEITH Tardieu et al. 2015                      ##########
##########                                                                       ##########
#___Santiago_Alvarez Prado_24/09/2015__________________________________________________________________________


#Agregar o sacar LA de acuerdo a lo que quiera estimar
ET0_ft_func <- function(manip,vpd,zalt,zvent,lat,J,Tmoy,Rs,hmoy,PPFD,t){
#__________________________________________________________________________
#########   1. Param?tres equation  ###################################
#  

#Jout= flujo transpiratorio (mm3 s-1 plant-1)
s=0.1849*Tmoy^2+1.005*Tmoy+50.87    	#pente de la courbe T (?C)-P(T) Pa ?C-1
s <-as.numeric(s)
ro=0.00002*Tmoy^2-0.0047*Tmoy+1.2921  #densit? air sec Kg m-3
gama = 0.0643*Tmoy+64.9  		       	  #Constante psychrom?trique Pa ?K-1
lambda=(-0.0024*Tmoy+2.5008)*1000000 	#Chaleur latente vaporisation J Kg-1
conversion = -0.1424*Tmoy+43.917		  #convertir gs de mol m-2 s-1 en m s-1
vpd1 = vpd * 1000                    #deficit de presion de vapor en Pa
#Constantes. Jugar con valores de gs y ga
Cp = 1012				                  	  #chaleur massique air J Kg-1?K-1


#### Rayonnement
phi <- lat*(pi/180) # lattitude (radian)
dr <- 1+0.033*cos(((2*pi)/(365))*J) #  inverse relative distance T-S
deltaMin <- 0.409*sin((((2*pi)/(365))*J)-1.39) #solar declination
b <- (2*pi*(J-81))/364
Sc <- 0.1645*sin(2*b)-0.1255*cos(b)-0.025*sin(b)
omega <- (pi/12)*((t + 0.06667*(0-3.8772300)+Sc)-12) # sunset hour angle per hour
w1 <- omega - ((pi*1)/24) 
#w1 <- ifelse(w1<0, 0, w1)
w2 <- omega + ((pi*1)/24)
#w2 <- ifelse(w2<0, 0, w2)
Ra <- ((12*60)/pi)*0.082*dr*((w2 -w1)*sin(phi)*sin(deltaMin) + (cos(phi)*cos(deltaMin)*(sin(w2)-sin(w1)))) # Hourly extraterrestrial radiation (MJ.m-2.hour-1)
Ra <- ifelse(Ra<0, 0, Ra)
Rso <- (0.75+(2*10^-5)*zalt)*Ra  # clear sky solar radiation(MJ.m-2.hour-1). 
albedo <- 0.23 # Grass = 0.20-0.25 segun FAO. Maize= 0.10-0.23 / 0.23 en estado vegetativo  
Rns <- (1-albedo)*Rs # net solar or shortwave radiation (MJ.m-2.hour-1)
Rns <- ifelse((is.na(Rns)==TRUE),0,Rns)
e0Tmoy <-	0.6108*exp((17.27*Tmoy)/(Tmoy+237.3))	# saturation vapour pressure Tmax
ea <-	(e0Tmoy*(hmoy/100)) # actual vapor pressure [Pa]
#RR <-ifelse(Rso==0,0.8, Rs/Rso)
RR <-Rs/Rso
RR1 <- ifelse(RR<0, 0, RR)
RR2 <- ifelse(RR1>1, 1, RR1)
Rnl <- ((2.043*10^-10)*((Tmoy+273.16)^4)*(0.34-(0.14*sqrt(ea)))*((1.35*(RR2))-0.35)) #  net longwave radiation (MJ.m-2.hour-1). Divido por 24 para llevarlo a hora
#Rn <- (Rns-Rnl)# net radiation (MJ.m-2.hour-1)
Rn <- (Rns-Rnl)*277.7778# net radiation (W.m-2)
#Rn <-as.numeric(Rn)
Rn <- ifelse((is.na(Rn)==TRUE),0,Rn)
#G <- ifelse(Ra==0,0.5*Rn1, 0.1*Rn1)
#Rn = Rn1-G
# Ejemplo: Rn=66.7245
#__________________________________________________________________________
#########   3. ET0 PM-FT 2015   ######################


ga= 0.016 # mean of all manip for 2011
# For the others :
# ZB13=0.012
# ZB12=0.012
# ZA13=0.010
# ZA16=0.03 
if(manip=="ZB13"){ga=0.012}
if(manip=="ZB12"){ga=0.012}  
if(manip=="ZA13"){ga=0.010}  
if(manip=="ZA16"){ga=0.030}  

gsmin <- 0.05 # De acuerdo con Allen et al., 2006 
gsmax <- 0.4  #  Esta en moles. De acuerdo con Allen et al., 2006 
gsd <- ifelse(PPFD>500, gsmax,gsmax*(PPFD/500))# original= 500/PPFD

#gsd <- as.numeric(gsd)
gs1 <- ifelse(PPFD==0, gsmin, gsd)
#Ejemplo: gs1=0.314  

gs <- gs1/conversion
#L1 <- L*0.2

A1 <- s*Rn + ro*Cp*ga*vpd1
A2 <- lambda*(s + gama*(1 + ga/gs))
Jout <- (A1/A2)*1000000#*L1 # mg/m? sec

return(Jout)
}