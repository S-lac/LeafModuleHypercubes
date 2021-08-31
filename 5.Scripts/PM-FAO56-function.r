##########                                                                       ##########
##########          CALCUL DE L'ET0 PENMAN-MONTEITH D'APRES LA FAO               ##########
##########                                                                       ##########
#___Emilie_MILLET_05/06/2013__________________________________________________________________________


#__________________________________________________________________________
#########   1. PARAMETRES DE LA LOCALISATION   ############################
#zalt <-  altitude du site (m)
#zvent <- hauteur de mesure du vent  (m)
#lat <- latitude du site en décimal si possible
#__________________________________________________________________________
#########   2. VARIABLES JOURNALIERES   ###################################
#J <- jour julien
#Tmin en °C
#Tmax  en °C
#Tmoy en °C
#Rhmax en %
#Rhmin en %
#Rs solar radiation (MJ.m-2.day-1)
#uz  windspeed (m.s-1)

ET0_func <- function(zalt,zvent,lat,J,Tmax,Tmin,Tmoy,Rhmax,Rhmin,Rs,uz){
  #__________________________________________________________________________
  #########   3. CALCUL DES PARAMETRES DE L'EQUATION   ######################
  #### 3.1. VPD
  P <- 101.3*((293-0.0065*zalt)/293)^5.26 # pression atm (kPa)
  gama <- ((1.013*10^-3)*P)/(0.622*2.45) # psychrometric constant (kPa.°C-1)
  deltaMaj <- (4098*(0.6108*exp((17.27*Tmoy)/(Tmoy+237.3))))/(Tmoy+237.3)^2 # slope of saturation pressure curve (kPa.°C-1)
  e0Tmin <- 0.6108*exp((17.27*Tmin)/(Tmin+237.3))	 #	saturation vapour pressure Tmin
  e0Tmax <-	0.6108*exp((17.27*Tmax)/(Tmax+237.3))	#saturation vapour pressure Tmax
  ea <-	((e0Tmin*(Rhmax/100))+(e0Tmax*(Rhmin/100)))/2	# actual vapour pressure
  vpd <- mean(c(e0Tmin,e0Tmax),na.rm=T)-ea		 # VPD
  #don$vpd[don$Jour==J] <- vpd
  #### 3.2. Rayonnement
  phi <- lat*(pi/180) # lattitude (radian)
  dr <- 1+0.33*cos(((2*pi)/(365))*J) #	inverse relative distance T-S
  deltaMin <- 0.409*sin((((2*pi)/(365))*J)-1.39) #solar declination
  omegas <- acos((-tan(phi))*(tan(deltaMin))) # sunset hour angle
  Ra <- ((24*60)/pi)*0.082*dr*((omegas*sin(phi)*sin(deltaMin))+(cos(phi)*cos(deltaMin)*sin(omegas)))  # extraterrestrial radiation (MJ.m-2.day-1)
  Rso <- (0.75+(2*10^-5)*zalt)*Ra  # clear sky solar radiation(MJ.m-2.day-1)
  albedo <- 0.23 # albedo of a hypothetic reference grass
  Rns <- (1-albedo)*Rs # net solar or shortwave radiation (MJ.m-2.day-1)
  Rnl <- (4.903*10^-9)*(((Tmin+273.16)^4+(Tmax+273.16)^4)/2)*(0.34-(0.14*sqrt(ea)))*((1.35*(Rs/Rso))-0.35) #  net longwave radiation (MJ.m-2.day-1)
  Rn <- Rns-Rnl # net radiation (MJ.m-2.day-1)
  #### 3.3. Soil heat flux
  G <- 0 #soil heat flux (MJ m-2 day-1)  DAILY
  #G <- 0.1*Rn #soil heat flux (MJ m-2 day-1)  HOURLY
  #### 3.4. Vent à 2m
  u2 <- uz*(4.87/(log((67.8*zvent)-5.42))) # vitesse à 2m
  #__________________________________________________________________________
  #########   4. ET0 PM d'apres FAO56   ######################
  A1 <- 0.408*deltaMaj*(Rn-G)
  A2 <- gama*(900/(Tmoy+273))*u2*vpd
  A3 <- deltaMaj+(gama*(1+0.34*u2))
  ET0 <- (A1+A2)/A3
  
  return(print(ET0))
}
