#
#   Suncourse, astronomical function
#
#
#location parameters (montpellier)
#
lat <- 43.61 * pi / 180
long <- 3.87 * pi / 180
#
#
# conversion hUTC vers heure solaire et reverse (ratp, grevet)
#
eqntime <- function(doy) {
  om = 0.017202 * (doy - 3.244)
  teta = om + 0.03344 * sin(om) * (1 + 0.021 * cos(om)) - 1.3526
  tphi = 0.91747 * sin(teta) / cos(teta)
  dphi = atan(tphi) - om + 1.3526
  dphi = sapply(dphi,function(x) {
    if ((x + 1) <= 0)
      x <- ((x + 1 + 1000 * pi) %% pi) -1
    x
  })
  dphi * 229.2 / 60# Equation of time (in hour) : discrepency between "sun at zenith" and "noon UTC" at Greeenwich
}
#
UTC2sun <- function(hUTC,doy) (hUTC + long / 15 - eqntime(doy)) %% 24
sun2UTC <- function(hsun,doy) (hsun - long / 15 + eqntime(doy)) %% 24
#
#declinaison du soleil
decsun <- function(doy) {
  alpha = 2 * pi * (doy - 1) / 365
  0.006918 - 0.399912 * cos(alpha) + 0.070257 * sin(alpha)
}
#
angh <- function(hsun) 2 * pi / 24 * (hsun - 12)
# elevation
thetasun <- function(doy,hsun) {
  dec <- decsun(doy)
  ah <- angh(hsun)
  costheta <- sin(lat) * sin(dec) + cos(lat) * cos(dec) * cos(ah)
  asin(costheta)
}
# angle azimutal  (sens positif Nord vers west entre nord et sun)
phisun <- function(doy,hsun) {
  dec <- decsun(doy)
  ah <- angh(hsun)
  D = sin(lat) * cos(ah) - cos(lat) * sin(dec) / cos(dec)
  pi + ifelse(ah < 0, -1, 1) * ifelse(D == 0, pi / 2, abs(atan(sin(ah) / D)))
} 
# Intensite relative
Irel <- function(doy,hsun) {
  dec <- decsun(doy)
  ah <- angh(hsun)
  costheta <- sin(lat) * sin(dec) + cos(lat) * cos(dec) * cos(ah)
  costheta
}
# suncourse
suncourse <- function(doy, dt=60) {
  h <- seq(dt/2,24*60,dt) / 60
  theta <- thetasun(doy,h)
  phi <- phisun(doy,h)
  I <- Irel(doy,h)
  sel <- theta > 0
  data.frame(hsun=h[sel], hUTC = sun2UTC(h[sel],doy),elev = theta[sel] / pi * 180, azim = phi[sel] / pi * 180, I = I[sel] / sum(I[sel]))
}
#duree du jour pour tous les jours de l'annee
djour <- sapply(1:365,function(x) nrow(suncourse(x,1))/60)
#
# Spitters estimation of RdRs (diffus/global).Rg is global irradiance (W.m-2), hUTC
#
spitters <- function(hUTC,Rg,doy) {
  hsun <- UTC2sun(hUTC,doy)
  costheta <- Irel(doy,hsun)
  Io <- 1370 * (1 + 0.033 * cos(2 * pi * (doy - 4) / 366))#eclairement (w/m2) a la limitte de l'atmosphere dans un plan perpendiculaire aux rayons du soleil, fonction du jour
  So <- Io * costheta
  RsRso <- Rg / So
  R <- 0.847 - 1.61 * costheta + 1.04 * costheta * costheta
  K <- (1.47 - R) / 1.66
  #
  res <- R
  res[RsRso <= 0.22] <- 1
  sel <- RsRso <= 0.35 & RsRso > 0.22 
  res[sel] <- 1 - 6.4 * (RsRso[sel] - 0.22)^2
  sel <- RsRso <= K & RsRso > 0.35
  res[sel] <- 1.47 - 1.66 * RsRso[sel]
  res
}
  

