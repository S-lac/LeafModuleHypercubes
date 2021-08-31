Leaf_dev_stages <- function(Nfinal, # Max number of leaves
                            Leaf_tip_emerg, # Number of tip at emerg
                            atip, # Phyllo
                            tt_firstligule_sinceapp, # thermal time of first lig since appear (30-50)
                            First_LiguloChrone) # Ligulochrone 
{
  
  
  ## Centered at emergence ! 
  btip = Leaf_tip_emerg  
  Lagmax = 5.4
  leaf_no_init_at_emerg = 4.5
  LIR = 2*atip
  k_bl = 1.412
  Nlim = 6.617
  alpha_tr = 0.5
  k_ll = 2.4
  tt_to_emerg = 0
  abl = k_bl * atip;
  tt_lim_tip = ( (Nlim - btip) / atip )   + tt_to_emerg;
  bbl = Nlim - (abl * tt_lim_tip);
  tt_ll1=tt_firstligule_sinceapp;
  a_ll1 = First_LiguloChrone
  a_ll2=k_ll*a_ll1;
  b_ll1=1-(a_ll1*tt_ll1);
  N_limll = alpha_tr * Nfinal ;
  b_ll2 = (((N_limll- b_ll1) * (a_ll1 - a_ll2))/a_ll1)+b_ll1;
  
  ## Initiation of variables: 
  Dev_Stages <- data.frame("Leaf"=NA,"tt_Ini"=NA,"tt_BegElong"=NA,"tt_Tip"=NA,"tt_EndElong"=NA,"tt_Lig"=NA)
  Leaves <- seq(1,Nfinal, 1)
  tt_ini <- rep(NA,Nfinal)
  tt_tip <- rep(NA,Nfinal)
  tt_startExp <- rep(NA,Nfinal)
  tt_endExp <- rep(NA,Nfinal)
  tt_lig <- rep(NA,Nfinal)
  
  for (leaf in Leaves)
  {
    
    ######################################
    
    ## Initiation : 
    if (leaf < leaf_no_init_at_emerg)
      {
      if (leaf < floor(leaf_no_init_at_emerg))
        {
        tt_ini[leaf] <- 0
        }
      else 
        {
        tt_ini[leaf] <- (1 - (leaf_no_init_at_emerg - floor(leaf_no_init_at_emerg))) * 1/LIR
        }
      }
    else 
      {
        tt_ini[leaf] <- tt_ini[leaf-1] + 1/LIR
      }
      
    ######################################
    
    ## Tip appearance : 
    if (leaf <=Leaf_tip_emerg)
    {
      tt_tip[leaf] <- 0
    }
    else 
    {
      tt_tip[leaf] <- ((leaf)-btip)/atip + tt_to_emerg
    }

    
    ######################################
    
    ## Begignning of elongation : 
    tt_bl = ((leaf)-bbl)/abl;
    if (tt_bl>=tt_tip[leaf])
    {
      tt_startExp[leaf]= tt_tip[leaf];
    }
    else 
    { 
      tt_startExp[leaf]= tt_bl + tt_to_emerg;
    }
    
    ######################################
    
    ## Ligulation : 
    if (leaf <= N_limll)
    {
      tt_lig[leaf]= ( (leaf-b_ll1) / a_ll1 ) + tt_to_emerg ;
    }
    else 
    {
      tt_lig[leaf]= ( (leaf-b_ll2) / a_ll2 ) + tt_to_emerg  ;
    }
    
    ######################################
    
    lag = Lagmax * leaf;
    ## End of elongation : 
    if (leaf <= N_limll)
    {
      tt_endExp[leaf]= ( (leaf-b_ll1) / a_ll1 - lag ) + tt_to_emerg  ;
    }
    else if (leaf <= Nfinal-1) {
      tt_endExp[leaf]= ( (leaf-b_ll2) / a_ll2 - lag ) + tt_to_emerg  ;
    }
    else 
    {
      tt_endExp[leaf]=  tt_endExp[leaf-1];
		}
    
  }
  
  Dev_Stages <- as.data.frame(cbind(Leaves,tt_ini,tt_tip,tt_startExp,tt_endExp,tt_lig))
  colnames(Dev_Stages) <- c("Leaf","nInitLeaves","nLeafTips","nLeafExp","nLeafFullyExp","nLigules")
  
  return(Dev_Stages)
}


WidthAtLowRad <- function( Nfinal ) # Max number of leaves 
{
  
  Parameters <- data.frame(W_min = 1.5 ,
                           L_rank = 0.2,
                           Rank_Max = 0.66,
                           W_max = 7.8,
                           Width_Curve = 0.29,
                           SlopeScale = 0.25)
  
  DataSim <- data.frame(Leaf = seq(1,nff,1),
                        SimWidth = rep(NA,nff))
  
  for (leaf in seq(1,nff,1))
  {
    
    SemiPlateau <- (round(Parameters$Width_Curve*nff)/3)
    
    # First non mature leaves : 
    if (leaf <= round(Parameters$L_rank*nff))
    {
      
      DataSim$SimWidth[leaf] <- Parameters$W_min
      
      # Upward :
    } else if (leaf <= round(Parameters$Rank_Max*nff)-SemiPlateau){
      
      Num <- (Parameters$W_max -Parameters$W_min)
      Den <- round(Parameters$Rank_Max*nff)-SemiPlateau-round(Parameters$L_rank*nff)
      UpSlope <- Num / Den
      
      DataSim$SimWidth[leaf] <- Parameters$W_min + (UpSlope * (leaf-round(Parameters$L_rank*nff)))
      
      # Plateau :
    } else if (leaf <= round(Parameters$Rank_Max*nff)+SemiPlateau-1){
      
      DataSim$SimWidth[leaf] <- Parameters$W_max
      
      # Downward :
    } else {
      
      LimSup <- round(Parameters$Rank_Max*nff)+SemiPlateau-1
      DownSlope <- UpSlope*Parameters$SlopeScale
      DataSim$SimWidth[leaf] <- DataSim$SimWidth[leaf-1] - (DownSlope * (leaf-LimSup))
      
    }
    
  }
  
  return(DataSim)
}


LERmax <- function(Nfinal) # Max number of leaves 
{
  
  Amax  <- 0.94;
  SIGMA <- 0.35 * Nfinal;
  ALPHA  <-  0.59 * Nfinal;
  DOWN  <-  (2*SIGMA^2);
  
  A6 = 1 / exp( (-(6-ALPHA)^2) / DOWN) ;
  
  DataSim <- data.frame(Leaf = seq(1,Nfinal,1),
                        LERmax = rep(NA,nff))  
  
  for (leaf in seq(1,nff,1))
  {
    
    DataSim$Leaf[leaf] <- leaf
    
    UP = (-((leaf+1)-ALPHA)^2) ;
    DataSim$LERmax[leaf] <- A6 * exp(UP / DOWN);
    

  
  }
  
  return(DataSim)
  
  
}
