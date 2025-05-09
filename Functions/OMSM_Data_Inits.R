#### OMSM_Datalist ####
## This function generates the list containing all the data the OMSM will need to run
OMSM_Datalist <- function(
    Data.sources, 
    Data.zoops, 
    Sources, 
    Tracers, 
    TDF_meta, 
    TDF_proto
) {
  N_T <- length(Tracers$all) # number of tracers
  N_S <- length(Sources) # number of organic matter sources
  N_Z <- nrow(Data.zoops) # number of zooplankton samples
  sd.anl <- 0.52
  sdTDF_scaling_factor <- 1
  sdY_scaling_factor <- 1
  
  Data.1 <- list(
    ## How many tracers will be used in this model?
    N_T = N_T,
    ## How many zooplankton samples are there?
    N_Z = N_Z,
    ## Where should the model look for zooplankton amino acid isotope data?
    Y   = as.matrix(Data.zoops[Tracers$all]),
    ## and the analytical uncertainty in those measurements
    sd_Y = as.matrix(Data.zoops[SDTracers$all])*sdY_scaling_factor,
    ## Which tracers will be used in the FWL equation and not see the data?
    i_T_FWL = which(Tracers$all %in% Tracers$FWL),
    ## Which tracers will be used in the MTS equation and not see the data?
    i_T_MTS = which(Tracers$all %in% Tracers$MTS),
    ## Which tracers will be used in the mixing model and see the data?
    i_T_mix = which(Tracers$all %in% Tracers$mix),
    ## Which tracers do not fractionate?
    i_T_non = which(Tracers$all %in% Tracers$non),
    ## Which tracers do fractionate?
    i_T_frac = which(Tracers$all %in% Tracers$frac),
    ## Which tracers have constant TDFs throughout the food web?
    i_T_const = which(Tracers$all %in% Tracers$constTDF),
    ## Which tracers have variable TDFs throughout the food web?
    i_T_var = which(Tracers$all %in% Tracers$varTDF),
    ## Where should the model look for TDF data and its associated uncertainty?
    # for metazoan TDFs
    TDF_m   = as.numeric(TDF_meta[Tracers$frac]),
    sdTDF_m = as.numeric(TDF_meta[SDTracers$frac])*sdTDF_scaling_factor,
    # for protozoan TDFs
    TDF_p   = as.numeric(TDF_proto[Tracers$frac]),
    sdTDF_p = as.numeric(TDF_proto[SDTracers$frac])*sdTDF_scaling_factor
  )
  
  ## Where should the model look for organic matter source data?
  ## making generalized names A-F for a max of six possible organic matter sources
  Sources.alpha <- c("X_A","X_B","X_C","X_D","X_E","X_F")
  SDSources.alpha <- c("sdX_A","sdX_B","sdX_C","sdX_D","sdX_E","sdX_F")
  nsams_alpha   <- c("N_A","N_B","N_C","N_D","N_E","N_F")
  ## Loop through all sources described in setup chunks above to add that data
  ## to the Data.1 list as individual matrices for each source
  for (i in 1:N_S) {
    sams_source <- Data.sources[Data.sources[["Source"]]==Sources[i],c(Tracers$all,SDTracers$all)]
    Data.1[[Sources.alpha[i]]] = as.matrix(sams_source[Tracers$all])
    Data.1[[SDSources.alpha[i]]] = as.double(colMeans(sams_source[SDTracers$all]))
    Data.1[[nsams_alpha[i]]]   = sum((Data.sources[["Source"]]==Sources[i])*1)
  }
  # Setting i_T_ vectors =0 if no tracers in that ccategory were included
  if(length(Data.1$i_T_non)==0){
    Data.1$i_T_non <- 0
  }
  if(length(Data.1$i_T_frac)==0){
    Data.1$i_T_frac <- 0
  }
  if(length(Data.1$i_T_const)==0){
    Data.1$i_T_const <- 0
  }
  if(length(Data.1$i_T_var)==0){
    Data.1$i_T_var <- 0
  }
  return(Data.1)
}




#### OMSM_Initlist ####
## Generates a list containing the initial conditions for the OMSM
OMSM_Initlist <- function(
    seed = 222,
    Nchains = 3,
    N_Z,
    N_S,
    N_T
) {
  # Define MCMC initial values. Edit Nchains or ranges if adjusting MCMC settings.
  set.seed(seed)  # For reproducibility
  Nchains <<- Nchains   # Number of chains
  inits_1 <- vector(mode = "list", length = Nchains)
  # names for organic matter sources
  means_alpha <- c("mean_A", "mean_B", "mean_C", "mean_D", "mean_E", "mean_F")
  # and their uncertainty
  sd_alpha <- c("sd_A", "sd_B", "sd_C", "sd_D", "sd_E", "sd_F")
  # This function will generate starting values for each chain and save them as a list
  fun_init_1 <- function(i) {
    inits <- list(
      pz = rdirichlet(n = N_Z, rep(1, N_S)),  # Mixing coefficients
      FWL = runif(n = N_Z, 1, 5),             # Food web length
      MTS = runif(n = N_Z, 0, 5),             # Metazoan trophic steps
      .RNG.seed = i + 1,
      .RNG.name = c("base::Super-Duper", "base::Wichmann-Hill", "base::Marsaglia-Multicarry")[i %% 3 + 1]
    )
    # We define starting values for sources in a separate loop so we can be flexible with the number of sources included in the model
    for (j in 1:N_S) {
      inits[[means_alpha[j]]] <- runif(N_T, -40, 20)
      inits[[sd_alpha[j]]] <- runif(N_T, 0, 2)
    }
    return(inits)
  }
  for (i in 1:Nchains) inits_1[[i]] <- fun_init_1(i)
  return(inits_1)
}