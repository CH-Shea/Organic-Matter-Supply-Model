# ## Simulating zooplankton data! ##
# use this function is d15NPhe is being treated as a fractionating tracer in the model
Sim_Zoop <- function(
    zoops.n = 21, # number of zooplankton to simulate
    incr = 0.2, # increment between zooplankton f values
    disperse = 0, # number of SDs by which to increase variance of zooplankton data beyond source means
    PTS, # a vector containing PTS values to simulate
    MTS, # a vector containing MTS values to simulate
    PTS.sd, # a vector containing SD of PTS values
    MTS.sd, # a vector containing SD of PTS values
    TDF_m, # trophic discrimination factors
    TDF_m.sd, # SD of trophic discrimination factors
    TDF_p, # trophic discrimination factors
    TDF_p.sd, # SD of trophic discrimination factors
    Sources, # names of organic matter sources
    Data.sources, # data for organic matter sources
    Tracers, # names of tracers
    SDTracers, # names of tracer SDs
    variables, # additional variables to include
    Use_Dirichlet = TRUE, # should the model usee Dirichlet to simulate zooplankton data?
    zoops.f, # if no Dirichlet, please supply matrix of f values
    seed = 222 # option to set seed
){
  
  # 
  # zoops.n = 50
  # PTS = c(0,1)
  # MTS = c(1,2)
  # PTS.sd = 0.5
  # MTS.sd = 0.5
  # TDF_m = TDF_meta[Tracers$frac] # trophic discrimination factors
  # TDF_m.sd = TDF_meta[SDTracers$frac] # SD of trophic discrimination factors
  # TDF_p = TDF_proto[Tracers$frac] # trophic discrimination factors
  # TDF_p.sd = TDF_proto[SDTracers$frac] # SD of trophic discrimination factors
  # Sources = Sources
  # Data.sources = Data.sources
  # Tracers = Tracers
  # SDTracers = SDTracers
  # # "Group" = "Group",
  # Use_Dirichlet = FALSE
  # # zoops.f = zoops.f
  # seed = 123
  
  ## Doing the mixing!
  
  ## Use Dirichlet to generate normally distributed mixing parameters
  ## Do not use Dirichlet to explicitly define mixing parameters
  if (Use_Dirichlet == TRUE) {
    ## the following lines are meant to simulate a handful of zooplankton samples
    ## via the Dirichlet distribution
    ## specify the average contribution of each organic matter source to the base 
    ## of the food web in the order they were originally described
    f_all <- rep(1, length(Sources))
    # c(
    # 40/100, # surface particles
    # 5/100, # large particles
    # 40/100, # small particles
    # 5/100  # submicron particles
    # )
    # increase the precision parameter to pull samples around the mean value 
    # decrease the precision parameter to push samples towards the edges of the distribution
    precision <- 0.5
    alpha   <- precision*f_all # shape parameters
    # now we use a Dirichlet distribution to produce a handful of zooplankton with
    # compositions distributed around the mean mixture composition
    zoops.n <- zoops.n # n zoops we'll simulate
    set.seed(seed) # for repeatability
    
    # simulated f using dirichlet
    zoops.f <- data.frame(rdirichlet(n=zoops.n, alpha=alpha)) 
    colnames(zoops.f) <- Sources # columns refer to fractional importance of each particle
    zoops.f$PTS <- runif(n=zoops.n, min = min(PTS), max = max(PTS))
    zoops.f$MTS <- runif(n=zoops.n, min = min(MTS), max = max(MTS))
    zoops.f$FWL <- zoops.f$PTS + zoops.f$MTS
    zoops.f <- zoops.f[order(zoops.f$FWL, zoops.f$PTS, zoops.f$MTS, # order by trophic level
                             -zoops.f[,1], -zoops.f[,2]), ] # then by the first two compositional columns
    row.names(zoops.f) <- seq(1,zoops.n)
  } else {
    ## The following lines allow the user to specify the relative contribution of each organic matter source to each zooplankton sample.
    # enter the % contribution of each organic matter source to each zooplankton sample graphically below: 
    # the columns refer to organic matter sources in the order they were originally specified. 
    # Each row is a zooplankton sample
    
    # Generate a small number of representative samples
    # zoops.f <- data.frame(
    #   t(matrix(nrow = length(Sources),
    #            c(98,01,01,
    #              01,98,01,
    #              01,01,98,
    #              49,49,04,
    #              49,02,49,
    #              02,49,49,
    #              33,33,34)*0.01))
    # )
    
    # Generate a larger number off representative samples by uniformly sampling compositional space
    zoops.f <-
      data.frame("Surface" = NA,
                 "Large" = NA,
                 "Small" = NA,
                 "PTS" = NA,
                 "MTS" = NA)
    # "total" = NA)
    # incr = 0.2
    index = 1
    for (A in seq(0,1,incr)) {
      fA = A
      for (B in seq(0,1-A,incr)) {
        fB = B
        fC = 1-fA-fB
        for(i in PTS) {
          for(j in MTS) {
            # total = fA+fB+fC
            newrow = data.frame("Surface" = fA, "Large" = fB, "Small" = fC, 
                                PTS = i, MTS = j) #, "total" = total)
            zoops.f <- rbind(zoops.f, newrow)
            index = index+1
          }
        }
      }
    }
    zoops.f <- zoops.f[-1,]
    zoops.f$FWL <- zoops.f$PTS + zoops.f$MTS
    
    
    colnames(zoops.f) <- c(Sources,"PTS","MTS","FWL") # columns refer to fractional importance of each particle
    zoops.n <- nrow(zoops.f) # rows refer to samples
    zoops.f <- zoops.f[order(zoops.f$FWL, zoops.f$PTS, zoops.f$MTS, # order by trophic level
                             -zoops.f[,1], -zoops.f[,2]), ] # then by the first two compositional columns
    row.names(zoops.f) <- seq(1,zoops.n)
  }
  
  # can visualize Dirichlet distribution by plotting these on a ternary diagram
  # must pool multiple sources per axis if more than 3 sources are present
  par(mar=rep(0, 4))
  TernaryPlot(alab="f(Small)", blab="f(Large)", clab="f(Surface)")
  TernaryPoints(
    data.frame(A=zoops.f$Small, 
               B=zoops.f$Large, 
               C=zoops.f$Surface))
  
  ## calculating the tracer values after mixing, before any trophic discrimination
  # first calculating the mean tracer value for each organic matter source group
  src.mn <- 
    aggregate(Data.sources[Tracers$all], # aggregate source data
              by=list(Group = Data.sources[["Group"]]), # by organic matter source group
              FUN = mean, na.rm=TRUE)[-1] # taking a mean
  # and the SD
  # we will calculate the standard deviation within each population
  src.SD <- 
    aggregate(Data.sources[Tracers$all], # aggregate source data
              by=list(Group = Data.sources[["Group"]]), # by organic matter source group
              FUN = sd, na.rm=TRUE)[-1] # calculating SD within the population
  # in real data the zooplankton delta values are often outside of the mean
  # values in organic matter sources, so we'll increase the variance of our
  # simulated samples by adding between group variance to the source data.
  src.mn <- src.mn + 
    sign((src.mn - matrix(rep(colMeans(src.mn),3), nrow = 3, byrow = TRUE)))*
    disperse * src.SD
  
  base.sim <- data.frame(matrix(ncol = length(c("PTS","MTS","FWL",Sources,Tracers$all)), nrow = zoops.n))
  colnames(base.sim) <- c("PTS","MTS","FWL",Sources,Tracers$all)
  base.sim[,variables] <- seq(1,zoops.n) # paste("sim ",seq(1,zoops.n))
  base.sim[,c("PTS","MTS",Sources)] <- zoops.f[,c("PTS","MTS",Sources)]
  base.sim$FWL <- base.sim$PTS + base.sim$MTS
  for (i in 1:zoops.n) {
    base.sim[i,Tracers$all] <- colSums(t(zoops.f[i,c(Sources)]) * src.mn)
  }
  
  ## Adding trophic discrimination
  
  # We'll define our some parameters we will need to account for trophic
  # discrimination and propagate error through
  FWL.sd <- MTS.sd
  TDF_m <-    # trophic discrimination factors
    TDF_meta[c(Tracers$frac)]
  TDF_m.sd <- # SD of trophic discrimination factors
    TDF_meta[c(SDTracers$frac)]
  TDF_p <-    # trophic discrimination factors
    TDF_proto[c(Tracers$frac)]
  TDF_p.sd <- # SD of trophic discrimination factors
    TDF_proto[c(SDTracers$frac)]
  
  # initialize empty data frame to store simulated zoop data
  Data.zoops <- 
    data.frame(matrix(ncol = length(c(variables,Tracers$all)), nrow = zoops.n))
  colnames(Data.zoops) <- c(variables,Tracers$all)
  # initialize empty data frame to store zooplankton data before adding process noise
  mn <- 
    data.frame(matrix(ncol = length(Tracers$all), nrow = 1))
  colnames(mn) <- Tracers$all
  var <- 
    data.frame(matrix(ncol = length(Tracers$all), nrow = 1))
  colnames(mn) <- Tracers$all
  
  # Now, in the below loop we can propagate uncertainty or we can use an
  # analytical uncertainty for the step where we pull the sample d15N values
  # from a normal. Propagating uncertainty gives a pretty large spread of
  # values while analytical uncertainty is typically <1 permil, often <0.5.
  # Which method we use really depends on the application of the simulated data
  # and since we are using this for model validation we will by default us the
  # lesser analytical uncertainty. Nonetheless, the code to calculate variance
  # is still included below.
  sd.anl <- 0.52
  
  ## loop to carry out trophic discrimination
  set.seed(987)
  for (iz in 1:nrow(base.sim)) { # zooplankton d15N loop
    # adding trophic discrimination for amino acids with constant TDFs
    mn[c(Tracers$constTDF)] <- 
      base.sim[iz,c(Tracers$constTDF)] + 
      TDF_m[c(Tracers$constTDF)]*base.sim[iz,"FWL"]
    # adding trophic discrimination for amino acids with variable TDFs
    mn[Tracers$varTDF] <- 
      base.sim[iz,Tracers$varTDF] + 
      TDF_p[Tracers$varTDF]*base.sim[iz,"PTS"] +
      TDF_m[Tracers$varTDF]*base.sim[iz,"MTS"]
    # calculating the variance and sd by propagating error through mixing and
    # trophic discrimination
    var[c(Tracers$constTDF)] <- 
      sqrt(
        sum(zoops.f[iz,]^2*sd.anl^2)^2 +
          base.sim[iz,"FWL"]^2*TDF_m.sd[c(SDTracers$constTDF)]^2 +
          TDF_m[c(Tracers$constTDF)]^2*FWL.sd^2
      )
    var[Tracers$varTDF] <- 
      sqrt(
        sum(zoops.f[iz,]^2*sd.anl^2)^2 +
          base.sim[iz,"PTS"]^2*TDF_p.sd[SDTracers$varTDF]^2 +
          TDF_p[Tracers$varTDF]^2*PTS.sd^2 +
          base.sim[iz,"MTS"]^2*TDF_m.sd[SDTracers$varTDF]^2 +
          TDF_m[Tracers$varTDF]^2*MTS.sd^2
      )
    # and drawing the zooplankton d15N from a normal with that mn and sd
    for (j in c(Tracers$frac)) {
      Data.zoops[iz,j] <- rnorm(1, as.numeric(mn[1,j]), 0.001)
    }
    var <- sum(zoops.f[iz,]^2*sd.anl^2)
    # for (k in Tracers$non) {
    #   Data.zoops[iz,k] <- rnorm(1, as.numeric(base.sim[iz,k]), 0.01)
    # }
  }
  Data.zoops[,variables] <- seq(1,zoops.n) # paste("sim ",seq(1,zoops.n))
  Data.zoops[,SDTracers$all] <- sd.anl
  
  # Let's just create one dataframe with both zooplankton compositional and 
  # d15N data in it
  zoops.sim <- cbind(Data.zoops,zoops.f)
  
  # and let's print the data frame
  # datatable(zoops.sim)
  Data.zoops <<- Data.zoops
  zoops.sim <<- zoops.sim
  base.sim <<- base.sim
  src.mn <<- src.mn
  src.SD <<- src.SD
  zoops.f <<- zoops.f
  zoops.n <<- zoops.n
  mn <<- mn
  variance <<- var
  MTS <<- MTS
  PTS <<- PTS
  FWL <<- MTS + PTS
  
}


# ## Simulating zooplankton data! ##
# use this function is d15NPhe is being treated as a conservative tracer in the model
Sim_Zoop_RealPhe <- function(
    zoops.n = 21, # number of zooplankton to simulate
    incr = 0.2, # increment between zooplankton f values
    disperse = 0, # number of SDs by which to increase variance of zooplankton data beyond source means
    PTS, # a vector containing PTS values to simulate
    MTS, # a vector containing MTS values to simulate
    PTS.sd, # a vector containing SD of PTS values
    MTS.sd, # a vector containing SD of PTS values
    TDF_m, # trophic discrimination factors
    TDF_m.sd, # SD of trophic discrimination factors
    TDF_p, # trophic discrimination factors
    TDF_p.sd, # SD of trophic discrimination factors
    Sources, # names of organic matter sources
    Data.sources, # data for organic matter sources
    Tracers, # names of tracers
    SDTracers, # names of tracer SDs
    variables, # additional variables to include
    Use_Dirichlet = TRUE, # should the model usee Dirichlet to simulate zooplankton data?
    zoops.f, # if no Dirichlet, please supply matrix of f values
    seed = 222 # option to set seed
){
  
  # 
  # zoops.n = 50
  # PTS = c(0,1)
  # MTS = c(1,2)
  # PTS.sd = 0.5
  # MTS.sd = 0.5
  # TDF_m = TDF_meta[Tracers$frac] # trophic discrimination factors
  # TDF_m.sd = TDF_meta[SDTracers$frac] # SD of trophic discrimination factors
  # TDF_p = TDF_proto[Tracers$frac] # trophic discrimination factors
  # TDF_p.sd = TDF_proto[SDTracers$frac] # SD of trophic discrimination factors
  # Sources = Sources
  # Data.sources = Data.sources
  # Tracers = Tracers
  # SDTracers = SDTracers
  # # "Group" = "Group",
  # Use_Dirichlet = FALSE
  # # zoops.f = zoops.f
  # seed = 123
  
  ## Doing the mixing!
  
  ## Use Dirichlet to generate normally distributed mixing parameters
  ## Do not use Dirichlet to explicitly define mixing parameters
  if (Use_Dirichlet == TRUE) {
    ## the following lines are meant to simulate a handful of zooplankton samples
    ## via the Dirichlet distribution
    ## specify the average contribution of each organic matter source to the base 
    ## of the food web in the order they were originally described
    f_all <- rep(1, length(Sources))
    # c(
    # 40/100, # surface particles
    # 5/100, # large particles
    # 40/100, # small particles
    # 5/100  # submicron particles
    # )
    # increase the precision parameter to pull samples around the mean value 
    # decrease the precision parameter to push samples towards the edges of the distribution
    precision <- 0.5
    alpha   <- precision*f_all # shape parameters
    # now we use a Dirichlet distribution to produce a handful of zooplankton with
    # compositions distributed around the mean mixture composition
    zoops.n <- zoops.n # n zoops we'll simulate
    set.seed(seed) # for repeatability
    
    # simulated f using dirichlet
    zoops.f <- data.frame(rdirichlet(n=zoops.n, alpha=alpha)) 
    colnames(zoops.f) <- Sources # columns refer to fractional importance of each particle
    zoops.f$PTS <- runif(n=zoops.n, min = min(PTS), max = max(PTS))
    zoops.f$MTS <- runif(n=zoops.n, min = min(MTS), max = max(MTS))
    zoops.f$FWL <- zoops.f$PTS + zoops.f$MTS
    zoops.f <- zoops.f[order(zoops.f$FWL, zoops.f$PTS, zoops.f$MTS, # order by trophic level
                             -zoops.f[,1], -zoops.f[,2]), ] # then by the first two compositional columns
    row.names(zoops.f) <- seq(1,zoops.n)
  } else {
    ## The following lines allow the user to specify the relative contribution of each organic matter source to each zooplankton sample.
    # enter the % contribution of each organic matter source to each zooplankton sample graphically below: 
    # the columns refer to organic matter sources in the order they were originally specified. 
    # Each row is a zooplankton sample
    
    # Generate a small number of representative samples
    # zoops.f <- data.frame(
    #   t(matrix(nrow = length(Sources),
    #            c(98,01,01,
    #              01,98,01,
    #              01,01,98,
    #              49,49,04,
    #              49,02,49,
    #              02,49,49,
    #              33,33,34)*0.01))
    # )
    
    # Generate a larger number off representative samples by uniformly sampling compositional space
    zoops.f <-
      data.frame("Surface" = NA,
                 "Large" = NA,
                 "Small" = NA,
                 "PTS" = NA,
                 "MTS" = NA)
    # "total" = NA)
    # incr = 0.2
    index = 1
    for (A in seq(0,1,incr)) {
      fA = A
      for (B in seq(0,1-A,incr)) {
        fB = B
        fC = 1-fA-fB
        for(i in PTS) {
          for(j in MTS) {
            # total = fA+fB+fC
            newrow = data.frame("Surface" = fA, "Large" = fB, "Small" = fC, 
                                PTS = i, MTS = j) #, "total" = total)
            zoops.f <- rbind(zoops.f, newrow)
            index = index+1
          }
        }
      }
    }
    zoops.f <- zoops.f[-1,]
    zoops.f$FWL <- zoops.f$PTS + zoops.f$MTS
    
    
    colnames(zoops.f) <- c(Sources,"PTS","MTS","FWL") # columns refer to fractional importance of each particle
    zoops.n <- nrow(zoops.f) # rows refer to samples
    zoops.f <- zoops.f[order(zoops.f$FWL, zoops.f$PTS, zoops.f$MTS, # order by trophic level
                             -zoops.f[,1], -zoops.f[,2]), ] # then by the first two compositional columns
    row.names(zoops.f) <- seq(1,zoops.n)
  }
  
  # can visualize Dirichlet distribution by plotting these on a ternary diagram
  # must pool multiple sources per axis if more than 3 sources are present
  par(mar=rep(0, 4))
  TernaryPlot(alab="f(Small)", blab="f(Large)", clab="f(Surface)")
  TernaryPoints(
    data.frame(A=zoops.f$Small, 
               B=zoops.f$Large, 
               C=zoops.f$Surface))
  
  ## calculating the tracer values after mixing, before any trophic discrimination
  # first calculating the mean tracer value for each organic matter source group
  src.mn <- 
    aggregate(Data.sources[Tracers$all], # aggregate source data
              by=list(Group = Data.sources[["Group"]]), # by organic matter source group
              FUN = mean, na.rm=TRUE)[-1] # taking a mean
  # and the SD
  # we will calculate the standard deviation within each population
  src.SD <- 
    aggregate(Data.sources[Tracers$all], # aggregate source data
              by=list(Group = Data.sources[["Group"]]), # by organic matter source group
              FUN = sd, na.rm=TRUE)[-1] # calculating SD within the population
  # in real data the zooplankton delta values are often outside of the mean
  # values in organic matter sources, so we'll increase the variance of our
  # simulated samples by adding between group variance to the source data.
  src.mn <- src.mn + 
    sign((src.mn - matrix(rep(colMeans(src.mn),3), nrow = 3, byrow = TRUE)))*
    disperse * src.SD
  
  base.sim <- data.frame(matrix(ncol = length(c("PTS","MTS","FWL",Sources,Tracers$all)), nrow = zoops.n))
  colnames(base.sim) <- c("PTS","MTS","FWL",Sources,Tracers$all)
  base.sim[,variables] <- seq(1,zoops.n) # paste("sim ",seq(1,zoops.n))
  base.sim[,c("PTS","MTS",Sources)] <- zoops.f[,c("PTS","MTS",Sources)]
  base.sim$FWL <- base.sim$PTS + base.sim$MTS
  for (i in 1:zoops.n) {
    base.sim[i,Tracers$all] <- colSums(t(zoops.f[i,c(Sources)]) * src.mn)
  }
  
  ## Adding trophic discrimination
  
  # We'll define our some parameters we will need to account for trophic
  # discrimination and propagate error through
  FWL.sd <- MTS.sd
  TDF_m <-    # trophic discrimination factors
    TDF_meta[c(Tracers$frac,"d15NPhe")]
  TDF_m.sd <- # SD of trophic discrimination factors
    TDF_meta[c(SDTracers$frac,"SDd15NPhe")]
  TDF_p <-    # trophic discrimination factors
    TDF_proto[c(Tracers$frac,"d15NPhe")]
  TDF_p.sd <- # SD of trophic discrimination factors
    TDF_proto[c(SDTracers$frac,"SDd15NPhe")]
  
  # initialize empty data frame to store simulated zoop data
  Data.zoops <- 
    data.frame(matrix(ncol = length(c(variables,Tracers$all)), nrow = zoops.n))
  colnames(Data.zoops) <- c(variables,Tracers$all)
  # initialize empty data frame to store zooplankton data before adding process noise
  mn <- 
    data.frame(matrix(ncol = length(Tracers$all), nrow = 1))
  colnames(mn) <- Tracers$all
  var <- 
    data.frame(matrix(ncol = length(Tracers$all), nrow = 1))
  colnames(mn) <- Tracers$all
  
  # Now, in the below loop we can propagate uncertainty or we can use an
  # analytical uncertainty for the step where we pull the sample d15N values
  # from a normal. Propagating uncertainty gives a pretty large spread of
  # values while analytical uncertainty is typically <1 permil, often <0.5.
  # Which method we use really depends on the application of the simulated data
  # and since we are using this for model validation we will by default us the
  # lesser analytical uncertainty. Nonetheless, the code to calculate variance
  # is still included below.
  sd.anl <- 0.52
  
  ## loop to carry out trophic discrimination
  set.seed(987)
  for (iz in 1:nrow(base.sim)) { # zooplankton d15N loop
    # adding trophic discrimination for amino acids with constant TDFs
    mn[c(Tracers$constTDF,"d15NPhe")] <- 
      base.sim[iz,c(Tracers$constTDF,"d15NPhe")] + 
      TDF_m[c(Tracers$constTDF,"d15NPhe")]*base.sim[iz,"FWL"]
    # adding trophic discrimination for amino acids with variable TDFs
    mn[Tracers$varTDF] <- 
      base.sim[iz,Tracers$varTDF] + 
      TDF_p[Tracers$varTDF]*base.sim[iz,"PTS"] +
      TDF_m[Tracers$varTDF]*base.sim[iz,"MTS"]
    # calculating the variance and sd by propagating error through mixing and
    # trophic discrimination
    var[c(Tracers$constTDF,"d15NPhe")] <- 
      sqrt(
        sum(zoops.f[iz,]^2*sd.anl^2)^2 +
          base.sim[iz,"FWL"]^2*TDF_m.sd[c(SDTracers$constTDF,"SDd15NPhe")]^2 +
          TDF_m[c(Tracers$constTDF,"d15NPhe")]^2*FWL.sd^2
      )
    var[Tracers$varTDF] <- 
      sqrt(
        sum(zoops.f[iz,]^2*sd.anl^2)^2 +
          base.sim[iz,"PTS"]^2*TDF_p.sd[SDTracers$varTDF]^2 +
          TDF_p[Tracers$varTDF]^2*PTS.sd^2 +
          base.sim[iz,"MTS"]^2*TDF_m.sd[SDTracers$varTDF]^2 +
          TDF_m[Tracers$varTDF]^2*MTS.sd^2
      )
    # and drawing the zooplankton d15N from a normal with that mn and sd
    for (j in c(Tracers$frac,"d15NPhe")) {
      Data.zoops[iz,j] <- rnorm(1, as.numeric(mn[1,j]), 0.001)
    }
    var <- sum(zoops.f[iz,]^2*sd.anl^2)
    # for (k in Tracers$non) {
    #   Data.zoops[iz,k] <- rnorm(1, as.numeric(base.sim[iz,k]), 0.01)
    # }
  }
  Data.zoops[,variables] <- seq(1,zoops.n) # paste("sim ",seq(1,zoops.n))
  Data.zoops[,SDTracers$all] <- sd.anl
  
  # Let's just create one dataframe with both zooplankton compositional and 
  # d15N data in it
  zoops.sim <- cbind(Data.zoops,zoops.f)
  
  # and let's print the data frame
  # datatable(zoops.sim)
  Data.zoops <<- Data.zoops
  zoops.sim <<- zoops.sim
  base.sim <<- base.sim
  src.mn <<- src.mn
  src.SD <<- src.SD
  zoops.f <<- zoops.f
  zoops.n <<- zoops.n
  mn <<- mn
  variance <<- var
  MTS <<- MTS
  PTS <<- PTS
  FWL <<- MTS + PTS
  
}

# ## Simulating zooplankton data! ##
# use this function if d15NPhe and d15NLys are being treated as a conservative tracer in the model
Sim_Zoop_RealPheLys <- function(
    zoops.n = 21, # number of zooplankton to simulate
    incr = 0.2, # increment between zooplankton f values
    disperse = 0, # number of SDs by which to increase variance of zooplankton data beyond source means
    PTS, # a vector containing PTS values to simulate
    MTS, # a vector containing MTS values to simulate
    PTS.sd, # a vector containing SD of PTS values
    MTS.sd, # a vector containing SD of PTS values
    TDF_m, # trophic discrimination factors
    TDF_m.sd, # SD of trophic discrimination factors
    TDF_p, # trophic discrimination factors
    TDF_p.sd, # SD of trophic discrimination factors
    Sources, # names of organic matter sources
    Data.sources, # data for organic matter sources
    Tracers, # names of tracers
    SDTracers, # names of tracer SDs
    variables, # additional variables to include
    Use_Dirichlet = TRUE, # should the model usee Dirichlet to simulate zooplankton data?
    zoops.f, # if no Dirichlet, please supply matrix of f values
    seed = 222 # option to set seed
){
  
  # 
  # zoops.n = 50
  # PTS = c(0,1)
  # MTS = c(1,2)
  # PTS.sd = 0.5
  # MTS.sd = 0.5
  # TDF_m = TDF_meta[Tracers$frac] # trophic discrimination factors
  # TDF_m.sd = TDF_meta[SDTracers$frac] # SD of trophic discrimination factors
  # TDF_p = TDF_proto[Tracers$frac] # trophic discrimination factors
  # TDF_p.sd = TDF_proto[SDTracers$frac] # SD of trophic discrimination factors
  # Sources = Sources
  # Data.sources = Data.sources
  # Tracers = Tracers
  # SDTracers = SDTracers
  # # "Group" = "Group",
  # Use_Dirichlet = FALSE
  # # zoops.f = zoops.f
  # seed = 123
  
  ## Doing the mixing!
  
  ## Use Dirichlet to generate normally distributed mixing parameters
  ## Do not use Dirichlet to explicitly define mixing parameters
  if (Use_Dirichlet == TRUE) {
    ## the following lines are meant to simulate a handful of zooplankton samples
    ## via the Dirichlet distribution
    ## specify the average contribution of each organic matter source to the base 
    ## of the food web in the order they were originally described
    f_all <- rep(1, length(Sources))
    # c(
    # 40/100, # surface particles
    # 5/100, # large particles
    # 40/100, # small particles
    # 5/100  # submicron particles
    # )
    # increase the precision parameter to pull samples around the mean value 
    # decrease the precision parameter to push samples towards the edges of the distribution
    precision <- 0.5
    alpha   <- precision*f_all # shape parameters
    # now we use a Dirichlet distribution to produce a handful of zooplankton with
    # compositions distributed around the mean mixture composition
    zoops.n <- zoops.n # n zoops we'll simulate
    set.seed(seed) # for repeatability
    
    # simulated f using dirichlet
    zoops.f <- data.frame(rdirichlet(n=zoops.n, alpha=alpha)) 
    colnames(zoops.f) <- Sources # columns refer to fractional importance of each particle
    zoops.f$PTS <- runif(n=zoops.n, min = min(PTS), max = max(PTS))
    zoops.f$MTS <- runif(n=zoops.n, min = min(MTS), max = max(MTS))
    zoops.f$FWL <- zoops.f$PTS + zoops.f$MTS
    zoops.f <- zoops.f[order(zoops.f$FWL, zoops.f$PTS, zoops.f$MTS, # order by trophic level
                             -zoops.f[,1], -zoops.f[,2]), ] # then by the first two compositional columns
    row.names(zoops.f) <- seq(1,zoops.n)
  } else {
    ## The following lines allow the user to specify the relative contribution of each organic matter source to each zooplankton sample.
    # enter the % contribution of each organic matter source to each zooplankton sample graphically below: 
    # the columns refer to organic matter sources in the order they were originally specified. 
    # Each row is a zooplankton sample
    
    # Generate a small number of representative samples
    # zoops.f <- data.frame(
    #   t(matrix(nrow = length(Sources),
    #            c(98,01,01,
    #              01,98,01,
    #              01,01,98,
    #              49,49,04,
    #              49,02,49,
    #              02,49,49,
    #              33,33,34)*0.01))
    # )
    
    # Generate a larger number off representative samples by uniformly sampling compositional space
    zoops.f <-
      data.frame("Surface" = NA,
                 "Large" = NA,
                 "Small" = NA,
                 "PTS" = NA,
                 "MTS" = NA)
    # "total" = NA)
    # incr = 0.2
    index = 1
    for (A in seq(0,1,incr)) {
      fA = A
      for (B in seq(0,1-A,incr)) {
        fB = B
        fC = 1-fA-fB
        for(i in PTS) {
          for(j in MTS) {
            # total = fA+fB+fC
            newrow = data.frame("Surface" = fA, "Large" = fB, "Small" = fC, 
                                PTS = i, MTS = j) #, "total" = total)
            zoops.f <- rbind(zoops.f, newrow)
            index = index+1
          }
        }
      }
    }
    zoops.f <- zoops.f[-1,]
    zoops.f$FWL <- zoops.f$PTS + zoops.f$MTS
    
    
    colnames(zoops.f) <- c(Sources,"PTS","MTS","FWL") # columns refer to fractional importance of each particle
    zoops.n <- nrow(zoops.f) # rows refer to samples
    zoops.f <- zoops.f[order(zoops.f$FWL, zoops.f$PTS, zoops.f$MTS, # order by trophic level
                             -zoops.f[,1], -zoops.f[,2]), ] # then by the first two compositional columns
    row.names(zoops.f) <- seq(1,zoops.n)
  }
  
  # can visualize Dirichlet distribution by plotting these on a ternary diagram
  # must pool multiple sources per axis if more than 3 sources are present
  par(mar=rep(0, 4))
  TernaryPlot(alab="f(Small)", blab="f(Large)", clab="f(Surface)")
  TernaryPoints(
    data.frame(A=zoops.f$Small, 
               B=zoops.f$Large, 
               C=zoops.f$Surface))
  
  ## calculating the tracer values after mixing, before any trophic discrimination
  # first calculating the mean tracer value for each organic matter source group
  src.mn <- 
    aggregate(Data.sources[Tracers$all], # aggregate source data
              by=list(Group = Data.sources[["Group"]]), # by organic matter source group
              FUN = mean, na.rm=TRUE)[-1] # taking a mean
  # and the SD
  # we will calculate the standard deviation within each population
  src.SD <- 
    aggregate(Data.sources[Tracers$all], # aggregate source data
              by=list(Group = Data.sources[["Group"]]), # by organic matter source group
              FUN = sd, na.rm=TRUE)[-1] # calculating SD within the population
  # in real data the zooplankton delta values are often outside of the mean
  # values in organic matter sources, so we'll increase the variance of our
  # simulated samples by adding between group variance to the source data.
  src.mn <- src.mn + 
    sign((src.mn - matrix(rep(colMeans(src.mn),3), nrow = 3, byrow = TRUE)))*
    disperse * src.SD
  
  base.sim <- data.frame(matrix(ncol = length(c("PTS","MTS","FWL",Sources,Tracers$all)), nrow = zoops.n))
  colnames(base.sim) <- c("PTS","MTS","FWL",Sources,Tracers$all)
  base.sim[,variables] <- seq(1,zoops.n) # paste("sim ",seq(1,zoops.n))
  base.sim[,c("PTS","MTS",Sources)] <- zoops.f[,c("PTS","MTS",Sources)]
  base.sim$FWL <- base.sim$PTS + base.sim$MTS
  for (i in 1:zoops.n) {
    base.sim[i,Tracers$all] <- colSums(t(zoops.f[i,c(Sources)]) * src.mn)
  }
  
  ## Adding trophic discrimination
  
  # We'll define our some parameters we will need to account for trophic
  # discrimination and propagate error through
  FWL.sd <- MTS.sd
  TDF_m <-    # trophic discrimination factors
    TDF_meta[c(Tracers$frac,"d15NPhe","d15NLys")]
  TDF_m.sd <- # SD of trophic discrimination factors
    TDF_meta[c(SDTracers$frac,"SDd15NPhe","SDd15NLys")]
  TDF_p <-    # trophic discrimination factors
    TDF_proto[c(Tracers$frac,"d15NPhe","d15NLys")]
  TDF_p.sd <- # SD of trophic discrimination factors
    TDF_proto[c(SDTracers$frac,"SDd15NPhe","SDd15NLys")]
  
  # initialize empty data frame to store simulated zoop data
  Data.zoops <- 
    data.frame(matrix(ncol = length(c(variables,Tracers$all)), nrow = zoops.n))
  colnames(Data.zoops) <- c(variables,Tracers$all)
  # initialize empty data frame to store zooplankton data before adding process noise
  mn <- 
    data.frame(matrix(ncol = length(Tracers$all), nrow = 1))
  colnames(mn) <- Tracers$all
  var <- 
    data.frame(matrix(ncol = length(Tracers$all), nrow = 1))
  colnames(mn) <- Tracers$all
  
  # Now, in the below loop we can propagate uncertainty or we can use an
  # analytical uncertainty for the step where we pull the sample d15N values
  # from a normal. Propagating uncertainty gives a pretty large spread of
  # values while analytical uncertainty is typically <1 permil, often <0.5.
  # Which method we use really depends on the application of the simulated data
  # and since we are using this for model validation we will by default us the
  # lesser analytical uncertainty. Nonetheless, the code to calculate variance
  # is still included below.
  sd.anl <- 0.52
  
  ## loop to carry out trophic discrimination
  set.seed(987)
  for (iz in 1:nrow(base.sim)) { # zooplankton d15N loop
    # adding trophic discrimination for amino acids with constant TDFs
    mn[c(Tracers$constTDF,"d15NPhe","d15NLys")] <- 
      base.sim[iz,c(Tracers$constTDF,"d15NPhe","d15NLys")] + 
      TDF_m[c(Tracers$constTDF,"d15NPhe","d15NLys")]*base.sim[iz,"FWL"]
    # adding trophic discrimination for amino acids with variable TDFs
    mn[Tracers$varTDF] <- 
      base.sim[iz,Tracers$varTDF] + 
      TDF_p[Tracers$varTDF]*base.sim[iz,"PTS"] +
      TDF_m[Tracers$varTDF]*base.sim[iz,"MTS"]
    # calculating the variance and sd by propagating error through mixing and
    # trophic discrimination
    var[c(Tracers$constTDF,"d15NPhe","d15NLys")] <- 
      sqrt(
        sum(zoops.f[iz,]^2*sd.anl^2)^2 +
          base.sim[iz,"FWL"]^2*TDF_m.sd[c(SDTracers$constTDF,"SDd15NPhe","SDd15NLys")]^2 +
          TDF_m[c(Tracers$constTDF,"d15NPhe","d15NLys")]^2*FWL.sd^2
      )
    var[Tracers$varTDF] <- 
      sqrt(
        sum(zoops.f[iz,]^2*sd.anl^2)^2 +
          base.sim[iz,"PTS"]^2*TDF_p.sd[SDTracers$varTDF]^2 +
          TDF_p[Tracers$varTDF]^2*PTS.sd^2 +
          base.sim[iz,"MTS"]^2*TDF_m.sd[SDTracers$varTDF]^2 +
          TDF_m[Tracers$varTDF]^2*MTS.sd^2
      )
    # and drawing the zooplankton d15N from a normal with that mn and sd
    for (j in c(Tracers$frac,"d15NPhe","d15NLys")) {
      Data.zoops[iz,j] <- rnorm(1, as.numeric(mn[1,j]), 0.001)
    }
    var <- sum(zoops.f[iz,]^2*sd.anl^2)
    # for (k in Tracers$non) {
    #   Data.zoops[iz,k] <- rnorm(1, as.numeric(base.sim[iz,k]), 0.01)
    # }
  }
  Data.zoops[,variables] <- seq(1,zoops.n) # paste("sim ",seq(1,zoops.n))
  Data.zoops[,SDTracers$all] <- sd.anl
  
  # Let's just create one dataframe with both zooplankton compositional and 
  # d15N data in it
  zoops.sim <- cbind(Data.zoops,zoops.f)
  
  # and let's print the data frame
  # datatable(zoops.sim)
  Data.zoops <<- Data.zoops
  zoops.sim <<- zoops.sim
  base.sim <<- base.sim
  src.mn <<- src.mn
  src.SD <<- src.SD
  zoops.f <<- zoops.f
  zoops.n <<- zoops.n
  mn <<- mn
  variance <<- var
  MTS <<- MTS
  PTS <<- PTS
  FWL <<- MTS + PTS
  
}

# ## Simulating zooplankton data! ##
# use this function if all SAAs are being treated as a conservative tracer in the model
Sim_Zoop_RealSAAs <- function(
    zoops.n = 21, # number of zooplankton to simulate
    incr = 0.2, # increment between zooplankton f values
    disperse = 0, # number of SDs by which to increase variance of zooplankton data beyond source means
    PTS, # a vector containing PTS values to simulate
    MTS, # a vector containing MTS values to simulate
    PTS.sd, # a vector containing SD of PTS values
    MTS.sd, # a vector containing SD of PTS values
    TDF_m, # trophic discrimination factors
    TDF_m.sd, # SD of trophic discrimination factors
    TDF_p, # trophic discrimination factors
    TDF_p.sd, # SD of trophic discrimination factors
    Sources, # names of organic matter sources
    Data.sources, # data for organic matter sources
    Tracers, # names of tracers
    SDTracers, # names of tracer SDs
    variables, # additional variables to include
    Use_Dirichlet = TRUE, # should the model usee Dirichlet to simulate zooplankton data?
    zoops.f, # if no Dirichlet, please supply matrix of f values
    seed = 222 # option to set seed
){
  
  # 
  # zoops.n = 50
  # PTS = c(0,1)
  # MTS = c(1,2)
  # PTS.sd = 0.5
  # MTS.sd = 0.5
  # TDF_m = TDF_meta[Tracers$frac] # trophic discrimination factors
  # TDF_m.sd = TDF_meta[SDTracers$frac] # SD of trophic discrimination factors
  # TDF_p = TDF_proto[Tracers$frac] # trophic discrimination factors
  # TDF_p.sd = TDF_proto[SDTracers$frac] # SD of trophic discrimination factors
  # Sources = Sources
  # Data.sources = Data.sources
  # Tracers = Tracers
  # SDTracers = SDTracers
  # # "Group" = "Group",
  # Use_Dirichlet = FALSE
  # # zoops.f = zoops.f
  # seed = 123
  
  ## Doing the mixing!
  
  ## Use Dirichlet to generate normally distributed mixing parameters
  ## Do not use Dirichlet to explicitly define mixing parameters
  if (Use_Dirichlet == TRUE) {
    ## the following lines are meant to simulate a handful of zooplankton samples
    ## via the Dirichlet distribution
    ## specify the average contribution of each organic matter source to the base 
    ## of the food web in the order they were originally described
    f_all <- rep(1, length(Sources))
    # c(
    # 40/100, # surface particles
    # 5/100, # large particles
    # 40/100, # small particles
    # 5/100  # submicron particles
    # )
    # increase the precision parameter to pull samples around the mean value 
    # decrease the precision parameter to push samples towards the edges of the distribution
    precision <- 0.5
    alpha   <- precision*f_all # shape parameters
    # now we use a Dirichlet distribution to produce a handful of zooplankton with
    # compositions distributed around the mean mixture composition
    zoops.n <- zoops.n # n zoops we'll simulate
    set.seed(seed) # for repeatability
    
    # simulated f using dirichlet
    zoops.f <- data.frame(rdirichlet(n=zoops.n, alpha=alpha)) 
    colnames(zoops.f) <- Sources # columns refer to fractional importance of each particle
    zoops.f$PTS <- runif(n=zoops.n, min = min(PTS), max = max(PTS))
    zoops.f$MTS <- runif(n=zoops.n, min = min(MTS), max = max(MTS))
    zoops.f$FWL <- zoops.f$PTS + zoops.f$MTS
    zoops.f <- zoops.f[order(zoops.f$FWL, zoops.f$PTS, zoops.f$MTS, # order by trophic level
                             -zoops.f[,1], -zoops.f[,2]), ] # then by the first two compositional columns
    row.names(zoops.f) <- seq(1,zoops.n)
  } else {
    ## The following lines allow the user to specify the relative contribution of each organic matter source to each zooplankton sample.
    # enter the % contribution of each organic matter source to each zooplankton sample graphically below: 
    # the columns refer to organic matter sources in the order they were originally specified. 
    # Each row is a zooplankton sample
    
    # Generate a small number of representative samples
    # zoops.f <- data.frame(
    #   t(matrix(nrow = length(Sources),
    #            c(98,01,01,
    #              01,98,01,
    #              01,01,98,
    #              49,49,04,
    #              49,02,49,
    #              02,49,49,
    #              33,33,34)*0.01))
    # )
    
    # Generate a larger number off representative samples by uniformly sampling compositional space
    zoops.f <-
      data.frame("Surface" = NA,
                 "Large" = NA,
                 "Small" = NA,
                 "PTS" = NA,
                 "MTS" = NA)
    # "total" = NA)
    # incr = 0.2
    index = 1
    for (A in seq(0,1,incr)) {
      fA = A
      for (B in seq(0,1-A,incr)) {
        fB = B
        fC = 1-fA-fB
        for(i in PTS) {
          for(j in MTS) {
            # total = fA+fB+fC
            newrow = data.frame("Surface" = fA, "Large" = fB, "Small" = fC, 
                                PTS = i, MTS = j) #, "total" = total)
            zoops.f <- rbind(zoops.f, newrow)
            index = index+1
          }
        }
      }
    }
    zoops.f <- zoops.f[-1,]
    zoops.f$FWL <- zoops.f$PTS + zoops.f$MTS
    
    
    colnames(zoops.f) <- c(Sources,"PTS","MTS","FWL") # columns refer to fractional importance of each particle
    zoops.n <- nrow(zoops.f) # rows refer to samples
    zoops.f <- zoops.f[order(zoops.f$FWL, zoops.f$PTS, zoops.f$MTS, # order by trophic level
                             -zoops.f[,1], -zoops.f[,2]), ] # then by the first two compositional columns
    row.names(zoops.f) <- seq(1,zoops.n)
  }
  
  # can visualize Dirichlet distribution by plotting these on a ternary diagram
  # must pool multiple sources per axis if more than 3 sources are present
  par(mar=rep(0, 4))
  TernaryPlot(alab="f(Small)", blab="f(Large)", clab="f(Surface)")
  TernaryPoints(
    data.frame(A=zoops.f$Small, 
               B=zoops.f$Large, 
               C=zoops.f$Surface))
  
  ## calculating the tracer values after mixing, before any trophic discrimination
  # first calculating the mean tracer value for each organic matter source group
  src.mn <- 
    aggregate(Data.sources[Tracers$all], # aggregate source data
              by=list(Group = Data.sources[["Group"]]), # by organic matter source group
              FUN = mean, na.rm=TRUE)[-1] # taking a mean
  # and the SD
  # we will calculate the standard deviation within each population
  src.SD <- 
    aggregate(Data.sources[Tracers$all], # aggregate source data
              by=list(Group = Data.sources[["Group"]]), # by organic matter source group
              FUN = sd, na.rm=TRUE)[-1] # calculating SD within the population
  # in real data the zooplankton delta values are often outside of the mean
  # values in organic matter sources, so we'll increase the variance of our
  # simulated samples by adding between group variance to the source data.
  src.mn <- src.mn + 
    sign((src.mn - matrix(rep(colMeans(src.mn),3), nrow = 3, byrow = TRUE)))*
    disperse * src.SD
  
  base.sim <- data.frame(matrix(ncol = length(c("PTS","MTS","FWL",Sources,Tracers$all)), nrow = zoops.n))
  colnames(base.sim) <- c("PTS","MTS","FWL",Sources,Tracers$all)
  base.sim[,variables] <- seq(1,zoops.n) # paste("sim ",seq(1,zoops.n))
  base.sim[,c("PTS","MTS",Sources)] <- zoops.f[,c("PTS","MTS",Sources)]
  base.sim$FWL <- base.sim$PTS + base.sim$MTS
  for (i in 1:zoops.n) {
    base.sim[i,Tracers$all] <- colSums(t(zoops.f[i,c(Sources)]) * src.mn)
  }
  
  ## Adding trophic discrimination
  
  # We'll define our some parameters we will need to account for trophic
  # discrimination and propagate error through
  FWL.sd <- MTS.sd
  TDF_m <-    # trophic discrimination factors
    TDF_meta[c(Tracers$frac,"d15NPhe","d15NLys","d15NGly","d15NSer")]
  TDF_m.sd <- # SD of trophic discrimination factors
    TDF_meta[c(SDTracers$frac,"SDd15NPhe","SDd15NLys")]
  TDF_p <-    # trophic discrimination factors
    TDF_proto[c(Tracers$frac,"d15NPhe","d15NLys","d15NGly","d15NSer")]
  TDF_p.sd <- # SD of trophic discrimination factors
    TDF_proto[c(SDTracers$frac,"SDd15NPhe","SDd15NLys")]
  
  # initialize empty data frame to store simulated zoop data
  Data.zoops <- 
    data.frame(matrix(ncol = length(c(variables,Tracers$all)), nrow = zoops.n))
  colnames(Data.zoops) <- c(variables,Tracers$all)
  # initialize empty data frame to store zooplankton data before adding process noise
  mn <- 
    data.frame(matrix(ncol = length(Tracers$all), nrow = 1))
  colnames(mn) <- Tracers$all
  var <- 
    data.frame(matrix(ncol = length(Tracers$all), nrow = 1))
  colnames(mn) <- Tracers$all
  
  # Now, in the below loop we can propagate uncertainty or we can use an
  # analytical uncertainty for the step where we pull the sample d15N values
  # from a normal. Propagating uncertainty gives a pretty large spread of
  # values while analytical uncertainty is typically <1 permil, often <0.5.
  # Which method we use really depends on the application of the simulated data
  # and since we are using this for model validation we will by default us the
  # lesser analytical uncertainty. Nonetheless, the code to calculate variance
  # is still included below.
  sd.anl <- 0.52
  
  ## loop to carry out trophic discrimination
  set.seed(987)
  for (iz in 1:nrow(base.sim)) { # zooplankton d15N loop
    # adding trophic discrimination for amino acids with constant TDFs
    mn[c(Tracers$constTDF,"d15NPhe","d15NLys","d15NGly","d15NSer")] <- 
      base.sim[iz,c(Tracers$constTDF,"d15NPhe","d15NLys","d15NGly","d15NSer")] + 
      TDF_m[c(Tracers$constTDF,"d15NPhe","d15NLys","d15NGly","d15NSer")]*base.sim[iz,"FWL"]
    # adding trophic discrimination for amino acids with variable TDFs
    mn[Tracers$varTDF] <- 
      base.sim[iz,Tracers$varTDF] + 
      TDF_p[Tracers$varTDF]*base.sim[iz,"PTS"] +
      TDF_m[Tracers$varTDF]*base.sim[iz,"MTS"]
    # calculating the variance and sd by propagating error through mixing and
    # trophic discrimination
    var[c(Tracers$constTDF,"d15NPhe","d15NLys","d15NGly","d15NSer")] <- 
      sqrt(
        sum(zoops.f[iz,]^2*sd.anl^2)^2 +
          base.sim[iz,"FWL"]^2*TDF_m.sd[c(SDTracers$constTDF,"SDd15NPhe","SDd15NLys")]^2 +
          TDF_m[c(Tracers$constTDF,"d15NPhe","d15NLys","d15NGly","d15NSer")]^2*FWL.sd^2
      )
    var[Tracers$varTDF] <- 
      sqrt(
        sum(zoops.f[iz,]^2*sd.anl^2)^2 +
          base.sim[iz,"PTS"]^2*TDF_p.sd[SDTracers$varTDF]^2 +
          TDF_p[Tracers$varTDF]^2*PTS.sd^2 +
          base.sim[iz,"MTS"]^2*TDF_m.sd[SDTracers$varTDF]^2 +
          TDF_m[Tracers$varTDF]^2*MTS.sd^2
      )
    # and drawing the zooplankton d15N from a normal with that mn and sd
    for (j in c(Tracers$frac,"d15NPhe","d15NLys","d15NGly","d15NSer")) {
      Data.zoops[iz,j] <- rnorm(1, as.numeric(mn[1,j]), 0.001)
    }
    var <- sum(zoops.f[iz,]^2*sd.anl^2)
    # for (k in Tracers$non) {
    #   Data.zoops[iz,k] <- rnorm(1, as.numeric(base.sim[iz,k]), 0.01)
    # }
  }
  Data.zoops[,variables] <- seq(1,zoops.n) # paste("sim ",seq(1,zoops.n))
  Data.zoops[,SDTracers$all] <- sd.anl
  
  # Let's just create one dataframe with both zooplankton compositional and 
  # d15N data in it
  zoops.sim <- cbind(Data.zoops,zoops.f)
  
  # and let's print the data frame
  # datatable(zoops.sim)
  Data.zoops <<- Data.zoops
  zoops.sim <<- zoops.sim
  base.sim <<- base.sim
  src.mn <<- src.mn
  src.SD <<- src.SD
  zoops.f <<- zoops.f
  zoops.n <<- zoops.n
  mn <<- mn
  variance <<- var
  MTS <<- MTS
  PTS <<- PTS
  FWL <<- MTS + PTS
  
}