#### Simulating zooplankton data! ####

Sim_Zoop <- function(
    nzoops = 21, # number of zooplankton to simulate
    incr = 0.2, # increment between zooplankton f values
    disperse = 0, # number of SDs by which to increase variance of zooplankton data beyond source means
    PTS, # a vector containing PTS values to simulate
    MTS, # a vector containing MTS values to simulate
    TDF_m, # trophic discrimination factors
    TDF_p, # trophic discrimination factors
    Sources, # names of organic matter sources
    Data.sources, # data for organic matter sources
    Tracers, # names of tracers
    Variables, # additional variables to include
    Random_Samples = TRUE, # should the model use Dirichlet to simulate zooplankton data?
    seed = 222 # option to set seed
){
  
  ## Use Dirichlet to generate normally distributed mixing parameters
  ## Do not use Dirichlet to explicitly define mixing parameters
  if (Random_Samples == TRUE) {
    ## the following lines are meant to simulate a handful of zooplankton samples
    ## via the Dirichlet distribution
    ## specify the average contribution of each organic matter source to the base 
    ## of the food web in the order they were originally described
    f_all <- rep(1, length(Sources))
    # increase the precision parameter to pull samples around the mean value 
    # decrease the precision parameter to push samples towards the edges of the distribution
    precision <- 0.5
    alpha   <- precision*f_all # shape parameters
    # now we use a Dirichlet distribution to produce a handful of zooplankton with
    # compositions distributed around the mean mixture composition
    set.seed(seed) # for repeatability
    # simulated f using dirichlet
    zoops.f <- data.frame(rdirichlet(n=nzoops, alpha=alpha)) 
    colnames(zoops.f) <- Sources # columns refer to fractional importance of each particle
    zoops.f$PTS <- runif(n=nzoops, min = min(PTS), max = max(PTS))
    zoops.f$MTS <- runif(n=nzoops, min = min(MTS), max = max(MTS))
    zoops.f$FWL <- zoops.f$PTS + zoops.f$MTS
    zoops.f <- zoops.f[order(zoops.f$FWL, zoops.f$PTS, zoops.f$MTS, # order by trophic level
                             -zoops.f[,1], -zoops.f[,2]), ] # then by the first two compositional columns
    row.names(zoops.f) <- seq(1,nzoops)
  } else {
    # OTHERWISE
    # Generate a larger number of representative samples by uniformly sampling compositional space
    zoops.f <-
      data.frame("Surface" = NA,
                 "Large" = NA,
                 "Small" = NA,
                 "PTS" = NA,
                 "MTS" = NA)
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
    nzoops <- nrow(zoops.f) # rows refer to samples
    zoops.f <- zoops.f[order(zoops.f$FWL, zoops.f$PTS, zoops.f$MTS, # order by trophic level
                             -zoops.f[,1], -zoops.f[,2]), ] # then by the first two compositional columns
    row.names(zoops.f) <- seq(1,nzoops)
  }
  
  # can visualize Dirichlet distribution by plotting these on a ternary diagram
  # must pool multiple sources per axis if more than 3 sources are present
  if(length(Sources) == 3){
    par(mar=rep(0, 4))
    TernaryPlot(alab=Sources[1], blab=Sources[2], clab=Sources[3])
    TernaryPoints(
      data.frame(A=zoops.f[Sources[1]], 
                 B=zoops.f[Sources[2]], 
                 C=zoops.f[Sources[3]]))
  }else{
    if(length(Sources) == 2){
      plot(
        x = zoops.f[[Sources[1]]],
        y = zoops.f[[Sources[2]]],
        xlab = Sources[1], 
        ylab = Sources[2])
    }
  }
  par(mar=rep(2, 4))
  hist(zoops.f$PTS, main = "Protistan Trophic Steps")
  hist(zoops.f$MTS, main = "Metazoan Trophic Steps")
  # print(datatable(round(zoops.f,2)))
  
  ## calculating the tracer values after mixing, before any trophic discrimination
  # calculating the mean tracer value for each organic matter source group
  src.mn <- 
    aggregate(Data.sources[Tracers$all], # aggregate source data
              by=list(Source = Data.sources[["Source"]]), # by organic matter source group
              FUN = mean, na.rm=TRUE)[-1] # taking a mean
  # we will calculate the standard deviation within each population
  src.SD <- 
    aggregate(Data.sources[Tracers$all], # aggregate source data
              by=list(Source = Data.sources[["Source"]]), # by organic matter source group
              FUN = sd, na.rm=TRUE)[-1] # calculating SD within the population
  
  # in real data the zooplankton delta values are often outside of the mean
  # values in organic matter sources, so we'll increase the variance of our
  # simulated samples by adding between group variance to the source data.
  src.mn <- src.mn + 
    sign((src.mn - matrix(rep(colMeans(src.mn),3), nrow = 3, byrow = TRUE)))*
    disperse * src.SD
  
  # now we'll calculate tracer values for each sample at the base of the food web
  base.sim <- data.frame(matrix(ncol = length(c("PTS","MTS","FWL",Sources,Tracers$all)), nrow = nzoops))
  colnames(base.sim) <- c("PTS","MTS","FWL",Sources,Tracers$all)
  base.sim[,Variables] <- seq(1,nzoops) # paste("sim ",seq(1,nzoops))
  base.sim[,c("PTS","MTS",Sources)] <- zoops.f[,c("PTS","MTS",Sources)]
  base.sim$FWL <- base.sim$PTS + base.sim$MTS
  for (i in 1:nzoops) {
    base.sim[i,Tracers$all] <- colSums(t(zoops.f[i,c(Sources)]) * src.mn)
  }
  
  ## Adding trophic discrimination
  
  # We'll define our some parameters we will need to account for trophic
  TDF_m <-    # trophic discrimination factors
    TDF_meta[c(Tracers$frac)]
  TDF_p <-    # trophic discrimination factors
    TDF_proto[c(Tracers$frac)]
  
  # initialize empty data frame to store simulated zoop data
  Data.zoops <- 
    data.frame(matrix(ncol = length(c(Variables,Tracers$all)), nrow = nzoops))
  colnames(Data.zoops) <- c(Variables,Tracers$all)
  # initialize empty data frame to store zooplankton data before adding process noise
  mn <- 
    data.frame(matrix(ncol = length(Tracers$all), nrow = 1))
  colnames(mn) <- Tracers$all
  
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
    Data.zoops[iz,Tracers$frac] <- 
      as.numeric(mn[1,Tracers$frac])
  }
  Data.zoops[,Variables] <- seq(1,nzoops) # paste("sim ",seq(1,nzoops))
  Data.zoops[,SDTracers$all] <- 0.5
  
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
  nzoops <<- nzoops
  mn <<- mn
  MTS <<- MTS
  PTS <<- PTS
  FWL <<- MTS + PTS
  
}




#### Simulating zooplankton data! ####
# use this function is d15NPhe is being treated as a conservative tracer in the model
Sim_Zoop_RealPhe <- function(
    nzoops = 21, # number of zooplankton to simulate
    incr = 0.2, # increment between zooplankton f values
    disperse = 0, # number of SDs by which to increase variance of zooplankton data beyond source means
    PTS, # a vector containing PTS values to simulate
    MTS, # a vector containing MTS values to simulate
    TDF_m, # trophic discrimination factors
    TDF_p, # trophic discrimination factors
    Sources, # names of organic matter sources
    Data.sources, # data for organic matter sources
    Tracers, # names of tracers
    Variables, # additional variables to include
    Random_Samples = TRUE, # should the model use Dirichlet to simulate zooplankton data?
    seed = 222 # option to set seed
){
  
  ## Use Dirichlet to generate normally distributed mixing parameters
  ## Do not use Dirichlet to explicitly define mixing parameters
  if (Random_Samples == TRUE) {
    ## the following lines are meant to simulate a handful of zooplankton samples
    ## via the Dirichlet distribution
    ## specify the average contribution of each organic matter source to the base 
    ## of the food web in the order they were originally described
    f_all <- rep(1, length(Sources))
    # increase the precision parameter to pull samples around the mean value 
    # decrease the precision parameter to push samples towards the edges of the distribution
    precision <- 0.5
    alpha   <- precision*f_all # shape parameters
    # now we use a Dirichlet distribution to produce a handful of zooplankton with
    # compositions distributed around the mean mixture composition
    set.seed(seed) # for repeatability
    # simulated f using dirichlet
    zoops.f <- data.frame(rdirichlet(n=nzoops, alpha=alpha)) 
    colnames(zoops.f) <- Sources # columns refer to fractional importance of each particle
    zoops.f$PTS <- runif(n=nzoops, min = min(PTS), max = max(PTS))
    zoops.f$MTS <- runif(n=nzoops, min = min(MTS), max = max(MTS))
    zoops.f$FWL <- zoops.f$PTS + zoops.f$MTS
    zoops.f <- zoops.f[order(zoops.f$FWL, zoops.f$PTS, zoops.f$MTS, # order by trophic level
                             -zoops.f[,1], -zoops.f[,2]), ] # then by the first two compositional columns
    row.names(zoops.f) <- seq(1,nzoops)
  } else {
    # OTHERWISE
    # Generate a larger number of representative samples by uniformly sampling compositional space
    zoops.f <-
      data.frame("Surface" = NA,
                 "Large" = NA,
                 "Small" = NA,
                 "PTS" = NA,
                 "MTS" = NA)
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
    nzoops <- nrow(zoops.f) # rows refer to samples
    zoops.f <- zoops.f[order(zoops.f$FWL, zoops.f$PTS, zoops.f$MTS, # order by trophic level
                             -zoops.f[,1], -zoops.f[,2]), ] # then by the first two compositional columns
    row.names(zoops.f) <- seq(1,nzoops)
  }
  
  # can visualize Dirichlet distribution by plotting these on a ternary diagram
  # must pool multiple sources per axis if more than 3 sources are present
  if(length(Sources) == 3){
    par(mar=rep(0, 4))
    TernaryPlot(alab=Sources[1], blab=Sources[2], clab=Sources[3])
    TernaryPoints(
      data.frame(A=zoops.f[Sources[1]], 
                 B=zoops.f[Sources[2]], 
                 C=zoops.f[Sources[3]]))
  }else{
    if(length(Sources) == 2){
      plot(
        x = zoops.f[[Sources[1]]],
        y = zoops.f[[Sources[2]]],
        xlab = Sources[1], 
        ylab = Sources[2])
    }
  }
  par(mar=rep(2, 4))
  hist(zoops.f$PTS, main = "Protistan Trophic Steps")
  hist(zoops.f$MTS, main = "Metazoan Trophic Steps")
  # print(datatable(round(zoops.f,2)))
  
  ## calculating the tracer values after mixing, before any trophic discrimination
  # calculating the mean tracer value for each organic matter source group
  src.mn <- 
    aggregate(Data.sources[Tracers$all], # aggregate source data
              by=list(Source = Data.sources[["Source"]]), # by organic matter source group
              FUN = mean, na.rm=TRUE)[-1] # taking a mean
  # we will calculate the standard deviation within each population
  src.SD <- 
    aggregate(Data.sources[Tracers$all], # aggregate source data
              by=list(Source = Data.sources[["Source"]]), # by organic matter source group
              FUN = sd, na.rm=TRUE)[-1] # calculating SD within the population
  
  # in real data the zooplankton delta values are often outside of the mean
  # values in organic matter sources, so we'll increase the variance of our
  # simulated samples by adding between group variance to the source data.
  src.mn <- src.mn + 
    sign((src.mn - matrix(rep(colMeans(src.mn),3), nrow = 3, byrow = TRUE)))*
    disperse * src.SD
  
  # now we'll calculate tracer values for each sample at the base of the food web
  base.sim <- data.frame(matrix(ncol = length(c("PTS","MTS","FWL",Sources,Tracers$all)), nrow = nzoops))
  colnames(base.sim) <- c("PTS","MTS","FWL",Sources,Tracers$all)
  base.sim[,Variables] <- seq(1,nzoops) # paste("sim ",seq(1,nzoops))
  base.sim[,c("PTS","MTS",Sources)] <- zoops.f[,c("PTS","MTS",Sources)]
  base.sim$FWL <- base.sim$PTS + base.sim$MTS
  for (i in 1:nzoops) {
    base.sim[i,Tracers$all] <- colSums(t(zoops.f[i,c(Sources)]) * src.mn)
  }
  
  ## Adding trophic discrimination
  
  # We'll define our some parameters we will need to account for trophic discrimination
  TDF_m <-    # trophic discrimination factors
    TDF_meta[c(Tracers$frac,"d15NPhe")]
  TDF_p <-    # trophic discrimination factors
    TDF_proto[c(Tracers$frac,"d15NPhe")]
  
  # initialize empty data frame to store simulated zoop data
  Data.zoops <- 
    data.frame(matrix(ncol = length(c(Variables,Tracers$all)), nrow = nzoops))
  colnames(Data.zoops) <- c(Variables,Tracers$all)
  # initialize empty data frame to store zooplankton data before adding process noise
  mn <- 
    data.frame(matrix(ncol = length(Tracers$all), nrow = 1))
  colnames(mn) <- Tracers$all
  
  ## loop to carry out trophic discrimination
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
    Data.zoops[iz,c(Tracers$frac,"d15NPhe")] <- 
      as.numeric(mn[1,c(Tracers$frac,"d15NPhe")])
  }
  Data.zoops[,Variables] <- seq(1,nzoops) # paste("sim ",seq(1,nzoops))
  Data.zoops[,SDTracers$all] <- 0.5
  
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
  nzoops <<- nzoops
  mn <<- mn
  MTS <<- MTS
  PTS <<- PTS
  FWL <<- MTS + PTS
  
}




#### Simulating zooplankton data! ####
# use this function if d15NPhe and d15NLys are being treated as a conservative tracer in the model
Sim_Zoop_RealPheLys <- function(
    nzoops = 21, # number of zooplankton to simulate
    incr = 0.2, # increment between zooplankton f values
    disperse = 0, # number of SDs by which to increase variance of zooplankton data beyond source means
    PTS, # a vector containing PTS values to simulate
    MTS, # a vector containing MTS values to simulate
    TDF_m, # trophic discrimination factors
    TDF_p, # trophic discrimination factors
    Sources, # names of organic matter sources
    Data.sources, # data for organic matter sources
    Tracers, # names of tracers
    Variables, # additional variables to include
    Random_Samples = TRUE, # should the model use Dirichlet to simulate zooplankton data?
    seed = 222 # option to set seed
){
  
  ## Use Dirichlet to generate normally distributed mixing parameters
  ## Do not use Dirichlet to explicitly define mixing parameters
  if (Random_Samples == TRUE) {
    ## the following lines are meant to simulate a handful of zooplankton samples
    ## via the Dirichlet distribution
    ## specify the average contribution of each organic matter source to the base 
    ## of the food web in the order they were originally described
    f_all <- rep(1, length(Sources))
    # increase the precision parameter to pull samples around the mean value 
    # decrease the precision parameter to push samples towards the edges of the distribution
    precision <- 0.5
    alpha   <- precision*f_all # shape parameters
    # now we use a Dirichlet distribution to produce a handful of zooplankton with
    # compositions distributed around the mean mixture composition
    set.seed(seed) # for repeatability
    # simulated f using dirichlet
    zoops.f <- data.frame(rdirichlet(n=nzoops, alpha=alpha)) 
    colnames(zoops.f) <- Sources # columns refer to fractional importance of each particle
    zoops.f$PTS <- runif(n=nzoops, min = min(PTS), max = max(PTS))
    zoops.f$MTS <- runif(n=nzoops, min = min(MTS), max = max(MTS))
    zoops.f$FWL <- zoops.f$PTS + zoops.f$MTS
    zoops.f <- zoops.f[order(zoops.f$FWL, zoops.f$PTS, zoops.f$MTS, # order by trophic level
                             -zoops.f[,1], -zoops.f[,2]), ] # then by the first two compositional columns
    row.names(zoops.f) <- seq(1,nzoops)
  } else {
    # OTHERWISE
    # Generate a larger number of representative samples by uniformly sampling compositional space
    zoops.f <-
      data.frame("Surface" = NA,
                 "Large" = NA,
                 "Small" = NA,
                 "PTS" = NA,
                 "MTS" = NA)
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
    nzoops <- nrow(zoops.f) # rows refer to samples
    zoops.f <- zoops.f[order(zoops.f$FWL, zoops.f$PTS, zoops.f$MTS, # order by trophic level
                             -zoops.f[,1], -zoops.f[,2]), ] # then by the first two compositional columns
    row.names(zoops.f) <- seq(1,nzoops)
  }
  
  # can visualize Dirichlet distribution by plotting these on a ternary diagram
  # must pool multiple sources per axis if more than 3 sources are present
  if(length(Sources) == 3){
    par(mar=rep(0, 4))
    TernaryPlot(alab=Sources[1], blab=Sources[2], clab=Sources[3])
    TernaryPoints(
      data.frame(A=zoops.f[Sources[1]], 
                 B=zoops.f[Sources[2]], 
                 C=zoops.f[Sources[3]]))
  }else{
    if(length(Sources) == 2){
      plot(
        x = zoops.f[[Sources[1]]],
        y = zoops.f[[Sources[2]]],
        xlab = Sources[1], 
        ylab = Sources[2])
    }
  }
  par(mar=rep(2, 4))
  hist(zoops.f$PTS, main = "Protistan Trophic Steps")
  hist(zoops.f$MTS, main = "Metazoan Trophic Steps")
  # print(datatable(round(zoops.f,2)))
  
  ## calculating the tracer values after mixing, before any trophic discrimination
  # calculating the mean tracer value for each organic matter source group
  src.mn <- 
    aggregate(Data.sources[Tracers$all], # aggregate source data
              by=list(Source = Data.sources[["Source"]]), # by organic matter source group
              FUN = mean, na.rm=TRUE)[-1] # taking a mean
  # we will calculate the standard deviation within each population
  src.SD <- 
    aggregate(Data.sources[Tracers$all], # aggregate source data
              by=list(Source = Data.sources[["Source"]]), # by organic matter source group
              FUN = sd, na.rm=TRUE)[-1] # calculating SD within the population
  
  # in real data the zooplankton delta values are often outside of the mean
  # values in organic matter sources, so we'll increase the variance of our
  # simulated samples by adding between group variance to the source data.
  src.mn <- src.mn + 
    sign((src.mn - matrix(rep(colMeans(src.mn),3), nrow = 3, byrow = TRUE)))*
    disperse * src.SD
  
  # now we'll calculate tracer values for each sample at the base of the food web
  base.sim <- data.frame(matrix(ncol = length(c("PTS","MTS","FWL",Sources,Tracers$all)), nrow = nzoops))
  colnames(base.sim) <- c("PTS","MTS","FWL",Sources,Tracers$all)
  base.sim[,Variables] <- seq(1,nzoops) # paste("sim ",seq(1,nzoops))
  base.sim[,c("PTS","MTS",Sources)] <- zoops.f[,c("PTS","MTS",Sources)]
  base.sim$FWL <- base.sim$PTS + base.sim$MTS
  for (i in 1:nzoops) {
    base.sim[i,Tracers$all] <- colSums(t(zoops.f[i,c(Sources)]) * src.mn)
  }
  
  ## Adding trophic discrimination
  
  # We'll define our some parameters we will need to account for trophic discrimination
  TDF_m <-    # trophic discrimination factors
    TDF_meta[c(Tracers$frac,"d15NPhe","d15NLys")]
  TDF_p <-    # trophic discrimination factors
    TDF_proto[c(Tracers$frac,"d15NPhe","d15NLys")]
  
  # initialize empty data frame to store simulated zoop data
  Data.zoops <- 
    data.frame(matrix(ncol = length(c(Variables,Tracers$all)), nrow = nzoops))
  colnames(Data.zoops) <- c(Variables,Tracers$all)
  # initialize empty data frame to store zooplankton data before adding process noise
  mn <- 
    data.frame(matrix(ncol = length(Tracers$all), nrow = 1))
  colnames(mn) <- Tracers$all
  
  ## loop to carry out trophic discrimination
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
    Data.zoops[iz,c(Tracers$frac,"d15NPhe","d15NLys")] <- 
      as.numeric(mn[1,c(Tracers$frac,"d15NPhe","d15NLys")])
  }
  Data.zoops[,Variables] <- seq(1,nzoops) # paste("sim ",seq(1,nzoops))
  Data.zoops[,SDTracers$all] <- 0.5
  
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
  nzoops <<- nzoops
  mn <<- mn
  MTS <<- MTS
  PTS <<- PTS
  FWL <<- MTS + PTS
  
}




#### Simulating zooplankton data! ####
# use this function if all SAAs are being treated as a conservative tracer in the model
Sim_Zoop_RealSAAs <- function(
    nzoops = 21, # number of zooplankton to simulate
    incr = 0.2, # increment between zooplankton f values
    disperse = 0, # number of SDs by which to increase variance of zooplankton data beyond source means
    PTS, # a vector containing PTS values to simulate
    MTS, # a vector containing MTS values to simulate
    TDF_m, # trophic discrimination factors
    TDF_p, # trophic discrimination factors
    Sources, # names of organic matter sources
    Data.sources, # data for organic matter sources
    Tracers, # names of tracers
    Variables, # additional variables to include
    Random_Samples = TRUE, # should the model use Dirichlet to simulate zooplankton data?
    seed = 222 # option to set seed
){
  
  ## Use Dirichlet to generate normally distributed mixing parameters
  ## Do not use Dirichlet to explicitly define mixing parameters
  if (Random_Samples == TRUE) {
    ## the following lines are meant to simulate a handful of zooplankton samples
    ## via the Dirichlet distribution
    ## specify the average contribution of each organic matter source to the base 
    ## of the food web in the order they were originally described
    f_all <- rep(1, length(Sources))
    # increase the precision parameter to pull samples around the mean value 
    # decrease the precision parameter to push samples towards the edges of the distribution
    precision <- 0.5
    alpha   <- precision*f_all # shape parameters
    # now we use a Dirichlet distribution to produce a handful of zooplankton with
    # compositions distributed around the mean mixture composition
    set.seed(seed) # for repeatability
    # simulated f using dirichlet
    zoops.f <- data.frame(rdirichlet(n=nzoops, alpha=alpha)) 
    colnames(zoops.f) <- Sources # columns refer to fractional importance of each particle
    zoops.f$PTS <- runif(n=nzoops, min = min(PTS), max = max(PTS))
    zoops.f$MTS <- runif(n=nzoops, min = min(MTS), max = max(MTS))
    zoops.f$FWL <- zoops.f$PTS + zoops.f$MTS
    zoops.f <- zoops.f[order(zoops.f$FWL, zoops.f$PTS, zoops.f$MTS, # order by trophic level
                             -zoops.f[,1], -zoops.f[,2]), ] # then by the first two compositional columns
    row.names(zoops.f) <- seq(1,nzoops)
  } else {
    # OTHERWISE
    # Generate a larger number of representative samples by uniformly sampling compositional space
    zoops.f <-
      data.frame("Surface" = NA,
                 "Large" = NA,
                 "Small" = NA,
                 "PTS" = NA,
                 "MTS" = NA)
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
    nzoops <- nrow(zoops.f) # rows refer to samples
    zoops.f <- zoops.f[order(zoops.f$FWL, zoops.f$PTS, zoops.f$MTS, # order by trophic level
                             -zoops.f[,1], -zoops.f[,2]), ] # then by the first two compositional columns
    row.names(zoops.f) <- seq(1,nzoops)
  }
  
  # can visualize Dirichlet distribution by plotting these on a ternary diagram
  # must pool multiple sources per axis if more than 3 sources are present
  if(length(Sources) == 3){
    par(mar=rep(0, 4))
    TernaryPlot(alab=Sources[1], blab=Sources[2], clab=Sources[3])
    TernaryPoints(
      data.frame(A=zoops.f[Sources[1]], 
                 B=zoops.f[Sources[2]], 
                 C=zoops.f[Sources[3]]))
  }else{
    if(length(Sources) == 2){
      plot(
        x = zoops.f[[Sources[1]]],
        y = zoops.f[[Sources[2]]],
        xlab = Sources[1], 
        ylab = Sources[2])
    }
  }
  par(mar=rep(2, 4))
  hist(zoops.f$PTS, main = "Protistan Trophic Steps")
  hist(zoops.f$MTS, main = "Metazoan Trophic Steps")
  # print(datatable(round(zoops.f,2)))
  
  ## calculating the tracer values after mixing, before any trophic discrimination
  # calculating the mean tracer value for each organic matter source group
  src.mn <- 
    aggregate(Data.sources[Tracers$all], # aggregate source data
              by=list(Source = Data.sources[["Source"]]), # by organic matter source group
              FUN = mean, na.rm=TRUE)[-1] # taking a mean
  # we will calculate the standard deviation within each population
  src.SD <- 
    aggregate(Data.sources[Tracers$all], # aggregate source data
              by=list(Source = Data.sources[["Source"]]), # by organic matter source group
              FUN = sd, na.rm=TRUE)[-1] # calculating SD within the population
  
  # in real data the zooplankton delta values are often outside of the mean
  # values in organic matter sources, so we'll increase the variance of our
  # simulated samples by adding between group variance to the source data.
  src.mn <- src.mn + 
    sign((src.mn - matrix(rep(colMeans(src.mn),3), nrow = 3, byrow = TRUE)))*
    disperse * src.SD
  
  # now we'll calculate tracer values for each sample at the base of the food web
  base.sim <- data.frame(matrix(ncol = length(c("PTS","MTS","FWL",Sources,Tracers$all)), nrow = nzoops))
  colnames(base.sim) <- c("PTS","MTS","FWL",Sources,Tracers$all)
  base.sim[,Variables] <- seq(1,nzoops) # paste("sim ",seq(1,nzoops))
  base.sim[,c("PTS","MTS",Sources)] <- zoops.f[,c("PTS","MTS",Sources)]
  base.sim$FWL <- base.sim$PTS + base.sim$MTS
  for (i in 1:nzoops) {
    base.sim[i,Tracers$all] <- colSums(t(zoops.f[i,c(Sources)]) * src.mn)
  }
  
  ## Adding trophic discrimination
  
  # We'll define our some parameters we will need to account for trophic discrimination
  TDF_m <-    # trophic discrimination factors
    TDF_meta[c(Tracers$frac,"d15NPhe","d15NLys","d15NGly","d15NSer")]
  TDF_p <-    # trophic discrimination factors
    TDF_proto[c(Tracers$frac,"d15NPhe","d15NLys","d15NGly","d15NSer")]
  
  # initialize empty data frame to store simulated zoop data
  Data.zoops <- 
    data.frame(matrix(ncol = length(c(Variables,Tracers$all)), nrow = nzoops))
  colnames(Data.zoops) <- c(Variables,Tracers$all)
  # initialize empty data frame to store zooplankton data before adding process noise
  mn <- 
    data.frame(matrix(ncol = length(Tracers$all), nrow = 1))
  colnames(mn) <- Tracers$all
  
  ## loop to carry out trophic discrimination
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
    Data.zoops[iz,c(Tracers$frac,"d15NPhe","d15NLys","d15NGly","d15NSer")] <- 
      as.numeric(mn[1,c(Tracers$frac,"d15NPhe","d15NLys","d15NGly","d15NSer")])
  }
  Data.zoops[,Variables] <- seq(1,nzoops) # paste("sim ",seq(1,nzoops))
  Data.zoops[,SDTracers$all] <- 0.5
  
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
  nzoops <<- nzoops
  mn <<- mn
  MTS <<- MTS
  PTS <<- PTS
  FWL <<- MTS + PTS
  
}