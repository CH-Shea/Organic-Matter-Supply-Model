## Simulating zooplankton data! ##
Sim_Zoop <- function(
    zoops.n,
    PTS,
    MTS,
    PTS.sd,
    MTS.sd,
    TDF_m, # trophic discrimination factors
    TDF_m.sd, # SD of trophic discrimination factors
    TDF_p, # trophic discrimination factors
    TDF_p.sd, # SD of trophic discrimination factors
    Sources,
    Data.sources,
    Tracers,
    SDTracers,
    variables,
    Use_Dirichlet,
    zoops.f,
    seed
){

# First, even is we are assuming conservative behavior in Phe, we want to simulate data that reflects fractionation in Phe
if(sum(which(Tracers$non == "d15NPhe")) > 0) {
  Tracers$constTDF <- c(Tracers$constTDF, "d15NPhe")
  SDTracers$constTDF <- c(SDTracers$constTDF, "SDd15NPhe")
  Tracers$non <- Tracers$non[-which(Tracers$non == "d15NPhe")]
  SDTracers$non <- SDTracers$non[-which(SDTracers$non == "SDd15NPhe")]
}
  
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
  zoops.n <- 10 # n zoops we'll simulate
  set.seed(seed) # for repeatability
  zoops.f <- data.frame(rdirichlet(n=zoops.n, alpha=alpha)) # simulated f
  colnames(zoops.f) <- Sources # columns refer to fractional importance of each particle
  
} else {
  ## The following lines allow the user to specify the relative contribution of each organic matter source to each zooplankton sample.
  # enter the % contribution of each organic matter source to each zooplankton sample graphically below: 
  # the columns refer to organic matter sources in the order they were originally specified. 
  # Each row is a zooplankton sample
  zoops.f <- data.frame(
    t(matrix(nrow = length(Sources),
             c(98,01,01,
               01,98,01,
               01,01,98,
               49,49,04,
               49,02,49,
               02,49,49,
               33,33,34)*0.01))
  )
  colnames(zoops.f) <- Sources # columns refer to fractional importance of each particle
  zoops.n <- nrow(zoops.f) # rows refer to samples
}

# can visualize Dirichlet distribution by plotting these on a ternary diagram
# must pool multiple sources per axis if more than 3 sources are present
par(mar=rep(0, 4))
TernaryPlot(alab="f(Small)", blab="f(Large)", clab="f(Surface)")
TernaryPoints(
  data.frame(A=zoops.f$Small, 
             B=zoops.f$Large, 
             C=zoops.f$Surface))

# calculating the tracer values after mixing, before any trophic discrimination
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
base.sim <- data.frame(matrix(ncol = length(c(variables,Tracers$all)), nrow = zoops.n))
colnames(base.sim) <- c(variables,Tracers$all)
for (i in 1:zoops.n) {
  base.sim[i,Tracers$all] <- colSums(t(zoops.f[i,]) * src.mn)
}
base.sim[,variables] <- seq(1,zoops.n) # paste("sim ",seq(1,zoops.n))

## Adding trophic discrimination

# We'll define our some parameters we will need to account for trophic
# discrimination and propagate error through
FWL <-    # total food web length
  PTS + MTS
FWL.sd <- 0.5
TDF_m <-    # trophic discrimination factors
  TDF_meta[Tracers$frac]
TDF_m.sd <- # SD of trophic discrimination factors
  TDF_meta[SDTracers$frac]
TDF_p <-    # trophic discrimination factors
  TDF_proto[Tracers$frac]
TDF_p.sd <- # SD of trophic discrimination factors
  TDF_proto[SDTracers$frac]

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
for (iz in 1:zoops.n) { # zooplankton d15N loop
  # setting mn_z and mn_base equal for conservative tracers
  if(length(Tracers$non) > 0) {
    mn[c(Tracers$non)] <- base.sim[iz,c(Tracers$non)]
  }
  # adding trophic discrimination for amino acids with constant TDFs
  if(length(Tracers$constTDF) > 0) {
    mn[c(Tracers$constTDF)] <- base.sim[iz,c(Tracers$constTDF)] + 
      TDF_m[c(Tracers$constTDF)]*FWL
  }
  # adding trophic discrimination for amino acids with variable TDFs
  if(length(Tracers$varTDF) > 0) {
    mn[Tracers$varTDF] <- base.sim[iz,Tracers$varTDF] + 
      TDF_p[Tracers$varTDF]*PTS +
      TDF_m[Tracers$varTDF]*MTS
  }
  # calculating the variance and sd by propagating error through mixing and
  # trophic discrimination
  var[c(Tracers$non)] <-
    sum(zoops.f[iz,]^2*sd.anl^2)^2
  var[c(Tracers$constTDF)] <- 
    sqrt(
      sum(zoops.f[iz,]^2*sd.anl^2)^2 +
        FWL^2*TDF_m.sd[c(SDTracers$constTDF)]^2 +
        TDF_m[c(Tracers$constTDF)]^2*FWL.sd^2
    )
  var[Tracers$varTDF] <- 
    sqrt(
      sum(zoops.f[iz,]^2*sd.anl^2)^2 +
        PTS^2*TDF_p.sd[SDTracers$varTDF]^2 +
        TDF_p[Tracers$varTDF]^2*PTS.sd^2 +
        MTS^2*TDF_m.sd[SDTracers$varTDF]^2 +
        TDF_m[Tracers$varTDF]^2*MTS.sd^2
    )
  # and drawing the zooplankton d15N from a normal with that mn and sd
  for (j in c(Tracers$all)) {
    Data.zoops[iz,j] <- rnorm(1, as.numeric(mn[1,j]), 0.01)
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
zoops.f <<- zoops.f
zoops.n <<- zoops.n
mn <<- mn
variance <<- var
MTS <<- MTS
PTS <<- PTS
FWL <<- FWL
}