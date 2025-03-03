OMSM_Zoop_Extract <- function(
    rjo_1,
    Sources,
    Data.zoops,
    Tracers.all,
    Tracers.frac,
    Tracers.varTDF
){
  
  ## Grabbing MCMC model output
  sams <<- as.matrix(rjo_1$mcmc)
  nsams <<- nrow(sams)
  nsources <<- length(Sources)
  nzoops <<- nrow(Data.zoops)
  ntracers <<- length(Tracers.all)
  nTracers.frac <<- length(Tracers.frac)
  nTracers.varTDF <<- length(Tracers.varTDF)
  
## 2. Zooplankton food web parameter posteriors (mixing coefficients, trophic positions, and δ15N values at the base of the food web) ##
sams <- as.matrix(rjo_1$mcmc)
# we'll make a list which contains samples and summary data for each posterior
posts.zoops <- list(
  f = list(
    samples = data.frame(matrix(nrow = nsams,ncol=0)),
    HDI95 = data.frame(matrix(nrow = 2,ncol=0)),
    HDI90 = data.frame(matrix(nrow = 2,ncol=0)),
    HDI75 = data.frame(matrix(nrow = 2,ncol=0)),
    HDI50 = data.frame(matrix(nrow = 2,ncol=0)),
    mean = data.frame(matrix(nrow = 1,ncol=0)),
    mode = data.frame(matrix(nrow = 1,ncol=0))
  ),
  trophic = list(
    samples = data.frame(matrix(nrow = nsams,ncol=0)),
    HDI95 = data.frame(matrix(nrow = 2,ncol=0)),
    HDI90 = data.frame(matrix(nrow = 2,ncol=0)),
    HDI75 = data.frame(matrix(nrow = 2,ncol=0)),
    HDI50 = data.frame(matrix(nrow = 2,ncol=0)),
    mean = data.frame(matrix(nrow = 1,ncol=0)),
    mode = data.frame(matrix(nrow = 1,ncol=0))
  ),
  base = list(
    samples = data.frame(matrix(nrow = nsams*nsources,ncol=0)),
    HDI95 = data.frame(matrix(nrow = nsources*2,ncol=0)),
    HDI90 = data.frame(matrix(nrow = nsources*2,ncol=0)),
    HDI75 = data.frame(matrix(nrow = nsources*2,ncol=0)),
    HDI50 = data.frame(matrix(nrow = nsources*2,ncol=0)),
    mean = data.frame(matrix(nrow = nsources,ncol=0)),
    mode = data.frame(matrix(nrow = nsources,ncol=0))
  ),
  zoop = list(
    samples = data.frame(matrix(nrow = nsams*nzoops,ncol=0)),
    HDI95 = data.frame(matrix(nrow = nzoops*2,ncol=0)),
    HDI90 = data.frame(matrix(nrow = nzoops*2,ncol=0)),
    HDI75 = data.frame(matrix(nrow = nzoops*2,ncol=0)),
    HDI50 = data.frame(matrix(nrow = nzoops*2,ncol=0)),
    mean = data.frame(matrix(nrow = nzoops,ncol=0)),
    mode = data.frame(matrix(nrow = nzoops,ncol=0))
  )
)

# grabbing f(source) posteriors
temp1<-data.frame(matrix(ncol=length(variables),nrow=nzoops*nsams))
temp2<-data.frame(matrix(ncol=length(variables),nrow=nzoops*2))
temp3<-data.frame(matrix(ncol=length(variables),nrow=nzoops))
colnames(temp1)<-colnames(temp2)<-colnames(temp3)<-variables
# initial columns give info on each zooplankton sample
# Data.zoops$Size <- factor(Data.zoops$Size, levels = c("0.2-0.5 mm","0.5-1.0 mm","1-2 mm","2-5 mm",">5"))
for (j in 1:nzoops) {
  temp1[((j-1)*nsams+1):(nsams*j),] <- Data.zoops[j,variables]
  temp2[((j-1)*2+1):(2*j),] <- Data.zoops[j,variables]
  temp3[j,] <- Data.zoops[j,variables]
}
temp1$Size <- factor(temp1$Size, levels = c("0.2-0.5 mm","0.5-1.0 mm","1-2 mm","2-5 mm",">5 mm"))
temp2$Size <- factor(temp2$Size, levels = c("0.2-0.5 mm","0.5-1.0 mm","1-2 mm","2-5 mm",">5 mm"))
temp3$Size <- factor(temp3$Size, levels = c("0.2-0.5 mm","0.5-1.0 mm","1-2 mm","2-5 mm",">5 mm"))
posts.zoops$f$samples <- posts.zoops$trophic$samples <- posts.zoops$base$samples <- posts.zoops$zoop$samples <- temp1
posts.zoops$f$HDI95 <- posts.zoops$trophic$HDI95 <- posts.zoops$base$HDI95 <- posts.zoops$zoop$HDI95 <-
  posts.zoops$f$HDI90 <- posts.zoops$trophic$HDI90 <- posts.zoops$base$HDI90 <- posts.zoops$zoop$HDI90 <-
  posts.zoops$f$HDI75 <- posts.zoops$trophic$HDI75 <- posts.zoops$base$HDI75 <- posts.zoops$zoop$HDI75 <-
  posts.zoops$f$HDI50 <- posts.zoops$trophic$HDI50 <- posts.zoops$base$HDI50 <- posts.zoops$zoop$HDI50 <- temp2
posts.zoops$f$mean <- posts.zoops$trophic$mean <- posts.zoops$base$mean <- posts.zoops$zoop$mean <-
  posts.zoops$f$mode <- posts.zoops$trophic$mode <- posts.zoops$base$mode <- posts.zoops$zoop$mode <- temp3

# next columns give respective f data for each source
for (j in 1:nsources) {
  temp1<-temp2<-temp3<-temp4<-temp5<-temp6<-temp7<-c()
  for (i in 1:nzoops) {
    temp1 <- c(temp1, sams[,paste("pz[",i,",",j,"]",sep = "")])
    temp2 <- c(temp2, HDIofMCMC(sams[,paste("pz[",i,",",j,"]",sep = "")], credMass=0.95))
    temp3 <- c(temp3, HDIofMCMC(sams[,paste("pz[",i,",",j,"]",sep = "")], credMass=0.90))
    temp4 <- c(temp4, HDIofMCMC(sams[,paste("pz[",i,",",j,"]",sep = "")], credMass=0.75))
    temp5 <- c(temp5, HDIofMCMC(sams[,paste("pz[",i,",",j,"]",sep = "")], credMass=0.50))
    temp6 <- c(temp6, mean(sams[,paste("pz[",i,",",j,"]",sep = "")]))
    temp7 <- c(temp7, post.mode(sams[,paste("pz[",i,",",j,"]",sep = "")]))
  }
  posts.zoops$f$samples[Sources[j]] <- temp1
  posts.zoops$f$HDI95[Sources[j]] <- temp2
  posts.zoops$f$HDI90[Sources[j]] <- temp3
  posts.zoops$f$HDI75[Sources[j]] <- temp4
  posts.zoops$f$HDI50[Sources[j]] <- temp5
  posts.zoops$f$mean[Sources[j]] <- temp6
  posts.zoops$f$mode[Sources[j]] <- temp7
}
# next columns give respective trophic parameters
trophs <- c("FWL","PTS","MTS")
for (j in 1:length(trophs)) {
  temp1<-temp2<-temp3<-temp4<-temp5<-temp6<-temp7<-c()
  for (i in 1:nzoops) {
    temp1 <- c(temp1, sams[,paste(trophs[j],"[",i,"]",sep = "")])
    temp2 <- c(temp2, HDIofMCMC(sams[,paste(trophs[j],"[",i,"]",sep = "")], credMass=0.95))
    temp3 <- c(temp3, HDIofMCMC(sams[,paste(trophs[j],"[",i,"]",sep = "")], credMass=0.90))
    temp4 <- c(temp4, HDIofMCMC(sams[,paste(trophs[j],"[",i,"]",sep = "")], credMass=0.75))
    temp5 <- c(temp5, HDIofMCMC(sams[,paste(trophs[j],"[",i,"]",sep = "")], credMass=0.50))
    temp6 <- c(temp6, mean(sams[,paste(trophs[j],"[",i,"]",sep = "")]))
    temp7 <- c(temp7, post.mode(sams[,paste(trophs[j],"[",i,"]",sep = "")]))
  }
  posts.zoops$trophic$samples[trophs[j]] <- temp1
  posts.zoops$trophic$HDI95[trophs[j]] <- temp2
  posts.zoops$trophic$HDI90[trophs[j]] <- temp3
  posts.zoops$trophic$HDI75[trophs[j]] <- temp4
  posts.zoops$trophic$HDI50[trophs[j]] <- temp5
  posts.zoops$trophic$mean[trophs[j]] <- temp6
  posts.zoops$trophic$mode[trophs[j]] <- temp7
}
# next columns give respective δ15N_base values
for (j in 1:ntracers) {
  temp1<-temp2<-temp3<-temp4<-temp5<-temp6<-temp7<-c()
  for (i in 1:nzoops) {
    temp1 <- c(temp1, sams[,paste("mean_b[",i,",",j,"]",sep = "")])
    temp2 <- c(temp2, HDIofMCMC(sams[,paste("mean_b[",i,",",j,"]",sep = "")], credMass=0.95))
    temp3 <- c(temp3, HDIofMCMC(sams[,paste("mean_b[",i,",",j,"]",sep = "")], credMass=0.90))
    temp4 <- c(temp4, HDIofMCMC(sams[,paste("mean_b[",i,",",j,"]",sep = "")], credMass=0.75))
    temp5 <- c(temp5, HDIofMCMC(sams[,paste("mean_b[",i,",",j,"]",sep = "")], credMass=0.50))
    temp6 <- c(temp6, mean(sams[,paste("mean_b[",i,",",j,"]",sep = "")]))
    temp7 <- c(temp7, post.mode(sams[,paste("mean_b[",i,",",j,"]",sep = "")]))
  }
  posts.zoops$base$samples[Tracers.all[j]] <- temp1
  posts.zoops$base$HDI95[Tracers.all[j]] <- temp2
  posts.zoops$base$HDI90[Tracers.all[j]] <- temp3
  posts.zoops$base$HDI75[Tracers.all[j]] <- temp4
  posts.zoops$base$HDI50[Tracers.all[j]] <- temp5
  posts.zoops$base$mean[Tracers.all[j]] <- temp6
  posts.zoops$base$mode[Tracers.all[j]] <- temp7
}
# next columns give respective δ15N_zoop values
for (j in 1:ntracers) {
  temp1<-temp2<-temp3<-temp4<-temp5<-temp6<-temp7<-c()
  for (i in 1:nzoops) {
    temp1 <- c(temp1, sams[,paste("mean_z[",i,",",j,"]",sep = "")])
    temp2 <- c(temp2, HDIofMCMC(sams[,paste("mean_z[",i,",",j,"]",sep = "")], credMass=0.95))
    temp3 <- c(temp3, HDIofMCMC(sams[,paste("mean_z[",i,",",j,"]",sep = "")], credMass=0.90))
    temp4 <- c(temp4, HDIofMCMC(sams[,paste("mean_z[",i,",",j,"]",sep = "")], credMass=0.75))
    temp5 <- c(temp5, HDIofMCMC(sams[,paste("mean_z[",i,",",j,"]",sep = "")], credMass=0.50))
    temp6 <- c(temp6, mean(sams[,paste("mean_z[",i,",",j,"]",sep = "")]))
    temp7 <- c(temp7, post.mode(sams[,paste("mean_z[",i,",",j,"]",sep = "")]))
  }
  posts.zoops$zoop$samples[Tracers.all[j]] <- temp1
  posts.zoops$zoop$HDI95[Tracers.all[j]] <- temp2
  posts.zoops$zoop$HDI90[Tracers.all[j]] <- temp3
  posts.zoops$zoop$HDI75[Tracers.all[j]] <- temp4
  posts.zoops$zoop$HDI50[Tracers.all[j]] <- temp5
  posts.zoops$zoop$mean[Tracers.all[j]] <- temp6
  posts.zoops$zoop$mode[Tracers.all[j]] <- temp7
}

# at times we will want thinned samlpes for plotting purposes
thin.by <- round(samsPerChain/100,0) # will select 1 out of every n samples
thin.f <- seq(1,nrow(posts.zoops$f$samples),thin.by)
posts.zoops$f$thin <- posts.zoops$f$samples[thin.f,]
thin.trophic <- seq(1,nrow(posts.zoops$trophic$samples),thin.by)
posts.zoops$trophic$thin <- posts.zoops$trophic$samples[thin.trophic,]
thin.base <- seq(1,nrow(posts.zoops$base$samples),thin.by)
posts.zoops$base$thin <- posts.zoops$base$samples[thin.base,]
thin.zoop <- seq(1,nrow(posts.zoops$zoop$samples),thin.by)
posts.zoops$zoop$thin <- posts.zoops$zoop$samples[thin.zoop,]

# at times we'll want this data in long format so we'll handle that here as well
posts.zoops.long <- list(
  f = list(
    samples = melt(posts.zoops$f$samples, 
                   id.vars = variables, value.name = "f", variable.name = "Source"),
    thin = melt(posts.zoops$f$thin, 
                id.vars = variables, value.name = "f", variable.name = "Source"),
    HDI95 = melt(posts.zoops$f$HDI95, 
                 id.vars = variables, value.name = "f", variable.name = "Source"),
    HDI90 = melt(posts.zoops$f$HDI90, 
                 id.vars = variables, value.name = "f", variable.name = "Source"),
    HDI75 = melt(posts.zoops$f$HDI75, 
                 id.vars = variables, value.name = "f", variable.name = "Source"),
    HDI50 = melt(posts.zoops$f$HDI50, 
                 id.vars = variables, value.name = "f", variable.name = "Source"),
    mean = melt(posts.zoops$f$mean, 
                id.vars = variables, value.name = "f", variable.name = "Source"),
    mode = melt(posts.zoops$f$mode, 
                id.vars = variables, value.name = "f", variable.name = "Source")
  ),
  trophic = list(
    samples = melt(posts.zoops$trophic$samples, 
                   id.vars = variables, value.name = "Value", variable.name = "Param"),
    thin = melt(posts.zoops$trophic$thin, 
                id.vars = variables, value.name = "Value", variable.name = "Param"),
    HDI95 = melt(posts.zoops$trophic$HDI95, 
                 id.vars = variables, value.name = "Value", variable.name = "Param"),
    HDI90 = melt(posts.zoops$trophic$HDI90, 
                 id.vars = variables, value.name = "Value", variable.name = "Param"),
    HDI75 = melt(posts.zoops$trophic$HDI75, 
                 id.vars = variables, value.name = "Value", variable.name = "Param"),
    HDI50 = melt(posts.zoops$trophic$HDI50, 
                 id.vars = variables, value.name = "Value", variable.name = "Param"),
    mean = melt(posts.zoops$trophic$mean, 
                id.vars = variables, value.name = "Value", variable.name = "Param"),
    mode = melt(posts.zoops$trophic$mode, 
                id.vars = variables, value.name = "Value", variable.name = "Param")
  ),
  base = list(
    samples = melt(posts.zoops$base$samples, 
                   id.vars = variables, value.name = "Value", variable.name = "Tracer"),
    thin = melt(posts.zoops$base$thin, 
                id.vars = variables, value.name = "Value", variable.name = "Tracer"),
    HDI95 = melt(posts.zoops$base$HDI95, 
                 id.vars = variables, value.name = "Value", variable.name = "Tracer"),
    HDI90 = melt(posts.zoops$base$HDI90, 
                 id.vars = variables, value.name = "Value", variable.name = "Tracer"),
    HDI75 = melt(posts.zoops$base$HDI75, 
                 id.vars = variables, value.name = "Value", variable.name = "Tracer"),
    HDI50 = melt(posts.zoops$base$HDI50, 
                 id.vars = variables, value.name = "Value", variable.name = "Tracer"),
    mean = melt(posts.zoops$base$mean, 
                id.vars = variables, value.name = "Value", variable.name = "Tracer"),
    mode = melt(posts.zoops$base$mode, 
                id.vars = variables, value.name = "Value", variable.name = "Tracer")
  ),
  zoop = list(
    samples = melt(posts.zoops$zoop$samples, 
                   id.vars = variables, value.name = "Value", variable.name = "Tracer"),
    thin = melt(posts.zoops$zoop$thin, 
                id.vars = variables, value.name = "Value", variable.name = "Tracer"),
    HDI95 = melt(posts.zoops$zoop$HDI95, 
                 id.vars = variables, value.name = "Value", variable.name = "Tracer"),
    HDI90 = melt(posts.zoops$zoop$HDI90, 
                 id.vars = variables, value.name = "Value", variable.name = "Tracer"),
    HDI75 = melt(posts.zoops$zoop$HDI75, 
                 id.vars = variables, value.name = "Value", variable.name = "Tracer"),
    HDI50 = melt(posts.zoops$zoop$HDI50, 
                 id.vars = variables, value.name = "Value", variable.name = "Tracer"),
    mean = melt(posts.zoops$zoop$mean, 
                id.vars = variables, value.name = "Value", variable.name = "Tracer"),
    mode = melt(posts.zoops$zoop$mode, 
                id.vars = variables, value.name = "Value", variable.name = "Tracer")
  )
)

# # last, we will at times want to pool two mixing parameters to fit on a ternary diagram
# # Large + Small particles
# temp1<-temp2<-temp3<-temp4<-temp5<-temp6<-temp7<- c()
# for (i in 1:nzoops) {
#     temp1 <- c(temp1, sams[,paste("pz[",i,",2]",sep = "")]+sams[,paste("pz[",i,",3]",sep = "")])
#     temp2 <- 
#       c(temp2, HDIofMCMC(sams[,paste("pz[",i,",2]",sep = "")]+
#                            sams[,paste("pz[",i,",3]",sep = "")], credMass=0.95))
#     temp3 <- 
#       c(temp3, HDIofMCMC(sams[,paste("pz[",i,",2]",sep = "")]+
#                            sams[,paste("pz[",i,",3]",sep = "")], credMass=0.90))
#     temp4 <- 
#       c(temp4, HDIofMCMC(sams[,paste("pz[",i,",2]",sep = "")]+
#                            sams[,paste("pz[",i,",3]",sep = "")], credMass=0.75))
#     temp5 <- 
#       c(temp5, HDIofMCMC(sams[,paste("pz[",i,",2]",sep = "")]+
#                            sams[,paste("pz[",i,",3]",sep = "")], credMass=0.50))
#     temp6 <- 
#       c(temp6, mean(sams[,paste("pz[",i,",2]",sep = "")]+
#                       sams[,paste("pz[",i,",3]",sep = "")]))
#     temp7 <- 
#       c(temp7, post.mode(sams[,paste("pz[",i,",2]",sep = "")]+
#                            sams[,paste("pz[",i,",3]",sep = "")]))
#   }
#   posts.zoops$f$samples$LargeSmall <- temp1
#   posts.zoops$f$HDI95$LargeSmall <- temp2
#   posts.zoops$f$HDI90$LargeSmall <- temp3
#   posts.zoops$f$HDI75$LargeSmall <- temp4
#   posts.zoops$f$HDI50$LargeSmall <- temp5
#   posts.zoops$f$mean$LargeSmall <- temp6
#   posts.zoops$f$mode$LargeSmall <- temp7
#   
# # Small + Submicron particles
# temp1<-temp2<-temp3<-temp4<-temp5<-temp6<-temp7<- c()
# for (i in 1:nzoops) {
#     temp1 <- c(temp1, sams[,paste("pz[",i,",4]",sep = "")]+
#                  sams[,paste("pz[",i,",3]",sep = "")])
#     temp2 <- 
#       c(temp2, HDIofMCMC(sams[,paste("pz[",i,",4]",sep = "")]+
#                            sams[,paste("pz[",i,",3]",sep = "")], credMass=0.95))
#     temp3 <- 
#       c(temp3, HDIofMCMC(sams[,paste("pz[",i,",4]",sep = "")]+
#                            sams[,paste("pz[",i,",3]",sep = "")], credMass=0.90))
#     temp4 <- 
#       c(temp4, HDIofMCMC(sams[,paste("pz[",i,",4]",sep = "")]+
#                            sams[,paste("pz[",i,",3]",sep = "")], credMass=0.75))
#     temp5 <- 
#       c(temp5, HDIofMCMC(sams[,paste("pz[",i,",4]",sep = "")]+
#                            sams[,paste("pz[",i,",3]",sep = "")], credMass=0.50))
#     temp6 <- 
#       c(temp6, mean(sams[,paste("pz[",i,",4]",sep = "")]+
#                       sams[,paste("pz[",i,",3]",sep = "")]))
#     temp7 <- 
#       c(temp7, post.mode(sams[,paste("pz[",i,",4]",sep = "")]+
#                            sams[,paste("pz[",i,",3]",sep = "")]))
#   }
#   posts.zoops$f$samples$SmallSub <- temp1
#   posts.zoops$f$HDI95$SmallSub <- temp2
#   posts.zoops$f$HDI90$SmallSub <- temp3
#   posts.zoops$f$HDI75$SmallSub <- temp4
#   posts.zoops$f$HDI50$SmallSub <- temp5
#   posts.zoops$f$mean$SmallSub <- temp6
#   posts.zoops$f$mode$SmallSub <- temp7

# # Large + Surface particles
# temp1<-temp2<-temp3<-temp4<-temp5<-temp6<-temp7<- c()
# for (i in 1:nzoops) {
#     temp1 <- c(temp1, sams[,paste("pz[",i,",1]",sep = "")]+
#                  sams[,paste("pz[",i,",2]",sep = "")])
#     temp2 <- 
#       c(temp2, HDIofMCMC(sams[,paste("pz[",i,",1]",sep = "")]+
#                            sams[,paste("pz[",i,",2]",sep = "")], credMass=0.95))
#     temp3 <- 
#       c(temp3, HDIofMCMC(sams[,paste("pz[",i,",1]",sep = "")]+
#                            sams[,paste("pz[",i,",2]",sep = "")], credMass=0.90))
#     temp4 <- 
#       c(temp4, HDIofMCMC(sams[,paste("pz[",i,",1]",sep = "")]+
#                            sams[,paste("pz[",i,",2]",sep = "")], credMass=0.75))
#     temp5 <- 
#       c(temp5, HDIofMCMC(sams[,paste("pz[",i,",1]",sep = "")]+
#                            sams[,paste("pz[",i,",2]",sep = "")], credMass=0.50))
#     temp6 <- 
#       c(temp6, mean(sams[,paste("pz[",i,",1]",sep = "")]+
#                       sams[,paste("pz[",i,",2]",sep = "")]))
#     temp7 <- 
#       c(temp7, post.mode(sams[,paste("pz[",i,",1]",sep = "")]+
#                            sams[,paste("pz[",i,",2]",sep = "")]))
#   }
#   posts.zoops$f$samples$SurfLarge <- temp1
#   posts.zoops$f$HDI95$SurfLarge <- temp2
#   posts.zoops$f$HDI90$SurfLarge <- temp3
#   posts.zoops$f$HDI75$SurfLarge <- temp4
#   posts.zoops$f$HDI50$SurfLarge <- temp5
#   posts.zoops$f$mean$SurfLarge <- temp6
#   posts.zoops$f$mode$SurfLarge <- temp7
#   

posts.zoops <<- posts.zoops
posts.zoops.long <<- posts.zoops.long

}
