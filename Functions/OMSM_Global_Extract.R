OMSM_Global_Extract <- function(
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

## 1. Global model parameter posteriors (TDFs, d15N(sources)) ##
# we'll make a list which contains samples and summary data for each set of posteriors
posts.global <- list(
  TDF_meta = list(
    samples = data.frame(matrix(nrow = nsams,ncol=0)),
    HDI95 = data.frame(matrix(nrow = 2,ncol=0)),
    HDI90 = data.frame(matrix(nrow = 2,ncol=0)),
    HDI75 = data.frame(matrix(nrow = 2,ncol=0)),
    HDI50 = data.frame(matrix(nrow = 2,ncol=0)),
    mean = data.frame(matrix(nrow = 1,ncol=0)),
    mode = data.frame(matrix(nrow = 1,ncol=0))
  ),
  TDF_proto = list(
    samples = data.frame(matrix(nrow = nsams,ncol=0)),
    HDI95 = data.frame(matrix(nrow = 2,ncol=0)),
    HDI90 = data.frame(matrix(nrow = 2,ncol=0)),
    HDI75 = data.frame(matrix(nrow = 2,ncol=0)),
    HDI50 = data.frame(matrix(nrow = 2,ncol=0)),
    mean = data.frame(matrix(nrow = 1,ncol=0)),
    mode = data.frame(matrix(nrow = 1,ncol=0))
  ),
  mean_Sources = list(
    samples = data.frame(matrix(nrow = nsams*nsources,ncol=0)),
    HDI95 = data.frame(matrix(nrow = 2*nsources,ncol=0)),
    HDI90 = data.frame(matrix(nrow = 2*nsources,ncol=0)),
    HDI75 = data.frame(matrix(nrow = 2*nsources,ncol=0)),
    HDI50 = data.frame(matrix(nrow = 2*nsources,ncol=0)),
    mean = data.frame(matrix(nrow = 1*nsources,ncol=0)),
    mode = data.frame(matrix(nrow = 1*nsources,ncol=0))
  )
)

# grabbing TDF posteriors
for (i in which(Tracers.all %in% Tracers.frac)) {
  # samples of TDF_meta
  posts.global$TDF_meta$samples[,Tracers.all[i]] <- 
    sams[,paste("TDF_meta[",i,"]",sep = "")]
  # HDI min and max of TDF_meta
  posts.global$TDF_meta$HDI95[,Tracers.all[i]] <- 
    HDIofMCMC(sams[,paste("TDF_meta[",i,"]",sep = "")], credMass=0.95)
  posts.global$TDF_meta$HDI90[,Tracers.all[i]] <- 
    HDIofMCMC(sams[,paste("TDF_meta[",i,"]",sep = "")], credMass=0.90)
  posts.global$TDF_meta$HDI75[,Tracers.all[i]] <- 
    HDIofMCMC(sams[,paste("TDF_meta[",i,"]",sep = "")], credMass=0.75)
  posts.global$TDF_meta$HDI50[,Tracers.all[i]] <- 
    HDIofMCMC(sams[,paste("TDF_meta[",i,"]",sep = "")], credMass=0.50)
  # mean of TDF_meta
  posts.global$TDF_meta$mean[,Tracers.all[i]] <- 
    mean(sams[,paste("TDF_meta[",i,"]",sep = "")])
  # mean of TDF_meta
  posts.global$TDF_meta$mode[,Tracers.all[i]] <- 
    post.mode(sams[,paste("TDF_meta[",i,"]",sep = "")])
}
for (i in which(Tracers.all %in% Tracers.varTDF)) {
  # samples of TDF proto
  posts.global$TDF_proto$samples[Tracers.all[i]] <- 
    sams[,paste("TDF_proto[",i,"]",sep = "")]
  # HDI min and max of TDF_proto
  posts.global$TDF_proto$HDI95[Tracers.all[i]] <- 
    HDIofMCMC(sams[,paste("TDF_proto[",i,"]",sep = "")], credMass = 0.95)
  posts.global$TDF_proto$HDI90[Tracers.all[i]] <- 
    HDIofMCMC(sams[,paste("TDF_proto[",i,"]",sep = "")], credMass = 0.90)
  posts.global$TDF_proto$HDI75[Tracers.all[i]] <- 
    HDIofMCMC(sams[,paste("TDF_proto[",i,"]",sep = "")], credMass = 0.75)
  posts.global$TDF_proto$HDI50[Tracers.all[i]] <- 
    HDIofMCMC(sams[,paste("TDF_proto[",i,"]",sep = "")], credMass = 0.50)
  # mean of TDF_proto
  posts.global$TDF_proto$mean[Tracers.all[i]] <- 
    mean(sams[,paste("TDF_proto[",i,"]",sep = "")])
  # mean of TDF_proto
  posts.global$TDF_proto$mode[Tracers.all[i]] <- 
    post.mode(sams[,paste("TDF_proto[",i,"]",sep = "")])
}
# grabbing d15N(source) posteriors
source_alpha <- c("A","B","C","D","E","F")
temp1<-temp2<-temp3<-c()
# 1st column tells which source
for (j in 1:length(Sources)) {
  temp1 <- c(temp1, rep(Sources[j],nsams))
  temp2 <- c(temp2, rep(Sources[j],2))
  temp3 <- c(temp3, Sources[j])
}
posts.global$mean_Sources$samples$Group <- temp1
posts.global$mean_Sources$HDI95$Group <-
  posts.global$mean_Sources$HDI90$Group <-
  posts.global$mean_Sources$HDI75$Group <-
  posts.global$mean_Sources$HDI50$Group <- temp2
posts.global$mean_Sources$mean$Group <-
  posts.global$mean_Sources$mode$Group <-temp3
# next columns give respective d15N data
for (i in 1:length(Tracers.all)) {
  temp1<-temp2<-temp3<-temp4<-temp5<-temp6<-temp7<-c()
  for (j in 1:length(Sources)) {
    temp1 <- c(temp1, sams[,paste("mean_",source_alpha[j],"[",i,"]",sep = "")])
    temp2 <- 
      c(temp2, HDIofMCMC(sams[,paste("mean_",source_alpha[j],"[",i,"]",sep = "")], credMass=0.95))
    temp3 <- 
      c(temp3, HDIofMCMC(sams[,paste("mean_",source_alpha[j],"[",i,"]",sep = "")], credMass=0.90))
    temp4 <- 
      c(temp4, HDIofMCMC(sams[,paste("mean_",source_alpha[j],"[",i,"]",sep = "")], credMass=0.75))
    temp5 <- 
      c(temp5, HDIofMCMC(sams[,paste("mean_",source_alpha[j],"[",i,"]",sep = "")], credMass=0.50))
    temp6 <- c(temp6, mean(sams[,paste("mean_",source_alpha[j],"[",i,"]",sep = "")]))
    temp7 <- c(temp7, post.mode(sams[,paste("mean_",source_alpha[j],"[",i,"]",sep = "")]))
  }
  posts.global$mean_Sources$samples[Tracers.all[i]] <- temp1
  posts.global$mean_Sources$HDI95[Tracers.all[i]] <- temp2
  posts.global$mean_Sources$HDI90[Tracers.all[i]] <- temp3
  posts.global$mean_Sources$HDI75[Tracers.all[i]] <- temp4
  posts.global$mean_Sources$HDI50[Tracers.all[i]] <- temp5
  posts.global$mean_Sources$mean[Tracers.all[i]] <- temp6
  posts.global$mean_Sources$mode[Tracers.all[i]] <- temp7
}

# at times we'll want this data in long format so we'll handle that here as well
posts.global.long <- list(
  TDF_meta = list(
    samples = melt(posts.global$TDF_meta$samples, value.name = "Value", variable.name = "Tracer"),
    HDI95 = melt(posts.global$TDF_meta$HDI95, value.name = "Value", variable.name = "Tracer"),
    HDI90 = melt(posts.global$TDF_meta$HDI90, value.name = "Value", variable.name = "Tracer"),
    HDI75 = melt(posts.global$TDF_meta$HDI75, value.name = "Value", variable.name = "Tracer"),
    HDI50 = melt(posts.global$TDF_meta$HDI50, value.name = "Value", variable.name = "Tracer"),
    mean = melt(posts.global$TDF_meta$mean, value.name = "Value", variable.name = "Tracer"),
    mode = melt(posts.global$TDF_meta$mode, value.name = "Value", variable.name = "Tracer")
  ),
  TDF_proto = list(
    samples = melt(posts.global$TDF_proto$samples, value.name = "Value", variable.name = "Tracer"),
    HDI95 = melt(posts.global$TDF_proto$HDI95, value.name = "Value", variable.name = "Tracer"),
    HDI90 = melt(posts.global$TDF_proto$HDI90, value.name = "Value", variable.name = "Tracer"),
    HDI75 = melt(posts.global$TDF_proto$HDI75, value.name = "Value", variable.name = "Tracer"),
    HDI50 = melt(posts.global$TDF_proto$HDI50, value.name = "Value", variable.name = "Tracer"),
    mean = melt(posts.global$TDF_proto$mean, value.name = "Value", variable.name = "Tracer"),
    mode = melt(posts.global$TDF_proto$mode, value.name = "Value", variable.name = "Tracer")
  ),
  mean_Sources = list(
    samples = melt(posts.global$mean_Sources$samples, 
                   id.vars = "Group", value.name = "Value", variable.name = "Tracer"),
    HDI95 = melt(posts.global$mean_Sources$HDI95, 
                 id.vars = "Group", value.name = "Value", variable.name = "Tracer"),
    HDI90 = melt(posts.global$mean_Sources$HDI90, 
                 id.vars = "Group", value.name = "Value", variable.name = "Tracer"),
    HDI75 = melt(posts.global$mean_Sources$HDI75, 
                 id.vars = "Group", value.name = "Value", variable.name = "Tracer"),
    HDI50 = melt(posts.global$mean_Sources$HDI50, 
                 id.vars = "Group", value.name = "Value", variable.name = "Tracer"),
    mean = melt(posts.global$mean_Sources$mean, 
                id.vars = "Group", value.name = "Value", variable.name = "Tracer"),
    mode = melt(posts.global$mean_Sources$mode, 
                id.vars = "Group", value.name = "Value", variable.name = "Tracer")
  )
)

posts.global <<- posts.global
posts.global.long <<- posts.global.long
}