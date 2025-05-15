#### plot_sources #####
## plots the tracer values for organic matter sources
plot_sources <- function(
    Data.sources,
    Tracers
) {
  # creating a "long" version of organic matter source data
  Sources.long <- melt(Data.sources, id.vars=c(Variables), measure.vars =c(Tracers$mix),
                       variable.name = "Tracer", value.name = "Value")
  output <-
    ggplot(Sources.long, aes(x = Tracer, y = Value, color = Source, group=Source))+
    geom_point(alpha=0.4, size = 2, position = position_dodge(width = 0.8))+
    labs(color="Organic Matter Source", shape="Food Web Base", fill="Zooplankton Samples")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust=1))
  return(output)
}




#### plot_sources_PCA #####
## plots PCA of tracer values for organic matter sources
plot_sources_PCA <- function(
    PCA, Data.sources
) {
  output <-
    # Plot PC1&2 results and project variable vectors and print
    fviz_pca_biplot(PCA, axes = c(1,2), # axes tell which components to plot
                    geom = "point", addEllipses = TRUE, repel=TRUE,
                    col.ind = Data.sources[["Source"]],
                    label = c("quali", "var")) + 
    labs(color = "Source", fill = "Source", shape = "Source")
  return(output)
}




#### plot_sources_LDA #####
## plots LDA of tracer values for organic matter sources
plot_sources_LDA <- function(
    class.train
) {
  # Generate biplot with training data only
  output1 <-
    ggplot(data = class.train, 
           aes(x = LD2, y = LD1, color = Source, shape = Source), size = 2) + 
    geom_point(alpha = 0.75, size=3) + 
    stat_ellipse(type = 't', alpha = 0.75) +
    labs(color = "Source", shape = "Source")
  if(length(Sources) > 3){
    output2 <- 
      ggplot(data = class.train, 
             aes(x = LD3, y = LD1, color = Source, shape = Source), size = 2) + 
      geom_point(alpha = 0.75, size=3) +  
      stat_ellipse(type = 't', alpha = 0.75)+
      labs(color = "Source", shape = "Source")
    
    return(
      ggarrange(ncol = 2, widths = c(1,1.2), 
                output1, output2,
                common.legend = TRUE)
    )
  }else{
    return(output1)
  }
}





#### plot_Source_Consumer_Sim ####
## plots tracer values and LDA (if possible) for sources consumers and the food web base
plot_Source_Consumer_Sim <- function(
    Data.sources,
    Data.zoops,
    base.sim,
    Tracers, 
    LDA.full
) {
  
  Sources.long <- melt(Data.sources, id.vars=c(Variables), measure.vars =c(Tracers$MTS, Tracers$FWL, Tracers$mix),
                       variable.name = "Tracer", value.name = "Value")
  Zoops.long <- melt(Data.zoops, id.vars = Variables, measure.vars = c(Tracers$MTS, Tracers$FWL, Tracers$mix),
                     variable.name = "Tracer", value.name = "Value")
  Base.long <- melt(base.sim, id.vars = Variables, measure.vars = c(Tracers$MTS, Tracers$FWL, Tracers$mix),
                    variable.name = "Tracer", value.name = "Value")
  
  AA_plots <- ggplot(Sources.long) +
    geom_point(aes(x = Tracer, y = Value, color = Source, shape = Source), alpha = 0.6,
               position = position_nudge(x = 0.1)) +
    geom_point(data = Zoops.long, aes(x = Tracer, y = Value, fill = "Zooplankton"),
               color = "black", shape = 17, position = position_nudge(x = -0.1)) +
    geom_point(data = Base.long, aes(x = Tracer, y = Value, fill = "Food Web Base"),
               color = "brown", shape = 16, position = position_nudge(x = 0)) +
    ylab(expression(delta^{15}*N*" (\u2030)")) +
    labs(color = "Organic Matter Source", shape = "Organic Matter Source", fill = "Simulated Data") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust = 1), axis.title.x = element_blank())+
    scale_x_discrete(labels = substr(c(Tracers$MTS, Tracers$FWL, Tracers$mix),5,7))
  
  if (runMV) {
    pred.zoops <- predict(LDA.full, Data.zoops)
    class.zoops <- data.frame(Predicted = pred.zoops$class, pred.zoops$x,
                              lapply(Data.zoops[Variables], as.character))
    
    pred.base <- predict(LDA.full, base.sim)
    class.base <- data.frame(Predicted = pred.base$class, pred.base$x,
                             lapply(Data.zoops[Variables], as.character))
    
    plot.mix.LD12 <- 
      ggplot(class.train, aes(x = LD2, y = LD1, color = Source, shape = Source)) +
      geom_point(alpha = 0.6) +
      stat_ellipse(type = "t", alpha = 0.6) +
      geom_point(data = class.base, aes(x = LD2, y = LD1, fill = "Food Web Base"),
                 color = "brown", shape = 16) +
      geom_point(data = class.zoops, aes(x = LD2, y = LD1, fill = "Zooplankton"),
                 color = "black", shape = 17) +
      theme(legend.position = "none")
    # Combine tracer and LDA plots.
    output <- 
      ggarrange(AA_plots, plot.mix.LD12, ncol = 2, widths = c(2, 1.5),
                common.legend = TRUE, legend = "right")
    return(output)
  } else return(AA_plots)
}




#### plot_sourcepost_sim ####
## plots the source tracer posteriors relative to the organic matter source data
plot_sourcepost_sim <- function(
    posts.long, 
    Data.sources, 
    Tracers
){
  
  Sources.long <- melt(Data.sources, id.vars=c(Variables), measure.vars =c(Tracers$MTS, Tracers$FWL, Tracers$mix),
                       variable.name = "Tracer", value.name = "Value")
  
  theme_set(theme_classic2()+
              theme(panel.grid.major.x = element_line(colour = "grey95"),
                    panel.grid.major.y = element_line(colour = "grey95")))
  
  ## plotting posterior δ15N values for organic matter sources compared to data
  output <- 
    ggplot()+
    geom_density(data = posts.long$source$samples[
      which(posts.long$source$samples$Tracer %in% Tracers$mix),],
      aes(x=Value, fill=Source), alpha=0.8, color="grey10", linewidth=0.5)+
    geom_point(data=Sources.long[
      which(Sources.long$Tracer %in% Tracers$mix),],
      aes(x=Value, y=0, fill=Source, shape = "Source Data"),
      color="grey10",size=2,stroke=0.75,
      show.legend = FALSE)+
    scale_shape_manual(values=c(24))+
    xlab("Value")+
    facet_wrap(~Tracer, scales = "free", nrow=2)+
    scale_x_continuous(n.breaks=4)+
    labs(fill="Organic Matter Source",color="95% HDI")+
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),plot.title = element_text(hjust = 0.5),
          legend.position="bottom")
  return(output)
}




#### plot_basepost_sim ####
## generates plots comparing the posteriors for tracer values at the base of the food web with simulated zooplankton and organic matter swource data
plot_basepost_sim <- function(
    posts, 
    posts.long, 
    Data.sources, 
    Data.zoops,
    base.sim,
    Tracers, 
    Variables,
    LDA.full
) {
  
  Sources.long <- melt(Data.sources, id.vars=c(Variables), measure.vars =c(Tracers$MTS, Tracers$FWL, Tracers$mix),
                       variable.name = "Tracer", value.name = "Value")
  Zoops.long <- melt(Data.zoops, id.vars = Variables, measure.vars = c(Tracers$MTS, Tracers$FWL, Tracers$mix),
                     variable.name = "Tracer", value.name = "Value")
  Base.long <- melt(base.sim, id.vars = Variables, measure.vars = c(Tracers$MTS, Tracers$FWL, Tracers$mix),
                    variable.name = "Tracer", value.name = "Value")
  
  ## Generating AA δ15N comparison plot
  posts.long$base$thin$Source <- as.factor(posts.long$base$thin$Source)
  
  AAplots <-
    ggplot(Sources.long)+
    geom_point(aes(x = Tracer, y = Value, color = Source),
               position = position_nudge(x = 0.2))+
    geom_point(data = Zoops.long,
               aes(x = Tracer, y = Value, fill="simulated zooplankton"), pch = "triangle",
               position = position_nudge(x = -0.2))+
    geom_point(data = posts.long$base$thin[which(posts.long$base$thin$Tracer %in% c(Tracers$MTS, Tracers$FWL, Tracers$mix)),],
               aes(x = Tracer, y = Value), color = "grey60", alpha=0.1, pch=16,
               position = position_dodge2(width=0.2))+
    geom_point(data = posts.long$base$mean[which(posts.long$base$mean$Tracer %in% c(Tracers$MTS, Tracers$FWL, Tracers$mix)),],
               aes(x = Tracer, y = Value, shape="posterior mean"), color = "brown", alpha=1,
               position = position_dodge2(width = 0.2))+
    # position = position_nudge(x = 0.25))+
    geom_point(data = Base.long,
               aes(x = Tracer, y = Value, shape="true value"), color = "red", alpha=1,
               position = position_dodge2(width = 0.2))+
    labs(color="Organic Matter Source", shape="Food Web Base", fill="Zooplankton Samples")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust=1))
  
  ## Generating LDA comparison plot
  if(runMV == TRUE){
    # project simulated base into LDA space
    pred.base <- predict(LDA.full, base.sim)
    class.base <- data.frame(Predicted = pred.base$class, pred.base$x,
                             lapply(Data.zoops[Variables], as.character))
    # poject zooplankton into LD space and predict classification
    pred.base.MCMC = predict(LDA.full, 
                             posts$base$mode)
    #class.zoops = data.frame('Class' = pred.zoops$class, pred.zoops$x)
    # generate data frame which stores predicted classification of particles indexed by size and depth
    class.base.MCMC = data.frame('Predicted' = pred.base.MCMC$class, pred.base.MCMC$x,
                                 lapply(posts$base$mode[Variables],as.character))
    # # Predicted probabilities of class membership
    # head(pred.base.MCMC$posterior,n=dim(pred.base.MCMC$posterior))
    # # linear discriminants
    # head(pred.base.MCMC$x,n=dim(pred.base.MCMC$x))
    
    pred.base.sams <- predict(LDA.full,
                              posts$base$thin[Tracers$mix])
    # generate data frame which stores predicted classification of particles indexed by size and depth
    class.base.sams = data.frame('Predicted' = pred.base.sams$class, pred.base.sams$x,
                                 lapply(posts$base$thin[Variables],as.character))
    # Predicted probabilities of class membership
    # head(pred.base.sams$posterior,n=dim(pred.base.sams$posterior))
    # linear discriminants
    # head(pred.base.sams$x,n=dim(pred.base.sams$x))
    
    plot.mix.LD12 <- 
      ggplot(data = class.train, 
             aes(x = LD1, y = LD2), size = 2) + 
      geom_point(aes(shape=Source), alpha = 1, color="grey20") +
      stat_ellipse(aes(group=Source),type = 't', alpha = 0.8)+
      stat_ellipse(data = class.base.sams,
                   aes(x = LD1, y = LD2, color=Source, fill=Source),
                   type = 't', alpha = 0.5, size=0.5, geom = "polygon")+
      geom_point(data = class.base,
                 aes(x = LD1, y = LD2, fill=Source, alpha = ""),
                 shape=23, color="black", size=3, stroke=1)+
      scale_alpha_manual(values = c(1,1,1,1))+
      geom_point(data = class.base.MCMC,
                 aes(x = LD1, y = LD2, fill=Source),
                 size=3.5, shape=21, stroke=1)+
      labs(fill = "Posterior 95% HDI and Mode",
           color = "Simulated Zooplankton Sample",
           alpha = "True Value",
           shape = "Organic Matter Source"
      )+
      guides(color="none")
    output2 <-
      ggarrange(AAplots, plot.mix.LD12, ncol = 2, widths = c(2, 1.5),
                common.legend = TRUE, legend = "left")
    return(output2)
  }else {
    return(AAplots)
  }
}




#### plot_zooppost_sim ####
## plots posterior zooplankton tracer values and compares with tracer values in zooplankton data.
plot_zooppost_sim <- function(
    posts.long, 
    Data.zoops, 
    Data.sources,
    base.sim,
    Tracers
) {
  
  Sources.long <- melt(Data.sources, id.vars=c(Variables), measure.vars =c(Tracers$MTS, Tracers$FWL, Tracers$mix),
                       variable.name = "Tracer", value.name = "Value")
  Zoops.long <- melt(Data.zoops, id.vars = Variables, measure.vars = c(Tracers$MTS, Tracers$FWL, Tracers$mix),
                     variable.name = "Tracer", value.name = "Value")
  Base.long <- melt(base.sim, id.vars = Variables, measure.vars = c(Tracers$MTS, Tracers$FWL, Tracers$mix),
                    variable.name = "Tracer", value.name = "Value")
  
  Ymax <- data.frame(
    PTS = max(posts$trophic$HDI95$PTS),
    MTS = max(posts$trophic$HDI95$MTS),
    FWL = max(posts$trophic$HDI95$FWL)
  )
  
  ## Generating AA δ15N comparison plot
  posts.long$zoop$thin$Source <- as.factor(posts.long$zoop$thin$Source)
  
  output <-
    ggplot(subset(Sources.long, Tracer == Tracers$mix))+
    geom_point(aes(x = Tracer, y = Value, color = Source),
               position = position_nudge(x = 0.2))+
    geom_point(data = subset(Zoops.long, Tracer == Tracers$mix),
               aes(x = Tracer, y = Value, fill="simulated zooplankton"), pch = "triangle",
               position = position_nudge(x = -0.2))+
    geom_point(data = posts.long$zoop$thin[which(posts.long$zoop$thin$Tracer %in% Tracers$mix),],
               aes(x = Tracer, y = Value), color = "grey60", alpha=0.1, pch=16,
               position = position_dodge2(width=0.1))+
    geom_point(data = posts.long$zoop$mean[which(posts.long$zoop$mean$Tracer %in% Tracers$mix),],
               aes(x = Tracer, y = Value, shape="posterior mean"), color = "brown", alpha=1,
               position = position_dodge2(width = 0.1))+
    labs(color="Organic Matter Source", shape="Zooplankton Posterior", fill="Zooplankton Samples")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust=1))
  return(output)
}




#### plot_FWLpost_sim ####
## plots the posterior for FWL relative to the true value for simulated samples
plot_FWLpost_sim <- function(
    posts, 
    zoops.f
) {
  
  nzoops <- nrow(zoops.f)
  
  Ymax <- data.frame(
    PTS = max(posts$trophic$HDI95$PTS),
    MTS = max(posts$trophic$HDI95$MTS),
    FWL = max(posts$trophic$HDI95$FWL)
  )
  
  output <-
    ggplot()+
    geom_line(data = posts$trophic$HDI95,
              aes(x = Source, y = FWL, group = Source),
              alpha = 0.5, size = 1, color = "steelblue4")+
    geom_line(data = posts$trophic$HDI90,
              aes(x = Source, y = FWL, group = Source),
              alpha = 0.7, size = 1.5, color = "steelblue4")+
    geom_line(data = posts$trophic$HDI75,
              aes(x = Source, y = FWL, group = Source),
              alpha = 0.9, size = 2, color = "steelblue4")+
    geom_line(data = posts$trophic$HDI50,
              aes(x = Source, y = FWL, group = Source),
              alpha = 1, size = 3, color = "steelblue4")+
    geom_text(data = posts$trophic$mode,
              aes(x = Source, y = FWL, group = Source, 
                  label = round(FWL,1)),
              size = 1.6, color = "lightskyblue1")+
    scale_linetype_manual(values=c(1,1))+
    geom_segment(
      aes(y=zoops.f$FWL, x=seq(1,nzoops)-0.5, xend=seq(1,nzoops)+0.5,
          color=""),
      lty=1,size=0.5,alpha=1)+
    scale_color_manual(values=c("goldenrod"))+
    labs(fill="Model Posterior", color="True Value",lty="95% HDI")+
    xlab("Simulated Zooplankton Sample") + ylab("FWL")+
    scale_x_continuous(breaks = seq(1,nzoops,2))+
    scale_y_continuous(breaks = c(-1.0,-0.5, 0.0, 0.5 , 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0), 
                       labels = c(-1.0, "" , 0.0,  "" , 1.0, "" , 2.0, "" , 3.0, "" , 4.0, "" , 5.0))+
    coord_cartesian(ylim = c(0,Ymax$FWL), xlim = c(0,nzoops+1), expand = FALSE)
  return(output)
}




#### plot_PTSpost_sim ####
## plots the posterior for PTS relative to the true value for simulated samples
plot_PTSpost_sim <- function(
    posts, 
    zoops.f
) {
  
  nzoops <- nrow(zoops.f)
  
  Ymax <- data.frame(
    PTS = max(posts$trophic$HDI95$PTS),
    MTS = max(posts$trophic$HDI95$MTS),
    FWL = max(posts$trophic$HDI95$FWL)
  )
  
  output <-
    ggplot()+
    geom_line(data = posts$trophic$HDI95,
              aes(x = Source, y = PTS, group = Source),
              alpha = 0.5, size = 1, color = "steelblue4")+
    geom_line(data = posts$trophic$HDI90,
              aes(x = Source, y = PTS, group = Source),
              alpha = 0.7, size = 1.5, color = "steelblue4")+
    geom_line(data = posts$trophic$HDI75,
              aes(x = Source, y = PTS, group = Source),
              alpha = 0.9, size = 2, color = "steelblue4")+
    geom_line(data = posts$trophic$HDI50,
              aes(x = Source, y = PTS, group = Source),
              alpha = 1, size = 3, color = "steelblue4")+
    geom_text(data = posts$trophic$mode,
              aes(x = Source, y = PTS, group = Source, 
                  label = round(PTS,1)),
              size = 1.6, color = "lightskyblue1")+
    scale_linetype_manual(values=c(1,1))+
    geom_segment(
      aes(y=zoops.f$PTS, x=seq(1,nzoops)-0.5, xend=seq(1,nzoops)+0.5,
          color=""),
      lty=1,size=0.5,alpha=1)+
    scale_color_manual(values=c("goldenrod"))+
    labs(fill="Model Posterior", color="True Value",lty="95% HDI")+
    xlab("Simulated Zooplankton Sample") + ylab("PTS")+
    scale_x_continuous(breaks = seq(1,nzoops,2))+
    scale_y_continuous(breaks = c(-1.0,-0.5, 0.0, 0.5 , 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0), 
                       labels = c(-1.0, "" , 0.0,  "" , 1.0, "" , 2.0, "" , 3.0, "" , 4.0, "" , 5.0))+
    coord_cartesian(ylim = c(0,Ymax$FWL), xlim = c(0,nzoops+1), expand = FALSE)
  return(output)
}




#### plot_MTSpost_sim ####
## plots the posterior for MTS relative to the true value for simulated samples
plot_MTSpost_sim <- function(
    posts, 
    zoops.f
) {
  
  nzoops <- nrow(zoops.f)
  
  Ymax <- data.frame(
    PTS = max(posts$trophic$HDI95$PTS),
    MTS = max(posts$trophic$HDI95$MTS),
    FWL = max(posts$trophic$HDI95$FWL)
  )
  
  output <-
    ggplot()+
    geom_line(data = posts$trophic$HDI95,
              aes(x = Source, y = MTS, group = Source),
              alpha = 0.5, size = 1, color = "steelblue4")+
    geom_line(data = posts$trophic$HDI90,
              aes(x = Source, y = MTS, group = Source),
              alpha = 0.7, size = 1.5, color = "steelblue4")+
    geom_line(data = posts$trophic$HDI75,
              aes(x = Source, y = MTS, group = Source),
              alpha = 0.9, size = 2, color = "steelblue4")+
    geom_line(data = posts$trophic$HDI50,
              aes(x = Source, y = MTS, group = Source),
              alpha = 1, size = 3, color = "steelblue4")+
    geom_text(data = posts$trophic$mode,
              aes(x = Source, y = MTS, group = Source, 
                  label = round(MTS,1)),
              size = 1.6, color = "lightskyblue1")+
    scale_linetype_manual(values=c(1,1))+
    geom_segment(
      aes(y=zoops.f$MTS, x=seq(1,nzoops)-0.5, xend=seq(1,nzoops)+0.5,
          color=""),
      lty=1,size=0.5,alpha=1)+
    scale_color_manual(values=c("goldenrod"))+
    labs(fill="Model Posterior", color="True Value",lty="95% HDI")+
    xlab("Simulated Zooplankton Sample") + ylab("MTS")+
    scale_x_continuous(breaks = seq(1,nzoops,2))+
    scale_y_continuous(breaks = c(-1.0,-0.5, 0.0, 0.5 , 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0), 
                       labels = c(-1.0, "" , 0.0,  "" , 1.0, "" , 2.0, "" , 3.0, "" , 4.0, "" , 5.0))+
    coord_cartesian(ylim = c(0,Ymax$FWL), xlim = c(0,nzoops+1), expand = FALSE)
  return(output)
}




#### plot_Trophic_TrueMod ####
plot_Trophic_TrueMod <- function(
    data.truemod.long,
    variable,
    Ymax
){
  
  output <-
    ggplot(data = subset(data.truemod.long, parameter == variable))+
    geom_point(
      aes(x = trophic_true, y = trophic_model),
      size = 1, stroke=1, shape = 1, color = "steelblue4"
    )+
    geom_smooth(
      aes(x = trophic_true, y = trophic_model),
      formula = y ~ x,
      color = "royalblue4", fill = "skyblue1",
      method = "lm"
      # method = "glm", method.args = list(family = binomial())
    )+
    geom_abline(slope = 1, intercept = 0, size = 0.5, color = "goldenrod")+
    xlab("True Value")+
    coord_cartesian(ylim = c(0,Ymax$FWL), xlim = c(0,Ymax$FWL), expand = FALSE)+
    scale_x_continuous(breaks = c(-1.0,-0.5, 0.0, 0.5 , 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0), 
                       labels = c(-1.0, "" , 0.0,  "" , 1.0, "" , 2.0, "" , 3.0, "" , 4.0, "" , 5.0))+
    scale_y_continuous(breaks = c(-1.0,-0.5, 0.0, 0.5 , 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0), 
                       labels = c(-1.0, "" , 0.0,  "" , 1.0, "" , 2.0, "" , 3.0, "" , 4.0, "" , 5.0))
  return(output)
}




#### plot_Trophic_TrueMod_Disc ####
plot_Trophic_TrueMod_Disc <- function(
    data.truemod.long,
    variable,
    Ymax
){
  
  output <-
    ggplot(data = subset(data.truemod.long, parameter == variable))+
    geom_point(
      aes(x = trophic_true, y = trophic_disc),
      size = 1, stroke=1, shape = 1, color = "steelblue4"
    )+
    geom_smooth(
      aes(x = trophic_true, y = trophic_disc),
      formula = y ~ x,
      color = "royalblue4", fill = "skyblue1",
      method = "lm"
      # method = "glm", method.args = list(family = binomial())
    )+
    geom_abline(slope = 0, intercept = 0, size = 0.5, color = "goldenrod")+
    xlab("True Value")+
    coord_cartesian(ylim = c(-1,1), xlim = c(0,Ymax$FWL), expand = FALSE)+
    scale_x_continuous(breaks = c(-1.0,-0.5, 0.0, 0.5 , 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0), 
                       labels = c(-1.0, "" , 0.0,  "" , 1.0, "" , 2.0, "" , 3.0, "" , 4.0, "" , 5.0))+
    scale_y_continuous(breaks = c(-1,-0.50,0,0.5,1), labels = c("-1","","0","","1"))
  return(output)
}




#### plotall_trophicpost_sim ####
## plots the posteriors for trphic parameters relative to the true values for simulated samples
plotall_trophicpost_sim <- function(
    posts, 
    zoops.f
) {
  
  Ymax <- data.frame(
    PTS = max(posts$trophic$HDI95$PTS),
    MTS = max(posts$trophic$HDI95$MTS),
    FWL = max(posts$trophic$HDI95$FWL)
  )
  
  nzoops <- nrow(zoops.f)
  
  plot.PTS <- plot_PTSpost_sim(posts, zoops.f)
  plot.MTS <- plot_MTSpost_sim(posts, zoops.f)
  plot.FWL <- plot_FWLpost_sim(posts, zoops.f)
  
  data.true <- zoops.f[c("PTS","MTS","FWL")]
  data.model <- posts$trophic$mode[c("PTS","MTS","FWL")]
  
  mix.disc <- data.model-data.true
  mix.disc.abs <- abs(mix.disc[c("PTS","MTS","FWL")])
  
  data.model.long <- 
    melt(data.model[c("PTS","MTS","FWL")], 
         value.name = "trophic_model",
         variable.name = "parameter"
    )
  data.true.long <- 
    melt(data.true[c("PTS","MTS","FWL")], 
         value.name = "trophic_true",
         variable.name = "parameter"
    )
  data.truemod.long <-
    cbind(
      data.true.long,
      data.model.long["trophic_model"]
    )
  data.truemod.long$trophic_disc <-
    data.truemod.long$trophic_model -
    data.truemod.long$trophic_true
  
  plot_list2 <- 
    lapply(
      c("PTS","MTS","FWL"), 
      plot_Trophic_TrueMod, 
      data.truemod.long = data.truemod.long,
      Ymax = Ymax)
  for (i in seq(1,2)) {
    plot_list2[[i]] <- plot_list2[[i]] + no.x.axis
  }
  for (i in seq(1,3)) {
    plot_list2[[i]] <- plot_list2[[i]] + theme(axis.title.y = element_blank())
  }
  plot_list3 <- 
    lapply(
      c("PTS","MTS","FWL"), 
      plot_Trophic_TrueMod_Disc, 
      data.truemod.long = data.truemod.long,
      Ymax = Ymax)
  for (i in seq(1,2)) {
    plot_list3[[i]] <- plot_list3[[i]] + no.x.axis
  }
  for (i in seq(1,3)) {
    plot_list3[[i]] <- plot_list3[[i]] + theme(axis.title.y = element_blank())
  }
  
  col1<-
    ggarrange(
      plot.PTS + no.x.axis, 
      plot.MTS + no.x.axis, 
      plot.FWL,
      ncol = 1, heights = c(2,2,2.6), common.legend = TRUE , legend="none")
  col2<-
    ggarrange(
      plotlist = plot_list2,
      ncol = 1, heights = c(2,2,2.6), common.legend = TRUE , legend="none")
  col3<-
    ggarrange(
      plotlist = plot_list3,
      ncol = 1, heights = c(2,2,2.6), common.legend = TRUE , legend="none")
  
  output <-
    ggarrange(
      ggarrange(
        ggplot()+
          ggtitle("Posterior Mode and HDI by Sample")+
          theme(plot.title = element_text(hjust = 0.6)),
        ggplot()+
          ggtitle("Modeled\nvs True")+
          theme(plot.title = element_text(hjust = 0.5)),
        ggplot()+
          ggtitle("Discrepancy")+
          theme(plot.title = element_text(hjust = 1)),
        ncol = 3, widths = c(4,1,1)
      ),
      ggarrange(
        col1, col2, col3,
        ncol = 3, widths = c(4,1,1)
      ),
      nrow = 2, heights = c(0.15,1)
    )
  return(output)
  
}




#### plot_Post_Error_Trophic ####
## plots modelled vs true regressions for trophic parameters
plot_Post_Error_Trophic <- function(
    zoops.f,
    posts
){
  data.true <- zoops.f[c("PTS","MTS","FWL")]
  data.model <- posts$trophic$mode[c("PTS","MTS","FWL")]
  
  
  data.true.long <-
    melt(
      data.true,
      value.name = "true_value",
      variable.name = "Parameter"
    )
  data.model.long <- 
    melt(
      data.model,
      value.name = "model_value",
      variable.name = "Parameter"
    )
  
  data.truemod.long <-
    cbind(
      data.true.long,
      data.model.long[,"model_value"])
  colnames(data.truemod.long) <-
    c("parameter","true_value","model_value")
  
  data.truemod.long$disc <-
    data.truemod.long$model_value -
    data.truemod.long$true_value
  output1 <- 
    ggplot(data = data.truemod.long)+
    geom_abline(slope = 1, intercept = 0, size = 1)+
    geom_point(
      aes(x = true_value, y = model_value,
          color = parameter), shape = 1,
      position = position_dodge(width = 0.2),
      size = 2
    )+
    geom_smooth(
      aes(x = true_value, y = model_value,
          color = parameter, fill = parameter),
      formula = y ~ x,
      method = "lm"
    )+
    xlab("True Value")+
    ylab("Modeled Trophic Parameter")+
    coord_cartesian(expand = FALSE)
  
  output2 <- 
    ggplot(data = data.truemod.long,)+
    geom_abline(slope = 0, intercept = 0, size = 1)+
    geom_point(
      aes(x = true_value, y = disc,
          color = parameter, group = parameter),
      position = position_dodge(width = 0.05),
      size = 2, stroke=1, shape = 1
    )+
    geom_smooth(
      aes(x = true_value, y = disc,
          color = parameter, fill = parameter),
      formula = y ~ x,
      method = "lm"
    )+
    xlab("True Value")+
    ylab("Modeled - True")+
    coord_cartesian(ylim = c(-1,1),expand = FALSE)
  
  
  return(list(output1, output2))
}




#### plot_fpost_sim ####
## Generates a single plot comparing mixing coefficient posteriors with true values for simiulated data
plot_fpost_sim = function(
    posts, 
    zoops.f, 
    variable
) {
  nzoops <- nrow(zoops.f)
  nsources <- length(Sources)
  
  posts$f$samples$Source <- as.character(posts$f$samples$Source)
  posts$f$mode[paste0(Sources,"round")] <- round(posts$f$mode[Sources]*100,0)
  variableround <- paste0(variable,"round")
  
  zoops.f$Samplemin <- seq(1,nzoops)-0.5
  zoops.f$Samplemax <- seq(1,nzoops)+0.5
  zoops.f$Color <- ""
  
  output <- 
    ggplot()+
    scale_fill_manual(values=c("steelblue1","steelblue4"))+
    geom_line(
      data = posts$f$HDI95,
      aes_string(x = "Source", y = variable, group = "Source"),
      alpha = 0.5, size = 1, color = "steelblue4")+
    geom_line(
      data = posts$f$HDI90,
      aes_string(x = "Source", y = variable, group = "Source"),
      alpha = 0.7, size = 1.5, color = "steelblue4")+
    geom_line(
      data = posts$f$HDI75,
      aes_string(x = "Source", y = variable, group = "Source"),
      alpha = 0.9, size = 2, color = "steelblue4")+
    geom_line(
      data = posts$f$HDI50,
      aes_string(x = "Source", y = variable, group = "Source"),
      alpha = 1, size = 3, color = "steelblue4")+
    geom_text(
      data = posts$f$mode,
      aes_string(x = "Source", y = variable, group = "Source",
                 label = variableround),
      size = 1.6, color = "lightskyblue1")+
    geom_segment(
      data = zoops.f,
      aes_string(y=variable, x="Samplemin", xend="Samplemax", color="Color"),
      lty=1,size=0.5,alpha=1)+
    scale_color_manual(values=c("goldenrod"))+
    labs(fill="Model Posterior", color="True Value", lty="95% HDI")+
    xlab("Simulated Zooplankton Sample") + ylab(paste0("f(",variable,")")) +
    scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1), labels = c("0","","0.5","","1.0"))+
    scale_x_continuous(breaks = seq(1,nzoops,2))+
    coord_cartesian(ylim = c(0,1), xlim = c(0,nzoops+1), expand = FALSE)
  return(output)
}




#### plot_fpost_TrueMod ####
plot_fpost_TrueMod <- function(
    data.truemod.long,
    variable
){
  
  output <-
    ggplot(data = subset(data.truemod.long, Source == variable))+
    geom_point(
      aes(x = f_true, y = f_model),
      size = 1, stroke=1, shape = 1, color = "steelblue4"
    )+
    geom_smooth(
      aes(x = f_true, y = f_model),
      formula = y ~ x,
      color = "royalblue4", fill = "skyblue1",
      # method = "lm"
      method = "glm", method.args = list(family = binomial())
    )+
    geom_abline(slope = 1, intercept = 0, size = 0.5, color = "goldenrod")+
    xlab("True Value")+
    coord_cartesian(ylim = c(0,1), xlim = c(0,1), expand = FALSE)+
    scale_x_continuous(breaks = c(0,0.5,1))+
    scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1), labels = c("0","","0.5","","1.0"))
  return(output)
}




#### plot_fpost_TrueMod_Disc ####
plot_fpost_TrueMod_Disc <- function(
    data.truemod.long,
    variable
){
  
  output <-
    ggplot(data = subset(data.truemod.long, Source == variable))+
    geom_point(
      aes(x = f_true, y = f_disc),
      size = 1, stroke=1, shape = 1, color = "steelblue4"
    )+
    geom_smooth(
      aes(x = f_true, y = f_disc),
      formula = y ~ x,
      color = "royalblue4", fill = "skyblue1",
      method = "lm"
      # method = "glm", method.args = list(family = binomial())
    )+
    geom_abline(slope = 0, intercept = 0, size = 0.5, color = "goldenrod")+
    xlab("True Value")+
    coord_cartesian(ylim = c(-1,1), xlim = c(0,1), expand = FALSE)+
    scale_x_continuous(breaks = c(0,0.5,1))+
    scale_y_continuous(breaks = c(-1,-0.50,0,0.5,1), labels = c("-1","","0","","1"))
  return(output)
}




#### plotall_fpost_sim ####
## generates a multipanel plot containing all of the true vs modelled mixing coefficient plots
plotall_fpost_sim <- function(
    posts, 
    zoops.f, 
    Sources
) {
  
  theme_set(theme_classic2()+
              theme(panel.grid.major.x = element_line(colour = "grey95"),
                    panel.grid.major.y = element_line(colour = "grey95")))
  
  plot_list1 <- 
    lapply(
      Sources, 
      plot_fpost_sim, 
      posts = posts, zoops.f = zoops.f)
  for (i in seq(1,nsources-1)) {
    plot_list1[[i]] <- plot_list1[[i]] + no.x.axis
  }
  
  data.true <- zoops.f[Sources]
  data.model <- posts$f$mode[Sources]
  
  mix.disc <- data.model-data.true
  mix.disc.abs <- abs(mix.disc[Sources])
  
  data.model.long <- 
    melt(data.model[Sources], 
         value.name = "f_model",
         variable.name = "Source"
    )
  data.true.long <- 
    melt(data.true[Sources], 
         value.name = "f_true",
         variable.name = "Source"
    )
  data.truemod.long <-
    cbind(
      data.true.long,
      data.model.long["f_model"]
    )
  data.truemod.long$f_disc <-
    data.truemod.long$f_model -
    data.truemod.long$f_true
  
  plot_list2 <- 
    lapply(
      Sources, 
      plot_fpost_TrueMod, 
      data.truemod.long = data.truemod.long)
  for (i in seq(1,nsources-1)) {
    plot_list2[[i]] <- plot_list2[[i]] + no.x.axis
  }
  for (i in seq(1,nsources)) {
    plot_list2[[i]] <- plot_list2[[i]] + theme(axis.title.y = element_blank())
  }
  plot_list3 <- 
    lapply(
      Sources, 
      plot_fpost_TrueMod_Disc, 
      data.truemod.long = data.truemod.long)
  for (i in seq(1,nsources-1)) {
    plot_list3[[i]] <- plot_list3[[i]] + no.x.axis
  }
  for (i in seq(1,nsources)) {
    plot_list3[[i]] <- plot_list3[[i]] + theme(axis.title.y = element_blank())
  }
  
  plot_heights <- c(rep(1,nsources-1),1+(0.5-nsources*0.07))
  
  nsource <- length(Sources)
  col1<-
    ggarrange(
      plotlist = plot_list1,
      nrow = nsource, ncol = 1, heights = plot_heights, legend = "none", common.legend = TRUE)
  col2<-
    ggarrange(
      plotlist = plot_list2,
      nrow = nsource, ncol = 1, heights = plot_heights, legend = "none", common.legend = TRUE)
  col3<-
    ggarrange(
      plotlist = plot_list3,
      nrow = nsource, ncol = 1, heights = plot_heights, legend = "none", common.legend = TRUE)
  output <-
    ggarrange(
      ggarrange(
        ggplot()+
          ggtitle("Posterior Mode and HDI by Sample")+
          theme(plot.title = element_text(hjust = 0.6)),
        ggplot()+
          ggtitle("Modeled\nvs True")+
          theme(plot.title = element_text(hjust = 0.5)),
        ggplot()+
          ggtitle("Discrepancy")+
          theme(plot.title = element_text(hjust = 1)),
        ncol = 3, widths = c(4,0.8,1)
      ),
      ggarrange(
        col1, col2, col3,
        ncol = 3, widths = c(4,1,1)
      ),
      nrow = 2, heights = c(0.17-(nsources*0.015),1)
    )
  return(output)
}




#### plot_Post_Error_f ####
## plots modelled vs true regressions for trophic parameters
plot_Post_Error_f <- function(
    zoops.f,
    posts,
    Sources
){
  
  data.true <- zoops.f[Sources]
  data.model <- posts$f$mode[Sources]
  
  mix.disc <- data.model-data.true
  mix.disc.abs <- abs(mix.disc[Sources])
  
  data.model.long <- 
    melt(data.model[Sources], 
         value.name = "f_model",
         variable.name = "Source"
    )
  data.true.long <- 
    melt(data.true[Sources], 
         value.name = "f_true",
         variable.name = "Source"
    )
  data.truemod.long <-
    cbind(
      data.true.long,
      data.model.long["f_model"]
    )
  data.truemod.long$f_disc <-
    data.truemod.long$f_model -
    data.truemod.long$f_true
  
  output1 <-
    ggplot(data = data.truemod.long,)+
    geom_abline(slope = 1, intercept = 0, size = 1)+
    geom_point(
      aes(x = f_true, y = f_model,
          color = Source, group = Source),
      position = position_dodge(width = 0.05),
      size = 2, stroke=1, shape = 1
    )+
    geom_smooth(
      aes(x = f_true, y = f_model,
          color = Source, fill = Source),
      formula = y ~ x,
      method = "glm", method.args = list(family = binomial())
    )+
    xlab("True Value")+
    ylab("Modeled Mixing Coefficient")+
    coord_cartesian(ylim = c(0,1), expand = FALSE)
  
  output2 <-
    ggplot(data = data.truemod.long,)+
    geom_abline(slope = 0, intercept = 0, size = 1)+
    geom_point(
      aes(x = f_true, y = f_disc,
          color = Source, group = Source),
      position = position_dodge(width = 0.05),
      size = 2, stroke=1, shape = 1
    )+
    geom_smooth(
      aes(x = f_true, y = f_disc,
          color = Source, fill = Source),
      formula = y ~ x,
      method = "lm"
    )+
    xlab("True Value")+
    ylab("Modeled - True")+
    coord_cartesian(ylim = c(-1,1),expand = FALSE)
  
  return(list(output1, output2))
}

