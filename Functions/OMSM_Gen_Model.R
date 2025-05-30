# The model is written in the BUGS language and stored as a text string. 
# To make the model more flexible, we will define the model in multiple strings, and the concatenate them together.
# We have multiple string options up through the mixing model portion of the code to accommodate different number of organic matter sources.
OMSM_Gen_Model <- function(
    N_S,
    Data.OMSM
){
  if(Data.OMSM$i_T_non[1]==0){Tracers_NoFrac <- FALSE}else{Tracers_NoFrac <- TRUE}
  if(Data.OMSM$i_T_const[1]==0){Tracers_ConstFrac <- FALSE}else{Tracers_ConstFrac <- TRUE}
  if(Data.OMSM$i_T_var[1]==0){Tracers_VarFrac <- FALSE}else{Tracers_VarFrac <- TRUE}
  
OMSM_Source_Data_and_Priors_2 <-
  "/* observe tracers in sources */
  # Each row will be a sample
  # Each column will be a tracer
    
  # Source A
  for (i in 1:N_A) { # Loop though all samples of this type
    for (j in 1:N_T) { # Loop through all tracers for each sample
      X_A[i,j] ~ dnorm(mean_A[j], 1/sd_A[j]^2)
    }
  }
  
  # Source B
  for (i in 1:N_B) { # Loop though all samples of this type
    for (j in 1:N_T) { # Loop through all tracers for each sample
      X_B[i,j] ~ dnorm(mean_B[j], 1/sd_B[j]^2)
    }
  }
  
  /* uniform priors for tracer means in prey/sources */
  for (j in 1:N_T) {
    mean_A[j] ~ dunif(-100,100)
    mean_B[j] ~ dunif(-100,100)
  }
  /* approximate reciprocal priors for tracer sd in prey/sources */
  for (j in 1:N_T) {
    sd_A[j] ~ dgamma(0.001, 0.001) T(1*10^-50, 2)
    sd_B[j] ~ dgamma(0.001, 0.001) T(1*10^-50, 2)
  }
  
  /* deterministic for TDFs */
  # metazoan TDFs
  for (j in i_T_frac) {
    TDF_meta[j] <- TDF_m[j]
    sdTDF_meta[j] <- sdTDF_m[j]
  }
  # protozoan TDFs
  for (j in i_T_frac) {
    TDF_proto[j] <- TDF_p[j]
    sdTDF_proto[j] <- sdTDF_p[j]
  }
  
  /* priors for the zooplankton */
    for (iz in 1:N_Z) { # N_Z is the number of zooplankton
      pz[iz, 1:3] ~ ddirch(c(1,1,1)) # uniform priors
    }"
OMSM_Source_Data_and_Priors_3 <-
  "/* observe tracers in sources */
  # Each row will be a sample
  # Each column will be a tracer
    
  # Source A
  for (i in 1:N_A) { # Loop though all samples of this type
    for (j in 1:N_T) { # Loop through all tracers for each sample
      X_A[i,j] ~ dnorm(mean_A[j], 1/sd_A[j]^2)
    }
  }
  
  # Source B
  for (i in 1:N_B) { # Loop though all samples of this type
    for (j in 1:N_T) { # Loop through all tracers for each sample
      X_B[i,j] ~ dnorm(mean_B[j], 1/sd_B[j]^2)
    }
  }
  
  # Source C
  for (i in 1:N_C) { # Loop though all samples of this type
    for (j in 1:N_T) { # Loop through all tracers for each sample
      X_C[i,j] ~ dnorm(mean_C[j], 1/sd_C[j]^2)
    }
  }
  
  /* uniform priors for tracer means in prey/sources */
  for (j in 1:N_T) {
    mean_A[j] ~ dunif(-100,100)
    mean_B[j] ~ dunif(-100,100)
    mean_C[j] ~ dunif(-100,100)
  }
  /* approximate reciprocal priors for tracer sd in prey/sources */
  for (j in 1:N_T) {
    sd_A[j] ~ dgamma(0.001, 0.001) T(1*10^-50, 2)
    sd_B[j] ~ dgamma(0.001, 0.001) T(1*10^-50, 2)
    sd_C[j] ~ dgamma(0.001, 0.001) T(1*10^-50, 2)
  }
  
  /* deterministic for TDFs */
  # metazoan TDFs
  for (j in i_T_frac) {
    TDF_meta[j] <- TDF_m[j]
    sdTDF_meta[j] <- sdTDF_m[j]
  }
  # protozoan TDFs
  for (j in i_T_frac) {
    TDF_proto[j] <- TDF_p[j]
    sdTDF_proto[j] <- sdTDF_p[j]
  }
  
  /* priors for the zooplankton */
    for (iz in 1:N_Z) { # N_Z is the number of zooplankton
      pz[iz, 1:3] ~ ddirch(c(1,1,1)) # uniform priors
    }"
OMSM_Source_Data_and_Priors_4 <-
  "/* observe tracers in sources */
  # Each row will be a sample
  # Each column will be a tracer
    
  # Source A
  for (i in 1:N_A) { # Loop though all samples of this type
    for (j in 1:N_T) { # Loop through all tracers for each sample
      X_A[i,j] ~ dnorm(mean_A[j], 1/sd_A[j]^2)
    }
  }
  
  # Source B
  for (i in 1:N_B) { # Loop though all samples of this type
    for (j in 1:N_T) { # Loop through all tracers for each sample
      X_B[i,j] ~ dnorm(mean_B[j], 1/sd_B[j]^2)
    }
  }
  
  # Source C
  for (i in 1:N_C) { # Loop though all samples of this type
    for (j in 1:N_T) { # Loop through all tracers for each sample
      X_C[i,j] ~ dnorm(mean_C[j], 1/sd_C[j]^2)
    }
  }
  
  # Source D
  for (i in 1:N_D) { # Loop though all samples of this type
    for (j in 1:N_T) { # Loop through all tracers for each sample
      X_D[i,j] ~ dnorm(mean_D[j], 1/sd_D[j]^2)
    }
  }
  
  /* uniform priors for tracer means in prey/sources */
  for (j in 1:N_T) {
    mean_A[j] ~ dunif(-100,100)
    mean_B[j] ~ dunif(-100,100)
    mean_C[j] ~ dunif(-100,100)
    mean_D[j] ~ dunif(-100,100)
  }
  /* approximate reciprocal priors for tracer sd in prey/sources */
  for (j in 1:N_T) {
    sd_A[j] ~ dgamma(0.001, 0.001) T(1*10^-50, 2)
    sd_B[j] ~ dgamma(0.001, 0.001) T(1*10^-50, 2)
    sd_C[j] ~ dgamma(0.001, 0.001) T(1*10^-50, 2)
    sd_D[j] ~ dgamma(0.001, 0.001) T(1*10^-50, 2)
  }
  
  /* deterministic for TDFs */
  # metazoan TDFs
  for (j in i_T_frac) {
    TDF_meta[j] <- TDF_m[j]
    sdTDF_meta[j] <- sdTDF_m[j]
  }
  # protozoan TDFs
  for (j in i_T_frac) {
    TDF_proto[j] <- TDF_p[j]
    sdTDF_proto[j] <- sdTDF_p[j]
  }
  
  /* priors for the zooplankton */
    for (iz in 1:N_Z) { # N_Z is the number of zooplankton
      pz[iz, 1:4] ~ ddirch(c(1,1,1,1)) # uniform priors
    }"
OMSM_Source_Data_and_Priors_5 <-
  "/* observe tracers in sources */
  # Each row will be a sample
  # Each column will be a tracer
    
  # Source A
  for (i in 1:N_A) { # Loop though all samples of this type
    for (j in 1:N_T) { # Loop through all tracers for each sample
      X_A[i,j] ~ dnorm(mean_A[j], 1/sd_A[j]^2)
    }
  }
  
  # Source B
  for (i in 1:N_B) { # Loop though all samples of this type
    for (j in 1:N_T) { # Loop through all tracers for each sample
      X_B[i,j] ~ dnorm(mean_B[j], 1/sd_B[j]^2)
    }
  }
  
  # Source C
  for (i in 1:N_C) { # Loop though all samples of this type
    for (j in 1:N_T) { # Loop through all tracers for each sample
      X_C[i,j] ~ dnorm(mean_C[j], 1/sd_C[j]^2)
    }
  }
  
  # Source D
  for (i in 1:N_D) { # Loop though all samples of this type
    for (j in 1:N_T) { # Loop through all tracers for each sample
      X_D[i,j] ~ dnorm(mean_D[j], 1/sd_D[j]^2)
    }
  }
  
  # Source E
  for (i in 1:N_E) { # Loop though all samples of this type
    for (j in 1:N_T) { # Loop through all tracers for each sample
      X_E[i,j] ~ dnorm(mean_E[j], 1/sd_E[j]^2)
    }
  }
  
  /* uniform priors for tracer means in prey/sources */
  for (j in 1:N_T) {
    mean_A[j] ~ dunif(-100,100)
    mean_B[j] ~ dunif(-100,100)
    mean_C[j] ~ dunif(-100,100)
    mean_D[j] ~ dunif(-100,100)
    mean_E[j] ~ dunif(-100,100)
  }
  /* approximate reciprocal priors for tracer sd in prey/sources */
  for (j in 1:N_T) {
    sd_A[j] ~ dgamma(0.001, 0.001) T(1*10^-50, 2)
    sd_B[j] ~ dgamma(0.001, 0.001) T(1*10^-50, 2)
    sd_C[j] ~ dgamma(0.001, 0.001) T(1*10^-50, 2)
    sd_D[j] ~ dgamma(0.001, 0.001) T(1*10^-50, 2)
    sd_E[j] ~ dgamma(0.001, 0.001) T(1*10^-50, 2)
  }
  
  /* deterministic for TDFs */
  # metazoan TDFs
  for (j in i_T_frac) {
    TDF_meta[j] <- TDF_m[j]
    sdTDF_meta[j] <- sdTDF_m[j]
  }
  # protozoan TDFs
  for (j in i_T_frac) {
    TDF_proto[j] <- TDF_p[j]
    sdTDF_proto[j] <- sdTDF_p[j]
  }
  
  /* priors for the zooplankton */
    for (iz in 1:N_Z) { # N_Z is the number of zooplankton
      pz[iz, 1:5] ~ ddirch(c(1,1,1,1,1)) # uniform priors
    }"
OMSM_Source_Data_and_Priors_6 <-
  "/* observe tracers in sources */
  # Each row will be a sample
  # Each column will be a tracer
    
  # Source A
  for (i in 1:N_A) { # Loop though all samples of this type
    for (j in 1:N_T) { # Loop through all tracers for each sample
      X_A[i,j] ~ dnorm(mean_A[j], 1/sd_A[j]^2)
    }
  }
  
  # Source B
  for (i in 1:N_B) { # Loop though all samples of this type
    for (j in 1:N_T) { # Loop through all tracers for each sample
      X_B[i,j] ~ dnorm(mean_B[j], 1/sd_B[j]^2)
    }
  }
  
  # Source C
  for (i in 1:N_C) { # Loop though all samples of this type
    for (j in 1:N_T) { # Loop through all tracers for each sample
      X_C[i,j] ~ dnorm(mean_C[j], 1/sd_C[j]^2)
    }
  }
  
  # Source D
  for (i in 1:N_D) { # Loop though all samples of this type
    for (j in 1:N_T) { # Loop through all tracers for each sample
      X_D[i,j] ~ dnorm(mean_D[j], 1/sd_D[j]^2)
    }
  }
  
  # Source E
  for (i in 1:N_E) { # Loop though all samples of this type
    for (j in 1:N_T) { # Loop through all tracers for each sample
      X_E[i,j] ~ dnorm(mean_E[j], 1/sd_E[j]^2)
    }
  }
  
  # Source F
  for (i in 1:N_F) { # Loop though all samples of this type
    for (j in 1:N_T) { # Loop through all tracers for each sample
      X_F[i,j] ~ dnorm(mean_F[j], 1/sd_F[j]^2)
    }
  }
  
  /* uniform priors for tracer means in prey/sources */
  for (j in 1:N_T) {
    mean_A[j] ~ dunif(-100,100)
    mean_B[j] ~ dunif(-100,100)
    mean_C[j] ~ dunif(-100,100)
    mean_D[j] ~ dunif(-100,100)
    mean_E[j] ~ dunif(-100,100)
    mean_F[j] ~ dunif(-100,100)
  }
  /* approximate reciprocal priors for tracer sd in prey/sources */
  for (j in 1:N_T) {
    sd_A[j] ~ dgamma(0.001, 0.001) T(1*10^-50, 2)
    sd_B[j] ~ dgamma(0.001, 0.001) T(1*10^-50, 2)
    sd_C[j] ~ dgamma(0.001, 0.001) T(1*10^-50, 2)
    sd_D[j] ~ dgamma(0.001, 0.001) T(1*10^-50, 2)
    sd_E[j] ~ dgamma(0.001, 0.001) T(1*10^-50, 2)
    sd_F[j] ~ dgamma(0.001, 0.001) T(1*10^-50, 2)
  }
  
  /* deterministic for TDFs */
  # metazoan TDFs
  for (j in i_T_frac) {
    TDF_meta[j] <- TDF_m[j]
    sdTDF_meta[j] <- sdTDF_m[j]
  }
  # protozoan TDFs
  for (j in i_T_frac) {
    TDF_proto[j] <- TDF_p[j]
    sdTDF_proto[j] <- sdTDF_p[j]
  }
  
  /* priors for the zooplankton */
    for (iz in 1:N_Z) { # N_Z is the number of zooplankton
      pz[iz, 1:6] ~ ddirch(c(1,1,1,1,1,1)) # uniform priors
    }"

OMSM_Mixing_2 <-
  "  /* calculate tracers at the base of the food web and then in zooplankton */
    for (iz in 1:N_Z) {
      
      ## Carrying out organic matter source mixing at the base of the food web
      for (j in 1:N_T) {
        mean_b[iz,j] <- pz[iz,1]*mean_A[j] + pz[iz,2]*mean_B[j]
      }"
OMSM_Mixing_3 <-
  "  /* calculate tracers at the base of the food web and then in zooplankton */
    for (iz in 1:N_Z) {
      
      ## Carrying out organic matter source mixing at the base of the food web
      for (j in 1:N_T) {
        mean_b[iz,j] <- pz[iz,1]*mean_A[j] + pz[iz,2]*mean_B[j] + pz[iz,3]*mean_C[j]
      }"
OMSM_Mixing_4 <-
  "  /* calculate tracers at the base of the food web and then in zooplankton */
    for (iz in 1:N_Z) {
      
      ## Carrying out organic matter source mixing at the base of the food web
      for (j in 1:N_T) {
        mean_b[iz,j] <- pz[iz,1]*mean_A[j] + pz[iz,2]*mean_B[j] + pz[iz,3]*mean_C[j] + pz[iz,4]*mean_D[j]
      }"
OMSM_Mixing_5 <-
  "  /* calculate tracers at the base of the food web and then in zooplankton */
    for (iz in 1:N_Z) {
      
      ## Carrying out organic matter source mixing at the base of the food web
      for (j in 1:N_T) {
        mean_b[iz,j] <- pz[iz,1]*mean_A[j] + pz[iz,2]*mean_B[j] + pz[iz,3]*mean_C[j] + 
                        pz[iz,4]*mean_D[j] + pz[iz,5]*mean_E[j]
      }"
OMSM_Mixing_6 <-
  "  /* calculate tracers at the base of the food web and then in zooplankton */
    for (iz in 1:N_Z) {
      
      ## Carrying out organic matter source mixing at the base of the food web
      for (j in 1:N_T) {
        mean_b[iz,j] <- pz[iz,1]*mean_A[j] + pz[iz,2]*mean_B[j] + pz[iz,3]*mean_C[j] + 
                        pz[iz,4]*mean_D[j] + pz[iz,5]*mean_E[j] + pz[iz,6]*mean_F[j]
      }"

OMSM_Discrimination_Conservative <-
  "      ## Non-fractionating AAs
      for (j in i_T_non) {
        mean_z[iz,j] <- mean_b[iz, j]
        va_z[iz,j]   <- sd_Y[iz,j]^2
      }"
OMSM_Discrimination_Constant <-
  "      ## AAs with constant TDFs
      for (j in i_T_const) {
        mean_z[iz,j]  <- mean_b[iz, j] + FWL[iz]*TDF_meta[j]
        va_z[iz,j]    <- sd_Y[iz,j]^2 + 
                         (FWL[iz]*TDF_meta[j] * (sdTDF_meta[j]*TDF_meta[j]^-1))^2
      }"
OMSM_Discrimination_Variable <-
"      ## AAs with variable TDFs
      for (j in i_T_var) {
        mean_z[iz,j]  <- mean_b[iz, j]  + PTS[iz]*TDF_proto[j] + MTS[iz]*TDF_meta[j]
        va_z[iz,j]    <- sd_Y[iz,j]^2 + 
                         (PTS[iz]*TDF_proto[j] * (sdTDF_proto[j]*TDF_proto[j]^-1))^2 +
                         (MTS[iz]*TDF_meta[j] * (sdTDF_meta[j]*TDF_meta[j]^-1))^2
      }"

OMSM_Trophic_Parameters <-
  "      ## Carrying out trophic enrichment from food web base to zooplankton
      ## calculating FWL from MTS and PTS
      mean_FWL[iz] <- (Y[iz,i_T_FWL] - mean_b[iz,i_T_FWL]) * TDF_meta[i_T_FWL]^-1
      va_FWL[iz] <- mean_FWL[iz] *
                    (sd_Y[iz,i_T_FWL] * Y[iz,i_T_FWL]^-1)^2 + 
                    (sdTDF_meta[i_T_FWL] * TDF_meta[i_T_FWL]^-1)^2
      mean_MTS[iz] <- (Y[iz,i_T_MTS] - mean_b[iz,i_T_MTS]) * TDF_meta[i_T_MTS]^-1
      va_MTS[iz] <- mean_MTS[iz] *
                    (sd_Y[iz,i_T_MTS] * Y[iz,i_T_MTS]^-1)^2 + 
                    (sdTDF_meta[i_T_MTS] * TDF_meta[i_T_MTS]^-1)^2
      FWL[iz] ~ dnorm(mean_FWL[iz], va_FWL[iz]^-1) T(1,)
      MTS[iz] ~ dnorm(mean_MTS[iz], va_MTS[iz]^-1) T(0,)
      PTS[iz] <- FWL[iz] - MTS[iz]
"

OMSM_Observing_Consumers <-
  "      ## tracer observations
      for (j in i_T_mix) {
        Y[iz,j] ~ dnorm(mean_z[iz,j] , 1/va_z[iz,j])
      }
    }
}"

# Pasting at the chunks together that we'll need for the model.
OMSM_Full <-
  paste(
    "model{",
    eval(as.name(paste(
      "OMSM_Source_Data_and_Priors",
      N_S,sep="_"))),
    eval(as.name(paste(
      "OMSM_Mixing",
      N_S,sep="_"))),
    OMSM_Trophic_Parameters,
    if(Data.OMSM$i_T_non[1]==0){"  "}
    else{OMSM_Discrimination_Conservative},
    if(Data.OMSM$i_T_const[1]==0){"  "}
    else{OMSM_Discrimination_Constant},
    if(Data.OMSM$i_T_var[1]==0){"  "}
    else{OMSM_Discrimination_Variable},
    OMSM_Observing_Consumers,
    sep = " \n "
  )

return(OMSM_Full)
}