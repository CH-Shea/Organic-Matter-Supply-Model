Post_Error_Quant_Trophic <- function(
    zoops.f,
    posts
){
  data.true <- zoops.f[c("PTS","MTS","FWL")]
  
  data.model <- posts$trophic$mode[c("PTS","MTS","FWL")]
  
  troph.disc <- data.model-data.true
  
  
  troph.disc.abs <- abs(troph.disc[c("PTS","MTS","FWL")])
  
  # troph.disc.sum <- aggregate(.~Model,
  #                             data=troph.disc.abs,
  #                             FUN=sum)
  # troph.disc.sum
  
  troph.disc.mean <-colMeans(troph.disc.abs)
  print("Mean error rates associated with each trohic parameter:")
  print(troph.disc.mean)
  
  for (i in 1:3) {
    print(paste("maximum error -", c("PTS","MTS","FWL")[i]))
    print(max(troph.disc.abs[[i]]))
  }
}




Post_Error_Quant_f <- function(
    zoops.f,
    posts
){
  data.true <- zoops.f[Sources]
  data.model <- posts$f$mode[Sources]
  
  mix.disc <- data.model-data.true
  mix.disc.abs <- abs(mix.disc[Sources])
  
  # mix.disc.sum <- colSums(mix.disc.abs)
  # print("Summed Absolute Discrepancy:")
  # mix.disc.sum
  
  mix.disc.mean <- colMeans(mix.disc.abs)
  print("Mean error rates associated with each mixing parameter:")
  print(mix.disc.mean)
  
  for (i in Sources) {
    print(paste("maximum error -", i, "particles:"))
    print(max(mix.disc.abs[[i]]))
  }
}