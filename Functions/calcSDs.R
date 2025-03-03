# defining function to propagate uncertainty for mean calculations
SDmean <- function(data, na.rm = FALSE)
{
  if (prod(!is.na(data)) == 1) {
    squareds = data^2
    sumofsquareds = sum(squareds)
    root = sqrt(sumofsquareds)
    SD = root/length(as.double(data))
    return(SD)
  } else { if (na.rm == TRUE) {
    data = na.omit(data)
    squareds = data^2
    sumofsquareds = sum(squareds)
    root = sqrt(sumofsquareds)
    SD = root/length(data)
    return(SD)
  } else {
    return(NA)
  }
  }
}

# defining a function to propagate error through summing calculations
SDsum <- function(data)
{
  if (is.numeric(data) == TRUE) {
    squareds = data^2
    sumofsquareds = sum(squareds)
    root = sqrt(sumofsquareds)
    SD = root
    return(SD)
  } else {
    return(NA)
  }
}