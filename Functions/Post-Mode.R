# the mode is more diagnostic than the mean for the posteriors, so we're going to go ahead and write a function to find the point of maximum probability for a PDF
post.mode <- function(post) {
  # first calculate the density function for the posterior PDF
  # this bins the data in n bins
  post.d <- density(post, n=512*2)
  # then find the index of the bin with the highest density
  bin.index <-with(post.d, which.max(y))
  # then find the variable value of the PDF in that bin
  with(post.d, x[bin.index])
}