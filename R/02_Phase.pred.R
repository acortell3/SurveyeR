### Function 2. Predict specific site
#' @title Phase.pred
#' @description It performs the Dirichlet process and returns the chronological
#' probability of an assemblage or site, given the control dataset
#' @param dated: A data.frame with the known data where the first column is the
#'  assigned period and successive columns are the ordered types
#' @param estimate: Data on which we want to estimate. Rownames are the names of
#' the sites and are used to compute the probability, so they must be correct.
#' Columns are number of each type per site
#' @param prior: Prior to be selected. Options are "Haldane", "Jeffreys", "Perk",
#' "BayesLaplace","Hierarchical". Default is "Perk".
#' @param site: The site on which we want to estimate. It must be spelled as in
#' the rownames
#' @returns a vector with the prediction for the selected site where each position
#' is the probability per period
#' @export

Phase.pred <- function(dated, estimate, prior = "Perk", site){

  ## Define periods

  ## Standardise name for first column
  colnames(dated)[1] <- "Period"

  ## Total dates sites
  tds <- nrow(dated)

  ## Total number of periods
  Ps <- length(levels(as.factor(dated$Period)))

  ## A vector where each value is the number of sites in that period
  prs <- c(rep(NA,Ps))
  for (i in 1:length(prs)){
    prs[i] <- nrow(subset(dated,Period == i))
  }

  # Defines P(m* = m[i] | D)
  est.p <- c(prs/tds)

  ## Switch to types

  ## Vector with arrow type to create data frames
  svaria <- colnames(dated)[-1]

  ## List with subsets of types per period
  post <- list()
  for (i in 1:length(prs)){
    post[[i]] <- subset(dated, Period == i, svaria)
  }

  ## Select prior
  if (prior == "Haldane"){
    teo_distribution <- 0
  } else if (prior == "Jeffreys"){
    teo_distribution <- 1/2
  } else if (prior == "Perk"){
    teo_distribution <- 1/length(svaria)
  } else if (prior == "BayesLaplace"){
    teo_distribution <- 1
  } else if (prior == "Hierarchical"){
    teo_distribution <- sqrt(2)/length(svaria)
  }

  ## Compute Dirichlet probabilities
  types <- list()
  y <- list()
  Dir.res <- list()

  for (j in 1:length(post)){
    for (i in 1:length(svaria)){
      types[[i]] <- sum(post[[j]][,i])
    }

    for (i in 1:length(svaria)){
      y[[i]] <- MASS::fractions(types[[i]] + teo_distribution)
    }

    Dir.res[[j]] <- unlist(y)
  }

  ## Select the site to observe
  y.star <- subset(estimate[,-1], rownames(estimate) == site)
  y.star <- y.star
  n.star <-sum(y.star)

  ## Compute the likelihood
  likelihood <- rep(NA,length(post))

  for (i in 1:length(post)){
    likelihood[i] <- extraDistr::ddirmnom(y.star, n.star, as.numeric(Dir.res[[i]]))
  }

  ## Results
  res <- likelihood*est.p /sum(likelihood*est.p)
  names(res) <- paste0("Period.", 1:length(res))
  return(res)

}
