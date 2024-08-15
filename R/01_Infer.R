### Function 1. Learning process
#' @title Infer
#' @description It provides the theoretical probability distribution of the
#' types given the control (or any) dataset and a non-informative prior
#' @param data: A data.frame with the known data where the first column is the
#'  assigned period and successive columns are the ordered types
#' @param periods: A factor where the name of each value is the site/level
#' and the value is the period to which the site/level corresponds.
#' @param prior: Prior to be selected. Options are "Haldane", "Jeffreys", "Perk",
#' "BayesLaplace","Hierarchical". Default is "Perk".
#' @returns A list where each element is the period. Each element/period is a data frame
#' with the curve for each type
#' @export

Infer <- function(data, prior = "Perk"){

  ## Standardise name for first column
  colnames(data)[1] <- "Period"

  # How many periods
  P <- length(levels(as.factor(data$Period)))

  # How many types
  Ty <- ncol(data)-1

  ## Vector with arrow type to create data frames
  svaria <- paste0("Type.",seq(1,Ty))

  ## Select prior
  if (prior == "Haldane"){
    teo_distribution <- 0
  } else if (prior == "Jeffreys"){
    teo_distribution <- 1/2
  } else if (prior == "Perk"){
    teo_distribution <- 1/Ty
  } else if (prior == "BayesLaplace"){
    teo_distribution <- 1
  } else if (prior == "Hierarchical"){
    teo_distribution <- sqrt(2)/Ty
  }

  # Utility to help with the beta distribution
  x <- seq(0,1,length=1000)

  # To store the results
  periods_train <- list()

  for (i in 1:P){
    newdata <- subset(data,Period == i, svaria)

    sum_newdata <- list()
    ys <- list()

    for (j in 1:Ty){
      ## Data frames with the sum of type [j] in period [i]
      sum_newdata[[j]] <- sum(newdata[,j])
      ## Data frame with the proportions of type [j] in period [i]
      ys[[j]] <- MASS::fractions(sum(newdata[,j])+teo_distribution)
    }

    # create a new data.frame with the sum of the rest of type[j] in period[i] (result in proportions)
    ybs <- list()
    ytot <- sum(unlist(ys))

    for (j in 1:Ty){
      ybs[[j]] <- ytot-ys[[j]]
    }

    # create a new data.frame with marginal distribution for type[j] in period[i]
    type <- list()
    for (j in 1:Ty){
      type[[j]] <- dbeta(x, ys[[j]], ybs[[j]])
    }

    # create a vector with the combination of all marginal distributions
    beta_dist <- as.data.frame(matrix(unlist(type), ncol = Ty))
    colnames(beta_dist) <- paste0("Type.", 1:Ty)

    periods_train[[i]] <- beta_dist
    names(periods_train)[[i]] <- paste("period_",i)

  }
  return(periods_train)
}
