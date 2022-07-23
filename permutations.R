plot_gbi_null_custom <- function(x, ylim=c(0, 1)){
  plot(x$mcmc)
  par(mfrow = c(which.min(c(4,length(x$observed))),1))
  for(i in 1:length(x$observed)){
    for(k in 1:length(x$mcmc)){
      pval_gr <- cumsum(x$mcmc[[k]][,i] > x$observed[i])/(1:length(x$mcmc[[k]][,i]))
      if(k == 1){
        plot(pval_gr, main = names(x$observed)[i], type = "l", col = k, xlab = "Iteration", ylab = "P(Random > Observed)", ylim = ylim)
        abline(h=0.05, lty=2)
      }else{
        points(pval_gr, type = "l", col = k)
      }
    }
  }
}

gbi_MCMC_custom <- function (data, ind_constraint = NULL, group_constraint = NULL,
          samples = 1000, thin = 100, burnin = 1000, chains = 2, FUN = NULL)
{
  if (!is.matrix(data) | any(!c(data) %in% c(0, 1)) | any(is.na(data))) {
    stop("Data must be a matrix containing only 1s and 0s")
  }
  N <- ncol(data)
  G <- nrow(data)
  if (is.null(group_constraint))
    group_constraint <- rep(1, G)
  if (is.null(ind_constraint))
    ind_constraint <- rep(1, N)
  if (!is.vector(group_constraint) | length(group_constraint) !=
      G | any(is.na(group_constraint)))
    stop("group_constraint must be a vector with length equal to the number of groups")
  if (!is.vector(ind_constraint) | length(ind_constraint) !=
      N | any(is.na(ind_constraint)))
    stop("ind_constraint must be a vector with length equal to the number of individuals")
  if (is.null(FUN)) {
    FUN <- function(gbi) {
      x <- get_numerator(gbi, data_format = "GBI", return = "vector")
      d <- get_denominator(gbi, data_format = "GBI", return = "vector")
      sri <- x/d
      res <- c(mean(sri, na.rm = T), sd(sri, na.rm = T),
               sd(sri, na.rm = T)/mean(sri, na.rm = T), mean(sri >
                                                               0, na.rm = T))
      names(res) <- c("Mean", "SD", "CV", "Non-zero")
      return(res)
    }
  }
  if (!is.function(FUN))
    stop("FUN must be a function")
  observed <- FUN(data)
  chain_res <- list()
  chain_res <- parallel::mclapply(1:chains, function(k) {
    gbi.p <- data
    res_matrix <- matrix(nrow = samples, ncol = length(observed))
    colnames(res_matrix) <- names(observed)
    for (i in 1:burnin) {
      if (i %% 1000 == 0) {
        cat(paste0("Chain ", k, " Burn-in iteration: ", i, "\n"))
      }
      for (j in 1:thin) {
        cols <- sample(N, 2)
        rows <- sample(G, 2)
        trial_matrix <- gbi.p[rows, cols]
        if (all(rowSums(trial_matrix) == 1) & all(colSums(trial_matrix) ==
                                                  1) & group_constraint[rows[1]] == group_constraint[rows[2]] &
            ind_constraint[cols[1]] == ind_constraint[cols[2]]) {
          trial_matrix <- ifelse(trial_matrix == 1,
                                 0, 1)
          gbi.p[rows, cols] <- trial_matrix
        }
      }
    }
    for (i in 1:samples) {
      if (i %% 1000 == 0) {
        cat(paste0("Chain ", k, " Sampling iteration: ", i, "\n"))
      }
      for (j in 1:thin) {
        cols <- sample(N, 2)
        rows <- sample(G, 2)
        trial_matrix <- gbi.p[rows, cols]
        if (all(rowSums(trial_matrix) == 1) & all(colSums(trial_matrix) ==
                                                  1) & group_constraint[rows[1]] == group_constraint[rows[2]] &
            ind_constraint[cols[1]] == ind_constraint[cols[2]]) {
          trial_matrix <- ifelse(trial_matrix == 1,
                                 0, 1)
          gbi.p[rows, cols] <- trial_matrix
        }
      }
      res_matrix[i, ] <- FUN(gbi.p)
    }
    return(coda::mcmc(res_matrix, thin = thin, start = thin, end = thin * samples))
  },
  mc.cores=10
  )
  chain_res <- coda::mcmc.list(chain_res)
  results <- list(FUN = FUN, observed = observed, mcmc = chain_res)
  class(results) <- "gbi_null"
  return(results)
}
