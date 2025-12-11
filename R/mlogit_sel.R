#' Lambda2 Function - Bivariate Normal Copula
#'
#' Computes the bivariate normal copula joint probability
#'
#' @param u First probability (must be in (0,1))
#' @param v Second probability (must be in (0,1))
#' @param rho Correlation parameter (must be in (-1,1))
#' @return Joint probability from bivariate normal copula
#' @export
Lambda2 <- function(u, v, rho) {
  # Bound inputs to avoid numerical issues
  u <- pmax(pmin(u, 1 - 1e-10), 1e-10)
  v <- pmax(pmin(v, 1 - 1e-10), 1e-10)
  rho <- pmin(pmax(rho, -0.999999), 0.999999)
  
  # Transform to standard normal quantiles
  #u_inv <- qnorm(u)
  #v_inv <- qnorm(v)
  
  # Bivariate normal CDF with correlation rho
  # Using pbivnorm or approximation
  normal_cop <- copula::normalCopula(param = rho, dim = 2)
  bvn_cdf <- copula::pCopula(c(u, v), normal_cop)

  return(bvn_cdf)
}

#' Multinomial Probabilities
#'
#' Computes multinomial logit probabilities
#'
#' @param X Matrix of covariates
#' @param beta_list List of coefficient vectors for each category (except baseline)
#' @return Matrix of probabilities
#' @export
multinomial_probs <- function(X, beta_list) {
  n <- nrow(X)
  q_minus_1 <- length(beta_list)
  q <- q_minus_1 + 1
  eta <- sapply(beta_list, function(b) X %*% b)
  eta <- cbind(eta, rep(0, n))
  eta_max <- apply(eta, 1, max)
  exp_eta <- exp(eta - eta_max)
  pi <- exp_eta / rowSums(exp_eta)
  return(pi)
}

#' Multinomial Logit Model with Sample Selection
#'
#' Estimates a multinomial logit model with sample selection correction using
#' a two-step approach. First estimates the selection equation, then estimates
#' the outcome equation accounting for selection bias through a copula-based
#' dependence structure.
#'
#' @param Y Numeric vector of outcomes (categorical variable)
#' @param S Binary vector indicating selection (1 = observed, 0 = not observed)
#' @param X Matrix of covariates for the outcome equation
#' @param Z Matrix of additional covariates for the selection equation
#' @param baseline_last Logical, if TRUE the last category is baseline (default TRUE)
#' @param init_beta Optional list of initial values for beta parameters
#' @param init_gamma Optional list of initial values for gamma parameters
#' @return A list containing:
#'   \item{delta}{Selection equation coefficients}
#'   \item{beta}{List of outcome equation coefficients for each category}
#'   \item{gamma}{List of dependence parameters for each category}
#'   \item{q}{Number of outcome categories}
#'   \item{convergence}{Convergence code from optimization (0 = success)}
#'   \item{logLik}{Log-likelihood value at convergence}
#'   \item{theta_hat}{Vector of all estimated parameters (beta, gamma)}
#'   \item{nu_hat}{Fitted selection probabilities}
#'   \item{Y_mapped}{Mapped outcome values}
#'   \item{X}{Covariate matrix (stored for bootstrap)}
#'   \item{S}{Selection vector (stored for bootstrap)}
#' @examples
#' \dontrun{
#' # Generate example data
#' set.seed(123)
#' n <- 500
#' X <- cbind(intercept = 1, x1 = rnorm(n), x2 = rnorm(n))
#' Z <- cbind(z1 = rnorm(n), z2 = rnorm(n))
#' W <- cbind(X, Z)
#' 
#' # Generate selection
#' delta_true <- c(0.5, 0.02, 0.1, -0.3, 0.05)
#' prob_select <- plogis(W %*% delta_true)
#' S <- rbinom(n, 1, prob_select)
#' 
#' # Generate outcomes
#' eta1 <- X %*% c(-1, 0.05, 0.2)
#' eta2 <- X %*% c(0.5, -0.03, 0.15)
#' probs <- exp(cbind(eta1, eta2, 0))
#' probs <- probs / rowSums(probs)
#' Y <- apply(probs, 1, function(p) sample(0:2, 1, prob = p))
#' Y[S == 0] <- NA
#' 
#' # Fit model
#' result <- mlogit_sel(Y = Y, S = S, X = X, Z = Z)
#' print(result)
#' }
#' @export
mlogit_sel <- function(Y, S, X, Z, 
                       baseline_last = TRUE,
                       init_beta = NULL,
                       init_gamma = NULL) {
  n <- length(Y)
  dx <- ncol(X)
  W <- cbind(X, Z)
  
  # Prepare outcome coding
  Y_obs_vals <- Y[S == 1]
  unique_y <- sort(unique(Y_obs_vals))
  q <- length(unique_y)
  
  if (baseline_last) {
    category_map <- setNames(1:q, unique_y)
    Y_mapped <- ifelse(S == 1, category_map[as.character(Y)], NA)
  } else {
    Y_mapped <- ifelse(S == 1, Y + 1, NA)
  }
  
  # STEP 1: Selection equation
  fit_sel <- stats::glm(S ~ W - 1, family = stats::binomial(link = "logit"))
  delta_hat <- stats::coef(fit_sel)
  nu_hat <- stats::fitted(fit_sel)  # P(S=1 | W)
  
  # STEP 2: Outcome & dependence parameters — handle initial values
  num_beta_params <- (q - 1) * dx
  
  # --- Handle init_beta ---
  if (is.null(init_beta)) {
    beta0 <- rep(0, num_beta_params)
  } else {
    if (!is.list(init_beta) || length(init_beta) != (q - 1)) {
      stop("'init_beta' must be a list of length ", q - 1)
    }
    if (!all(sapply(init_beta, function(b) length(b) == dx))) {
      stop("Each element of 'init_beta' must be a numeric vector of length ", dx)
    }
    beta0 <- unlist(init_beta)
  }
  
  # --- Handle init_gamma ---
  if (is.null(init_gamma)) {
    gamma0 <- rep(0, num_beta_params)
  } else {
    if (!is.list(init_gamma) || length(init_gamma) != (q - 1)) {
      stop("'init_gamma' must be a list of length ", q - 1)
    }
    if (!all(sapply(init_gamma, function(g) length(g) == dx))) {
      stop("Each element of 'init_gamma' must be a numeric vector of length ", dx)
    }
    gamma0 <- unlist(init_gamma)
  }
  
  init_params <- c(beta0, gamma0)
  
  # Likelihood function
  negloglik_step2 <- function(params, X, nu_hat, Y_mapped, S, q, dx) {
    n <- nrow(X)
    num_beta <- (q - 1) * dx
    beta_flat <- params[1:num_beta]
    gamma_flat <- params[(num_beta + 1):length(params)]
    
    beta_list <- split(beta_flat, rep(1:(q - 1), each = dx))
    gamma_list <- split(gamma_flat, rep(1:(q - 1), each = dx))
    
    pi_mat <- multinomial_probs(X, beta_list)
    rho_mat <- sapply(gamma_list, function(g) tanh(X %*% g))
    rho_mat <- cbind(rho_mat, rep(0, n))
    
    sel_idx <- which(S == 1)
    Y_sel <- Y_mapped[sel_idx]
    pi_sel <- pi_mat[sel_idx, , drop = FALSE]
    nu_sel <- nu_hat[sel_idx]
    rho_sel <- rho_mat[sel_idx, , drop = FALSE]
    
    loglik <- 0
    m <- length(Y_sel)
    
    for (i in 1:m) {
      k <- Y_sel[i]
      u_ik <- pi_sel[i, k]
      v_i  <- nu_sel[i]
      rho_ik <- rho_sel[i, k]
      
      joint_prob <- Lambda2(u_ik, v_i, rho_ik)
      
      if (joint_prob <= 0 || is.na(joint_prob)) {
        return(1e10)
      }
      loglik <- loglik + log(joint_prob)
    }
    return(-loglik)
  }
  
  # Optimization with logging
  cat("Starting Step 2: Estimating outcome and dependence parameters...\n")
  cat("  - Number of parameters (beta + gamma):", length(init_params), "\n")
  cat("  - Initial negative log-likelihood: ", 
      round(negloglik_step2(init_params, X, nu_hat, Y_mapped, S, q, dx), 2), 
      "\n", sep = "")
  
  opt <- stats::optim(
    par = init_params,
    fn = negloglik_step2,
    X = X,
    nu_hat = nu_hat,
    Y_mapped = Y_mapped,
    S = S,
    q = q,
    dx = dx,
    method = "BFGS",
    control = list(trace = 1, maxit = 8000)
  )
  
  cat("\nStep 2 optimization completed.\n")
  cat("  - Convergence code:", opt$convergence, "\n")
  cat("  - Final negative log-likelihood:", round(opt$value, 4), "\n")
  if (opt$convergence == 0) {
    cat("  - Status: Successful convergence\n")
  } else {
    cat("  - Status: Warning - check convergence\n")
  }
  
  if (opt$convergence != 0) warning("Optimization may not have converged")
  
  # Extract results
  num_beta <- (q - 1) * dx
  beta_flat <- opt$par[1:num_beta]
  gamma_flat <- opt$par[(num_beta + 1):length(opt$par)]
  
  beta_list <- split(beta_flat, rep(1:(q - 1), each = dx))
  gamma_list <- split(gamma_flat, rep(1:(q - 1), each = dx))
  
  result <- list(
    delta = delta_hat,
    beta = beta_list,
    gamma = gamma_list,
    q = q,
    convergence = opt$convergence,
    logLik = -opt$value,
    theta_hat = opt$par,
    nu_hat = nu_hat,
    Y_mapped = Y_mapped,
    X = X,
    S = S,
    dx = dx
  )
  
  class(result) <- c("mlogit_sel", "list")
  return(result)
}

#' Multiplier Bootstrap for Multinomial Logit with Sample Selection
#'
#' Implements multiplier bootstrap inference for the two-step multinomial logit
#' model with sample selection. Computes bootstrap standard errors and confidence
#' intervals using the multiplier bootstrap procedure.
#'
#' @param model_fit Output from mlogit_sel function
#' @param B Number of bootstrap replications (default 1000)
#' @param confidence_level Confidence level for intervals (default 0.95)
#' @param regularization Regularization parameter lambda as proportion of max eigenvalue (default 1e-6)
#' @param verbose Print progress messages (default TRUE)
#' @return A list containing:
#'   \item{theta_hat}{Original parameter estimates}
#'   \item{se}{Bootstrap standard errors}
#'   \item{ci_lower}{Lower confidence interval bounds}
#'   \item{ci_upper}{Upper confidence interval bounds}
#'   \item{bootstrap_draws}{Matrix of bootstrap parameter draws (B x d_theta)}
#'   \item{vcov}{Estimated variance-covariance matrix}
#'   \item{param_names}{Names of parameters}
#' @examples
#' \dontrun{
#' # Fit model
#' result <- mlogit_sel(Y = Y, S = S, X = X, Z = Z)
#' 
#' # Compute bootstrap inference
#' boot_result <- bootstrap_inference(result, B = 1000)
#' 
#' # View results
#' print(boot_result$se)
#' print(boot_result$ci_lower)
#' print(boot_result$ci_upper)
#' }
#' @export
bootstrap_inference <- function(model_fit, 
                                B = 1000, 
                                confidence_level = 0.95,
                                regularization = 1e-6,
                                verbose = TRUE) {
  
  if (verbose) {
    cat("\n============================================================================\n")
    cat("Multiplier Bootstrap Inference\n")
    cat("============================================================================\n\n")
  }
  
  # Extract fitted objects
  theta_hat <- model_fit$theta_hat
  X <- model_fit$X
  S <- model_fit$S
  Y_mapped <- model_fit$Y_mapped
  nu_hat <- model_fit$nu_hat
  q <- model_fit$q
  dx <- model_fit$dx
  
  n <- nrow(X)
  d_theta <- length(theta_hat)
  num_beta <- (q - 1) * dx
  
  # Create parameter names
  beta_names <- paste0("beta", rep(1:(q-1), each = dx), "_", 
                       rep(1:dx, times = q-1))
  gamma_names <- paste0("gamma", rep(1:(q-1), each = dx), "_", 
                        rep(1:dx, times = q-1))
  param_names <- c(beta_names, gamma_names)
  
  if (verbose) {
    cat("Computing individual log-likelihood contributions...\n")
  }
  
  # Define individual log-likelihood function
  loglik_i <- function(i, params) {
    if (S[i] == 0) return(0)
    
    beta_flat <- params[1:num_beta]
    gamma_flat <- params[(num_beta + 1):length(params)]
    
    beta_list <- split(beta_flat, rep(1:(q - 1), each = dx))
    gamma_list <- split(gamma_flat, rep(1:(q - 1), each = dx))
    
    # Compute probabilities for this observation
    pi_vec <- multinomial_probs(matrix(X[i, ], nrow = 1), beta_list)[1, ]
    rho_vec <- sapply(gamma_list, function(g) tanh(sum(X[i, ] * g)))
    rho_vec <- c(rho_vec, 0)
    
    k <- Y_mapped[i]
    u_k <- pi_vec[k]
    v_i <- nu_hat[i]
    rho_k <- rho_vec[k]
    
    joint_prob <- Lambda2(u_k, v_i, rho_k)
    
    if (joint_prob <= 0 || is.na(joint_prob)) {
      return(-1e10)
    }
    
    return(log(joint_prob))
  }
  
  # Step 1: Compute score vectors for each observation
  if (verbose) {
    cat("Step 1: Computing score vectors for each observation...\n")
  }
  
  scores <- matrix(0, nrow = n, ncol = d_theta)
  
  for (i in 1:n) {
    if (S[i] == 1) {
      scores[i, ] <- numDeriv::grad(function(theta) loglik_i(i, theta), theta_hat)
    }
  }
  
  # Step 2: Estimate observed information matrix
  if (verbose) {
    cat("Step 2: Computing observed information matrix...\n")
  }
  
  # Compute Hessian for each observation
  hessian_sum <- matrix(0, nrow = d_theta, ncol = d_theta)
  
  for (i in 1:n) {
    if (S[i] == 1) {
      H_i <- numDeriv::hessian(function(theta) loglik_i(i, theta), theta_hat)
      hessian_sum <- hessian_sum + H_i
    }
  }
  
  I_n <- -hessian_sum / n
  
  # Regularization
  eigenvals <- eigen(I_n, symmetric = TRUE, only.values = TRUE)$values
  lambda_max <- max(eigenvals)
  lambda <- regularization * lambda_max
  
  I_n_reg <- I_n + lambda * diag(d_theta)
  
  # Invert to get covariance estimate
  Sigma_hat <- tryCatch({
    solve(I_n_reg)
  }, error = function(e) {
    warning("Information matrix inversion failed. Using pseudoinverse.")
    MASS::ginv(I_n_reg)
  })
  
  if (verbose) {
    cat("  - Regularization parameter lambda:", format(lambda, scientific = TRUE), "\n")
    cat("  - Matrix condition number:", format(kappa(I_n_reg), scientific = TRUE), "\n\n")
  }
  
  # Step 3: Bootstrap
  if (verbose) {
    cat("Step 3: Running", B, "bootstrap replications...\n")
  }
  
  bootstrap_draws <- matrix(0, nrow = B, ncol = d_theta)
  
  for (b in 1:B) {
    if (verbose && b %% 100 == 0) {
      cat("  - Replication", b, "/", B, "\n")
    }
    
    # Draw multipliers
    G <- stats::rnorm(n)
    
    # Compute perturbed score average
    s_bar_b <- colMeans(G * scores)
    
    # Form bootstrap draw
    bootstrap_draws[b, ] <- theta_hat + Sigma_hat %*% s_bar_b
  }
  
  # Step 4: Construct confidence intervals
  if (verbose) {
    cat("\nStep 4: Computing confidence intervals...\n")
  }
  
  alpha <- 1 - confidence_level
  
  se <- apply(bootstrap_draws, 2, stats::sd)
  ci_lower <- apply(bootstrap_draws, 2, stats::quantile, probs = alpha/2)
  ci_upper <- apply(bootstrap_draws, 2, stats::quantile, probs = 1 - alpha/2)
  
  if (verbose) {
    cat("\n============================================================================\n")
    cat("Bootstrap Inference Completed\n")
    cat("============================================================================\n\n")
    cat("Summary Statistics:\n")
    cat("  - Number of parameters:", d_theta, "\n")
    cat("  - Bootstrap replications:", B, "\n")
    cat("  - Confidence level:", confidence_level * 100, "%\n\n")
  }
  
  # Create results table
  results_table <- data.frame(
    Parameter = param_names,
    Estimate = theta_hat,
    Std.Error = se,
    CI.Lower = ci_lower,
    CI.Upper = ci_upper
  )
  
  if (verbose) {
    cat("Parameter Estimates with", confidence_level * 100, "% Confidence Intervals:\n")
    print(results_table, row.names = FALSE, digits = 4)
    cat("\n")
  }
  
  result <- list(
    theta_hat = theta_hat,
    se = se,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    bootstrap_draws = bootstrap_draws,
    vcov = Sigma_hat,
    param_names = param_names,
    results_table = results_table,
    confidence_level = confidence_level
  )
  
  class(result) <- c("bootstrap_inference", "list")
  return(result)
}

#' Summary Method for Bootstrap Inference
#'
#' Print a summary of bootstrap inference results
#'
#' @param object Output from bootstrap_inference function
#' @param ... Additional arguments (not used)
#' @export
summary.bootstrap_inference <- function(object, ...) {
  cat("\nBootstrap Inference Summary\n")
  cat("===========================\n\n")
  cat("Confidence Level:", object$confidence_level * 100, "%\n")
  cat("Number of Parameters:", length(object$theta_hat), "\n")
  cat("Number of Bootstrap Draws:", nrow(object$bootstrap_draws), "\n\n")
  
  print(object$results_table, row.names = FALSE, digits = 4)
  
  invisible(object)
}

