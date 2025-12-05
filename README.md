# mlogitsel

Multinomial Logit Model with Sample Selection and Bootstrap Inference

## Installation on GitHub
# devtools::install_github("onilboussim/mlogitsel")
```

## Usage

### Basic Model Estimation

```r
library(mlogitsel)

# Generate example data
set.seed(123)
n <- 500
X <- cbind(intercept = 1, x1 = rnorm(n), x2 = rnorm(n))
Z <- cbind(z1 = rnorm(n), z2 = rnorm(n))
W <- cbind(X, Z)

# Generate selection
delta_true <- c(0.5, 0.02, 0.1, -0.3, 0.05)
prob_select <- plogis(W %*% delta_true)
S <- rbinom(n, 1, prob_select)

# Generate outcomes
eta1 <- X %*% c(-1, 0.05, 0.2)
eta2 <- X %*% c(0.5, -0.03, 0.15)
probs <- exp(cbind(eta1, eta2, 0))
probs <- probs / rowSums(probs)
Y <- apply(probs, 1, function(p) sample(0:2, 1, prob = p))
Y[S == 0] <- NA

# Fit model
result <- mlogit_sel(Y = Y, S = S, X = X, Z = Z)

# View results
print(result)
```

### Bootstrap Inference

```r
# Compute bootstrap standard errors and confidence intervals
boot_result <- bootstrap_inference(result, B = 1000, confidence_level = 0.95)

# View summary
print(boot_result$results_table)

# Access components
boot_result$se            # Standard errors
boot_result$ci_lower      # Lower CI bounds
boot_result$ci_upper      # Upper CI bounds
boot_result$vcov          # Variance-covariance matrix
```

## Features

- Two-step estimation for multinomial logit with sample selection
- Copula-based dependence structure
- Multiplier bootstrap for inference
- Standard errors and confidence intervals
- Comprehensive documentation


