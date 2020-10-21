library("abtest")
library("rstan")
library("bridgesampling")

##  See https://arxiv.org/abs/1905.02068

###  Is employee retention improved by training
##  
##   consultant claims at 6 months employee retention
##   is improved by 15% on avg
##   95% CI = (.025, .0275)
##   
##
##   1000 observations (500 in each group)
##   at 6 months
##   w/o training = 249  = p-hat = .498
##   w   training = 269  = p-hat = .538
##
##   delta = 4%
##


data("seqdata")

#-------------------------------------------------------------------------------
# Prior
#-------------------------------------------------------------------------------


###  Is employee retention improved by training
##  
##   consultant claims at 6 months employee retention
##   is improved by 15% on avg
##   95% CI = (.025, .0275)
##   



# elicit prior
prior_par <- elicit_prior(
  q = c(0.025, 0.15, 0.275),
  prob = c(.025, .5, .975),
  what = "arisk"
)

# plot prior
cairo_ps(
  paste0("plots/prior_arisk.eps"),
  width = 530 / 72,
  height = 350 / 72
)
plot_prior(prior_par, what = "arisk")
#dev.off()

cairo_ps(
  paste0("plots/prior_logor.eps"),
  width = 530 / 72,
  height = 350 / 72
)
plot_prior(prior_par)
#dev.off()

cairo_ps(
  paste0("plots/prior_p1p2.eps"),
  width = 530 / 72,
  height = 350 / 72
)
plot_prior(prior_par, what = "p1p2")
#dev.off()

#-------------------------------------------------------------------------------
# A/B test
#-------------------------------------------------------------------------------

set.seed(1)
ab <- ab_test(data = seqdata, prior_par = prior_par)
print(ab)

#-------------------------------------------------------------------------------
# Posterior
#-------------------------------------------------------------------------------

# plot posterior probabilities
cairo_ps(
  paste0("plots/post_probs.eps"),
  width = 500 / 72,
  height = 200 / 72
)
prob_wheel(ab)
#dev.off()

# plot sequential analysis
cairo_ps(
  paste0("plots/sequential.eps"),
  width = 530 / 72,
  height = 400 / 72
)
plot_sequential(ab, thin = 4)
#dev.off()

# plot posterior
cairo_ps(
  paste0("plots/posterior_arisk.eps"),
  width = 530 / 72,
  height = 400 / 72
)
plot_posterior(ab, what = "arisk")
#dev.off()

cairo_ps(
  paste0("plots/posterior_logor.eps"),
  width = 530 / 72,
  height = 400 / 72
)
plot_posterior(ab, what = "logor")
dev.off()

cairo_ps(
  paste0("plots/posterior_p1p2.eps"),
  width = 530 / 72,
  height = 400 / 72
)
plot_posterior(ab, what = "p1p2")
#dev.off()

#-------------------------------------------------------------------------------
# Simulation (Laplace vs. Bridge Sampling)
#-------------------------------------------------------------------------------

H1code <-
  "data {
 int<lower = 1> n1;
 int<lower = 1> n2;
 int<lower = 0, upper = n1> y1;
 int<lower = 0, upper = n2> y2;
 real mu_beta;
 real mu_psi;
 real<lower = 0> sigma_beta;
 real<lower = 0> sigma_psi;
}
parameters {
 real psi;
 real beta;
}
transformed parameters {
 real eta1 = beta - psi/2;
 real eta2 = beta + psi/2;
 real<lower = 0, upper = 1> p1 = inv_logit(eta1);
 real<lower = 0, upper = 1> p2 = inv_logit(eta2);
}
model {
 // likelihood
 target += y1 * log(p1) + (n1 - y1) * log(1 - p1);
 target += y2 * log(p2) + (n2 - y2) * log(1 - p2);

 // priors
 target += normal_lpdf(beta | mu_beta, sigma_beta);
 target += normal_lpdf(psi | mu_psi, sigma_psi);
}"
stanmodelH1 <- stan_model(model_code = H1code)

H0code <-
  "data {
 int<lower = 1> n1;
 int<lower = 1> n2;
 int<lower = 0, upper = n1> y1;
 int<lower = 0, upper = n2> y2;
 real mu_beta;
 real<lower = 0> sigma_beta;
}
parameters {
 real beta;
}
transformed parameters {
 real eta1 = beta;
 real eta2 = beta;
 real<lower = 0, upper = 1> p1 = inv_logit(eta1);
 real<lower = 0, upper = 1> p2 = inv_logit(eta2);
}
model {
 // likelihood
 target += y1 * log(p1) + (n1 - y1) * log(1 - p1);
 target += y2 * log(p2) + (n2 - y2) * log(1 - p2);

 // prior
 target += normal_lpdf(beta | mu_beta, sigma_beta);
}"
stanmodelH0 <- stan_model(model_code = H0code)

# function for computing log BF10
logbf10stan <- function(data,
                        prior_par = list(
                          mu_psi = 0,
                          sigma_psi = 1,
                          mu_beta = 0,
                          sigma_beta = 1
                        ),
                        warmup = 500,
                        nchains = 4,
                        nperchain = 10000,
                        adapt_delta = .99,
                        method = "normal",
                        silent = FALSE) {
  
  # H1
  stanfitH1 <- sampling(
    stanmodelH1,
    data = c(
      data,
      list(
        mu_beta = prior_par$mu_beta,
        sigma_beta = prior_par$sigma_beta,
        mu_psi = prior_par$mu_psi,
        sigma_psi = prior_par$sigma_psi
      )
    ),
    warmup = warmup,
    chains = nchains,
    iter = warmup + nperchain,
    control = list(adapt_delta = adapt_delta)
  )
  bridgeH1 <-
    bridge_sampler(stanfitH1,
                   method = method,
                   silent = silent)
  
  # H0
  stanfitH0 <- sampling(
    stanmodelH0,
    data = c(data, list(
      mu_beta = prior_par$mu_beta,
      sigma_beta = prior_par$sigma_beta
    )),
    warmup = warmup,
    chains = nchains,
    iter = warmup + nperchain,
    control = list(adapt_delta = adapt_delta)
  )
  bridgeH0 <-
    bridge_sampler(stanfitH0,
                   method = method,
                   silent = silent)
  
  logbf10 <- bf(bridgeH1, bridgeH0, log = TRUE)
  
  return(logbf10)
  
}

# simulation settings
n1 <- c(5, 10, 20, 50, 100)
n2 <- n1

prop1 <- 1:4 / 5
prop2 <- prop1

prior_par <- list(mu_psi = 0,
                  sigma_psi = 1,
                  mu_beta = 0,
                  sigma_beta = 1)

# run simulation
iter <- 1
r <- list()

set.seed(1)

for (i in seq_along(n1)) {
  for (j in seq_along(n2)) {
    for (k in seq_along(prop1)) {
      for (l in seq_along(prop2)) {
        
        print(iter)
        
        data <- list(y1 = prop1[k] * n1[i],
                     y2 = prop2[l] * n2[j],
                     n1 = n1[i], n2 = n2[j])
        ab <- ab_test(data = data, prior_par = prior_par)
        laplace_logbf10 <- ab$logbf$bf10
        bridge_logbf10 <- logbf10stan(data = data, prior_par)
        r[[iter]] <- list(data = data, laplace_logbf10 = laplace_logbf10,
                          bridge_logbf10 = bridge_logbf10)
        save(r, file = "laplace_bridge.Rdata")
        iter <- iter + 1
        
      }
    }
  }
}

# plotting function
plot_panel <- function(n1, n2, cex = 1.8, cex.axis = 1.8, lwd = 2,
                       cex.lab = 1.7, cex.main = 1.8) {
  
  index <- sapply(r, function(x) x$data$n1 == n1 && x$data$n2 == n2)
  laplace_logbf10 <- sapply(r[index], function(x) x$laplace_logbf10)
  bridge_logbf10 <- sapply(r[index], function(x) x$bridge_logbf10$bf)
  
  xticks <- pretty(laplace_logbf10)
  yticks <- pretty(bridge_logbf10)
  xlim <- range(xticks)
  ylim <- range(yticks)
  
  plot(1, type = "n", axes = FALSE, xlim = xlim, ylim = ylim,
       xlab = "", ylab = "")
  lines(c(max(c(xlim[1], ylim[1])), min(c(xlim[2], ylim[2]))),
        c(max(c(xlim[1], ylim[1])), min(c(xlim[2], ylim[2]))), lwd = lwd)
  points(laplace_logbf10, bridge_logbf10, cex = cex, pch = 21, bg = "grey")
  axis(1, at = xticks, cex.axis = cex.axis)
  axis(2, at = yticks, cex.axis = cex.axis, las = 1)
  mtext("Laplace Log BF10", 1, cex = cex.lab, line = 3.4)
  mtext("Bridge Log BF10", 2, cex = cex.lab, line = 3.8)
  
  mtext(paste0("n1 = ", n1, ", n2 = ", n2), 3, cex = cex.main)
  
}

# multi-panel plot of results
cairo_ps("plots/laplace_bridge.eps",
         height = 5 * 310 / 72,
         width = 5 * 310 / 72)
op <- par(mfrow = c(5, 5), mar = c(5, 6, 4, 4))
for (i in seq_along(n1)) {
  for (j in seq_along(n2)) {
    plot_panel(n1 = n1[i], n2 = n2[j])
  }
}
par(op)
dev.off()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Appendix
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#-------------------------------------------------------------------------------
# Prior
#-------------------------------------------------------------------------------

# default prior for Bayesian A/B test
cairo_ps(
  paste0("plots_appendix/prior_arisk_default.eps"),
  width = 530 / 72,
  height = 350 / 72
)
plot_prior(what = "arisk")
dev.off()

cairo_ps(
  paste0("plots_appendix/prior_logor_default.eps"),
  width = 530 / 72,
  height = 350 / 72
)
plot_prior(what = "logor")
dev.off()

cairo_ps(
  paste0("plots_appendix/prior_p1p2_default.eps"),
  width = 530 / 72,
  height = 350 / 72
)
plot_prior(what = "p1p2")
dev.off()

#-------------------------------------------------------------------------------
# A/B test
#-------------------------------------------------------------------------------

# conduct default A/B test
set.seed(1)
ab_default <- ab_test(data = seqdata)
print(ab_default)

#-------------------------------------------------------------------------------
# Posterior
#-------------------------------------------------------------------------------

### default

# plot posterior probabilities
cairo_ps(
  paste0("plots_appendix/post_probs_default.eps"),
  width = 500 / 72,
  height = 200 / 72
)
prob_wheel(ab_default)
dev.off()

# plot sequential Bayesian A/B test results
cairo_ps(
  paste0("plots_appendix/sequential_default.eps"),
  width = 530 / 72,
  height = 400 / 72
)
plot_sequential(ab_default, thin = 4)
dev.off()

# plot posterior
cairo_ps(
  paste0("plots_appendix/posterior_arisk_default.eps"),
  width = 530 / 72,
  height = 400 / 72
)
plot_posterior(ab_default, what = "arisk")
dev.off()

cairo_ps(
  paste0("plots_appendix/posterior_logor_default.eps"),
  width = 530 / 72,
  height = 400 / 72
)
plot_posterior(ab_default, what = "logor")
dev.off()

cairo_ps(
  paste0("plots_appendix/posterior_p1p2_default.eps"),
  width = 530 / 72,
  height = 400 / 72
)
plot_posterior(ab_default, what = "p1p2")
#dev.off()