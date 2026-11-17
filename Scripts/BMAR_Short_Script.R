
#---------------------------------------------------#
#      Bayesian meta-analysis and meta-regression   #
#---------------------------------------------------#

  

if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}


pacman::p_load(
  bayesmeta, bayestestR, bayesplot, binom, brms,
  DescTools, dplyr,
  extraDistr,
  flextable, forcats,
  ggdist, ggplot2, ggridges, glue, gt, gtsummary,
  here, 
  kableExtra, 
  loo, 
  meta, metafor,
  patchwork, posterior, 
  readxl, rstan, rstantools,
  scales, shinystan, stringr, 
  tictoc, tidybayes, tidyverse, tiff, 
  webshot2
)

# --------------------------------------------------
#                      Data
# --------------------------------------------------

ds <- read_xlsx(here("Data", "ncc_data.xlsx"))


## **Data exploration**

colnames(ds)

head(ds)

summary(ds)

glimpse(ds)


ds$n <- as.integer(ds$n) 
ds$y <- as.integer(ds$y) 

glimpse(ds)


df1 <- ds %>% 
  select(Study,
         Study_period,
         n, 
         y)


df1 <- df1 %>% 
  mutate(
    pi = binom.confint(y,n, method = "wilson")$mean,
    lb = binom.confint(y,n, method = "wilson")$lower,
    ub = binom.confint(y,n, method = "wilson")$upper,
    vi = (pi*(1-pi))/n
  ) 

## Another way of calculating proportion and variance 
# df1 <- escalc(
#   measure = "PR",
#   data = df1,
#   ni = df1$n,
#   xi = df1$y,
#   var.names = c("Raw_p","Raw_var")
# )

ggplot(df1,
       aes(x = pi,
           y = Study)) +
  geom_point(colour = "orange", size = 2) +
  geom_errorbar(
    aes(xmin = lb, xmax = ub)
  )+
  theme_minimal(base_size = 8) +
  labs(
    title = "Forest plot without pooled estimate",
    x = "Prevalence (95% CI)"
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank()
  )


plot(density(df1$pi),
     main = "Distribution of the reported proportions",
     xlab = "Study-specific prevalence",
     xlim = c(0, 1)
)

# plot(density(df1$Raw_p),
#      main = "Distribution of the reported proportions",
#      xlab = "Study-specific prevalence",
#      xlim = c(0, 1)
# )


# --------------------------------------------------
#       Frequentist meta-analysis
# --------------------------------------------------

mod_raw <- metaprop(
  event = y,
  n = n,
  studlab = Study,
  data = df1,
  sm = "PRAW",
  method = "Inverse",
  method.tau = "DL",
  method.ci = "NAsm"
)

summary(mod_raw)


forest(
  mod_raw,
  common = FALSE,
  print.tau2 = TRUE,
  print.Q = TRUE,
  print.pval.Q = TRUE,
  print.I2 = TRUE,
  leftcols = c("studlab", "event", "n"),
  leftlabs = c("Study", "Cases", "Total"),
  xlab = "Prevalence",
  digits = 2,
  fontsize = 15,
  spacing = 0.9,
  col.square = "navy",
  col.square.lines = "navy",
  col.diamond = "maroon",
  col.diamond.lines = "maroon"
)


mod_logit <- metaprop(
  event = y,
  n = n,
  studlab = Study,
  data = df1,
  sm = "PLO",
  method = "Inverse",
  method.tau = "DL",
  method.ci = "NAsm"
)

summary(mod_logit)



forest(
  mod_logit,
  common = FALSE,
  print.tau2 = TRUE,
  print.Q = TRUE,
  print.pval.Q = TRUE,
  print.I2 = TRUE,
  leftcols = c("studlab", "event", "n"),
  leftlabs = c("Study", "Cases", "Total"),
  xlab = "Prevalence",
  digits = 2,
  fontsize = 15,
  spacing = 0.9,
  col.square = "navy",
  col.square.lines = "navy",
  col.diamond = "maroon",
  col.diamond.lines = "maroon"
)


mod_gllm <- metaprop(
  event = y,
  n = n,
  studlab = Study,
  data = df1,
  sm = "PLOGIT",
  method = "GLMM",
  common = FALSE,
  random = TRUE
)

summary(mod_gllm)



forest(
  mod_gllm,
  common = FALSE,
  print.tau2 = TRUE,
  print.Q = TRUE,
  print.pval.Q = TRUE,
  print.I2 = TRUE,
  leftcols = c("studlab", "event", "n"),
  leftlabs = c("Study", "Cases", "Total"),
  xlab = "Prevalence",
  digits = 2,
  fontsize = 15,
  spacing = 0.9,
  col.square = "navy",
  col.square.lines = "navy",
  col.diamond = "maroon",
  col.diamond.lines = "maroon"
)


# --------------------------------------------------
#                    Priors
# --------------------------------------------------

set.seed(123)
n_samples <- 10000

y <- rt(n_samples, df = 3) * 2.5
x <- plogis(y)
hist(x, breaks = 50, main = "Histogram of inv_logit(y) where y ~ Student_t(3, 0, 2.5)",
     xlab = "inv_logit(y)", col = "lightblue", border = "white")



y <- runif(n_samples, min = 0, max = 1)
x <- plogis(y)
hist(x, breaks = 50, main = "Histogram of inv_logit(y) where y ~ Uniform(0, 1)",
     xlab = "inv_logit(y)", col = "lightblue", border = "white")



y <- rnorm(n_samples, mean = 0, sd = 1)
x <- plogis(y)
hist(x, breaks = 50, main = "Histogram of inv_logit(y) where y ~ Normal(0, 1)",
     xlab = "inv_logit(y)", col = "lightblue", border = "white")


y <- rnorm(n_samples, mean = 0, sd = 1.5)
x <- plogis(y)
hist(x, breaks = 50, main = "Histogram of inv_logit(y) where y ~ Normal(0, 1.5)",
     xlab = "inv_logit(y)", col = "lightblue", border = "white")



y <- rbeta(n_samples, shape1 = 1, shape2 = 1)
x <- plogis(y)
hist(x, breaks = 50, main = "Histogram of inv_logit(y) where y ~ Beta(1, 1)",
     xlab = "inv_logit(y)", col = "lightblue", border = "white")


y <- rlogis(n_samples, location = 0, scale = 1)
x <- plogis(y)
hist(x, breaks = 50, main = "Histogram of inv_logit(y) where y ~ logistic(0, 1)",
     xlab = "inv_logit(y)", col = "lightblue", border = "white")



y <- rcauchy(n_samples, location = 0, scale = 0.5)
y_pos <- y[y > 0]
y_trim <- y_pos[y_pos < quantile(y_pos, 0.99)]
hist(y_trim, breaks = 100, freq = FALSE,
     main = "Half-Cauchy(0, 0.5) Prior for sd",
     xlab = "sd", col = "skyblue", border = "white")
lines(density(y_trim), col = "darkblue", lwd = 2)


set.seed(123)
N_samples <- 1000
N_studies <- nrow(df1)

# Overall intercept: Logistic(0, 1)
mu_samples <- rlogis(
  N_samples,
  location = 0,
  scale = 1
)

# mu_samples <- rnorm(
#   N_samples,
#   -1.386294,
#   1
# )

# Between-study heterogeneity: Half-Cauchy(0, 0.5)
tau_samples <- abs(
  rcauchy(
    N_samples,
    location = 0,
    scale = 0.5
  )
)

# --------------------------------------------------
#     Prior predictive check
# --------------------------------------------------

prior_pred <- map_dfr(seq_len(N_samples), function(i) {
  
  u <- rnorm(
    N_studies,
    mean = 0,
    sd = tau_samples[i]
  )
  
  
  p <- plogis(mu_samples[i] + u)
  
  
  y_pred <- rbinom(
    N_studies,
    size = df1$n,
    prob = p
  )
  
  tibble(
    Study = df1$Study,
    iter = i,
    mu = mu_samples[i],
    tau = tau_samples[i],
    u = u,
    p_pred = p,
    y_pred = y_pred,
    prevalence_pred = y_pred / df1$n
  )
})



ggplot() +
  geom_density(
    data = prior_pred,
    aes(x = prevalence_pred),
    fill = "steelblue",
    alpha = 0.4
  ) +
  geom_density(
    data = df1,
    aes(x = y / n),
    color = "red",
    linewidth = 1
  ) +
  labs(
    title = "Observed vs prior predictive prevalence",
    x = "Prevalence",
    y = "Density"
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    labels = label_percent()
  ) +
  theme_minimal()


# --------------------------------------------------
#          Prior only model
# --------------------------------------------------

prior_set <- c(prior(logistic(0,1), class = "b", coef = "Intercept"),
               prior(cauchy(0,0.5), class="sd"))


prior_model <- brm(
  y | trials(n) ~ 0 + Intercept + (1 | Study),
  data = df1,
  family = binomial(link = "logit"),
  prior = prior_set,
  sample_prior = "only",
  chains = 4,
  iter = 2000,
  warmup = 1000,
  cores = 4,
  seed = 123
)

saveRDS(prior_model, here("Outputs", file = "prior_model.rds"))

## --------results

prior_model <- readRDS (here::here("Outputs","prior_model.rds"))

# Simulate replicated case counts
yrep <- posterior_predict(prior_model)


prior_pred <- as.data.frame(as.table(yrep)) %>%
  as_tibble() %>%
  mutate(
    Study = df1$Study[as.integer(Var2)],
    n = df1$n[as.integer(Var2)],
    prevalence_pred = Freq / n,
    iter = as.integer(Var1)
  )



ggplot() +
  geom_density(
    data = prior_pred,
    aes(x = prevalence_pred),
    fill = "steelblue",
    alpha = 0.4
  ) +
  geom_density(
    data = df1,
    aes(x = y / n),
    color = "red",
    linewidth = 1
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    labels = label_percent()
  ) +
  labs(
    title = "Observed vs prior predictive prevalence",
    x = "Prevalence",
    y = "Density"
  ) +
  theme_minimal()


## We can look at the estimates for both parameters separately

prior_draws <- posterior::as_draws_df(prior_model) %>%
  as.data.frame()

prior_parameters <- prior_draws %>%
  transmute(
    prevalence = plogis(b_Intercept),
    tau = sd_Study__Intercept
  )

p1 <- ggplot(prior_parameters, aes(x = prevalence)) +
  geom_density(fill = "steelblue", alpha = 0.5) +
  scale_x_continuous(
    limits = c(0, 1),
    labels = scales::label_percent()
  ) +
  labs(
    title = "Prior distribution of Ps",
    x = "Overall prevalence",
    y = "Density"
  ) +
  theme_minimal()

# Plot 2: Between-study heterogeneity prior
p2 <- ggplot(prior_parameters, aes(x = tau)) +
  geom_density(fill = "orange", alpha = 0.5) +
  coord_cartesian(xlim = c(0, 3)) +
  labs(
    title = "Prior distribution of Tau",
    x = expression("Between-study SD (" * tau * ")"),
    y = "Density"
  ) +
  theme_minimal()


combined_plot <- p1 + p2 +
  patchwork::plot_layout(ncol = 2)


print(combined_plot)

# --------------------------------------------------
#  Bayesian meta-analysis model
# --------------------------------------------------

prior_set <- c(prior(logistic(0,1), class = Intercept),
               prior(cauchy(0,0.5), class=sd))


b_ma_mod <-brm(
  y|trials(n)~(1|Study),
  data=df1,
  prior=prior_set,
  family=binomial(),
  chains = 4, 
  iter=10000,
  warmup = 5000,
  control = list(adapt_delta = 0.80),
  seed = 123)


saveRDS(b_ma_mod, here("Outputs", file = "b_ma_mod.rds"))


# --------------------------------------------------
#     Model diagnostics
# --------------------------------------------------

b_ma_mod <- readRDS (here::here("Outputs","b_ma_mod.rds"))

prior_summary(b_ma_mod)
plot(b_ma_mod,warmup=1000)
pp_check(b_ma_mod) 
# launch_shinystan(b_ma_mod)

# --------------------------------------------------
#     Results of BMA
# --------------------------------------------------

summary(b_ma_mod)
fixef(b_ma_mod) %>% plogis()



posterior_draw <- as_draws_array(b_ma_mod)
posterior <- posterior_draw %>% 
  rename_variables(
    Proportion = "b_Intercept",
    tau = "sd_Study__Intercept"
  ) %>% 
  mutate_variables(
    Proportion = plogis(Proportion)
  )
summary_results <- as_draws_df(posterior) %>%
  summarise(
    across(
      c(Proportion, tau),
      list(
        Median = ~ median(.x),
        Lower_95 = ~ quantile(.x, 0.025),
        Upper_95 = ~ quantile(.x, 0.975)
      ),
      .names = "{.col}_{.fn}"
    )
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = c("Parameter", ".value"),
    names_pattern = "^(Proportion|tau)_(Median|Lower_95|Upper_95)$"
  )

summary_results

# posterior_summary(
#   b_ma_mod,
#   robust = TRUE, ## for median
#   probs = c(0.025, 0.975)
# )

# --------------------------------------------------
#     Sensitivity analysis
# --------------------------------------------------


prior_set <- c(prior(normal(0,1.5), class = Intercept), ## Only change here
               prior(cauchy(0,0.5), class=sd))

b_ma_mod_1 <-brm(                        ## Give new model name
  y|trials(n)~(1|Study),
  data=df1,
  prior=prior_set,
  family=binomial(),
  chains = 4, 
  iter=10000,
  warmup = 5000,
  control = list(adapt_delta = 0.80),
  seed = 123)
saveRDS(b_ma_mod_1, here("Outputs", file = "ma_se_norm_p.rds"))



b_ma_mod_1 <- readRDS (here::here("Outputs", "ma_se_norm_p.rds"))

summary(b_ma_mod_1)
fixef(b_ma_mod_1) %>% plogis()
pp_check(b_ma_mod_1)
prior_summary(b_ma_mod_1)



log(0.33/(1-0.33)) 

y <- rnorm(n_samples, mean = -0.7081851, sd = 1)
x <- plogis(y)
hist(x, breaks = 50, main = "Histogram of inv_logit(y) where y ~ Normal(0.33, 1)",
     xlab = "inv_logit(y)", col = "lightblue", border = "white")



prior_set <- c(prior(normal(-0.7081851,1), class = Intercept), 
               prior(cauchy(0,0.5), class=sd))

b_ma_mod_2 <-brm(                        
  y|trials(n)~(1|Study),
  data=df1,
  prior=prior_set,
  family=binomial(),
  chains = 4, 
  iter=10000,
  warmup = 5000,
  control = list(adapt_delta = 0.80),
  seed = 123)
saveRDS(b_ma_mod_2, here("Outputs", file = "ma_se_norm_0.33.rds"))


b_ma_mod_2 <- readRDS (here::here("Outputs", "ma_se_norm_0.33.rds"))

summary(b_ma_mod_2)
fixef(b_ma_mod_2) %>% plogis()
pp_check(b_ma_mod_2)
prior_summary(b_ma_mod_2)


y <- rnorm(n_samples, mean = log(0.33/(1-0.33)), sd = 0.1)
x <- plogis(y)
hist(x, breaks = 50, main = "Histogram of inv_logit(y) where y ~ Normal(0.33, 0.1)",
     xlab = "inv_logit(y)", col = "lightblue", border = "white")



prior_set <- c(prior(normal(-0.7081851,0.1), class = Intercept), 
               prior(cauchy(0,0.5), class=sd))

b_ma_mod_3 <-brm(                        
  y|trials(n)~(1|Study),
  data=df1,
  prior=prior_set,
  family=binomial(),
  chains = 4, 
  iter=10000,
  warmup = 5000,
  control = list(adapt_delta = 0.80),
  seed = 123)

saveRDS(b_ma_mod_3, here("Outputs", file = "ma_se_norm_0.33_0.1.rds"))


b_ma_mod_3 <- readRDS (here::here("Outputs", "ma_se_norm_0.33_0.1.rds"))

summary(b_ma_mod_3)
fixef(b_ma_mod_3) %>% plogis()
pp_check(b_ma_mod_3)
prior_summary(b_ma_mod_3)

# --------------------------------------------------
#     Visualization of results
# --------------------------------------------------


if (requireNamespace("hexbin", quietly = TRUE)) { 
  mcmc_hex(posterior, pars = c("Proportion", "tau")) 
}




study.draws <- spread_draws(b_ma_mod, 
                            r_Study[Study,], 
                            b_Intercept) %>%
  mutate(b_Intercept = r_Study + b_Intercept) 


pooled.effect.draws <- spread_draws(b_ma_mod, 
                                    b_Intercept) %>%
  mutate(Study = "Pooled Effect") 


forest.data <- bind_rows(study.draws, 
                         pooled.effect.draws) %>%
  ungroup() %>%
  mutate(b_Intercept = plogis(b_Intercept)) %>% 
  mutate(Study = str_replace_all(Study, "[.]", " ")) %>% 
  mutate(Study = reorder(Study, 
                         b_Intercept, 
                         decreasing = F)) 


forest.data$Study <- factor(forest.data$Study)


forest.data.summary <- group_by(forest.data, 
                                Study) %>%
  median_qi(b_Intercept) 

p <- ggplot(aes(b_Intercept, 
           relevel(Study, 
                   "Pooled Effect", 
                   after = Inf)), 
       data = forest.data) +
  geom_vline(xintercept = plogis(fixef(b_ma_mod)[1, 1]), 
             color = "#E64B35", 
             size = 1) + 
  geom_vline(xintercept = plogis(fixef(b_ma_mod)[1, 3:4]), 
             color = "grey", 
             linetype = 2) + 
  geom_density_ridges(
    aes(fill = Study),
    rel_min_height = 0.01,
    col = NA,
    scale = 1,
    alpha = 1
  ) +
  scale_fill_manual(
    values = c(
      "Pooled Effect" = "#E64B35",
      setNames(
        rep("#4292C6", length(setdiff(
          levels(forest.data$Study), "Pooled Effect"
        ))),
        setdiff(levels(forest.data$Study), "Pooled Effect")
      )
    ),
    guide = "none"
  )  + 
  stat_pointinterval(.width = 0.95, 
                     linewidth = 0.5) + 
  geom_text(data = mutate_if(forest.data.summary, 
                             is.numeric, 
                             round, 2), 
            aes(label = glue("{b_Intercept} [{.lower}, {.upper}]"),x = Inf), 
            hjust = "inward") + 
  labs(x = "Estimated Proportion", 
       y = NULL, 
       title = "Estimated Proportion \n(95%BCI)", 
       subtitle = "Study") + 
  theme_minimal() + 
  theme(
    text = element_text(family = "Times"), 
    axis.title.x = element_text(size = 14, family = "Times", colour = "black"), 
    axis.text.x = element_text(size = 10, family = "Times", colour = "black"), 
    axis.title.y = element_text(size = 10, family = "Times", colour = "black"), 
    axis.text.y = element_text(size = 10, family = "Times", colour = "black"), 
    plot.title = element_text(hjust = 1, size = 12, family = "Times", face = "bold"), 
    plot.subtitle = element_text(hjust = -0.1, size = 12, family = "Times", face = "bold"), #
    plot.margin = margin(20, 20, 40, 20), 
    legend.text = element_text(size = 10, family = "Times"), 
    legend.title = element_text(size = 12, family = "Times"), 
    axis.line = element_line(size = 0.5), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank() 
  )
p


# ggsave(
#   here("Outputs", "stan_forest_plot.pdf"),
#   plot = p,
#   width = 10,
#   height = 8,
#   unit = "in"
# )

# --------------------------------------------------
#     Model in Stan
# --------------------------------------------------
stan_model <- "
  data{
    int<lower = 1> N; 
    array[N] int<lower = 0> n;
    array[N] int<lower = 0> y;
  }
  
 parameters{
   real mu;
   real <lower = 0> sigma;
   vector[N] theta;
 } 
 
 model{
   mu ~ logistic(0,1);
   sigma ~ cauchy(0,0.5);
   
   theta ~ normal(mu,sigma);
   
   y ~ binomial_logit(n,theta);
   
 }
 
 generated quantities {
   real pooled_prevalence;
   vector[N] prevalence;
   
   pooled_prevalence = inv_logit(mu);
   prevalence = inv_logit(theta);
 }
"

stan_data <- list(
  N = nrow(df1),
  n = as.integer(df1$n),
  y = as.integer(df1$y)
)


stan_ma_mod <- stan(
  model_code = stan_model,
  data = stan_data,
  chains = 4,
  iter = 10000,
  warmup = 5000,
  seed = 123
)


saveRDS(stan_ma_mod, here("Outputs", file = "stan_ma_mod.rds"))



stan_ma_mod <- readRDS (here::here("Outputs","stan_ma_mod.rds"))

summary(
  stan_ma_mod ,
  pars = c("mu", "sigma", "pooled_prevalence")
)$summary


#### **Trace plot**
mcmc_trace(
  as.array(stan_ma_mod),
  pars = c("mu", "sigma")
)


#### **Posterior distribution**
mcmc_hist(
  as.array(stan_ma_mod),
  pars = c("mu", "sigma")
)


#### **Check divergent transitions**
sampler_params <- get_sampler_params(
  stan_ma_mod,
  inc_warmup = FALSE
)

sum(sapply(
  sampler_params,
  function(x) sum(x[, "divergent__"])
))

### **Results**
stan_post_draws <- rstan::extract(stan_ma_mod)

# stan_ma_results <- data.frame(
#   Study = df1$Study,
#   Median = apply(stan_post_draws$prevalence, 2, median),
#   Lower = apply(stan_post_draws$prevalence, 2, quantile, 0.025),
#   Upper = apply(stan_post_draws$prevalence, 2, quantile, 0.975)
# )
# 
# stan_ma_results

quantile(
  stan_post_draws$pooled_prevalence,
  probs = c(0.025, 0.50, 0.975)
)


# --------------------------------------------------
#     Sub-group analysis
# --------------------------------------------------

b_ma_mod_periodA <-brm(
  y|trials(n)~(1|Study),
  data=df1 %>% filter(Study_period == "A"),
  prior=c(prior(logistic(0,1), class = Intercept),
          prior(cauchy(0,0.5), class=sd)),
  family=binomial(),
  chains = 4, 
  iter=10000,
  warmup = 5000,
  control = list(adapt_delta = 0.80),
  seed = 123)

saveRDS(b_ma_mod_periodA, here("Outputs", file = "b_ma_mod_periodA.rds"))


b_ma_mod_periodB <-brm(
  y|trials(n)~(1|Study),
  data=df1 %>% filter(Study_period == "B"),
  prior=c(prior(logistic(0,1), class = Intercept),
          prior(cauchy(0,0.5), class=sd)),
  family=binomial(),
  chains = 4, 
  iter=10000,
  warmup = 5000,
  control = list(adapt_delta = 0.80),
  seed = 123)


saveRDS(b_ma_mod_periodB, here("Outputs", file = "b_ma_mod_periodB.rds"))



b_ma_mod_periodA <- readRDS (here::here("Outputs","b_ma_mod_periodA.rds"))
b_ma_mod_periodB <- readRDS (here::here("Outputs","b_ma_mod_periodB.rds"))

draws_A <- as_draws_df(b_ma_mod_periodA)
draws_B <- as_draws_df(b_ma_mod_periodB)

prev_A <- plogis(draws_A$b_Intercept)
prev_B <- plogis(draws_B$b_Intercept)


results <- tibble(
  Period = c("Period A", "Period B"),
  Median = c(median(prev_A), median(prev_B)),
  Lower = c(
    quantile(prev_A, 0.025),
    quantile(prev_B, 0.025)
  ),
  Upper = c(
    quantile(prev_A, 0.975),
    quantile(prev_B, 0.975)
  )
) %>%
  mutate(
    across(c(Median, Lower, Upper), ~ .x * 100)
  )

print(results)



# Difference in prevalence (Period B - Period A)
difference <- prev_B - prev_A

difference_results <- tibble(
  Median_difference = median(difference) * 100,
  Lower_95_CrI = quantile(difference, 0.025) * 100,
  Upper_95_CrI = quantile(difference, 0.975) * 100,
  Probability_increase = mean(difference > 0),
  Probability_decrease = mean(difference < 0)
)

print(difference_results)


plot_data <- results %>%
  mutate(Period = factor(Period, levels = c("Period A", "Period B")))

p1 <- ggplot(
  plot_data,
  aes(x = Period, y = Median)
) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = Lower, ymax = Upper),
    width = 0.15,
    linewidth = 0.8
  ) +
  labs(
    title = "Pooled NCC prevalence by study period",
    x = NULL,
    y = "Prevalence (%)"
  ) +
  theme_minimal()


p2 <- ggplot(
  tibble(difference = difference * 100),
  aes(x = difference)
) +
  geom_density(
    fill = "steelblue",
    alpha = 0.5
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  labs(
    title = "Posterior distribution of the difference",
    x = "Difference in prevalence (B - A), percentage points",
    y = "Density"
  ) +
  theme_minimal()

# Combine plots
p1 + p2 + plot_layout(ncol = 2)

## What is the probability that prevalence was reduced in period B?
prob_B_lower <- mean(prev_B < prev_A)
prob_B_lower

# --------------------------------------------------
#     Meta-regression
# --------------------------------------------------
### **Data set for meta-regression**


df2 <- ds %>% 
  select(Study,
         Study_period,
         Region,
         n, 
         y)
unique(df2$Study_period)
unique(df2$Region)


table(df2$Study_period)
table(df2$Region)



df2 <- df2 %>% 
  mutate(
    Period_cat = case_when(
      Study_period == "A" ~ 0,
      Study_period == "B" ~ 1
    ),
    Region_cat = case_when(
      Region == "AFR" ~ 0,
      Region == "AMR" ~ 1,
      Region == "SEAR" ~ 2
    )
  ) %>% 
  mutate(
    Period_cat = as.factor(Period_cat),
    Region_cat = as.factor(Region_cat)
  ) %>% 
  mutate(
    Period_cat = relevel(Period_cat, ref = "0"),
    Region_cat = relevel(Region_cat, ref = "0")
  )


# --------------------------------------------------
#     Frequentist meta-regression
# --------------------------------------------------


frq_mr_p <- rma.glmm(
  data = df2,
  mods = ~ Period_cat,
  xi = y,
  ni = n,
  measure = "PLO",
  method = "ML"
  
)
summary(frq_mr_p)


results <- data.frame(
  Parameter = c(
    "Intercept (Period A)",
    "Period B vs Period A"
  ),
  Estimate = as.numeric(coef(frq_mr_p)),
  SE = frq_mr_p$se,
  OR = exp(as.numeric(coef(frq_mr_p))),
  CI_lower = exp(frq_mr_p$ci.lb),
  CI_upper = exp(frq_mr_p$ci.ub),
  p_value = frq_mr_p$pval
)

results


# --------------------------------------------------
#     Prior predictive check
# --------------------------------------------------

prior_set <- c(prior(logistic(0,1), class = Intercept),
               prior(normal(0,1), class = b),
               prior(cauchy(0,0.5), class=sd))


mr_prior_model <- brm(
  y | trials(n) ~ Study_period + (1 | Study),
  data = df2,
  family = binomial(),
  prior = prior_set,
  sample_prior = "only",
  chains = 4,
  iter = 4000,
  warmup = 1000
)

saveRDS(mr_prior_model, here("Outputs", file = "mr_prior_model.rds"))



mr_prior_model <- readRDS(
  here::here("Outputs", "mr_prior_model.rds")
)


yrep <- posterior_predict(mr_prior_model)

prior_pred <- as.data.frame(as.table(yrep)) %>%
  as_tibble() %>%
  mutate(
    Study = df2$Study[as.integer(Var2)],
    Study_period = df2$Study_period[as.integer(Var2)],
    n = df2$n[as.integer(Var2)],
    prevalence_pred = Freq / n,
    iter = as.integer(Var1)
  )

ggplot() +
  geom_density(
    data = prior_pred,
    aes(x = prevalence_pred),
    fill = "steelblue",
    alpha = 0.4
  ) +
  geom_density(
    data = df2,
    aes(x = y / n),
    color = "red",
    linewidth = 1
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    labels = scales::label_percent()
  ) +
  labs(
    title = "Prior predictive check: Meta-regression",
    subtitle = "Observed vs prior-predicted prevalence",
    x = "Prevalence",
    y = "Density"
  ) +
  theme_minimal()



ggplot() +
  geom_density(
    data = prior_pred,
    aes(x = prevalence_pred),
    fill = "steelblue",
    alpha = 0.4
  ) +
  geom_density(
    data = df2,
    aes(x = y / n),
    color = "red",
    linewidth = 1
  ) +
  facet_wrap(~Study_period) +
  scale_x_continuous(
    limits = c(0, 1),
    labels = scales::label_percent()
  ) +
  labs(
    title = "Prior predictive check by study period",
    x = "Prevalence",
    y = "Density"
  ) +
  theme_minimal()



prior_draws <- as_draws_df(mr_prior_model)
prior_beta <- prior_draws %>%
  mutate(
    p_A = plogis(b_Intercept),
    p_B = plogis(b_Intercept + b_Study_periodB),
    beta = b_Study_periodB,
    OR = exp(beta),
    prevalence_difference = p_B - p_A
  )

ggplot(prior_beta, aes(x = OR)) +
  geom_density(fill = "steelblue", alpha = 0.4) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_x_log10() +
  labs(
    title = "Prior implications for the odds ratio",
    x = "Odds ratio (log scale)",
    y = "Density"
  ) +
  theme_minimal()

ggplot(prior_beta, aes(x = prevalence_difference)) +
  geom_density(fill = "steelblue", alpha = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_x_continuous(
    labels = scales::label_percent()
  ) +
  labs(
    title = "Prior predictive implications of study period",
    x = "Difference in prevalence (B - A)",
    y = "Density"
  ) +
  theme_minimal()

prior_beta %>%
  summarise(
    across(
      c(beta, OR, prevalence_difference),
      list(
        median = median,
        lower = ~quantile(.x, 0.025),
        upper = ~quantile(.x, 0.975)
      )
    )
  )


# --------------------------------------------------
#     Bayesian meta-regression: period
# --------------------------------------------------

mr_mod_p_raw <- brm(
  y | trials(n) ~ Period_cat + (1 | Study),
  data = df2,
  family = binomial(link = "logit"),
  prior = c(prior(logistic(0,1), class = Intercept),
            prior(normal(0,1), class = b),
            prior(cauchy(0,0.5), class=sd)),
  chains = 4,
  cores = 4,
  iter = 10000,
  seed = 123)

saveRDS(mr_mod_p_raw, here("Outputs", file = "mr_mod_p_raw.rds"))


mr_mod_p_raw <- readRDS(
  here::here("Outputs", "mr_mod_p_raw.rds")
)

summary(mr_mod_p_raw)
fixef(mr_mod_p_raw) %>% exp()
loo(mr_mod_p_raw)
pp_check(mr_mod_p_raw)
plot(mr_mod_p_raw)
prior_summary(mr_mod_p_raw)


mr_mod_p <- brm(
  y | trials(n) ~ Period_cat + (1 | Study),
  data = df2,
  family = binomial(link = "logit"),
  prior = c(prior(logistic(0,1), class = Intercept),
            prior(normal(0,1), class = b),
            prior(cauchy(0,0.5), class=sd)),
  chains = 4,
  cores = 4,
  iter = 10000,
  save_pars = save_pars(all = TRUE), 
  seed = 123)

mr_mod_p <- add_criterion(mr_mod_p, 
                          "loo", 
                          moment_match = TRUE)

saveRDS(mr_mod_p, here("Outputs", file = "mr_mod_p.rds"))



mr_mod_p <- readRDS(
  here::here("Outputs", "mr_mod_p.rds")
)

summary(mr_mod_p)
fixef(mr_mod_p) %>% exp()
loo(mr_mod_p)
pp_check(mr_mod_p)
plot(mr_mod_p)
prior_summary(mr_mod_p)


## Model comparison
loo_p_raw <- loo(mr_mod_p_raw)
loo_p <- loo(mr_mod_p)
loo_compare(loo_p_raw,loo_p)


loo_p_raw <- loo(mr_mod_p_raw, moment_match = TRUE)
loo_p <- loo(mr_mod_p, moment_match = TRUE)
loo_compare(loo_p_raw,loo_p)

# --------------------------------------------------
#     Bayesian meta-regression: region
# --------------------------------------------------

mr_mod_r <- brm(
  y | trials(n) ~ Region_cat + (1 | Study),
  data = df2,
  family = binomial(link = "logit"),
  prior = c(prior(logistic(0,1), class = Intercept),
            prior(normal(0,1), class = b),
            prior(cauchy(0,0.5), class=sd)),
  chains = 4,
  cores = 4,
  iter = 10000,
  save_pars = save_pars(all = TRUE), 
  seed = 123)

mr_mod_r <- add_criterion(mr_mod_r, 
                           "loo", 
                           moment_match = TRUE)

saveRDS(mr_mod_r, here("Outputs", file = "mr_mod_r.rds"))


mr_mod_r <- readRDS(
  here::here("Outputs", "mr_mod_r.rds")
)

summary(mr_mod_r)
fixef(mr_mod_r) %>% exp()
loo(mr_mod_r)
pp_check(mr_mod_r)
plot(mr_mod_r)
prior_summary(mr_mod_r)


loo_p <- loo(mr_mod_p)
loo_r <- loo(mr_mod_r)
loo_compare(loo_p,loo_r)

# --------------------------------------------------
#     Bayesian meta-regression: period and region
# --------------------------------------------------

mr_mod_pr <- brm(
  y | trials(n) ~ Period_cat + Region_cat + (1 | Study),
  data = df2,
  family = binomial(link = "logit"),
  prior = c(prior(logistic(0,1), class = Intercept),
            prior(normal(0,1), class = b),
            prior(cauchy(0,0.5), class=sd)),
  chains = 4,
  cores = 4,
  iter = 10000,
  save_pars = save_pars(all = TRUE), 
  seed = 123)

mr_mod_pr <- add_criterion(mr_mod_pr, 
                          "loo", 
                          moment_match = TRUE)

saveRDS(mr_mod_pr, here("Outputs", file = "mr_mod_pr.rds"))


mr_mod_pr <- readRDS(
  here::here("Outputs", "mr_mod_pr.rds")
)

summary(mr_mod_pr)
fixef(mr_mod_pr) %>% exp()
loo(mr_mod_pr)
pp_check(mr_mod_pr)
plot(mr_mod_pr)
prior_summary(mr_mod_pr)
# pareto_k_values(loo(mr_mod_pr))
# plot(loo(mr_mod_pr))

loo_p <- loo(mr_mod_p)
loo_pr <- loo(mr_mod_pr)
loo_compare(loo_p,loo_pr)


mr_mod_pr_i <- brm(
  y | trials(n) ~ Period_cat * Region_cat + (1 | Study),
  data = df2,
  family = binomial(link = "logit"),
  prior = c(prior(logistic(0,1), class = Intercept),
            prior(normal(0,1), class = b),
            prior(cauchy(0,0.5), class=sd)),
  chains = 4,
  cores = 4,
  iter = 10000,
  save_pars = save_pars(all = TRUE), 
  seed = 123)

mr_mod_pr_i <- add_criterion(mr_mod_pr_i, 
                           "loo", 
                           moment_match = TRUE)

saveRDS(mr_mod_pr_i, here("Outputs", file = "mr_mod_pr_i.rds"))



mr_mod_pr_i <- readRDS(
  here::here("Outputs", "mr_mod_pr_i.rds")
)

prior_summary(mr_mod_pr_i)
summary(mr_mod_pr_i)
fixef(mr_mod_pr_i) %>% exp()
loo(mr_mod_pr_i)
pp_check(mr_mod_pr_i)
plot(mr_mod_pr_i)
# pareto_k_values(loo(mr_mod_pr_i))
# plot(loo(mr_mod_pr_i))

# --------------------------------------------------
#     Bayesian meta-regression: Random intercept and slope
# --------------------------------------------------

gen_ds <- read_xlsx(here("Data", "ncc_data_by_gender.xlsx"))

df3 <- gen_ds %>% 
  mutate(Gender_cat = case_when(
    Gender == "Male" ~ 0,
    Gender == "Female" ~ 1
  ),
    Gender_cat = factor(Gender_cat),
    Gender_cat = relevel(Gender_cat, ref = "0")
  )


rr_mr_g <- brm(
  y | trials(n) ~ Gender_cat + (1 + Gender_cat | Study),
  data = df3,
  family = binomial(link = "logit"),
  prior = c(
    prior(logistic(0, 1), class = Intercept),
    prior(normal(0, 1), class = b),
    prior(exponential(2), class = sd, group = Study, coef = Intercept),
    prior(exponential(2), class = sd, group = Study, coef = Gender_cat1),
    prior(lkj(2), class = cor)
  ),
  chains = 4,
  cores = 4,
  iter = 10000,
  save_pars = save_pars(all = TRUE),
  seed = 123
)



saveRDS(rr_mr_g, here("Outputs", file = "rr_mr_g.rds"))


rr_mr_g <- readRDS(
  here::here("Outputs", "rr_mr_g.rds")
)

prior_summary(rr_mr_g)
summary(rr_mr_g)
fixef(rr_mr_g) %>% exp()
loo(rr_mr_g)
pp_check(rr_mr_g)
plot(rr_mr_g)

# --------------------------------------------------
#     Bayesian meta-analysis of measures of association
# --------------------------------------------------

## Data

data(dat.bcg, package = "metafor")

dat <- escalc(
  measure = "OR",
  ai = tpos,
  bi = tneg,
  ci = cpos,
  di = cneg,
  data = dat.bcg
)


dat %>%
  select(trial, author, year, yi, vi) %>%
  print()

plot(density(dat$yi),
     main = "Distribution of ORs")


## **Model**


priors <- c(
  prior(normal(0, 1.5), class = "Intercept"),
  prior(cauchy(0, 0.5), class = "sd")
)


bcg_bayes <- brm(
  yi | se(sqrt(vi)) ~ 1 + (1 | trial),
  data = dat,
  family = gaussian(),
  prior = priors,
  chains = 4,
  cores = 4,
  iter = 10000,
  warmup = 5000,
  seed = 1234,
  control = list(
    adapt_delta = 0.95,
    max_treedepth = 12
  )
)

saveRDS(bcg_bayes, here("Outputs", file = "bcg_bayes.rds"))


## **Model diagnostics**


bcg_bayes <- readRDS(
  here::here("Outputs", "bcg_bayes.rds")
)

plot(bcg_bayes)
rhat(bcg_bayes)
neff_ratio(bcg_bayes)
pp_check(bcg_bayes)

## **Summary**


summary(bcg_bayes)
fixef(bcg_bayes)
VarCorr(bcg_bayes)


draws <- as_draws_df(bcg_bayes)
pooled_OR <- exp(draws$b_Intercept)

quantile(
  pooled_OR,
  probs = c(0.025, 0.50, 0.975)
)

# Posterior probabilities
## Probability that BCG vaccination reduces TB odds
mean(pooled_OR < 1)

## Probability that BCG vaccination increases TB odds
mean(pooled_OR > 1)


## **Visualization**

#### **Posterior distribution**

plot(
  density(pooled_OR),
  main = "Posterior distribution of the pooled OR",
  xlab = "Odds ratio",
  lwd = 2
)

abline(v = 1, lty = 2, col = "red")


#### **Forest plots**

study_effects <- dat %>%
  mutate(
    observed_OR = exp(yi),
    observed_lower = exp(yi - 1.96 * sqrt(vi)),
    observed_upper = exp(yi + 1.96 * sqrt(vi)),
    study = paste(author, year)
  )


pooled_OR <- quantile(
  exp(draws$b_Intercept),
  probs = c(0.025, 0.50, 0.975)
)


posterior_study <- posterior_linpred(
  bcg_bayes,
  re_formula = NULL,
  transform = FALSE
)

bayes_effects <- data.frame(
  trial = dat$trial,
  median = apply(posterior_study, 2, median),
  lower = apply(posterior_study, 2, quantile, probs = 0.025),
  upper = apply(posterior_study, 2, quantile, probs = 0.975)
) %>%
  mutate(
    bayes_OR = exp(median),
    bayes_lower = exp(lower),
    bayes_upper = exp(upper)
  )


forest_data <- study_effects %>%
  left_join(bayes_effects, by = "trial") %>%
  mutate(
    study = factor(study, levels = study[order(observed_OR)])
  )


ggplot(forest_data, aes(x = observed_OR, y = study)) +
  geom_errorbar(
    aes(xmin = observed_lower, xmax = observed_upper),
    orientation = "y",
    width = 0.2
  ) +
  geom_point(size = 2.5) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    color = "red"
  ) +
  geom_vline(
    xintercept = pooled_OR[2],
    color = "blue",
    linewidth = 0.8
  ) +
  annotate(
    "errorbar",
    xmin = pooled_OR[1],
    xmax = pooled_OR[3],
    y = 0.5,
    width = 0.2,
    color = "blue",
    linewidth = 1
  ) +
  annotate(
    "point",
    x = pooled_OR[2],
    y = 0.5,
    shape = 18,
    size = 3,
    color = "blue"
  ) +
  scale_x_log10() +
  labs(
    title = "Bayesian Random-Effects Meta-Analysis of BCG Vaccination",
    x = "Odds ratio (log scale)",
    y = "Study"
  ) +
  theme_minimal()




## **Shrinkage**

ggplot(forest_data, aes(y = as.numeric(study))) +
  

  geom_errorbar(
    aes(
      xmin = observed_lower,
      xmax = observed_upper,
      y = as.numeric(study) + 0.15
    ),
    orientation = "y",
    width = 0.08,
    linewidth = 0.8,
    color = "orange"
  ) +
  geom_point(
    aes(
      x = observed_OR,
      y = as.numeric(study) + 0.15
    ),
    shape = 1,
    size = 2.8,
    color = "black"
  ) +
  

  geom_errorbar(
    aes(
      xmin = bayes_lower,
      xmax = bayes_upper,
      y = as.numeric(study) - 0.15
    ),
    orientation = "y",
    width = 0.08,
    linewidth = 0.8,
    color = "blue"
  ) +
  geom_point(
    aes(
      x = bayes_OR,
      y = as.numeric(study) - 0.15
    ),
    shape = 18,
    size = 2.8,
    color = "blue"
  ) +
  
  # Null effect
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    color = "red"
  ) +
  
  scale_x_log10() +
  scale_y_continuous(
    breaks = seq_along(levels(forest_data$study)),
    labels = levels(forest_data$study),
    expand = expansion(add = 0.5)
  ) +
  labs(
    title = "Bayesian Shrinkage in Random-Effects Meta-Analysis",
    x = "Odds ratio (log scale)",
    y = "Study"
  ) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank()
  )

##----------------------- END-------------------------------------
