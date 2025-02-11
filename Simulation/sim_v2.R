
library(glmmTMB)
library(broom.mixed)
library(tidyverse)
library(IPMbook)
# Define mean and temporal variability (SD) of the demographic parameters
mean.sj <- 0.3 # Mean juvenile survival (probability scale)
sd.sj.e <- 0.05# Uncertainty of mean juv. survival as SD on natural scale
sd.sj.t <- 0.20 # Temporal variability on the logit scale
mean.sa <- 0.54 # Mean adult survival (probability scale)
sd.sa.e <- 0.05 # Uncertainty of mean ad. survival as SD on natural scale
sd.sa.t <- 0.20 # Temporal variability on the logit scale
mean.f1 <- 2 # Mean productivity of 1y old females
sd.f1.t <- 0.3 # Temporal variability on the natural scale
mean.fa <- 2 # Mean productivity of adult females
sd.fa.e <- 0.03 # Uncertainty of mean productivity as SD on natural scale
sd.fa.t <- 0.3 # Temporal variability on the natural scale
sd.N.e <- 1 # Uncertainty of population size
# Define the number of years with predictions and burn-in length
T <- 5 # Length of Markov chain
u <- 1 # Length of burn-in period
n_sim <- 50 # Number of simulations per area





sim_lirype <- function(n_areas, n_sim, beta){
  
  
  ## Create dataset for storing simulation results
  sim_dat <- tibble(
    Nt = numeric(),
    Nobs = numeric(),
    Na = numeric(),
    Nj = numeric(),
    sa = numeric(),
    sa_obs = numeric(),
    sj = numeric(),
    sj_obs = numeric(),
    fa = numeric(),
    f1 = numeric(),
    r = numeric(),
    robs = numeric(),
    time = numeric(),
    n_sim = numeric(),
    population = factor(),
    latitude = numeric(),
    harvest_rate = factor(),
    n_areas = numeric()
  )

  # Define harvest rate and number of areas
  harvest_rate <- c(0, 0.15, 0.30)
  #harvest_rate_areas <- matrix(sample(harvest_rate, 3), nrow = n_areas, ncol = length(harvest_rate)) # samples harvest rate per year per area, without replacement (each area will experience all harvest rates)
  #harvest_rate_matrix <- cbind(matrix(rep(NA, u), nrow = n_areas, ncol = u, byrow = T), harvest_rate_areas, rep(0, ncol = u)) # samples harvest rate per year per area, without replacement (each area will experience all harvest rates)
  
for (k in 1:n_areas) {
  for (i in 1:n_sim) {
    
    ## Generate areas with different area and densities of ptarmigan
    
    hunting_areas <- tibble(
      area = rnorm(1, 100, 10),
      density = rnorm(1, 20, 4),
      N_ad = (area * density) * 0.4,
      N_juv = (area * density) * 0.6,
      K_density = density * runif(1, 0.8, 4),
      latitude = runif(1, 40, 70)
    )
    ## Draw harvest rates
    harvest_rate_sim <- c(0, sample(harvest_rate, 3), 0)
    
    
    # Generate demographic values from normal distributions
    sj <- plogis(rnorm(T, qlogis(mean.sj), sd.sj.t))
    sa <- plogis(rnorm(T, qlogis(mean.sa), sd.sa.t))
    fa <- f1 <- rnorm(T, mean.fa, sd.fa.t)
    sj_obs <- sa_obs <- c() # vector for survival with observation error
   
    # Define population matrix and initial stage-specific population sizes
    N <- matrix(NA, nrow = 2, ncol = T + 1)
    Nobs <- matrix(NA, nrow = 2, ncol = T + 1) # for storing N with observation error
    N[, 1] <- c(hunting_areas$N_juv, hunting_areas$N_ad)
    
    # Project population forwards
    r <- robs <- numeric(T)
    
    ## Burn-in
    for (t in 1:u) {
      A <- matrix(c(sj[t] * f1[t], sj[t] * fa[t], sa[t], sa[t]), ncol = 2, byrow = TRUE)
      N[, t + 1] <- A %*% N[, t]
      Nobs[, t+1] <- round(rnorm(1, N[, t+1], sd.N.e))
      r[t] <- log(sum(N[, t + 1])) - log(sum(N[, t])) # Annual population growth rate
      robs[t] <- log(sum(Nobs[, t + 1])) - log(sum(Nobs[, t])) # Annual population growth rate with observation error
    }

    # Simulating with harvest
    for (t in (u + 1):T) {
      sj[t] <- sj[t] * (1 - harvest_rate_sim[t]) + (beta * hunting_areas$latitude)
      sa[t] <- sa[t] * (1 - harvest_rate_sim[t]) + (beta * hunting_areas$latitude)
      sj_obs[t] <- rbeta2(1, sj[t], sd.sj.e) ## Adding observation error to sj
      sa_obs[t] <- rbeta2(1, sa[t], sd.sa.e) ## Adding observation error to sa
      
      A <- matrix(c(sj[t] * f1[t], sj[t] * fa[t], sa[t], sa[t]), ncol = 2, byrow = TRUE)
      N[, t + 1] <- A %*% N[, t]
      Nobs[, t+1] <- round(rnorm(1, N[, t+1], sd.N.e))
      r[t] <- log(sum(N[, t + 1])) - log(sum(N[, t])) # Annual population growth rate
      robs[t] <- log(sum(Nobs[, t + 1])) - log(sum(Nobs[, t])) # Annual population growth rate with observation error
      # Add data to dataframe
      sim_dat <- sim_dat %>%
        add_row(
          Nt = sum(N[, t]),
          Nobs = sum(Nobs[, t]),
          Na = N[2, t],
          Nj = N[1, t],
          sa = sa[t],
          sa_obs = sa_obs[t],
          sj = sj[t],
          sj_obs = sj_obs[t],
          f1 = f1[t],
          fa = fa[t],
          r = exp(r[t]),
          robs = exp(robs[t]),
          time = t - u,
          n_sim = i,
          population = as.factor(k),
          harvest_rate = as.factor(harvest_rate_sim[t]),
          latitude = hunting_areas$latitude
        )
      
    }
  }
}
 sim_dat$n_areas <- n_areas
 sim_dat$harvest_rate <- relevel(sim_dat$harvest_rate, ref = "0")
  return(sim_dat)
}

sim_dat <- sim_lirype(n_areas = 10, n_sim = 2, beta = 0.003)
sim_dat

mod_sim <- lm((sa_obs) ~ harvest_rate + latitude  , data = filter(sim_dat, n_sim == 2))
summary(mod_sim)

sim_dat %>% 
  group_by(n_sim) %>%
  reframe(compute_slope_ci(cur_data()), .groups = "drop")
compute_slope_ci(sim_dat)
## Functions for computing effect and CI


# Function to compute slope, p-value, and confidence interval - growth rate
compute_slope_ci_r <- function(data) {
  #data <- data %>% mutate(lead_density = lead(density))  
  
  model <- lm((robs) ~ harvest_rate + latitude, data = data)
  
  # Extract coefficients, confidence intervals, and p-values 
  results <- broom.mixed::tidy(model, conf.int = TRUE)
  
  # Select only the rows for harvest_rate and latitude
  results <- results %>%
    filter(term %in% c("harvest_rate0.15", "harvest_rate0.3", "latitude")) %>%
    select(term, estimate, conf.low, conf.high, p.value)
  
  return(results)
}

# Function to compute slope, p-value, and confidence interval - survival
compute_slope_ci_s <- function(data) {
  #data <- data %>% mutate(lead_density = lead(density))  
  
  model <- glm((sa_obs) ~ harvest_rate + latitude, data = data, family = quasibinomial)
  
  # Extract coefficients, confidence intervals, and p-values 
  results <- broom::tidy(model, conf.int = TRUE)
  
  # Select only the rows for harvest_rate and latitude
  results <- results %>%
    filter(term %in% c("harvest_rate0.15", "harvest_rate0.3", "latitude")) %>%
    select(term, estimate, conf.low, conf.high, p.value)
  
  return(results)
}


#saveRDS(results, file = "Simulation/Results/lirype.RDS")
#results <- readRDS("Simulation/Results/lirype.RDS")


n_areas <- seq(15, 35, 1)
results_list <- list()


for (i in 1:length(n_areas)){
  results_list[[i]] <- sim_lirype(n_areas = n_areas[i],
                                       n_sim = 200,
                                  beta = 0.003) %>% 
    mutate(n_areas = n_areas[i])}

#saveRDS(results_list, file = "Simulation/Results/lirype.RDS") # Save simulation results
results <- readRDS("Simulation/Results/lirype.RDS")

results_s <- bind_rows(results_list) %>%
  filter(time != 4) %>% 
  group_by(n_sim, n_areas) %>%
  reframe(compute_slope_ci_s(cur_data()), .groups = "drop") %>% 
  ungroup() %>% 
  group_by(term, n_areas) %>% 
  summarize(mean_p = mean(p.value),
            prop_p = sum(p.value < 0.05)/length(p.value))
  
saveRDS(results_s, "Report/data/lirype_survival.RDS")
print(results_s, n = Inf)

## Plotting proportion of simulations with p-values above 0.95 - Survival

## Changing names for facet_wrap
facet_labels <- c(
  "harvest_rate0.15" = "Harvest rate = 15 %",
  "harvest_rate0.3" = "Harvest rate = 30 %",
  "latitude" = "Climatic covariate"
)

lirype_s <- results_s %>% 
  filter(term != "(Intercept)") %>% 
  ggplot(aes(x = n_areas, y = prop_p)) + 
  geom_line()+
  geom_hline(yintercept = 0.95, color = "red")+
  facet_wrap(~term, labeller = as_labeller(facet_labels))+
  theme_minimal(base_size = 20)+
  labs(y = "Proportion of simulations with significant effect",
       x = "Number of hunting areas in the simulation",
       title = "Effects of hunting on survival")
lirype_s

## Plotting proportion of simulations with p-values above 0.95 - growth rate
results_r <- bind_rows(results_list) %>%
  #filter(time != 4) %>% 
  group_by(n_sim, n_areas) %>%
  reframe(compute_slope_ci_r(cur_data()), .groups = "drop") %>% 
  ungroup() %>% 
  group_by(term, n_areas) %>% 
  summarize(mean_p = mean(p.value),
            prop_p = sum(p.value < 0.05)/length(p.value))

results_r
saveRDS(results_r, "Report/data/lirype_growth_rate.RDS")
## Changing names for facet_wrap
facet_labels <- c(
  "harvest_rate0.15" = "Harvest rate = 15 %",
  "harvest_rate0.3" = "Harvest rate = 30 %",
  "latitude" = "Climatic covariate"
)

lirype_r <- results_r %>% 
  filter(term != "(Intercept)") %>% 
  ggplot(aes(x = n_areas, y = prop_p)) + 
  geom_line()+
  geom_hline(yintercept = 0.95, color = "red")+
  facet_wrap(~term, labeller = as_labeller(facet_labels))+
  theme_minimal(base_size = 20)+
  labs(y = "Proportion of simulations with significant effect",
       x = "Number of hunting areas in the simulation",
       title = "Effects of hunting on population growth rate")
lirype_r
