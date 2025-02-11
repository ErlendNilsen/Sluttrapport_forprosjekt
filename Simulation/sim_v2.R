# Define mean and temporal variability (SD) of the demographic parameters
mean.sj <- 0.3 # Mean juvenile survival (probability scale)
sd.sj.t <- 0.25 # Temporal variability on the logit scale
mean.sa <- 0.54 # Mean adult survival (probability scale)
sd.sa.t <- 0.07 # Temporal variability on the logit scale
mean.f1 <- 2 # Mean productivity of 1y old females
sd.f1.t <- 0.3 # Temporal variability on the natural scale
mean.fa <- 2 # Mean productivity of adult females
sd.fa.t <- 0.3 # Temporal variability on the natural scale
# Define the number of years with predictions and burn-in length
T <- 5 # Length of Markov chain
u <- 1 # Length of burn-in period
n_sim <- 50 # Number of simulations per area

harvest_rate_areas


## Generate areas with different area and densities of ptarmigan

hunting_areas <- tibble(
  area = runif(n_areas, 50, 150),
  density = runif(n_areas, 5, 25),
  N_ad = (area * density) * 0.25,
  N_juv = (area * density) * 0.75,
  K_density = density * runif(1, 0.8, 4),
  latitude = runif(n_areas, 40, 70)
)
hunting_areas

## Create dataset for storing simulation results
sim_dat <- tibble(
  Na = numeric(),
  Nj = numeric(),
  sa = numeric(),
  sj = numeric(),
  fa = numeric(),
  f1 = numeric(),
  r = numeric(),
  # density = numeric(),
  time = numeric(),
  n_sim = numeric(),
  population = factor(),
  # area = numeric(),
  # latitude = numeric(),
  harvest_rate = factor(),
  n_areas = numeric()
)

sim_lirype <- function(n_areas, n_sim){

  # Define harvest rate and number of areas
  harvest_rate <- c(0, 0.15, 0.30)
  harvest_rate_areas <- matrix(sample(harvest_rate, 3), nrow = n_areas, ncol = length(harvest_rate)) # samples harvest rate per year per area, without replacement (each area will experience all harvest rates)
  harvest_rate_matrix <- cbind(matrix(rep(NA, u), nrow = n_areas, ncol = u, byrow = T), harvest_rate_areas, rep(0, ncol = u)) # samples harvest rate per year per area, without replacement (each area will experience all harvest rates)
  
for (k in 1:n_areas) {
  for (i in 1:n_sim) {
    # Generate demographic values from normal distributions
    sj <- plogis(rnorm(T, qlogis(mean.sj), sd.sj.t))
    sa <- plogis(rnorm(T, qlogis(mean.sa), sd.sa.t))
    f1 <- rnorm(T, mean.f1, sd.f1.t)
    fa <- rnorm(T, mean.fa, sd.fa.t)
    # Define population matrix and initial stage-specific population sizes
    N <- matrix(NA, nrow = 2, ncol = T + 1)
    N[, 1] <- c(hunting_areas$N_juv[k], hunting_areas$N_ad[k])
    # Project population forwards
    r <- numeric(T)
    ## Burn-in
    for (t in 1:u) {
      A <- matrix(c(sj[t] * f1[t], sj[t] * fa[t], sa[t], sa[t]), ncol = 2, byrow = TRUE)
      N[, t + 1] <- A %*% N[, t]
      r[t] <- log(sum(N[, t + 1])) - log(sum(N[, t])) # Annual population growth rate
    }

    # Simulating with harvest
    for (t in (u + 1):T) {
      sj[t] <- sj[t] * (1 - harvest_rate_matrix[1, t])
      sa[t] <- sa[t] * (1 - harvest_rate_matrix[1, t])
      A <- matrix(c(sj[t] * f1[t], sj[t] * fa[t], sa[t], sa[t]), ncol = 2, byrow = TRUE)
      N[, t + 1] <- A %*% N[, t]
      r[t] <- log(sum(N[, t + 1])) - log(sum(N[, t])) # Annual population growth rate

      # Add data to dataframe
      sim_dat <- sim_dat %>%
        add_row(
          Na = N[2, t],
          Nj = N[1, t],
          # density,
          sa = sa[t],
          sj = sj[t],
          f1 = f1[t],
          fa = fa[t],
          r = exp(r[t]),
          time = t - u,
          n_sim = i,
          population = as.factor(k),
          harvest_rate = as.factor(harvest_rate_matrix[k, t])
          # area = sub_dat$area,
          # latitude = sub_dat$latitude
        )
      
    }
  }
}
 sim_dat$n_areas <- n_areas
  return(sim_dat)
}

sim_dat <- sim_lirype(n_areas = 21, n_sim = 2)

sim_dat

mod_sim <- lm(sa ~ harvest_rate + population, data = filter(sim_dat, n_sim == 1))
summary(mod_sim)

## Functions for computing effect and CI


# Function to compute slope, p-value, and confidence interval
compute_slope_ci <- function(data) {
  #data <- data %>% mutate(lead_density = lead(density))  
  
  model <- lm(sa ~ harvest_rate + population, data = data)
  
  # Extract coefficients, confidence intervals, and p-values 
  results <- broom::tidy(model, conf.int = TRUE)
  
  # Select only the rows for harvest_rate and latitude
  results <- results %>%
    filter(term %in% c("(Intercept)", "harvest_rate0.15", "harvest_rate0.3")) %>%
    select(term, estimate, conf.low, conf.high, p.value)
  
  return(results)
}

#saveRDS(results, file = "Simulation/Results/lirype.RDS")
#results <- readRDS("Simulation/Results/lirype.RDS")


n_areas <- seq(10,30,1)
results_list <- list()


for (i in 1:length(n_areas)){
  results_list[[i]] <- sim_lirype(n_areas = n_areas[i],
                                       n_sim = 100) %>% 
    mutate(n_areas = n_areas[i])}

results <- bind_rows(results_list) %>%
  group_by(n_sim, n_areas) %>%
  summarise(compute_slope_ci(cur_data()), .groups = "drop") %>% 
  ungroup() %>% 
  group_by(term, n_areas) %>% 
  summarize(mean_p = mean(p.value),
            prop_p = sum(p.value < 0.05)/length(p.value))
  

results

## Plotting proportion of simulations with p-values above 0.95

## Changing names for facet_wrap
facet_labels <- c(
  "harvest_rate0.15" = "Høstingsrate = 15 %",
  "harvest_rate0.3" = "Høstingsrate = 30 %",
  "latitude" = "Kovariat (f.eks klimaeffekt)"
)

lirype_plot <- results %>% 
  filter(term != "(Intercept)") %>% 
  ggplot(aes(x = n_areas, y = prop_p)) + 
  geom_line()+
  geom_hline(yintercept = 0.95, color = "red")+
  facet_wrap(~term, labeller = as_labeller(facet_labels))+
  theme_minimal(base_size = 20)+
  labs(y = "Andel simuleringer med p < 0.05",
       x = "Antall områder i simuleringen",
       title = "Effekt av høsting og kovariat på tetthet - Lirype")
lirype_plot
