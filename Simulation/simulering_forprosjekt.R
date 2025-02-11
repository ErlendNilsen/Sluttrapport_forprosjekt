library(tidyverse)
library(broom)
### Formålet med denne simuleringa er å beregne utvalgsstørrelse for å kunne estimere effekter av høsting og klima gjennom eit høstingseksperiment
#### Eksperimentet går over 3 år - første år er høstingsrate (0, 15 eller 30%) trekt tilfeldig,


harvest_rate <- c(0, 0.15, 0.30)
n_areas <- 100 # number of experimental areas/units
harvest_rate_areas <- matrix(sample(harvest_rate, 3), nrow = n_areas, ncol = length(harvest_rate)) # samples harvest rate per year per area, without replacement (each area will experience all harvest rates)

## Generate areas with different area and densities of ptarmigan

hunting_areas <- tibble(area = runif(n_areas, 50, 150), 
                        density = runif(n_areas, 5, 25),
                        N_ad = (area*density)*0.25,
                        N_juv = (area*density)*0.75,
                        K_density = density * runif(1, 0.8, 4),
                        latitude = runif(n_areas, 40, 70))
hunting_areas

## Function for beta parameters for survival

estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

estBetaParams(0.3, 0.01)

plot(density(rbeta(1000, 6, 14)))

### Funksjon for simulering ----
rype_simulering <- function(harvest_rates,
                            n_areas,
                            n_generations,
                            n_sim,
                            beta_latitude) {
  extinction_time <- c() # vector for saving extinction time
  sim_dat <- tibble(
    N_ad = numeric(),
    N_juv = numeric(),
    p_ad = numeric(),
    p_juv = numeric(),
    m_ad = numeric(),
    density = numeric(),
    time = numeric(),
    n_sim = numeric(),
    #harvest = numeric(),
    extinction_generation = numeric(),
    population = numeric(),
    area = numeric(),
    latitude = numeric(),
    harvest_rate = numeric()
  ) # dataframe for saving data

  for (k in 1:n_areas) {
    # Subsetting data for takseringsområde
    sub_dat <- hunting_areas[k,]
    harvest_rates_sub <- c(0, harvest_rates[k,], 0) ## Adds a 0 in the start because first time-step is used to initialize the population structure - adds a trailing zero because the effect of harvest in time t is on density in time t+1

    # Simulate ptarmigan population dynamics
    for (i in 1:n_sim) {
      extinct <- FALSE
      harvest <- 0
      density <- sub_dat$density # starting density year 1
      p_ad <- p_juv <- m_ad <-  NA


      # Population specific simulation parameters


      
      K <- sub_dat$K_density * sub_dat$area
      N <- list() # List for storing simulated population sizes
      extinction_threshold <- 0.1 * sub_dat$area # Population size below which extinction is assumed, equals density of 0.1/m^2

      # Initial population vector using mean density from estimates as starting point
      N[[1]] <- matrix(c(sub_dat$N_juv, sub_dat$N_ad),
                       ncol = 1)
      N_ad <- N[[1]][2]
      N_juv <- N[[1]][1]

      for (t in 1:(n_generations)) {
        # Vital rates for post-breeding

        m_ad[t] <- rpois(1, lambda = 4) / 2 ## Mean & SD from Indre Salten, number of chicks per pair / 2
        m_juv <- m_ad[t] ## Assume 1yr reproductive rate is equal to adults
        p_juv[t] <- rbeta(1, 6, 14) ## Survival from fledgling to next census. Mean = 0.3, sigma^2 = 0.01
        p_ad[t] <- rbeta(1, 26.2872, 22.3928) ## adult survival from non-harvested ptarmigan (Sandercock et al. 2011): Mean = 0.54. Using variance = 0.005
        
        p_juv[t] <- (p_juv[t] * (1-harvest_rates_sub[t])) + (beta_latitude * (mean(hunting_areas$latitude) - sub_dat$latitude)) ## add effect of latitude on survival, latitude mean-centered
        p_ad[t] <- (p_ad[t] * (1-harvest_rates_sub[t])) + (beta_latitude * (mean(hunting_areas$latitude) - sub_dat$latitude))
        
        f_juv <- m_juv * p_juv[t] ## fertility juv
        f_ad <- m_ad[t] * p_ad[t] ## fertility adult


        # Projection matrix
        A <- matrix(c(f_juv, f_ad, p_juv[t], p_ad[t]), 
                    nrow = 2, byrow = T)
        #D <- (K - sum(N[[t - 1]])) / K
        #I <- diag(2)
        #N[[t]] <- N[[t - 1]] + D * ((A - I) %*% N[[t - 1]]) ## Adding density dependence
        # N[[t]] <- N[[t-1]] +   ((A-I) %*% N[[t-1]])
        N[[t+1]] <- A %*% N[[t]]
        N_ad[t+1] <- ifelse(N[[t+1]][2] < 0, 0, N[[t+1]][2]) ## If negative replace with 0
        N_juv[t+1] <- ifelse(N[[t+1]][1] < 0, 0, N[[t+1]][1])
        density[t+1] <- (N_ad[t+1] + N_juv[t+1])/sub_dat$area

          
        # Check for extinction
        if (sum(N[[t+1]]) < extinction_threshold) {
          extinct <- TRUE
          extinction_generation <- t
          break
        } else {
          extinction_generation <- 0
        }
      }
      # Add data to dataframe
      sim_dat <- sim_dat %>%
        add_row(N_ad = N[[t]][2],
          N_juv,
          density,
          p_ad,
          p_juv,
          m_ad,
          time = seq(1, t),
          n_sim = i,
          extinction_generation,
          population = k,
          harvest_rate = harvest_rates_sub,
          area = sub_dat$area,
          latitude = sub_dat$latitude
        )
    }
  }
  print(sim_dat)
  #saveRDS(sim_dat, paste0(
#    "C:/Users/christoffer.hilde/OneDrive - NINA/Prosjekter/Statskog rypeforvaltning/Data/Forprosjekt/Simuleringer/sim_dat.RDS"
#  ))
}


test_sim <- rype_simulering(harvest_rates = harvest_rate_areas,
                n_generations = 2,
                n_areas = 1,
                n_sim = 1,
                beta_latitude = 0.000) %>% 
  mutate(harvest_rate = as.factor(harvest_rate),
         population = as.factor(population))

library(glmmTMB)

summary(glmmTMB((density) ~ (harvest_rate)  + latitude + m_ad + (1|population) , data = filter(test_sim, n_sim == 1)))
mod <- (glm((p_ad) ~ (harvest_rate)  + latitude , data = filter(test_sim, n_sim == 1), family = quasibinomial(link = "logit")))
summary(mod)
ggpredict(mod)

## Functions for computing effect and CI


# Function to compute slope, p-value, and confidence interval
compute_slope_ci <- function(data) {
  data <- data %>% mutate(lead_density = lead(density))  
  
  model_density <- lm(lead_density ~ harvest_rate + m_ad + latitude, data = data)
  
  # Extract coefficients, confidence intervals, and p-values 
  results <- broom::tidy(model, conf.int = TRUE)
  
  # Select only the rows for harvest_rate and latitude
  results <- results %>%
    filter(term %in% c("(Intercept)", "harvest_rate0.15", "harvest_rate0.3", "latitude")) %>%
    select(term, estimate, conf.low, conf.high, p.value)
  
  return(results)
}

# Apply function per simulation
results <- test_sim %>%
  group_by(n_sim) %>%
  summarise(compute_slope_ci(cur_data()), .groups = "drop")

# View output
head(results)

## summary per sim

results %>% 
  group_by(term) %>% 
  summarize(mean_p = mean(p.value),
            prop_p = sum(p.value < 0.05)/length(p.value))


#### Simulating with number of areas from 10 - 30
n_areas <- seq(20,25,1)
results_list <- list()
for (i in 1:length(n_areas)){
results_list[[i]] <- rype_simulering(harvest_rates = harvest_rate_areas,
                n_generations = 4,
                n_areas = n_areas[i],
                n_sim = 100,
                beta_latitude = 0.005) %>% 
  mutate(harvest_rate = as.factor(harvest_rate),
         population = as.factor(population),
         n_areas = n_areas[i])}

results <- bind_rows(results_list) %>%
  group_by(n_sim, n_areas) %>%
  summarise(compute_slope_ci(cur_data()), .groups = "drop")

results

## Summary per simulation

results <- results %>% 
  group_by(term, n_areas) %>% 
  summarize(mean_p = mean(p.value),
            prop_p = sum(p.value < 0.05)/length(p.value))
results
saveRDS(results, file = "Simulation/Results/lirype.RDS")
results <- readRDS("Simulation/Results/lirype.RDS")
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
svg("Simulation/Results/lirype_resultat.svg",
    width = 10,
    height = 10)
lirype_plot
dev.off()
####################################################
#################### Capercaillie ##################
####################################################

## Survival from Åhlen et al. 2013 - P_ad = 0.68 (CI = 0.62 - 0.75)


## Function for simulating capercaillie harvest

harvest_rate_CC <- c(0, 0.10, 0.20)
n_areas_CC <- 100 # number of experimental areas/units
harvest_rate_areas_CC <- matrix(sample(harvest_rate_CC, 3), nrow = n_areas, ncol = length(harvest_rate_CC)) # samples harvest rate per year per area, without replacement (each area will experience all harvest rates)

## Generate areas with different area and densities of ptarmigan

hunting_areas_CC <- tibble(area = runif(n_areas, 50, 150), 
                        density_ad = runif(n_areas, 4, 15), 
                        density_juv = runif(n_areas, 0, 15), 
                        K_density = (density_ad + density_juv) * runif(1, 0.5, 4),
                        latitude = runif(n_areas, 40, 70))
hunting_areas_CC


### Funksjon for simulering ----
CC_simulering <- function(harvest_rates,
                            n_areas,
                            n_generations,
                            n_sim,
                            beta_latitude) {
  extinction_time <- c() # vector for saving extinction time
  sim_dat <- tibble(
    N_ad = numeric(),
    N_juv = numeric(),
    p_ad = numeric(),
    p_juv = numeric(),
    density = numeric(),
    time = numeric(),
    n_sim = numeric(),
    harvest = numeric(),
    extinction_generation = numeric(),
    population = numeric(),
    area = numeric(),
    latitude = numeric(),
    harvest_rate = numeric()
  ) # dataframe for saving data
  
  for (k in 1:n_areas) {
    # Subsetting data for takseringsområde
    sub_dat <- hunting_areas[k,]
    harvest_rates_sub <- c(0, harvest_rates[k,], 0) ## Adds a 0 in the start because first time-step is used to initialize the population structure - adds a trailing zero because the effect of harvest in time t is on density in time t+1
    
    # Simulate ptarmigan population dynamics
    for (i in 1:n_sim) {
      extinct <- FALSE
      harvest <- 0
      density <- sub_dat$density_ad + sub_dat$density_juv # starting density year 1
      p_ad <- p_juv <- NA
      
      
      # Population specific simulation parameters
      
      
      
      K <- sub_dat$K_density * sub_dat$area
      N <- list() # List for storing simulated population sizes
      extinction_threshold <- 0.1 * sub_dat$area # Population size below which extinction is assumed, equals density of 0.1/m^2
      
      # Initial population vector using mean density from estimates as starting point
      N[[1]] <- c(round(sub_dat$area * sub_dat$density_juv), round(sub_dat$area * sub_dat$density_ad))
      N_ad <- N[[1]][2]
      N_juv <- N[[1]][1]
      
      for (t in 2:(n_generations+1)) {
        # Vital rates for post-breeding
        
        m_ad <- rpois(1, lambda = 2.5)
        m_juv <- m_ad ## Assume 1yr reproductive rate is equal to adults
        p_juv[t] <- rbeta(1, 18.8, 28.2) ## Survival from fledgling to next census. Mean = 0.3, sigma^2 = 0.01
        p_ad[t] <- rbeta(1, 28.9136, 13.6064) ## adult survival from Åhlgren et al. 2011 : Mean = 0.68 Using variance = 0.005
        
        p_juv[t] <-p_juv[t] + beta_latitude * (mean(hunting_areas$latitude) - sub_dat$latitude) ## add effect of latitude on survival, latitude mean-centered
        p_ad[t] <- p_ad[t] + beta_latitude * (mean(hunting_areas$latitude) - sub_dat$latitude)
        
        f_juv <- m_juv * p_juv[t] ## fertility juv
        f_ad <- m_ad * p_ad[t] ## fertility adult
        
        
        # Projection matrix
        A <- matrix(c(f_juv, f_ad, p_juv[t], p_ad[t]), nrow = 2, byrow = T)
        #D <- (K - sum(N[[t - 1]])) / K
        I <- diag(2)
        #N[[t]] <- N[[t - 1]] + D * ((A - I) %*% N[[t - 1]]) ## Adding density dependence
        # N[[t]] <- N[[t-1]] +   ((A-I) %*% N[[t-1]])
        N[[t]] <- A %*%N[[t-1]]
        N_ad[t] <- ifelse(N[[t]][2] < 0, 0, N[[t]][2]) ## If negative replace with 0
        N_juv[t] <- ifelse(N[[t]][1] < 0, 0, N[[t]][1])
        density[t] <- (N_ad[t] + N_juv[t])/sub_dat$area
        
        # Apply harvest and update population - assume equal harvest probability for adults and juveniles
        
        harvest[t] <- (N_ad[t] + N_juv[t]) * harvest_rates_sub[t]
        
        #age_ratio <- N_ad[t] / (N_juv[t] + N_ad[t])
        #N[[t]][2] <- N_ad[t] - (harvest[t] * age_ratio)
        #N[[t]][1] <- N_juv[t] - (harvest[t] * (1 - age_ratio))
        
        N[[t]][2] <- N_ad[t] * (1 - harvest_rates_sub[t])
        N[[t]][1] <- N_juv[t] * (1 - harvest_rates_sub[t])
        
        
        # Check for extinction
        if (sum(N[[t]]) < extinction_threshold) {
          extinct <- TRUE
          extinction_generation <- t
          break
        } else {
          extinction_generation <- 0
        }
      }
      # Add data to dataframe
      sim_dat <- sim_dat %>%
        add_row(N_ad,
                N_juv,
                density,
                p_ad,
                p_juv,
                time = seq(1, t),
                n_sim = i,
                harvest,
                extinction_generation,
                population = k,
                harvest_rate = harvest_rates_sub,
                area = sub_dat$area,
                latitude = sub_dat$latitude
        )
    }
  }
  print(sim_dat)
  #saveRDS(sim_dat, paste0(
  #    "C:/Users/christoffer.hilde/OneDrive - NINA/Prosjekter/Statskog rypeforvaltning/Data/Forprosjekt/Simuleringer/sim_dat.RDS"
  #  ))
}

test_sim_CC <- CC_simulering(harvest_rates = harvest_rate_areas_CC,
                            n_generations = 4,
                            n_areas = 30,
                            n_sim = 50,
                            beta_latitude = 0.005) %>% 
  mutate(harvest_rate = as.factor(harvest_rate),
         population = as.factor(population))


summary(lm(lead(density) ~ (harvest_rate)  + latitude  , data = filter(test_sim_CC, n_sim == 1)))

#Function to compute slope, p-value, and confidence interval for Capercaillie
compute_slope_ci_CC <- function(data) {
  data <- data %>% mutate(lead_density = lead(density))  
  
  model <- lm(lead_density ~ harvest_rate + latitude, data = data)
  
  # Extract coefficients, confidence intervals, and p-values 
  results <- broom::tidy(model, conf.int = TRUE)
  
  # Select only the rows for harvest_rate and latitude
  results <- results %>%
    filter(term %in% c("(Intercept)", "harvest_rate0.1", "harvest_rate0.2", "latitude")) %>%
    select(term, estimate, conf.low, conf.high, p.value)
  
  return(results)
}


results_CC <- test_sim_CC %>%
  group_by(n_sim) %>%
  summarise(compute_slope_ci_CC(cur_data()), .groups = "drop")

results_CC %>% 
  group_by(term) %>% 
  summarize(mean_p = mean(p.value),
            prop_p = sum(p.value < 0.05)/length(p.value))


#### Simulating with number of areas from 10 - 30
n_areas <- c(seq(10,30, 1))
results_list_CC <- list()
for (i in 1:length(n_areas)){
  results_list_CC[[i]] <- CC_simulering(harvest_rates = harvest_rate_areas_CC,
                                       n_generations = 4,
                                       n_areas = n_areas[i],
                                       n_sim = 100,
                                       beta_latitude = 0.005) %>% 
    mutate(harvest_rate = as.factor(harvest_rate),
           population = as.factor(population),
           n_areas = n_areas[i])}

#
results_CC <- bind_rows(results_list_CC) %>%
  group_by(n_sim, n_areas) %>%
  summarise(compute_slope_ci_CC(cur_data()), .groups = "drop")

results_CC

## Summary per simulation

results_CC <- results_CC %>% 
  group_by(term, n_areas) %>% 
  summarize(mean_p = mean(p.value),
            prop_p = sum(p.value < 0.05)/length(p.value))

saveRDS(results_CC, file = "Simulation/Results/skogsfugl.RDS")

## Plotting proportion of simulations with p-values above 0.95

facet_labels <- c(
  "harvest_rate0.1" = "Høstingsrate = 10 %",
  "harvest_rate0.2" = "Høstingsrate = 20 %",
  "latitude" = "Kovariat (f.eks klimaeffekt)"
)

skogsfugl_plot <- results_CC %>% 
  filter(term != "(Intercept)") %>% 
  ggplot(aes(x = n_areas, y = prop_p)) + 
  geom_line()+
  geom_hline(yintercept = 0.95, color = "red")+
  facet_wrap(~term, labeller = as_labeller(facet_labels))+
  theme_minimal(base_size = 20)+
  labs(y = "Andel simuleringer med p < 0.05",
       x = "Antall områder i simuleringen",
       title = "Effekt av høsting og kovariat på tetthet - Skogsfugl")

skogsfugl_plot

svg("Simulation/Results/skogsfugl_resultat.svg",
    width = 10,
    height = 10)
skogsfugl_plot
dev.off()



