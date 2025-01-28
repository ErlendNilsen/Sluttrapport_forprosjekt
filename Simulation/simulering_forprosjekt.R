library(tidyverse)

### Formålet med denne simuleringa er å beregne utvalgsstørrelse for å kunne estimere effekter av høsting og klima gjennom eit høstingseksperiment
#### Eksperimentet går over 3 år - første år er høstingsrate (0, 15 eller 30%) trekt tilfeldig,


harvest_rate <- c(0, 0.15, 0.30)
n_areas <- 10 # number of experimental areas/units
harvest_rate_areas <- matrix(sample(harvest_rate, 3), nrow = n_areas, ncol = length(harvest_rate)) # samples harvest rate per year per area, without replacement (each area will experience all harvest rates)

## Generate areas with different area and densities of ptarmigan

hunting_areas <- tibble(area = runif(n_areas, 50, 150), density_ad = runif(n_areas, 5, 15), density_juv = runif(n_areas, 0, 20), K_density = (density_ad + density_juv) * runif(1, 0.5, 2))
hunting_areas

### Funksjon for simulering ----
rype_simulering <- function(harvest_rates,
                            n_areas,
                            n_generations,
                            n_sim) {
  extinction_time <- c() # vector for saving extinction time
  sim_dat <- tibble(
    N_ad = numeric(),
    N_juv = numeric(),
    time = numeric(),
    n_sim = numeric(),
    harvest = numeric(),
    extinction_generation = numeric(),
    population = numeric(),
    area = numeric(),
    #K_density = numeric(),
    harvest_rate = numeric()
  ) # dataframe for saving data

  for (k in 1:n_areas) {
    # Subsetting data for takseringsområde
    sub_dat <- hunting_areas[k,]
    harvest_rates_sub <- c(0, harvest_rates[k,]) ## Adds a 0 in the start because first time-step is used to initialize the population structure

    # Simulate ptarmigan population dynamics
    for (i in 1:n_sim) {
      extinct <- FALSE
      harvest <- 0


      # Population specific simulation parameters


      
      K <- sub_dat$K_density * sub_dat$area
      N <- list() # List for storing simulated population sizes
      extinction_threshold <- 0.5 * sub_dat$area # Population size below which extinction is assumed, 350 equals density of 0.5/m^2

      # Initial population vector using mean density from estimates as starting point
      N[[1]] <- c(round(sub_dat$area * sub_dat$density_juv), round(sub_dat$area * sub_dat$density_ad))
      N_ad <- N[[1]][2]
      N_juv <- N[[1]][1]

      for (t in 2:n_generations) {
        # Vital rates for post-breeding

        m_ad <- rpois(1, lambda = 4) / 2 ## Mean & SD from Indre Salten, number of chicks per pair / 2
        m_juv <- m_ad ## Assume 1yr reproductive rate is equal to adults
        p_juv <- rbeta(1, 2.85, 6.65) ## Survival from fledgling to next census. Survival = 0.3, sigma^2 = 0.02
        p_ad <- rbeta(1, 6.1668, 5.2532) ## adult survival from non-harvested ptarmigan (Sandercock et al. 2011): 0.54. Using variance = 0.02
        f_juv <- m_juv * p_juv ## fertility juv
        f_ad <- m_ad * p_ad ## fertility adult


        # Projection matrix
        A <- matrix(c(f_juv, f_ad, p_juv, p_ad), nrow = 2, byrow = T)
        D <- (K - sum(N[[t - 1]])) / K
        I <- diag(2)
        N[[t]] <- N[[t - 1]] + D * ((A - I) %*% N[[t - 1]]) ## Adding density dependence
        # N[[t]] <- N[[t-1]] +   ((A-I) %*% N[[t-1]])
        # N[[t]] <- A %*%N[[t-1]]
        N_ad[t] <- ifelse(N[[t]][2] < 0, 0, N[[t]][2]) ## If negative replace with 0
        N_juv[t] <- ifelse(N[[t]][1] < 0, 0, N[[t]][1])


        harvest[t] <- (N_ad[t] + N_juv[t]) * harvest_rates_sub[t]



        # Apply harvest and update population - assume equal harvest probability for adults and juveniles
        age_ratio <- N_ad[t] / (N_juv[t] + N_ad[t])
        N[[t]][2] <- N_ad[t] - (harvest[t] * age_ratio)
        N[[t]][1] <- N_juv[t] - (harvest[t] * (1 - age_ratio))

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
          time = seq(1, t),
          n_sim = i,
          harvest,
          extinction_generation,
          population = k,
          harvest_rate = harvest_rates_sub,
          area = sub_dat$area
        )
    }
  }
  print(sim_dat)
#  saveRDS(sim_dat, paste0(
#    "C:/Users/christoffer.hilde/OneDrive - NINA/Prosjekter/Statskog rypeforvaltning/Data/Forprosjekt/Simuleringer/sim_dat.RDS"
#  ))
}


rype_simulering(harvest_rates = harvest_rate_areas,
                n_generations = 4,
                n_areas = 4,
                n_sim = 1)
                