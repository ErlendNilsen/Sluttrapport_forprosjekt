

sp_data<-included_full_details |>
  separate_rows(species_if_multiple_seperate_with_a_comma_use_full_latin_names_e_g_lagopus_lagopus_lagopus_muta, sep = ", ") |>
  mutate(
    species_corrected = map_chr(
      species_if_multiple_seperate_with_a_comma_use_full_latin_names_e_g_lagopus_lagopus_lagopus_muta,
      ~ {
        res <- tryCatch(name_backbone(name = .x), error = function(e) NULL)
        if (!is.null(res) && "species" %in% names(res) && !is.na(res$species)) {
          res$species # Return corrected name
        } else {
          .x # Return original if not found
        }
      }
    )
  ) 


sp_EBV<-sp_data |> 
  separate_rows(outcome_variable_ebv_name_not_exhaustive_indicate_name_that_is_not_included_in_original_framework, sep = ", ") |>
  group_by(outcome_variable_ebv_name_not_exhaustive_indicate_name_that_is_not_included_in_original_framework) |>
  mutate(outcome_variable_ebv_name_not_exhaustive_indicate_name_that_is_not_included_in_original_framework = case_when(
    outcome_variable_ebv_name_not_exhaustive_indicate_name_that_is_not_included_in_original_framework %in% c("Behaviour", "Display behaviour/ Gobbling", "Foraging behaviour", "Gobbling activity") ~ "Behaviour",
    TRUE ~ outcome_variable_ebv_name_not_exhaustive_indicate_name_that_is_not_included_in_original_framework
  )) 
# Aggregate data to get counts
heatmap_data <- sp_EBV |> 
  count(outcome_variable_ebv_name_not_exhaustive_indicate_name_that_is_not_included_in_original_framework, species_corrected) |> 
  drop_na()

# Rename columns for easier use in ggplot
colnames(heatmap_data) <- c("EBVname", "Species", "Count")
ggplot(heatmap_data, aes(x = EBVname, y = Species, fill = Species)) +
  geom_tile(aes(alpha=Count) )+
  scale_fill_manual(values = species_colours) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels for readability
  labs(title = "Species Distribution by EBV name",
       x = "EBV name",
       y = "Species",
       fill="")

################################################################################

included_full_details |>
  separate_rows(species_if_multiple_seperate_with_a_comma_use_full_latin_names_e_g_lagopus_lagopus_lagopus_muta, sep = ", ") |>
  mutate(
    species_corrected = map_chr(
      species_if_multiple_seperate_with_a_comma_use_full_latin_names_e_g_lagopus_lagopus_lagopus_muta,
      ~ {
        res <- tryCatch(name_backbone(name = .x), error = function(e) NULL)
        if (!is.null(res) && "species" %in% names(res) && !is.na(res$species)) {
          res$species # Return corrected name
        } else {
          .x # Return original if not found
        }
      }
    )
  ) |>
  group_by(species_corrected) |>
  tally() |>
  drop_na() |>
  ggplot(aes(reorder(species_corrected, n), n)) +
  geom_col(fill=species_colours[3]) + # Changed from geom_histogram(stat="identity") to geom_col()
  labs(
    x = "Species",
    y = "Number of Responses"
  ) +
  theme_minimal() +
  coord_flip() + 
  theme(text = element_text(size = 15))


################################################################################

country_data<-included_full_details |>
  separate_rows(study_area_location_country_e_g_sweden_if_multiple_seperate_with_a_comma, sep = ", ") |> # Split at comma and create new rows
  group_by(study_area_location_country_e_g_sweden_if_multiple_seperate_with_a_comma) |>
  mutate(study_area_locat
         ion_country_e_g_sweden_if_multiple_seperate_with_a_comma = case_when(
           study_area_location_country_e_g_sweden_if_multiple_seperate_with_a_comma %in% c(
             "Alaska", "California", "Colorado", "Florida", "Georgia","Georiga", "Georgia (USA)", "Idaho", "Illinois", "Indiana", "Kansas", "Kentucky", "Louisiana", "Maine", "Michigan", "Mississippi", "Missouri", "Montana", "Nevada", "New York", "North Carolina", "Ohio",
             "Oklahoma", "Oregon", "Pennsylvania", "Poland", "Rhode Island", "Rhode island", "South Caolina", "South Caraolina", "South Carolina", "South Dakota", "Tennessee", "Texas", "Texas,Texas", "US", "Utah",
             "Virgina", "Virginia", "West Virginia", "Wisconsin", "Wisonsin"
           ) ~ "USA",
           TRUE ~ study_area_location_country_e_g_sweden_if_multiple_seperate_with_a_comma
         ))

country_sp_data<-country_data |> 
  separate_rows(species_if_multiple_seperate_with_a_comma_use_full_latin_names_e_g_lagopus_lagopus_lagopus_muta, sep = ", ") |>
  mutate(
    species_corrected = map_chr(
      species_if_multiple_seperate_with_a_comma_use_full_latin_names_e_g_lagopus_lagopus_lagopus_muta,
      ~ {
        res <- tryCatch(name_backbone(name = .x), error = function(e) NULL)
        if (!is.null(res) && "species" %in% names(res) && !is.na(res$species)) {
          res$species # Return corrected name
        } else {
          .x # Return original if not found
        }
      }
    )
  ) 

# Aggregate data to get counts
heatmap_data <- country_sp_data %>%
  count(study_area_location_country_e_g_sweden_if_multiple_seperate_with_a_comma, species_corrected)|> 
  drop_na()

colnames(heatmap_data) <- c("Country", "Species", "Count")


ggplot(heatmap_data, aes(x = Country, y = Species, size = Count)) +
  geom_point(colour = species_colours[3], alpha = 0.8)+
  scale_fill_manual(values = species_colours) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels for readability
  labs(x = "Country",
       y = "Species",
       fill = "") +
  theme(legend.position = "none", legend.title = element_blank()) +
  theme(text = element_text(size = 15))

```

Different species were also subject to research focusing on different outcome variables (@fig-species_ebv). \[More text...\]

```{r}
#| echo: false
#| warning: false
#| message: false
#| label: fig-species_ebv
#| fig-cap: "The distribution of EBV classes by species"

sp_data<-included_full_details |>
  separate_rows(species_if_multiple_seperate_with_a_comma_use_full_latin_names_e_g_lagopus_lagopus_lagopus_muta, sep = ", ") |>
  mutate(
    species_corrected = map_chr(
      species_if_multiple_seperate_with_a_comma_use_full_latin_names_e_g_lagopus_lagopus_lagopus_muta,
      ~ {
        res <- tryCatch(name_backbone(name = .x), error = function(e) NULL)
        if (!is.null(res) && "species" %in% names(res) && !is.na(res$species)) {
          res$species # Return corrected name
        } else {
          .x # Return original if not found
        }
      }
    )
  ) 


sp_EBV<-sp_data |> 
  separate_rows(outcome_variable_ebv_name_not_exhaustive_indicate_name_that_is_not_included_in_original_framework, sep = ", ") |>
  group_by(outcome_variable_ebv_name_not_exhaustive_indicate_name_that_is_not_included_in_original_framework) |>
  mutate(outcome_variable_ebv_name_not_exhaustive_indicate_name_that_is_not_included_in_original_framework = case_when(
    outcome_variable_ebv_name_not_exhaustive_indicate_name_that_is_not_included_in_original_framework %in% c("Behaviour", "Display behaviour/ Gobbling", "Foraging behaviour", "Gobbling activity") ~ "Behaviour",
    TRUE ~ outcome_variable_ebv_name_not_exhaustive_indicate_name_that_is_not_included_in_original_framework
  )) 

# Aggregate data to get counts
heatmap_data <- sp_EBV |>
  mutate(outcome_variable_ebv_name_not_exhaustive_indicate_name_that_is_not_included_in_original_framework = case_when(
    outcome_variable_ebv_name_not_exhaustive_indicate_name_that_is_not_included_in_original_framework %in% c("Population density *", "Species abundance", "Population growth rate") ~ "Population density/abundance/growth*",
    TRUE ~ outcome_variable_ebv_name_not_exhaustive_indicate_name_that_is_not_included_in_original_framework
  )) |>
  count(outcome_variable_ebv_name_not_exhaustive_indicate_name_that_is_not_included_in_original_framework, species_corrected) |> 
  drop_na()

# Rename columns for easier use in ggplot
colnames(heatmap_data) <- c("EBVname", "Species", "Count")

ggplot(heatmap_data, aes(x = EBVname, y = Species, size = Count)) +
  geom_point(colour = species_colours[3], alpha = 0.8)+
  scale_fill_manual(values = species_colours) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels for readability
  labs(x = "EBV name",
       y = "Species",
       fill="") +
  theme(legend.position = "none", legend.title = element_blank()) +
  theme(text = element_text(size = 15))





