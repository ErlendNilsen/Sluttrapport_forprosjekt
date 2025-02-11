

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