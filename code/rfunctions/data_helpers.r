# Function to calculate DayRelativeToNearestHCT for infections data
calculate_day_relative_to_nearest_hct <- function(patient_id, event_day, hct_meta) {
    patient_hcts <- hct_meta %>%
        filter(PatientID == patient_id) %>%
        pull(TimepointOfTransplant)
    nearest_hct <- patient_hcts[which.min(abs(event_day - patient_hcts))]
    event_day - nearest_hct
}

# Define a custom theme for plots
custom_theme <- theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.75), 
    axis.text = element_text(color = "black", size = 10)
)

# Function to prepare taxa counts for a given taxonomic level. This will join the counts with the patient metadata.
#' @param taxlevel A string indicating the taxonomic level to filter by (e.g., "asv", "genus", "species").
#' @return A data frame containing the prepared taxa counts with metadata.
prep_taxa_counts <- function(taxlevel) {
   
   # Check caps to use taxlevel for filtering
   taxlevel_filter <- if (taxlevel == "asv") {
        toupper(taxlevel)
    } else {
        tools::toTitleCase(tolower(taxlevel))
    }
    
    # Transplant metadata
    tblhctmeta <- suppressMessages(read_csv("data/tblhctmeta.csv"))
    
    # Load stool sample metadata
    tblASVsamples <- suppressMessages(
        read_csv("data/tblASVsamples.csv",
            col_types = cols(
                .default = col_guess(),  # Automatically guess other columns
                Pool = col_character())) %>%  # Force Pool column to be character
            filter(PatientID %in% unique(tblhctmeta$PatientID)) %>%
            select(SampleID, PatientID, Consistency, Timepoint) %>%
            rowwise() %>%
            mutate(DayRelativeToNearestHCT = calculate_day_relative_to_nearest_hct(PatientID, Timepoint, tblhctmeta)) %>%
            ungroup() %>%
            distinct(SampleID, .keep_all = TRUE)
    )
    
    # Load clinical metadata and score patients who developed infections
    load("analyses/processed_data/tblpatientday_clinical.RData")
    patient_infection_status <- tblpatientday_clinical %>%
        group_by(PatientID) %>%
        summarise(across(ends_with("_infection"), ~ any(.x == 1), .names = "{col}"))
    
    # Load and transform the counts for the specified taxonomic level
    tblrel_taxlevel_wide <- suppressMessages(
        read_csv(paste0("data/tblcounts_", taxlevel, "_wide.csv")) %>%
            pivot_longer(
                cols = -1,
                names_to = "Sample",
                values_to = "Count"
            ) %>%
            pivot_wider(names_from = all_of(taxlevel_filter), values_from = "Count") %>%
            mutate(count_sum = rowSums(across(2:ncol(.)))) %>%
            mutate(across(2:(ncol(.)-1), ~ .x / count_sum, .names = "{col}_abund")) %>%
            select(Sample, ends_with("_abund"))
    )
    
    # Join counts with stool sample metadata
    tblrel_taxlevel_meta <- tblrel_taxlevel_wide %>%
            left_join(tblASVsamples, by = c("Sample" = "SampleID")) %>%
            filter(!is.na(PatientID)) %>%
            left_join(patient_infection_status, by = "PatientID") %>%
            mutate(across(ends_with("_infection"), ~ replace_na(.x, FALSE)))
    
    return(tblrel_taxlevel_meta)
}


# Function to plot microbial composition timeline for a given PatientID
plot_patient_timeline_asv <- function(patient_id, tblrel_asv_meta) {
    
    # Load taxonomy data with colors and filter to tax level
    taxlevel_color <- read_csv("data/tblASVtaxonomy_silva132_v4v5_filter.csv") %>%
        select(ASV, HexColor, ColorOrder) %>%
        distinct(ASV, .keep_all = TRUE)
    
    # Filter data for the given PatientID and time points between -5 and 20
    patient_data <- tblrel_asv_meta %>%
        filter(PatientID == patient_id, DayRelativeToNearestHCT >= -30, DayRelativeToNearestHCT <= 100) %>%
        select(DayRelativeToNearestHCT, ends_with("_abund")) %>%
        pivot_longer(
            cols = ends_with("_abund"),
            names_to = "Taxa",
            values_to = "RelativeAbundance"
        ) %>%
        mutate(Taxa = str_remove(Taxa, "_abund"))
    
    # Merge patient data with taxonomy colors and order by ColorOrder
    patient_data <- patient_data %>%
        left_join(taxlevel_color, by = c("Taxa" = "ASV")) %>%
        arrange(ColorOrder)
    
    # Plot microbial composition for the given PatientID with assigned colors and order
    ggplot(patient_data, aes(x = DayRelativeToNearestHCT, y = RelativeAbundance, fill = reorder(Taxa, -ColorOrder))) +
        geom_bar(stat = "identity", position = "stack") +
        scale_fill_manual(values = setNames(patient_data$HexColor, patient_data$Taxa)) +
        labs(
            title = paste("Microbial Composition of Patient", patient_id),
            x = "Days Relative to Nearest HCT",
            y = "16S relative abundance"
        ) +
        custom_theme +
        theme(legend.position = "none")
}

# Fucntion to plot drug administration timeline for a given PatientID
# This function assumes that tbldrug is a data frame with columns: PatientID, Factor, StartDayRelativeToNearestHCT, StopDayRelativeToNearestHCT
plot_patient_drug_timeline <- function(patient_id, tbldrug, tblhctmeta, tblrel_asv_meta) {
  
  # Get the bounds of the X axis from tblrel_asv_meta
  x_bounds <- tblrel_asv_meta %>%
    filter(PatientID == patient_id) %>%
    summarise(
      min_day = min(DayRelativeToNearestHCT, na.rm = TRUE),
      max_day = max(DayRelativeToNearestHCT, na.rm = TRUE)
    )
  
  # Filter tblhctmeta to include only patients with a single HCT
  single_hct_patients <- tblhctmeta %>%
    group_by(PatientID) %>%
    filter(n() == 1) %>%
    pull(PatientID)
  
  # Check if the patient is in the single HCT list
  if (!(patient_id %in% single_hct_patients)) {
    message("The patient had multiple HCTs and cannot plot the timeline.")
    return(NULL)
  }
  
  # Filter tbldrug for the given patient
  patient_drug_data <- tbldrug %>%
    filter(PatientID == patient_id)
  
  # Check if there is ASV data for the patient
  if (nrow(x_bounds) == 0 || is.na(x_bounds$min_day) || is.na(x_bounds$max_day)) {
    message("No ASV data available for the patient to define X axis bounds.")
    return(NULL)
  }
  
  # Create the dumbbell plot using geom_segment
  plot <- ggplot(patient_drug_data, aes(y = Factor, x = StartDayRelativeToNearestHCT, xend = StopDayRelativeToNearestHCT)) +
    geom_segment(aes(x = StartDayRelativeToNearestHCT, xend = StopDayRelativeToNearestHCT, 
                     y = Factor, yend = Factor), linewidth = 3.5, color = "gray75") +
    scale_x_continuous(limits = c(x_bounds$min_day, x_bounds$max_day)) +
    scale_y_discrete(drop = TRUE, limits = function(y) {
      # Retain Y categories with data within the X axis limits
      intersect(y, patient_drug_data %>%
           filter(StartDayRelativeToNearestHCT <= x_bounds$max_day & StopDayRelativeToNearestHCT >= x_bounds$min_day) %>%
           pull(Factor) %>% unique())
    }) +
    geom_segment(data = patient_drug_data %>%
                   mutate(
                     StartDayRelativeToNearestHCT = pmax(StartDayRelativeToNearestHCT, x_bounds$min_day),
                     StopDayRelativeToNearestHCT = pmin(StopDayRelativeToNearestHCT, x_bounds$max_day)
                   ) %>%
                   filter(StartDayRelativeToNearestHCT <= StopDayRelativeToNearestHCT),
                 aes(x = StartDayRelativeToNearestHCT, xend = StopDayRelativeToNearestHCT, 
                     y = Factor, yend = Factor), linewidth = 3.5, color = "gray75") +
    labs(
      title = paste("Drug Administration Timeline for Patient", patient_id),
      x = "Days Relative to Nearest HCT",
      y = "Drug",
      color = "Drug"
    ) +
    custom_theme
  
  return(plot)
}

