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
            select(-c("gut metagenome", "metagenome", "mouse gut metagenome", "<not present>")) %>%
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
