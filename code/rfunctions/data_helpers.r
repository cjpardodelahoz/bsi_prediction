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