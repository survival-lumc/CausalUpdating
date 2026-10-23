library(tidyverse)

### ----------------------------------------------------------------------------
### ------------------------- Load data ----------------------------------------
### ----------------------------------------------------------------------------

# Load first scenario
newevents <- 25 
oldevents <- 25
load(paste("Results/summ", 
           newevents, "n", oldevents, "o.rda", sep = ""))
updating_periods <- cbind(periods, newevents, oldevents)
summ <- bind_cols(summm1, data.frame(oldevents))

# Load second scenario
newevents <- 25 
oldevents <- 55
load(paste("Results/summ", 
           newevents, "n", oldevents, "o.rda", sep = ""))
updating_periods <- rbind(updating_periods, cbind(periods, newevents, oldevents))
summ <- rbind(summ, bind_cols(summm1, data.frame(oldevents)))

# Load third scenario
newevents <- 25 
oldevents <- 85
load(paste("Results/summ", 
           newevents, "n", oldevents, "o.rda", sep = ""))
updating_periods <- rbind(updating_periods, cbind(periods, newevents, oldevents))
summ <- rbind(summ, bind_cols(summm1, data.frame(oldevents)))

# Load fourth scenario
newevents <- 25 
oldevents <- 115
load(paste("Results/summ", 
           newevents, "n", oldevents, "o.rda", sep = ""))
updating_periods <- rbind(updating_periods, cbind(periods, newevents, oldevents))
summ <- rbind(summ, bind_cols(summm1, data.frame(oldevents)))
summ <- add_column(summ, Method = rep(rownames(summm1), 4), .before = 1)

updating_periods$type <- case_when(
  updating_periods$type == "Factual update" ~ "Data interval for factual update",
  updating_periods$type == "Interventional update" ~ "Data interval for interventional update"
  )

### ----------------------------------------------------------------------------
### -------------------- Plot of the periods -----------------------------------
### ----------------------------------------------------------------------------

cols <- c("Data interval for interventional update" = "#BC8E7D", "Data interval for factual update" = "#2F3D70")
labels <- c(
  "25" = "Updating every 25 new events using 25 old events",
  "55" = "Updating every 25 new events using 55 old events",
  "85" = "Updating every 25 new events using 85 old events",
  "115" = "Updating every 25 new events using 115 old events"
)

gap_start <- as.Date("2020-05-05")
gap_end   <- as.Date("2020-07-27")
gap_length <- as.numeric(gap_end - gap_start)
gap_fraction <- 0.4

compress_date <- function(x) {
  as.numeric(x) -
    ifelse(
      x >= gap_end,
      gap_length * (1 - gap_fraction),
      ifelse(
        x >= gap_start,
        (as.numeric(x) - as.numeric(gap_start)) * (1 - gap_fraction),
        0
      )
    )
}

c_InDate <- compress_date(InDate)
c_OutDate <- compress_date(OutDate)
c_startdate <- compress_date(startdate)
c_enddate <- compress_date(enddate)

updating_periods <- updating_periods |>
  mutate(
    start_plot = compress_date(start),
    end_plot   = compress_date(end),
    cross_gap  = start < gap_start & end > gap_end
  ) |>
  group_by(newevents, oldevents) |>
  mutate(end_pred_plot = c(end_plot[-1], compress_date(enddate))) |>
  ungroup()

cross_gap <- updating_periods %>% filter(cross_gap) # Lines that cross the gap
no_gap <- updating_periods %>% filter(!cross_gap)   # Lines that don't cross the gap

ggplot(no_gap)  +
  geom_segment(
    aes(x = start_plot, xend = end_plot, y = period, yend = period, colour = type),
    linewidth = 0.8,
    lineend = "round"
  ) +
  geom_segment(    # Solid part before the gap
    data = cross_gap,
    aes(x = start_plot, xend = compress_date(gap_start), y = period, yend = period, colour = type),
    linewidth = 0.8,
    lineend = "round"
  ) +
  geom_segment(     # Dotted part representing June–August
    data = cross_gap,
    aes(x = compress_date(gap_start), xend = compress_date(gap_end), y = period, yend = period, colour = type),
    linewidth = 0.8,
    linetype = "11"
  ) +
  geom_segment(      # Solid part after the gap
    data = cross_gap,
    aes(x = compress_date(gap_end), xend = end_plot, y = period, yend = period, colour = type),
    linewidth = 0.8,
    lineend = "round"
  ) +
  geom_segment(
    aes(x = end_plot + 1, xend = end_pred_plot, y = period, yend = period, linetype = "Model deployment period"),
    linewidth = 0.5, color = "#555555", 
    arrow = arrow(length = unit(0.13, "cm"))
  ) +
  geom_point(aes(x = start_plot, y = period, colour = type), size = 2) +
  geom_point(aes(x = end_plot, y = period, colour = type), size = 2) +
  geom_point(data = cross_gap, aes(x = start_plot, y = period, colour = type), size = 2) +
  geom_point(data = cross_gap, aes(x = end_plot, y = period, colour = type), size = 2) +
  scale_colour_manual(
    name = "type",
    values = cols,
    aesthetics = c("colour", "fill")
  ) +
  scale_linetype_manual(
    name = NULL,
    values = c("Model deployment period" = "solid")
  )+
  geom_vline(
    xintercept = c(c_InDate, c_OutDate),
    linetype = "longdash",
    color = "darkred",
    linewidth = 0.5
  ) +
  geom_vline(
    xintercept = c(c_startdate, c_enddate),
    linetype = "solid",
    color = "#48275B",
    linewidth = 0.5
  ) +
  annotate("rect", 
           xmin = c_InDate, xmax = c_OutDate, ymin = -Inf, ymax = Inf, 
           alpha = 0.1, fill = "darkred") +
  annotate("rect", 
           xmin = c_startdate, xmax = c_enddate, ymin = -Inf, ymax = Inf, 
           alpha = 0.1, fill = "#48275B") +
  annotate("text", x = c_InDate + (c_OutDate - c_InDate)/2, y = 1, color = "darkred",
           label = "Treatment adoption \n phase", size = 3.5) +
  annotate("text", x = c_enddate - (c_enddate - c_startdate)/2, y = 6.6, color = "#48275B",
           label = "Validation data", size = 3.5) +
  scale_x_continuous(
    breaks = compress_date(as.Date(c("2020-05-01","2020-08-01", "2020-09-01", "2020-10-01",
                                     "2020-11-01", "2020-12-01","2021-01-01"))),
    labels = c("May\n2020", "Aug\n2020", "Sep\n2020", "Oct\n2020", "Nov\n2020",
               "Dec\n2020", "Jan\n2021")) +
  scale_y_discrete(expand = expansion(mult = c(0.15, 0.18))) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(size = 13, hjust = 0),
    plot.margin = margin(1, 12, 1, 5),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.margin = margin(0, 0, 0, 0),
    panel.spacing.x = unit(1.7, "lines"),
    axis.text = element_text(size = 9)
  ) +
  facet_wrap(~ oldevents, nrow = 2, labeller = as_labeller(labels))


### ----------------------------------------------------------------------------
### ------------------------ Tables --------------------------------------------
### ----------------------------------------------------------------------------
library(knitr)
library(kableExtra)

summ2 <- summ |> select(-c("Brier.score"))

# Table 1 for the main paper
# Part A
summ_wide <- summ2 %>%
  filter(oldevents %in% c("115")) %>%
  select(-c("oldevents")) %>%
  filter(!(Method %in% c("Interventional Refit (groups)", "Interventional Recalibration (groups)"))) %>% 
  mutate(across(c("Scaled.Brier", "Unos.c.index", "AUCt"), ~ { 
    # Extract the first numeric value 
    first_value <- as.numeric(sub(".*?(-?[0-9.]+).*", "\\1", .x))
    # Find the maximum within the column
    max_value <- max(first_value, na.rm = TRUE) 
    # Format the cell
    formatted <- paste0( "\\makecell{", gsub(" ", "\\\\\\\\ \\\\text{", .x), "}}" )
    # Bold if first value is the maximum 
    ifelse( first_value == max_value, paste0("\\textbf{", formatted, "}"), formatted ) 
    } ) ) %>%
  mutate(across(c("Calibration.intercept"), ~ { 
    # Extract the first numeric value 
    first_value <- as.numeric(sub(".*?(-?[0-9.]+).*", "\\1", .x))
    # Find the maximum within the column
    best_value <- min(abs(first_value), na.rm = TRUE) 
    # Format the cell
    formatted <- paste0( "\\makecell{", gsub(" ", "\\\\\\\\ \\\\text{", .x), "}}" )
    # Bold if first value is the maximum 
    ifelse( abs(first_value) == best_value, paste0("\\textbf{", formatted, "}"), formatted ) 
    } ) ) %>%
  mutate(across(c("OE.ratio", "Calibration.slope"), ~ { 
    # Extract the first numeric value 
    first_value <- as.numeric(sub(".*?(-?[0-9.]+).*", "\\1", .x))
    # Find the maximum within the column
    best_value <- min(abs(log(first_value)), na.rm = TRUE)
    # Format the cell
    formatted <- paste0( "\\makecell{", gsub(" ", "\\\\\\\\ \\\\text{", .x), "}}" )
    # Bold if first value is the maximum 
    ifelse( abs(log(first_value)) == best_value, paste0("\\textbf{", formatted, "}"), formatted ) 
    } ) ) %>%
  mutate(
    Method = ifelse(
      grepl(" ", Method),
      paste0("\\makecell{", gsub(" ", "\\\\\\\\", Method), "}"),
      Method
    )
  )

tab <- kable(
  summ_wide,
  format = "latex",
  booktabs = TRUE,
  escape = FALSE,
  col.names = c(
    " ", "\\makecell{Scaled \\\\ Brier}", "\\makecell{Uno's \\\\ c-index}", 
    "\\makecell{AUCt}", "\\makecell{OE \\\\ ratio}", "\\makecell{Calib. \\\\ intercept}",
    "\\makecell{Calib. \\\\ slope}"
  ),
  align = rep("c", 7)
) 

for (i in 1:4) {
  tab <- tab %>%
    row_spec(i, extra_latex_after = "\\addlinespace[6pt]")
}

tab



# Table 1 for the appendix paper
# Part A
summ_wide <- summ2 %>%
  filter(oldevents %in% c("25", "55")) %>%
  pivot_wider(
    id_cols = Method,
    names_from = oldevents,
    values_from = -c("Method", "oldevents"),
    names_glue = "{.value}_{oldevents}",
    names_vary = "slowest" # This gives the correct variables ordering
  ) %>% 
  mutate(across(matches("^Scaled.Brier|Unos.c.index|AUCt"), ~ { 
    # Extract the first numeric value 
    first_value <- as.numeric(sub(".*?(-?[0-9.]+).*", "\\1", .x))
    # Find the maximum within the column
    max_value <- max(first_value, na.rm = TRUE) 
    # Format the cell
    formatted <- paste0( "\\makecell{", gsub(" ", "\\\\\\\\ \\\\text{", .x), "}}" )
    # Bold if first value is the maximum 
    ifelse( first_value == max_value, paste0("\\textbf{", formatted, "}"), formatted ) 
  } ) ) %>%
  mutate(across(starts_with("Calibration.intercept"), ~ { 
    # Extract the first numeric value 
    first_value <- as.numeric(sub(".*?(-?[0-9.]+).*", "\\1", .x))
    # Find the maximum within the column
    best_value <- min(abs(first_value), na.rm = TRUE) 
    # Format the cell
    formatted <- paste0( "\\makecell{", gsub(" ", "\\\\\\\\ \\\\text{", .x), "}}" )
    # Bold if first value is the maximum 
    ifelse( abs(first_value) == best_value, paste0("\\textbf{", formatted, "}"), formatted ) 
  } ) ) %>%
  mutate(across(matches("^OE.ratio|Calibration.slope"), ~ { 
    # Extract the first numeric value 
    first_value <- as.numeric(sub(".*?(-?[0-9.]+).*", "\\1", .x))
    # Find the maximum within the column
    best_value <- min(abs(log(first_value)), na.rm = TRUE)
    # Format the cell
    formatted <- paste0( "\\makecell{", gsub(" ", "\\\\\\\\ \\\\text{", .x), "}}" )
    # Bold if first value is the maximum 
    ifelse( abs(log(first_value)) == best_value, paste0("\\textbf{", formatted, "}"), formatted ) 
  } ) ) %>%
  mutate(
    Method = ifelse(
      grepl(" ", Method),
      paste0("\\makecell{", gsub(" ", "\\\\\\\\", Method), "}"),
      Method
    )
  )

tab <- kable(
  summ_wide,
  format = "latex",
  booktabs = TRUE,
  escape = FALSE,
  col.names = c(
    " ",
    "\\makecell{Scaled \\\\ Brier}", "\\makecell{Uno's \\\\ c-index}", 
    "\\makecell{AUCt}", "\\makecell{OE \\\\ ratio}", "\\makecell{Calib. \\\\ intercept}",
    "\\makecell{Calib. \\\\ slope}",
    "\\makecell{Scaled \\\\ Brier}", "\\makecell{Uno's \\\\ c-index}", 
    "\\makecell{AUCt}", "\\makecell{OE \\\\ ratio}", "\\makecell{Calib. \\\\ intercept}",
    "\\makecell{Calib. \\\\ slope}"
  ),
  align = rep("c", 13)
) %>%
  add_header_above(
    c(
      " " = 1,
      "25 new events and 25 old events" = 6,
      "25 new events and 55 old events" = 6
    )
  )

for (i in 1:6) {
  tab <- tab %>%
    row_spec(i, extra_latex_after = "\\addlinespace[6pt]")
}

tab


# Part B
summ_wide <- summ2 %>%
  filter(oldevents %in% c("85", "115")) %>%
  pivot_wider(
    id_cols = Method,
    names_from = oldevents,
    values_from = -c("Method", "oldevents"),
    names_glue = "{.value}_{oldevents}",
    names_vary = "slowest" # This gives the correct variables ordering
  ) %>% 
  mutate(across(matches("^Scaled.Brier|Unos.c.index|AUCt"), ~ { 
    # Extract the first numeric value 
    first_value <- as.numeric(sub(".*?(-?[0-9.]+).*", "\\1", .x))
    # Find the maximum within the column
    max_value <- max(first_value, na.rm = TRUE) 
    # Format the cell
    formatted <- paste0( "\\makecell{", gsub(" ", "\\\\\\\\ \\\\text{", .x), "}}" )
    # Bold if first value is the maximum 
    ifelse( first_value == max_value, paste0("\\textbf{", formatted, "}"), formatted ) 
  } ) ) %>%
  mutate(across(starts_with("Calibration.intercept"), ~ { 
    # Extract the first numeric value 
    first_value <- as.numeric(sub(".*?(-?[0-9.]+).*", "\\1", .x))
    # Find the maximum within the column
    best_value <- min(abs(first_value), na.rm = TRUE) 
    # Format the cell
    formatted <- paste0( "\\makecell{", gsub(" ", "\\\\\\\\ \\\\text{", .x), "}}" )
    # Bold if first value is the maximum 
    ifelse( abs(first_value) == best_value, paste0("\\textbf{", formatted, "}"), formatted ) 
  } ) ) %>%
  mutate(across(matches("^OE.ratio|Calibration.slope"), ~ { 
    # Extract the first numeric value 
    first_value <- as.numeric(sub(".*?(-?[0-9.]+).*", "\\1", .x))
    # Find the maximum within the column
    best_value <- min(abs(log(first_value)), na.rm = TRUE)
    # Format the cell
    formatted <- paste0( "\\makecell{", gsub(" ", "\\\\\\\\ \\\\text{", .x), "}}" )
    # Bold if first value is the maximum 
    ifelse( abs(log(first_value)) == best_value, paste0("\\textbf{", formatted, "}"), formatted ) 
  } ) ) %>%
  mutate(
    Method = ifelse(
      grepl(" ", Method),
      paste0("\\makecell{", gsub(" ", "\\\\\\\\", Method), "}"),
      Method
    )
  )

tab <- kable(
  summ_wide,
  format = "latex",
  booktabs = TRUE,
  escape = FALSE,
  col.names = c(
    " ",
    "\\makecell{Scaled \\\\ Brier}", "\\makecell{Uno's \\\\ c-index}", 
    "\\makecell{AUCt}", "\\makecell{OE \\\\ ratio}", "\\makecell{Calib. \\\\ intercept}",
    "\\makecell{Calib. \\\\ slope}",
    "\\makecell{Scaled \\\\ Brier}", "\\makecell{Uno's \\\\ c-index}", 
    "\\makecell{AUCt}", "\\makecell{OE \\\\ ratio}", "\\makecell{Calib. \\\\ intercept}",
    "\\makecell{Calib. \\\\ slope}"
  ),
  align = rep("c", 13)
) %>%
  add_header_above(
    c(
      " " = 1,
      "25 new events and 85 old events" = 6,
      "25 new events and 115 old events" = 6
    )
  )

for (i in 1:6) {
  tab <- tab %>%
    row_spec(i, extra_latex_after = "\\addlinespace[6pt]")
}

tab
