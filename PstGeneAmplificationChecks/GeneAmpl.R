library(tidyverse)
library(viridis)

get_accession <- function(filename) {
  split_filename <- strsplit(filename, "_")[[1]]
  paste(split_filename[1:(length(split_filename) - 2)], collapse = "_")
}

files <- list.files(path = "./data", pattern = "*.txt", full.names = TRUE)

df_list <- lapply(files, function(file) {
  df <- read_csv(file, col_names = c("Gene", "PercAmb", "NAmbiguous", "NTotal", "Notes"))
  df$Accession <- get_accession(basename(file))
  df$PercAmb <- as.numeric(df$PercAmb)
  df$Coverage <- abs(df$PercAmb - 1)
  df
})

df_combined <- bind_rows(df_list)

df_combined <- df_combined |> filter(str_detect(Gene, "^PST"))
df_combined$binned_Coverage <- cut(df_combined$Coverage, breaks = seq(0,1,by=0.1))

df_aggregated <- df_combined |>
  filter(Coverage >= 0.10) |>
  group_by(NTotal, binned_Coverage) |>
  summarise(n=n(), Coverage=mean(Coverage), .groups="drop")

df_aggregated <- df_aggregated |> mutate(NTotal = as.integer(NTotal))

min_NTotal <- min(df_aggregated$NTotal)
max_NTotal <- max(df_aggregated$NTotal)
complete_NTotal <- tibble(NTotal = seq(min_NTotal, max_NTotal))

max_n <- 80
max_coverage <- 1
scaling_factor <- 100/80

barplt <- ggplot(df_aggregated, aes(x = NTotal, y = n)) +
  geom_bar(stat = "identity", aes(fill = forcats::fct_rev(binned_Coverage)), position = "stack", width = 3) +
  scale_fill_manual(values = hcl.colors(9, "Temps"), 
                    labels = c("[90-100]", "[80-90]", "[70-80]", "[60-70]", "[50-60]", "[40-50]", "[30-40]", "[20-30]", "[10-20]"),
                    guide = guide_legend(nrow = 1, label.position = "bottom", title.position = "top")) +
  stat_smooth(aes(y = Coverage * max_n, group = 1), size = 1, color = "plum4", fill = "plum") +
  labs(x = "Gene Length", y = "Number of Genes", fill = "Amplification Bin [%]") +
  theme_minimal() +
  scale_y_continuous(sec.axis = sec_axis(~ . * scaling_factor,
                                         name = "Average Amplification (%)",
                                         breaks = seq(0, 100, by = 20)),
                     limits = c(0, 80),
                     expand = c(0, 0)) +
  scale_x_continuous(limits = c(1000, 2210), breaks = seq(1000, 2400, by = 200), expand = c(0, 0)) +
  theme(legend.position = "bottom",
        legend.box = "horizontal",
        legend.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 8),
        legend.spacing.x = unit(0.5, 'cm'),
        axis.title.y.right = element_text(color = "plum4"), 
        axis.text.y.right = element_text(color = "plum4"))
barplt