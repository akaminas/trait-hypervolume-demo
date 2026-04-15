
`hypervolume_penguins.R`:
source("hypervolume_penguins.R")

```r
# Trait-based diversity with kernel hypervolumes
# Alexandros Kaminas

rm(list = ls())

library(hypervolume)
library(palmerpenguins)
library(dplyr)
library(ggplot2)

set.seed(42)

dir.create("figures", showWarnings = FALSE)

# --------------------------------------------------
# 1. Load and prepare data
# --------------------------------------------------

traits_raw <- penguins |>
  as_tibble() |>
  select(species, bill_length_mm, bill_depth_mm, flipper_length_mm) |>
  filter(
    !is.na(species),
    !is.na(bill_length_mm),
    !is.na(bill_depth_mm),
    !is.na(flipper_length_mm)
  ) |>
  filter(species %in% c("Adelie", "Gentoo"))

# Standardise trait axes before hypervolume estimation
traits_scaled <- traits_raw |>
  mutate(
    bill_length_z = as.numeric(scale(bill_length_mm)),
    bill_depth_z = as.numeric(scale(bill_depth_mm)),
    flipper_length_z = as.numeric(scale(flipper_length_mm))
  )

trait_cols <- c("bill_length_z", "bill_depth_z", "flipper_length_z")

adelie_traits <- traits_scaled |>
  filter(species == "Adelie") |>
  select(all_of(trait_cols))

gentoo_traits <- traits_scaled |>
  filter(species == "Gentoo") |>
  select(all_of(trait_cols))

# --------------------------------------------------
# 2. Estimate bandwidths
# --------------------------------------------------

bw_adelie <- estimate_bandwidth(adelie_traits, method = "silverman")
bw_gentoo <- estimate_bandwidth(gentoo_traits, method = "silverman")

cat("Bandwidths\n")
cat("Adelie:\n")
print(bw_adelie)
cat("\nGentoo:\n")
print(bw_gentoo)
cat("\n")

# --------------------------------------------------
# 3. Construct Gaussian kernel hypervolumes
# --------------------------------------------------
# samples.per.point is set below the default to keep the example lightweight.
# Increase it for more stable estimates in a more serious analysis.

hv_adelie <- hypervolume_gaussian(
  data = adelie_traits,
  name = "Adelie",
  kde.bandwidth = bw_adelie,
  samples.per.point = 200,
  quantile.requested = 0.95,
  quantile.requested.type = "probability",
  verbose = FALSE
)

hv_gentoo <- hypervolume_gaussian(
  data = gentoo_traits,
  name = "Gentoo",
  kde.bandwidth = bw_gentoo,
  samples.per.point = 200,
  quantile.requested = 0.95,
  quantile.requested.type = "probability",
  verbose = FALSE
)

cat("Hypervolume summaries\n")
print(hv_adelie)
print(hv_gentoo)
cat("\n")

# --------------------------------------------------
# 4. Pairwise set operations and overlap statistics
# --------------------------------------------------

hv_set <- hypervolume_set(
  hv1 = hv_adelie,
  hv2 = hv_gentoo,
  check.memory = FALSE,
  verbose = FALSE
)

overlap_stats <- hypervolume_overlap_statistics(hv_set)

cat("Overlap statistics\n")
print(overlap_stats)
cat("\n")

# --------------------------------------------------
# 5. Save a raw trait-space plot
# --------------------------------------------------

p_raw <- ggplot(
  traits_scaled,
  aes(x = bill_length_z, y = flipper_length_z, shape = species)
) +
  geom_point(size = 2.5, alpha = 0.8) +
  labs(
    title = "Trait space of two penguin species",
    subtitle = "Traits are z-standardised before hypervolume estimation",
    x = "Bill length (z-score)",
    y = "Flipper length (z-score)",
    shape = "Species"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  filename = "figures/trait_space_raw.png",
  plot = p_raw,
  width = 7,
  height = 5,
  dpi = 300
)

# --------------------------------------------------
# 6. Save a hypervolume pairs plot
# --------------------------------------------------

hv_list <- hypervolume_join(hv_adelie, hv_gentoo)

png(
  filename = "figures/hypervolume_pairs.png",
  width = 1600,
  height = 1600,
  res = 220
)

plot(
  hv_list,
  show.3d = FALSE,
  show.random = TRUE,
  show.density = FALSE,
  show.data = TRUE,
  show.contour = TRUE,
  contour.type = "kde",
  names = c("Bill length", "Bill depth", "Flipper length"),
  colors = c("#1f78b4", "#33a02c"),
  cex.random = 0.4,
  cex.data = 0.8,
  show.legend = TRUE,
  verbose = FALSE
)

dev.off()

# --------------------------------------------------
# 7. Export a concise summary table
# --------------------------------------------------

summary_tbl <- tibble(
  metric = names(overlap_stats),
  value = as.numeric(overlap_stats)
)

write.csv(summary_tbl, "figures/overlap_statistics.csv", row.names = FALSE)

cat("Saved files:\n")
cat("- figures/trait_space_raw.png\n")
cat("- figures/hypervolume_pairs.png\n")
cat("- figures/overlap_statistics.csv\n")
