library(sf)
library(dplyr)
library(readr)

gpkg  <- "buchanan_contours.gpkg"
layer <- "contours"
csv   <- "edit_depths_template.csv"   # the CSV you edited

# 1) Read contours
cnt <- sf::st_read(gpkg, layer = layer, quiet = TRUE)

# 2) Ensure a stable row id (must match the order used when you exported the CSV)
if (!"id_row" %in% names(cnt)) {
  cnt$id_row <- seq_len(nrow(cnt))
}

template <- cnt %>% st_drop_geometry() %>%
  transmute(id_row, current_DEPTH_M = DEPTH_M, new_depth_ft = NA_real_)
# write_csv(template, "edit_depths_template.csv")
message("Open edit_depths_template.csv, fill new_depth_ft (in FEET), save, then run the UPDATE block below.")


# 3) Read the CSV and normalize names/types
depths <- readr::read_csv(csv, show_col_types = FALSE)
names(depths) <- tolower(gsub("[^a-z0-9]+","_", names(depths)))  # safe names

stopifnot(all(c("id_row","new_depth_ft") %in% names(depths)))
depths$id_row      <- as.integer(depths$id_row)
depths$new_depth_ft <- suppressWarnings(as.numeric(depths$new_depth_ft))

# 4) Make sure DEPTH_M exists & is numeric (do this OUTSIDE mutate)
if (!"DEPTH_M" %in% names(cnt)) cnt$DEPTH_M <- NA_real_
cnt$DEPTH_M <- suppressWarnings(as.numeric(cnt$DEPTH_M))

# 5) Join and update (feet → meters only where provided)
cnt <- cnt |>
  dplyr::left_join(dplyr::select(depths, id_row, new_depth_ft), by = "id_row") |>
  dplyr::mutate(DEPTH_M = ifelse(!is.na(new_depth_ft), new_depth_ft * 0.3048, DEPTH_M)) |>
  dplyr::select(-new_depth_ft)

# 6) Overwrite the layer with corrected depths
sf::st_write(cnt, gpkg, layer = layer, delete_layer = TRUE, quiet = TRUE)
cat("Updated DEPTH_M written back to", gpkg, "layer", layer, "\n")