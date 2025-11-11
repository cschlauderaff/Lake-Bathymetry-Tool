library(tmap); library(terra)

t_srs <- "EPSG:26915"
r_geo <- rast("buchanan_georef.tif")      # your georeferenced scan
cnt <- st_read(out_gpkg, layer = out_layer, quiet = TRUE)

# midpoints for labels
midpts <- st_line_sample(cnt, sample = 0.5) |> st_cast("POINT") |> st_as_sf()
midpts$id_row <- seq_len(nrow(cnt))

tmap_mode("plot")

# pick the right layer call for the raster (see #2 below)
if (terra::nlyr(r_geo) >= 3) {
  base <- tm_shape(r_geo) + tm_rgb(alpha = 0.6)
} else {
  base <- tm_shape(r_geo) + tm_raster(alpha = 0.6, palette = "Greys", legend.show = FALSE)
}

idx_map <- base +
  tm_shape(cnt) + tm_lines(lwd = 1) +
  tm_shape(midpts) + tm_text("id_row", size = 0.7, col = "black") +
  tm_layout(title = "Contour line index (id_row)")

tmap_save(idx_map, "buchanan_contours_index.png", dpi = 300, width = 8, height = 6, units = "in")