# ============================================================
# Buchanan Lake — Depth Heatmap from Digitized Contours
# ============================================================
# Outputs:
#   - buchanan_depth_heatmap.tif   (GeoTIFF raster, depth meters)
#   - buchanan_depth_heatmap.png   (PNG figure)
#
# install.packages(c("sf","terra","dplyr","gstat","tmap","viridisLite"))
library(sf)
library(terra)
library(dplyr)
library(gstat)
library(tmap)
library(viridisLite)

sf::sf_use_s2(FALSE)  # planar ops in UTM (avoids s2 quirks)

# ----------------------------
# User paths / params
# ----------------------------
hydro_gpkg     <- "water_dnr_hydrography_uncompressed.gpkg"  # official lake polygon
contours_gpkg  <- "buchanan_contours.gpkg"                   # your traced contours (DEPTH_M, LINESTRING)
contours_layer <- NULL  # set if your contours are not the first layer

crs_proj   <- "EPSG:26915"   # NAD83 / UTM zone 15N (meters)
dow_num    <- 56020900L      # Buchanan Lake
grid_res_m <- 10             # output raster resolution (m)
step_m     <- 25             # spacing for points along lines & shoreline (m)
idw_power  <- 2.5            # IDW exponent
do_smooth  <- TRUE           # small smoothing for nicer gradient

out_tif <- "buchanan_depth_heatmap.tif"
out_png <- "buchanan_depth_heatmap.png"

# ----------------------------
# 1) Read Buchanan polygon from hydrography (dowlknum)
# ----------------------------
stopifnot(file.exists(hydro_gpkg))
ow <- sf::st_read(hydro_gpkg, quiet = TRUE)

dow_field <- if ("dowlknum" %in% names(ow)) "dowlknum" else if ("DOWLKNUM" %in% names(ow)) "DOWLKNUM" else stop("No DOW field found.")
island_field <- if ("island" %in% names(ow)) "island" else if ("ISLAND" %in% names(ow)) "ISLAND" else NULL

ow <- ow |>
  mutate(DOW = suppressWarnings(as.integer(.data[[dow_field]])),
         ISL = if (!is.null(island_field)) suppressWarnings(as.integer(.data[[island_field]])) else NA_integer_) |>
  filter(DOW == dow_num, is.na(ISL) | ISL == 0)

stopifnot(nrow(ow) > 0)

lake_poly <- ow |>
  st_make_valid() |>
  st_collection_extract("POLYGON") |>
  st_union() |>
  st_transform(crs_proj)

# ----------------------------
# 2) Read your contours (DEPTH_M, LINESTRING)
# ----------------------------
stopifnot(file.exists(contours_gpkg))
if (is.null(contours_layer)) {
  contours_layer <- sf::st_layers(contours_gpkg)$name[1]
}
cnt <- sf::st_read(contours_gpkg, layer = contours_layer, quiet = TRUE)
stopifnot("DEPTH_M" %in% names(cnt))

cnt <- cnt |>
  st_transform(crs_proj) |>
  mutate(DEPTH_M = as.numeric(DEPTH_M)) |>
  filter(is.finite(DEPTH_M))

gt <- unique(as.character(st_geometry_type(cnt, by_geometry = TRUE)))
if (!any(gt %in% c("LINESTRING","MULTILINESTRING"))) {
  stop("Contours must be LINESTRING/MULTILINESTRING.")
}
# explode to LINESTRINGs (handles MULTILINESTRING cleanly)
cnt_ls <- st_cast(cnt, "LINESTRING", warn = FALSE)

# ----------------------------
# 3) Sample points along contours & along shoreline (depth=0)
# ----------------------------
# Helper: sample ~every step_m along each LINESTRING, carrying a depth value
set.seed(42)  # reproducible random sampling along lines

sample_points_from_lines <- function(lines_sf, depth_vec, step_m, crs) {
  g <- st_geometry(lines_sf)
  pts <- vector("list", length(g))
  n_out <- 0L
  for (i in seq_along(g)) {
    # random sampling breaks up regular artifacts along contours
    s <- st_line_sample(g[i], density = 1/step_m, type = "random")
    if (length(s) == 0) next
    p <- st_cast(s, "POINT")
    n <- length(p)
    n_out <- n_out + 1L
    pts[[n_out]] <- st_sf(DEPTH_M = rep(depth_vec[i], n), geometry = p)
  }
  if (n_out == 0L) return(st_sf(DEPTH_M = numeric(0), geometry = st_sfc(crs = crs)))
  out <- do.call(rbind, pts)
  st_set_crs(out, crs)
}

# Contour points (carry their depth)
pts_cnt <- sample_points_from_lines(cnt_ls, cnt_ls$DEPTH_M, step_m, st_crs(crs_proj))

# Shoreline points at 0 m (ensure LINESTRING input to st_line_sample)
shore_ml <- st_boundary(lake_poly)                          # MULTILINESTRING
shore_ls <- st_cast(shore_ml, "LINESTRING", warn = FALSE)   # explode to LINESTRINGs
# sample along each linestring and combine
shore_pts_list <- lapply(st_geometry(shore_ls), function(g) {
  s <- st_line_sample(g, density = 1/step_m)
  if (length(s) == 0) return(NULL)
  st_cast(s, "POINT")
})
shore_pts_geom <- do.call(c, shore_pts_list)
shore_pts <- st_sf(DEPTH_M = rep(0, length(shore_pts_geom)),
                   geometry = st_sfc(shore_pts_geom, crs = st_crs(crs_proj)))

# Combine & keep only points inside a small buffer of the lake
pts_all <- rbind(pts_cnt, shore_pts)
buf <- st_buffer(lake_poly, 5)  # 5 m
inside <- st_within(pts_all, buf, sparse = FALSE)[, 1]
pts_all <- pts_all[inside, ]
stopifnot(nrow(pts_all) > 0)

# ----------------------------
# 4) Interpolate to a raster (IDW), clip & (optionally) smooth
# ----------------------------
# Template with the lake's bounding box; no values yet (that's OK)
r_tmpl <- rast(ext(vect(lake_poly)), crs = crs_proj, resolution = grid_res_m)

# Optional: tighten to lake bbox, but still don't mask yet
r_tmpl <- crop(r_tmpl, vect(lake_poly))

# Grid cell centers for prediction
xy <- terra::xyFromCell(r_tmpl, 1:terra::ncell(r_tmpl))
grd_sf <- st_as_sf(data.frame(x = xy[,1], y = xy[,2]),
                   coords = c("x","y"), crs = crs_proj)

# IDW (gstat works directly with sf)
idw_res <- gstat::idw(DEPTH_M ~ 1,
                      locations = pts_all,
                      newdata   = grd_sf,
                      idp = idw_power)

pred <- as.numeric(idw_res$var1.pred)

# Put predictions back into the template raster
r_idw <- r_tmpl
values(r_idw) <- pred

# Now it's safe to post-process
r_idw <- clamp(r_idw, lower = 0)
r_idw <- mask(r_idw, vect(lake_poly))

# gentle de-speckle: two 5x5 passes
if (do_smooth) {
  w5 <- matrix(1, 5, 5) / 25
  r_idw <- focal(r_idw, w = w5, fun = sum, na.rm = TRUE)
  r_idw <- focal(r_idw, w = w5, fun = sum, na.rm = TRUE)
  r_idw <- mask(r_idw, vect(lake_poly))
}

writeRaster(r_idw, out_tif, overwrite = TRUE)

r_plot <- r_idw
r_plot_hr <- terra::rast(ext(r_plot), crs = crs(r_plot),
                         resolution = terra::res(r_plot) / 2)
r_plot <- terra::resample(r_plot, r_plot_hr, method = "bilinear")

# Legend breaks (clean ticks; adjust step if you prefer 1.0 m)
zmax  <- ceiling(terra::global(r_plot, "max", na.rm = TRUE)[[1]] * 2) / 2
brks  <- seq(0, zmax, by = 0.5)
pal   <- rev(viridisLite::viridis(length(brks) - 1))

tmap_mode("plot")

b0 <- sf::st_bbox(lake_poly)
dx <- as.numeric(b0["xmax"] - b0["xmin"])
dy <- as.numeric(b0["ymax"] - b0["ymin"])
pad_x <- 0.01 * dx          # ~2% on each side (tight → legend feels closer)
pad_y <- 0.12 * dy          # ~12% above/below for title & scale bar
b1 <- b0
b1["xmin"] <- b0["xmin"] - pad_x
b1["xmax"] <- b0["xmax"] + pad_x
b1["ymin"] <- b0["ymin"] - pad_y
b1["ymax"] <- b0["ymax"] + pad_y

map_heat_final <- 
  tm_shape(r_plot, bbox = b1) +
  tm_raster(breaks = brks, palette = pal, title = "Depth (m)",
            alpha = 0.98, legend.is.portrait = TRUE) +
  tm_shape(lake_poly) + tm_borders(lwd = 0.7, col = "#2c3e50") + 
  tm_layout(
    asp = 0,                                # use page area efficiently
    main.title = NULL,
    # Legend outside on the right (tall, readable)
    legend.outside = FALSE,
    legend.position = c(0.76, 0.85),
    legend.just = c("left", "center"),
    legend.title.size = 1.0,
    legend.text.size  = 0.9,
    legend.bg.color = "white",
    legend.frame = TRUE,
    # Compact margins (space for legend but no giant white half-page)
    outer.margins = c(0.02, 0.01, 0.03, 0.01),   # bottom, left, top, right
    inner.margins = c(0.01, 0.01, 0.02, 0.01),
    frame = FALSE
  ) +
  tm_credits("Buchanan Lake — Depth (m)",
             position = c(0.35, 0.985),
             just = c("center", "top"),
             size = 1.6, bg.color = NA) +
  tm_scale_bar(position = c(0.35, 0.06), text.size = 0.7)

tmap::tmap_save(
  map_heat_final,
  filename = "buchanan_depth_heatmap_layout.png",
  dpi = 600, width = 9, height = 6.8, units = "in"
)

cat("Wrote:\n - ", out_tif, "\n - ", out_png, "\n", sep = "")
