# ============================================================
# Buchanan Lake – Digitize Depth Contours (R only, no QGIS)
# ============================================================
# What you get:
#   1) buchanan_georef.tif        – georeferenced scan (UTM 15N)
#   2) buchanan_contours.gpkg     – your traced contour lines (DEPTH_M column)
#   3) buchanan_contours_map.png  – pretty figure (scan + colorized lines)
#
# Required packages (install once):
# install.packages(c("sf","terra","tiff","leaflet","leafem","mapedit","tmap","dplyr","readr"))
# ------------------------------------------------------------

library(sf)
library(terra)
library(tiff)
library(leaflet)
library(leafem)
library(mapedit)
library(tmap)
library(dplyr)
library(readr)
library(raster)

# --------------------
# User settings
# --------------------
img_path <- "c2366010(1).tif"     # your scanned depth-map image
t_srs    <- "EPSG:26915"          # NAD83 / UTM zone 15N (meters)
out_tif  <- "buchanan_georef.tif"
out_gpkg <- "buchanan_contours.gpkg"
out_layer <- "contours"
out_png  <- "buchanan_contours_map.png"

# Optional: polygon to mask the scan to the lake footprint (makes the final map cleaner)
# Set to your Buchanan polygon if you have it; otherwise leave NULL to skip.
lake_poly_gpkg  <- NULL   # e.g. "DNR_Hydrography_OpenWater.gpkg"
lake_poly_layer <- NULL   # e.g. the layer name inside that gpkg
lake_dow        <- 56020900L

# Optional: perform a tight crop around the map frame BEFORE georeferencing
DO_MANUAL_CROP <- TRUE

# --------------------
# Helper: display TIFF upright at native resolution and collect clicks
# --------------------
click_points_on_scan <- function(img_file, prompt, min_pts = 4) {
  a <- tiff::readTIFF(img_file, native = FALSE, convert=TRUE)    # matrix/array; [0..1]
  if (length(dim(a)) == 3) {
    a <- 0.2989*a[,,1] + 0.5870*a[,,2] + 0.1140*a[,,3]  # grayscale
  }
  h <- nrow(a); w <- ncol(a)
  dev.new(width = 9, height = 9)
  par(xaxs = "i", yaxs = "i", mar = c(2,2,2,2))
  plot(NA, xlim = c(0, w), ylim = c(h, 0), xlab = "pixel", ylab = "line", asp = 1)
  rasterImage(a, 0, h, w, 0, interpolate = FALSE)       # crisp, upright (origin top-left)
  title(prompt)
  pts <- locator(type = "p")
  stopifnot(length(pts$x) >= min_pts)
  data.frame(x = pts$x, y = pts$y)
}

# --------------------
# 0) (Optional) Crop the scan to remove legend/margins
# --------------------
src_for_gcps <- img_path
if (DO_MANUAL_CROP) {
  message("Crop: Click top-left then bottom-right corners of the map frame (Esc to finish).")
  crop_clicks <- click_points_on_scan(img_path,
                                      prompt = "Crop box: click top-left, then bottom-right", min_pts = 2)
  if (nrow(crop_clicks) >= 2) {
    x1 <- min(crop_clicks$x); x2 <- max(crop_clicks$x)
    y1 <- min(crop_clicks$y); y2 <- max(crop_clicks$y)
    w  <- ceiling(x2 - x1);    h  <- ceiling(y2 - y1)
    # Use GDAL translate -srcwin to crop by pixel window
    cropped <- "buchanan_crop.tif"
    sf::gdal_utils("translate", img_path, cropped,
                   options = c("-srcwin", floor(x1), floor(y1), w, h))
    src_for_gcps <- cropped
    message(paste("Wrote:", cropped))
  }
}

# --------------------
# 1) Pick GCPs on the (cropped) scan (upright, full-res)
# --------------------
gcps_px <- click_points_on_scan(src_for_gcps,
                                prompt = "Click 6–12 GCPs on the SCAN (distinct features); Esc when done", min_pts = 4)

# --------------------
# 2) Pick the SAME points (same order) on a basemap (WGS84)
# --------------------
map_init <- leaflet() |> addTiles() |> setView(-95.5553, 46.4496, 13)
message("Now click the SAME features in the SAME order on the basemap. Click 'Done' when finished.")
gcp_ll <- mapedit::drawFeatures(map_init)  # sf POINTS; EPSG:4326
gcp_ll <- st_transform(gcp_ll, 4326)
stopifnot(nrow(gcp_ll) == nrow(gcps_px))
lonlat  <- st_coordinates(gcp_ll)

# Build -gcp vector for GDAL (pixel, line, lon, lat)
gcp_vec <- unlist(lapply(seq_len(nrow(gcps_px)), function(i) {
  c("-gcp", gcps_px$x[i], gcps_px$y[i], lonlat[i, "X"], lonlat[i, "Y"])
}))

# --------------------
# 3) Georeference with thin-plate spline & transparency
# --------------------
vrt_gcp <- tempfile(fileext = ".vrt")
sf::gdal_utils(
  util = "translate",
  source = src_for_gcps,                 # use your cropped file if you made one; else the original
  destination = vrt_gcp,
  options = c(gcp_vec, "-a_srs", "EPSG:4326", "-of", "VRT")
)

# 3b) Warp to the target CRS into another VRT (TPS + cubic; add alpha in VRT safely):
vrt_warp <- tempfile(fileext = ".vrt")
sf::gdal_utils(
  util = "warp",
  source = vrt_gcp,
  destination = vrt_warp,
  options = c("-t_srs", t_srs, "-tps", "-r", "cubic", "-dstalpha", "-overwrite", "-of", "VRT")
)

# 3c) Translate the warped VRT to GeoTIFF (optionally keep ALPHA band):
sf::gdal_utils(
  util = "translate",
  source = vrt_warp,
  destination = out_tif,
  options = c("-of", "GTiff", "-co", "COMPRESS=DEFLATE", "-co", "PREDICTOR=2", "-co", "ALPHA=YES")
)

message(paste("Wrote georeferenced raster:", out_tif))

# --------------------
# 4) Optional: mask the georef raster to the lake polygon (cleaner figure)
# --------------------
r_geo <- rast(out_tif)
if (!is.null(lake_poly_gpkg) && !is.null(lake_poly_layer) && file.exists(lake_poly_gpkg)) {
  lp <- st_read(lake_poly_gpkg, layer = lake_poly_layer, quiet = TRUE)
  if ("DOWLKNUM" %in% names(lp)) {
    lp <- lp |>
      mutate(DOWLKNUM = suppressWarnings(as.integer(DOWLKNUM))) |>
      filter(DOWLKNUM == lake_dow)
  }
  if (nrow(lp) > 0) {
    lp <- st_transform(lp, t_srs)
    r_geo <- mask(crop(r_geo, vect(lp)), vect(lp))
  }
}

# --------------------
# 5) Trace contours on TOP of the scan in Leaflet
#    (overlay the TIFF in WGS84 and draw polylines on it)
# --------------------
r_wgs84 <- terra::project(r_geo, "EPSG:4326")

# 2) Convert SpatRaster -> RasterLayer for leaflet::addRasterImage()
r_raster <- raster::raster(r_wgs84)

# 3) Palette domain from actual values
rng <- range(raster::values(r_raster), na.rm = TRUE)
pal <- leaflet::colorNumeric("Greys", domain = rng, na.color = NA)

# 4) Map centered on the raster’s extent
e <- raster::extent(r_raster)

mv <- leaflet() %>%
  addTiles() %>%
  # IMPORTANT: project = FALSE (it's already EPSG:4326)
  addRasterImage(r_raster, colors = pal, opacity = 0.75, project = FALSE) %>%
  fitBounds(lng1 = e@xmin, lat1 = e@ymin, lng2 = e@xmax, lat2 = e@ymax)

message("Trace contour lines ON the scan (Polyline tool). Double-click to finish each; click Done.")
cnt_sf <- mapedit::drawFeatures(mv)   # returns EPSG:4326
cnt_sf <- st_transform(cnt_sf, t_srs) # back to UTM 15N

# Add depth field (meters). Enter depths now, or fill later in a GIS/table editor.
cnt_sf$DEPTH_M <- NA_real_

# OPTIONAL: quick console entry of depths (feet) in drawing order, e.g. "5,10,15,20"
ans <- readline("Enter depths (feet) for each line, comma-separated (or press Enter to skip): ")
if (nzchar(ans)) {
  v <- as.numeric(strsplit(ans, ",")[[1]])
  if (length(v) == nrow(cnt_sf)) {
    cnt_sf$DEPTH_M <- v * 0.3048
  } else {
    message("Count mismatch; leaving DEPTH_M as NA. You can edit later.")
  }
}

# Save contours
st_write(cnt_sf, out_gpkg, layer = out_layer,
         delete_dsn = !file.exists(out_gpkg), quiet = TRUE)
message(paste("Saved contours to", out_gpkg, "layer:", out_layer))
