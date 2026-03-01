library(sf)
library(dplyr)
library(stringr)
library(spdep)
library(mgcv)
library(tidyr)
library(ggplot2)
library(purrr)
library(scales)

## individual betas

df1 <- read.csv("C:/Users/A02369659/Documents/12-2020/Marine Project/chapter 3/Spatial and WTP/ZIPcodesdistances/ols_travel2.csv", stringsAsFactors = FALSE)
colnames(df1)
nrow(df1)


df<- read.csv("C:/Users/A02369659/Documents/12-2020/Marine Project/New Data/beta_i2.csv", stringsAsFactors = FALSE)
colnames(df)
nrow(df)
4452/12

# Keep only the first occurrence of each unique "pid" (not for 441 441)
df <- df %>%
  distinct(pid, .keep_all = TRUE)
nrow(df)
nrow(df1)

df  <- merge(df1, df[c("cluc.1","cluc.2","cluc.3","cluc.4","cluc.5","pid")], by="pid", all.x=TRUE)
colnames(df)
nrow(df)
getwd()
write.csv(df,"ols_b1_2.csv")


library(sf)
library(dplyr)
library(stringr)

#### Setup + helper to standardize ZIPs
standardize_zip5 <- function(df) {
  nm  <- names(df)
  lnm <- tolower(nm)
  candidates <- c("zip5","zip","zip_code","zcta5","zcta","zcta5ce10","name")
  hit <- candidates[candidates %in% lnm][1]
  if (is.na(hit)) stop("No ZIP-like column found. Available names: ", paste(nm, collapse=", "), call.=FALSE)
  
  col <- nm[match(hit, lnm)]         # original-cased column name
  z   <- as.character(df[[col]])
  
  # Try to extract an exact 5-digit ZIP; if missing, take first run of digits and pad
  zip5 <- str_extract(z, "\\b\\d{5}\\b")
  miss <- is.na(zip5)
  if (any(miss)) {
    tmp <- str_extract(z[miss], "\\d+")
    zip5[miss] <- tmp
  }
  zip5 <- str_pad(substr(zip5, 1, 5), 5, pad = "0")
  
  df$Zip5 <- zip5
  df
}

# 2) Read polygons + make valid + standardize Zip
zip_poly <- st_read("C:/Users/A02369659/Documents/12-2020/Marine Project/chapter 3/Spatial and WTP/ZIPcodesdistances/OR_zip/OR_zipCodes_m.shp", quiet = TRUE) |>
  st_make_valid() |> # fixes invalid geometries (common in shapefiles).
  standardize_zip5() # adds the clean Zip5 to the polygon layer.
head(zip_poly)
nrow(zip_poly) # 429 ZIP polygons

#install.packages("janitor")      # run once
library(janitor)

# 3) Read the table + standardize Zip 441x371
people <- read.csv("C:/Users/A02369659/Documents/12-2020/Marine Project/chapter 3/Spatial and WTP/ZIPcodesdistances/ols_b1_2.csv", stringsAsFactors = FALSE) |>
  clean_names() |>
  standardize_zip5()
head(people)
people$zip
people$zip5
nrow(people)
colnames(people)


keep <- c("Zip5","pid","cluc_1","cluc_2","cluc_3","cluc_4","cluc_5", "total_travel_time","total_kilometers", "mr1","mr2","mr3","mr4","mr5","age35","bach","inc75","latino1","gender1","emp1","urban")

people_keep <- people |>
  select(any_of(keep)) |> #keep only that list if it exists.
  rename_with(~paste0("p_", .), -Zip5)# add p_ in front of all columns except Zip5
anyDuplicated(names(people_keep))
colnames(people_keep)
names(people_keep)
nrow(people_keep)
head(people_keep)

# Join people → polygons
poly_people <- left_join(zip_poly, people_keep, by = "Zip5")
class(poly_people)
colnames(poly_people)
nrow(poly_people) #A left_join from polygons (429 rows) to people (371 rows) duplicates polygons for ZIPs that match multiple people.
glimpse(poly_people)
#plot(poly_people)

###### here is the true map!!
library(sf); library(dplyr); library(ggplot2)

# 1) Put ZIP polygons in a metric CRS
zip_m   <- st_transform(zip_poly, 5070)      # put polygons in a metric CRS
colnames(zip_m)
zip_pts <- st_point_on_surface(zip_m["Zip5"])  # one point per ZIP. makes one point per polygon (so 429 points, one for each ZIP, placed inside the polygon). Now have a geometry for each ZIP. Your shapefile has some ZIPs that are multipart (same ZIP code broken into several polygons). When sf makes one point “on the surface” of each feature, it has to decide which attributes go with that point. It’s warning: “I’m assuming all parts of this feature have the same attributes.”

# 2) Match each person to their ZIP point (NO dplyr join on sf)
idx <- match(as.character(people_keep$Zip5), zip_pts$Zip5) # Take each person’s ZIP (people_keep$Zip5) Find which ZIP point it corresponds to (match(...))
if (anyNA(idx)) {
  cat(sum(is.na(idx)), "respondents have ZIPs missing from the shapefile:\n")
  print(setdiff(people_keep$Zip5, zip_pts$Zip5))
}

resp_sf <- st_sf(people_keep, geometry = st_geometry(zip_pts)[idx], crs = st_crs(zip_pts)) #Pull that point geometry (st_geometry(zip_pts)[idx]) Build an sf object with the people data and the ZIP point geometry
nrow(resp_sf)  # should be 371 or 441 for 441x441. each one now spatial (a POINT). So avoided the polygon duplication mess.
glimpse(resp_sf)

##gpt
# --- packages
library(sf)
library(dplyr)
library(ggplot2)
library(spdep)
library(mgcv)
library(gstat)
library(automap)
library(purrr)
library(tidyr)
library(scales)





# ============================================================
# 0) Prep data (points) and choose the 5 attribute columns
# ============================================================
resp_m <- st_transform(resp_sf, 5070)

numify <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "NaN", ".", "Inf", "-Inf")] <- NA
  suppressWarnings(as.numeric(x))
}

coef_cols <- grep("^p_cluc[_.]\\d+$", names(resp_m), value = TRUE)
stopifnot(length(coef_cols) == 5)

resp_m <- resp_m %>%
  mutate(across(all_of(coef_cols), numify))

# Name them as in the paper (Figure 5 / Table 8 / Table 9)
# IMPORTANT: if your order differs, fix it here
attr_names <- c("SQ", "MRS", "Job", "Carbon", "Cost")
names_map  <- setNames(attr_names, coef_cols)

# ============================================================
# 1) Table 8: Global Moran’s I (analytical + permutation)
# ============================================================
set.seed(12345)
resp_w <- resp_m %>% st_jitter(amount = 10)
coords <- st_coordinates(resp_w)

make_knn_connected <- function(coords, k_start = 6, k_max = 30) {
  k <- k_start
  repeat {
    nb <- knn2nb(knearneigh(coords, k = k), sym = TRUE)
    if (n.comp.nb(nb)$nc == 1 || k >= k_max) break
    k <- k + 1
  }
  lw <- nb2listw(nb, style = "W")
  attr(lw, "k_used") <- k
  lw
}

lw_full <- make_knn_connected(coords, k_start = 6)

moran_tbl_analytical <- function(sfpoints, cols, listw, nm_map) {
  map_dfr(cols, function(v) {
    x  <- sfpoints[[v]]
    ok <- is.finite(x)

    lw_sub <- spdep::subset.listw(listw, subset = ok,
                                  zero.policy = TRUE, adjust.n = TRUE)

    mt <- spdep::moran.test(x[ok], lw_sub,
                            randomisation = TRUE, zero.policy = TRUE)

    I   <- as.numeric(mt$estimate[["Moran I statistic"]])
    E_I <- as.numeric(mt$estimate[["Expectation"]])
    V_I <- as.numeric(mt$estimate[["Variance"]])

    tibble(
      attribute = nm_map[[v]],
      I      = I,
      E_I    = E_I,
      Var_I  = V_I,
      z      = (I - E_I)/sqrt(V_I),
      p_value = mt$p.value,
      k_used  = attr(listw, "k_used"),
      n_used  = sum(ok)
    )
  })
}

moran_tbl_perm <- function(sfpoints, cols, listw, nm_map, nsim = 999, base_seed = 4242) {
  map2_dfr(cols, seq_along(cols), function(v, i) {
    x  <- sfpoints[[v]]
    ok <- is.finite(x)

    lw_sub <- spdep::subset.listw(listw, subset = ok,
                                  zero.policy = TRUE, adjust.n = TRUE)

    set.seed(base_seed + i)
    mt <- spdep::moran.mc(x[ok], lw_sub, nsim = nsim,
                          zero.policy = TRUE, alternative = "greater")

    tibble(
      attribute = nm_map[[v]],
      I_perm = as.numeric(mt$statistic),
      p_perm = mt$p.value
    )
  })
}

tab8_ana  <- moran_tbl_analytical(resp_w, coef_cols, lw_full, names_map)
tab8_perm <- moran_tbl_perm(resp_w, coef_cols, lw_full, names_map, nsim = 999)

Table8 <- left_join(tab8_ana, tab8_perm, by = "attribute") %>%
  select(attribute, I, E_I, Var_I, z, p_value, I_perm, p_perm, n_used, k_used)

Table8
# write.csv(Table8, "Table8_MoranI_attributes.csv", row.names = FALSE)

# ============================================================
# 2) Table 9: Spatial GAM diagnostics + incremental fit vs null
# ============================================================
xy  <- st_coordinates(resp_m)
dat <- resp_m %>% st_drop_geometry() %>% mutate(x = xy[,1], y = xy[,2])

fit_spatial <- function(col, k = 80) {
  gam(as.formula(paste0(col, " ~ s(x, y, k = ", k, ")")),
      data = dat, family = gaussian(), method = "REML")
}

fit_null <- function(col) {
  gam(as.formula(paste0(col, " ~ 1")),
      data = dat, family = gaussian(), method = "REML")
}

mods      <- setNames(lapply(coef_cols, fit_spatial, k = 80), attr_names)
mods_null <- setNames(lapply(coef_cols, fit_null), attr_names)

extract_table9 <- function(m, m0, nm) {
  sm <- summary(m)
  st <- sm$s.table[1, , drop = FALSE]   # row for s(x,y)

  at <- anova(m0, m, test = "F")

  tibble(
    attribute = nm,
    edf     = as.numeric(st[,"edf"]),
    F       = as.numeric(st[,"F"]),
    p       = as.numeric(st[,"p-value"]),
    devexpl = sm$dev.expl,
    r2_adj  = sm$r.sq,
    dfdiff  = as.numeric(at$Df[2]),
    devdiff = as.numeric(at$Deviance[2]),
    Finc    = as.numeric(at$F[2]),
    pinc    = as.numeric(at$`Pr(>F)`[2])
  )
}

Table9 <- map2_dfr(names(mods), seq_along(mods), function(nm, i){
  extract_table9(mods[[nm]], mods_null[[nm]], nm)
})

Table9
# write.csv(Table9, "Table9_GAM_diagnostics_attributes.csv", row.names = FALSE)

# ============================================================
# 3) Figure 5: Smooth spatial variation maps for the 5 attributes
# ============================================================
state_m <- st_union(st_transform(zip_poly, 5070))

grid_poly <- st_make_grid(state_m, cellsize = 10000, what = "polygons") %>%
  st_as_sf() %>%
  st_intersection(state_m)

grid_cent <- st_centroid(grid_poly)
gc        <- st_coordinates(grid_cent)
newdat    <- data.frame(x = gc[,1], y = gc[,2])

pred_mat <- sapply(mods, function(m) predict(m, newdata = newdat, type = "response"))
pred_df  <- as.data.frame(pred_mat)
pred_sf  <- cbind(grid_poly, pred_df)

# If you want ONE common scale across facets, standardize to z-scores:
pred_z <- as.data.frame(scale(pred_df))
names(pred_z) <- paste0("z_", names(pred_z))
pred_sf_z <- cbind(grid_poly, pred_z)

surf <- pred_sf_z %>%
  pivot_longer(cols = starts_with("z_"),
               names_to = "attribute",
               values_to = "coef_z") %>%
  mutate(attribute = sub("^z_", "", attribute))

surf$attribute <- factor(surf$attribute, levels = c("SQ","MRS","Job","Carbon","Cost"))

L <- max(abs(surf$coef_z), na.rm = TRUE)

ggplot(surf) +
  geom_sf(aes(fill = coef_z), color = NA) +
  geom_sf(data = st_transform(zip_poly, 5070),
          fill = NA, color = "grey60", linewidth = 0.15) +
  scale_fill_gradient2(low = "white", mid = "white", high = "firebrick",
                       midpoint = 0, limits = c(-L, L),
                       oob = scales::squish, name = "Coefficient") +
  coord_sf() +
  theme_void(base_size = 12) +
  theme(legend.position = "right",
        strip.text = element_text(face = "bold")) +
  facet_wrap(~ attribute, ncol = 2)
