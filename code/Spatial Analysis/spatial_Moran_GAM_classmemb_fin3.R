## results generate slighlty different local morans because distincts jitters random every time, mran.mc montecarlo and because of jitter also W. 

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
#head(zip_poly)
#nrow(zip_poly) # 429 ZIP polygons

#install.packages("janitor")      # run once
library(janitor)

# 3) Read the table + standardize Zip 441x371
people <- read.csv("C:/Users/A02369659/Documents/12-2020/Marine Project/chapter 3/Spatial and WTP/ZIPcodesdistances/ols_travel2_covlcm.csv", stringsAsFactors = FALSE) |>
clean_names() |>
standardize_zip5()


keep <- c("Zip5","pid","cluc_1","cluc_2","cluc_3","cluc_4","cluc_5", "total_travel_time","total_kilometers", "mr1","mr2","mr3","mr4","mr5","age35","bach","inc75","latino1","gender1","emp1","urban")

people_keep <- people |>
  select(any_of(keep)) |> #keep only that list if it exists.
  rename_with(~paste0("p_", .), -Zip5)# add p_ in front of all columns except Zip5


# Join people table → polygons
poly_people <- left_join(zip_poly, people_keep, by = "Zip5")
class(poly_people)
colnames(poly_people)
nrow(poly_people) #A left_join from polygons (429 rows) to people (371 rows) duplicates polygons for ZIPs that match multiple people.

###### Spatial ANalysis
library(sf); library(dplyr); library(ggplot2)

# 1) Put ZIP polygons in a metric CRS
zip_m   <- st_transform(zip_poly, 5070)      # put polygons in a metric CRS
colnames(zip_m)
zip_pts <- st_point_on_surface(zip_m["Zip5"])  # one point per ZIP. makes one point per polygon (so 429 points, one for each ZIP, placed inside the polygon). Now have a geometry for each ZIP. Shapefile has some ZIPs that are multipart (same ZIP code broken into several polygons). When sf makes one point “on the surface” of each feature, it has to decide which attributes go with that point. It’s warning: “I’m assuming all parts of this feature have the same attributes.”

# 2) Match each person to their ZIP point (NO dplyr join on sf)
idx <- match(as.character(people_keep$Zip5), zip_pts$Zip5) # Take each person’s ZIP (people_keep$Zip5) Find which ZIP point it corresponds to (match(...))
if (anyNA(idx)) {
  cat(sum(is.na(idx)), "respondents have ZIPs missing from the shapefile:\n")
  print(setdiff(people_keep$Zip5, zip_pts$Zip5))
}

resp_sf <- st_sf(people_keep, geometry = st_geometry(zip_pts)[idx], crs = st_crs(zip_pts)) #Pull that point geometry (st_geometry(zip_pts)[idx]) Build an sf object with the people data and the ZIP point geometry
#nrow(resp_sf)  # should be 371. each one now spatial (a POINT). So avoided the polygon duplication mess.

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

# --- 0) Data & CRS (meters) - Make spatial weights and run Moran’s I

resp_m  <- st_transform(resp_sf,   5070)
colnames(resp_m)

#########

# grab the class-prob columns (works with "clu_1" or "clu.1")
prob_cols <- grep("^p_cluc[_.]\\d+$", names(resp_m), value = TRUE) # rom all column names, give me the ones that start with p_cluc_ or p_cluc. and end with digits
stopifnot(length(prob_cols) > 0)
resp_m <- resp_m %>%
  mutate(across(all_of(prob_cols), ~ as.numeric(trimws(.x))))  # blanks -> NA are fine. forced them to numeric.

#######

# --- 1) Moran's I table (Table 8 style) -------------------------------
# inverse-square, row-standardized spatial weights (k-NN → radius → 1/d^2)

###### 1st weights with 3 sub-grohs
resp_w  <- st_transform(resp_m, 5070) %>% st_jitter(amount = 10) # move points a tiny bit (10 m) so that points that were exactly on top of each other don’t break k-NN.
coords  <- st_coordinates(resp_w) # extract x,y matrix to feed to spdep.

# Then the helper:
make_knn_connected <- function(coords, k_start = 4, k_max = 30) {
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

# try k-nearest neighbors with k = 4; check if the neighbor graph is fully connected (n.comp.nb(...) == 1); if not, increase k → 5 → 6 → … until it becomes one connected graph or you hit k_max; convert to a spatial weights list (nb2listw); tag the object with the k actually used.

lw_full <- make_knn_connected(coords, k_start = 6)
#started at 6, and it told tells that “neighbor object has 3 sub-graphs”, so it kept increasing until it found a k that connects everything (it ended with k = 8 in the Moran’s output). knn not fully connected when k =6. 3 separate islands (sub-graphs). But the helper keeps increasing N and fixed this with knn=8.

# Moran’s helper:
moran_tbl <- function(sfpoints, prob_cols, listw, nsim = 999) {
  map_dfr(prob_cols, function(v){
    x  <- sfpoints[[v]]
    ok <- is.finite(x)
    if (!any(ok)) return(tibble(class = sub("^p_cluc[_.]","class ", v),
                                I = NA_real_, p_perm = NA_real_, n_used = 0L,
                                k_used = attr(listw, "k_used")))
    lw_sub <- spdep::subset.listw(listw, subset = ok, zero.policy = TRUE, adjust.n = TRUE)
    mt <- moran.mc(x[ok], lw_sub, nsim = nsim, zero.policy = TRUE, alternative = "greater")
    tibble(
      class  = sub("^p_cluc[_.]","class ", v),
      I      = as.numeric(mt$statistic),
      p_perm = mt$p.value,
      n_used = sum(ok),
      k_used = attr(listw, "k_used")
    )
  })
}

# loop over each class-prob column; drop rows where that column is NA; subset the weights to those non-NA rows (this is the subtle bit most people forget!); run a permutation Moran’s I (moran.mc) with 999 sims; return a tibble with Moran’s I, p-value, how many obs were used, and k.

#Global Moran's I → one number for the whole study area → “is there spatial autocorrelation overall?”
mi_tbl <- moran_tbl(resp_w, prob_cols, lw_full, nsim = 999)
print(mi_tbl)

# I ran the Monte Carlo (permutation) version of global Moran’s I, once per variable (p_cluc_1 … p_cluc_5). That’s why I got a tibble with 5 rows, not 366 row.

### class 1: positive spatial autocorrelation, significant (p = .001) class 2,5: small but significant classes 3–4: basically no spatial pattern.

# Analytical Moran (not permutation)
library(spdep)
library(dplyr)
library(purrr)
library(spdep)

moran_tbl_analytical <- function(sfpoints, prob_cols, listw) {
  map_dfr(prob_cols, function(v) {
    x  <- sfpoints[[v]]
    ok <- is.finite(x)
    
    if (!any(ok)) {
      return(tibble(
        class   = sub("^p_cluc[_.]", "class ", v),
        I       = NA_real_,
        E_I     = NA_real_,
        Var_I   = NA_real_,
        z       = NA_real_,
        p_value = NA_real_,
        k_used  = attr(listw, "k_used")
      ))
    }
    
    # match weights to the non-missing rows
    lw_sub <- spdep::subset.listw(listw, subset = ok,
                                  zero.policy = TRUE, adjust.n = TRUE)
    
    mt <- spdep::moran.test(x[ok], lw_sub,
                            randomisation = TRUE,  # classic null
                            zero.policy   = TRUE)
    
    tibble(
      class   = sub("^p_cluc[_.]", "class ", v),
      I       = as.numeric(mt$estimate[["Moran I statistic"]]),
      E_I     = as.numeric(mt$estimate[["Expectation"]]),
      Var_I   = as.numeric(mt$estimate[["Variance"]]),
      z       = as.numeric((mt$estimate[["Moran I statistic"]] -
                              mt$estimate[["Expectation"]]) /
                             sqrt(mt$estimate[["Variance"]])),
      p_value = mt$p.value,
      k_used  = attr(listw, "k_used")
    )
  })
}


# I already have these:
# resp_w  <- ...  (your jittered points)
# prob_cols <- grep("^p_cluc[_.]\\d+$", names(resp_w), value = TRUE)
# lw_full <- make_knn_connected(...)

mi_ana_tbl <- moran_tbl_analytical(resp_w, prob_cols, lw_full)
mi_ana_tbl

moran_tbl_perm <- function(sfpoints, prob_cols, listw, nsim = 999) {
  map_dfr(prob_cols, function(v){
    x  <- sfpoints[[v]]
    ok <- is.finite(x)
    if (!any(ok)) {
      return(tibble(
        class   = sub("^p_cluc[_.]","class ", v),
        I_perm  = NA_real_,
        p_perm  = NA_real_,
        n_used  = 0L,
        k_used  = attr(listw, "k_used")
      ))
    }
    lw_sub <- spdep::subset.listw(listw, subset = ok, zero.policy = TRUE, adjust.n = TRUE)
    mt <- spdep::moran.mc(x[ok], lw_sub, nsim = nsim,
                          zero.policy = TRUE, alternative = "greater")
    tibble(
      class   = sub("^p_cluc[_.]","class ", v),
      I_perm  = as.numeric(mt$statistic),
      p_perm  = mt$p.value,
      n_used  = sum(ok),
      k_used  = attr(listw, "k_used")
    )
  })
}

ana  <- moran_tbl_analytical(resp_w, prob_cols, lw_full)
perm <- moran_tbl_perm(resp_w, prob_cols, lw_full, nsim = 999)

mi_all <- left_join(ana, perm, by = c("class","k_used"))
mi_all

## Moran's I until here → explain that results can slightly be different in the preamble. 

##### MAPS #####

####################
## ============================================================
## PACKAGES
## ============================================================
library(sf)
library(dplyr)
library(ggplot2)
library(mgcv)
library(tidyr)
library(scales)

## ============================================================
## 0) DATA TO METRIC CRS
## ============================================================
state_m <- st_union(st_transform(zip_poly, 5070)) #dissolves all ZIP polygons into one state outline polygon (a single geometry). clip the grid.
resp_m  <- st_transform(resp_sf, 5070) #reprojects ZIP polygons into EPSG:5070 (meters)

## ============================================================
## 1) BUILD WORKING DATA FRAME
## ============================================================
xy  <- st_coordinates(resp_m) #Extracts an n x 2 matrix of point coordinates: columns are x and y (in meters).
dat <- resp_m |>
  st_drop_geometry() |> #turns the sf object into a plain data frame (keeps attributes, drops geometry).
  mutate(x = xy[,1], y = xy[,2]) #adds numeric coordinate columns so GAM can model s(x,y)

prob_cols <- grep("^p_cluc[_.]\\d+$", names(dat), value = TRUE) #select outcome variables
stopifnot(length(prob_cols) > 0) # aborts if none are found (so later code doesn’t silently fail).

dat <- dat |>
  mutate(across(all_of(prob_cols),
                ~ suppressWarnings(as.numeric(trimws(.x))))) #celan data

## ============================================================
## 2) PREDICTION GRID (10 km)
## ============================================================
grid_poly <- st_make_grid(state_m, cellsize = 10000, what = "polygons") |> #makes a square grid covering state_m’s bounding box, with 10,000 m = 10 km cells.
  st_as_sf() |> # converts grid to an sf object.
  st_intersection(state_m) # st_intersection(state_m)

grid_cent <- st_centroid(grid_poly) # computes the centroid of each grid cell polygon.
gc        <- st_coordinates(grid_cent) #extracts centroid x/y.
newdat    <- data.frame(x = gc[,1], y = gc[,2]) #the prediction dataset (one row per grid cell), with x/y only.

## ============================================================
## 3) FIT EACH CLASS (UNORDERED GAMs) - Because fit separate GAMs independently for each class probability, rather than a joint multinomial model.
## ============================================================

#fit the GAMS
fit_one <- function(col, k = 80) { #Defines a function that fits one surface for one class-prob column (e.g., p_cluc_1)
  fml <- as.formula(paste0(col, " ~ s(x, y, k = ", k, ")")) #join in one string. Builds a formula like: p_cluc_1 ~ s(x, y, k = 80) s(x,y) is a 2D smooth over space. k=80 is the basis dimension (upper bound on wiggliness / complexity).
  m <- gam(fml, data = dat, family = quasibinomial(), method = "REML") # fits the model for that class.family = quasibinomial():uses a logit-type mean function but allows overdispersion (variance not forced to be exactly binomial). works when your “probabilities” aren’t literally 0/1 Bernoulli outcomes.chooses the smoothing penalty via REML
  p <- predict(m, newdata = newdat, type = "response") # Predicts at the grid centroids. type="response" returns predictions on the probability scale (0–1), not the linear predictor.
  pmin(pmax(p, 0), 1) #Clips any small numerical overshoots below 0 or above 1 back into [0,1]
}

pred_list <- lapply(prob_cols, fit_one) # Fits one GAM per class probability column.Output: list of predicted vectors, each length = number of grid cells.
pred_mat  <- do.call(cbind, pred_list)   # cbind them into a matrix: rows = grid cells cols = classes (e.g., 5 columns)

## ============================================================
## 4) RENORMALIZE TO SIMPLEX
## ============================================================
row_sums <- rowSums(pred_mat) #For each grid cell, sum predicted probabilities across classes.
row_sums[row_sums == 0] <- 1 #If a row sum is 0 (rare but possible numerically), set to 1 to avoid division-by-zero.
pred_norm <- pred_mat / row_sums #Divide each row by its sum → now class probs at each grid cell sum to 1.

colnames(pred_norm) <- sub("^p_cluc[_.]", "class_", prob_cols) #Renames columns from p_cluc_1 → class_1.
pi_df <- as.data.frame(pred_norm) # Converts matrix to data frame.

## ============================================================
## 5) PROBABILITY MAPS (optional, same as before)
## ============================================================
surf <- cbind(grid_poly, pi_df) |> #attaches predicted probs to each grid cell polygon.
  pivot_longer(starts_with("class_"), #
               names_to = "class",
               values_to = "prob")

pal <- c("#FFF200","#FFB300","#FF8F00","#E65100","#B71C1C")

ggplot() +
  geom_sf(data = surf, aes(fill = prob), color = NA) +
  geom_sf(data = st_transform(zip_poly, 5070),
          fill = NA, color = "grey60", linewidth = 0.15) +
  scale_fill_gradientn(colours = pal,
                       limits = c(0,1),
                       breaks = seq(0,1,by=.125),
                       labels = percent) +
  coord_sf() +
  theme_minimal(base_size = 12) +
  theme(axis.text = element_blank(),
        panel.grid = element_blank(),
        legend.title = element_text(size = 10)) +
  labs(title = "Unordered GAMs, renormalized to sum to 1",
       fill  = "Pr(class)") +
  facet_wrap(~ class, ncol = 2)

## ============================================================
## 6) DOMINANT MAP (style like your screenshot)
## ============================================================
# give id before pivot
grid_poly_id <- grid_poly |>
  mutate(cell_id = row_number())

surf_norm <- cbind(grid_poly_id, pi_df) |>
  pivot_longer(
    starts_with("class_"),
    names_to  = "class",
    values_to = "prob"
  )

dominant_norm <- surf_norm |>
  group_by(cell_id) |>
  slice_max(prob, n = 1, with_ties = FALSE) |>
  ungroup()

# same colors you liked
dom_cols <- c(
  "class_1" = "red",
  "class_2" = "steelblue",
  "class_3" = "forestgreen",
  "class_4" = "purple",
  "class_5" = "goldenrod"
)

ggplot(dominant_norm) +
  geom_sf(aes(fill = class), color = NA) +                    # filled regions
  geom_sf(data = st_transform(zip_poly, 5070),                # ZIP outlines
          fill = NA, color = "grey60", linewidth = 0.15) +
  scale_fill_manual(values = dom_cols,
                    name   = "Dominant class",
                    drop   = FALSE) +
  coord_sf() +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_line(color = "grey90", linewidth = 0.25),
    axis.text  = element_blank()
  ) +
  labs(title = "Dominant latent class (hierarchical, sums to 1)")

##
# Fit and store models instead of only returning predictions
fit_one_model <- function(col, k = 80, fam = quasibinomial()) {
  fml <- as.formula(paste0(col, " ~ s(x, y, k = ", k, ")"))
  gam(fml, data = dat, family = fam, method = "REML")
}
mods <- lapply(prob_cols, fit_one_model)
names(mods) <- sub("^p_cluc[_.]", "class_", prob_cols)

# Extract edf, F, p-value, deviance explained
extract_summ <- function(m){
  sm <- summary(m)
  st <- sm$s.table[1, , drop = FALSE]   # s(x,y) row
  data.frame(
    edf   = st[,"edf"],
    F     = st[,"F"],
    pval  = st[,"p-value"],
    dev_expl = sm$dev.expl,
    r2_adj   = sm$r.sq
  )
}
summ_tbl <- do.call(rbind, lapply(mods, extract_summ))
summ_tbl
## ============================================================
## ============================================================
##  B) INCREMENTAL TESTS (NULL vs SPATIAL)  — COEFFICIENTS
## ============================================================

stopifnot(exists("mods"))              # your spatial GAMs list
stopifnot(exists("coef_cols"))         # the columns you fit

# pretty class labels like "Class 1", "Class 2", ...
pretty_class <- function(nm) paste("Class", sub(".*_(\\d+)$", "\\1", nm))

# 1) GAM diagnostics for each spatial model
extract_summ <- function(m, nm){
  sm <- summary(m)
  st <- sm$s.table[1, , drop = FALSE]     # row for s(x,y)
  data.frame(
    class    = pretty_class(nm),
    edf      = st[,"edf"],
    F        = st[,"F"],
    pval     = st[,"p-value"],
    dev_expl = sm$dev.expl,
    r2_adj   = sm$r.sq,
    row.names = NULL
  )
}
summ_tbl <- do.call(rbind, Map(extract_summ, mods, names(mods)))

# 2) Null (intercept-only) models with gaussian family
mods_null <- setNames(lapply(coef_cols, function(col)
  gam(as.formula(paste0(col, " ~ 1")),
      data = dat, family = gaussian(), method = "REML")),
  names(mods)
)

# 3) Incremental fit: null vs spatial
inc_tbl <- do.call(rbind, Map(function(m0, m1, nm){
  at <- anova(m0, m1, test = "F")
  data.frame(
    class    = pretty_class(nm),
    df_diff  = at$Df[2],
    dev_diff = at$Deviance[2],
    F_inc    = at$F[2],
    p_inc    = at$`Pr(>F)`[2],
    row.names = NULL
  )
}, mods_null, mods, names(mods)))

# 4) Merge to final table
final_tbl <- dplyr::left_join(summ_tbl, inc_tbl, by = "class")

# Optional: format like your example (scientific for p-values)
fmt_e <- function(x) ifelse(is.na(x), NA, formatC(x, format = "e", digits = 2))
final_tbl_out <- final_tbl |>
  dplyr::mutate(
    pval  = fmt_e(pval),
    p_inc = fmt_e(p_inc)
  )

final_tbl_out

