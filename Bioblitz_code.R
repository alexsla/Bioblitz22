### NEW BIOBLITZ CODE ###
###   Alex Slavenko   ###

library(tidyverse)
library(sf)
library(terra)
library(eks)
library(scatterpie)
library(ggridges)
library(ggforce)
library(stringr)
library(ggVennDiagram)
library(lubridate)
library(emmeans)
library(ggsignif)
library(ggfortify)
library(scales)
source("R/bb_functions.R")

## LOAD DATA
# Load the clean dataset of biodiversity records, and the cypher file
dat_all <- read_csv("data/Records v1 7Mar23.csv")

cypher <- read_csv("data/Morphospecies cypher v1 7Mar23.csv")

synonyms <- read_csv("data/Synonyms v1 7Mar23.csv")

dat_all <- dat_all %>%
  left_join(synonyms, by = c("morphospecies" = "synonym")) %>%
  mutate(
    morphospecies = case_when(is.na(accepted) ~ morphospecies,
                              T ~ accepted),
    std_date = parse_date_time(std_date, c("dmy", "ymd")),
    record_period = case_when(
      year(std_date) == 2022 & month(std_date) == 9 ~ "bb22",
      year(std_date) < 1992 ~ "historical",
      T ~ "recent"
    ),
    occ = 1,
    year = year(std_date)
  ) %>%
  dplyr::select(-accepted) %>%
  left_join(cypher) %>%
  drop_na(maingroup) %>%
  filter(!maingroup %in% c("Lampreys", "Protozoan")) %>%
  mutate(std_date = as.character(std_date),
         X = row_number())

councils_full <- unique(dat_all$council)

# generate discrete colour palette
bb_palette <- tibble(
  maingroup = sort(unique(dat_all$maingroup)),
  col = c(
    "#8C468C",
    "#ED1B24",
    "#404254",
    "#51340B",
    "#2A764F",
    "#EDC919",
    "#4C86C6",
    "#4563AC",
    "#DBA527",
    "#51B747",
    "#00008B",
    "#F7861C"
  )
)

# generate gradient colour palettes
taxa <-
  c(
    "Arachnids",
    "Birds",
    "Cartilaginous fishes",
    "Chromists",
    "Frogs",
    "Fungi",
    "Insects",
    "Mammals",
    "Other invertebrates",
    "Plants",
    "Ray-finned fishes",
    "Reptiles"
  )
palettes <-
  list(
    Arachnids = c("#ddcedd", "#c9acc9", "#b58ab5", "#a168a0", "#8c468c"),
    Birds = c("#f0d4d5", "#faadae", "#fc8583", "#f85a55", "#ed1b24"),
    `Cartilaginous fishes` = c("#d1d2d6", "#adaeb6", "#8c8e99", "#666877", "#404254"),
    Chromists = c("#dad4cc", "#c2b8aa", "#a0907a", "#7c6648", "#51340b"),
    Frogs = c("#bcdccb", "#97c2ab", "#74a88b", "#508f6d", "#2a764f"),
    Fungi = c("#f2edd4", "#f1e5ac", "#f0dc83", "#eed357", "#edc919"),
    Insects = c("#d6e1ee", "#b4cae4", "#93b3da", "#719cd0", "#4c86c6"),
    Mammals = c("#d6e1ee", "#adc1de", "#87a2ce", "#6582be", "#4563ac"),
    `Other invertebrates` = c("#ece5d3", "#e8d6aa", "#e4c681", "#e0b657", "#dba527"),
    Plants = c("#d8ecd6", "#b8dfb3", "#97d290", "#76c56c", "#51b747"),
    `Ray-finned fishes` = c("#bdbde1", "#8c8ccb", "#5959b4", "#2b2b9f", "#00008b"),
    Reptiles = c("#fbebdc", "#fdd3ac", "#fdba7d", "#fba150", "#f7861c")
  )

# load map of Victoria and reproject to Lambert equal-area
map_VIC <- read_sf("SpatialData/AD_LGA_AREA_POLYGON.shp") %>%
  filter(STATE == "VIC")

# filter based on LGAs
map_LGA <- map_VIC %>%
  st_transform(crs = 3112) %>%
  filter(NAME %in% str_to_upper(councils_full)) %>%
  rename(council = NAME) %>%
  mutate(
    council = case_when(
      council == "BAYSIDE" ~ "BSD",
      council == "CARDINIA" ~ "CAR",
      council == "CASEY" ~ "CAS",
      council == "FRANKSTON" ~ "FRA",
      council == "GREATER DANDENONG" ~ "GDA",
      council == "KINGSTON" ~ "KIN",
      council == "KNOX" ~ "KNO",
      council == "MONASH" ~ "MON",
      council == "MORNINGTON PENINSULA" ~ "MPE",
      council == "YARRA RANGES" ~  "YRA"
    )
  )

LGA_coords <- st_centroid(map_LGA) %>%
  st_coordinates() %>%
  as_tibble() %>%
  add_column(council = map_LGA$council)

# create grid at 500m resolution
# grid_LGA <- st_make_grid(map_LGA,
#                          cellsize = c(500, 500))  %>%
#   st_as_sf() %>%
#   filter(st_intersects(st_geometry(.), summarise(map_LGA), sparse = F)[,1]) %>%
#   mutate(FID = row_number())
#
# write_sf(grid_LGA, "SpatialData/grid_LGA.shp")

grid_LGA <- read_sf("SpatialData/grid_LGA.shp")

# filter records outside of LGAs
dat_all <- dat_all %>%
  select(-council) %>%
  left_join(
    st_as_sf(
      .,
      coords = c("decimalLongitude", "decimalLatitude"),
      crs = "WGS84"
    ) %>%
      st_transform(crs = 3112) %>%
      st_join(map_LGA) %>%
      as_tibble() %>%
      select(X, council),
    by = "X"
  ) %>%
  filter(!is.na(council))

# convert data to grid form
occ_GRID <- dat_all %>%
  filter(record_period %in% c("recent", "bb22")) %>%
  select(decimalLatitude,
         decimalLongitude,
         recordID,
         record_period,
         maingroup) %>%
  st_as_sf(
    coords = c("decimalLongitude", "decimalLatitude"),
    crs = st_crs(map_VIC)
  ) %>%
  st_transform(crs = 3112) %>%
  st_intersection(grid_LGA) %>%
  as_tibble() %>%
  select(-geometry) %>%
  group_by(FID, maingroup, record_period) %>%
  summarise(occ = n())

# load map of green spaces and reproject
greenspaces <- read_sf("data/VPA_Draft_Open_Space_Data.shp") %>%
  filter(OS_TYPE %in% c("Private open space", "Public open space")) %>%
  mutate(
    council = case_when(
      LGA == "BAYSIDE" ~ "BSD",
      LGA == "CARDINIA" ~ "CAR",
      LGA == "CASEY" ~ "CAS",
      LGA == "FRANKSTON" ~ "FRA",
      LGA == "GREATER DANDENONG" ~ "GDA",
      LGA == "KINGSTON" ~ "KIN",
      LGA == "KNOX" ~ "KNO",
      LGA == "MONASH" ~ "MON",
      LGA == "MORNINGTON" ~ "MPE",
      LGA == "YARRA RANGES" ~  "YRA"
    )
  ) %>%
  filter(!is.na(council)) %>%
  st_transform(st_crs(map_LGA))

# counts per green space
occ_greenspace <- dat_all %>%
  filter(record_period %in% c("recent", "bb22")) %>%
  select(decimalLatitude,
         decimalLongitude,
         recordID,
         record_period,
         maingroup) %>%
  st_as_sf(
    coords = c("decimalLongitude", "decimalLatitude"),
    crs = st_crs(map_VIC)
  ) %>%
  st_transform(crs = 3112) %>%
  st_intersection(greenspaces) %>%
  as_tibble() %>%
  select(-geometry) %>%
  group_by(FID, maingroup, record_period) %>%
  summarise(occ = n())

sp_greenspace <- dat_all %>%
  filter(record_period %in% c("recent", "bb22")) %>%
  select(
    decimalLatitude,
    decimalLongitude,
    recordID,
    record_period,
    maingroup,
    morphospecies
  ) %>%
  st_as_sf(
    coords = c("decimalLongitude", "decimalLatitude"),
    crs = st_crs(map_VIC)
  ) %>%
  st_transform(crs = 3112) %>%
  st_intersection(greenspaces) %>%
  as_tibble() %>%
  select(FID, maingroup, record_period, morphospecies) %>%
  unique() %>%
  group_by(FID, maingroup, record_period) %>%
  summarise(occ = n())

# # load participant data
# participants <- read_csv("data/BB21_participants_FOR_COUNCILS_v1 13July2022.csv")

## GENERATE DELIVERABLES: TOTAL AND PER COUNCIL
councils <- c("ALL", unique(dat_all$council))

for (i in councils) {
  if (i == "ALL") {
    dat <- dat_all
    grid_dat <- grid_LGA
    green_dat <- greenspaces
    LGA_dat <- map_LGA
  } else {
    dat <- dat_all %>% filter(council == i)
    grid_dat <- grid_LGA %>%
      filter(st_intersects(st_geometry(.), map_LGA %>% filter(council == i), sparse = F)[, 1])
    green_dat <- greenspaces %>% filter(council == i)
    LGA_dat <- map_LGA %>% filter(council == i)
  }
  
  # LEGEND MAP
  p_legend <- map_VIC %>%
    ggplot() +
    geom_sf(fill = "grey95", colour = "grey75") +
    geom_sf(
      data = greenspaces,
      fill = "darkolivegreen",
      colour = NA,
      alpha = .8
    ) +
    geom_sf(
      data = LGA_dat,
      fill = NA,
      colour = "black",
      size = 2
    ) +
    theme_bw() +
    coord_sf(crs = 3112,
             xlim = st_bbox(LGA_dat)[c(1, 3)],
             ylim = st_bbox(LGA_dat)[c(2, 4)]) +
    theme(
      panel.background = element_rect(fill = 'lightblue'),
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  ggsave(
    paste0("outputs/Maps/legend_map_", i, ".png"),
    plot = p_legend,
    width = 2000,
    height = 1500,
    units = "px",
    dpi = 300
  )
  
  # TABLE 1: species list, Nrecords total/period, FFG/EPBC, origin, common name, etc...
  Splist.Nrecords(dat,
                  full_dat = dat_all,
                  cypher = cypher,
                  path = "outputs/Table 1",
                  i)
  
  ## TABLE 2: emerging stories [lazarus, new, new_introduced]
  emerging(
    file = paste0("outputs/Table 1/Splist_Nrecords_", i, ".csv"),
    path = "outputs/Table 2",
    i
  )
  
  ## TABLE 3: species richness by main group, by period, unique species
  SpRichness(
    file = paste0("outputs/Table 1/Splist_Nrecords_", i, ".csv"),
    path = "outputs/Table 3",
    group = maingroup,
    i
  )
  
  ## TABLE 4: species richness by subgroup, by period, unique species
  SpRichness(
    file = paste0("outputs/Table 1/Splist_Nrecords_", i, ".csv"),
    path = "outputs/Table 4",
    group = subgroup,
    i
  )
  
  ## TABLE 7: most common species from BB22
  Common.sp(
    file = paste0("outputs/Table 1/Splist_Nrecords_", i, ".csv"),
    path = "outputs/Table 7",
    i
  )
  
  ## FIGURE 1
  Common.sp.p(
    file = paste0("outputs/Table 7/Common_species_", i, ".csv"),
    path = "outputs/Figure 1",
    i,
    palette = bb_palette
  )
  
  ## TABLE 8: most common group from BB22
  Common.group(
    file = paste0("outputs/Table 3/Species_richness_by_maingroup_", i, ".csv"),
    path = "outputs/Table 8",
    group = maingroup,
    i
  )
  
  ## FIGURE 2
  Common.group.p(
    file = paste0("outputs/Table 8/Common_maingroup_", i, ".csv"),
    path = "outputs/Figure 2",
    i,
    group = maingroup,
    count = Nrecords_bioblitz,
    palette = bb_palette,
    yaxis = T
  )
  
  Common.group.p(
    file = paste0("outputs/Table 8/Common_maingroup_", i, ".csv"),
    path = "outputs/Figure 2",
    i,
    group = maingroup,
    count = Nuniquesp_bioblitz,
    palette = bb_palette,
    yaxis = T
  )
  
  ## TABLE 9: most common subgroups from BB22
  Common.group(
    file = paste0("outputs/Table 4/Species_richness_by_subgroup_", i, ".csv"),
    path = "outputs/Table 9",
    group = subgroup,
    i
  )
  
  ## FIGURE 3
  Common.group.p(
    file = paste0("outputs/Table 9/Common_subgroup_", i, ".csv"),
    path = "outputs/Figure 3",
    i,
    group = subgroup,
    count = Nrecords_bioblitz,
    palette = bb_palette,
    yaxis = T
  )
  
  Common.group.p(
    file = paste0("outputs/Table 9/Common_subgroup_", i, ".csv"),
    path = "outputs/Figure 3",
    i,
    group = subgroup,
    count = Nuniquesp_bioblitz,
    palette = bb_palette,
    yaxis = T
  )
  
  ## TABLE 10: BB22 summaries by greenspace
  Nrecords.GS(
    dat = green_dat,
    occ_dat = occ_greenspace,
    sp = sp_greenspace,
    path = "outputs/Table 10",
    i
  )
  
  # ## TABLE 11: BB22 summary of participants
  # Nrecords.PAR(dat,
  #              path = "outputs/Table 11",
  #              i)
  #
  # ## FIGURE 4
  # Common.group.p(file = paste0("outputs/Table 11/BBparticipants_", i, ".csv"),
  #                path = "outputs/Figure 4",
  #                i,
  #                group = Participant_username,
  #                count = Nrecords_bb21,
  #                palette = bb_palette,
  #                percent = T)
  #
  # ## TABLE 12: BB22 summary for groups of participants
  # BBparticipants_type(dat,
  #                     participants %>% select(-c(Nrecords_BB, Perc_records_BB)),
  #                     path = "outputs/Table 12",
  #                     i)
  #
  # ## FIGURE 5
  # Participant.p(file = paste0("outputs/Table 12/BBparticipants_type_", i, ".csv"),
  #               path = "outputs/Figure 5",
  #               i)
  #
  # ## TABLE 13: BB22 summary of record bins and n participants
  # BBparticipants_bin(dat,
  #                    path = "outputs/Table 13",
  #                    i)
  #
  # ## TABLE 14: BB22 summary for (taxa) main groups
  # BBmaingroup(dat,
  #             path = "outputs/Table 14",
  #             i)
  
  ## FIGURE X: cumulative known species
  Trend.sp(
    file = paste0("outputs/Table 1/Splist_Nrecords_", i, ".csv"),
    path = "outputs/Figure X",
    i,
    palette = bb_palette
  )
  
  ## FIGURE X2: no. observations
  Trend.obs(dat,
            path = "outputs/Figure X2",
            i,
            palette = bb_palette)
  
  ## MAPS
  Map.p(
    dat,
    grid_dat,
    occ_dat = occ_GRID,
    basemap = map_LGA,
    path = "outputs/Maps/Total",
    i,
    period = "recent"
  )
  
  # Map.p(dat,
  #       grid_dat,
  #       occ_dat = occ_GRID,
  #       basemap = map_LGA,
  #       path = "outputs/Maps/Total",
  #       i,
  #       period = "recent",
  #       type = "Ndensity")
  
  Map.p(
    dat,
    grid_dat,
    occ_dat = occ_GRID,
    basemap = map_LGA,
    path = "outputs/Maps/Total",
    i,
    period = "bb22"
  )
  
  # Map.p(dat,
  #       grid_dat,
  #       occ_dat = occ_GRID,
  #       basemap = map_LGA,
  #       path = "outputs/Maps/Total",
  #       i,
  #       period = "bb22",
  #       type = "Ndensity")
  
  Map.gs(
    dat,
    green_dat,
    occ_dat = occ_greenspace,
    basemap = map_LGA,
    path = "outputs/Maps/Total",
    i,
    period = "bb22"
  )
  
  Map.gs(
    dat,
    green_dat,
    occ_dat = occ_greenspace,
    basemap = map_LGA,
    path = "outputs/Maps/Total",
    i,
    period = "bb22",
    type = "Greenspace"
  )
  
  # Produce maps per maingroup
  for (k in taxa) {
    Map.p(
      dat %>% filter(maingroup == k),
      grid_dat,
      occ_dat = occ_GRID %>% filter(maingroup == k),
      basemap = map_LGA,
      path = paste0("outputs/Maps/", k),
      i,
      period = "recent",
      palette = palettes[[k]]
    )
    
    # Map.p(dat %>% filter(maingroup == k),
    #       grid_dat,
    #       occ_dat = occ_GRID %>% filter(maingroup == k),
    #       basemap = map_LGA,
    #       path = paste0("outputs/Maps/", k),
    #       i,
    #       period = "recent",
    #       type = "Ndensity",
    #       palette = palettes[[k]])
    
    Map.p(
      dat %>% filter(maingroup == k),
      grid_dat,
      occ_dat = occ_GRID %>% filter(maingroup == k),
      basemap = map_LGA,
      path = paste0("outputs/Maps/", k),
      i,
      period = "bb22",
      palette = palettes[[k]]
    )
    
    # Map.p(dat %>% filter(maingroup == k),
    #       grid_dat,
    #       occ_dat = occ_GRID %>% filter(maingroup == k),
    #       basemap = map_LGA,
    #       path = paste0("outputs/Maps/", k),
    #       i,
    #       period = "bb22",
    #       type = "Ndensity",
    #       palette = palettes[[k]])
    
    Map.gs(
      dat %>% filter(maingroup == k),
      green_dat,
      occ_dat = occ_greenspace %>% filter(maingroup == k),
      basemap = map_LGA,
      path = paste0("outputs/Maps/", k),
      i,
      period = "bb22",
      palette = palettes[[k]]
    )
    
    Map.gs(
      dat %>% filter(maingroup == k),
      green_dat,
      occ_dat = occ_greenspace %>% filter(maingroup == k),
      basemap = map_LGA,
      path = paste0("outputs/Maps/", k),
      i,
      period = "bb22",
      type = "Greenspace",
      palette = palettes[[k]]
    )
  }
}

## General outputs
Common.group_all <-
  lapply(dir("outputs/Table 8", full.names = T), read_csv)

Common.group_Nrecords <- bind_rows(lapply(1:length(councils),
                                          function(x)
                                            Common.group_all[[x]] %>%
                                            select(maingroup, Nrecords_bioblitz) %>%
                                            mutate(council = councils[[x]])))  %>%
  spread(maingroup, Nrecords_bioblitz) %>%
  mutate_if(is.numeric, funs(replace_na(., 0))) %>%
  left_join(LGA_coords %>% bind_rows(tibble(
    X = min(LGA_coords$X - 5000),
    Y = max(LGA_coords$Y + 5000),
    council = "ALL"
  )))

Common.group_Nsp <- bind_rows(lapply(1:length(councils),
                                     function(x)
                                       Common.group_all[[x]] %>%
                                       select(maingroup, Nsp_bioblitz) %>%
                                       mutate(council = councils[[x]])))  %>%
  spread(maingroup, Nsp_bioblitz) %>%
  mutate_if(is.numeric, funs(replace_na(., 0))) %>%
  left_join(LGA_coords %>% bind_rows(tibble(
    X = min(LGA_coords$X - 5000),
    Y = max(LGA_coords$Y + 5000),
    council = "ALL"
  )))

no_zero_groups <-
  colSums(Common.group_Nrecords[, bb_palette$maingroup])
no_zero_groups <- no_zero_groups[which(no_zero_groups > 0)]

plot_pie_map(
  p_map = map_LGA %>%
    ggplot() +
    geom_sf(
      fill = "white",
      colour = "dark grey",
      inherit.aes = F
    ) +
    geom_scatterpie(
      aes(
        x = X,
        y = Y,
        group = council,
        r = r / 2
      ),
      data = Common.group_Nrecords %>%
        mutate(r = rowSums(Common.group_Nrecords[, bb_palette$maingroup])),
      cols = bb_palette$maingroup,
      color = NA
    ) +
    theme_bw() +
    scale_fill_manual(
      values = as.character(
        bb_palette %>%
          filter(maingroup %in% names(no_zero_groups)) %>%
          select(col) %>%
          as_vector()
      ),
      name = "Group"
    ) +
    geom_scatterpie_legend(
      rowSums(Common.group_Nrecords[, bb_palette$maingroup]) / 3,
      x = max(LGA_coords$X) + 20000,
      y = min(LGA_coords$Y) - 5000,
      labeller = function(x)
        x * 2
    ) +
    labs(title = "Bioblitz 2022 No. observations",
         x = "Longitude",
         y = "Latitude"),
  dat = Common.group_Nrecords,
  palette = bb_palette,
  scale.fact = .5
)

ggsave(
  paste0("outputs/Summary/Common_groups_Nrecords.jpg"),
  width = 2000,
  height = 1500,
  units = "px",
  dpi = 300
)

plot_pie_map(
  p_map = map_LGA %>%
    ggplot() +
    geom_sf(
      fill = "white",
      colour = "dark grey",
      inherit.aes = F
    ) +
    geom_scatterpie(
      aes(
        x = X,
        y = Y,
        group = council,
        r = r * 4
      ),
      data = Common.group_Nsp %>%
        mutate(r = rowSums(Common.group_Nsp[, bb_palette$maingroup])),
      cols = bb_palette$maingroup,
      color = NA
    ) +
    theme_bw() +
    scale_fill_manual(
      values = as.character(
        bb_palette %>%
          filter(maingroup %in% names(no_zero_groups)) %>%
          select(col) %>%
          as_vector()
      ),
      name = "Group"
    ) +
    geom_scatterpie_legend(
      rowSums(Common.group_Nsp[, bb_palette$maingroup]) * 4,
      x = max(LGA_coords$X) + 20000,
      y = min(LGA_coords$Y) - 5000,
      labeller = function(x)
        round(x / 4,-1)
    ) +
    labs(title = "Bioblitz 2022 No. species",
         x = "Longitude",
         y = "Latitude"),
  dat = Common.group_Nsp,
  bb_palette,
  scale.fact = 4
)

ggsave(
  paste0("outputs/Summary/Common_groups_Nsp.jpg"),
  width = 2000,
  height = 1500,
  units = "px",
  dpi = 300
)

groups_no_ones <-
  dat_all %>% filter(record_period == "bb22") %>% group_by(maingroup, morphospecies) %>% count() %>% ungroup() %>% count(maingroup) %>% filter(n > 1) %>% pull(maingroup)

dat_all %>% filter(record_period == "bb22") %>% group_by(maingroup, morphospecies) %>% count() %>% filter(maingroup %in% groups_no_ones) %>%
  ggplot(aes(
    y = reorder(maingroup, n, FUN = median, decreasing = T),
    x = n,
    fill = maingroup,
    colour = maingroup
  )) +
  stat_density_ridges(
    quantile_lines = TRUE,
    quantiles = 2,
    alpha = .7,
    jittered_points = TRUE,
    position = position_points_jitter(width = 0.05, height = 0),
    point_shape = '|',
    point_size = 3,
    point_alpha = 1
  ) +
  scale_x_continuous(trans = "log10") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = as.character(
    bb_palette %>%
      filter(maingroup %in% groups_no_ones) %>%
      select(col) %>%
      as_vector()
  ),
  name = "Group") +
  scale_colour_manual(values = as.character(
    bb_palette %>%
      filter(maingroup %in% groups_no_ones) %>%
      select(col) %>%
      as_vector()
  ),
  name = "Group") +
  labs(title = "No. observations per species",
       x = "Observations",
       y = "Group")

ggsave(
  paste0("outputs/Summary/Nrecords_per_sp.jpg"),
  width = 2000,
  height = 1500,
  units = "px",
  dpi = 300
)

dat_all %>%
  filter(!record_period == "historical") %>%
  group_by(year, maingroup) %>%
  count(morphospecies) %>%
  summarise(
    upper = quantile(n, probs = .95),
    lower = quantile(n, probs = .05),
    n = mean(n)
  ) %>%
  ggplot(aes(x = year, colour = maingroup, fill = maingroup)) +
  geom_line(aes(y = n)) +
  geom_ribbon(aes(ymax = upper, ymin = lower),
              colour = NA,
              alpha = .25) +
  scale_colour_manual(values = bb_palette$col, name = "Group") +
  scale_fill_manual(values = bb_palette$col, name = "Group") +
  theme_bw() +
  facet_wrap(vars(maingroup), scales = "free_y", ncol = 3) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Year",
       y = "Mean no. obersvations per species")

ggsave(
  paste0("outputs/Summary/Nrecords_per_sp_per_year.jpg"),
  width = 2000,
  height = 1500,
  units = "px",
  dpi = 300
)

dat_all %>%
  filter(!record_period == "historical") %>%
  select(year, maingroup, morphospecies) %>%
  unique() %>%
  group_by(year, maingroup) %>%
  count() %>%
  summarise(
    upper = quantile(n, probs = .95),
    lower = quantile(n, probs = .05),
    n = mean(n)
  ) %>%
  ggplot(aes(x = year, y = n, colour = maingroup)) +
  geom_point() +
  geom_smooth(se = F) +
  scale_colour_manual(values = bb_palette$col, name = "Group") +
  theme_bw() +
  facet_wrap(vars(maingroup), scales = "free_y", ncol = 3) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Year",
       y = "No. species observed per year")

ggsave(
  paste0("outputs/Summary/Nsp_per_year.jpg"),
  width = 2000,
  height = 1500,
  units = "px",
  dpi = 300
)

dat_2022 <- dat_all %>%
  filter(year == 2022) %>%
  mutate(eventDate = ymd(std_date))

dat_2022 %>%
  group_by(eventDate, maingroup, record_period) %>%
  count() %>%
  ggplot(aes(
    x = eventDate,
    y = n,
    colour = maingroup,
    fill = maingroup,
    linetype = factor(record_period, labels = c("Bioblitz", "Baseline"))
  )) +
  geom_point(alpha = .1) +
  geom_smooth() +
  scale_colour_manual(values = bb_palette$col) +
  scale_fill_manual(values = bb_palette$col) +
  theme_bw() +
  facet_wrap(vars(maingroup), scales = "free_y", ncol = 3) +
  labs(
    y = "No. observations per day",
    x = "Date",
    colour = "Group",
    fill = "Group",
    linetype = "Record period"
  )

ggsave(
  paste0("outputs/Summary/Nrecords_2022.jpg"),
  width = 2000,
  height = 1500,
  units = "px",
  dpi = 300
)

dat_2022 %>%
  group_by(eventDate, maingroup, record_period, morphospecies) %>%
  count() %>%
  select(-n) %>%
  group_by(eventDate, maingroup, record_period) %>%
  count() %>%
  ggplot(aes(
    x = eventDate,
    y = n,
    colour = maingroup,
    fill = maingroup,
    linetype = factor(record_period, labels = c("Bioblitz", "Baseline"))
  )) +
  geom_point(alpha = .1) +
  geom_smooth() +
  scale_colour_manual(values = bb_palette$col) +
  scale_fill_manual(values = bb_palette$col) +
  theme_bw() +
  facet_wrap(vars(maingroup), scales = "free_y", ncol = 3) +
  labs(
    y = "No. species per day",
    x = "Date",
    colour = "Group",
    fill = "Group",
    linetype = "Record period"
  )

ggsave(
  paste0("outputs/Summary/Nsp_2022.jpg"),
  width = 2000,
  height = 1500,
  units = "px",
  dpi = 300
)

dat_september <- dat_all %>%
  filter(!record_period == "historical") %>%
  mutate(
    eventDate = ymd(std_date),
    month = month(eventDate),
    day = day(eventDate)
  ) %>%
  filter(month == 9)

dat_september %>%
  mutate(date = dmy(paste(day, month, "2022", sep = "-"))) %>%
  group_by(date, maingroup, record_period, year) %>%
  count() %>%
  ggplot(aes(
    x = factor(record_period, labels = c("Bioblitz", "Baseline")),
    y = n,
    fill = maingroup
  )) +
  geom_boxplot() +
  scale_fill_manual(values = bb_palette$col) +
  theme_bw() +
  facet_wrap(vars(maingroup), scales = "free_y", ncol = 3) +
  labs(y = "No. observations per day",
       x = "Record period",
       fill = "Group")

ggsave(
  paste0("outputs/Summary/Nrecords_aov.jpg"),
  width = 2000,
  height = 1500,
  units = "px",
  dpi = 300
)

sep_av <- aov(
  n ~ record_period * maingroup,
  dat = dat_september %>%
    mutate(date = dmy(paste(
      day, month, "2022", sep = "-"
    ))) %>%
    group_by(date, maingroup, record_period, year) %>%
    count()
)

summary(sep_av)

sep_res <-
  lsmeans(sep_av, pairwise ~ record_period * maingroup, adjust = "bonferroni")

sep_cont <-
  summary(sep_res[[2]])[!is.na(summary(sep_res[[2]])[, 4]), ] %>%
  filter(
    str_count(contrast, "Arachnids") == 2 |
      str_count(contrast, "Birds") == 2 |
      str_count(contrast, "Cartilaginous fishes") == 2 |
      str_count(contrast, "Chromists") == 2 |
      str_count(contrast, "Frogs") == 2 |
      str_count(contrast, "Fungi") == 2 |
      str_count(contrast, "Insects") == 2 |
      str_count(contrast, "Mammals") == 2 |
      str_count(contrast, "Other invertebrates") == 2 |
      str_count(contrast, "Plants") == 2 |
      str_count(contrast, "Ray-finned fishes") == 2 |
      str_count(contrast, "Reptiles") == 2
  )

sep_cont %>% write_csv("outputs/Summary/Nrecords_aov.csv")

dat_september %>%
  mutate(date = dmy(paste(day, month, "2022", sep = "-"))) %>%
  group_by(date, maingroup, record_period, year, morphospecies) %>%
  count() %>%
  select(-n) %>%
  group_by(date, maingroup, record_period, year) %>%
  count() %>%
  ggplot(aes(
    x = factor(record_period, labels = c("Bioblitz", "Baseline")),
    y = n,
    fill = maingroup
  )) +
  geom_boxplot() +
  scale_fill_manual(values = bb_palette$col) +
  theme_bw() +
  facet_wrap(vars(maingroup), scales = "free_y", ncol = 3) +
  labs(y = "No. sp per day",
       x = "Record period",
       fill = "Group")

ggsave(
  paste0("outputs/Summary/Nsp_aov.jpg"),
  width = 2000,
  height = 1500,
  units = "px",
  dpi = 300
)

sep_av <- aov(
  n ~ record_period * maingroup,
  dat = dat_september %>%
    mutate(date = dmy(paste(
      day, month, "2022", sep = "-"
    ))) %>%
    group_by(date, maingroup, record_period, year, morphospecies) %>%
    count() %>%
    select(-n) %>%
    group_by(date, maingroup, record_period, year) %>%
    count()
)

summary(sep_av)

sep_res <-
  lsmeans(sep_av, pairwise ~ record_period * maingroup, adjust = "bonferroni")

sep_cont <-
  summary(sep_res[[2]])[!is.na(summary(sep_res[[2]])[, 4]), ] %>%
  filter(
    str_count(contrast, "Arachnids") == 2 |
      str_count(contrast, "Birds") == 2 |
      str_count(contrast, "Cartilaginous fishes") == 2 |
      str_count(contrast, "Chromists") == 2 |
      str_count(contrast, "Frogs") == 2 |
      str_count(contrast, "Fungi") == 2 |
      str_count(contrast, "Insects") == 2 |
      str_count(contrast, "Mammals") == 2 |
      str_count(contrast, "Other invertebrates") == 2 |
      str_count(contrast, "Plants") == 2 |
      str_count(contrast, "Ray-finned fishes") == 2 |
      str_count(contrast, "Reptiles") == 2
  )

sep_cont %>% write_csv("outputs/Summary/Nsp_aov.csv")

## emerging story plot
all_emerging <- read_csv("outputs/Table 2/emerging_ALL.csv") %>%
  rowwise() %>%
  mutate(New = sum(Nnewnative, Nnewintro)) %>%
  select(-c(Nnew, Pnewnative, Pnewintro, Nnewnative, Nnewintro)) %>%
  rename(Rediscovered = Nlazarus) %>%
  full_join(
    read_csv("outputs/Table 1/Splist_Nrecords_ALL.csv") %>% mutate(extinct = case_when(
      is.na(first_record_after1991) ~ "Extinct"
    )) %>% drop_na(extinct) %>% count(maingroup) %>% rename(Extinct = n)
  ) %>%
  full_join(
    read_csv("outputs/Table 1/Splist_Nrecords_ALL.csv") %>% mutate(extant = case_when(
      is.na(new) &
        is.na(lazarus) &
        is.na(new_introduced) &
        !is.na(first_record_after1991) ~ "Extant"
    )) %>% drop_na(extant) %>% count(maingroup) %>% rename(Extant = n)
  ) %>%
  full_join(
    read_csv("outputs/Table 1/Splist_Nrecords_ALL.csv") %>% mutate(known = case_when(
      is.na(new) &
        is.na(new_introduced) ~ "Known biodversity"
    )) %>% drop_na(known) %>% count(maingroup) %>% rename(`Known biodiversity` = n)
  ) %>% full_join(
    read_csv("outputs/Table 1/Splist_Nrecords_ALL.csv") %>% mutate(
      extinct = case_when(is.na(first_record_after1991) ~ "Extinct"),
      updated = case_when(is.na(extinct) ~ "Updated biodiversity")
    ) %>% drop_na(updated) %>% count(maingroup) %>% rename(`Updated biodiversity` = n)
  ) %>%
  full_join(
    read_csv("outputs/Table 1/Splist_Nrecords_ALL.csv") %>% mutate(
      refound = case_when(
        is.na(new) &
          is.na(lazarus) &
          is.na(new_introduced) &
          !is.na(first_record_after1991) &
          most_recent_record == 2022 ~ "Re-found"
      )
    ) %>% drop_na(refound) %>% count(maingroup) %>% rename(`Re-found` = n)
  ) %>%
  full_join(
    read_csv("outputs/Table 1/Splist_Nrecords_ALL.csv") %>% mutate(
      nrefound = case_when(
        is.na(new) &
          is.na(lazarus) &
          is.na(new_introduced) &
          !is.na(first_record_after1991) &
          most_recent_record < 2022 ~ "Not re-found"
      )
    ) %>% drop_na(nrefound) %>% count(maingroup) %>% rename(`Not re-found` = n)
  ) %>%
  full_join(
    read_csv("outputs/Table 1/Splist_Nrecords_ALL.csv") %>% mutate(pextinct = case_when(
      is.na(first_record_after1991) |
        lazarus == "Lazarus" ~ "Possibly extinct"
    )) %>% drop_na(pextinct) %>% count(maingroup) %>% rename(`Possibly extinct` = n)
  ) %>%
  gather(story, n,-maingroup)

origin_emerging <-
  read_csv("outputs/Table 1/Splist_Nrecords_ALL.csv") %>% drop_na(lazarus) %>% count(Origin) %>% rename(Rediscovered = n) %>%
  full_join(
    read_csv("outputs/Table 1/Splist_Nrecords_ALL.csv") %>% drop_na(new) %>% count(Origin) %>% rename(New = n)
  ) %>%
  full_join(
    read_csv("outputs/Table 1/Splist_Nrecords_ALL.csv") %>% mutate(extinct = case_when(
      is.na(first_record_after1991) ~ "Extinct"
    )) %>% drop_na(extinct) %>% count(Origin) %>% rename(Extinct = n)
  ) %>%
  full_join(
    read_csv("outputs/Table 1/Splist_Nrecords_ALL.csv") %>% mutate(extant = case_when(
      is.na(new) &
        is.na(lazarus) &
        is.na(new_introduced) &
        !is.na(first_record_after1991) ~ "Extant"
    )) %>% drop_na(extant) %>% count(Origin) %>% rename(Extant = n)
  ) %>%
  full_join(
    read_csv("outputs/Table 1/Splist_Nrecords_ALL.csv") %>% mutate(known = case_when(
      is.na(new) &
        is.na(new_introduced) ~ "Known biodversity"
    )) %>% drop_na(known) %>% count(Origin) %>% rename(`Known biodiversity` = n)
  ) %>% full_join(
    read_csv("outputs/Table 1/Splist_Nrecords_ALL.csv") %>% mutate(
      extinct = case_when(is.na(first_record_after1991) ~ "Extinct"),
      updated = case_when(is.na(extinct) ~ "Updated biodiversity")
    ) %>% drop_na(updated) %>% count(Origin) %>% rename(`Updated biodiversity` = n)
  ) %>%
  full_join(
    read_csv("outputs/Table 1/Splist_Nrecords_ALL.csv") %>% mutate(
      refound = case_when(
        is.na(new) &
          is.na(lazarus) &
          is.na(new_introduced) &
          !is.na(first_record_after1991) &
          most_recent_record == 2022 ~ "Re-found"
      )
    ) %>% drop_na(refound) %>% count(Origin) %>% rename(`Re-found` = n)
  ) %>%
  full_join(
    read_csv("outputs/Table 1/Splist_Nrecords_ALL.csv") %>% mutate(
      nrefound = case_when(
        is.na(new) &
          is.na(lazarus) &
          is.na(new_introduced) &
          !is.na(first_record_after1991) &
          most_recent_record < 2022 ~ "Not re-found"
      )
    ) %>% drop_na(nrefound) %>% count(Origin) %>% rename(`Not re-found` = n)
  ) %>%
  full_join(
    read_csv("outputs/Table 1/Splist_Nrecords_ALL.csv") %>% mutate(pextinct = case_when(
      is.na(first_record_after1991) |
        lazarus == "Lazarus" ~ "Possibly extinct"
    )) %>% drop_na(pextinct) %>% count(Origin) %>% rename(`Possibly extinct` = n)
  ) %>%
  gather(story, n,-Origin)

stories <- unique(all_emerging$story)

for (i in stories) {
  story <- all_emerging %>% filter(story == i, n > 0)
  story_origin <-
    origin_emerging %>% filter(story == i, n > 0) %>% drop_na(Origin)
  
  p_story <- story %>%
    ggplot() +
    geom_arc_bar(
      aes(
        x0 = 0,
        y0 = 0,
        r0 = 0.6,
        r = 1,
        amount = n,
        fill = maingroup
      ),
      stat = "pie",
      colour = NA
    ) +
    geom_text(x = 0,
              y = 0,
              aes(label = sum(n)),
              size = 20) +
    scale_fill_manual(
      values = filter(bb_palette, maingroup %in% story$maingroup)$col,
      name = "Group"
    ) +
    coord_fixed() +
    theme_void() +
    theme(legend.position = "none")
  
  if (i %in% c("Known biodiversity", "Updated biodiversity")) {
    library(grImport2)
    library(rphylopic)
    
    phylopics <-
      list(
        Arachnids = png::readPNG("phylopic/Arachnids.png"),
        Birds = png::readPNG("phylopic/Birds.png"),
        `Cartilaginous fishes` = png::readPNG("phylopic/Cartilaginous fishes.png"),
        Chromists = png::readPNG("phylopic/Chromists.png"),
        Frogs = png::readPNG("phylopic/Frogs.png"),
        Fungi = png::readPNG("phylopic/Fungi.png"),
        Insects = png::readPNG("phylopic/Insects.png"),
        Mammals = png::readPNG("phylopic/Mammals.png"),
        `Other invertebrates` = png::readPNG("phylopic/Other invertebrates.png"),
        Plants = png::readPNG("phylopic/Plants.png"),
        `Ray-finned fishes` = png::readPNG("phylopic/Ray-finned fishes.png"),
        Reptiles = png::readPNG("phylopic/Reptiles.png")
      )
    
    coords_ALL <- data.frame(X = 0, Y = 0)
    
    radius_ALL <- 1.05
    
    angles_ALL <- calc_angles(story$n)
    
    angles_ALL <- (angles_ALL + replace_na(lag(angles_ALL), 0)) / 2
    
    palette <- bb_palette %>%
      mutate(
        angle = angles_ALL,
        x = coords_ALL$X + radius_ALL * sin(pi * 2 * angle / 360),
        y = coords_ALL$Y + radius_ALL * cos(pi * 2 * angle / 360)
      )
    
    for (k in story$maingroup) {
      p_story <- p_story +
        add_phylopic(
          phylopics[[k]],
          x = as.numeric(filter(palette, maingroup == k)[, 4]),
          y = as.numeric(filter(palette, maingroup == k)[, 5]),
          ysize = 0.1,
          color = as.character(filter(palette, maingroup == k)[, 2]),
          alpha = 1
        )
    }
  }
  
  p_story_origin <- story_origin %>%
    ggplot() +
    geom_arc_bar(
      aes(
        x0 = 0,
        y0 = 0,
        r0 = 0.6,
        r = 1,
        amount = n,
        fill = Origin
      ),
      stat = "pie",
      colour = NA
    ) +
    geom_text(x = 0,
              y = 0,
              aes(label = sum(n)),
              size = 20) +
    scale_fill_manual(values = c("deepskyblue4", "darkorchid4"),
                      name = "Origin") +
    coord_fixed() +
    theme_void() +
    theme(legend.position = "none")
  
  ggsave(
    paste0("outputs/Summary/emerging_", i, ".pdf"),
    plot = p_story,
    width = 2000,
    height = 1500,
    units = "px",
    dpi = 300
  )
  
  ggsave(
    paste0("outputs/Summary/origin_", i, ".pdf"),
    plot = p_story_origin,
    width = 2000,
    height = 1500,
    units = "px",
    dpi = 300
  )
}

# ## participants summary plots
# participants %>%
#   mutate(Participant_type2 = case_when(Participant_type == "BB organiser" ~ BBorganiser_type, T ~ Participant_type)) %>%
#   ggplot(aes(x = factor(Participant_type, levels = c("General public", "BB organiser")),
#              y = Nrecords_BB,
#              fill = factor(Participant_type2, levels = c("General public", "Council staff", "Members of supporting associations")))) +
#   geom_boxplot() +
#   scale_y_continuous(trans = "log10") +
#   theme_bw() +
#   scale_fill_viridis_d(option = "mako") +
#   labs(x = "",
#        y = "Number records / participant",
#        fill = "") +
#   theme(legend.position = "bottom")
#
# ggsave(paste0("outputs/Summary/Nrecords_by_participant_type.jpg"),
#        width = 2000,
#        height = 1500,
#        units = "px",
#        dpi = 300)
#
# dALL.cleanup.groups_corrected %>%
#   filter(record_period == "bb22") %>%
#   left_join(participants, by = c("identifiedBy" = "Participant_username")) %>%
#   mutate(Participant_type2 = case_when(Participant_type == "BB organiser" ~ BBorganiser_type, T ~ Participant_type)) %>%
#   select(identifiedBy, Participant_type2, Participant_type, morphospecies) %>%
#   drop_na() %>%
#   unique() %>%
#   group_by(identifiedBy, Participant_type2, Participant_type) %>%
#   count(identifiedBy) %>%
#   ggplot(aes(x = factor(Participant_type, levels = c("General public", "BB organiser")),
#              y = n,
#              fill = factor(Participant_type2, levels = c("General public", "Council staff", "Members of supporting associations")))) +
#   geom_boxplot() +
#   scale_y_continuous(trans = "log10") +
#   theme_bw() +
#   scale_fill_viridis_d(option = "mako") +
#   labs(x = "",
#        y = "Number species / participant",
#        fill = "") +
#   theme(legend.position = "bottom")
#
# ggsave(paste0("outputs/Summary/Nsp_by_participant_type.jpg"),
#        width = 2000,
#        height = 1500,
#        units = "px",
#        dpi = 300)
#
# participants %>%
#   mutate(Participant_type2 = case_when(Participant_type == "BB organiser" ~ BBorganiser_type, T ~ Participant_type)) %>%
#   arrange(desc(Nrecords_BB)) %>%
#   mutate(Nrecords_tot = cumsum(Nrecords_BB), Participant_username = reorder(Participant_username, Nrecords_tot)) %>%
#   ggplot(aes(x = Participant_username,
#              y = Nrecords_tot,
#              fill = factor(Participant_type2, levels = c("General public", "Council staff", "Members of supporting associations")))) +
#   geom_bar(stat = "identity", colour = NA) +
#   theme_bw() +
#   scale_fill_viridis_d(option = "mako") +
#   labs(x = "",
#        y = "Accumulated number of records",
#        fill = "") +
#   theme(axis.ticks.x = element_blank(),
#         axis.text.x = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         legend.position = "bottom")
#
# ggsave(paste0("outputs/Summary/Nrecords_by_participants.jpg"),
#        width = 2000,
#        height = 1500,
#        units = "px",
#        dpi = 300)
#
# ## assess significance of overlap in sampled species between participant types (are different types of participants recording similar species?) - if null rejected overlap is lower than expected by chance, suggesting examined participant type locates different species
#
# sig_overlap <- list()
# for (k in c("Total", names(palettes))){
#   if(k == "Total") {
#     venn_data <- lapply(c("General public", "Council staff", "Members of supporting associations"),
#                         function(x)
#                           dALL.cleanup.groups_corrected %>%
#                           filter(record_period == "bb22") %>%
#                           left_join(participants, by = c("identifiedBy" = "Participant_username")) %>%
#                           mutate(Participant_type2 = case_when(Participant_type == "BB organiser" ~ BBorganiser_type, T ~ Participant_type)) %>%
#                           select(identifiedBy, Participant_type2, Participant_type, morphospecies) %>%
#                           drop_na() %>%
#                           unique() %>%
#                           filter(Participant_type2 == x) %>%
#                           pull(morphospecies) %>%
#                           unique())
#   } else {
#     venn_data <- lapply(c("General public", "Council staff", "Members of supporting associations"),
#                         function(x)
#                           dALL.cleanup.groups_corrected %>%
#                           filter(record_period == "bb22", maingroup == k) %>%
#                           left_join(participants, by = c("identifiedBy" = "Participant_username")) %>%
#                           mutate(Participant_type2 = case_when(Participant_type == "BB organiser" ~ BBorganiser_type, T ~ Participant_type)) %>%
#                           select(identifiedBy, Participant_type2, Participant_type, morphospecies) %>%
#                           drop_na() %>%
#                           unique() %>%
#                           filter(Participant_type2 == x) %>%
#                           pull(morphospecies) %>%
#                           unique())
#   }
#
#   names(venn_data) <- c("General public", "BB organiser - Council", "BB organiser - Supporting")
#
#   venn_p <- ggVennDiagram(venn_data,
#                           label_percent_digit = 1,
#                           set_size = 3) +
#     scale_fill_viridis_c(option = "plasma") +
#     scale_colour_viridis_d(option = "mako") +
#     labs(fill = "No. species") +
#     theme(legend.position = "bottom") +
#     scale_x_continuous(expand = expansion(mult = .2))
#
#   venn_p
#
#   ggsave(paste0("outputs/Summary/Nsp_by_participants_Venn_", k, ",.jpg"),
#          width = 2000,
#          height = 1500,
#          units = "px",
#          dpi = 300)
#
#   venn_groups <- as.numeric(str_remove(str_extract(layer_data(venn_p, 4)$label, ".*\n"), "\n"))
#
#   # calculate significance of overlap
#   p1 <- calc_enrich(venn_groups,
#                     group = "General public",
#                     plot = T)
#
#   p2 <- calc_enrich(venn_groups,
#                     group = "BB organiser - Council",
#                     plot = T)
#
#   p3 <- calc_enrich(venn_groups,
#                     group = "BB organiser - Supporting",
#                     plot = T)
#
#   cowplot::plot_grid(p1$plot, p2$plot, p3$plot,
#                      ncol = 1)
#
#   ggsave(paste0("outputs/Summary/Nsp_by_participants_overlap_", k, ",.jpg"),
#          width = 2000,
#          height = 1500,
#          units = "px",
#          dpi = 300)
#
#   sig_overlap[[k]] <- tibble(participant_type = c("General public", "BB organiser - Council", "BB organiser - Supporting"),
#                              maingroup = rep(k, 3),
#                              p = c(p1$p, p2$p, p3$p)) %>%
#     mutate(sig = case_when(p <= 0.001 ~ "***",
#                            p <= 0.01 ~ "**",
#                            p <= 0.05 ~ "*",
#                            p > 0.05 ~ "n.s."))
# }
#
# bind_rows(sig_overlap) %>% write_csv("outputs/Summary/Nsp_overlap_hypergeometric_tests.csv")


# update events
dat_events <- read_csv("data/LGAs_Events_dataset.csv") %>%
  mutate(
    lga = case_when(
      lga == "bayside" ~ "BSD",
      lga == "cardinia" ~ "CAR",
      lga == "casey" ~ "CAS",
      lga == "frankston" ~ "FRA",
      lga == "gd" ~ "GDA",
      lga == "kingston" ~ "KIN",
      lga == "knox" ~ "KNO",
      lga == "monash" ~ "MON",
      lga == "mps" ~ "MPE",
      lga == "yr" ~  "YRA"
    )
  )

dat_events_unique <- dat_events %>%
  select(event_year, event_description, start_date, end_date) %>%
  unique()

dat_all_events <- dat_all

dat_all_events$event_description <- "Baseline"

for (i in 1:nrow(dat_events_unique)) {
  event_interval <-
    interval(dmy(dat_events_unique[i, 3]), dmy(dat_events_unique[i, 4]))
  event_year <- as.character(dat_events_unique[i, 1])
  event_d <- as.character(dat_events_unique[i, 2])
  
  dat_all_events <- dat_all_events %>%
    mutate(
      record_period = case_when(
        ymd(std_date) %within% event_interval ~ event_year,
        T ~ record_period
      ),
      event_description = case_when(
        ymd(std_date) %within% event_interval ~ event_d,
        T ~ event_description
      )
    )
}

select(lga, event_year) %>%
  drop_na(lga) %>%
  mutate(event_year = factor(
    event_year,
    levels = c(
      "bb_20",
      "cnc_21",
      "bb_21",
      "gsb_21",
      "cnc_22",
      "bb_22",
      "gsb_22"
    )
  )) %>%
  group_by(lga) %>%
  arrange(event_year) %>%
  mutate(event_number = row_number()) %>%
  ungroup() %>%
  mutate(event_number = paste0("event_", letters[event_number])) %>%
  rename(council = lga,
         record_period = event_year)

dat_events_lga <- dat_events_lga %>%
  full_join(
    dat_events %>%
      select(lga, event_year) %>%
      rename(council = lga,
             record_period = event_year) %>%
      expand(council, record_period)
  ) %>%
  drop_na(council) %>%
  mutate(record_period = factor(
    record_period,
    levels = c(
      "bb_20",
      "inter_bb20_cnc21",
      "cnc_21",
      "inter_cnc21_bb21",
      "bb_21",
      "inter_bb21_gsb21",
      "gsb_21",
      "inter_gsb21_cnc22",
      "cnc_22",
      "inter_cnc22_bb22",
      "bb_22",
      "inter_bb22_gsb22",
      "gsb_22"
    )
  )) %>%
  arrange(council, record_period) %>%
  left_join(dat_events_lga) %>%
  group_by(council, event_number) %>%
  ungroup()

council_events <- list()
for (i in unique(dat_events_lga$council)) {
  council_events[[i]] <- dat_events_lga %>%
    filter(council == i)
  
  index = 1
  for (k in 1:nrow(council_events[[i]])) {
    if (is.na(council_events[[i]][k, 3])) {
      council_events[[i]][k, 3] <- paste0("interperiod_", letters[index])
    } else {
      index = index + 1
    }
  }
}

dat_events_lga <- bind_rows(council_events)

dat_all_events <- dat_all_events %>%
  left_join(dat_events_lga)

dat_all_events %>% write_csv("data/Records w Events v1 3Apr23.csv")


# create emerging stories summaries for each event and interperiod:
events_emerging <- list()
for (i in 1:nrow(dat_events_unique)) {
  event_i <- as.character(dat_events_unique[i, 1])
  
  event_year <- year(dmy(dat_events_unique[i, 3]))
  
  dat_event <- dat_all_events %>%
    filter(ymd(std_date) <= dmy(dat_events_unique[i, 4])) %>%
    mutate(
      record_period = case_when(
        record_period == event_i ~ "event",
        year < event_year - 30 ~ "historical",
        T ~ "recent"
      )
    )
  
  dat_event <- dat_event %>%
    group_by(morphospecies, record_period, council) %>%
    summarise(Nrecords_total = sum(occ)) %>%
    spread(record_period, Nrecords_total) %>%
    mutate_if(is.numeric, replace_na, 0) %>%
    mutate(
      total = sum(event, historical, recent, na.rm = T),
      Perc_records_total = total / nrow(dat_event) * 100,
      Perc_records_historical = historical / total * 100,
      Perc_records_recent = recent / total * 100,
      Perc_records_bioblitz = event / total * 100,
      distribution = case_when(
        Perc_records_bioblitz == 0 ~ 1,
        Perc_records_bioblitz == 100 ~ 2,
        T ~ 3
      )
    ) %>%
    ungroup() %>%
    left_join(dat_event %>%
                group_by(morphospecies) %>%
                summarise(most_recent_record = max(year))) %>%
    left_join(
      dat_event %>%
        filter(record_period %in% c("recent", "event")) %>%
        group_by(morphospecies) %>%
        summarise(first_record_after1991 = min(year))
    ) %>%
    left_join(cypher, by = c("morphospecies")) %>%
    mutate(
      lazarus = case_when(historical > 0 &
                            recent == 0 & event > 0 ~ "Lazarus"),
      new = case_when(historical == 0 &
                        recent == 0 & event > 0 ~ "New"),
      new_introduced = case_when(
        historical == 0 &
          recent == 0 &
          event > 0 & origin == "Introduced" ~ "New introduced"
      )
    ) %>%
    rename(
      CommonName = cname,
      Origin = origin,
      EPBC = epbc,
      FFGAct = ffga,
      Nrecords_total = total,
      Nrecords_historical = historical,
      Nrecords_recent = recent,
      Nrecords_bioblitz = event
    ) %>%
    mutate(CommonName = str_replace_all(str_to_title(str_replace_all(
      CommonName, "-", "_"
    )), "_", "-")) %>%
    select(
      council,
      maingroup,
      morphospecies,
      CommonName,
      Origin,
      family,
      genus,
      subgroup,
      EPBC,
      FFGAct,
      Nrecords_total,
      Perc_records_total,
      Nrecords_historical,
      Perc_records_historical,
      Nrecords_recent,
      Perc_records_recent,
      Nrecords_bioblitz,
      Perc_records_bioblitz,
      most_recent_record,
      first_record_after1991,
      distribution,
      lazarus,
      new,
      new_introduced
    )
  
  events_emerging[[i]] <-
    dat_event %>%
    mutate(
      new = case_when(new == "New" ~ 1,
                      T ~ 0),
      lazarus = case_when(lazarus == "Lazarus" ~ 1,
                          T ~ 0),
      new_introduced = case_when(new_introduced == "New introduced" ~ 1,
                                 T ~ 0)
    ) %>%
    group_by(council, maingroup) %>%
    summarise(
      Nnew = sum(new),
      Nnewintro = sum(new_introduced),
      Nlazarus = sum(lazarus)
    ) %>%
    mutate(
      Nnewnative = Nnew - Nnewintro,
      Pnewintro = (Nnewintro / Nnew) * 100,
      Pnewnative = (Nnewnative / Nnew) * 100
    ) %>%
    select(council,
           maingroup,
           Nnew,
           Nnewnative,
           Pnewnative,
           Nnewintro,
           Pnewintro,
           Nlazarus)  %>%
    rowwise() %>%
    mutate(New = sum(Nnewnative, Nnewintro)) %>%
    select(-c(Nnew, Pnewnative, Pnewintro, Nnewnative, Nnewintro)) %>%
    rename(Rediscovered = Nlazarus) %>%
    full_join(
      dat_event %>% mutate(extinct = case_when(
        is.na(first_record_after1991) ~ "Extinct"
      )) %>% drop_na(extinct) %>% count(council, maingroup) %>% rename(Extinct = n)
    ) %>%
    full_join(
      dat_event %>% mutate(extant = case_when(
        is.na(new) &
          is.na(lazarus) &
          is.na(new_introduced) &
          !is.na(first_record_after1991) ~ "Extant"
      )) %>% drop_na(extant) %>% count(council, maingroup) %>% rename(Extant = n)
    ) %>%
    full_join(
      dat_event %>% mutate(known = case_when(
        is.na(new) &
          is.na(new_introduced) ~ "Known biodversity"
      )) %>% drop_na(known) %>% count(council, maingroup) %>% rename(`Known biodiversity` = n)
    ) %>% full_join(
      dat_event %>% mutate(
        extinct = case_when(is.na(first_record_after1991) ~ "Extinct"),
        updated = case_when(is.na(extinct) ~ "Updated biodiversity")
      ) %>% drop_na(updated) %>% count(council, maingroup) %>% rename(`Updated biodiversity` = n)
    ) %>%
    full_join(
      dat_event %>% mutate(
        refound = case_when(
          is.na(new) &
            is.na(lazarus) &
            is.na(new_introduced) &
            !is.na(first_record_after1991) &
            most_recent_record == event_year ~ "Re-found"
        )
      ) %>% drop_na(refound) %>% count(council, maingroup) %>% rename(`Re-found` = n)
    ) %>%
    full_join(
      dat_event %>% mutate(
        nrefound = case_when(
          is.na(new) &
            is.na(lazarus) &
            is.na(new_introduced) &
            !is.na(first_record_after1991) &
            most_recent_record < event_year ~ "Not re-found"
        )
      ) %>% drop_na(nrefound) %>% count(council, maingroup) %>% rename(`Not re-found` = n)
    ) %>%
    full_join(
      dat_event %>% mutate(pextinct = case_when(
        is.na(first_record_after1991) |
          lazarus == "Lazarus" ~ "Possibly extinct"
      )) %>% drop_na(pextinct) %>% count(council, maingroup) %>% rename(`Possibly extinct` = n)
    ) %>%
    gather(story, n,-c(council, maingroup)) %>%
    mutate(event_year = event_i)
}

bind_rows(events_emerging) %>% write_csv("outputs/Summary/events_emerging.csv")


p_events <- bind_rows(events_emerging) %>%
  group_by(council, story, event_year) %>%
  summarise(n = sum(n, na.rm = T)) %>%
  filter(story %in% c("Extant", "Updated biodiversity")) %>%
  spread(story, n) %>%
  mutate(change = `Updated biodiversity` - Extant) %>%
  select(-c(Extant, `Updated biodiversity`)) %>%
  left_join(
    dat_events %>%
      rename(council = lga) %>%
      select(council, event_year, event_description)
  ) %>%
  rename(participation = event_description) %>%
  mutate(participation = case_when(str_detect(event_year, "inter") ~ "Interperiod", T ~ "Event")) %>%
  left_join(
    dat_events_unique %>%
      rowwise() %>%
      mutate(period = length(seq(
        dmy(start_date), dmy(end_date), by = "day"
      ))) %>%
      select(event_year, period)
  ) %>%
  ggplot(aes(
    x = factor(
      event_year,
      levels = c(
        "bb_20",
        "inter_bb20_cnc21",
        "cnc_21",
        "inter_cnc21_bb21",
        "bb_21",
        "inter_bb21_gsb21",
        "gsb_21",
        "inter_gsb21_cnc22",
        "cnc_22",
        "inter_cnc22_bb22",
        "bb_22"
      ),
      labels = c(
        "BioBlitz 2020",
        "Interperiod 1",
        "City Nature Challenge 2021",
        "Interperiod 2",
        "BioBlitz 2021",
        "Interperiod 3",
        "Great Southern Bioblitz 2021",
        "Interperiod 4",
        "City Nature Challenge 2022",
        "Interperiod 5",
        "BioBlitz 2022"
      )
    ),
    y = change / period,
    fill = participation
  )) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  labs(x = "Timeline", y = "Rate of biodiversity accumulation (new species per day)") +
  scale_fill_viridis_d(option = "turbo",
                       name = "",
                       direction = -1)

ggsave(
  "outputs/Summary/events_timeline.png",
  plot = p_events,
  width = 2000,
  height = 1500,
  units = "px",
  dpi = 300,
  scale = 1.5
)

bind_rows(events_emerging) %>%
  group_by(council, story, event_year) %>%
  summarise(n = sum(n, na.rm = T)) %>%
  filter(story %in% c("Extant", "Updated biodiversity")) %>%
  spread(story, n) %>%
  mutate(change = `Updated biodiversity` - Extant) %>%
  select(-c(Extant, `Updated biodiversity`)) %>%
  left_join(
    dat_events %>%
      rename(council = lga) %>%
      select(council, event_year, event_description)
  ) %>%
  rename(participation = event_description) %>%
  mutate(participation = case_when(is.na(participation) ~ "no", T ~ "yes")) %>%
  left_join(
    dat_events_unique %>%
      rowwise() %>%
      mutate(period = length(seq(
        dmy(start_date), dmy(end_date), by = "day"
      ))) %>%
      select(event_year, period)
  ) %>%
  group_by(event_year) %>%
  summarise(
    mean = mean(change),
    median = median(change),
    IQR_25 = quantile(change)[2],
    IQR_75 = quantile(change)[4],
    min = min(change),
    max = max(change),
    n = n()
  ) %>%
  write_csv("outputs/Summary/events_timeline_stats.csv")

p_events_group_free <- bind_rows(events_emerging) %>%
  group_by(maingroup, council, story, event_year) %>%
  summarise(n = sum(n, na.rm = T)) %>%
  filter(story %in% c("Extant", "Updated biodiversity")) %>%
  spread(story, n) %>%
  mutate(change = `Updated biodiversity` - Extant) %>%
  select(-c(Extant, `Updated biodiversity`)) %>%
  left_join(
    dat_events %>%
      rename(council = lga) %>%
      select(council, event_year, event_description)
  ) %>%
  rename(participation = event_description) %>%
  mutate(participation = case_when(str_detect(event_year, "inter") ~ "Interperiod", T ~ "Event")) %>%
  left_join(
    dat_events_unique %>%
      rowwise() %>%
      mutate(period = length(seq(
        dmy(start_date), dmy(end_date), by = "day"
      ))) %>%
      select(event_year, period)
  ) %>%
  ggplot(aes(
    x = factor(
      event_year,
      levels = c(
        "bb_20",
        "inter_bb20_cnc21",
        "cnc_21",
        "inter_cnc21_bb21",
        "bb_21",
        "inter_bb21_gsb21",
        "gsb_21",
        "inter_gsb21_cnc22",
        "cnc_22",
        "inter_cnc22_bb22",
        "bb_22"
      ),
      labels = c(
        "BioBlitz 2020",
        "Interperiod 1",
        "City Nature Challenge 2021",
        "Interperiod 2",
        "BioBlitz 2021",
        "Interperiod 3",
        "Great Southern Bioblitz 2021",
        "Interperiod 4",
        "City Nature Challenge 2022",
        "Interperiod 5",
        "BioBlitz 2022"
      )
    ),
    y = change / period,
    fill = participation
  )) +
  geom_boxplot() +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    panel.grid.major.x = element_line(color = c("grey65", "grey90"))
  ) +
  labs(x = "Timeline", y = "Rate of biodiversity accumulation (new species per day)") +
  scale_fill_viridis_d(option = "turbo",
                       name = "",
                       direction = -1) +
  facet_wrap( ~ maingroup, scales = "free_y")

ggsave(
  "outputs/Summary/events_timeline_group_free.png",
  plot = p_events_group_free,
  width = 2000,
  height = 1500,
  units = "px",
  dpi = 300,
  scale = 1.5
)

p_events_group <- bind_rows(events_emerging) %>%
  group_by(maingroup, council, story, event_year) %>%
  summarise(n = sum(n, na.rm = T)) %>%
  filter(story %in% c("Extant", "Updated biodiversity")) %>%
  spread(story, n) %>%
  mutate(change = `Updated biodiversity` - Extant) %>%
  select(-c(Extant, `Updated biodiversity`)) %>%
  left_join(
    dat_events %>%
      rename(council = lga) %>%
      select(council, event_year, event_description)
  ) %>%
  rename(participation = event_description) %>%
  mutate(participation = case_when(str_detect(event_year, "inter") ~ "Interperiod", T ~ "Event")) %>%
  left_join(
    dat_events_unique %>%
      rowwise() %>%
      mutate(period = length(seq(
        dmy(start_date), dmy(end_date), by = "day"
      ))) %>%
      select(event_year, period)
  ) %>%
  ggplot(aes(
    x = factor(
      event_year,
      levels = c(
        "bb_20",
        "inter_bb20_cnc21",
        "cnc_21",
        "inter_cnc21_bb21",
        "bb_21",
        "inter_bb21_gsb21",
        "gsb_21",
        "inter_gsb21_cnc22",
        "cnc_22",
        "inter_cnc22_bb22",
        "bb_22"
      ),
      labels = c(
        "BioBlitz 2020",
        "Interperiod 1",
        "City Nature Challenge 2021",
        "Interperiod 2",
        "BioBlitz 2021",
        "Interperiod 3",
        "Great Southern Bioblitz 2021",
        "Interperiod 4",
        "City Nature Challenge 2022",
        "Interperiod 5",
        "BioBlitz 2022"
      )
    ),
    y = change / period,
    fill = participation
  )) +
  geom_boxplot() +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    panel.grid.major.x = element_line(color = c("grey65", "grey90"))
  ) +
  labs(x = "Timeline", y = "Rate of biodiversity accumulation (new species per day)") +
  scale_fill_viridis_d(option = "turbo",
                       name = "",
                       direction = -1) +
  facet_wrap( ~ maingroup)

ggsave(
  "outputs/Summary/events_timeline_group.png",
  plot = p_events_group,
  width = 2000,
  height = 1500,
  units = "px",
  dpi = 300,
  scale = 1.5
)

p_events_linear <- bind_rows(events_emerging) %>%
  group_by(council, story, event_year) %>%
  summarise(n = sum(n, na.rm = T)) %>%
  filter(story %in% c("Extant", "Updated biodiversity")) %>%
  spread(story, n) %>%
  mutate(change = `Updated biodiversity` - Extant) %>%
  select(-c(Extant, `Updated biodiversity`)) %>%
  left_join(
    dat_events %>%
      rename(council = lga) %>%
      select(council, event_year, event_description)
  ) %>%
  rename(participation = event_description) %>%
  mutate(participation = case_when(is.na(participation) ~ "no", T ~ "yes")) %>%
  left_join(
    dat_events_unique %>%
      rowwise() %>%
      mutate(period = length(seq(
        dmy(start_date), dmy(end_date), by = "day"
      ))) %>%
      select(event_year, period)
  ) %>%
  ggplot(aes(x = period, y = change, colour = participation)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_y_continuous(breaks = c(0, 50, 100, 150, 200, 250, 300, 350)) +
  scale_x_continuous(breaks = c(0, 50, 100, 150, 200, 250)) +
  labs(x = "Number of days", y = "Change in known biodiversity") +
  scale_colour_viridis_d(option = "turbo", name = "During event?")

ggsave(
  "outputs/Summary/events_ancova_linear.png",
  plot = p_events_linear,
  width = 2000,
  height = 1500,
  units = "px",
  dpi = 300,
  scale = 1.5
)

log_to_linear <- function(x) {
  10^x
}

linear_to_log <- function(x) {
  log10(x)
}

trans <- trans_new("log10_to_linear", transform = log_to_linear, inverse = linear_to_log)

p_events_linear_pred <- expand.grid(period = log10(c(1:30)), participation = c("no", "yes")) %>%
  mutate(y = predict(lm(change ~ period * participation, data = dat_events_emerging), ., interval = "confidence")[,1],
         ymin = predict(lm(change ~ period * participation, data = dat_events_emerging), ., interval = "confidence")[,2],
         ymax = predict(lm(change ~ period * participation, data = dat_events_emerging), ., interval = "confidence")[,3]) %>%
  ggplot(aes(x = period, colour = participation, fill = participation)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = .15, colour = NA) +
  geom_line(aes(y = y), size = 1) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = "Number of days", y = "Change in known biodiversity") +
  scale_colour_viridis_d(option = "turbo", name = "During event?") +
  scale_fill_viridis_d(option = "turbo", name = "During event?") +
  scale_x_continuous(trans = trans, breaks = log10(c(1, 10, 20, 30)), labels = c(1, 10, 20, 30)) +
  scale_y_continuous(trans = trans, breaks = log10(c(1, 25, 50, 75, 100)), labels = c(1, 25, 50, 75, 100)) +
  scale_colour_viridis_d(option = "turbo", name = "During event?") +
  scale_fill_viridis_d(option = "turbo", name = "During event?")

ggsave(
  "outputs/Summary/events_ancova_linear_pred.png",
  plot = p_events_linear_pred,
  width = 2000,
  height = 1500,
  units = "px",
  dpi = 300,
  scale = 1.5
)

p_events_log10 <- bind_rows(events_emerging) %>%
  group_by(council, story, event_year) %>%
  summarise(n = sum(n, na.rm = T)) %>%
  filter(story %in% c("Extant", "Updated biodiversity")) %>%
  spread(story, n) %>%
  mutate(change = `Updated biodiversity` - Extant) %>%
  select(-c(Extant, `Updated biodiversity`)) %>%
  left_join(
    dat_events %>%
      rename(council = lga) %>%
      select(council, event_year, event_description)
  ) %>%
  rename(participation = event_description) %>%
  mutate(participation = case_when(is.na(participation) ~ "no", T ~ "yes")) %>%
  left_join(
    dat_events_unique %>%
      rowwise() %>%
      mutate(period = length(seq(
        dmy(start_date), dmy(end_date), by = "day"
      ))) %>%
      select(event_year, period)
  ) %>%
  ggplot(aes(x = period, y = change, colour = participation)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_y_continuous(trans = "log10",
                     breaks = c(0, 10, 50, 100, 150, 200, 250, 300, 350)) +
  scale_x_continuous(trans = "log10",
                     breaks = c(0, 10, 50, 100, 150, 200, 250)) +
  labs(x = "Number of days", y = "Change in known biodiversity") +
  scale_colour_viridis_d(option = "turbo", name = "During event?")

ggsave(
  "outputs/Summary/events_ancova_log10.png",
  plot = p_events_log10,
  width = 2000,
  height = 1500,
  units = "px",
  dpi = 300,
  scale = 1.5
)

p_events_log105 <- bind_rows(events_emerging) %>%
  group_by(council, story, event_year) %>%
  summarise(n = sum(n, na.rm = T)) %>%
  filter(story %in% c("Extant", "Updated biodiversity")) %>%
  spread(story, n) %>%
  mutate(change = `Updated biodiversity` - Extant) %>%
  select(-c(Extant, `Updated biodiversity`)) %>%
  left_join(
    dat_events %>%
      rename(council = lga) %>%
      select(council, event_year, event_description)
  ) %>%
  rename(participation = event_description) %>%
  mutate(participation = case_when(is.na(participation) ~ "no", T ~ "yes")) %>%
  left_join(
    dat_events_unique %>%
      rowwise() %>%
      mutate(period = length(seq(
        dmy(start_date), dmy(end_date), by = "day"
      ))) %>%
      select(event_year, period)
  ) %>%
  ggplot(aes(
    x = log(period, base = 5),
    y = change,
    colour = participation
  )) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_y_continuous(trans = "log10",
                     breaks = c(0, 10, 50, 100, 150, 200, 250, 300, 350)) +
  scale_x_continuous(breaks = log(c(0, 10, 50, 100, 150, 200, 250), base = 5),
                     labels = c(0, 10, 50, 100, 150, 200, 250)) +
  labs(x = "Number of days", y = "Change in known biodiversity") +
  scale_colour_viridis_d(option = "turbo", name = "During event?")

ggsave(
  "outputs/Summary/events_ancova_log105.png",
  plot = p_events_log105,
  width = 2000,
  height = 1500,
  units = "px",
  dpi = 300,
  scale = 1.5
)

dat_events_emerging <- bind_rows(events_emerging) %>%
  group_by(council, story, event_year) %>%
  summarise(n = sum(n)) %>%
  filter(story %in% c("Extant", "Updated biodiversity")) %>%
  spread(story, n) %>% mutate(change = `Updated biodiversity` - Extant) %>%
  select(-c(Extant, `Updated biodiversity`)) %>%
  left_join(
    dat_events %>%
      rename(council = lga) %>%
      select(council, event_year, event_description)
  ) %>%
  rename(participation = event_description) %>%
  mutate(participation = case_when(is.na(participation) ~ "no", T ~ "yes")) %>%
  left_join(
    dat_events_unique %>%
      rowwise() %>%
      mutate(period = length(seq(
        dmy(start_date), dmy(end_date), by = "day"
      ))) %>%
      select(event_year, period)
  ) %>%
  mutate(change = log10(change),
         period = log10(period)) %>%
  ungroup()

broom::tidy(lm(change ~ period * participation, data = dat_events_emerging)) %>%
  mutate_if(is.numeric, round, 3) %>%
  write_csv("outputs/Summary/events_model_coefficients.csv")

broom::glance(lm(change ~ period * participation, data = dat_events_emerging)) %>%
  mutate_if(is.numeric, round, 3) %>%
  write_csv("outputs/Summary/events_model_summary.csv")

# species accumulation curves:
cum_sp <- list()
for (i in bb_palette$maingroup) {
  cum_sp[[i]] <- dat_all_events %>%
    filter(maingroup == i) %>%
    group_by(morphospecies) %>%
    arrange(std_date) %>%
    slice(1L) %>%
    ungroup() %>%
    arrange(std_date) %>%
    group_by(std_date) %>%
    count() %>%
    ungroup(std_date) %>%
    full_join(tibble(std_date = seq(
      min(dat_all_events$std_date),
      max(dat_all_events$std_date),
      by = "day"
    ))) %>%
    arrange(std_date) %>%
    mutate(n = replace_na(n, 0),
           maingroup = i,
           n_sum = cumsum(n))
}
cum_sp <- bind_rows(cum_sp)

scaled_cum_sp <- list()
for (i in bb_palette$maingroup) {
  group_cumsum <- cum_sp %>%
    filter(std_date >= dmy("01012020")) %>%
    filter(maingroup == i) %>%
    pull(n_sum)
  
  scaled_cum_sp[[i]] <- cum_sp %>%
    filter(std_date >= dmy("01012020")) %>%
    filter(maingroup == i) %>%
    mutate(n_sum = group_cumsum) %>%
    ungroup() %>%
    add_row(
      maingroup = i,
      std_date = dmy("01012020"),
      n = 1,
      n_sum = min(group_cumsum)
    ) %>%
    add_row(
      maingroup = i,
      std_date = dmy("30102022"),
      n = 1,
      n_sum = max(group_cumsum)
    )
}

p_cumsum <- cum_sp %>%
  ungroup() %>%
  group_by(std_date) %>%
  summarise(n = sum(n)) %>%
  mutate(n_sum = cumsum(n)) %>%
  filter(std_date >= dmy("01012020")) %>%
  ggplot() +
  geom_rect(
    aes(
      xmin = start_date,
      xmax = end_date,
      ymin = min(
        cum_sp %>%
          ungroup() %>%
          group_by(std_date) %>%
          summarise(n = sum(n)) %>%
          mutate(n_sum = cumsum(n)) %>%
          filter(std_date >= dmy("01012020")) %>%
          pull(n_sum)
      ),
      ymax = max(
        cum_sp %>%
          ungroup() %>%
          group_by(std_date) %>%
          summarise(n = sum(n)) %>%
          mutate(n_sum = cumsum(n)) %>%
          filter(std_date >= dmy("01012020")) %>%
          pull(n_sum)
      )
    ),
    data = dat_events %>%
      select(event_year, start_date, end_date) %>% unique() %>%
      filter(!str_detect(event_year, "inter")) %>%
      mutate(
        start_date = dmy(start_date),
        end_date = dmy(end_date)
      ) %>%
      filter(start_date < dmy("01102022")),
    fill = "darkred",
    alpha = .4
  ) +
  geom_text(
    aes(
      x = x,
      y = y,
      label = str_to_upper(event_year)
    ),
    data = dat_events %>%
      select(event_year, start_date, end_date) %>% unique() %>%
      filter(!str_detect(event_year, "inter")) %>%
      mutate(
        start_date = dmy(start_date),
        end_date = dmy(end_date)
      ) %>%
      filter(start_date < dmy("01102022")) %>%
      mutate(
        x = start_date + (end_date - start_date) / 2,
        y = c(15680, 15680, 15680, 15680, 15680, 15680)
      ),
    angle = 45,
    hjust = 0.5
  ) +
  geom_line(aes(x = std_date, y = n_sum), linewidth = 1.5) +
  theme_bw() +
  labs(y = "Cumulative number of species",
       x = "Date")

ggsave(
  "outputs/Summary/sp_cumsum.png",
  plot = p_cumsum,
  width = 2000,
  height = 1500,
  units = "px",
  dpi = 300,
  scale = 1.5
)

p_cumsum_group <- cum_sp %>%
  filter(std_date >= dmy("01012020")) %>%
  ggplot() +
  geom_rect(
    aes(
      xmin = start_date,
      xmax = end_date,
      ymin = min_sum,
      ymax = max_sum
    ),
    data = expand_grid(
      left_join(
        bind_rows(scaled_cum_sp) %>%
          group_by(maingroup) %>%
          arrange(desc(std_date)) %>%
          slice(1L) %>%
          select(-c(n, std_date)) %>%
          rename(max_sum = n_sum),
        bind_rows(scaled_cum_sp) %>%
          group_by(maingroup) %>%
          arrange(std_date) %>%
          slice(1L) %>%
          select(-c(n, std_date)) %>%
          rename(min_sum = n_sum)
      ),
      dat_events %>%
        select(event_year, start_date, end_date) %>% unique() %>%
        filter(!str_detect(event_year, "inter")) %>%
        mutate(
          start_date = dmy(start_date),
          end_date = dmy(end_date)
        ) %>%
        filter(start_date < dmy("01102022"))
    ),
    fill = "darkred",
    alpha = .4
  ) +
  geom_line(
    aes(x = std_date, y = n_sum, colour = maingroup),
    data = bind_rows(scaled_cum_sp),
    linewidth = .8
  ) +
  theme_bw() +
  labs(y = "Cumulative number of species",
       x = "Date",
       colour = "Group") +
  scale_colour_manual(values = bb_palette$col) +
  facet_wrap( ~ maingroup, scales = "free_y")

ggsave(
  "outputs/Summary/sp_cumsum_group.png",
  plot = p_cumsum_group,
  width = 2000,
  height = 1500,
  units = "px",
  dpi = 300,
  scale = 1.5
)

left_join(
  cum_sp %>%
    left_join(
      dat_events %>%
        select(event_year, start_date, end_date) %>% unique() %>%
        filter(!str_detect(event_year, "inter")) %>%
        mutate(std_date = dmy(start_date)) %>%
        select(event_year, std_date)
    ) %>%
    drop_na(event_year) %>%
    mutate(start_sp = n_sum) %>%
    select(maingroup, event_year, start_sp),
  cum_sp %>%
    left_join(
      dat_events %>%
        select(event_year, start_date, end_date) %>% unique() %>%
        filter(!str_detect(event_year, "inter")) %>%
        mutate(std_date = dmy(end_date)) %>%
        select(event_year, std_date)
    ) %>%
    drop_na(event_year) %>%
    mutate(end_sp = n_sum) %>%
    select(maingroup, event_year, end_sp)
) %>%
  write_csv("outputs/Summary/sp_events.csv")

### MADDIE'S RESEARCH QUESTIONS ###
Records_PlanZone_Events <-
  read_csv("data/Records_PlanZone_Events.csv")

# zone_chisq <- Records_PlanZone_Events %>%
#   group_by(council, zone) %>%
#   summarise(records_count = sum(records_count)) %>%
#   pivot_wider(names_from = council,
#               values_from = records_count) %>%
#   mutate_if(is.numeric, .funs = function(x) replace_na(x, 0)) %>%
#   ungroup() %>%
#   select(-zone) %>%
#   rowSums() %>%
#   chisq.test()
#
# zone_chisq

stars <- function(p) {
  case_when(p <= 0.001 ~ "***",
            p <= 0.01  ~ "**",
            p <= 0.05  ~ "*",
            TRUE       ~ "ns")
}

zone_anova <- Records_PlanZone_Events %>%
  group_by(council, record_period) %>%
  mutate(records_sum = sum(records_count)) %>%
  ungroup() %>%
  mutate(records_prop = records_count / records_sum) %>%
  lm(records_prop ~ zone, data = .)

zone_anova %>%
  broom::glance() %>%
  mutate_if(is.numeric, round, 3) %>%
  write_csv("outputs/Summary/zone_model_summary.csv")

TukeyHSD(aov(zone_anova), conf.level = .95) %>%
  broom::tidy() %>%
  mutate_if(is.numeric, round, 3) %>%
  write_csv("outputs/Summary/zone_model_comparisons.csv")

p_zone <- Records_PlanZone_Events %>%
  group_by(council, record_period) %>%
  mutate(records_sum = sum(records_count)) %>%
  ungroup() %>%
  mutate(records_prop = records_count / records_sum) %>%
  ggplot(aes(x = zone, y = records_prop, fill = zone)) +
  geom_boxplot() +
  theme_bw() +
  scale_fill_viridis_d() +
  labs(x = "Landuse Zone",
       y = "Proportion of records",
       fill = "Landuse Zone")

ggsave(
  "outputs/Summary/zone_boxplot.png",
  plot = p_zone,
  width = 2000,
  height = 1500,
  units = "px",
  dpi = 300,
  scale = 1.5
)

Records_PlanZone_Events %>%
  group_by(council, record_period) %>%
  mutate(records_sum = sum(records_count)) %>%
  ungroup() %>%
  mutate(records_prop = records_count / records_sum) %>%
  group_by(zone) %>%
  summarise(
    mean = mean(records_prop),
    median = median(records_prop),
    IQR_25 = quantile(records_prop)[2],
    IQR_75 = quantile(records_prop)[4],
    min = min(records_prop),
    max = max(records_prop),
    n = n()
  ) %>%
  write_csv("outputs/Summary/zone_boxplot_stats.csv")

p_zone2 <- Records_PlanZone_Events %>%
  group_by(council, zone) %>%
  summarise(records_count = sum(records_count)) %>%
  ggplot(aes(x = council, y = records_count, fill = zone)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_viridis_d() +
  labs(x = "LGA",
       y = "Proportion of records",
       fill = "Landuse Zone")

ggsave(
  "outputs/Summary/zone_barplot.png",
  plot = p_zone2,
  width = 2000,
  height = 1500,
  units = "px",
  dpi = 300,
  scale = 1.5
)

Records_GS <- read_csv("data/Records_GS.csv") %>%
  mutate(Remnant_veg_area_ha = as.numeric(units::set_units(units::set_units(Remnant_veg_area_m2, "m2"), "ha"))) %>%
  select(-Remnant_veg_area_m2)

gs_glm <-
  glm(
    Num_record ~ log(GS_Area_ha) + GS_ACCESS + GS_category + log(Surrounding_Population) + Remnant_veg_area_ha,
    data = Records_GS %>%
      filter(!GS_category == "Civic squares and promenades",
             GS_Area_ha > 0) %>%
      mutate(Remnant_veg_area_ha = Remnant_veg_area_ha /
               GS_Area_ha),
    family = "poisson"
  )

gs_glm %>%
  broom::tidy() %>%
  mutate_if(is.numeric, round, 3) %>%
  write_csv("outputs/Summary/gs_model_coefficients.csv")

gs_glm %>%
  broom::glance() %>%
  mutate(deviance.explained = with(summary(gs_glm), 1 - deviance / null.deviance)) %>%
  mutate_if(is.numeric, round, 3) %>%
  write_csv("outputs/Summary/gs_model_summary.csv")

p_area <-
  sjPlot::plot_model(
    gs_glm,
    type = "eff",
    terms = c("GS_Area_ha", "GS_ACCESS", "GS_category")
  ) +
  theme_bw() +
  scale_colour_viridis_d(option = "cividis") +
  scale_fill_viridis_d(option = "cividis") +
  labs(
    x = "Area of greenspace (ha)",
    y = "Predicted number of records",
    title = "",
    fill = "Greenspace access",
    colour = "Greenspace access"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

ggsave(
  "outputs/Summary/gs_area.png",
  plot = p_area,
  width = 2000,
  height = 1500,
  units = "px",
  dpi = 300,
  scale = 1.5
)

# Extract the data for the fourth facet (Natural and semi-natural open space)
p_area_natural <- ggplot_build(p_area)

p_area_natural <- subset(p_area_natural$data[[2]], PANEL == 4) %>%
  ggplot(aes(x = x)) +
  geom_ribbon(aes(ymax = ymax, ymin = ymin, fill = fill), alpha = .15) +
  geom_line(aes(y = y, colour = colour), size = 1, data = subset(p_area_natural$data[[1]], PANEL == 4)) +
  theme_bw() +
  scale_colour_viridis_d(option = "cividis", labels = c("Closed", "Open")) +
  scale_fill_viridis_d(option = "cividis", labels = c("Closed", "Open")) +
  labs(
    x = "Area of greenspace (ha)",
    y = "Predicted number of records",
    title = "Natural and semi-natural open space",
    fill = "Greenspace access",
    colour = "Greenspace access"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

ggsave(
  "outputs/Summary/gs_area_natural.png",
  plot = p_area_natural,
  width = 2000,
  height = 1500,
  units = "px",
  dpi = 300,
  scale = 1.5
)

p_pop <-
  sjPlot::plot_model(
    gs_glm,
    type = "eff",
    terms = c("Surrounding_Population", "GS_ACCESS", "GS_category")
  ) +
  theme_bw() +
  scale_colour_viridis_d(option = "cividis") +
  scale_fill_viridis_d(option = "cividis") +
  labs(
    x = "2021 Population of Statistical Area 2 greenspace falls into",
    y = "Predicted number of records",
    title = "",
    fill = "Greenspace access",
    colour = "Greenspace access"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

ggsave(
  "outputs/Summary/gs_pop.png",
  plot = p_pop,
  width = 2000,
  height = 1500,
  units = "px",
  dpi = 300,
  scale = 1.5
)

# Extract the data for the fourth facet (Natural and semi-natural open space)
p_pop_natural <- ggplot_build(p_pop)

p_pop_natural <- subset(p_pop_natural$data[[2]], PANEL == 4) %>%
  ggplot(aes(x = x)) +
  geom_ribbon(aes(ymax = ymax, ymin = ymin, fill = fill), alpha = .15) +
  geom_line(aes(y = y, colour = colour), size = 1, data = subset(p_pop_natural$data[[1]], PANEL == 4)) +
  theme_bw() +
  scale_colour_viridis_d(option = "cividis", labels = c("Closed", "Open")) +
  scale_fill_viridis_d(option = "cividis", labels = c("Closed", "Open")) +
  labs(
    x = "2021 Population of Statistical Area 2 greenspace falls into",
    y = "Predicted number of records",
    title = "Natural and semi-natural open space",
    fill = "Greenspace access",
    colour = "Greenspace access"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

ggsave(
  "outputs/Summary/gs_pop_natural.png",
  plot = p_pop_natural,
  width = 2000,
  height = 1500,
  units = "px",
  dpi = 300,
  scale = 1.5
)

p_veg <-
  sjPlot::plot_model(
    gs_glm,
    type = "eff",
    terms = c("Remnant_veg_area_ha", "GS_ACCESS", "GS_category")
  ) +
  theme_bw() +
  scale_colour_viridis_d(option = "cividis") +
  scale_fill_viridis_d(option = "cividis") +
  scale_x_continuous(limits = c(0, 1)) +
  labs(
    x = "Proportion of remnant vegetation within greenspace",
    y = "Predicted number of records",
    title = "",
    fill = "Greenspace access",
    colour = "Greenspace access"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

ggsave(
  "outputs/Summary/gs_veg.png",
  plot = p_veg,
  width = 2000,
  height = 1500,
  units = "px",
  dpi = 300,
  scale = 1.5
)

# Extract the data for the fourth facet (Natural and semi-natural open space)
p_veg_natural <- ggplot_build(p_veg)

p_veg_natural <- subset(p_veg_natural$data[[2]], PANEL == 4) %>%
  ggplot(aes(x = x)) +
  geom_ribbon(aes(ymax = ymax, ymin = ymin, fill = fill), alpha = .15) +
  geom_line(aes(y = y, colour = colour), size = 1, data = subset(p_veg_natural$data[[1]], PANEL == 4)) +
  theme_bw() +
  scale_colour_viridis_d(option = "cividis", labels = c("Closed", "Open")) +
  scale_fill_viridis_d(option = "cividis", labels = c("Closed", "Open")) +
  labs(
    x = "Proportion of remnant vegetation within greenspace",
    y = "Predicted number of records",
    title = "Natural and semi-natural open space",
    fill = "Greenspace access",
    colour = "Greenspace access"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

ggsave(
  "outputs/Summary/gs_veg_natural.png",
  plot = p_veg_natural,
  width = 2000,
  height = 1500,
  units = "px",
  dpi = 300,
  scale = 1.5
)

p_cat <-
  sjPlot::plot_model(gs_glm,
                     type = "eff",
                     terms = c("GS_category", "GS_ACCESS")) +
  theme_bw() +
  scale_colour_viridis_d(option = "cividis") +
  labs(
    colour = "Greenspace access",
    y = "Predicted number of records",
    title = "",
    x = ""
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  "outputs/Summary/gs_cat.png",
  plot = p_cat,
  width = 2000,
  height = 1500,
  units = "px",
  dpi = 300,
  scale = 1.5
)
