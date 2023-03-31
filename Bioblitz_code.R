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
source("R/bb_functions.R")

## LOAD DATA
# Load the clean dataset of biodiversity records, and the cypher file 
dat_all <- read_csv("data/Records v1 7Mar23.csv")

cypher <- read_csv("data/Morphospecies cypher v1 7Mar23.csv")

synonyms <- read_csv("data/Synonyms v1 7Mar23.csv")

dat_all <- dat_all %>%
  left_join(synonyms, by = c("morphospecies" = "synonym")) %>%
  mutate(morphospecies = case_when(is.na(accepted) ~ morphospecies,
                                   T ~ accepted),
         std_date = parse_date_time(std_date, c("dmy", "ymd")),
         record_period = case_when(year(std_date) == 2022 & month(std_date) == 9 ~ "bb22",
                                   year(std_date) < 1992 ~ "historical",
                                   T ~ "recent"),
         occ = 1,
         year = year(std_date)) %>%
  dplyr::select(-accepted) %>%
  left_join(cypher) %>%
  drop_na(maingroup) %>%
  filter(!maingroup %in% c("Lampreys", "Protozoan")) %>%
  mutate(std_date = as.character(std_date),
         X = row_number())

councils_full <- unique(dat_all$council)

# generate discrete colour palette
bb_palette <- tibble(maingroup = sort(unique(dat_all$maingroup)),
                     col = c("#8C468C",
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
                             "#F7861C"))

# generate gradient colour palettes
taxa <- c("Arachnids", "Birds", "Cartilaginous fishes", "Chromists", "Frogs", "Fungi", "Insects", "Mammals", "Other invertebrates", "Plants", "Ray-finned fishes", "Reptiles")
palettes <- list(Arachnids = c("#ddcedd", "#c9acc9", "#b58ab5", "#a168a0", "#8c468c"),
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
                 Reptiles = c("#fbebdc", "#fdd3ac", "#fdba7d", "#fba150", "#f7861c"))

# load map of Victoria and reproject to Lambert equal-area
map_VIC <- read_sf("SpatialData/AD_LGA_AREA_POLYGON.shp") %>%
  filter(STATE == "VIC")

# filter based on LGAs
map_LGA <- map_VIC %>%
  st_transform(crs = 3112) %>%
  filter(NAME %in% str_to_upper(councils_full)) %>%
  rename(council = NAME) %>%
  mutate(council = case_when(council == "BAYSIDE" ~ "BSD",
                             council == "CARDINIA" ~ "CAR",
                             council == "CASEY" ~ "CAS",
                             council == "FRANKSTON" ~ "FRA",
                             council == "GREATER DANDENONG" ~ "GDA",
                             council == "KINGSTON" ~ "KIN",
                             council == "KNOX" ~ "KNO",
                             council == "MONASH" ~ "MON",
                             council == "MORNINGTON PENINSULA" ~ "MPE",
                             council == "YARRA RANGES" ~  "YRA"))

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
  left_join(st_as_sf(., coords = c("decimalLongitude", "decimalLatitude"), crs = "WGS84") %>%
              st_transform(crs = 3112) %>%
              st_join(map_LGA) %>%
              as_tibble() %>%
              select(X, council), by = "X") %>%
  filter(!is.na(council))

# convert data to grid form
occ_GRID <- dat_all %>%
  filter(record_period %in% c("recent", "bb22")) %>%
  select(decimalLatitude, decimalLongitude, recordID, record_period, maingroup) %>%
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = st_crs(map_VIC)) %>%
  st_transform(crs = 3112) %>%
  st_intersection(grid_LGA) %>%
  as_tibble() %>%
  select(-geometry) %>%
  group_by(FID, maingroup, record_period) %>%
  summarise(occ = n())

# load map of green spaces and reproject
greenspaces <- read_sf("data/VPA_Draft_Open_Space_Data.shp") %>%
  filter(OS_TYPE %in% c("Private open space", "Public open space")) %>%
  mutate(council = case_when(LGA == "BAYSIDE" ~ "BSD",
                             LGA == "CARDINIA" ~ "CAR",
                             LGA == "CASEY" ~ "CAS",
                             LGA == "FRANKSTON" ~ "FRA",
                             LGA == "GREATER DANDENONG" ~ "GDA",
                             LGA == "KINGSTON" ~ "KIN",
                             LGA == "KNOX" ~ "KNO",
                             LGA == "MONASH" ~ "MON",
                             LGA == "MORNINGTON" ~ "MPE",
                             LGA == "YARRA RANGES" ~  "YRA")) %>%
  filter(!is.na(council)) %>%
  st_transform(st_crs(map_LGA))

# counts per green space
occ_greenspace <- dat_all %>%
  filter(record_period %in% c("recent", "bb22")) %>%
  select(decimalLatitude, decimalLongitude, recordID, record_period, maingroup) %>%
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = st_crs(map_VIC)) %>%
  st_transform(crs = 3112) %>%
  st_intersection(greenspaces) %>%
  as_tibble() %>%
  select(-geometry) %>%
  group_by(FID, maingroup, record_period) %>%
  summarise(occ = n())

sp_greenspace <- dat_all %>%
  filter(record_period %in% c("recent", "bb22")) %>%
  select(decimalLatitude, decimalLongitude, recordID, record_period, maingroup, morphospecies) %>%
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = st_crs(map_VIC)) %>%
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

for (i in councils){
  if(i == "ALL") {
    dat <- dat_all
    grid_dat <- grid_LGA
    green_dat <- greenspaces
    LGA_dat <- map_LGA
  } else {
    dat <- dat_all %>% filter(council == i)
    grid_dat <- grid_LGA %>%
      filter(st_intersects(st_geometry(.), map_LGA %>% filter(council == i), sparse = F)[,1])
    green_dat <- greenspaces %>% filter(council == i)
    LGA_dat <- map_LGA %>% filter(council == i)
  }
  
  # LEGEND MAP
  p_legend <- map_VIC %>%
    ggplot() +
    geom_sf(fill = "grey95", colour = "grey75") +
    geom_sf(data = greenspaces, fill = "darkolivegreen", colour = NA, alpha = .8) +
    geom_sf(data = LGA_dat, fill = NA, colour = "black", size = 2) +
    theme_bw() +
    coord_sf(crs = 3112,
             xlim = st_bbox(LGA_dat)[c(1,3)],
             ylim = st_bbox(LGA_dat)[c(2,4)]) +
    theme(panel.background = element_rect(fill = 'lightblue'),
          panel.grid = element_blank())
  
  ggsave(paste0("outputs/Maps/legend_map_", i,".png"),
         plot = p_legend,
         width = 2000,
         height = 1500, 
         units = "px",
         dpi = 300)
  
  # TABLE 1: species list, Nrecords total/period, FFG/EPBC, origin, common name, etc...
  Splist.Nrecords(dat,
                  full_dat = dat_all,
                  cypher = cypher,
                  path = "outputs/Table 1",
                  i)
  
  ## TABLE 2: emerging stories [lazarus, new, new_introduced]
  emerging(file = paste0("outputs/Table 1/Splist_Nrecords_", i, ".csv"),
           path = "outputs/Table 2",
           i)
  
  ## TABLE 3: species richness by main group, by period, unique species
  SpRichness(file = paste0("outputs/Table 1/Splist_Nrecords_", i, ".csv"),
             path = "outputs/Table 3",
             group = maingroup,
             i)
  
  ## TABLE 4: species richness by subgroup, by period, unique species
  SpRichness(file = paste0("outputs/Table 1/Splist_Nrecords_", i, ".csv"),
             path = "outputs/Table 4",
             group = subgroup,
             i)
  
  ## TABLE 7: most common species from BB22
  Common.sp(file = paste0("outputs/Table 1/Splist_Nrecords_", i, ".csv"),
            path = "outputs/Table 7",
            i)
  
  ## FIGURE 1
  Common.sp.p(file = paste0("outputs/Table 7/Common_species_", i, ".csv"),
              path = "outputs/Figure 1",
              i,
              palette = bb_palette)
  
  ## TABLE 8: most common group from BB22
  Common.group(file = paste0("outputs/Table 3/Species_richness_by_maingroup_", i, ".csv"),
               path = "outputs/Table 8",
               group = maingroup,
               i)
  
  ## FIGURE 2
  Common.group.p(file = paste0("outputs/Table 8/Common_maingroup_", i, ".csv"),
                 path = "outputs/Figure 2",
                 i,
                 group = maingroup,
                 count = Nrecords_bioblitz,
                 palette = bb_palette,
                 yaxis = T)
  
  Common.group.p(file = paste0("outputs/Table 8/Common_maingroup_", i, ".csv"),
                 path = "outputs/Figure 2",
                 i,
                 group = maingroup,
                 count = Nuniquesp_bioblitz,
                 palette = bb_palette,
                 yaxis = T)
  
  ## TABLE 9: most common subgroups from BB22
  Common.group(file = paste0("outputs/Table 4/Species_richness_by_subgroup_", i, ".csv"),
               path = "outputs/Table 9",
               group = subgroup,
               i)
  
  ## FIGURE 3
  Common.group.p(file = paste0("outputs/Table 9/Common_subgroup_", i, ".csv"),
                 path = "outputs/Figure 3",
                 i,
                 group = subgroup,
                 count = Nrecords_bioblitz,
                 palette = bb_palette,
                 yaxis = T)
  
  Common.group.p(file = paste0("outputs/Table 9/Common_subgroup_", i, ".csv"),
                 path = "outputs/Figure 3",
                 i,
                 group = subgroup,
                 count = Nuniquesp_bioblitz,
                 palette = bb_palette,
                 yaxis = T)
  
  ## TABLE 10: BB22 summaries by greenspace
  Nrecords.GS(dat = green_dat,
              occ_dat = occ_greenspace,
              sp = sp_greenspace,
              path = "outputs/Table 10",
              i)
  
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
  Trend.sp(file = paste0("outputs/Table 1/Splist_Nrecords_", i, ".csv"),
           path = "outputs/Figure X",
           i,
           palette = bb_palette)
  
  ## FIGURE X2: no. observations
  Trend.obs(dat,
            path = "outputs/Figure X2",
            i,
            palette = bb_palette)
  
  ## MAPS
  Map.p(dat,
        grid_dat,
        occ_dat = occ_GRID,
        basemap = map_LGA,
        path = "outputs/Maps/Total",
        i,
        period = "recent")
  
  Map.p(dat,
        grid_dat,
        occ_dat = occ_GRID,
        basemap = map_LGA,
        path = "outputs/Maps/Total",
        i,
        period = "recent",
        type = "Ndensity")
  
  Map.p(dat,
        grid_dat,
        occ_dat = occ_GRID,
        basemap = map_LGA,
        path = "outputs/Maps/Total",
        i,
        period = "bb22")
  
  Map.p(dat,
        grid_dat,
        occ_dat = occ_GRID,
        basemap = map_LGA,
        path = "outputs/Maps/Total",
        i,
        period = "bb22",
        type = "Ndensity")
  
  Map.gs(dat,
         green_dat,
         occ_dat = occ_greenspace,
         basemap = map_LGA,
         path = "outputs/Maps/Total",
         i,
         period = "bb22")
  
  Map.gs(dat,
         green_dat,
         occ_dat = occ_greenspace,
         basemap = map_LGA,
         path = "outputs/Maps/Total",
         i,
         period = "bb22",
         type = "Greenspace")
  
  # Produce maps per maingroup
  for (k in taxa){
    Map.p(dat %>% filter(maingroup == k),
          grid_dat,
          occ_dat = occ_GRID %>% filter(maingroup == k),
          basemap = map_LGA,
          path = paste0("outputs/Maps/", k),
          i,
          period = "recent",
          palette = palettes[[k]])
    
    Map.p(dat %>% filter(maingroup == k),
          grid_dat,
          occ_dat = occ_GRID %>% filter(maingroup == k),
          basemap = map_LGA,
          path = paste0("outputs/Maps/", k),
          i,
          period = "recent",
          type = "Ndensity",
          palette = palettes[[k]])
    
    Map.p(dat %>% filter(maingroup == k),
          grid_dat,
          occ_dat = occ_GRID %>% filter(maingroup == k),
          basemap = map_LGA,
          path = paste0("outputs/Maps/", k),
          i,
          period = "bb22",
          palette = palettes[[k]])
    
    Map.p(dat %>% filter(maingroup == k),
          grid_dat,
          occ_dat = occ_GRID %>% filter(maingroup == k),
          basemap = map_LGA,
          path = paste0("outputs/Maps/", k),
          i,
          period = "bb22",
          type = "Ndensity",
          palette = palettes[[k]])
    
    Map.gs(dat %>% filter(maingroup == k),
           green_dat,
           occ_dat = occ_greenspace %>% filter(maingroup == k),
           basemap = map_LGA,
           path = paste0("outputs/Maps/", k),
           i,
           period = "bb22",
           palette = palettes[[k]])
    
    Map.gs(dat %>% filter(maingroup == k),
           green_dat,
           occ_dat = occ_greenspace %>% filter(maingroup == k),
           basemap = map_LGA,
           path = paste0("outputs/Maps/", k),
           i,
           period = "bb22",
           type = "Greenspace",
           palette = palettes[[k]])
  }
}

## General outputs
Common.group_all <- lapply(dir("outputs/Table 8", full.names = T), read_csv)

Common.group_Nrecords <- bind_rows(lapply(1:length(councils),
                                          function(x)
                                            Common.group_all[[x]] %>%
                                            select(maingroup, Nrecords_bioblitz) %>%
                                            mutate(council = councils[[x]])))  %>%
  spread(maingroup, Nrecords_bioblitz) %>%
  mutate_if(is.numeric, funs(replace_na(., 0))) %>%
  left_join(LGA_coords %>% bind_rows(tibble(X = min(LGA_coords$X-5000),
                                            Y = max(LGA_coords$Y+5000),
                                            council = "ALL")))

Common.group_Nsp <- bind_rows(lapply(1:length(councils),
                                     function(x)
                                       Common.group_all[[x]] %>%
                                       select(maingroup, Nsp_bioblitz) %>%
                                       mutate(council = councils[[x]])))  %>%
  spread(maingroup, Nsp_bioblitz) %>%
  mutate_if(is.numeric, funs(replace_na(., 0))) %>%
  left_join(LGA_coords %>% bind_rows(tibble(X = min(LGA_coords$X-5000),
                                            Y = max(LGA_coords$Y+5000),
                                            council = "ALL")))

no_zero_groups <- colSums(Common.group_Nrecords[, bb_palette$maingroup])
no_zero_groups <- no_zero_groups[which(no_zero_groups > 0)]

plot_pie_map(p_map = map_LGA %>%
               ggplot() +
               geom_sf(fill = "white", colour = "dark grey", inherit.aes = F) +
               geom_scatterpie(aes(x = X, y = Y, group = council, r = r/2),
                               data = Common.group_Nrecords %>%
                                 mutate(r = rowSums(Common.group_Nrecords[, bb_palette$maingroup])),
                               cols = bb_palette$maingroup, color = NA) +
               theme_bw() +
               scale_fill_manual(values = as.character(bb_palette %>%
                                                         filter(maingroup %in% names(no_zero_groups)) %>%
                                                         select(col) %>%
                                                         as_vector()),
                                 name = "Group") +
               geom_scatterpie_legend(rowSums(Common.group_Nrecords[, bb_palette$maingroup])/3,
                                      x = max(LGA_coords$X)+20000,
                                      y = min(LGA_coords$Y)-5000,
                                      labeller = function(x) x*2) +
               labs(title = "Bioblitz 2022 No. observations",
                    x = "Longitude",
                    y = "Latitude"),
             dat = Common.group_Nrecords,
             palette = bb_palette,
             scale.fact = .5)

ggsave(paste0("outputs/Summary/Common_groups_Nrecords.jpg"), 
       width = 2000,
       height = 1500, 
       units = "px",
       dpi = 300)

plot_pie_map(p_map = map_LGA %>%
               ggplot() +
               geom_sf(fill = "white", colour = "dark grey", inherit.aes = F) +
               geom_scatterpie(aes(x = X, y = Y, group = council, r = r*4),
                               data = Common.group_Nsp %>%
                                 mutate(r = rowSums(Common.group_Nsp[, bb_palette$maingroup])),
                               cols = bb_palette$maingroup, color = NA) +
               theme_bw() +
               scale_fill_manual(values = as.character(bb_palette %>%
                                                         filter(maingroup %in% names(no_zero_groups)) %>%
                                                         select(col) %>%
                                                         as_vector()),
                                 name = "Group") +
               geom_scatterpie_legend(rowSums(Common.group_Nsp[, bb_palette$maingroup])*4,
                                      x = max(LGA_coords$X)+20000,
                                      y = min(LGA_coords$Y)-5000,
                                      labeller = function(x) round(x/4, - 1)) +
               labs(title = "Bioblitz 2022 No. species",
                    x = "Longitude",
                    y = "Latitude"),
             dat = Common.group_Nsp,
             bb_palette,
             scale.fact = 4)

ggsave(paste0("outputs/Summary/Common_groups_Nsp.jpg"), 
       width = 2000,
       height = 1500, 
       units = "px",
       dpi = 300)

groups_no_ones <- dat_all %>% filter(record_period == "bb22") %>% group_by(maingroup, morphospecies) %>% count() %>% ungroup() %>% count(maingroup) %>% filter(n > 1) %>% pull(maingroup)

dat_all %>% filter(record_period == "bb22") %>% group_by(maingroup, morphospecies) %>% count() %>% filter(maingroup %in% groups_no_ones) %>%
  ggplot(aes(y = reorder(maingroup, n, FUN = median, decreasing = T), x = n, fill = maingroup, colour = maingroup)) +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2, alpha = .7,
                      jittered_points = TRUE,
                      position = position_points_jitter(width = 0.05, height = 0),
                      point_shape = '|', point_size = 3, point_alpha = 1) +
  scale_x_continuous(trans = "log10") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = as.character(bb_palette %>%
                                            filter(maingroup %in% groups_no_ones) %>%
                                            select(col) %>%
                                            as_vector()),
                    name = "Group") +
  scale_colour_manual(values = as.character(bb_palette %>%
                                              filter(maingroup %in% groups_no_ones) %>%
                                              select(col) %>%
                                              as_vector()),
                      name = "Group") +
  labs(title = "No. observations per species",
       x = "Observations",
       y = "Group")

ggsave(paste0("outputs/Summary/Nrecords_per_sp.jpg"), 
       width = 2000,
       height = 1500, 
       units = "px",
       dpi = 300)

dat_all %>%
  filter(!record_period == "historical") %>%
  group_by(year, maingroup) %>%
  count(morphospecies) %>%
  summarise(upper = quantile(n, probs = .95),
            lower = quantile(n, probs = .05),
            n = mean(n)) %>%
  ggplot(aes(x = year, colour = maingroup, fill = maingroup)) +
  geom_line(aes(y = n)) +
  geom_ribbon(aes(ymax = upper, ymin = lower), colour = NA, alpha = .25) +
  scale_colour_manual(values = bb_palette$col, name = "Group") +
  scale_fill_manual(values = bb_palette$col, name = "Group") +
  theme_bw() +
  facet_wrap(vars(maingroup), scales = "free_y", ncol = 3) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Year",
       y = "Mean no. obersvations per species")

ggsave(paste0("outputs/Summary/Nrecords_per_sp_per_year.jpg"), 
       width = 2000,
       height = 1500, 
       units = "px",
       dpi = 300)

dat_all %>%
  filter(!record_period == "historical") %>%
  select(year, maingroup, morphospecies) %>%
  unique() %>%
  group_by(year, maingroup) %>%
  count() %>%
  summarise(upper = quantile(n, probs = .95),
            lower = quantile(n, probs = .05),
            n = mean(n)) %>%
  ggplot(aes(x = year, y = n, colour = maingroup)) +
  geom_point() +
  geom_smooth(se = F) +
  scale_colour_manual(values = bb_palette$col, name = "Group") +
  theme_bw() +
  facet_wrap(vars(maingroup), scales = "free_y", ncol = 3) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Year",
       y = "No. species observed per year")

ggsave(paste0("outputs/Summary/Nsp_per_year.jpg"), 
       width = 2000,
       height = 1500, 
       units = "px",
       dpi = 300)

dat_2022 <- dat_all %>%
  filter(year == 2022) %>%
  mutate(eventDate = ymd(std_date))

dat_2022 %>%
  group_by(eventDate, maingroup, record_period) %>%
  count() %>%
  ggplot(aes(x = eventDate, y = n, colour = maingroup, fill = maingroup, linetype = factor(record_period, labels = c("Bioblitz", "Baseline")))) +
  geom_point(alpha = .1) +
  geom_smooth() +
  scale_colour_manual(values = bb_palette$col) +
  scale_fill_manual(values = bb_palette$col) +
  theme_bw() +
  facet_wrap(vars(maingroup), scales = "free_y", ncol = 3) +
  labs(y = "No. observations per day",
       x = "Date",
       colour = "Group",
       fill = "Group",
       linetype = "Record period")

ggsave(paste0("outputs/Summary/Nrecords_2022.jpg"), 
       width = 2000,
       height = 1500, 
       units = "px",
       dpi = 300)

dat_2022 %>%
  group_by(eventDate, maingroup, record_period, morphospecies) %>%
  count() %>%
  select(-n) %>%
  group_by(eventDate, maingroup, record_period) %>%
  count() %>%
  ggplot(aes(x = eventDate, y = n, colour = maingroup, fill = maingroup, linetype = factor(record_period, labels = c("Bioblitz", "Baseline")))) +
  geom_point(alpha = .1) +
  geom_smooth() +
  scale_colour_manual(values = bb_palette$col) +
  scale_fill_manual(values = bb_palette$col) +
  theme_bw() +
  facet_wrap(vars(maingroup), scales = "free_y", ncol = 3) +
  labs(y = "No. species per day",
       x = "Date",
       colour = "Group",
       fill = "Group",
       linetype = "Record period")

ggsave(paste0("outputs/Summary/Nsp_2022.jpg"), 
       width = 2000,
       height = 1500, 
       units = "px",
       dpi = 300)

dat_september <- dat_all %>%
  filter(!record_period == "historical") %>%
  mutate(eventDate = ymd(std_date),
         month = month(eventDate),
         day = day(eventDate)) %>%
  filter(month == 9)

dat_september %>%
  mutate(date = dmy(paste(day, month, "2022", sep = "-"))) %>%
  group_by(date, maingroup, record_period, year) %>%
  count() %>%
  ggplot(aes(x = factor(record_period, labels = c("Bioblitz", "Baseline")), y = n, fill = maingroup)) +
  geom_boxplot() +
  scale_fill_manual(values = bb_palette$col) +
  theme_bw() +
  facet_wrap(vars(maingroup), scales = "free_y", ncol = 3) +
  labs(y = "No. observations per day",
       x = "Record period",
       fill = "Group")

ggsave(paste0("outputs/Summary/Nrecords_aov.jpg"), 
       width = 2000,
       height = 1500, 
       units = "px",
       dpi = 300)

sep_av <- aov(n ~ record_period * maingroup,
              dat = dat_september %>%
                mutate(date = dmy(paste(day, month, "2022", sep = "-"))) %>%
                group_by(date, maingroup, record_period, year) %>%
                count())

summary(sep_av)

sep_res <- lsmeans(sep_av, pairwise ~ record_period * maingroup, adjust = "bonferroni")

sep_cont <- summary(sep_res[[2]])[!is.na(summary(sep_res[[2]])[,4]),] %>%
  filter(str_count(contrast, "Arachnids") == 2 |
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
           str_count(contrast, "Reptiles") == 2)

sep_cont %>% write_csv("outputs/Summary/Nrecords_aov.csv")

dat_september %>%
  mutate(date = dmy(paste(day, month, "2022", sep = "-"))) %>%
  group_by(date, maingroup, record_period, year, morphospecies) %>%
  count() %>%
  select(-n) %>%
  group_by(date, maingroup, record_period, year) %>%
  count() %>%
  ggplot(aes(x = factor(record_period, labels = c("Bioblitz", "Baseline")), y = n, fill = maingroup)) +
  geom_boxplot() +
  scale_fill_manual(values = bb_palette$col) +
  theme_bw() +
  facet_wrap(vars(maingroup), scales = "free_y", ncol = 3) +
  labs(y = "No. sp per day",
       x = "Record period",
       fill = "Group")

ggsave(paste0("outputs/Summary/Nsp_aov.jpg"), 
       width = 2000,
       height = 1500, 
       units = "px",
       dpi = 300)

sep_av <- aov(n ~ record_period * maingroup,
              dat = dat_september %>%
                mutate(date = dmy(paste(day, month, "2022", sep = "-"))) %>%
                group_by(date, maingroup, record_period, year, morphospecies) %>%
                count() %>%
                select(-n) %>%
                group_by(date, maingroup, record_period, year) %>%
                count())

summary(sep_av)

sep_res <- lsmeans(sep_av, pairwise ~ record_period * maingroup, adjust = "bonferroni")

sep_cont <- summary(sep_res[[2]])[!is.na(summary(sep_res[[2]])[,4]),] %>%
  filter(str_count(contrast, "Arachnids") == 2 |
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
           str_count(contrast, "Reptiles") == 2)

sep_cont %>% write_csv("outputs/Summary/Nsp_aov.csv")

## emerging story plot
all_emerging <- read_csv("outputs/Table 2/emerging_ALL.csv") %>%
  rowwise() %>%
  mutate(New = sum(Nnewnative, Nnewintro)) %>%
  select(-c(Nnew, Pnewnative, Pnewintro, Nnewnative, Nnewintro)) %>%
  rename(Rediscovered = Nlazarus) %>%
  full_join(read_csv("outputs/Table 1/Splist_Nrecords_ALL.csv") %>% mutate(extinct = case_when(is.na(first_record_after1991) ~ "Extinct")) %>% drop_na(extinct) %>% count(maingroup) %>% rename(Extinct = n)) %>%
  full_join(read_csv("outputs/Table 1/Splist_Nrecords_ALL.csv") %>% mutate(extant = case_when(is.na(new) & is.na(lazarus) & is.na(new_introduced) & !is.na(first_record_after1991) ~ "Extant")) %>% drop_na(extant) %>% count(maingroup) %>% rename(Extant = n)) %>%
  full_join(read_csv("outputs/Table 1/Splist_Nrecords_ALL.csv") %>% mutate(known = case_when(is.na(new) & is.na(new_introduced) ~ "Known biodversity")) %>% drop_na(known) %>% count(maingroup) %>% rename(`Known biodiversity` = n)) %>% full_join(read_csv("outputs/Table 1/Splist_Nrecords_ALL.csv") %>% mutate(extinct = case_when(is.na(first_record_after1991) ~ "Extinct"), updated = case_when(is.na(extinct) ~ "Updated biodiversity")) %>% drop_na(updated) %>% count(maingroup) %>% rename(`Updated biodiversity` = n)) %>%
  full_join(read_csv("outputs/Table 1/Splist_Nrecords_ALL.csv") %>% mutate(refound = case_when(is.na(new) & is.na(lazarus) & is.na(new_introduced) & !is.na(first_record_after1991) & most_recent_record == 2022 ~ "Re-found")) %>% drop_na(refound) %>% count(maingroup) %>% rename(`Re-found` = n)) %>%
  full_join(read_csv("outputs/Table 1/Splist_Nrecords_ALL.csv") %>% mutate(nrefound = case_when(is.na(new) & is.na(lazarus) & is.na(new_introduced) & !is.na(first_record_after1991) & most_recent_record < 2022 ~ "Not re-found")) %>% drop_na(nrefound) %>% count(maingroup) %>% rename(`Not re-found` = n)) %>%
  full_join(read_csv("outputs/Table 1/Splist_Nrecords_ALL.csv") %>% mutate(pextinct = case_when(is.na(first_record_after1991) | lazarus == "Lazarus" ~ "Possibly extinct")) %>% drop_na(pextinct) %>% count(maingroup) %>% rename(`Possibly extinct` = n)) %>%
  gather(story, n, -maingroup)

origin_emerging <- read_csv("outputs/Table 1/Splist_Nrecords_ALL.csv") %>% drop_na(lazarus) %>% count(Origin) %>% rename(Rediscovered = n) %>%
  full_join(read_csv("outputs/Table 1/Splist_Nrecords_ALL.csv") %>% drop_na(new) %>% count(Origin) %>% rename(New = n)) %>%
  full_join(read_csv("outputs/Table 1/Splist_Nrecords_ALL.csv") %>% mutate(extinct = case_when(is.na(first_record_after1991) ~ "Extinct")) %>% drop_na(extinct) %>% count(Origin) %>% rename(Extinct = n)) %>%
  full_join(read_csv("outputs/Table 1/Splist_Nrecords_ALL.csv") %>% mutate(extant = case_when(is.na(new) & is.na(lazarus) & is.na(new_introduced) & !is.na(first_record_after1991) ~ "Extant")) %>% drop_na(extant) %>% count(Origin) %>% rename(Extant = n)) %>%
  full_join(read_csv("outputs/Table 1/Splist_Nrecords_ALL.csv") %>% mutate(known = case_when(is.na(new) & is.na(new_introduced) ~ "Known biodversity")) %>% drop_na(known) %>% count(Origin) %>% rename(`Known biodiversity` = n)) %>% full_join(read_csv("outputs/Table 1/Splist_Nrecords_ALL.csv") %>% mutate(extinct = case_when(is.na(first_record_after1991) ~ "Extinct"), updated = case_when(is.na(extinct) ~ "Updated biodiversity")) %>% drop_na(updated) %>% count(Origin) %>% rename(`Updated biodiversity` = n)) %>%
  full_join(read_csv("outputs/Table 1/Splist_Nrecords_ALL.csv") %>% mutate(refound = case_when(is.na(new) & is.na(lazarus) & is.na(new_introduced) & !is.na(first_record_after1991) & most_recent_record == 2022 ~ "Re-found")) %>% drop_na(refound) %>% count(Origin) %>% rename(`Re-found` = n)) %>%
  full_join(read_csv("outputs/Table 1/Splist_Nrecords_ALL.csv") %>% mutate(nrefound = case_when(is.na(new) & is.na(lazarus) & is.na(new_introduced) & !is.na(first_record_after1991) & most_recent_record < 2022 ~ "Not re-found")) %>% drop_na(nrefound) %>% count(Origin) %>% rename(`Not re-found` = n)) %>%
  full_join(read_csv("outputs/Table 1/Splist_Nrecords_ALL.csv") %>% mutate(pextinct = case_when(is.na(first_record_after1991) | lazarus == "Lazarus" ~ "Possibly extinct")) %>% drop_na(pextinct) %>% count(Origin) %>% rename(`Possibly extinct` = n)) %>%
  gather(story, n, -Origin)

stories <- unique(all_emerging$story)

for (i in stories){
  story <- all_emerging %>% filter(story == i, n > 0)
  story_origin <- origin_emerging %>% filter(story == i, n > 0) %>% drop_na(Origin)
  
  p_story <- story %>%
    ggplot() +
    geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0.6, r = 1,
                     amount = n,
                     fill = maingroup), stat = "pie", colour = NA) +
    geom_text(x = 0, y = 0, aes(label = sum(n)), size = 20) +
    scale_fill_manual(values = filter(bb_palette, maingroup %in% story$maingroup)$col, name = "Group") +
    coord_fixed() +
    theme_void() +
    theme(legend.position = "none")
  
  if(i %in% c("Known biodiversity", "Updated biodiversity")) {
    library(grImport2)
    library(rphylopic)
    
    phylopics <- list(Arachnids = png::readPNG("phylopic/Arachnids.png"),
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
                      Reptiles = png::readPNG("phylopic/Reptiles.png"))
    
    coords_ALL <- data.frame(X = 0, Y = 0)
    
    radius_ALL <- 1.05
    
    angles_ALL <- calc_angles(story$n)
    
    angles_ALL <- (angles_ALL + replace_na(lag(angles_ALL), 0))/2
    
    palette <- bb_palette %>%
      mutate(angle = angles_ALL,
             x = coords_ALL$X + radius_ALL * sin(pi * 2 * angle/360),
             y = coords_ALL$Y + radius_ALL * cos(pi * 2 * angle/360))
    
    for (k in story$maingroup){
      p_story <- p_story +
        add_phylopic(phylopics[[k]],
                     x = as.numeric(filter(palette, maingroup == k)[,4]),
                     y = as.numeric(filter(palette, maingroup == k)[,5]),
                     ysize = 0.1,
                     color = as.character(filter(palette, maingroup == k)[,2]),
                     alpha = 1)
    }
  }
  
  p_story_origin <- story_origin %>%
    ggplot() +
    geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0.6, r = 1,
                     amount = n,
                     fill = Origin), stat = "pie", colour = NA) +
    geom_text(x = 0, y = 0, aes(label = sum(n)), size = 20) +
    scale_fill_manual(values = c("deepskyblue4", "darkorchid4"), name = "Origin") +
    coord_fixed() +
    theme_void() +
    theme(legend.position = "none")
  
  ggsave(paste0("outputs/Summary/emerging_", i,".pdf"),
         plot = p_story,
         width = 2000,
         height = 1500, 
         units = "px",
         dpi = 300)
  
  ggsave(paste0("outputs/Summary/origin_", i,".pdf"),
         plot = p_story_origin,
         width = 2000,
         height = 1500, 
         units = "px",
         dpi = 300)
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
  mutate(lga = case_when(lga == "bayside" ~ "BSD",
                         lga == "cardinia" ~ "CAR",
                         lga == "casey" ~ "CAS",
                         lga == "frankston" ~ "FRA",
                         lga == "gd" ~ "GDA",
                         lga == "kingston" ~ "KIN",
                         lga == "knox" ~ "KNO",
                         lga == "monash" ~ "MON",
                         lga == "mps" ~ "MPE",
                         lga == "yr" ~  "YRA"))

dat_events_unique <- dat_events %>%
  select(event_year, event_description, start_date, end_date) %>%
  unique()

dat_all$event_description <- "Baseline"

for(i in 1:nrow(dat_events_unique)){
  event_interval <- interval(dmy(dat_events_unique[i,3]), dmy(dat_events_unique[i,4]))
  event_year <- as.character(dat_events_unique[i,1])
  event_d <- as.character(dat_events_unique[i,2])
  
  dat_all <- dat_all %>%
    mutate(record_period = case_when(ymd(std_date) %within% event_interval ~ event_year,
                                     T ~ record_period),
           event_description = case_when(ymd(std_date) %within% event_interval ~ event_d,
                                         T ~ event_description)) 
}

dat_events_lga <- dat_events %>%
  select(lga, event_year) %>%
  drop_na(lga) %>%
  mutate(event_year = factor(event_year, levels = c("gsb_20", "cnc_21", "bb_21", "gsb_21", "cnc_22", "bb_22", "gsb_22"))) %>%
  group_by(lga) %>%
  arrange(event_year) %>%
  mutate(event_number = row_number()) %>%
  ungroup() %>%
  mutate(event_number = paste0("event_", letters[event_number])) %>%
  rename(council = lga,
         record_period = event_year)

dat_events_lga <- dat_events_lga %>%
  full_join(dat_events %>%
              select(lga, event_year) %>%
              rename(council = lga,
                     record_period = event_year) %>%
              expand(council, record_period)) %>%
  drop_na(council) %>%
  mutate(record_period = factor(record_period, levels = c("gsb_20", "inter_gsb20_cnc21", "cnc_21", "inter_cnc21_bb21", "bb_21", "inter_bb21_gsb21", "gsb_21", "inter_gsb21_cnc22", "cnc_22", "inter_cnc22_bb22", "bb_22", "inter_bb22_gsb22", "gsb_22"))) %>%
  arrange(council, record_period) %>%
  left_join(dat_events_lga) %>%
  group_by(council, event_number) %>%
  ungroup()

council_events <- list()
for (i in unique(dat_events_lga$council)){
  council_events[[i]] <- dat_events_lga %>%
    filter(council == i)
  
  index = 1
  for (k in 1:nrow(council_events[[i]])){
    if (is.na(council_events[[i]][k, 3])) {
      council_events[[i]][k, 3] <- paste0("interperiod_", letters[index])
    } else {
      index = index + 1
    }
  }
}

dat_events_lga <- bind_rows(council_events)

dat_all <- dat_all %>%
  left_join(dat_events_lga)

dat_all %>% write_csv("data/Records w Events v1 31Mar23.csv")
