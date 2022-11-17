### NEW BIOBLITZ CODE ###

library(tidyverse)
library(sf)
library(terra)
library(eks)
library(scatterpie)

## LOAD DATA
# Load the clean dataset of biodiversity records, and the cypher file 
dALL.cleanup.groups_corrected <- read_csv("data/ALLrecords_cleanup_groups_corrected_for_CESAR_15Nov2022.csv")

councils_full <- unique(dALL.cleanup.groups_corrected$council)

dALL.cleanup.groups_corrected <- dALL.cleanup.groups_corrected %>%
  mutate(council = case_when(council == "Cardinia" ~ "CAR",
                             council == "Casey" ~ "CAS",
                             council == "Frankston" ~ "FRA",
                             council == "Greater Dandenong" ~ "GDA",
                             council == "Kingston" ~ "KIN",
                             council == "Knox" ~ "KNO",
                             council == "Monash" ~ "MON",
                             council == "Mornington Peninsula" ~ "MPE",
                             council == "Yarra Ranges" ~  "YRA"))

# fix some issues with duplicate groupings
double_gen <- dALL.cleanup.groups_corrected %>%
  filter(morphospecies %in% (dALL.cleanup.groups_corrected %>%
                               select(morphospecies, genus) %>%
                               unique() %>%
                               group_by(morphospecies) %>%
                               count() %>%
                               filter(n > 1))$morphospecies)

dALL.cleanup.groups_corrected[which(dALL.cleanup.groups_corrected$morphospecies %in% double_gen$morphospecies), "genus"] <- sapply(str_split(double_gen$morphospecies, " "), "[[", 1)

dALL.cleanup.groups_corrected <- dALL.cleanup.groups_corrected %>%
  mutate(family = case_when(morphospecies == "Aira elegantissima" ~ "Poaceae",
                            morphospecies == "Ascocoryne sarcoides" ~ "Gelatinodiscaceae",
                            morphospecies == "Hemimycena lactea" ~ "Mycenaceae",
                            morphospecies == "Lejeunea drummondii" ~ "Lejeuneaceae",
                            morphospecies == "Megalurus gramineus" ~ "Megaluridae",
                            morphospecies == "Pieris rapae" ~ "Pieridae",
                            morphospecies == "Trachycarpus fortunei" ~ "Arecaceae",
                            morphospecies == "Xanthorrhoea australis" ~ "Xanthorrhoeaceae",
                            morphospecies == "Xanthorrhoea minor" ~ "Xanthorrhoeaceae",
                            T ~ family),
         maingroup = case_when(morphospecies == "Aira elegantissima" ~ "Plants",
                               morphospecies == "Pieris rapae" ~ "Insects",
                               morphospecies == "Serpula himantioides" ~ "Fungi",
                               T ~ maingroup),
         subgroup = case_when(morphospecies %in% c("Adiantum aethiopicum", "Anogramma leptophylla", "Cheilanthes austrotenuifolia", "Pellaea falcata", "Pteris tremula") ~ "Ferns",
                              morphospecies == "Aira elegantissima" ~ "Monocots",
                              morphospecies %in% c("Aploactisoma milesii", "Gymnapistes marmoratus", "Lepidotrigla papilio") ~ "Lionfishes and Sculpins",
                              family %in% c("Ardeiae", "Phalacrocoracidae", "Threskiornithidae") ~ "Waterbirds",
                              morphospecies %in% c("Borago officinalis", "Cynoglossum australe", "Cynoglossum australe", "Echium plantagineum", "Hackelia latifolia", "Hibbertia acicularis", "Hibbertia australis", "Hibbertia empetrifolia", "Hibbertia fasciculata", "Hibbertia riparia", "Hibbertia sericea", "Myosotis scorpioides", "Myosotis sylvatica", "Pentaglottis sempervirens") ~ "Eudicots",
                              family %in% c("Gobiidae", "Callionymidae", "Percichthyidae") ~ "Perciformes",
                              morphospecies == "Ikeda" ~ "Spoon worms",
                              morphospecies == "Lejeunea drummondii" ~ "Liverworts",
                              morphospecies == "Morus serrator" ~ "Seabirds",
                              morphospecies == "Pieris rapae" ~ "Butterflies",
                              morphospecies == "Serpula himantioides" ~ "Club fungi",
                              morphospecies == "Trachycarpus fortunei" ~ "Monocots"))
# For future trackability: Esti saved the clean dataset of biodiversity records from 'Biodiversity blitz 2021 R ENVIR v15 13Jul22.RData': write.csv(dALL.cleanup.groups_corrected, file="ALLrecords_cleanup_groups_corrected_for_CESAR_15Nov2022.csv")
vba <- read_csv("data/AH Cypher v14 15Feb22.csv") 

# generate colour palette
bb_palette <- tibble(maingroup = sort(unique(dALL.cleanup.groups_corrected$maingroup)),
                     col = c("black",
                             "darkorchid",
                             "slategrey",
                             "olivedrab3",
                             "darkgrey",
                             "darkred",
                             "bisque1",
                             "darkorange2",
                             "darksalmon",
                             "darkolivegreen",
                             "darkblue",
                             "deeppink3"))

# load map of Victoria and reproject to Lambert equal-area
map_VIC <- read_sf("SpatialData/AD_LGA_AREA_POLYGON.shp") %>%
  filter(STATE == "VIC")

# filter based on LGAs
map_LGA <- map_VIC %>%
  st_transform(crs = 3112) %>%
  filter(NAME %in% str_to_upper(councils_full)) %>%
  rename(council = NAME) %>%
  mutate(council = case_when(council == "CARDINIA" ~ "CAR",
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
#                          cellsize = 500)  %>%
#   st_as_sf() %>%
#   filter(st_intersects(st_geometry(.), summarise(map_LGA), sparse = F)[,1])
# 
# write_sf(grid_LGA, "SpatialData/grid_LGA.shp")

grid_LGA <- read_sf("SpatialData/grid_LGA.shp")

# convert data to grid form
occ_GRID <- dALL.cleanup.groups_corrected %>%
  filter(record_period %in% c("recent", "bb21")) %>%
  select(locality, decimalLatitude, decimalLongitude, identifiedBy, recordID, record_period, keep, maingroup) %>%
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = st_crs(map_VIC)) %>%
  st_transform(crs = 3112) %>%
  st_intersection(grid_LGA) %>%
  as_tibble() %>%
  select(-geometry) %>%
  group_by(FID, record_period) %>%
  summarise(occ = n())

## GENERATE DELIVERABLES: TOTAL AND PER COUNCIL
councils <- c("ALL", unique(dALL.cleanup.groups_corrected$council))

for (i in councils){
  if(i == "ALL") {
    dat <- dALL.cleanup.groups_corrected
    grid_dat <- grid_LGA
  } else {
    dat <- dALL.cleanup.groups_corrected %>% filter(council == i)
    grid_dat <- grid_LGA %>%
      filter(st_intersects(st_geometry(.), map_LGA %>% filter(council == i), sparse = F)[,1])
  }
  
  # TABLE 1: species list, Nrecords total/period, FFG/EPBC, origin, common name, etc...
  Splist.Nrecords <- dat %>%
    group_by(morphospecies, maingroup, genus, family, subgroup, record_period) %>%
    summarise(Nrecords_total = sum(occ)) %>%
    spread(record_period, Nrecords_total) %>%
    mutate_if(is.numeric, replace_na, 0) %>%
    mutate(total = sum(bb21, historical, recent, na.rm = T),
           Perc_records_total = total/nrow(dALL.cleanup.groups_corrected)*100,
           Perc_records_historical = historical/total*100,
           Perc_records_recent = recent/total*100,
           Perc_records_bioblitz = bb21/total*100,
           distribution = case_when(Perc_records_bioblitz == 0 ~ 1,
                                    Perc_records_bioblitz == 100 ~ 2,
                                    T ~ 3)) %>%
    ungroup() %>%
    left_join(dALL.cleanup.groups_corrected %>%
                group_by(morphospecies) %>%
                summarise(most_recent_record = max(year))) %>%
    left_join(dALL.cleanup.groups_corrected %>%
                filter(record_period %in% c("recent", "bb21")) %>%
                group_by(morphospecies) %>%
                summarise(first_record_after1991 = min(year))) %>%
    left_join(vba %>% select(-type), by = c("morphospecies" = "sname")) %>%
    mutate(origin = case_when(origin %in% c("Native", "native", "Native but some stands may be alien") ~ "Native",
                              origin %in% c("Introduced", "Naturalised alien", "Introduced but doubt it ever established a population in victoria", "introduced") ~  "Introduced"),
           lazarus = case_when(historical > 0 & recent == 0 & bb21 > 0 ~ "Lazarus"),
           new = case_when(historical == 0 & recent == 0 & bb21 > 0 ~ "New"),
           new_introduced = case_when(historical == 0 & recent == 0 & bb21 > 0 & origin == "Introduced" ~ "New introduced")) %>%
    rename(CommonName = cname,
           Origin = origin,
           EPBC = epbc,
           FFGAct = ffg,
           Nrecords_total = total,
           Nrecords_historical = historical,
           Nrecords_recent = recent,
           Nrecords_bioblitz = bb21) %>%
    select(maingroup, morphospecies, CommonName, Origin, family, genus, subgroup, EPBC, FFGAct, Nrecords_total, Perc_records_total, Nrecords_historical, Perc_records_historical, Nrecords_recent, Perc_records_recent, Nrecords_bioblitz, Perc_records_bioblitz, most_recent_record, first_record_after1991, distribution, lazarus, new, new_introduced)
  
  write_csv(Splist.Nrecords, paste0("outputs/Table 1/Splist_Nrecords_", i, ".csv"), na = "")
  
  ## TABLE 2: emerging stories [lazarus, new, new_introduced]
  emerging <- Splist.Nrecords %>%
    mutate(new = case_when(new == "New" ~ 1,
                           T ~ 0),
           lazarus = case_when(lazarus == "Lazarus" ~ 1,
                               T ~ 0),
           new_introduced = case_when(new_introduced == "New introduced" ~ 1,
                                      T ~ 0)) %>%
    group_by(maingroup) %>%
    summarise(Nnew = sum(new),
              Nnewintro = sum(new_introduced),
              Nlazarus = sum(lazarus)) %>%
    mutate(Nnewnative = Nnew - Nnewintro,
           Pnewintro = (Nnewintro/Nnew)*100,
           Pnewnative = (Nnewnative/Nnew)*100) %>%
    select(maingroup, Nnew, Nnewnative, Pnewnative, Nnewintro, Pnewintro, Nlazarus)
  
  write_csv(emerging, paste0("outputs/Table 2/emerging_", i, ".csv"), na = "")
  
  ## TABLE 3: species richness by main group, by period, unique species
  SpRichness.Groups <- Splist.Nrecords %>%
    group_by(maingroup) %>%
    summarise(Nrecords_total_from1991 = sum(Nrecords_recent + Nrecords_bioblitz),
              Nrecords_recent = sum(Nrecords_recent),
              Nrecords_bioblitz = sum(Nrecords_bioblitz)) %>%
    left_join(Splist.Nrecords %>%
                filter(Nrecords_recent > 0 | Nrecords_bioblitz > 0) %>%
                group_by(maingroup) %>%
                summarise(Nsp_total_from1991 = length(morphospecies))) %>%
    left_join(Splist.Nrecords %>%
                filter(Nrecords_recent > 0) %>%
                group_by(maingroup) %>%
                summarise(Nsp_recent = length(morphospecies))) %>%
    left_join(Splist.Nrecords %>%
                filter(Nrecords_bioblitz > 0) %>%
                group_by(maingroup) %>%
                summarise(Nsp_bioblitz = length(morphospecies))) %>%
    left_join(Splist.Nrecords %>%
                filter(Nrecords_recent > 0 & Nrecords_bioblitz == 0) %>%
                group_by(maingroup) %>%
                summarise(Nuniquesp_recent = length(morphospecies))) %>%
    left_join(Splist.Nrecords %>%
                filter(Nrecords_bioblitz > 0 & Nrecords_recent == 0) %>%
                group_by(maingroup) %>%
                summarise(Nuniquesp_bioblitz = length(morphospecies)))
  
  write_csv(SpRichness.Groups, paste0("outputs/Table 3/Species_richness_by_group_", i, ".csv"), na = "")
  
  ## TABLE 4: species richness by subgroup, by period, unique species
  SpRichness.Subgroups <- Splist.Nrecords %>%
    group_by(maingroup, subgroup) %>%
    summarise(Nrecords_total_from1991 = sum(Nrecords_recent + Nrecords_bioblitz),
              Nrecords_recent = sum(Nrecords_recent),
              Nrecords_bioblitz = sum(Nrecords_bioblitz)) %>%
    left_join(Splist.Nrecords %>%
                filter(Nrecords_recent > 0 | Nrecords_bioblitz > 0) %>%
                group_by(maingroup, subgroup) %>%
                summarise(Nsp_total_from1991 = length(morphospecies))) %>%
    left_join(Splist.Nrecords %>%
                filter(Nrecords_recent > 0) %>%
                group_by(maingroup, subgroup) %>%
                summarise(Nsp_recent = length(morphospecies))) %>%
    left_join(Splist.Nrecords %>%
                filter(Nrecords_bioblitz > 0) %>%
                group_by(maingroup, subgroup) %>%
                summarise(Nsp_bioblitz = length(morphospecies))) %>%
    left_join(Splist.Nrecords %>%
                filter(Nrecords_recent > 0 & Nrecords_bioblitz == 0) %>%
                group_by(maingroup, subgroup) %>%
                summarise(Nuniquesp_recent = length(morphospecies))) %>%
    left_join(Splist.Nrecords %>%
                filter(Nrecords_bioblitz > 0 & Nrecords_recent == 0) %>%
                group_by(maingroup, subgroup) %>%
                summarise(Nuniquesp_bioblitz = length(morphospecies))) %>%
    drop_na(subgroup) %>%
    arrange(subgroup)
  
  write_csv(SpRichness.Groups, paste0("outputs/Table 4/Species_richness_by_subgroup_", i, ".csv"), na = "")
  
  ## TABLE 7: most common species from BB21
  Common.sp <- Splist.Nrecords %>%
    mutate(Nrecords_total_from1991 = Nrecords_recent + Nrecords_bioblitz,
           Perc_records_bioblitz = (Nrecords_bioblitz/Nrecords_total_from1991)*100) %>%
    arrange(desc(Nrecords_bioblitz)) %>%
    slice(1:25) %>%
    select(morphospecies, CommonName, genus, family, subgroup, maingroup, Origin, Nrecords_bioblitz, Nrecords_total_from1991, Perc_records_bioblitz)
  
  write_csv(Common.sp, paste0("outputs/Table 7/Common_species_", i, ".csv"), na = "")
  
  ## FIGURE 1
  Common.sp %>%
    mutate(CommonName = case_when(is.na(CommonName) ~ morphospecies,
                                  T ~ CommonName)) %>%
    ggplot(aes(x = reorder(CommonName, -Nrecords_bioblitz), y = Nrecords_bioblitz, fill = maingroup)) +
    geom_bar(stat = "identity") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values = as.character(bb_palette %>%
                                              filter(maingroup %in% Common.sp$maingroup) %>%
                                              select(col) %>%
                                              as_vector()),
                      name = "Group") +
    labs(x = "Species",
         y = "No. observations",
         title = "Bioblitz 2021")
  
  ggsave(paste0("outputs/Figure 1/Common_species_", i, ".jpg"), 
         width = 2000,
         height = 1500, 
         units = "px",
         dpi = 300)
  
  ## TABLE 8: most common group from BB21
  Common.group <- SpRichness.Groups %>%
    mutate(Perc_records_bioblitz = (Nrecords_bioblitz/Nrecords_total_from1991)*100,
           Nuniquesp_bioblitz = replace_na(Nuniquesp_bioblitz, 0),
           Perc_newsp_bioblitz = (Nuniquesp_bioblitz/Nsp_total_from1991)*100,
           Perc_records_bioblitz = replace_na(Perc_records_bioblitz, 0),
           Perc_newsp_bioblitz = replace_na(Perc_newsp_bioblitz, 0)) %>%
    arrange(desc(Nrecords_bioblitz)) %>%
    select(maingroup, Nrecords_total_from1991, Nrecords_bioblitz, Perc_records_bioblitz, Nuniquesp_bioblitz, Perc_newsp_bioblitz, Nsp_bioblitz)
  
  write_csv(Common.group, paste0("outputs/Table 8/Common_group_", i, ".csv"), na = "")
  
  ## FIGURE X: cumulative known species
  Splist.Nrecords %>%
    group_by(first_record_after1991, maingroup) %>%
    count() %>%
    mutate(first_record_after1991 = case_when(is.na(first_record_after1991) ~ 1990,
                                              T ~ first_record_after1991)) %>%
    arrange(first_record_after1991) %>%
    ungroup() %>%
    group_by(maingroup) %>%
    mutate(n = cumsum(n)) %>%
    ungroup() %>%
    spread(maingroup, n) %>%
    mutate_all(funs(case_when(first_record_after1991 == 2021 ~ max(., na.rm = T), 
                              TRUE ~ .))) %>%
    gather(maingroup, n, -first_record_after1991) %>%
    drop_na() %>%
    ggplot(aes(x = first_record_after1991, y = n, color = maingroup)) +
    geom_step() +
    theme_bw() +
    scale_color_manual(values = as.character(bb_palette %>%
                                               filter(maingroup %in% Splist.Nrecords$maingroup) %>%
                                               select(col) %>%
                                               as_vector()),
                       name = "Group") +
    facet_wrap(facets = vars(maingroup), ncol = 3, scales = "free_y") +
    labs(x = "Year", y = "No. species")
  
  ggsave(paste0("outputs/Figure X/Cumulative_species_", i, ".jpg"), 
         width = 2000,
         height = 1500, 
         units = "px",
         dpi = 300)
  
  ## FIGURE X2: no. observations
  dat %>%
    filter(!record_period == "historical") %>%
    group_by(maingroup, year) %>%
    count() %>%
    spread(maingroup, n) %>%
    gather(maingroup, n, -year) %>%
    mutate(n = replace_na(n, 0)) %>%
    ggplot(aes(x = year, y = n, color = maingroup)) +
    geom_line() +
    theme_bw() +
    scale_color_manual(values = as.character(bb_palette %>%
                                               filter(maingroup %in% filter(dat, !record_period == "historical")$maingroup) %>%
                                               select(col) %>%
                                               as_vector()),
                       name = "Group") +
    facet_wrap(facets = vars(maingroup), ncol = 3, scales = "free_y") +
    labs(x = "Year", y = "No. observations")
  
  ggsave(paste0("outputs/Figure X2/Observations_", i, ".jpg"), 
         width = 2000,
         height = 1500, 
         units = "px",
         dpi = 300)
  
  ## MAPS
  grid_dat %>%
    left_join(occ_GRID %>% filter(record_period == "recent")) %>%
    ggplot() +
    geom_sf(data = map_LGA, fill = "white", colour = "dark grey", inherit.aes = F) +
    geom_sf(aes(fill = occ), colour = NA, alpha = .8) +
    scale_fill_viridis_c(na.value = NA, option = "magma", name = "Nrecords") +
    theme_bw() +
    coord_sf(xlim = st_bbox(grid_dat)[c(1,3)],
             ylim = st_bbox(grid_dat)[c(2,4)]) +
    ggtitle("Recent")
  
  ggsave(paste0("outputs/Maps/Nrecords_", i, "_recent.jpg"), 
         width = 2000,
         height = 1500, 
         units = "px",
         dpi = 300)
  
  grid_dat %>%
    left_join(occ_GRID %>% filter(record_period == "bb21")) %>%
    ggplot() +
    geom_sf(data = map_LGA, fill = "white", colour = "dark grey", inherit.aes = F) +
    geom_sf(aes(fill = occ), colour = NA, alpha = .8) +
    scale_fill_viridis_c(na.value = NA, option = "magma", name = "Nrecords") +
    theme_bw() +
    coord_sf(xlim = st_bbox(grid_dat)[c(1,3)],
             ylim = st_bbox(grid_dat)[c(2,4)]) +
    ggtitle("Bioblitz 2021")
  
  ggsave(paste0("outputs/Maps/Nrecords_", i, "_bioblitz.jpg"), 
         width = 2000,
         height = 1500, 
         units = "px",
         dpi = 300)
  
  kde_recent <- st_kde(dat %>%
                         filter(record_period == "recent") %>%
                         select(locality, decimalLatitude, decimalLongitude, identifiedBy, recordID, record_period, keep, maingroup) %>%
                         st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = st_crs(map_VIC)) %>%
                         st_transform(crs = 3112))
  
  kde_recent %>%
    ggplot() +
    geom_sf(data = map_LGA, fill = "white", colour = "dark grey", inherit.aes = F) +
    geom_sf(data = st_get_contour(kde_recent, cont = c(5, 25, 50, 75, 95)),
            aes(fill = label_percent(contlabel)), colour = NA, alpha = .8) +
    scale_fill_viridis_d(na.value = NA, option = "magma", name = "Nrecords") +
    theme_bw() +
    coord_sf(xlim = st_bbox(grid_dat)[c(1,3)],
             ylim = st_bbox(grid_dat)[c(2,4)]) +
    ggtitle("Recent")
  
  ggsave(paste0("outputs/Maps/Ndensity_", i, "_recent.jpg"), 
         width = 2000,
         height = 1500, 
         units = "px",
         dpi = 300)
  
  kde_bioblitz <- st_kde(dat %>%
                           filter(record_period == "bb21") %>%
                           select(locality, decimalLatitude, decimalLongitude, identifiedBy, recordID, record_period, keep, maingroup) %>%
                           st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = st_crs(map_VIC)) %>%
                           st_transform(crs = 3112))
  
  kde_bioblitz %>%
    ggplot() +
    geom_sf(data = map_LGA, fill = "white", colour = "dark grey", inherit.aes = F) +
    geom_sf(data = st_get_contour(kde_bioblitz, cont = c(5, 25, 50, 75, 95)),
            aes(fill = label_percent(contlabel)), colour = NA, alpha = .8) +
    scale_fill_viridis_d(na.value = NA, option = "magma", name = "Nrecords") +
    theme_bw() +
    coord_sf(xlim = st_bbox(grid_dat)[c(1,3)],
             ylim = st_bbox(grid_dat)[c(2,4)]) +
    ggtitle("Bioblitz 2021")
  
  ggsave(paste0("outputs/Maps/Ndensity_", i, "_bioblitz.jpg"), 
         width = 2000,
         height = 1500, 
         units = "px",
         dpi = 300)
}

## GENERATE FANCY LGA + COMMON SPECIES PIECHART MAPS
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

map_LGA %>%
  ggplot() +
  geom_sf(fill = "white", colour = "dark grey", inherit.aes = F) +
  geom_scatterpie(aes(x = X, y = Y, group = council, r = r/1.5),
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
                         labeller = function(x) x*1.5) +
  ggtitle("Bioblitz 2021 No. observations")

ggsave(paste0("outputs/Maps/Common_groups_Nrecords.jpg"), 
       width = 2000,
       height = 1500, 
       units = "px",
       dpi = 300)

map_LGA %>%
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
  ggtitle("Bioblitz 2021 No. species")

ggsave(paste0("outputs/Maps/Common_groups_Nsp.jpg"), 
       width = 2000,
       height = 1500, 
       units = "px",
       dpi = 300)
