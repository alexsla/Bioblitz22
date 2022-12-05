#### SUMMARY FUNCTIONS ####

### Function to generate summary table of number of records per species
Splist.Nrecords <- function(dat, full_dat, cypher, path, i){
  dat <-
    dat %>%
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
    mutate(CommonName = str_replace_all(str_to_title(str_replace_all(CommonName, "-", "_")), "_", "-")) %>%
    select(maingroup, morphospecies, CommonName, Origin, family, genus, subgroup, EPBC, FFGAct, Nrecords_total, Perc_records_total, Nrecords_historical, Perc_records_historical, Nrecords_recent, Perc_records_recent, Nrecords_bioblitz, Perc_records_bioblitz, most_recent_record, first_record_after1991, distribution, lazarus, new, new_introduced)
  
  write_csv(dat, paste0(path, "/Splist_Nrecords_", i, ".csv"), na = "")
}

### Function to generate emerging stories (lazarus, new species) per main group
emerging <- function(file, path, i) {
  dat <-
    read_csv(file) %>%
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
  
  write_csv(dat, paste0(path, "/emerging_", i, ".csv"), na = "")
}

### Function to calculate species richness per grouping variable
SpRichness <- function(file, path, group, i) {
  if(deparse(substitute(group)) == "maingroup")
    dat <- read_csv(file) %>%
      group_by({{group}}) else
        dat <- read_csv(file) %>%
          group_by(maingroup, {{group}})  
      
      dat <-
        dat %>%
        summarise(Nrecords_total_from1991 = sum(Nrecords_recent + Nrecords_bioblitz),
                  Nrecords_recent = sum(Nrecords_recent),
                  Nrecords_bioblitz = sum(Nrecords_bioblitz)) %>%
        left_join(dat %>%
                    filter(Nrecords_recent > 0 | Nrecords_bioblitz > 0) %>%
                    group_by({{group}}) %>%
                    summarise(Nsp_total_from1991 = length(morphospecies))) %>%
        left_join(dat %>%
                    filter(Nrecords_recent > 0) %>%
                    group_by({{group}}) %>%
                    summarise(Nsp_recent = length(morphospecies))) %>%
        left_join(dat %>%
                    filter(Nrecords_bioblitz > 0) %>%
                    group_by({{group}}) %>%
                    summarise(Nsp_bioblitz = length(morphospecies))) %>%
        left_join(dat %>%
                    filter(Nrecords_recent > 0 & Nrecords_bioblitz == 0) %>%
                    group_by({{group}}) %>%
                    summarise(Nuniquesp_recent = length(morphospecies))) %>%
        left_join(dat %>%
                    filter(Nrecords_bioblitz > 0 & Nrecords_recent == 0) %>%
                    group_by({{group}}) %>%
                    summarise(Nuniquesp_bioblitz = length(morphospecies))) %>%
        drop_na({{group}})
      
      write_csv(dat, paste0(path, "/Species_richness_by_", deparse(substitute(group)), "_", i, ".csv"), na = "")
}

### Function to find most common species
Common.sp <- function(file, path, i){
  dat <-
    read_csv(file) %>%
    mutate(Nrecords_total_from1991 = Nrecords_recent + Nrecords_bioblitz,
           Perc_records_bioblitz = (Nrecords_bioblitz/Nrecords_total_from1991)*100) %>%
    arrange(desc(Nrecords_bioblitz)) %>%
    slice(1:25) %>%
    select(morphospecies, CommonName, genus, family, subgroup, maingroup, Origin, Nrecords_bioblitz, Nrecords_total_from1991, Perc_records_bioblitz)
  
  write_csv(dat, paste0(path, "/Common_species_", i, ".csv"), na = "")
}

### Function to find most common group per grouping variable
Common.group <- function(file, path, group, i){
  dat <- read_csv(file) %>%
    ungroup()
  
  dat <- dat %>%
    mutate(Perc_records_bioblitz = (Nrecords_bioblitz/Nrecords_total_from1991)*100,
           Nuniquesp_bioblitz = replace_na(Nuniquesp_bioblitz, 0),
           Perc_newsp_bioblitz = (Nuniquesp_bioblitz/Nsp_total_from1991)*100,
           Perc_records_bioblitz = replace_na(Perc_records_bioblitz, 0),
           Perc_newsp_bioblitz = replace_na(Perc_newsp_bioblitz, 0)) %>%
    arrange(desc(Nrecords_bioblitz))
  
  if(deparse(substitute(group)) == "maingroup")
    dat <- dat %>%
    select({{group}}, Nrecords_total_from1991, Nrecords_bioblitz, Perc_records_bioblitz, Nuniquesp_bioblitz, Perc_newsp_bioblitz, Nsp_bioblitz) else
      dat <- dat %>%
    select(maingroup, {{group}}, Nrecords_total_from1991, Nrecords_bioblitz, Perc_records_bioblitz, Nuniquesp_bioblitz, Perc_newsp_bioblitz, Nsp_bioblitz)
  
  write_csv(dat, paste0(path, "/Common_", deparse(substitute(group)), "_", i, ".csv"), na = "")
}

### Function to generate summary table by green space
Nrecords.GS <- function(dat, occ_dat, sp, path, i){
  dat <-
    dat %>%
    as_tibble() %>%
    select(FID, PARK_NAME) %>%
    left_join(occ_dat %>% filter(record_period == "bb21")) %>%
    drop_na(maingroup) %>%
    mutate(maingroup = case_when(maingroup == "Arachnids" ~ "Nrecords_ARA",
                                 maingroup == "Birds" ~ "Nrecords_BIR",
                                 maingroup == "Frogs" ~ "Nrecords_FRO",
                                 maingroup == "Fungi" ~ "Nrecords_FUN",
                                 maingroup == "Insects" ~ "Nrecords_INS",
                                 maingroup == "Other invertebrates" ~ "Nrecords_INV",
                                 maingroup == "Mammals" ~ "Nrecords_MAM",
                                 maingroup == "Plants" ~ "Nrecords_PLA",
                                 maingroup == "Reptiles" ~ "Nrecords_REP")) %>%
    spread(maingroup, occ) %>%
    rename(Greenspace_name = PARK_NAME) %>%
    select(-FID) %>%
    rowwise() %>%
    mutate(Nrecords_total = sum(across(where(is.numeric)))) %>%
    left_join(dat %>%
                as_tibble() %>%
                select(FID, PARK_NAME) %>%
                left_join(sp %>% filter(record_period == "bb21")) %>%
                drop_na(maingroup) %>%
                mutate(maingroup = case_when(maingroup == "Arachnids" ~ "Nsp_ARA",
                                             maingroup == "Birds" ~ "Nsp_BIR",
                                             maingroup == "Frogs" ~ "Nsp_FRO",
                                             maingroup == "Fungi" ~ "Nsp_FUN",
                                             maingroup == "Insects" ~ "Nsp_INS",
                                             maingroup == "Other invertebrates" ~ "Nsp_INV",
                                             maingroup == "Mammals" ~ "Nsp_MAM",
                                             maingroup == "Plants" ~ "Nsp_PLA",
                                             maingroup == "Reptiles" ~ "Nsp_REP")) %>%
                spread(maingroup, occ) %>%
                rename(Greenspace_name = PARK_NAME) %>%
                select(-FID) %>%
                rowwise() %>%
                mutate(Nsp_total = sum(across(where(is.numeric)))))
  
  write_csv(dat, paste0(path, "/", i, "bb_Nrecords_GS_groups.csv.csv"), na = "") 
}

### Function to generate summary table of participants
Nrecords.PAR <- function(dat, path, i){
  dat <-
    dat %>%
    filter(record_period == "bb21") %>%
    group_by(identifiedBy) %>%
    summarise(Nrecords_bb21 = sum(occ)) %>%
    rename(Participant_username = identifiedBy) %>%
    mutate(Perc_records_bb21 = (Nrecords_bb21/sum(Nrecords_bb21))*100) %>%
    arrange(desc(Nrecords_bb21))
  
  write_csv(dat, paste0(path, "/BBparticipants_", i, ".csv"), na = "")
}

### Function to generate summary table of binned participants
BBparticipants_bin <- function(dat, path, i){
  dat <-
    dat %>%
    filter(record_period == "bb21") %>%
    group_by(identifiedBy) %>%
    summarise(Nrecords_bb21 = sum(occ)) %>%
    count(Nrecords_bin = cut(Nrecords_bb21, breaks = c(0, 5, 50, 200, max(Nrecords_bb21)+1),
                             labels = c("1 to 5",
                                        "6 to 50",
                                        "51 to 200",
                                        paste0("200 to ", max(Nrecords_bb21)+1)))) %>%
    rename(Nparticipants = n)
  
  write_csv(dat, paste0(path, "/BBparticipants_bins_", i, ".csv"), na = "")
}

### Function to generate summary table of main groups
BBmaingroup <- function(dat, path, i){
  dat <-
    dat %>%
    filter(record_period == "bb21") %>%
    mutate(maingroup = factor(maingroup, levels =c("Arachnids","Birds","Frogs","Fungi","Insects","Mammals","Other invertebrates","Plants","Ray-finned fishes","Reptiles"))) %>%
    group_by(maingroup)
  
  dat <- left_join(dat %>%
                     select(maingroup, identifiedBy) %>%
                     unique() %>%
                     count() %>%
                     rename(Nparticipants_bb21 = n),
                   dat %>%
                     count() %>%
                     rename(Nrecords_bb21 = n))
  
  write_csv(dat, paste0(path, "/BBmaingroup_", i, ".csv"), na = "")
}



#### PLOTTING FUNCTIONS ####

### Function to plot common species barplot
Common.sp.p <- function(file, path, i, palette){
  dat <- read_csv(file)
  
  ggsave(paste0(path, "/Common_species_", i, ".jpg"),
         plot = dat %>%
           mutate(CommonName = case_when(is.na(CommonName) ~ morphospecies,
                                         T ~ CommonName)) %>%
           ggplot(aes(x = reorder(CommonName, -Nrecords_bioblitz), y = Nrecords_bioblitz, fill = maingroup)) +
           geom_bar(stat = "identity") +
           theme_bw() +
           theme(axis.text.x = element_text(angle = 45, hjust = 1),
                 plot.margin = margin(t = .5, r = .5, b = .5, l = 1.5, unit = "cm")) +
           scale_fill_manual(values = as.character(bb_palette %>%
                                                     filter(maingroup %in% dat$maingroup) %>%
                                                     select(col) %>%
                                                     as_vector()),
                             name = "Group") +
           labs(x = "Species",
                y = "No. observations",
                title = "Bioblitz 2021"),
         width = 2000,
         height = 1500, 
         units = "px",
         dpi = 300)
}

### Function to plot common group lollipop plot
Common.group.p <- function(file, path, i, group, count, palette, percent = FALSE, yaxis = FALSE){
  dat <- read_csv(file) %>%
    mutate_if(is.numeric, .funs = function(x) replace_na(x, 0)) %>%
    arrange(desc({{count}}))
  
  sum.count <- sum(dat[,deparse(substitute(count))])
  if(sum.count == 0) print(paste0("All values of ", deparse(substitute(count)), " equal 0!"))
  else {
    if(nrow(dat) > 25)
      dat <-
        dat %>%
        slice(1:25L)
    
    p <- dat %>%
      ggplot(aes(x = reorder({{group}}, {{count}}), y = {{count}})) +
      geom_segment(aes(x = reorder({{group}}, {{count}}), xend = reorder({{group}}, {{count}}), y = 0, yend = {{count}}), color = "gray", lwd = 1) +
      coord_flip() +
      theme_bw() +
      labs(x = "",
           y = "")
    
    if("maingroup" %in% names(dat))
      p <- p +
      geom_point(aes(fill = maingroup), shape = 21, size = 3) +
      scale_fill_manual(values = as.character(palette %>%
                                                filter(maingroup %in% dat$maingroup) %>%
                                                select(col) %>%
                                                as_vector()),
                        name = "Group")
    else
      p <- p +
      geom_point(shape = 21, fill = "blue", size = 3)
    
    if(!percent)
      p <- p + geom_text(aes(label = {{count}}), nudge_y = max(dat[,deparse(substitute(count))])/10, color = "black", size = 3) else
        p <- p + geom_text(aes(label = scales::percent({{count}}/sum.count, accuracy = 0.01)), nudge_y = max(dat[,deparse(substitute(count))])/10, color = "black", size = 3)
    
    if(!yaxis)
      p <- p + theme(axis.text.y = element_blank())
    
    ggsave(paste0(path, "/Common_", deparse(substitute(group)),"_", deparse(substitute(count)), "_", i, ".jpg"), 
           plot = p,
           width = 2000,
           height = 1500, 
           units = "px",
           dpi = 300)
  }
}

### Function to generate cumulative species trend over time plot
Trend.sp <- function(file, path, i, palette){
  dat <- read_csv(file)
  
  p <- dat %>%
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
    scale_color_manual(values = as.character(palette %>%
                                               filter(maingroup %in% dat$maingroup) %>%
                                               select(col) %>%
                                               as_vector()),
                       name = "Group") +
    facet_wrap(facets = vars(maingroup), ncol = 3, scales = "free_y") +
    labs(x = "Year", y = "No. species")
  
  ggsave(paste0(path, "/Cumulative_species_", i, ".jpg"),
         plot = p,
         width = 2000,
         height = 1500, 
         units = "px",
         dpi = 300)
}

### Function to generate observation trend over time plot
Trend.obs <- function(dat, path, i, palette){
  p <- 
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
    scale_color_manual(values = as.character(palette %>%
                                               filter(maingroup %in% filter(dat, !record_period == "historical")$maingroup) %>%
                                               select(col) %>%
                                               as_vector()),
                       name = "Group") +
    facet_wrap(facets = vars(maingroup), ncol = 3, scales = "free_y") +
    labs(x = "Year", y = "No. observations")
  
  ggsave(paste0(path, "/Observations_", i, ".jpg"),
         plot = p,
         width = 2000,
         height = 1500, 
         units = "px",
         dpi = 300)
}

### Function to generate heatmap of observations (either gridded or density)
Map.p <- function(dat, grid_dat, occ_dat, basemap, path, i, period, type = "Nrecords"){
  if(type == "Nrecords") {
    p <- grid_dat %>%
      left_join(occ_dat %>%
                  filter(record_period == period) %>%
                  group_by(FID) %>%
                  summarise(occ = sum(occ))) %>%
      ggplot() +
      geom_sf(data = map_LGA, fill = "white", colour = "dark grey", inherit.aes = F) +
      geom_sf(aes(fill = occ), colour = NA, alpha = .8) +
      scale_fill_viridis_c(na.value = NA, option = "plasma", name = type) +
      theme_bw() +
      coord_sf(xlim = st_bbox(grid_dat)[c(1,3)],
               ylim = st_bbox(grid_dat)[c(2,4)]) +
      ggtitle(str_to_title(period))
  } else if(type == "Ndensity") {
    kde <- st_kde(dat %>%
                    filter(record_period == period) %>%
                    select(locality, decimalLatitude, decimalLongitude, identifiedBy, recordID, record_period, keep, maingroup) %>%
                    st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = 7844) %>%
                    st_transform(crs = 3112))
    
    p <- kde %>%
      ggplot() +
      geom_sf(data = map_LGA, fill = "white", colour = "dark grey", inherit.aes = F) +
      geom_sf(data = kde$sf %>%
                filter(contlabel %in% c(5, 25, 50, 75, 95)) %>%
                mutate(contlabel = factor(contlabel, levels = c(95, 75, 50, 25, 5))),
              aes(fill = label_percent(contlabel)), colour = NA, alpha = .8) +
      scale_fill_viridis_d(na.value = NA, option = "plasma", name = type) +
      theme_bw() +
      coord_sf(xlim = st_bbox(grid_dat)[c(1,3)],
               ylim = st_bbox(grid_dat)[c(2,4)]) +
      ggtitle(str_to_title(period))
  } else stop("type should be `Nrecords` or `Ndensity`")
  
  ggsave(paste0(path, "/", type, "_", i, "_", period, ".jpg"),
         plot = p,
         width = 2000,
         height = 1500,
         units = "px",
         dpi = 300)
}

## Function to generate greenspace map (either heatmap or point locality)
Map.gs <- function(dat, green_dat, occ_dat, basemap, path, i, period, type = "Records"){
  
  if(type == "Records") {
    p <- 
      green_dat %>%
      ggplot() +
      geom_sf(data = map_LGA, fill = "white", colour = "dark grey", inherit.aes = F) +
      geom_sf(fill = "darkolivegreen4", colour = NA, alpha = .8) +
      geom_sf(data = dat %>%
                filter(record_period == period) %>%
                st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = "WGS84"),
              colour = "blue", size = .5) +
      theme_bw() +
      coord_sf(xlim = st_bbox(grid_dat)[c(1,3)],
               ylim = st_bbox(grid_dat)[c(2,4)]) +
      ggtitle("Spatial location of records")
  } else if (type == "Greenspace") {
    p <- 
      green_dat %>%
      left_join(occ_greenspace %>%
                  filter(record_period == period) %>%
                  group_by(FID) %>%
                  summarise(occ = sum(occ))) %>%
      ggplot() +
      geom_sf(data = map_LGA, fill = "white", colour = "dark grey", inherit.aes = F) +
      geom_sf(aes(fill = occ), colour = "light grey", size = .1, alpha = .8) +
      scale_fill_viridis_c(na.value = NA, option = "plasma", name = "Nrecords") +
      theme_bw() +
      coord_sf(xlim = st_bbox(grid_dat)[c(1,3)],
               ylim = st_bbox(grid_dat)[c(2,4)]) +
      ggtitle("Records by urban greenspace")
  } else stop("type should be `Records` or `Greenspace`")
  
  ggsave(paste0(path, "/", type, "_map_", i, "_", period, ".jpg"),
         plot = p,
         width = 2000,
         height = 1500, 
         units = "px",
         dpi = 300)
}



#### GENERAL FUNCTIONS ####

### Function to calculate angles for per-group symbols around pie chart
calc_angles <- function(data){
  data <- as.numeric(data)
  angle_each <- 360/sum(data)
  angle_group <- data * angle_each 
  return(cumsum(angle_group))
}

### Function to draw piechart map with symbols
plot_pie_map <- function(p_map, dat, bb_palette, scale.fact){
  library(grImport2)
  library(rphylopic)
  
  # rsvg::rsvg_png("phylopic/Arachnids.svg", "phylopic/Arachnids.png", height = 128)
  # rsvg::rsvg_png("phylopic/Birds.svg", "phylopic/Birds.png", height = 128)
  # rsvg::rsvg_png("phylopic/Cartilaginous fishes.svg", "phylopic/Cartilaginous fishes.png", height = 128)
  # rsvg::rsvg_png("phylopic/Frogs.svg", "phylopic/Frogs.png", height = 128)
  # rsvg::rsvg_png("phylopic/Fungi.svg", "phylopic/Fungi.png", height = 128)
  # rsvg::rsvg_png("phylopic/Arachnids.svg", "phylopic/Arachnids.png", height = 128)
  # rsvg::rsvg_png("phylopic/Insects.svg", "phylopic/Insects.png", height = 128)
  # rsvg::rsvg_png("phylopic/Lampreys.svg", "phylopic/Lampreys.png", height = 128)
  # rsvg::rsvg_png("phylopic/Mammals.svg", "phylopic/Mammals.png", height = 128)
  # rsvg::rsvg_png("phylopic/Other invertebrates.svg", "phylopic/Other invertebrates.png", height = 128)
  # rsvg::rsvg_png("phylopic/Plants.svg", "phylopic/Plants.png", height = 128)
  # rsvg::rsvg_png("phylopic/Ray-finned fishes.svg", "phylopic/Ray-finned fishes.png", height = 128)
  # rsvg::rsvg_png("phylopic/Reptiles.svg", "phylopic/Reptiles.png", height = 128)
  
  phylopics <- list(Arachnids = png::readPNG("phylopic/Arachnids.png"),
                    Birds = png::readPNG("phylopic/Birds.png"),
                    `Cartilaginous fishes` = png::readPNG("phylopic/Cartilaginous fishes.png"),
                    Frogs = png::readPNG("phylopic/Frogs.png"),
                    Fungi = png::readPNG("phylopic/Fungi.png"),
                    Insects = png::readPNG("phylopic/Insects.png"),
                    Lampreys = png::readPNG("phylopic/Lampreys.png"),
                    Mammals = png::readPNG("phylopic/Mammals.png"),
                    `Other invertebrates` = png::readPNG("phylopic/Other invertebrates.png"),
                    Plants = png::readPNG("phylopic/Plants.png"),
                    `Ray-finned fishes` = png::readPNG("phylopic/Ray-finned fishes.png"),
                    Reptiles = png::readPNG("phylopic/Reptiles.png"))
  
  coords_ALL <- dat[1, c("X", "Y")]
  
  radius_ALL <- as.numeric(p_map$layers[[2]]$data[1, 4])*scale.fact + 3000
  
  angles_ALL <- calc_angles(dat[1, bb_palette$maingroup])
  
  angles_ALL <- (angles_ALL + replace_na(lag(angles_ALL), 0))/2
  
  bb_palette <- bb_palette %>%
    mutate(angle = angles_ALL,
           x = coords_ALL$X + radius_ALL * sin(pi * 2 * angle/360),
           y = coords_ALL$Y + radius_ALL * cos(pi * 2 * angle/360))
  
  groups <- colSums(dat[, bb_palette$maingroup])
  groups <- names(groups[which(groups > 0)])
  
  for (i in groups){
    p_map <- p_map +
      add_phylopic(phylopics[[i]],
                   x = as.numeric(filter(bb_palette, maingroup == i)[,4]),
                   y = as.numeric(filter(bb_palette, maingroup == i)[,5]),
                   ysize = 3000,
                   color = as.character(filter(bb_palette, maingroup == i)[,2]),
                   alpha = 1)
  }
  
  return(p_map)
}