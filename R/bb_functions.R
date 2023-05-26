#### SUMMARY FUNCTIONS ####

### Function to generate summary table of number of records per species
Splist.Nrecords <- function(dat, full_dat, cypher, path, i){
  dat <-
    dat %>%
    group_by(morphospecies, record_period) %>%
    summarise(Nrecords_total = sum(occ)) %>%
    spread(record_period, Nrecords_total) %>%
    mutate_if(is.numeric, replace_na, 0) %>%
    mutate(total = sum(bb22, historical, recent, na.rm = T),
           Perc_records_total = total/nrow(dat)*100,
           Perc_records_historical = historical/total*100,
           Perc_records_recent = recent/total*100,
           Perc_records_bioblitz = bb22/total*100,
           distribution = case_when(Perc_records_bioblitz == 0 ~ 1,
                                    Perc_records_bioblitz == 100 ~ 2,
                                    T ~ 3)) %>%
    ungroup() %>%
    left_join(dat %>%
                group_by(morphospecies) %>%
                summarise(most_recent_record = max(year))) %>%
    left_join(dat %>%
                filter(record_period %in% c("recent", "bb22")) %>%
                group_by(morphospecies) %>%
                summarise(first_record_after1991 = min(year))) %>%
    left_join(cypher, by = c("morphospecies")) %>%
    mutate(lazarus = case_when(historical > 0 & recent == 0 & bb22 > 0 ~ "Lazarus"),
           new = case_when(historical == 0 & recent == 0 & bb22 > 0 ~ "New"),
           new_introduced = case_when(historical == 0 & recent == 0 & bb22 > 0 & origin == "Introduced" ~ "New introduced")) %>%
    rename(CommonName = cname,
           Origin = origin,
           EPBC = epbc,
           FFGAct = ffga,
           Nrecords_total = total,
           Nrecords_historical = historical,
           Nrecords_recent = recent,
           Nrecords_bioblitz = bb22) %>%
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
        summarise(Nrecords_total_from1992 = sum(Nrecords_recent + Nrecords_bioblitz),
                  Nrecords_recent = sum(Nrecords_recent),
                  Nrecords_bioblitz = sum(Nrecords_bioblitz)) %>%
        left_join(dat %>%
                    filter(Nrecords_recent > 0 | Nrecords_bioblitz > 0) %>%
                    group_by({{group}}) %>%
                    summarise(Nsp_total_from1992 = length(morphospecies))) %>%
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
    mutate(Nrecords_total_from1992 = Nrecords_recent + Nrecords_bioblitz,
           Perc_records_bioblitz = (Nrecords_bioblitz/Nrecords_total_from1992)*100) %>%
    arrange(desc(Nrecords_bioblitz)) %>%
    slice(1:25) %>%
    select(morphospecies, CommonName, genus, family, subgroup, maingroup, Origin, Nrecords_bioblitz, Nrecords_total_from1992, Perc_records_bioblitz)
  
  write_csv(dat, paste0(path, "/Common_species_", i, ".csv"), na = "")
}

### Function to find most common group per grouping variable
Common.group <- function(file, path, group, i){
  dat <- read_csv(file) %>%
    ungroup()
  
  dat <- dat %>%
    mutate(Perc_records_bioblitz = (Nrecords_bioblitz/Nrecords_total_from1992)*100,
           Nuniquesp_bioblitz = replace_na(Nuniquesp_bioblitz, 0),
           Perc_newsp_bioblitz = (Nuniquesp_bioblitz/Nsp_total_from1992)*100,
           Perc_records_bioblitz = replace_na(Perc_records_bioblitz, 0),
           Perc_newsp_bioblitz = replace_na(Perc_newsp_bioblitz, 0)) %>%
    arrange(desc(Nrecords_bioblitz))
  
  if(deparse(substitute(group)) == "maingroup")
    dat <- dat %>%
    select({{group}}, Nrecords_total_from1992, Nrecords_bioblitz, Perc_records_bioblitz, Nuniquesp_bioblitz, Perc_newsp_bioblitz, Nsp_bioblitz) else
      dat <- dat %>%
    select(maingroup, {{group}}, Nrecords_total_from1992, Nrecords_bioblitz, Perc_records_bioblitz, Nuniquesp_bioblitz, Perc_newsp_bioblitz, Nsp_bioblitz)
  
  write_csv(dat, paste0(path, "/Common_", deparse(substitute(group)), "_", i, ".csv"), na = "")
}

### Function to generate summary table by green space
Nrecords.GS <- function(dat, occ_dat, sp, path, i){
  dat <-
    dat %>%
    as_tibble() %>%
    select(FID, PARK_NAME) %>%
    drop_na() %>%
    left_join(occ_dat %>% filter(record_period == "bb22")) %>%
    drop_na(maingroup, PARK_NAME) %>%
    group_by(PARK_NAME, maingroup, record_period) %>%
    summarise(occ = sum(occ)) %>%
    mutate(maingroup = case_when(maingroup == "Arachnids" ~ "Nrecords_ARA",
                                 maingroup == "Birds" ~ "Nrecords_BIR",
                                 maingroup == "Cartilaginous fishes" ~ "Nrecords_CAF",
                                 maingroup == "Chromists" ~ "Nrecords_CHR",
                                 maingroup == "Frogs" ~ "Nrecords_FRO",
                                 maingroup == "Fungi" ~ "Nrecords_FUN",
                                 maingroup == "Insects" ~ "Nrecords_INS",
                                 maingroup == "Other invertebrates" ~ "Nrecords_INV",
                                 maingroup == "Mammals" ~ "Nrecords_MAM",
                                 maingroup == "Plants" ~ "Nrecords_PLA",
                                 maingroup == "Ray-finned fishes" ~ "Nrecords_RAF",
                                 maingroup == "Reptiles" ~ "Nrecords_REP")) %>%
    spread(maingroup, occ) %>%
    rename(Greenspace_name = PARK_NAME) %>%
    mutate_if(is.numeric, function(x) replace_na(x, 0)) %>%
    rowwise() %>%
    mutate(Nrecords_total = sum(across(where(is.numeric)))) %>%
    left_join(dat %>%
                as_tibble() %>%
                select(FID, PARK_NAME) %>%
                drop_na() %>%
                left_join(occ_dat %>% filter(record_period == "bb22")) %>%
                drop_na(maingroup, PARK_NAME) %>%
                group_by(PARK_NAME, maingroup, record_period) %>%
                summarise(occ = sum(occ)) %>%
                mutate(maingroup = case_when(maingroup == "Arachnids" ~ "Nsp_ARA",
                                             maingroup == "Birds" ~ "Nsp_BIR",
                                             maingroup == "Cartilaginous fishes" ~ "Nsp_CAF",
                                             maingroup == "Chromists" ~ "Nsp_CHR",
                                             maingroup == "Frogs" ~ "Nsp_FRO",
                                             maingroup == "Fungi" ~ "Nsp_FUN",
                                             maingroup == "Insects" ~ "Nsp_INS",
                                             maingroup == "Other invertebrates" ~ "Nsp_INV",
                                             maingroup == "Mammals" ~ "Nsp_MAM",
                                             maingroup == "Plants" ~ "Nsp_PLA",
                                             maingroup == "Ray-finned fishes" ~ "Nsp_RAF",
                                             maingroup == "Reptiles" ~ "Nsp_REP")) %>%
                spread(maingroup, occ) %>%
                rename(Greenspace_name = PARK_NAME) %>%
                mutate_if(is.numeric, function(x) replace_na(x, 0)) %>%
                mutate_if(is.numeric, function(x) ifelse(x > 0, 1, 0)) %>%
                rowwise() %>%
                mutate(Nsp_total = sum(across(where(is.numeric)))))
  
  write_csv(dat, paste0(path, "/", i, "bb_Nrecords_GS_groups.csv.csv"), na = "") 
}

### Function to generate summary table of participants
Nrecords.PAR <- function(dat, path, i){
  dat <-
    dat %>%
    filter(record_period == "bb22") %>%
    group_by(identifiedBy) %>%
    summarise(Nrecords_bb22 = sum(occ)) %>%
    rename(Participant_username = identifiedBy) %>%
    mutate(Perc_records_bb22 = (Nrecords_bb22/sum(Nrecords_bb22))*100) %>%
    arrange(desc(Nrecords_bb22))
  
  write_csv(dat, paste0(path, "/BBparticipants_", i, ".csv"), na = "")
}

### Function to generate summary table of participant types
BBparticipants_type <- function(dat, participants, path, i){
  dat <-
    dat %>%
    filter(record_period == "bb22") %>%
    left_join(participants, by = c("identifiedBy" = "Participant_username"))
  
  dat_Nrecords <- bind_rows(dat %>%
                              count(Participant_type) %>%
                              drop_na(),
                            dat %>%
                              count(BBorganiser_type) %>%
                              drop_na()) %>%
    rename(Nrecords_BB = n)
  
  dat_Nsp <- bind_rows(dat %>%
                         select(Participant_type, morphospecies) %>%
                         unique() %>%
                         count(Participant_type) %>%
                         drop_na(),
                       dat %>%
                         select(BBorganiser_type, morphospecies) %>%
                         unique() %>%
                         count(BBorganiser_type) %>%
                         drop_na()) %>%
    rename(Nsp_BB = n)
  
  dat <- full_join(dat_Nrecords, dat_Nsp) %>%
    select(Participant_type, BBorganiser_type, Nrecords_BB, Nsp_BB)
  
  write_csv(dat, paste0(path, "/BBparticipants_type_", i, ".csv"), na = "")
}

### Function to generate summary table of binned participants
BBparticipants_bin <- function(dat, path, i){
  dat <-
    dat %>%
    filter(record_period == "bb22") %>%
    group_by(identifiedBy) %>%
    summarise(Nrecords_bb22 = sum(occ)) %>%
    count(Nrecords_bin = cut(Nrecords_bb22, breaks = c(0, 5, 50, 200, max(Nrecords_bb22)+1),
                             labels = c("1 to 5",
                                        "6 to 50",
                                        "51 to 200",
                                        paste0("200 to ", max(Nrecords_bb22)+1)))) %>%
    rename(Nparticipants = n)
  
  write_csv(dat, paste0(path, "/BBparticipants_bins_", i, ".csv"), na = "")
}

### Function to generate summary table of main groups
BBmaingroup <- function(dat, path, i){
  dat <-
    dat %>%
    filter(record_period == "bb22") %>%
    mutate(maingroup = factor(maingroup, levels =c("Arachnids","Birds","Frogs","Fungi","Insects","Mammals","Other invertebrates","Plants","Ray-finned fishes","Reptiles"))) %>%
    group_by(maingroup)
  
  dat <- left_join(dat %>%
                     select(maingroup, identifiedBy) %>%
                     unique() %>%
                     count() %>%
                     rename(Nparticipants_bb22 = n),
                   dat %>%
                     count() %>%
                     rename(Nrecords_bb22 = n))
  
  write_csv(dat, paste0(path, "/BBmaingroup_", i, ".csv"), na = "")
}



#### PLOTTING FUNCTIONS ####

### Function to plot common species barplot
Common.sp.p <- function(file, path, i, palette){
  dat <- read_csv(file)
  
  dat$CommonName[which(dat$CommonName %in% dat$CommonName[duplicated(dat$CommonName)])] <- NA
  
  ggsave(paste0(path, "/Common_species_", i, ".jpg"),
         plot = dat %>%
           mutate(CommonName = case_when(is.na(CommonName) ~ glue::glue("<i>{morphospecies}</i>"),
                                         T ~ CommonName)
                  ) %>%
           ggplot(aes(x = reorder(CommonName, -Nrecords_bioblitz), y = Nrecords_bioblitz, fill = maingroup)) +
           geom_bar(stat = "identity") +
           theme_bw() +
           theme(axis.text.x = ggtext::element_markdown(angle = 45, hjust = 1),
                 plot.margin = margin(t = .5, r = .5, b = .5, l = 1.5, unit = "cm")) +
           scale_fill_manual(values = as.character(bb_palette %>%
                                                     filter(maingroup %in% dat$maingroup) %>%
                                                     select(col) %>%
                                                     as_vector()),
                             name = "Group") +
           labs(x = "Species",
                y = "No. observations",
                title = "Bioblitz 2022"),
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

### Function to plot participant type lollipop plot
Participant.p <- function(file, path, i){
  dat <- read_csv(file) %>%
    mutate(y = case_when(Participant_type == "General public" ~ Participant_type,
                         BBorganiser_type == "Council staff" ~ "BB organiser - Council",
                         BBorganiser_type == "Members of supporting associations" ~ "BB organiser - Supporting")) %>%
    select(y, Nrecords_BB, Nsp_BB) %>%
    drop_na() %>%
    gather(type, count, -y) %>%
    group_by(type) %>%
    mutate(perc = (count / sum(count)) * 100,
           type = factor(type,
                         levels = c("Nrecords_BB", "Nsp_BB"),
                         labels = c("Observations", "Species"))) %>%
    arrange(desc(perc))
  
  p <- dat %>%
    ggplot(aes(x = y, y = perc, colour = type)) +
    geom_linerange(aes(x = y, ymin = 0, ymax = perc), lwd = 1, position = position_dodge2(width = 1, reverse = T)) +
    coord_flip() +
    theme_bw() +
    ylim(0, 100) +
    labs(x = "",
         y = "% of records",
         colour = "") +
    scale_colour_manual(values = c("deepskyblue2", "red")) +
    geom_point(shape = 21, fill = "purple", size = 3, position = position_dodge2(width = 1, reverse = T)) + 
    geom_text(aes(label = count, group = type, y = perc + 5),
              colour = "black", size = 3,
              position = position_dodge2(width = 1, reverse = T)) +
    theme(legend.position = "bottom")
  
  ggsave(paste0(path, "/BBparticipant_type_", deparse(substitute(count)), "_", i, ".jpg"), 
         plot = p,
         width = 2000,
         height = 1500, 
         units = "px",
         dpi = 300)
}

### Function to generate cumulative species trend over time plot
Trend.sp <- function(file, path, i, palette){
  dat <- read_csv(file)
  
  p <- dat %>%
    group_by(first_record_after1991, maingroup) %>%
    count() %>%
    mutate(first_record_after1991 = case_when(is.na(first_record_after1991) ~ 1991,
                                              T ~ first_record_after1991)) %>%
    arrange(first_record_after1991) %>%
    ungroup() %>%
    group_by(maingroup) %>%
    mutate(n = cumsum(n)) %>%
    ungroup() %>%
    spread(maingroup, n) %>%
    mutate_all(funs(case_when(first_record_after1991 == 2022 ~ max(., na.rm = T), 
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
Map.p <- function(dat, grid_dat, occ_dat, basemap, path, i, period, type = "Nrecords", palette = NULL){
  if(type == "Nrecords") {
    p <- grid_dat %>%
      left_join(occ_dat %>%
                  filter(record_period == period) %>%
                  group_by(FID) %>%
                  summarise(occ = sum(occ))) %>%
      ggplot() +
      geom_sf(data = map_VIC, fill = "grey95", colour = "grey75") +
      geom_sf(data = map_LGA, fill = NA, colour = "black", size = 2) +
      geom_sf(aes(fill = occ), colour = NA, alpha = .8) +
      theme_bw() +
      coord_sf(xlim = st_bbox(grid_dat)[c(1,3)],
               ylim = st_bbox(grid_dat)[c(2,4)]) +
      theme(panel.background = element_rect(fill = 'lightblue'),
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1))
    
    if(is.null(palette))
      p <- p + scale_fill_viridis_c(trans = "log10", na.value = NA, option = "plasma", name = "No. of records") else
        p <- p + scale_fill_gradientn(trans = "log10", na.value = NA, colours = palette, name = "No. of records")
      
      ggsave(paste0(path, "/", type, "_", i, "_", period, ".jpg"),
             plot = p,
             width = 2000,
             height = 1500,
             units = "px",
             dpi = 300)
  } else if(type == "Ndensity") {
    kde <- tryCatch(st_kde(dat %>%
                             filter(record_period == period) %>%
                             select(decimalLatitude, decimalLongitude, identifiedBy, recordID, record_period, maingroup) %>%
                             st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = 7844) %>%
                             st_transform(crs = 3112)),
                    error = function(e) e)
    
    if(!inherits(kde, "error")){
      p <- kde %>%
        ggplot() +
        geom_sf(data = map_VIC, fill = "grey95", colour = "grey75") +
        geom_sf(data = map_LGA, fill = NA, colour = "black", size = 2) +
        geom_sf(data = kde$sf %>%
                  filter(contlabel %in% c(5, 25, 50, 75, 95)) %>%
                  mutate(contlabel = factor(contlabel, levels = c(95, 75, 50, 25, 5))),
                aes(fill = label_percent(contlabel)), colour = NA, alpha = .8) +
        theme_bw() +
        coord_sf(xlim = st_bbox(grid_dat)[c(1,3)],
                 ylim = st_bbox(grid_dat)[c(2,4)]) +
        ggtitle(str_to_title(period)) +
        theme(panel.background = element_rect(fill = 'lightblue'),
              panel.grid = element_blank(),
              axis.text.x = element_text(angle = 45, hjust = 1))
      
      if(is.null(palette))
        p <- p + scale_fill_viridis_d(na.value = NA, option = "plasma", name = "Density of records") else
          p <- p + scale_fill_manual(na.value = NA, values = palette, name = "Density of records")
        
        ggsave(paste0(path, "/", type, "_", i, "_", period, ".jpg"),
               plot = p,
               width = 2000,
               height = 1500,
               units = "px",
               dpi = 300)
    }
  } else stop("type should be `Nrecords` or `Ndensity`")
}

## Function to generate greenspace map (either heatmap or point locality)
Map.gs <- function(dat, green_dat, occ_dat, basemap, path, i, period, type = "Records", palette = NULL){
  
  if(type == "Records") {
    if(is.null(palette))
      col_rec = "blue" else
        col_rec = tail(palette, 1)
    
    p <- 
      green_dat %>%
      ggplot() +
      geom_sf(data = map_VIC, fill = "grey95", colour = "grey75") +
      geom_sf(data = greenspaces, fill = "darkolivegreen", colour = NA, alpha = .8) +
      geom_sf(data = map_LGA, fill = NA, colour = "black", size = 2) +
      geom_sf(data = dat %>%
                filter(record_period == period) %>%
                st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = "WGS84"),
              fill = col_rec, shape = 21, colour = "black", size = 1, stroke = .5) +
      theme_bw() +
      coord_sf(xlim = st_bbox(grid_dat)[c(1,3)],
               ylim = st_bbox(grid_dat)[c(2,4)]) +
      theme(panel.background = element_rect(fill = 'lightblue'),
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1))
  } else if (type == "Greenspace") {
    p <- 
      green_dat %>%
      left_join(occ_dat %>%
                  filter(record_period == period) %>%
                  group_by(FID) %>%
                  summarise(occ = sum(occ))) %>%
      ggplot() +
      geom_sf(data = map_VIC, fill = "grey95", colour = "grey75") +
      geom_sf(aes(fill = occ), colour = "light grey", size = .1, alpha = .8) +
      geom_sf(data = map_LGA, fill = NA, colour = "black", size = 2) +
      theme_bw() +
      coord_sf(xlim = st_bbox(grid_dat)[c(1,3)],
               ylim = st_bbox(grid_dat)[c(2,4)]) +
      theme(panel.background = element_rect(fill = 'lightblue'),
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1))
    
    if(is.null(palette))
      p <- p + scale_fill_viridis_c(na.value = NA, option = "plasma", name = "No. of records") else
        p <- p + scale_fill_gradientn(na.value = NA, colours = palette, name = "No. of records")
      
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
plot_pie_map <- function(p_map, dat, palette, scale.fact){
  require(grImport2)
  require(rphylopic)
  
  # rsvg::rsvg_png("phylopic/Arachnids.svg", "phylopic/Arachnids.png", height = 128)
  # rsvg::rsvg_png("phylopic/Birds.svg", "phylopic/Birds.png", height = 128)
  # rsvg::rsvg_png("phylopic/Cartilaginous fishes.svg", "phylopic/Cartilaginous fishes.png", height = 128)
  # rsvg::rsvg_png("phylopic/Chromists.svg", "phylopic/Chromists.png", height = 128)
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
                    Chromists = png::readPNG("phylopic/Chromists.png"),
                    Frogs = png::readPNG("phylopic/Frogs.png"),
                    Fungi = png::readPNG("phylopic/Fungi.png"),
                    Insects = png::readPNG("phylopic/Insects.png"),
                    Mammals = png::readPNG("phylopic/Mammals.png"),
                    `Other invertebrates` = png::readPNG("phylopic/Other invertebrates.png"),
                    Plants = png::readPNG("phylopic/Plants.png"),
                    `Ray-finned fishes` = png::readPNG("phylopic/Ray-finned fishes.png"),
                    Reptiles = png::readPNG("phylopic/Reptiles.png"))
  
  coords_ALL <- dat[1, c("X", "Y")]
  
  radius_ALL <- as.numeric(p_map$layers[[2]]$data[1, 4])*scale.fact + 3000
  
  angles_ALL <- calc_angles(dat[1, palette$maingroup])
  
  angles_ALL <- (angles_ALL + replace_na(lag(angles_ALL), 0))/2
  
  palette <- palette %>%
    mutate(angle = angles_ALL,
           x = coords_ALL$X + radius_ALL * sin(pi * 2 * angle/360),
           y = coords_ALL$Y + radius_ALL * cos(pi * 2 * angle/360))
  
  groups <- colSums(dat[, palette$maingroup])
  groups <- names(groups[which(groups > 0)])
  
  for (i in groups){
    p_map <- p_map +
      add_phylopic(phylopics[[i]],
                   x = as.numeric(filter(palette, maingroup == i)[,4]),
                   y = as.numeric(filter(palette, maingroup == i)[,5]),
                   ysize = 3000,
                   color = as.character(filter(palette, maingroup == i)[,2]),
                   alpha = 1)
  }
  
  return(p_map)
}

### Function to calculate "enrichment" of species in category
calc_enrich <- function(venn_groups, group, lower.tail = T, plot = F){
  sp_ABC <- venn_groups[7]
  
  if(group == "General public") {
    sp_A <- venn_groups[1]
    sp_B <- venn_groups[2]
    sp_C <- venn_groups[3]
    sp_AB <- venn_groups[4]
    sp_BC <- venn_groups[6]
    sp_AC <- venn_groups[5]
  }
  
  if(group == "BB organiser - Council") {
    sp_A <- venn_groups[2]
    sp_B <- venn_groups[3]
    sp_C <- venn_groups[1]
    sp_AB <- venn_groups[6]
    sp_BC <- venn_groups[5]
    sp_AC <- venn_groups[4]
  }
  
  if(group == "BB organiser - Supporting") {
    sp_A <- venn_groups[3]
    sp_B <- venn_groups[1]
    sp_C <- venn_groups[2]
    sp_AB <- venn_groups[5]
    sp_BC <- venn_groups[4]
    sp_AC <- venn_groups[6]
  }
  
  overlap = sp_AB + sp_AC + sp_ABC # size of overlap
  group1 = sp_A + overlap # number of species in group A
  group2 = sum(venn_groups) - sp_A # number of species in groups B + C
  total = sum(venn_groups) # total number of species
  
  if(lower.tail)
    p <- phyper(overlap, group2, total - group2, group1, lower.tail = T) else
      p <- phyper(overlap - 1, group2, total - group2, group1, lower.tail = F)
  
  if(plot){
    d <- dhyper(0:sum(venn_groups), group2, total - group2, group1)
    
    d_plot_lims <- c(min(seq(0, 100, length = length(d))[min(which(d > 0))], (overlap/sum(venn_groups))*100),
                     max(seq(0, 100, length = length(d))[max(which(d > 0))], (overlap/sum(venn_groups))*100))
    
    pd <- tibble(Overlap = 0:sum(venn_groups),
                 Density = d) %>%
      mutate(Overlap = (Overlap / max(Overlap)) * 100) %>%
      ggplot(aes(x = Overlap, y = Density)) +
      geom_area(data = filter(tibble(Overlap = 0:sum(venn_groups),
                                     Density = d),
                              Overlap <= overlap) %>%
                  mutate(Overlap = (Overlap / sum(venn_groups)) * 100),
                fill = "red", alpha = .2, colour = NA) +
      geom_line(colour = "grey") +
      geom_vline(aes(xintercept = (overlap/sum(venn_groups))*100), colour = "deepskyblue2") +
      theme_bw() +
      scale_x_continuous(limits = c(max(d_plot_lims[1] - 15, 0),
                                    min(d_plot_lims[2] + 15, 100))) +
      labs(x = "% Overlap",
           title = paste("Probability density of overlap: ", group, sep = ""))
    
    print(pd)
    
    return(list(p = p, plot = pd))
  } else
    return(p)
}
