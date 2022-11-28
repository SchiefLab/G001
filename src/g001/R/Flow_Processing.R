#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE);
suppressPackageStartupMessages({
  #' load all required libraries.
  require(CytoML)
  require(flowWorkspace)
  require(cytoqc)
  require(tidyr)
  require(ggcyto)
  require(assertthat)
  require(dplyr)
  require(tidyr)
  require(ggrepel)
  require(stringr)
  require(purrr)
  require(lubridate)
})
if (!length(args) %in% 4:5) {
  stop("invalid number of arguments to flow processing script.")
}

source('src/g001/R/Flow_Functions.R')

key = args[1]
if (!key %in% c('vrc','fhrc')) {
  stop('First argument must be vrc or fhrc')
}
manifest = args[2]
flow_path = args[3]
override = ifelse(args[4] == 'yes', TRUE, FALSE)
output_path = file.path(
  ifelse(is.na(args[5]), file.path('..', 'flow_results'), args[5]),
  key
  )



# Saving output and messages to log file
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)

if (!interactive()) {
  # Adding timestamp to log so it stops overriding itself
  log_path <- file.path(output_path, paste0(key, '_log_', format(Sys.time(), "%Y%m%d_%H%M%S"),'.txt'))
  log_file <- file(log_path)
  sink(log_file,  split = T)
  sink(log_file,  type = 'message')
}

# Reading Manifest ----------------------------------------------------------------------------

print(paste("Processing",key,"flow data"))
print(paste("Reading",key,"manifest"))
read_manifest <- function(x){
  read.csv(x, na.strings = c(NA,''))
}

flow_manifest <- as.data.frame(read_manifest(manifest))
#check the columns
if (!all(
  c(
    "PTID",
    "Visit",
    "EXPERIMENT_NAME",
    "File",
    "FileType",
    "INX_Plate",
    "INX_Population",
    "KO_Probe",
    "IgG_Gate",
    "Tissue_State",
    "INX_Number_Of_Cells",
    "Note"
  ) %in% colnames(flow_manifest)
)) {
  stop(paste("Not all expected columns found in the", key, "flow manifest"))
}

#standardize the columns
flow_manifest$PTID <- as.character(flow_manifest$PTID)
flow_manifest$Visit <- as.character(flow_manifest$Visit)
df <- as.data.frame(stringr::str_split_fixed(flow_manifest$Visit,"_",2), stringsAsFactors = FALSE)
colnames(df) <- c("Visit","Tube")
flow_manifest$Visit <- df$Visit
flow_manifest$Tube <- as.character(df$Tube)
flow_manifest$EXPERIMENT_NAME <- as.character(flow_manifest$EXPERIMENT_NAME)
flow_manifest$File <- as.character(flow_manifest$File)
flow_manifest$FileType <- as.character(flow_manifest$FileType)
flow_manifest$INX_Plate <- as.character(flow_manifest$INX_Plate)
flow_manifest$INX_Population <- as.character(flow_manifest$INX_Population)
flow_manifest$KO_Probe <- as.character(flow_manifest$KO_Probe)
flow_manifest$IgG_Gate <- as.character(flow_manifest$IgG_Gate)
flow_manifest$Tissue_State <- as.character(flow_manifest$Tissue_State)
flow_manifest$INX_Number_Of_Cells <- as.numeric(flow_manifest$INX_Number_Of_Cells)

#Check for dup filenames but where one row is missing ptid/visit
dup_files <- flow_manifest %>% dplyr::count(File) %>% filter(n == 2) %>% pull('File')
flow_manifest <- flow_manifest %>% filter(!(File %in% dup_files & is.na(PTID)))


#Check that the only files with no visit are the xml files
#If not then infer the visit from the file name and record that change
files_with_no_visit <-
  flow_manifest %>% filter(is.na(Visit), FileType != ".xml")
if (nrow(files_with_no_visit) > 0) {
  cat("The following files have no visit in the manifest:\n")
  inferred_visits <- files_with_no_visit %>% rowwise() %>% do({
    vinf <-
      stringr::str_split(.$File, "_")[[1]][grepl("V\\d\\dA*", stringr::str_split(.$File, "_")[[1]])]
    cat(paste(.$File, " inferred ", vinf, "\n"))

    data.frame(File = .$File, Visit = vinf)
  })

  for (i in 1:nrow(inferred_visits)) {
    flow_manifest[flow_manifest$File == inferred_visits$File[i], "Visit"] = inferred_visits$Visit[i]
  }
}


# Fixing bad tissue
bad_ids_visits <-
  flow_manifest %>%
  filter((Visit %in% c("V02","V06","V07","V08","V10") & Tissue_State == 'Fresh') |
           (Visit %in% c("V05","V07A","V09") & Tissue_State == 'Frozen')) %>%
  distinct(PTID, Visit) %>%
  mutate(PTID_Visit = paste0(PTID , '(', Visit, ')')) %>%
  pull(PTID_Visit)

if (length(bad_ids_visits) > 0) {
  cat("The following samples have incorrect tissue state in the manifest (will infer based on visit number):\n")
  cat(paste0(bad_ids_visits, '\n'))

  flow_manifest <-
    flow_manifest %>%
    mutate(Tissue_State = case_when(
      Visit %in% c("V02","V06","V07","V08","V10") & Tissue_State == 'Fresh' ~ 'Frozen',
      Visit %in% c("V05","V07A","V09") & Tissue_State == 'Frozen' ~ 'Fresh',
      TRUE ~ Tissue_State)
    )
}

# Filling in some missing PTIDs for csv files
missing_csv_ptids <- flow_manifest %>%
  filter(is.na(PTID) & FileType == '.csv') %>%
  distinct(EXPERIMENT_NAME, Visit) %>%
  mutate(EXPERIMENT_NAME_Visit = paste0(EXPERIMENT_NAME , '(', Visit, ')')) %>%
  arrange(EXPERIMENT_NAME_Visit) %>%
  pull(EXPERIMENT_NAME_Visit)

if (length(missing_csv_ptids) > 0) {
  cat("The following csv files have missing PTID values (will infer based on other experiment files in manifest):\n")
  cat(paste0(missing_csv_ptids, '\n'))

  flow_manifest <- flow_manifest %>%
    group_by(EXPERIMENT_NAME) %>%
    mutate(PTID = ifelse(is.na(PTID) & FileType == '.csv',
                         unique(na.omit(PTID)), PTID)) %>%
    ungroup()
}

#Check that the only files with no index plate are the index fcs files
files_with_no_inx_plate <- flow_manifest %>% filter(is.na(INX_Plate), grepl("INX",File))
if (nrow(files_with_no_inx_plate) > 0) {
  cat("The following files have no INX Plate despite being INX files in the manifest:\n")
  inferred_visits <- files_with_no_inx_plate %>% rowwise() %>% do({
    cat(paste(.$File,"\n"))
    data.frame()
  })
}

#List incomplete experiments
incomplete_experiments <- flow_manifest %>% group_by(EXPERIMENT_NAME) %>% count(FileType) %>%
  spread(FileType,n) %>% ungroup() %>% replace(is.na(.), 0) %>% filter(`.csv` == 0 | `.xml` == 0 | `.fcs` == 0)
if (nrow(incomplete_experiments) > 0) {
  cat("The following experiments are incomplete, are missing csv, xml, or fcs files\n")
  print(incomplete_experiments$EXPERIMENT_NAME)
}

#summarize INX_Population
cat("INX_Population has the following entries\n")
print(table(flow_manifest$INX_Population))
if (nlevels(factor(flow_manifest$INX_Population)) > 2) {
  cat("Too many levels for INX_Population")
  stop("quitting")
}

#summarize INX_Plate
cat("INX_Plate has the following entries\n")
print(table(flow_manifest$INX_Plate))
if (!all(grepl("P\\d\\d",levels(factor(flow_manifest$INX_Plate))))) {
  cat("Some invalid INX_Plate entries.")
  stop("quitting")
}

#summarize KO_Probe
cat("KO_Probe has the following entries\n")
print(table(flow_manifest$KO_Probe))
if (nlevels(factor(flow_manifest$KO_Probe)) > 2 | !all(c("BB515", "BV570") %in% levels(factor(flow_manifest$KO_Probe)))) {
  cat("Too many or invalid levels for KO_Probe")
  stop("quitting")
}


#summarize IgG_Gate
cat("IgG_Gate has the following entries\n")
print(table(flow_manifest$IgG_Gate))
if (nlevels(factor(flow_manifest$IgG_Gate)) > 2 | !all(c("IgD-", "IgG+") %in% levels(factor(flow_manifest$IgG_Gate)))) {
  cat("Too many or invalid levels for IgG_Gate")
  stop("quitting")
}


#summarize Tissue_State
cat("Tissue_State has the following entries\n")
print(table(flow_manifest$Tissue_State))
if (nlevels(factor(flow_manifest$Tissue_State)) > 2 | !all(c("Fresh", "Frozen") %in% levels(factor(flow_manifest$Tissue_State)))) {
  cat("Too many or invalid levels for Tissue_State")
  stop("quitting")
}



## Next check for missing files.
print(paste0("searching ",flow_path))
files <- list.files(path = flow_path, recursive = TRUE, full.names = TRUE)
if (length(files) == 0) {
  cat("No Files in", flow_path)
  stop("quitting")
}


files_df <- tibble(Path = files, basename = basename(files)) %>%
  mutate(EXPERIMENT_NAME = basename(dirname(files)))
if (any(!(flow_manifest$File %in% files_df$basename))) {
  cat("Some files in manifest not found (might be running a subset of data):\n")
  print(flow_manifest$File[!(flow_manifest$File %in% files_df$basename)])
}



# Adding path to manifest and writing out and saving
flow_manifest_w_paths <- flow_manifest %>%
  rename(Path_short = Path) %>%
  left_join(files_df, by = c('File' = 'basename', 'EXPERIMENT_NAME')) %>%
  # dropping junk file type
  filter(FileType %in% c('.fcs','.xml','.csv'))
if (any(duplicated(flow_manifest_w_paths$File))) {
  cat("Warning: The following file names in the manifest links to multiple files:",
      paste0(unique(flow_manifest_w_paths$File[duplicated(flow_manifest_w_paths$File)]), collapse = ', '), '\n')
}

# Need to do loose matching for files that didn't get correct EXPERIMENT_NAME through file path
# This should only impact csv files
if (any(is.na(flow_manifest_w_paths$Path))) {
  cat(sum(is.na(flow_manifest_w_paths$Path)),
      "file(s) could not link with EXPERIMENT_NAME so trying linking only based on File name\n")
  loose_linking <- flow_manifest_w_paths %>%
    filter(is.na(Path)) %>%
    select(-Path) %>%
    left_join(files_df %>%  select(-EXPERIMENT_NAME), by = c('File' = 'basename')) %>%
    rename(Path_Loose_Matching = Path)
  if (any(duplicated(loose_linking$File))) {
    cat("The following file names in the manifest links to multiple files using loose matching with only File name: ",
        paste0(unique(loose_linking$File[duplicated(loose_linking$File)]), collapse = ', \n'))
    stop("quitting")
  }
  # Adding loose linking path into manifest
  flow_manifest_w_paths$Path[is.na(flow_manifest_w_paths$Path)] <-
    loose_linking$Path_Loose_Matching
}

# Only want to run selected data we have paths for
flow_manifest_w_paths <- flow_manifest_w_paths %>%
  filter(!is.na(Path))
if (nrow(flow_manifest_w_paths) == 0)
  stop('No files that have a valid path')


# Only want to add updated entries, but not override
if (file.exists(file.path(output_path, paste0(key, '_manifest_w_paths.csv')))) {
  flow_manifest_old <- readr::read_csv(file = file.path(output_path, paste0(key, '_manifest_w_paths.csv'))) %>%
    mutate(Tube = as.character(Tube))

  flow_manifest_out <- bind_rows(
    flow_manifest_old %>% anti_join(flow_manifest_w_paths,
                                    by = c("PTID", "EXPERIMENT_NAME", "File")),
    flow_manifest_w_paths
  )
} else {
  flow_manifest_out <- flow_manifest_w_paths
}

if (!interactive()) {
  readr::write_csv(x = flow_manifest_out, file = file.path(output_path, paste0(key, '_manifest_w_paths.csv')))
}


# Getting number of experiments (total and to run)
# Only want manifest entries that are listed as OK in the manifest
if (!any(names(flow_manifest_w_paths) == 'Additional.Notes'))
  flow_manifest_w_paths$`Additional.Notes` = NA

manifest_to_run <- flow_manifest_w_paths %>%
  filter(Note %in% c('OK', 'OK.') |
           (Note %>% str_detect("EXPERIMENT_NAME doesn't match") &
              !Note %>% str_detect("Experiment folder missing csv, fcs, or xml files.") &
              `Additional.Notes` == 'RE-EXPORTED') |
           `Additional.Notes` == 'Confirmed Correct') %>%
  group_by(EXPERIMENT_NAME) %>%
  mutate(All_3_File_Types = if_else(any(FileType == '.xml') &
                                      any(FileType == '.fcs') &
                                      any(FileType == '.csv'),
                                    TRUE, FALSE)) %>%
  # Need to have csv and fcs files for each exp, visit , and ptid
  group_by(EXPERIMENT_NAME, Visit, PTID) %>%
  mutate(All_2_File_Types_visit_ptid = if_else(any(FileType == '.fcs') &
                                                 any(FileType == '.csv'),
                                               TRUE, FALSE)) %>%
  ungroup()

# Need csv, xml, and fcs file types to run experiment
if (any(!manifest_to_run$All_3_File_Types)) {
  cat("The following experiments are missing at least one type of file (csv, xml, or fcs):\n",
      manifest_to_run %>%
        filter(!All_3_File_Types) %>% select(EXPERIMENT_NAME) %>%
        distinct() %>% pull() %>% paste(collapse = ', \n'),'\n'
  )
  manifest_to_run <- manifest_to_run %>% filter(All_3_File_Types)
}

# Need csv, and fcs file types for each experiment, visit, and id
if (manifest_to_run %>% filter(!All_2_File_Types_visit_ptid & FileType != '.xml') %>% nrow() > 0) {
  cat("The following experiments are missing at least one type of file (csv, xml, or fcs):\n",
      manifest_to_run %>%
        filter(!All_2_File_Types_visit_ptid & FileType != '.xml') %>%
        unite('id', EXPERIMENT_NAME, Visit, PTID, sep = '/') %>%
        pluck('id') %>% paste(collapse = ', \n')
  )
  manifest_to_run <- manifest_to_run %>%
    group_by(EXPERIMENT_NAME, Visit, PTID) %>%
    filter(all(All_2_File_Types_visit_ptid | FileType == '.xml')) %>%
    ungroup()
}


# Getting some values of what to run
total_exp_n <- flow_manifest_w_paths %>% distinct(EXPERIMENT_NAME) %>% count()
run_exp <- unique(manifest_to_run$EXPERIMENT_NAME)
run_exp_n <- length(run_exp)

# Pasting key variables that make up an analysis (experiment, visit, ptid)
exp_id_vis <- flow_manifest_w_paths %>%
  filter(FileType == '.fcs') %>%
  select(EXPERIMENT_NAME, Visit, PTID) %>%
  distinct()

run_exp_id_vis <- manifest_to_run %>%
  filter(FileType == '.fcs') %>%
  select(EXPERIMENT_NAME, Visit, PTID) %>%
  distinct()

cat(paste0(total_exp_n, " total experiments found in manifest, of which ", run_exp_n,
           " can be processed (i.e. non control experiments)\n"))

cat(paste0(nrow(exp_id_vis),
           " total experiments/visit/id combinations found in manifest, of which ",
           nrow(run_exp_id_vis),
           " can be processed (i.e. all available files)\n"))









####### Processing by experiment, visit, and ptid --------------------------------------------------------------------
overall_start_time <- Sys.time()

for (i in 1:nrow(run_exp_id_vis)) {

  cat('\n\n\n')
  timestamp()
  start_time <- Sys.time()


  # Getting experiment details
  current_exp = run_exp_id_vis$EXPERIMENT_NAME[i]
  current_visit = run_exp_id_vis$Visit[i]
  current_ptid = run_exp_id_vis$PTID[i]

  cat(paste0("Processing Experiment: ", current_exp,
             ", PTID: ", current_ptid, ', Visit: ', current_visit,
             '\n(', i,' of ',nrow(run_exp_id_vis), ')\n'))


  current_manifest <- manifest_to_run %>%
    filter(EXPERIMENT_NAME == current_exp &
             # Sometimes xml doesn't list visit since it covers multiple visits
             (Visit == current_visit | FileType == '.xml') &
             # Sometimes xml doesn't list PTID since it covers multiple
             (PTID == current_ptid | FileType == '.xml')
    )

  # Getting experiment details
  current_tissue = unique(current_manifest$Tissue_State %>% na.omit())
  # Using visit info to map to tissue if missing
  if (length(current_tissue) == 0) {
    current_tissue = case_when(current_visit %in% c('V05','V07A','V09') ~ 'Fresh',
                               current_visit %in% c('V02','V06','V07', 'V08','V09','V10')~ 'Frozen')
    cat('Warning: Using visit information (', current_visit,
        ') to extrapolate tissue status: ', current_tissue, '\n',
        sep = '')
  }

  output_dir <- file.path(output_path,
                          current_ptid,
                          current_visit,
                          current_tissue)

  output_prefix <- paste(current_ptid,
                         current_visit,
                         current_tissue,
                         sep = '_')

  if (!override & dir.exists(output_dir) & length(dir(output_dir)) > 0) {
    cat('"overide" = F and output directory not empty, so skipping')
    next()
  }



  # Getting specific file types
  fcs_files <- current_manifest$Path %>% str_subset("fcs$")
  if (length(fcs_files) == 0) {
    stop(current_exp, "has no FCS files!")
  }
  index_files <- fcs_files[grepl("INX",fcs_files)]
  experiment_files <- fcs_files[!grepl("INX",fcs_files)]
  xml_files <- current_manifest$Path %>% str_subset("xml$")
  csv_files <- current_manifest$Path %>% str_subset("csv$")



  # Getting Time of experiments -----------------------------------------------------------------

  # fcs_times <- read.FCSheader(fcs_files, keyword = c('$DATE','$BTIM','$ETIM')) %>%
  #   map_dfr(as.list) %>%
  #   mutate(diva_name = basename(fcs_files)) %>%
  #   rename_all(~str_replace(., '\\$', '')) %>%
  #   mutate(
  #     BDATE_TIME =  dmy_hms(paste0(DATE, ' ', BTIM)),
  #     EDATE_TIME =  dmy_hms(paste0(DATE, ' ', ETIM)),
  #     File_Interval = interval(BDATE_TIME, EDATE_TIME)
  #   )


  # linking the experiment and index names ------------------------------------------------------
  inx_to_link <- tibble(inx_name = basename(index_files)) %>%
    separate(inx_name, into = c("PTID1", "PTID2", "INX","VISIT", "TUBE", "TUBE2"),
             remove = FALSE, extra = 'drop', sep = "_") %>%
    select(inx_name, TUBE, TUBE2, VISIT)
  # left_join(fcs_times, by = c('inx_name' = 'diva_name'))

  exp_to_link <- tibble(exp_name = basename(experiment_files)) %>%
    separate(exp_name, into = c("PTID1", "PTID2","VISIT", "TUBE", "TUBE2"),
             remove = FALSE, extra = 'drop', sep = "_") %>%
    mutate(TUBE2 = TUBE2 %>% str_replace('.fcs', '')) %>%
    select(exp_name, TUBE, TUBE2, VISIT)
  # left_join(fcs_times, by = c('exp_name' = 'diva_name'))

  # Need to make a single correction to fix linking
  if (current_exp == '190705_AC05_DL')
    exp_to_link$TUBE[exp_to_link$exp_name == 'PubID_151_V08_01_002_016.fcs'] = 'XXX'
  if (current_exp == '190716_AB07_DL')
    exp_to_link$TUBE[exp_to_link$exp_name == 'PubID_005_V08_01_001_013.fcs'] = 'XXX'
  if (current_exp == '190717_AC05_DL') {
    exp_to_link$TUBE[exp_to_link$exp_name == 'PubID_001_V02_03_001_010.fcs'] = '04'
    inx_to_link$TUBE[inx_to_link$inx_name == 'PubID_001_INX_V02_03_002_011.fcs'] = '04'
  }

  # Only match on tubes if we know they can match
  if (all(inx_to_link$TUBE %in% exp_to_link$TUBE)) {
    fcs_files_linked <- full_join(exp_to_link, inx_to_link, by = c('TUBE', 'VISIT'))
  } else if (all(inx_to_link$TUBE2 %in% exp_to_link$TUBE2)) {
    #Try different tube placement
    fcs_files_linked <- full_join(exp_to_link, inx_to_link, by = c('TUBE2', 'VISIT'))
  } else {
    fcs_files_linked <- full_join(exp_to_link, inx_to_link, by = 'VISIT')
  }
  if (any(duplicated(na.omit(fcs_files_linked$inx_name)))) {
    stop('FCS Index file name linking failed')
  }


  # Open Diva Workspace -------------------------------------------------------------------------

  # read the workspace
  diva <- open_diva_xml(xml_files)


  diva_fcs_names <- diva_get_sample_groups(diva) %>%
    filter(!specimen %>% str_detect('Controls') & name %in% basename(fcs_files))
  # All experiment and index fcs files listed in xml should exist in folder
  if (length(fcs_files) != nrow(diva_fcs_names))
    cat('Warning: The following fcs files in xml file but not present:\n',
        paste0(basename(fcs_files)[!basename(fcs_files) %in% diva_fcs_names$name], collapse = '\n'),
        '\n', sep = '')



  # Exceptions for certain experiments
  scale_param <- if_else(current_exp == '190128_AC05_DL', "gate", "tube")

  # Listing experiments that can not use global worksheet, and which FCS file should be used
  need_local_gates <- c(
    '181113_AC03_DL' = 'PubID_051_V05_02_002.fcs',
    '191113_AB08_DL' = 'PubID_080_V07A_01_001.fcs',
    '190815_AP01_SLA' = basename(index_files[1]),
    '190828_AP01_AM' = basename(index_files[1])
  )




  if (!current_exp %in% names(need_local_gates)) {
    # Using global worksheet
    imported_workspace_global <-
      try(diva_to_gatingset(diva,
                            name = current_ptid,
                            subset = c(unique(fcs_files_linked$exp_name),
                                       na.omit(fcs_files_linked$inx_name)),
                            worksheet = 'global',
                            scale_level = scale_param,
                            emptyValue = FALSE,
                            truncate_max_range = FALSE)
      )
  } else {
    # Need to apply a local gating hierarchy to the rest of data
    gate_exp_global <-
      try(diva_to_gatingset(diva,
                            name = current_ptid,
                            subset = need_local_gates[current_exp == names(need_local_gates)],
                            worksheet = 'normal',
                            scale = scale_param,
                            emptyValue = FALSE,
                            truncate_max_range = FALSE)
      )
    imported_workspace_global <- gh_apply_to_new_fcs(gate_exp_global[[1]],
                                                     fcs_files,
                                                     emptyValue = FALSE,
                                                     truncate_max_range = FALSE)
  }
  if (!inherits(imported_workspace_global, "try-error")) {
    imported_workspace_global <- manually_adj_gates(
      ws_in = imported_workspace_global,
      IgD_color_in = ifelse(key == 'fhrc', 'V780-A', 'V800-A'),
      visit_in = current_visit)

    experimental_samples_global <- imported_workspace_global[unique(fcs_files_linked$exp_name)]
    if (any(!is.na(fcs_files_linked$inx_name)))
      index_samples_global <- imported_workspace_global[na.omit(unique(fcs_files_linked$inx_name))]

  } else {
    cat("Warning: failed_workspace_import\n")
    next;
  }




  # Need to add repeated to name if a repeated study
  if (all(fcs_files_linked$exp_name %>% str_detect('Repeat')))
    output_prefix <- paste0(output_prefix, '_Repeat')

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # Create plots for each fcs file to show the gating.
  write("Creating gating plots.", stdout())
  plots_global <- list()
  for (s in sampleNames(experimental_samples_global)) {
    write(s, stdout())
    suppressMessages(
      suppressWarnings(
        plots_global[[s]] <- autoplot(experimental_samples_global[[s]], bins = 64)
      )
    )
    pdf(file = file.path(
      output_dir,
      paste(output_prefix,
            stringr::str_split(s,"_")[[1]][3],
            "globalWorksheet.pdf", sep = "_")
    ))
    suppressWarnings(
      print(plots_global[[s]])
    )
    dev.off()
  }


  # Experiment Level Stats ----------------------------------------------------------------------


  exp_time <- purrr::map_df(unique(fcs_files_linked$exp_name), get_exp_time)

  # get the stats
  towrite_global <-
    gs_pop_get_count_fast(experimental_samples_global) %>%
    dplyr::filter(!grepl("INX", name)) %>%
    dplyr::mutate(tosplit = name) %>%
    separate(tosplit, into = c("PTID1","PTID2", "VISIT", "TUBE"), sep = "_", extra = 'drop') %>%
    dplyr::mutate(PTID = paste0(PTID1,'_',PTID2),
                  proportion = Count / ParentCount,
                  proportion = ifelse(is.nan(proportion), 0, proportion)) %>%
    dplyr::select(-PTID1,PTID2) %>%
    full_join(exp_time, by = 'name')


  #' Write a per-file summary
  write("Writing per file summary table.", stdout())
  write.csv(towrite_global,
            file = file.path(
              output_dir,
              paste(output_prefix,
                    "perfileTable.csv",
                    sep = "_")
            ),
            row.names = TRUE)

  # Transform so that we can merge with diva output
  out_summary <- towrite_global %>%
    ungroup %>%
    filter(Population %>% str_detect('CD19\\+CD20\\+/IgD-(/|$)', negate = TRUE)) %>%
    mutate(full_path = Population, full_parent = Parent,
           Population = basename(Population),
           Parent = basename(Parent),
           "Tube Name" = paste(gsub("V","",VISIT), TUBE, sep = "_"))

  # Need to recode for specific experiment (repeated)
  if (current_exp == "190913_AD01_SLA")
    out_summary$`Tube Name` = paste0(out_summary$`Tube Name`, '_Repeat_1')
  # Need to make a single correction to fix linking
  if (current_exp == '190705_AC05_DL')
    out_summary$`Tube Name`[out_summary$name == 'PubID_151_V08_01_002_016.fcs'] = '08_01_002'
  if (current_exp == '190716_AB07_DL')
    out_summary$`Tube Name`[out_summary$name == 'PubID_005_V08_01_001_013.fcs'] = '08_01_001'
  if (current_exp == '190717_AC05_DL')
    out_summary$`Tube Name`[out_summary$name == 'PubID_001_V02_03_001_010.fcs'] = '02_04'


  # Summarized Counts for Exp -------------------------------------------------------------------

  # Summarize files across tubes
  ptid_summary <-
    towrite_global %>% ungroup %>% group_by(PTID, VISIT, Population, Parent) %>%
    select(PTID, VISIT, Population, Parent, Count, ParentCount, proportion) %>%
    summarize(
      Count = sum(Count),
      ParentCount = sum(ParentCount),
      proportion = ifelse(is.na(sum(Count)/sum(ParentCount)),0,sum(Count)/sum(ParentCount)),
      `.groups` = "drop"
    )

  # Write out a CSV for the PTID and visit.
  write("Writing PTID summary table.", stdout())
  write.csv(ptid_summary,
            file = file.path(
              output_dir,
              paste(output_prefix,
                    "PTIDSummary.csv",
                    sep = "_")
            ),
            row.names = TRUE)



  # Index Level Stats ---------------------------------------------------------------------------

  #Skip in cases of no index files
  if (any(!is.na(fcs_files_linked$inx_name))) {

    no_inx_miss_linked <- fcs_files_linked %>% drop_na(inx_name)

    # make a gating set of the index sorted samples
    index_samples_global <- index_samples_global[no_inx_miss_linked$inx_name]


    # get the stats
    towrite_global_by_index <-
      purrr::map2_dfr(no_inx_miss_linked$exp_name,
                      no_inx_miss_linked$inx_name,
                      run_inx_stats_fun) %>%
      dplyr::filter(!grepl("INX", name)) %>%
      dplyr::mutate(tosplit = name) %>%
      separate(tosplit, into = c("PTID1","PTID2", "VISIT", "TUBE"), sep = "_", extra = 'drop') %>%
      dplyr::mutate(PTID = paste0(PTID1,'_',PTID2),
                    proportion = Count / ParentCount,
                    proportion = ifelse(is.nan(proportion), 0, proportion)) %>%
      left_join(current_manifest %>% select(PTID:File, INX_Population:INX_Number_Of_Cells),
                by = c('PTID', 'VISIT' = 'Visit', 'INX_File' = 'File')) %>%
      # Only want INX_Number_Of_Cells for sorting populations
      mutate(
        INX_Number_Of_Cells = case_when(
          INX_Population == 'Antigen-Specific' &
            Population %>% str_detect('eODGT8Double\\+$') ~ INX_Number_Of_Cells,
          INX_Population == 'Bulk' & IgG_Gate == 'IgG+' &
            Population %>% str_detect('IgD\\-IgG\\+/CD20hiCD38hi$|CD27hiCD38hi/IgD\\-IgG\\+$') ~ INX_Number_Of_Cells,
          INX_Population == 'Bulk' & IgG_Gate == 'IgD-' &
            Population %>% str_detect('CD27hiCD38hi/IgD\\-$') ~ INX_Number_Of_Cells
        ),
        Percent_Sorted = round(INX_Number_Of_Cells / Count * 100, 2)
      ) %>%
      dplyr::select(-PTID1,PTID2) %>%
      arrange(Min_Time, INX_Population)


    # Checking the sorting counts vs flow counts (flow should be slightly > num sorted)
    sort_link_check <- towrite_global_by_index %>%
      group_by(name, INX_File, INX_Population) %>%
      summarise(n = sum(!is.na(INX_Number_Of_Cells)),
                `.groups` = "drop")

    if (any(sort_link_check$n != 1))
      cat('Warning: At least one index file could not be linked to proper population\n')


    cat("The following counts and number sorted for each index file (Count, Num Sorted, Percent):\n",
        towrite_global_by_index %>%
          drop_na(INX_Number_Of_Cells) %>%
          unite(info, name, INX_File, INX_Population, Count, INX_Number_Of_Cells, Percent_Sorted, sep = '; ') %>%
          transmute(info = paste0(info, '\n')) %>%
          pull()
    )

    if (any(na.omit(towrite_global_by_index$Percent_Sorted) > 100))
      cat('Warning: At least one case where INX_Number_Of_Cells > Count\n')



    #' Write a per-file summary
    write("Writing per index file summary table.", stdout())
    write.csv(towrite_global_by_index,
              file = file.path(
                output_dir,
                paste(output_prefix,
                      "perINXfileTable.csv",
                      sep = "_")
              ),
              row.names = TRUE,
              na = '')


    # Summarized Counts for Bulk/A.S. -------------------------------------------------------------------


    # Summarize files across bulk and antigen specific types
    type_summary <-
      towrite_global_by_index %>% ungroup %>% group_by(PTID, VISIT, INX_Population, Population, Parent) %>%
      select(PTID, VISIT, INX_Population, Population, Parent, Count, ParentCount, proportion, INX_Number_Of_Cells) %>%
      summarize(
        Count = sum(Count),
        ParentCount = sum(ParentCount),
        proportion = ifelse(is.na(sum(Count)/sum(ParentCount)),0,sum(Count)/sum(ParentCount)),
        INX_Number_Of_Cells = sum(INX_Number_Of_Cells),
        `.groups` = "drop"
      ) %>%
      mutate(Percent_Sorted = round(INX_Number_Of_Cells / Count * 100, 2))


    # Write out a CSV for the PTID and visit.
    write("Writing Type (Bulk/Antigen Specific) summary table.", stdout())
    write.csv(type_summary,
              file = file.path(
                output_dir,
                paste(output_prefix,
                      "TypeSummary.csv",
                      sep = "_")
              ),
              row.names = TRUE,
              na = '')


  } else {
    cat('No index plots found for visit: ', current_visit)
  }




  # Reading in Diva CSV File(s) -----------------------------------------------------------------

  diva_csv_results <-
    map_dfr(csv_files, function(x) {
      xx <- try(read_diva_table(x), silent = TRUE)
      if (inherits(xx, "try-error")) {
        print(paste0("failed_diva_csv_read:", x))
        next;
      } else {
        xx
      }
    })

  # Need to recode if mislabeled visit 7/7A
  if (current_visit == "V07")
    diva_csv_results$`Tube Name` <- gsub("07A_","07_",diva_csv_results$`Tube Name`)
  if (current_visit == "V07A")
    diva_csv_results$`Tube Name` <- gsub("07_","07A_",diva_csv_results$`Tube Name`)
  # Need to add "_01" if missing that part of tube name
  bad_tube_index <- diva_csv_results$`Tube Name` %>% str_detect('_', negate = TRUE)
  if (any(bad_tube_index))
    diva_csv_results$`Tube Name`[bad_tube_index] <-
    paste0(diva_csv_results$`Tube Name`[bad_tube_index], '_01')




  # Comparing Diva (csv) and FCS at Experiment Level ----------------------------------------------------------------------------
  # If there are missing csv files need to add a warning and skip that comparison
  missing_csv_tubes <- bind_rows(
    anti_join(out_summary %>% distinct(`Tube Name`),
              diva_csv_results %>% distinct(`Tube Name`),
              by = 'Tube Name'),
    anti_join(diva_csv_results %>% distinct(`Tube Name`),
              out_summary %>% distinct(`Tube Name`),
              by = 'Tube Name')
  )
  if (nrow(missing_csv_tubes) > 0) {
    cat('Warning: csv file not present for the following Visit/Tube:\n',
        paste0(missing_csv_tubes$`Tube Name`, collapse = '\n'), '\n')
    out_summary_to_compare <- out_summary %>%
      filter(!`Tube Name` %in% missing_csv_tubes)
  } else {
    out_summary_to_compare <- out_summary
  }


  # combining diva and fcs results for comparisons and listing unmatched values
  comparison_results <- full_join(
    out_summary_to_compare %>%
      filter(!`Tube Name` %in% (missing_csv_tubes %>% pull())),
    diva_csv_results %>%
      filter(!`Tube Name` %in% (missing_csv_tubes %>% pull())),
    by = c("PTID" = "Specimen Name", "Tube Name","Population", "Parent")
  ) %>%
    mutate(count_diff = abs(Count - count_diva)) %>%
    select(PTID, VISIT, `Tube Name`, name,
           Population, Parent, Min_Time, Max_Time,
           Count, ParentCount, proportion,
           count_diva, proportion_diva, count_diff)

  if (any(is.na(comparison_results$proportion)))
    cat("The following populations are not available in the fcs files:\n",
        comparison_results %>%
          filter(is.na(proportion)) %>%
          transmute(paste0(Parent,'/',Population)) %>%
          distinct() %>% pull() %>% paste(collapse = ', \n'), '\n'
    )
  if (any(is.na(comparison_results$proportion_diva)))
    cat("The following populations are not available in the diva csv files:\n",
        comparison_results %>%
          filter(is.na(proportion_diva)) %>%
          transmute(paste0(Parent,'/',Population)) %>%
          distinct() %>% pull() %>% paste(collapse = ', \n'), '\n'
    )


  # Getting spearman between diva and fcs
  spearman_label <-  comparison_results %>%
    group_by(`Tube Name`) %>% do(rho = {
      cor(log(.$proportion),
          log(.$proportion_diva),
          method = "spearman",
          use = 'complete.obs')}) %>%
    unnest(rho) %>% mutate(rho_label = paste("rho = ", signif(rho,3)))


  comparison_results_out <- full_join(
    comparison_results,
    spearman_label %>% select(`Tube Name`,  correlation = rho),
    by = 'Tube Name'
  )


  # Write out a CSV comparison with diva
  write("Writing concordance with diva table.", stdout())
  write.csv(comparison_results_out,
            file = file.path(
              output_dir,
              paste(output_prefix,
                    "concordance.csv",
                    sep = "_")
            ),
            row.names = TRUE,
            na = '')



  # Concordance plot with diva -----------------------------------------------------------------

  # Create a concordance plot of the extracted vs the DiVA proportions.
  write("Generating VISC/DIVA concordance plot.", stdout())


  pdf(file = file.path(
    output_dir,
    paste(output_prefix,
          "concordancePlot.pdf",
          sep = "_")
  ))

  concordance_plot <- comparison_results %>%
    filter(!is.na(proportion)) %>%
    ggplot() + geom_point() +
    aes(x = proportion * 100, y = proportion_diva) +
    geom_smooth(method = "lm", formula = y ~ x) +
    ylab("% of Parent Population (Diva)") +
    xlab("% of Parent Population (VISC)") +
    facet_wrap(~`Tube Name`) +
    ggtitle(paste(
      "Agreement between DiVa and VISC for ",
      paste(current_ptid, current_visit, current_tissue, sep = "_"),
      "\n",
      paste("Subject", current_ptid,
            "Visit", current_visit,
            "Tissue", current_tissue)
    )) +
    ggrepel::geom_label_repel(aes(label = Population), force = 5) +
    geom_text(data = spearman_label,
              aes(label = rho_label),
              x = 15,
              y = 90)
  print(concordance_plot)
  dev.off()



  # Clearing temp files -------------------------------------------------------------------------

  gs_cleanup_temp(imported_workspace_global)

  print(Sys.time() - start_time)

}

cat('Overall Run Time:')
print(Sys.time() - overall_start_time)

