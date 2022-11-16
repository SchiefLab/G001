#' Traverses the processed flow data an collates it into a single file.
args = commandArgs(trailingOnly = TRUE);

library(data.table)
library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(stringr)
library(janitor)
library(wCorr)

flow_path <- args[1]

# Loading in manifest -------------------------------------------------------------------------
fhcrc_manifest <- read_csv(file.path(flow_path, "fhrc/fhcrc_manifest.csv"), col_select = -1)
vrc_manifest <- read_csv(file.path(flow_path, "vrc/vrc_manifest.csv"), col_select = -1) %>% 
  mutate(Tube = as.character(Tube))

# make full manifest
full_manifest <- bind_rows(fhcrc_manifest, vrc_manifest)


# Reading and parsing out filenames -------------------------------------------------------------------------


#' Getting info based on filename
#'
#' @param path_prefix_in
#' @param pattern_in
#'
#' @return tibble of path, ptid, visit, tisse, and csv name
#'
get_filenames_fun <-
  function(path_prefix_in = "flow_results",
           pattern_in = "perfileTable\\.csv") {
    tibble(
      filepath = list.files(
        file.path(path_prefix_in),
        pattern = pattern_in,
        recursive = TRUE,
        full.names = TRUE
      ) %>%
        str_subset('Archive', negate = TRUE)
    ) %>%
      mutate(name_to_split = str_replace(filepath, paste0(path_prefix_in, "/"), '')) %>%
      # parse into metadata
      separate(
        name_to_split,
        into = c('Site', "PTID", "VISIT", "Tissue", "csvname"),
        sep = "/",

        remove = TRUE
      ) %>%
      # Dealing with Repeated experiments
      mutate(nrep = if_else(csvname %>% str_detect('Repeat'),
                            2, 1))
  }

perfile_files <- get_filenames_fun(flow_path, "perfileTable\\.csv")
PTID_files <- get_filenames_fun(flow_path, "PTIDSummary\\.csv")
inx_files <- get_filenames_fun(flow_path, "perINXfileTable\\.csv")
type_files <- get_filenames_fun(flow_path, "TypeSummary\\.csv")
QC_files <- get_filenames_fun(flow_path, "concordance\\.csv")


# Fixing Gate Names ---------------------------------------------------------------------------

fix_gates_fun <- function(results_in) {

  results_in %>%
    mutate(
      # CXCRhiDP1hi -> CXCR5hiPD1hi
      Population = Population %>% str_replace('CXCRhiDP1hi', 'CXCR5hiPD1hi'),
      Parent = Parent %>% str_replace('CXCRhiDP1hi', 'CXCR5hiPD1hi'),

      # eODKO11 -> eODKO11-
      Population = Population %>% str_replace('eODKO11(?!-)', 'eODKO11-'),
      Parent = Parent %>% str_replace('eODKO11(?!-)', 'eODKO11-'),

      # typo in live gate name; Singlets/CD14- -> Singlets/Live|CD14\\-
      Population = Population %>% str_replace('Singlets/CD14\\-', 'Singlets/Live|CD14-'),
      Parent = Parent %>% str_replace('Singlets/CD14\\-', 'Singlets/Live|CD14-'),

      # For Plasmablast need to make CD19+CD20lo to match naming convention
      Population = if_else(VISIT == 'V07A',
                           Population %>% str_replace('Singlets/Live\\|Dump\\-/CD19\\+CD20\\+', 'Singlets/Live|Dump-/CD19+CD20lo'),
                           Population),
      Parent = if_else(VISIT == 'V07A',
                       Parent %>% str_replace('Singlets/Live\\|Dump\\-/CD19\\+CD20\\+', 'Singlets/Live|Dump-/CD19+CD20lo'),
                       Parent),

      # VRC FNA samples mislabeled GT8++KO- as GT8++KO+
      Population = Population %>% str_replace('GT8\\+\\+noKO/GT8\\+\\+KO\\+', 'GT8++noKO/GT8++KO-')
    ) %>%
    filter(
      # Keep in the BB515vsBV711 gate for now to send to lab
      # Population %>% str_detect('BB515vsBV711', negate = TRUE)
    )
}

# Reading in and combining flow data files ------------------------------------------------------------------------

read_in_flow_results_fun <-
  function(file_names_in,
           col_types_in,
           merge_in_index_file_info = FALSE) {
    temp_results <-
      full_join(
        file_names_in,
        map_dfr(file_names_in$filepath,
                ~ read_csv(.x, col_types = col_types_in, col_select = -1) %>%
                  # Dealing with Repeated experiments
                  mutate(nrep = if_else(.x %>% str_detect('Repeat'), 2, 1))),
        by = c('PTID', 'VISIT', 'nrep')
      ) %>%
      select(-filepath) %>%
      # Fixing gating names
      fix_gates_fun()


    if (merge_in_index_file_info) {
      temp_results  %>%
        left_join(
          full_manifest %>% select(-Tube),
          by = c(
            'PTID',
            'VISIT' = 'Visit' ,
            'EXPERIMENT_NAME',
            'INX_File' = 'File',
            'INX_Population',
            'KO_Probe',
            'IgG_Gate',
            'Tissue_State',
            'INX_Number_Of_Cells'
          )
        ) %>%
        select(-Path, -Note, -`Additional Notes`, -FileType)

    } else {
      temp_results
    }
  }

# Index level results -------------------------------------------------------------------------

# setting the columns expecting to see in the results csv files
index_level_cols <- cols(
  PTID = col_character(),
  VISIT = col_character(),
  INX_Population = col_character(),
  Population = col_character(),
  Parent = col_character(),
  Count = col_double(),
  ParentCount = col_double(),
  proportion = col_double(),
  INX_Number_Of_Cells = col_double(),
  Percent_Sorted = col_double()
)


# data by index file
flow_data_by_inx_files <- read_in_flow_results_fun(inx_files, index_level_cols, merge_in_index_file_info = TRUE)
#data by type (antigen specific/bulk)
flow_data_by_type <- read_in_flow_results_fun(type_files, index_level_cols)


# File level results -------------------------------------------------------------------------


# setting the columns expecting to see in the results csv files
file_level_cols <- cols(
  PTID = col_character(),
  VISIT = col_character(),
  Population = col_character(),
  Parent = col_character(),
  Count = col_double(),
  ParentCount = col_double(),
  proportion = col_double()
)

# data by fcs file
flow_data_by_file <- read_in_flow_results_fun(perfile_files, file_level_cols)
# data by sample
flow_data_by_PTID <- read_in_flow_results_fun(PTID_files, file_level_cols)


# QC results -------------------------------------------------------------------------

# setting the columns expecting to see in the results csv files
qc_level_cols <- cols(
  PTID = col_character(),
  VISIT = col_character(),
  `Tube Name` = col_character(),
  name = col_character(),
  Population = col_character(),
  Parent = col_character(),
  Min_Time = col_double(),
  Max_Time = col_double(),
  Count = col_double(),
  ParentCount = col_double(),
  proportion = col_double(),
  count_diva = col_double(),
  proportion_diva = col_double(),
  count_diff = col_double(),
  correlation = col_double()
)


# data by index file
qc_combined <- full_join(
  QC_files,
  map_dfr(QC_files$filepath,
          ~ read_csv(.x, col_types = qc_level_cols, col_select = -1) %>%
            # Dealing with Repeated experiments
            mutate(nrep = if_else(.x %>% str_detect('Repeat'), 2, 1)) %>%
            group_by(`Tube Name`) %>%
            mutate(cor_drop_low = cor(proportion[ParentCount > 5] , proportion_diva[ParentCount > 5],
                                      method = 'spearman', use = 'complete.obs'),
                   weighted_cor = weightedCorr(proportion[!is.na(ParentCount)], proportion_diva[!is.na(ParentCount)],
                                               method = "Spearman", weights = log10(ParentCount[!is.na(ParentCount)] + 1)),
                   pearson_cor = cor(proportion , proportion_diva,
                                     method = 'pearson', use = 'complete.obs'),
                   pearson_weighted_cor = weightedCorr(proportion[!is.na(ParentCount)], proportion_diva[!is.na(ParentCount)],
                                                       method = "Pearson", weights = log10(ParentCount[!is.na(ParentCount)] + 1))
            ) %>%
            ungroup()),
  by = c('PTID', 'VISIT', 'nrep')
) %>%
  mutate(prop_diff = abs(proportion * 100 - proportion_diva)) %>%
  select(-contains('cor'), everything(), -filepath)




# Outputting Results --------------------------------------------------------------------------

#string detect needs fhcrc and vrc specisl

dir.create(file.path(flow_path, 'Combined_Results'), recursive = TRUE, showWarnings = FALSE)

write_csv(flow_data_by_inx_files %>% filter(Site == 'fhrc'), file.path(flow_path, 'Combined_Results','FHCRC_Flow_Results_by_Index_File.csv'), na = '')
write_csv(flow_data_by_type %>% filter(Site == 'fhrc'), file.path(flow_path, 'Combined_Results','FHCRC_Flow_Results_by_Bulk_Antigen_Specific.csv'), na = '')

write_csv(flow_data_by_file %>% filter(Site == 'fhrc'), file.path(flow_path, 'Combined_Results','FHCRC_Flow_Results_by_Experiment.csv'), na = '')
write_csv(flow_data_by_PTID %>% filter(Site == 'fhrc'), file.path(flow_path, 'Combined_Results','FHCRC_Flow_Results_by_PTID_Visit.csv'), na = '')

write_csv(qc_combined %>% filter(Site == 'fhrc'), file.path(flow_path, 'Combined_Results','FHCRC_QC_Concordance.csv'), na = '')



write_csv(flow_data_by_inx_files %>% filter(Site == 'vrc'), file.path(flow_path, 'Combined_Results','VRC_Flow_Results_by_Index_File.csv'), na = '')
write_csv(flow_data_by_type %>% filter(Site == 'vrc'), file.path(flow_path, 'Combined_Results','VRC_Flow_Results_by_Bulk_Antigen_Specific.csv'), na = '')

write_csv(flow_data_by_file %>% filter(Site == 'vrc'), file.path(flow_path, 'Combined_Results','VRC_Flow_Results_by_Experiment.csv'), na = '')
write_csv(flow_data_by_PTID %>% filter(Site == 'vrc'), file.path(flow_path, 'Combined_Results','VRC_Flow_Results_by_PTID_Visit.csv'), na = '')

write_csv(qc_combined %>% filter(Site == 'vrc'), file.path(flow_path, 'Combined_Results','VRC_QC_Concordance.csv'), na = '')





# Getting Plate Level Data --------------------------------------------------------------------


flow_data_by_plate <- flow_data_by_inx_files %>%
  filter(!is.na(INX_Plate)) %>%
  group_by(PTID, VISIT, nrep, INX_Population, INX_Plate, Population, Parent) %>%
  select(PTID, VISIT, nrep, INX_Population, INX_Plate, Population, Parent, Count, ParentCount, proportion, INX_Number_Of_Cells) %>%
  summarize(
    Count = sum(Count),
    ParentCount = sum(ParentCount),
    proportion = ifelse(is.na(sum(Count)/sum(ParentCount)),0,sum(Count)/sum(ParentCount)),
    INX_Number_Of_Cells = sum(INX_Number_Of_Cells),
    `.groups` = "drop"
  ) %>%
  mutate(Percent_Sorted = round(INX_Number_Of_Cells / Count * 100, 2))





# Final Flow Data -----------------------------------------------------------------------------


pbmc_pop <- '/Lymphocytes/Singlets'
bcell_pop <- c(
  '/Lymphocytes/Singlets/Live|CD14-CD3-',
  '/Lymphocytes/Singlets/Live|Dump-/CD19+CD20+',
  '/Lymphocytes/Singlets/Live|Dump-/CD19+CD20lo'
)
igd_igg_pop <- c(
  '/Lymphocytes/Singlets/Live|CD14-CD3-/CD19+CD20+/IgD-IgG+',
  '/Lymphocytes/Singlets/Live|Dump-/CD19+CD20lo/CD27hiCD38hi/IgD-IgG+',
  '/Lymphocytes/Singlets/Live|Dump-/CD19+CD20+/IgD-IgG+'
)
igd_pop <- c(
  '/Lymphocytes/Singlets/Live|CD14-CD3-/CD19+CD20+/IgD-',
  '/Lymphocytes/Singlets/Live|Dump-/CD19+CD20lo/CD27hiCD38hi/IgD-',
  '/Lymphocytes/Singlets/Live|Dump-/CD19+CD20+/IgD-'
)

bulk_pop <- c(
  '/Lymphocytes/Singlets/Live|CD14-CD3-/CD19+CD20+/IgD-IgG+/CD20hiCD38hi',
  '/Lymphocytes/Singlets/Live|Dump-/CD19+CD20lo/CD27hiCD38hi/IgD-',
  '/Lymphocytes/Singlets/Live|Dump-/CD19+CD20lo/CD27hiCD38hi/IgD-IgG+'
)
ant_spec_pop <- c(
  '/Lymphocytes/Singlets/Live|CD14-CD3-/CD19+CD20+/IgD-IgG+/CD20hiCD38hi/eODKO11-/eODGT8Double+',
  '/Lymphocytes/Singlets/Live|Dump-/CD19+CD20+/IgD-IgG+/eODKO11-/eODGT8Double+',
  '/Lymphocytes/Singlets/Live|Dump-/CD19+CD20lo/CD27hiCD38hi/IgD-/eODKO11-/eODGT8Double+',
  '/Lymphocytes/Singlets/Live|Dump-/CD19+CD20lo/CD27hiCD38hi/IgD-IgG+/eODKO11-/eODGT8Double+'
)



final_flow_data <- flow_data_by_PTID %>%
  # Dropping first run for visits with multiple rep V02, V06, V07
  group_by(PTID,VISIT) %>%
  filter(nrep == max(nrep)) %>%
  # Getting PBMCs
  group_by(PTID, VISIT, Tissue) %>%
  summarise(
    root_Population = 'root',
    root_count = ParentCount[Parent == 'root'],

    PBMC_Population = pbmc_pop,
    PBMC_Count = Count[Population == pbmc_pop],

    BCell_Population = Population[Population %in% bcell_pop],
    BCell_Count = Count[Population %in% bcell_pop],

    IgD_IgG_Population =  ifelse(any(Population %in% igd_igg_pop),
                                     Population[Population %in% igd_igg_pop],
                                 NA_character_),
    IgD_IgG_Count = ifelse(any(Population %in% igd_igg_pop),
                           Count[Population %in% igd_igg_pop],
                                 NA_real_),

    IgD_Population =  ifelse(any(Population %in% igd_pop),
                                 Population[Population %in% igd_pop],
                                 NA_character_),
    IgD_Count = ifelse(any(Population %in% igd_pop),
                           Count[Population %in% igd_pop],
                           NA_real_),

    # only fresh has a bulk level
    Bulk_Population = ifelse(all(Tissue == 'Fresh'), Population[Population %in% bulk_pop], NA_character_),
    Bulk_Count = ifelse(all(Tissue == 'Fresh'), Count[Population %in% bulk_pop], NA_real_),

    IgG_KOneg_Population =  Population[Population %in% dirname(ant_spec_pop)],
    IgG_KOneg_count =  Count[Population %in% dirname(ant_spec_pop)],

    AG_Specific_Population = Population[Population %in% ant_spec_pop],
    AG_Specific_Count = Count[Population %in% ant_spec_pop],

    # Getting populations Kristen Requested
    IgD_IgG_GT8pospos_noKO_Population =
      ifelse(any(Population %>%
                   str_detect('(?=.*IgD-IgG\\+)(?=.*GT8\\+\\+noKO$)')),
             Population[Population %>%
                          str_detect('(?=.*IgD-IgG\\+)(?=.*GT8\\+\\+noKO$)')],
             NA_character_),
    IgD_IgG_GT8pospos_noKO_Count =
      ifelse(any(Population %>%
                   str_detect('(?=.*IgD-IgG\\+)(?=.*GT8\\+\\+noKO$)')),
             Count[Population %>%
                     str_detect('(?=.*IgD-IgG\\+)(?=.*GT8\\+\\+noKO$)')],
             NA_real_),

    IgD_IgG_GT8pospos_KOneg_Population =
      ifelse(any(Population %>% str_detect('(?=.*IgD-IgG\\+)(?=.*GT8\\+\\+KO\\-$)')),
             Population[Population %>% str_detect('(?=.*IgD-IgG\\+)(?=.*GT8\\+\\+KO\\-$)')],
             NA_character_),
    IgD_IgG_GT8pospos_KOneg_Count =
      ifelse(any(Population %>% str_detect('(?=.*IgD-IgG\\+)(?=.*GT8\\+\\+KO\\-$)')),
             Count[Population %>% str_detect('(?=.*IgD-IgG\\+)(?=.*GT8\\+\\+KO\\-$)')],
             NA_real_),


    IgD_GT8pospos_noKO_Population =
      ifelse(any(Population %>%
                   str_detect('(?=.*IgD-/)(?=.*GT8\\+\\+noKO$)')),
             Population[Population %>%
                          str_detect('(?=.*IgD-/)(?=.*GT8\\+\\+noKO$)')],
             NA_character_),
    IgD_GT8pospos_noKO_Count =
      ifelse(any(Population %>%
                   str_detect('(?=.*IgD-/)(?=.*GT8\\+\\+noKO$)')),
             Count[Population %>%
                     str_detect('(?=.*IgD-/)(?=.*GT8\\+\\+noKO$)')],
             NA_real_),

    IgD_GT8pospos_KOneg_Population =
      ifelse(any(Population %>% str_detect('(?=.*IgD-/)(?=.*GT8\\+\\+KO\\-$)')),
             Population[Population %>% str_detect('(?=.*IgD-/)(?=.*GT8\\+\\+KO\\-$)')],
             NA_character_),
    IgD_GT8pospos_KOneg_Count =
      ifelse(any(Population %>% str_detect('(?=.*IgD-/)(?=.*GT8\\+\\+KO\\-$)')),
             Count[Population %>% str_detect('(?=.*IgD-/)(?=.*GT8\\+\\+KO\\-$)')],
             NA_real_),

    `.groups` = 'drop'
  )


# Sample Swap
if (length(args) > 1) {
  # Applying Sample Swap fix: 2 PTIDs were swapped for V08 and V10
  # From David Leggat email 11/06/2020
  swap_info <- read_csv(args[2])
  final_flow_data <- final_flow_data %>%
    full_join(swap_info,
              by = c('PTID' = 'Bad_PTID', 'VISIT')) %>%
    mutate(PTID = if_else(is.na(Correct_PTID), PTID, Correct_PTID)) %>%
    select(-Correct_PTID)
}


write_csv(final_flow_data,
          file.path(flow_path, 'Combined_Results',
                    'Wide_Flow_Data_to_Merge.csv'),
          na = '')





# Final Flow Data by Bulk/Ant Specific

final_flow_data_by_type <- flow_data_by_type %>%
  # Dropping first run for pt  V02, V06, V07
  filter(!(PTID == 'PubID_028' & VISIT %in% c('V02', 'V06', 'V07') & nrep == 1)) %>%
  # Getting PBMCs
  group_by(PTID, VISIT, Tissue, INX_Population) %>%
  summarise(
    root_Population = 'root',
    root_count = ParentCount[Parent == 'root'],

    PBMC_Population = pbmc_pop,
    PBMC_Count = Count[Population == pbmc_pop],

    BCell_Population = Population[Population %in% bcell_pop],
    BCell_Count = Count[Population %in% bcell_pop],

    IgD_IgG_Population =  ifelse(any(Population %in% igd_igg_pop),
                                 Population[Population %in% igd_igg_pop],
                                 NA_character_),
    IgD_IgG_Count = ifelse(any(Population %in% igd_igg_pop),
                           Count[Population %in% igd_igg_pop],
                           NA_real_),

    IgD_Population =  ifelse(any(Population %in% igd_pop),
                             Population[Population %in% igd_pop],
                             NA_character_),
    IgD_Count = ifelse(any(Population %in% igd_pop),
                       Count[Population %in% igd_pop],
                       NA_real_),

    # only fresh has a bulk level
    Bulk_Population = ifelse(all(Tissue == 'Fresh'), Population[Population %in% bulk_pop], NA_character_),
    Bulk_Count = ifelse(all(Tissue == 'Fresh'), Count[Population %in% bulk_pop], NA_real_),

    IgG_KOneg_Population =  Population[Population %in% dirname(ant_spec_pop)],
    IgG_KOneg_count =  Count[Population %in% dirname(ant_spec_pop)],

    AG_Specific_Population = Population[Population %in% ant_spec_pop],
    AG_Specific_Count = Count[Population %in% ant_spec_pop],

    # Getting populations Kristen Requested
    IgD_IgG_GT8pospos_noKO_Population =
      ifelse(any(Population %>%
                   str_detect('(?=.*IgD-IgG\\+)(?=.*GT8\\+\\+noKO$)')),
             Population[Population %>%
                          str_detect('(?=.*IgD-IgG\\+)(?=.*GT8\\+\\+noKO$)')],
             NA_character_),
    IgD_IgG_GT8pospos_noKO_Count =
      ifelse(any(Population %>%
                   str_detect('(?=.*IgD-IgG\\+)(?=.*GT8\\+\\+noKO$)')),
             Count[Population %>%
                     str_detect('(?=.*IgD-IgG\\+)(?=.*GT8\\+\\+noKO$)')],
             NA_real_),

    IgD_IgG_GT8pospos_KOneg_Population =
      ifelse(any(Population %>% str_detect('(?=.*IgD-IgG\\+)(?=.*GT8\\+\\+KO\\-$)')),
             Population[Population %>% str_detect('(?=.*IgD-IgG\\+)(?=.*GT8\\+\\+KO\\-$)')],
             NA_character_),
    IgD_IgG_GT8pospos_KOneg_Count =
      ifelse(any(Population %>% str_detect('(?=.*IgD-IgG\\+)(?=.*GT8\\+\\+KO\\-$)')),
             Count[Population %>% str_detect('(?=.*IgD-IgG\\+)(?=.*GT8\\+\\+KO\\-$)')],
             NA_real_),


    IgD_GT8pospos_noKO_Population =
      ifelse(any(Population %>%
                   str_detect('(?=.*IgD-/)(?=.*GT8\\+\\+noKO$)')),
             Population[Population %>%
                          str_detect('(?=.*IgD-/)(?=.*GT8\\+\\+noKO$)')],
             NA_character_),
    IgD_GT8pospos_noKO_Count =
      ifelse(any(Population %>%
                   str_detect('(?=.*IgD-/)(?=.*GT8\\+\\+noKO$)')),
             Count[Population %>%
                     str_detect('(?=.*IgD-/)(?=.*GT8\\+\\+noKO$)')],
             NA_real_),

    IgD_GT8pospos_KOneg_Population =
      ifelse(any(Population %>% str_detect('(?=.*IgD-/)(?=.*GT8\\+\\+KO\\-$)')),
             Population[Population %>% str_detect('(?=.*IgD-/)(?=.*GT8\\+\\+KO\\-$)')],
             NA_character_),
    IgD_GT8pospos_KOneg_Count =
      ifelse(any(Population %>% str_detect('(?=.*IgD-/)(?=.*GT8\\+\\+KO\\-$)')),
             Count[Population %>% str_detect('(?=.*IgD-/)(?=.*GT8\\+\\+KO\\-$)')],
             NA_real_),

    `.groups` = 'drop'
  )


# Sample Swap
if (length(args) > 1) {
  # Applying Sample Swap fix: 2 PTIDs were swapped for V08 and V10
  # From David Leggat email 11/06/2020
  final_flow_data_by_type <- final_flow_data_by_type %>%
    full_join(swap_info,
              by = c('PTID' = 'Bad_PTID', 'VISIT')) %>%
    mutate(PTID = if_else(is.na(Correct_PTID), PTID, Correct_PTID)) %>%
    select(-Correct_PTID)
}


write_csv(final_flow_data_by_type,
          file.path(flow_path, 'Combined_Results',
                    'Wide_Flow_Data_by_Type_to_Merge.csv'),
          na = '')



