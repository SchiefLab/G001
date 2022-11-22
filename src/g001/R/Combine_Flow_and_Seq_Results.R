#' Combine processed flow and sequence data
args = commandArgs(trailingOnly = TRUE);

library(MIMOSA, warn.conflicts = FALSE)
library(data.table, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
library(purrr, warn.conflicts = FALSE)
library(readr, warn.conflicts = FALSE)
library(stringr, warn.conflicts = FALSE)
library(forcats, warn.conflicts = FALSE)
library(readxl, warn.conflicts = FALSE)
library(arrow, warn.conflicts = FALSE)
library(here, warn.conflicts = FALSE)
library(dplyr, warn.conflicts = FALSE)

fhcrc_manifest_path <- args[1]
vrc_manifest_path <- args[2]
seq_path <- args[3]
collated_output_path <- args[4]
output_path <- args[5]


# Sequence Data Loading and Setup -------------------------------------------------------------


# Getting
paired_candidates <- read_feather(file.path(seq_path, 'pairing_candidates.feather'))
#paired_candidates <- readr::read_csv(paste0(seq_path, 'pairing_candidates.csv'))


unpaired_data_org <- paired_candidates %>%
  filter(dropped_for == 'unpaired')

unpaired_data <- unpaired_data_org %>%
  mutate(
    Tissue = tissue %>% str_to_title(),
    Plate_Type = plate_type %>% str_to_title(),
    pos_call = case_when(
      chain == 'HEAVY' ~ v_call %>% str_detect('IGHV1-2'),
      chain %in% c('KAPPA', 'LAMBDA') ~ nchar(cdr3_aa) == 5
    )
  ) %>%
  # dropping missing calls for now
  filter(!is.na(pos_call)) %>%
  #if multiple chains Heavy must be negative or both Light to be considered a single negative call
  group_by(Subject = pub_id, Dose_Group = dose_group, Timepoint = timepoint, 
           Plate_Type, Tissue, Plate = plate, Well = well) %>%
  summarise(pos_call = min(pos_call[chain == 'Heavy'],
                           max(pos_call[chain != 'Heavy'])
  ),
  `.groups` = "drop") %>%
  # can only include negative values in unpaired data
  filter(!pos_call)




samples_in_paired_long <- paired_candidates %>%
  mutate(Plate_Type = paste0(
    str_replace(plate_type %>% str_to_title(),
                'Probe_specific', 'Ag_Spec'),
    '_Seq_Manifest')) %>%
  distinct(Subject = pub_id, Timepoint = timepoint, Plate_Type) %>%
  mutate(in_paired = TRUE)


samples_in_paired <- samples_in_paired_long %>%
  pivot_wider(names_from = Plate_Type,
              values_from = in_paired)

#seq_data_org <- readr::read_csv(file.path(seq_path, 'personalized.csv'))
seq_data_org <- read_feather(file.path(seq_path, 'personalize.feather'))


seq_data <- seq_data_org %>%
  mutate(
    # Tissue_heavy and Tissue_light should always match
    Tissue = tissue %>% str_to_title(),
    # Plate_Type_heavy and Plate_Type_light should always match
    Plate_Type = plate_type %>% str_to_title()
  ) %>%
  filter(!is.na(is_vrc01_class))


seq_results_long <- bind_rows(
  seq_data %>%
    mutate(is_paired = 'paired'),
) %>%
  group_by(Subject = pub_id, Dose_Group = dose_group, Timepoint = timepoint, 
           Plate_Type, Tissue) %>%
  dplyr::summarise(
    sequenced = dplyr::n(),
    n_VRC01 = sum(is_vrc01_class),

    n_hc_vh1_2_paired = sum(`has_vh1-2`, na.rm = TRUE),
    n_lc_5cdrl3_paired = sum(has_5_len_lcdr3, na.rm = TRUE),

    `.groups` = 'drop'
  ) %>%
  ungroup() %>%
  mutate(Plate_Type = str_replace(Plate_Type, 'Probe_specific', 'Ag_Spec'))



seq_results <- seq_results_long %>%
  pivot_wider(names_from = Plate_Type,
              values_from = c(sequenced, n_VRC01,
                              n_hc_vh1_2_paired, n_lc_5cdrl3_paired
              ))



# Flow Data ----------------------------------------------------------------------


flow_data <- read_csv(file.path(
  collated_output_path,
  'Wide_Flow_Data_to_Merge.csv'), show_col_types = FALSE) %>%
  rename(Visit = VISIT)


flow_data_by_type <- read_csv(file.path(
  collated_output_path,
  'Wide_Flow_Data_by_Type_to_Merge.csv'), show_col_types = FALSE) %>%
  rename(Visit = VISIT)


# Manifests -----------------------------------------------------------------------------------

# Need all manifest index plate counts for QC
fhcrc_flow_manifest <- read_csv(file.path(fhcrc_manifest_path), col_select = -1, show_col_types = FALSE)
vrc_flow_manifest <- read_csv(file.path(vrc_manifest_path), col_select = -1, show_col_types = FALSE) %>% 
  mutate(Tube = as.character(Tube))

flow_manifest_long <- bind_rows(fhcrc_flow_manifest, vrc_flow_manifest) %>%
  filter(!is.na(INX_Number_Of_Cells)) %>%
  group_by(PTID, Visit, Tissue_State, INX_Population) %>%
  summarise(
    inx_count_flow = sum(INX_Number_Of_Cells),
    `.groups` = 'drop'
  ) %>%
  ungroup() %>%
  rename(Tissue = Tissue_State, plate_type = INX_Population) %>%
  mutate(plate_type = str_replace(plate_type, 'Antigen-Specific', 'Ag_Spec'))


flow_manifest <- flow_manifest_long %>%
  pivot_wider(names_from = plate_type,
              values_from = c(inx_count_flow),
              names_prefix = 'Flow_Manifest_')



# Treatment and PubIDs ----------------------------------------------------

#Getting treatment
treat_path <- 'data/groups/G001treatment.csv'

treat_info <- read_csv(treat_path, col_select = -1, show_col_types = FALSE) %>% 
  rename(Site = site)


# Merging All ---------------------------------------------------------------------------------


combined_results <-
  full_join(flow_data,
            seq_results,
            by = c('PTID' = 'Subject', 'Visit' = 'Timepoint', 'Tissue')) %>%
  left_join(flow_manifest,
            by = c('PTID', 'Visit', 'Tissue')) %>%
  left_join(samples_in_paired,
            by = c('PTID' = 'Subject', 'Visit' = 'Timepoint')) %>%
  full_join(treat_info, by = c('PTID' = 'Volunteer Study ID')) %>%
  mutate(
    Dose = Group %>% fct_recode(`Low Dose` = 'Group 1',
                                `High Dose` = 'Group 2'),
    Visit = Visit %>% factor(),

    across(ends_with('_Ag_Spec'),
           ~ifelse(AG_Specific_Count == 0, 0, .)),
    across(ends_with('_Bulk'),
           ~ifelse(Bulk_Count == 0, 0, .)),

    # making seq indicator variable
    Seq_Attempted = case_when(
      AG_Specific_Count == 0  ~ 'No Seq Expected',
      sequenced_Ag_Spec > 0  ~ 'VRC01 Calls Made',
      is.na(Ag_Spec_Seq_Manifest)  ~ 'No Seq Performed',
      !is.na(Ag_Spec_Seq_Manifest)  ~ 'Seq Data but no VRC01 Calls'
    ),
    Bulk_Seq_Attempted = case_when(
      Bulk_Count == 0  ~ 'No Seq Expected (Bulk)',
      sequenced_Bulk > 0  ~ 'VRC01 Calls Made (Bulk)',
      is.na(Bulk_Seq_Manifest)  ~ 'No Seq Performed (Bulk)',
      !is.na(Bulk_Seq_Manifest)  ~ 'Seq Data but no VRC01 Calls (Bulk)'
    ),
    n_not_VRC01_Ag_Spec = sequenced_Ag_Spec - n_VRC01_Ag_Spec,
  ) %>%
  # Bringing in BL values (V02)
  group_by(PTID, Site, Dose) %>%
  mutate(
    AG_Specific_Count_BL = AG_Specific_Count[Visit == 'V02'],
    PBMC_Count_BL = PBMC_Count[Visit == 'V02'],
    IgD_IgG_Count_BL = IgD_IgG_Count[Visit == 'V02'],
    IgD_Count_BL = IgD_Count[Visit == 'V02'],
    n_VRC01_Ag_Spec_BL = n_VRC01_Ag_Spec[Visit == 'V02'],
    n_not_VRC01_Ag_Spec_BL = n_not_VRC01_Ag_Spec[Visit == 'V02'],
    sequenced_Ag_Spec_BL = sequenced_Ag_Spec[Visit == 'V02']
  ) %>%
  ungroup() %>%

  select(PTID, Site, Dose, Treatment, everything(),
         -Dose_Group, -Group) %>%
  arrange(desc(Dose))


write_csv(combined_results,
          file.path('data', 'combined_results', 'Wide_Flow_Seq_Data.csv'), na = '')




# Creating Results Data -----------------------------------------------------------------------

final_results <- combined_results %>%
  mutate(
    #using IgD for PB
    IgD_or_IgD_IgG_Count = ifelse(Visit %in% c('V07A') & !is.na(IgD_Count),
                                  IgD_Count,IgD_IgG_Count),

    BCell_of_PBMC = 100 * BCell_Count / PBMC_Count,

    IgD_IgG_of_PBMC = 100 * IgD_or_IgD_IgG_Count / PBMC_Count,

    IgD_IgG_of_BCell = 100 * IgD_or_IgD_IgG_Count / BCell_Count,

    KOneg_of_IgG = 100 * IgG_KOneg_count / IgD_or_IgD_IgG_Count,

    IgD_IgG_GC_Count = ifelse(Visit %in% c('V05','V09'),
                              Bulk_Count,
                              NA),

    KOneg_of_GC = 100 * IgG_KOneg_count / IgD_IgG_GC_Count,

    Ag_spec_of_KOneg = 100 * AG_Specific_Count / IgG_KOneg_count,
    Ag_spec_of_IgG = 100 * AG_Specific_Count / IgD_or_IgD_IgG_Count,
    Ag_spec_of_IgG_BL = 100 * AG_Specific_Count_BL /
      ifelse(Visit %in% c('V07A'), IgD_Count_BL, IgD_IgG_Count_BL),
    # FNA GC epitope specific
    Ag_spec_of_GC = 100 * AG_Specific_Count / IgD_IgG_GC_Count,
    Ag_spec_of_BCell = 100 * AG_Specific_Count / BCell_Count,
    Ag_spec_of_PBMC = 100 * AG_Specific_Count / PBMC_Count,

    Ag_spec_seq_of_Ag_spec = 100 * sequenced_Ag_Spec / AG_Specific_Count,
    VRC01_of_Ag_spec_seq = 100 * n_VRC01_Ag_Spec / sequenced_Ag_Spec,
    VRC01_of_Ag_spec_seq_bl_bound =
      ifelse(Visit == 'V02' & is.na(n_VRC01_Ag_Spec),
             100,
             VRC01_of_Ag_spec_seq),

    VRC01_of_Ag_spec_seq_BL = 100 * n_VRC01_Ag_Spec_BL / sequenced_Ag_Spec_BL,
    not_VRC01_of_Ag_spec_seq_BL = 100 * n_not_VRC01_Ag_Spec_BL / sequenced_Ag_Spec_BL,

    n_not_VRC01_VH1_2 = n_hc_vh1_2_paired_Ag_Spec - n_VRC01_Ag_Spec,
    n_not_VRC01_CDRL3_5amino = n_lc_5cdrl3_paired_Ag_Spec - n_VRC01_Ag_Spec,
    not_VRC01_of_Ag_spec_seq = 100 * n_not_VRC01_Ag_Spec / sequenced_Ag_Spec,
    not_VRC01_of_Ag_spec_seq_bl_bound =
      ifelse(Visit == 'V02' & is.na(n_not_VRC01_Ag_Spec),
             100,
             not_VRC01_of_Ag_spec_seq),
    not_VRC01_VH1_2_of_seq = 100 * n_not_VRC01_VH1_2 / sequenced_Ag_Spec,
    not_VRC01_CDRL3_5amino_of_seq = 100 * n_not_VRC01_CDRL3_5amino / sequenced_Ag_Spec,


    VRC01_of_BCell =
      if_else(AG_Specific_Count == 0,  0,
              100 * (Ag_spec_of_BCell / 100) * (VRC01_of_Ag_spec_seq / 100)),
    # Need to remove cases where No Seq Performed
    VRC01_of_BCell = if_else(Seq_Attempted == 'No Seq Performed',
                             NA_real_,
                             VRC01_of_BCell),

    not_VRC01_of_BCell =
      if_else(AG_Specific_Count == 0,  0,
              100 * (Ag_spec_of_BCell / 100) * (not_VRC01_of_Ag_spec_seq / 100)),
    # Need to remove cases where No Seq Performed
    not_VRC01_of_BCell = if_else(Seq_Attempted == 'No Seq Performed',
                                 NA_real_,
                                 not_VRC01_of_BCell),


    VRC01_of_IgG_seq_bl_bound =
      if_else(AG_Specific_Count == 0, 0,
              100 * (Ag_spec_of_IgG / 100) * (VRC01_of_Ag_spec_seq_bl_bound / 100),
      ),
    # Need to remove cases where No Seq Performed
    VRC01_of_IgG_seq_bl_bound = if_else(Seq_Attempted == 'No Seq Performed',
                                        NA_real_,
                                        VRC01_of_IgG_seq_bl_bound),
    VRC01_of_IgG =
      if_else(is.na(VRC01_of_Ag_spec_seq) & AG_Specific_Count > 0, NA_real_,
              VRC01_of_IgG_seq_bl_bound),

    VRC01_of_GC =
      if_else(AG_Specific_Count == 0, 0,
              100 * (Ag_spec_of_GC / 100) * (VRC01_of_Ag_spec_seq / 100)
      ),

    not_VRC01_of_IgG_seq_bl_bound =
      if_else(AG_Specific_Count == 0, 0,
              100 * (Ag_spec_of_IgG / 100) * (not_VRC01_of_Ag_spec_seq_bl_bound / 100),
      ),
    not_VRC01_of_IgG_seq_bl_bound = if_else(Seq_Attempted == 'No Seq Performed',
                                            NA_real_,
                                            not_VRC01_of_IgG_seq_bl_bound),
    not_VRC01_of_IgG =
      if_else(is.na(not_VRC01_of_Ag_spec_seq) & AG_Specific_Count > 0, NA_real_,
              not_VRC01_of_IgG_seq_bl_bound),
    not_VRC01_of_GC =
      if_else(AG_Specific_Count == 0, 0,
              100 * (Ag_spec_of_GC / 100) * (not_VRC01_of_Ag_spec_seq / 100)
      ),

    # sup vars
    sequenced_Ag_Spec_sup =  ifelse(
      Seq_Attempted %in% c("No Seq Performed", "Seq Data but no VRC01 Calls"),
      0,sequenced_Ag_Spec
    ),
    VRC01_of_BCell_sup =
      ifelse(
        Seq_Attempted %in% c("No Seq Performed", "Seq Data but no VRC01 Calls"),
        0, VRC01_of_BCell
      ),
    VRC01_of_IgG_sup =
      ifelse(
        Seq_Attempted %in% c("No Seq Performed", "Seq Data but no VRC01 Calls"),
        0, VRC01_of_IgG
      ),
    VRC01_of_GC_sup =
      ifelse(
        Seq_Attempted %in% c("No Seq Performed", "Seq Data but no VRC01 Calls"),
        0, VRC01_of_GC
      ),
    not_VRC01_of_BCell_sup =
      ifelse(
        Seq_Attempted %in% c("No Seq Performed", "Seq Data but no VRC01 Calls"),
        0, not_VRC01_of_BCell
      ),
    not_VRC01_of_IgG_sup =
      ifelse(
        Seq_Attempted %in% c("No Seq Performed", "Seq Data but no VRC01 Calls"),
        0, not_VRC01_of_IgG
      ),
    not_VRC01_of_GC_sup =
      ifelse(
        Seq_Attempted %in% c("No Seq Performed", "Seq Data but no VRC01 Calls"),
        0, not_VRC01_of_GC
      ),




    # Getting new gate flow estimates
    GT8pospos_noKO_of_BCell = 100 * ifelse(Visit %in% c('V07A') & !is.na(IgD_GT8pospos_noKO_Count),
                                           IgD_GT8pospos_noKO_Count / BCell_Count,
                                           IgD_IgG_GT8pospos_noKO_Count / BCell_Count),
    GT8pospos_noKO_of_IgG =  100 * ifelse(Visit %in% c('V07A') & !is.na(IgD_GT8pospos_noKO_Count),
                                          IgD_GT8pospos_noKO_Count / IgD_Count,
                                          IgD_IgG_GT8pospos_noKO_Count / IgD_IgG_Count),
    GT8pospos_noKO_of_GC = ifelse(Visit %in% c('V05','V09'),
                                  100 * IgD_IgG_GT8pospos_noKO_Count / Bulk_Count,
                                  NA),
    GT8pospos_KOneg_of_GT8pospos_noKO = 100 * ifelse(Visit %in% c('V07A') & !is.na(IgD_GT8pospos_KOneg_Count),
                                                     IgD_GT8pospos_KOneg_Count / IgD_GT8pospos_noKO_Count,
                                                     IgD_IgG_GT8pospos_KOneg_Count / IgD_IgG_GT8pospos_noKO_Count),
    GT8pospos_KOneg_of_IgG = 100 * ifelse(Visit %in% c('V07A') & !is.na(IgD_GT8pospos_KOneg_Count),
                                          IgD_GT8pospos_KOneg_Count / IgD_Count,
                                          IgD_IgG_GT8pospos_KOneg_Count / IgD_IgG_Count),

    IgG_of_IgD = 100 * IgD_IgG_Count / IgD_Count,
    IgG_of_IgD_GT8pospos_noKO = 100 * IgD_IgG_GT8pospos_noKO_Count / IgD_GT8pospos_noKO_Count,
    IgG_of_IgD_GT8pospos_KOneg = 100 * IgD_IgG_GT8pospos_KOneg_Count / IgD_GT8pospos_KOneg_Count,


    VRC01_of_GT8pospos =
      if_else(AG_Specific_Count == 0,  0,
              100 * pmin(AG_Specific_Count / ifelse(Visit %in% c('V07A'),
                                                    IgD_GT8pospos_noKO_Count,
                                                    IgD_IgG_GT8pospos_noKO_Count),
                         1) * (VRC01_of_Ag_spec_seq_bl_bound / 100)),
    # Need to remove cases where No Seq Performed
    VRC01_of_GT8pospos = if_else(Seq_Attempted == 'No Seq Performed',
                                 NA_real_,
                                 VRC01_of_GT8pospos),
    VRC01_of_GT8pospos_sup =
      ifelse(
        Seq_Attempted %in% c("No Seq Performed", "Seq Data but no VRC01 Calls"),
        0, VRC01_of_GT8pospos
      ),


    # Getting BL Value for response calls
    VRC01_of_IgG_BL = 100 * (Ag_spec_of_IgG_BL / 100) *
      if_else(is.na(sequenced_Ag_Spec_BL) | sequenced_Ag_Spec_BL == 0,
              1, VRC01_of_Ag_spec_seq_BL / 100),
    not_VRC01_of_IgG_BL = 100 * (Ag_spec_of_IgG_BL / 100) *
      if_else(is.na(sequenced_Ag_Spec_BL) | sequenced_Ag_Spec_BL == 0,
              1, VRC01_of_Ag_spec_seq_BL / 100),

    # response def
    Response = case_when(
      Tissue == 'Fresh' ~ as.numeric(n_VRC01_Ag_Spec > 0),
      n_VRC01_Ag_Spec > 0 & VRC01_of_IgG > VRC01_of_IgG_BL ~ 1,
      !is.na(VRC01_of_IgG) & VRC01_of_IgG <= VRC01_of_IgG_BL ~ 0,
      TRUE ~ NA_real_
    ),
    Response_sup = ifelse(
      Seq_Attempted %in% c("No Seq Performed", "Seq Data but no VRC01 Calls"),
      0, Response
    ),

    Response_notVRC01 = case_when(
      Tissue == 'Fresh' ~ as.numeric(n_not_VRC01_Ag_Spec > 0),
      n_not_VRC01_Ag_Spec > 0 & not_VRC01_of_IgG > not_VRC01_of_IgG_BL ~ 1,
      !is.na(not_VRC01_of_IgG) & not_VRC01_of_IgG <= not_VRC01_of_IgG_BL ~ 0,
      TRUE ~ NA_real_
    ),
    Response_sup_notVRC01 = ifelse(
      Seq_Attempted %in% c("No Seq Performed", "Seq Data but no VRC01 Calls"),
      0, Response_notVRC01
    ),



  ) %>%
  select(
    PTID,
    Site,
    Treatment,
    Visit,
    Tissue,

    `Number of PBMCs sorted` = root_count,
    `Number of PBMCs sorted (Lymphocytes/Singlets)` = PBMC_Count,
    `Number of B cells` = BCell_Count,
    `Percent of PBMCs (Lymphocytes/Singlets) that are B cells` = BCell_of_PBMC,
    `Number of IgD- B cells` = IgD_Count,
    `Number of IgD-IgG+ B cells` = IgD_IgG_Count,
    `Percent of IgD- B cells that are IgG+` = IgG_of_IgD,
    `Percent of B cells that are IgG+ B cells` = IgD_IgG_of_BCell,
    `Number of IgG+ GC B cells` = IgD_IgG_GC_Count,
    `Number of IgG+ B cells that are KO-` = IgG_KOneg_count,
    `Percent of IgG+ B cells that are KO-` = KOneg_of_IgG,
    `Number of epitope-specific (KO-GT8++) IgG+ B cells` = AG_Specific_Count,
    `Percent of IgG+ KO- B cells that are GT8++` = Ag_spec_of_KOneg,
    `Percent of B cells that are epitope-specific (KO-GT8++)` = Ag_spec_of_BCell,
    `Percent of IgG+ B cells that are epitope-specific (KO-GT8++)` = Ag_spec_of_IgG,
    `Percent of IgG+ GC B cells that are epitope-specific (KO-GT8++)` = Ag_spec_of_GC,

    `Number of epitope-specific (KO-GT8++) IgG+ B cells that have BCR heavy and light chains sequenced` = sequenced_Ag_Spec,
    `Number of epitope-specific (KO-GT8++) IgG+ B cells that have BCR heavy and light chains sequenced (missing seq to 0)` = sequenced_Ag_Spec_sup,
    `Percent of epitope-specific (KO-GT8++) IgG+ B cells that have BCR heavy and light chains sequenced` = Ag_spec_seq_of_Ag_spec,

    `Number of epitope-specific (KO-GT8++) sequenced IgG BCRs that are VRC01-class` = n_VRC01_Ag_Spec,
    `Percent of epitope-specific (KO-GT8++) sequenced IgG BCRs that are VRC01-class` = VRC01_of_Ag_spec_seq,

    `Number of epitope-specific (KO-GT8++) sequenced IgG BCRs that are not VRC01-class` = n_not_VRC01_Ag_Spec,

    `Number epitope-specific (KO-GT8++) sequenced IgG BCRs that are not VRC01-class but have a VH1-2 heavy chain` = n_not_VRC01_VH1_2,
    `Number epitope-specific (KO-GT8++) sequenced IgG BCRs that are not VRC01-class but have a 5-aa CDRL3` = n_not_VRC01_CDRL3_5amino,
    `Percent of epitope-specific (KO-GT8++) sequenced IgG BCRs that are not VRC01-class` = not_VRC01_of_Ag_spec_seq,
    `Percent of epitope-specific (KO-GT8++) sequenced IgG BCRs that are not VRC01-class but have a VH1-2 heavy chain` = not_VRC01_VH1_2_of_seq,
    `Percent of epitope-specific (KO-GT8++) sequenced IgG BCRs that are not VRC01-class but have a 5-aa CDRL3` = not_VRC01_CDRL3_5amino_of_seq,

    `Percent of B cells detected as VRC01-class` = VRC01_of_BCell,
    `Percent of B cells detected as VRC01-class (missing seq to 0)` = VRC01_of_BCell_sup,
    `Percent of IgG+ B cells detected as VRC01-class` = VRC01_of_IgG,
    `Percent of IgG+ B cells detected as VRC01-class (Baseline Upper Bound)` = VRC01_of_IgG_seq_bl_bound,
    `Percent of IgG+ B cells detected as VRC01-class (missing seq to 0)` = VRC01_of_IgG_sup,

    `Percent of IgG+ GC B cells detected as VRC01-class` = VRC01_of_GC,
    `Percent of IgG+ GC B cells detected as VRC01-class (missing seq to 0)` = VRC01_of_GC_sup,

    `Percent of B cells that are not VRC01-class` = not_VRC01_of_BCell,
    `Percent of B cells that are not VRC01-class (missing seq to 0)` = not_VRC01_of_BCell_sup,
    `Percent of IgG+ B cells that are not VRC01-class` = not_VRC01_of_IgG,
    `Percent of IgG+ B cells that are not VRC01-class (Baseline Upper Bound)` = not_VRC01_of_IgG_seq_bl_bound,
    `Percent of IgG+ B cells that are not VRC01-class (missing seq to 0)` = not_VRC01_of_IgG_sup,

    `Percent of IgG+ GC B cells that are not VRC01-class` = not_VRC01_of_GC,
    `Percent of IgG+ GC B cells that are not VRC01-class (missing seq to 0)` = not_VRC01_of_GC_sup,

    `Number of IgD- B cells that are GT8++ (without regard to KO binding status)` = IgD_GT8pospos_noKO_Count,
    `Number of IgD-IgG+ B cells that are GT8++ (without regard to KO binding status)` = IgD_IgG_GT8pospos_noKO_Count,
    `Percent of GT8++ IgD- B cells that are IgG+` = IgG_of_IgD_GT8pospos_noKO,
    `Percent of GT8++ IgG+ B cells detected as VRC01-class` = VRC01_of_GT8pospos,
    `Percent of GT8++ IgG+ B cells detected as VRC01-class (missing seq to 0)` = VRC01_of_GT8pospos_sup,
    `Percent of B cells that are GT8++ (without regard to KO binding status)` = GT8pospos_noKO_of_BCell,
    `Percent of IgG+ B cells that are GT8++ (without regard to KO binding status)` = GT8pospos_noKO_of_IgG,
    `Percent of IgG+ GC B cells that are GT8++ (without regard to KO binding status)` = GT8pospos_noKO_of_GC,
    `Number of IgD- B cells that are GT8++KO-` = IgD_GT8pospos_KOneg_Count,
    `Number of IgD-IgG+ B cells that are GT8++KO-` = IgD_IgG_GT8pospos_KOneg_Count,
    `Percent of GT8++KO- IgD- B cells that are IgG+` = IgG_of_IgD_GT8pospos_KOneg,
    `Percent of IgG+ B cells that are GT8++KO-` = GT8pospos_KOneg_of_IgG,
    `Percent of GT8++IgG+ B cells that are KO-` = GT8pospos_KOneg_of_GT8pospos_noKO,

    Response,
    `Response (missing seq to 0)` = Response_sup,
    `Response not VRC01` = Response_notVRC01,
    `Response not VRC01 (missing seq to 0)` = Response_sup_notVRC01,

    `Sequence Performed` = Seq_Attempted,
    contains('Population')

  ) %>%
  mutate_at(vars(matches("Percent")), ~if_else(is.nan(.), NA_real_,.)) %>%
  arrange(PTID)



write_csv(final_results,
          file.path('data', 'combined_results', 'Flow_Seq_Results.csv'),
          na = '')





# Adding GT8 and Ep Specific calls ----------------------------------------

# Ep Spec

ep_spec_response_data <- final_results %>%
  select(PTID, Visit, Treatment,
         ep_spec = `Number of epitope-specific (KO-GT8++) IgG+ B cells`,
         igg_bcell = `Number of IgD-IgG+ B cells`
  ) %>%
  group_by(PTID) %>%
  mutate(
    ep_spec_bl = ep_spec[Visit == 'V02'],
    igg_bcell_bl = igg_bcell[Visit == 'V02']
  )

ep_spec_fisher_response_info <- ep_spec_response_data %>%
  filter(!Visit %in% c('V02','V05','V07A','V09')) %>%
  group_by(PTID, Visit, Treatment) %>%
  mutate(
    response_p_Ep_Specific = fisher.test(matrix(
      c(ep_spec, ep_spec_bl,
        igg_bcell - ep_spec, igg_bcell_bl - ep_spec_bl),
      nrow = 2
    ), alternative = "greater")$p.value
  )


E <- ConstructMIMOSAExpressionSet(
  ep_spec_response_data %>%
    filter(!Visit %in% c('V05','V07A','V09')),
  reference = Visit %in% "V02",
  measure.columns = c("ep_spec", "igg_bcell"),
  other.annotations = c("Visit", "PTID"),
  default.cast.formula = component ~ PTID + Visit + Treatment,
  .variables = .(PTID, Treatment),
  featureCols = 1,
  ref.append.replace = "_REF")


result <- MIMOSA(igg_bcell + ep_spec ~ PTID + Treatment | Visit,
                 data = E,
                 ref = RefTreat %in% 'Reference',
                 subset = RefTreat %in% 'Treatment',
                 method = "mcmc",
                 burn = 5000,
                 iter = 20000, seed = 538157028)



MIMOSA_response_prob <- getZ(result)[, 'Pr.response']


ep_spec_MIMOSA_info <- bind_cols(
  map_dfr(result, ~.x@result@phenoData@data) %>%
    filter(RefTreat == 'Treatment') %>%
    as_tibble() %>%
    select(-RefTreat),
  MIMOSA_response_prob_Ep_Specific = MIMOSA_response_prob,
  Response_Ep_Specific = as.numeric(MIMOSA_response_prob > .99)
)


# GT8


gt8_spec_response_data <- final_results %>%
  select(PTID, Visit, Treatment,
         gt8_spec = `Number of IgD-IgG+ B cells that are GT8++ (without regard to KO binding status)`,
         igg_bcell = `Number of IgD-IgG+ B cells`
  ) %>%
  group_by(PTID) %>%
  mutate(
    gt8_spec_bl = gt8_spec[Visit == 'V02'],
    igg_bcell_bl = igg_bcell[Visit == 'V02']
  ) %>%
  ungroup() %>%
  filter(!is.na(gt8_spec))

gt8_fisher_response_info <- gt8_spec_response_data %>%
  filter(!Visit %in% c('V02','V05','V07A','V09')) %>%
  group_by(PTID, Visit, Treatment) %>%
  mutate(
    response_p_GT8_Specific = fisher.test(matrix(
      c(gt8_spec, gt8_spec_bl,
        igg_bcell - gt8_spec, igg_bcell_bl - gt8_spec_bl),
      nrow = 2
    ), alternative = "greater")$p.value
  )


E <- ConstructMIMOSAExpressionSet(
  gt8_spec_response_data %>%
    filter(!Visit %in% c('V05','V07A','V09')),
  reference = Visit %in% "V02",
  measure.columns = c("gt8_spec", "igg_bcell"),
  other.annotations = c("Visit", "PTID"),
  default.cast.formula = component ~ PTID + Visit + Treatment,
  .variables = .(PTID, Treatment),
  featureCols = 1,
  ref.append.replace = "_REF")


result <- MIMOSA(igg_bcell + gt8_spec ~ PTID + Treatment | Visit,
                 data = E,
                 ref = RefTreat %in% 'Reference',
                 subset = RefTreat %in% 'Treatment',
                 method = "mcmc",
                 burn = 5000,
                 iter = 20000, seed = 538187028)



MIMOSA_response_prob <- getZ(result)[, 'Pr.response']


gt8_MIMOSA_info <- bind_cols(
  map_dfr(result, ~.x@result@phenoData@data) %>%
    filter(RefTreat == 'Treatment') %>%
    as_tibble() %>%
    select(-RefTreat),
  MIMOSA_response_prob_GT8_Specific = MIMOSA_response_prob,
  Response_GT8 = as.numeric(MIMOSA_response_prob > .99)
)



final_results_with_calls <- final_results %>%
  full_join(ep_spec_MIMOSA_info %>%
              select(PTID, Visit, Treatment, Response_Ep_Specific),
            by = c('PTID', 'Visit', 'Treatment')) %>%
  full_join(gt8_MIMOSA_info %>%
              select(PTID, Visit, Treatment, Response_GT8),
            by = c('PTID', 'Visit', 'Treatment')) %>%
  mutate(    # Overriding response calls for FNA and PB based on 0.01%
    Response_Ep_Specific =
      case_when(Visit %in% c('V05','V09') ~
                  (`Percent of IgG+ GC B cells that are epitope-specific (KO-GT8++)` > 0.1) %>%
                  as.numeric(),
                Visit %in% c('V07A') ~
                  (`Percent of IgG+ B cells that are epitope-specific (KO-GT8++)` > 0.1) %>%
                  as.numeric(),
                TRUE ~ Response_Ep_Specific),
    Response_GT8 =
      case_when(Visit %in% c('V05','V09') ~
                  (`Percent of IgG+ GC B cells that are GT8++ (without regard to KO binding status)` > 0.1) %>%
                  as.numeric(),
                Visit %in% c('V07A') ~
                  (`Percent of IgG+ B cells that are GT8++ (without regard to KO binding status)` > 0.1) %>%
                  as.numeric(),
                TRUE ~ Response_GT8)) %>%
    select(-contains('Population'), everything(), contains('Population'))

write_csv(final_results_with_calls,
          file.path('data', 'combined_results', 'Flow_Seq_Results_with_calls.csv'),
          na = '')












# # Data Summary -----------------------------------------------------------------------------------------
# 
# # Reading in specimen info
# specimen_info <- read_csv(file.path('Data_Completeness','G001_ptid_info_&_Spec_Collection_29May2020.csv'))
# 
# 
# data_summary_data <- combined_results %>%
#   mutate(
#     Visit_Tissue = case_when(
#       Visit %in% c('V05','V09') ~ paste0(Visit, ' (FNA)'),
#       Visit %in% c('V07A') ~ paste0(Visit, ' (PB)'),
#       TRUE ~ paste0(Visit, ' (Frozen)'),
#     ) %>% factor()
#   )
# 
# data_summary_ptid <-   data_summary_data %>%
#   full_join(specimen_info,
#             by = c('PTID', 'Site', 'Dose', 'Visit_Tissue' = 'name')) %>%
#   mutate(
#     Flow_Seq_Data = case_when(
#       is.na(Seq_Attempted) ~ Flow_Missing_Reason,
#       # AG_Specific_Count < sequenced_Ag_Spec ~ 'Ag Spec < # Seq',
#       Seq_Attempted %in% c('No Seq Expected',
#                            'VRC01 Calls Made') ~ 'X',
#       Seq_Attempted == 'Seq Data but no VRC01 Calls' ~ 'No VRC01 Performed',
#       Seq_Attempted == 'No Seq Performed' ~ 'No Seq Data'
#     )
#   ) %>%
#   pivot_wider(id_cols = c('PTID', 'Site'),
#               names_from = Visit_Tissue,
#               names_sort = TRUE,
#               values_from = Flow_Seq_Data) %>%
#   full_join(combined_results %>% distinct(PTID, Treatment), by = 'PTID') %>%
#   select(PTID, Site, Treatment, everything()) %>%
#   arrange(PTID)
# 
# 
# write_csv(data_summary_ptid, file.path('data', 'combined_results', 'flow_and_seq_data_summary_by_ptid.csv'))
# 
# 
# 










# By Type Results ---------------------------------------------------------


combined_results_by_type <-
  left_join(flow_data_by_type %>%
              # Dropping super low count cases
              filter(PBMC_Count > 100),
            seq_results_long %>%
              mutate(Plate_Type =
                       str_replace(Plate_Type, 'Ag_Spec', 'Antigen-Specific')),
            by = c('PTID' = 'Subject', 'Visit' = 'Timepoint', 'Tissue',
                   'INX_Population' = 'Plate_Type')) %>%
  rename('Plate_Type' = 'INX_Population') %>%
  left_join(flow_manifest_long %>%
              mutate(Plate_Type =
                       str_replace(plate_type, 'Ag_Spec', 'Antigen-Specific')),
            by = c('PTID', 'Visit', 'Tissue', 'Plate_Type')) %>%
  left_join(samples_in_paired_long %>%
              mutate(Plate_Type =
                       str_replace(Plate_Type, 'Ag_Spec', 'Antigen-Specific') %>%
                       str_replace('_Seq_Manifest', '')),
            by = c('PTID' = 'Subject', 'Visit' = 'Timepoint', 'Plate_Type')) %>%
  full_join(treat_info, by = c('PTID' = 'Volunteer Study ID')) %>%
  mutate(
    Dose = Dose_Group %>% fct_recode(`Low Dose` = 'low_dose',
                                `High Dose` = 'high_dose'),

    # making seq indicator variable
    Seq_Attempted = case_when(
      (Plate_Type == 'Antigen-Specific' & AG_Specific_Count == 0) |
        (Plate_Type == 'Bulk' & Bulk_Count == 0) ~ 'No Seq Expected',
      sequenced > 0 ~ 'VRC01 Calls Made',
      is.na(in_paired) ~ 'No Seq Performed',
      !is.na(in_paired) ~
        'Seq Data but no VRC01 Calls'
    ),
    n_not_VRC01 = sequenced - n_VRC01) %>%
  # Bringing in BL values (V02)
  group_by(PTID, Site, Dose, Plate_Type) %>%
  mutate(
    AG_Specific_Count_BL = ifelse(Plate_Type == 'Antigen-Specific',
                                  AG_Specific_Count[Visit == 'V02'],
                                  NA),
    PBMC_Count_BL = ifelse(Plate_Type == 'Antigen-Specific',
                           PBMC_Count[Visit == 'V02'],
                           NA),
    IgD_Count_BL = ifelse(Plate_Type == 'Antigen-Specific',
                          IgD_Count[Visit == 'V02'],
                          NA),
    IgD_IgG_Count_BL = ifelse(Plate_Type == 'Antigen-Specific',
                          IgD_IgG_Count[Visit == 'V02'],
                          NA),
    n_VRC01_BL = ifelse(Plate_Type == 'Antigen-Specific',
                        n_VRC01[Visit == 'V02'],
                        NA),
    n_not_VRC01_BL = ifelse(Plate_Type == 'Antigen-Specific',
                            n_not_VRC01[Visit == 'V02'],
                            NA),
    sequenced_BL = ifelse(Plate_Type == 'Antigen-Specific',
                          sequenced[Visit == 'V02'],
                          NA)
  ) %>%
  ungroup() %>%
  select(PTID, Site, Dose, Treatment, everything(),
         -Dose_Group, -Group, -plate_type) %>%
  arrange(desc(Dose))


write_csv(combined_results_by_type,
          file.path('data', 'combined_results', 'Wide_Flow_Seq_Data_by_Plate_Type.csv'),
          na = '')




final_results_by_type <- combined_results_by_type %>%
  mutate(

    #using IgD for PB
    IgD_or_IgD_IgG_Count = ifelse(Visit %in% c('V07A') & !is.na(IgD_Count),
                                  IgD_Count,IgD_IgG_Count),

    BCell_of_PBMC = 100 * BCell_Count / PBMC_Count,
    IgD_IgG_of_PBMC = 100 * IgD_or_IgD_IgG_Count / PBMC_Count,
    IgD_IgG_of_BCell = 100 * IgD_or_IgD_IgG_Count / BCell_Count,
    KOneg_of_IgG = 100 * IgG_KOneg_count / IgD_or_IgD_IgG_Count,
    IgD_IgG_GC_Count = ifelse(Visit %in% c('V05','V09'),
                          Bulk_Count,
                          NA),
    KOneg_of_GC = 100 * IgG_KOneg_count / IgD_IgG_GC_Count,

    Ag_spec_of_KOneg = 100 * AG_Specific_Count / IgG_KOneg_count,
    Ag_spec_of_IgG = 100 * AG_Specific_Count / IgD_or_IgD_IgG_Count,
    Ag_spec_of_IgG_BL = 100 * AG_Specific_Count_BL /
      ifelse(Visit %in% c('V07A'), IgD_Count_BL, IgD_IgG_Count_BL),
    Ag_spec_of_GC = 100 * AG_Specific_Count / IgD_IgG_GC_Count,
    Ag_spec_of_BCell = 100 * AG_Specific_Count / BCell_Count,
    Ag_spec_of_PBMC = 100 * AG_Specific_Count / PBMC_Count,

    Bulk_of_IgG = 100 * Bulk_Count / IgD_or_IgD_IgG_Count,
    Bulk_of_BCell = 100 * Bulk_Count / BCell_Count,
    Bulk_of_PBMC = 100 * Bulk_Count / PBMC_Count,

    seq_of_Ag_spec = 100 * sequenced /
      ifelse(Plate_Type == 'Antigen-Specific', AG_Specific_Count, Bulk_Count),
    VRC01_of_Ag_spec_seq = 100 * n_VRC01 / sequenced,
    VRC01_of_Ag_spec_seq_bl_bound =
      ifelse(Visit == 'V02' & is.na(n_VRC01),
             100,
             VRC01_of_Ag_spec_seq),

    VRC01_of_Ag_spec_seq_BL = 100 * n_VRC01_BL / sequenced_BL,
    not_VRC01_of_Ag_spec_seq_BL = 100 * n_not_VRC01_BL / sequenced_BL,

    n_not_VRC01 = sequenced - n_VRC01,
    n_not_VRC01_VH1_2 = n_hc_vh1_2_paired - n_VRC01,
    n_not_VRC01_CDRL3_5amino = n_lc_5cdrl3_paired - n_VRC01,
    not_VRC01_of_Ag_spec_seq = 100 * n_not_VRC01 / sequenced,
    not_VRC01_of_Ag_spec_seq_bl_bound =
      ifelse(Visit == 'V02' & is.na(n_not_VRC01),
             100,
             not_VRC01_of_Ag_spec_seq),
    not_VRC01_VH1_2_of_seq = 100 * n_not_VRC01_VH1_2 / sequenced,
    not_VRC01_CDRL3_5amino_of_seq = 100 * n_not_VRC01_CDRL3_5amino / sequenced,


    VRC01_of_IgG_seq_bl_bound =
      if_else((Plate_Type == 'Antigen-Specific' & AG_Specific_Count == 0) |
                (Plate_Type == 'Bulk' & Bulk_Count == 0) , 0,
              100 *
                (ifelse(Plate_Type == 'Antigen-Specific',
                        Ag_spec_of_IgG,
                        Bulk_of_IgG) / 100) *
                (VRC01_of_Ag_spec_seq_bl_bound / 100),
      ),
    # Need to remove cases where No Seq Performed
    VRC01_of_IgG_seq_bl_bound = if_else(Seq_Attempted == 'No Seq Performed',
                                        NA_real_,
                                        VRC01_of_IgG_seq_bl_bound),

    VRC01_of_IgG =
      if_else(is.na(VRC01_of_Ag_spec_seq)  & AG_Specific_Count > 0, NA_real_,
              VRC01_of_IgG_seq_bl_bound),

    VRC01_of_GC =
      if_else(AG_Specific_Count == 0, 0,
              100 * (Ag_spec_of_GC / 100) * (VRC01_of_Ag_spec_seq / 100)
      ),

    VRC01_of_BCell =
      if_else((Plate_Type == 'Antigen-Specific' & AG_Specific_Count == 0) |
                (Plate_Type == 'Bulk' & Bulk_Count == 0) , 0,
              100 *
                (ifelse(Plate_Type == 'Antigen-Specific',
                        Ag_spec_of_BCell,
                        Bulk_of_BCell) / 100) *
                (VRC01_of_Ag_spec_seq / 100),
      ),

    # Need to remove cases where No Seq Performed
    VRC01_of_BCell = if_else(Seq_Attempted == 'No Seq Performed',
                             NA_real_,
                             VRC01_of_BCell),



    not_VRC01_of_IgG_seq_bl_bound =
      if_else((Plate_Type == 'Antigen-Specific' & AG_Specific_Count == 0) |
                (Plate_Type == 'Bulk' & Bulk_Count == 0) , 0,
              100 *
                (ifelse(Plate_Type == 'Antigen-Specific',
                        Ag_spec_of_IgG,
                        Bulk_of_IgG) / 100) *
                (not_VRC01_of_Ag_spec_seq_bl_bound / 100),
      ),
    # Need to remove cases where No Seq Performed
    not_VRC01_of_IgG_seq_bl_bound = if_else(Seq_Attempted == 'No Seq Performed',
                                            NA_real_,
                                            not_VRC01_of_IgG_seq_bl_bound),

    not_VRC01_of_IgG =
      if_else(is.na(not_VRC01_of_Ag_spec_seq)  & AG_Specific_Count > 0, NA_real_,
              not_VRC01_of_IgG_seq_bl_bound),


    not_VRC01_of_BCell =
      if_else((Plate_Type == 'Antigen-Specific' & AG_Specific_Count == 0) |
                (Plate_Type == 'Bulk' & Bulk_Count == 0) , 0,
              100 *
                (ifelse(Plate_Type == 'Antigen-Specific',
                        Ag_spec_of_BCell,
                        Bulk_of_BCell) / 100) *
                (not_VRC01_of_Ag_spec_seq / 100),
      ),

    # Need to remove cases where No Seq Performed
    not_VRC01_of_BCell = if_else(Seq_Attempted == 'No Seq Performed',
                                 NA_real_,
                                 not_VRC01_of_BCell),
    not_VRC01_of_GC =
      if_else(AG_Specific_Count == 0, 0,
              100 * (Ag_spec_of_GC / 100) * (not_VRC01_of_Ag_spec_seq / 100)
      ),

    # sup vars
    sequenced_sup =  ifelse(
      Seq_Attempted %in% c("No Seq Performed", "Seq Data but no VRC01 Calls"),
      0,sequenced
    ),
    VRC01_of_BCell_sup =
      ifelse(
        Seq_Attempted %in% c("No Seq Performed", "Seq Data but no VRC01 Calls"),
        0, VRC01_of_BCell
      ),
    VRC01_of_IgG_sup =
      ifelse(
        Seq_Attempted %in% c("No Seq Performed", "Seq Data but no VRC01 Calls"),
        0, VRC01_of_IgG
      ),
    VRC01_of_GC_sup =
      ifelse(
        Seq_Attempted %in% c("No Seq Performed", "Seq Data but no VRC01 Calls"),
        0, VRC01_of_GC
      ),
    not_VRC01_of_BCell_sup =
      ifelse(
        Seq_Attempted %in% c("No Seq Performed", "Seq Data but no VRC01 Calls"),
        0, not_VRC01_of_BCell
      ),
    not_VRC01_of_IgG_sup =
      ifelse(
        Seq_Attempted %in% c("No Seq Performed", "Seq Data but no VRC01 Calls"),
        0, not_VRC01_of_IgG
      ),
    not_VRC01_of_IgG_sup =
      ifelse(
        Seq_Attempted %in% c("No Seq Performed", "Seq Data but no VRC01 Calls"),
        0, not_VRC01_of_IgG
      ),
    not_VRC01_of_GC_sup =
      ifelse(
        Seq_Attempted %in% c("No Seq Performed", "Seq Data but no VRC01 Calls"),
        0, not_VRC01_of_GC
      ),


    # Getting new gate flow estimates
    GT8pospos_noKO_of_BCell = 100 * ifelse(Visit %in% c('V07A') & !is.na(IgD_GT8pospos_noKO_Count),
                                           IgD_GT8pospos_noKO_Count / BCell_Count,
                                           IgD_IgG_GT8pospos_noKO_Count / BCell_Count),
    GT8pospos_noKO_of_IgG =  100 * ifelse(Visit %in% c('V07A') & !is.na(IgD_GT8pospos_noKO_Count),
                                          IgD_GT8pospos_noKO_Count / IgD_Count,
                                          IgD_IgG_GT8pospos_noKO_Count / IgD_IgG_Count),
    GT8pospos_noKO_of_GC = ifelse(Visit %in% c('V05','V09'),
                                  100 * IgD_IgG_GT8pospos_noKO_Count / Bulk_Count,
                                  NA),
    GT8pospos_KOneg_of_GT8pospos_noKO = 100 * ifelse(Visit %in% c('V07A') & !is.na(IgD_GT8pospos_KOneg_Count),
                                                     IgD_GT8pospos_KOneg_Count / IgD_GT8pospos_noKO_Count,
                                                     IgD_IgG_GT8pospos_KOneg_Count / IgD_IgG_GT8pospos_noKO_Count),
    GT8pospos_KOneg_of_IgG = 100 * ifelse(Visit %in% c('V07A') & !is.na(IgD_GT8pospos_KOneg_Count),
                                          IgD_GT8pospos_KOneg_Count / IgD_Count,
                                          IgD_IgG_GT8pospos_KOneg_Count / IgD_IgG_Count),

    IgG_of_IgD = 100 * IgD_IgG_Count / IgD_Count,
    IgG_of_IgD_GT8pospos_noKO = 100 * IgD_IgG_GT8pospos_noKO_Count / IgD_GT8pospos_noKO_Count,
    IgG_of_IgD_GT8pospos_KOneg = 100 * IgD_IgG_GT8pospos_KOneg_Count / IgD_GT8pospos_KOneg_Count,

    VRC01_of_GT8pospos =
      if_else(AG_Specific_Count == 0,  0,
              100 * pmin(AG_Specific_Count / ifelse(Visit %in% c('V07A'),
                                                    IgD_GT8pospos_noKO_Count,
                                                    IgD_IgG_GT8pospos_noKO_Count),
                         1) * (VRC01_of_Ag_spec_seq_bl_bound / 100)),
    # Need to remove cases where No Seq Performed
    VRC01_of_GT8pospos = if_else(Seq_Attempted == 'No Seq Performed',
                                 NA_real_,
                                 VRC01_of_GT8pospos),
    VRC01_of_GT8pospos_sup =
      ifelse(
        Seq_Attempted %in% c("No Seq Performed", "Seq Data but no VRC01 Calls"),
        0, VRC01_of_GT8pospos
      ),


    # Getting BL Value for response calls
    VRC01_of_IgG_BL = 100 * (Ag_spec_of_IgG_BL / 100) *
      if_else(is.na(sequenced_BL) | sequenced_BL == 0,
              1, VRC01_of_Ag_spec_seq_BL / 100),
    not_VRC01_of_IgG_BL = 100 * (not_VRC01_of_Ag_spec_seq_BL / 100) *
      if_else(is.na(sequenced_BL) | sequenced_BL == 0,
              1, VRC01_of_Ag_spec_seq_BL / 100),

    # response def
    Response = case_when(
      Tissue == 'Fresh' ~ as.numeric(n_VRC01 > 0),
      n_VRC01 > 0 & VRC01_of_IgG > VRC01_of_IgG_BL ~ 1,
      !is.na(VRC01_of_IgG) & VRC01_of_IgG <= VRC01_of_IgG_BL ~ 0,
      TRUE ~ NA_real_
    ),
    Response_sup = ifelse(
      Seq_Attempted %in% c("No Seq Performed", "Seq Data but no VRC01 Calls"),
      0, Response
    ),

    Response_notVRC01 = case_when(
      Tissue == 'Fresh' ~ as.numeric(n_not_VRC01 > 0),
      n_not_VRC01 > 0 & not_VRC01_of_IgG > not_VRC01_of_IgG_BL ~ 1,
      !is.na(not_VRC01_of_IgG) & not_VRC01_of_IgG <= not_VRC01_of_IgG_BL ~ 0,
      TRUE ~ NA_real_
    ),
    Response_sup_notVRC01 = ifelse(
      Seq_Attempted %in% c("No Seq Performed", "Seq Data but no VRC01 Calls"),
      0, Response_notVRC01
    ),


  ) %>%
  select(
    PTID,
    Site,
    Treatment,
    Visit,
    Tissue,
    `Plate Type` = Plate_Type,

    `Number of PBMCs sorted` = root_count,
    `Number of PBMCs sorted (Lymphocytes/Singlets)` = PBMC_Count,
    `Number of B cells` = BCell_Count,
    `Percent of PBMCs (Lymphocytes/Singlets) that are B cells` = BCell_of_PBMC,
    `Number of IgD- B cells` = IgD_Count,
    `Number of IgD-IgG+ B cells` = IgD_IgG_Count,
    `Percent of IgD- B cells that are IgG+` = IgG_of_IgD,
    `Number of Bulk IgG+ B cells` = Bulk_Count,
    `Percent of B cells that are IgG+ B cells` = IgD_IgG_of_BCell,
    `Number of IgG+ GC B cells` = IgD_IgG_GC_Count,
    `Number of IgG+ B cells that are KO-` = IgG_KOneg_count,
    `Percent of IgG+ B cells that are KO-` = KOneg_of_IgG,
    `Number of epitope-specific (KO-GT8++) IgG+ B cells` = AG_Specific_Count,
    `Percent of IgG+ KO- B cells that are GT8++` = Ag_spec_of_KOneg,
    `Percent of B cells that are epitope-specific (KO-GT8++)` = Ag_spec_of_BCell,
    `Percent of IgG+ B cells that are epitope-specific (KO-GT8++)` = Ag_spec_of_IgG,
    `Percent of IgG+ GC B cells that are epitope-specific (KO-GT8++)` = Ag_spec_of_GC,

    `Number of epitope-specific (KO-GT8++)/Bulk IgG+ B cells that have BCR heavy and light chains sequenced` = sequenced,
    `Number of epitope-specific (KO-GT8++)/Bulk IgG+ B cells that have BCR heavy and light chains sequenced (missing seq to 0)` = sequenced_sup,
    `Percent of epitope-specific (KO-GT8++)/Bulk IgG+ B cells that have BCR heavy and light chains sequenced` = seq_of_Ag_spec,

    `Number of epitope-specific (KO-GT8++)/Bulk sequenced IgG BCRs that are VRC01-class` = n_VRC01,
    `Percent of epitope-specific (KO-GT8++) sequenced IgG BCRs that are VRC01-class` = VRC01_of_Ag_spec_seq,

    `Number of epitope-specific (KO-GT8++)/Bulk sequenced IgG BCRs that are not VRC01-class` = n_not_VRC01,
    `Number epitope-specific (KO-GT8++)/Bulk sequenced IgG BCRs that are not VRC01-class but have a VH1-2 heavy chain` = n_not_VRC01_VH1_2,
    `Number epitope-specific (KO-GT8++)/Bulk sequenced IgG BCRs that are not VRC01-class but have a 5-aa CDRL3` = n_not_VRC01_CDRL3_5amino,
    `Percent of epitope-specific (KO-GT8++)/Bulk sequenced IgG BCRs that are not VRC01-class` = not_VRC01_of_Ag_spec_seq,
    `Percent of epitope-specific (KO-GT8++)/Bulk sequenced IgG BCRs that are not VRC01-class but have a VH1-2 heavy chain` = not_VRC01_VH1_2_of_seq,
    `Percent of epitope-specific (KO-GT8++)/Bulk sequenced IgG BCRs that are not VRC01-class but have a 5-aa CDRL3` = not_VRC01_CDRL3_5amino_of_seq,

    `Percent of B cells detected as VRC01-class` = VRC01_of_BCell,
    `Percent of B cells detected as VRC01-class (missing seq to 0)` = VRC01_of_BCell_sup,
    `Percent of IgG+ B cells detected as VRC01-class` = VRC01_of_IgG,
    `Percent of IgG+ B cells detected as VRC01-class (Baseline Upper Bound)` = VRC01_of_IgG_seq_bl_bound,
    `Percent of IgG+ B cells detected as VRC01-class (missing seq to 0)` = VRC01_of_IgG_sup,

    `Percent of IgG+ GC B cells detected as VRC01-class` = VRC01_of_GC,
    `Percent of IgG+ GC B cells detected as VRC01-class (missing seq to 0)` = VRC01_of_GC_sup,

    `Percent of B cells that are not VRC01-class` = not_VRC01_of_BCell,
    `Percent of B cells that are not VRC01-class (missing seq to 0)` = not_VRC01_of_BCell_sup,
    `Percent of IgG+ B cells that are not VRC01-class` = not_VRC01_of_IgG,
    `Percent of IgG+ B cells that are not VRC01-class (Baseline Upper Bound)` = not_VRC01_of_IgG_seq_bl_bound,
    `Percent of IgG+ B cells that are not VRC01-class (missing seq to 0)` = not_VRC01_of_IgG_sup,

    `Percent of IgG+ GC B cells that are not VRC01-class` = not_VRC01_of_GC,
    `Percent of IgG+ GC B cells that are not VRC01-class (missing seq to 0)` = not_VRC01_of_GC_sup,

    `Number of IgD- B cells that are GT8++ (without regard to KO binding status)` = IgD_GT8pospos_noKO_Count,
    `Number of IgD-IgG+ B cells that are GT8++ (without regard to KO binding status)` = IgD_IgG_GT8pospos_noKO_Count,
    `Percent of GT8++ IgD- B cells that are IgG+` = IgG_of_IgD_GT8pospos_noKO,
    `Percent of GT8++ IgG+ B cells detected as VRC01-class` = VRC01_of_GT8pospos,
    `Percent of GT8++ IgG+ B cells detected as VRC01-class (missing seq to 0)` = VRC01_of_GT8pospos_sup,
    `Percent of B cells that are GT8++ (without regard to KO binding status)` = GT8pospos_noKO_of_BCell,
    `Percent of IgG+ B cells that are GT8++ (without regard to KO binding status)` = GT8pospos_noKO_of_IgG,
    `Percent of IgG+ GC B cells that are GT8++ (without regard to KO binding status)` = GT8pospos_noKO_of_GC,
    `Number of IgD- B cells that are GT8++KO-` = IgD_GT8pospos_KOneg_Count,
    `Number of IgD-IgG+ B cells that are GT8++KO-` = IgD_IgG_GT8pospos_KOneg_Count,
    `Percent of GT8++KO- IgD- B cells that are IgG+` = IgG_of_IgD_GT8pospos_KOneg,
    `Percent of IgG+ B cells that are GT8++KO-` = GT8pospos_KOneg_of_IgG,
    `Percent of GT8++IgG+ B cells that are KO-` = GT8pospos_KOneg_of_GT8pospos_noKO,

    Response,
    `Response (missing seq to 0)` = Response_sup,
    `Response not VRC01` = Response_notVRC01,
    `Response not VRC01 (missing seq to 0)` = Response_sup_notVRC01,

    `Sequence Performed` = Seq_Attempted,
    contains('Population')

  ) %>%
  mutate_at(vars(matches("Percent")), ~if_else(is.nan(.), NA_real_,.)) %>%
  arrange(PTID)


write_csv(final_results_by_type,
          file.path(output_path, 'Flow_Seq_Results_by_Plate_Type.csv'),
          na = '')




# # Reading in specimen info
# 
# 
# 
# data_summary_data_by_type <- combined_results_by_type %>%
#   mutate(
#     Visit_Tissue = case_when(
#       Visit %in% c('V05','V09') ~ paste0(Visit, ' (FNA)'),
#       Visit %in% c('V07A') ~ paste0(Visit, ' (PB)'),
#       TRUE ~ paste0(Visit, ' (Frozen)'),
#     ) %>% factor()
#   )
# 
# data_summary_ptid_by_type <-   data_summary_data_by_type %>%
#   filter(Plate_Type == 'Antigen-Specific') %>%
#   full_join(specimen_info,
#             by = c('PTID', 'Site', 'Dose', 'Visit_Tissue' = 'name')) %>%
#   mutate(
#     Flow_Seq_Data = case_when(
#       is.na(Seq_Attempted) & !is.na(Flow_Missing_Reason) ~ Flow_Missing_Reason,
#       # AG_Specific_Count < sequenced_Ag_Spec ~ 'Ag Spec < # Seq',
#       Seq_Attempted %in% c('No Seq Expected',
#                            'VRC01 Calls Made') ~ 'X',
#       Seq_Attempted == 'Seq Data but no VRC01 Calls' ~ 'No VRC01 Performed',
#       Seq_Attempted == 'No Seq Performed' ~ 'No Seq Data',
#       TRUE ~ 'No Flow Data During Sorting'
#     )
#   ) %>%
#   pivot_wider(id_cols = c('PTID', 'Site'),
#               names_from = Visit_Tissue,
#               names_sort = TRUE,
#               values_from = Flow_Seq_Data) %>%
#   full_join(combined_results %>% distinct(PTID, Treatment), by = 'PTID') %>%
#   select(PTID, Site, Treatment, unique(data_summary_data_by_type$Visit_Tissue)) %>%
#   arrange(PTID)
# 
# write_csv(data_summary_ptid_by_type,
#           file.path('data', 'combined_results',
#                     'flow_and_seq_data_summary_by_ptid_ag_spec_plate.csv'))

