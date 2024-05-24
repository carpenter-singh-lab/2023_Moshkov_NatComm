# Script to extract assay metadata from internal resources

library(tidyverse)

# https://imaging-platform.s3.amazonaws.com/projects/2018_01_09_PUMA_CBTS/workspace/analysis/cdrp/assay_readout.csv
assay_readout <- 
  read_csv('~/work/projects/2018_01_09_PUMA_CBTS/workspace/analysis/cdrp/assay_readout.csv')

# https://imaging-platform.s3.amazonaws.com/projects/2018_01_09_PUMA_CBTS/workspace/analysis/cdrp/Broad_MLPCN_pubchem-bioassay-aids.csv
bioassay_aids <- 
  read_csv('~/work/projects/2018_01_09_PUMA_CBTS/workspace/analysis/cdrp/Broad_MLPCN_pubchem-bioassay-aids.csv')

# https://raw.githubusercontent.com/CaicedoLab/2023_Moshkov_NatComm/20fd97178796148c2b464ddc2c025a91654a6433/assay_data/assay_metadata.csv
assay_metadata <- 
  read_csv("~/Documents/GitHub/2023_Moshkov_NatComm/assay_data/assay_metadata.csv")

assay_metadata_expanded <-
  assay_readout %>%
  select(-RESULT_VALUE, -ASSAY_UMOL_CONC, -CPD_ID, -n_compounds) %>%
  distinct() %>%
  unite(PUMA_ASSAY_ID, c(ASSAY_ID, ASSAY_OBS_ID), remove = TRUE)

# remove emails when writing output
assay_metadata_expanded <- 
  assay_metadata_expanded %>%
  mutate(ASSAY_DESC = str_remove_all(ASSAY_DESC, "\\([A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\\.[A-Za-z]{2,}\\)"))
  
assay_metadata_expanded <-
  assay_metadata_expanded %>%
  inner_join(assay_metadata %>% select(PUMA_ASSAY_ID, ASSAY_TYPE, READOUT_TYPE),
             by = "PUMA_ASSAY_ID") %>%
  separate(PUMA_ASSAY_ID, c("x", "y"), remove = FALSE) %>%
  arrange(x, y) %>%
  select(-x, -y)

# list of assays submitted by Broad Institute during the MLPCN screening campaigns 
# This helps map assays to the submissions in Pubchem, using the `AID` column
# For RegIDs that only list the seven character string for the assay code, e.g, `2013-01`, that is the summary assay for that project.

bioassay_aids <- 
  bioassay_aids %>%
  mutate(ASSAY_DB_ID = str_c("CBIP:", str_extract(RegID, "^\\d{4}(?:-\\d{2})?"))) %>%
  inner_join(assay_metadata_expanded %>% distinct(ASSAY_DB_ID), by = "ASSAY_DB_ID")

bioassay_aids %>% write_csv("bioassay_aids.csv")

assay_metadata_expanded %>% write_csv("assay_metadata_expanded.csv")
