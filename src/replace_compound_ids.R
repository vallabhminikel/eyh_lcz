setwd('~/d/sci/src/eyh_lcz')
moa = read_csv('data/moa_annotations.csv')
moa %>%
  distinct(nvp, inchi_key) -> nvp_inchi_map
write_tsv(nvp_inchi_map, 'ignore/nvp_inchi_map.tsv')
