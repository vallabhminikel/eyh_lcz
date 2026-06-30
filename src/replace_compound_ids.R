setwd('~/d/sci/src/eyh_lcz')
moa = read_csv('data/moa_annotations.csv')
moa %>%
  distinct(nvp, inchi_key) -> nvp_inchi_map
write_tsv(nvp_inchi_map, 'ignore/nvp_inchi_map.tsv')

# Sanitize a picks table: replace proprietary NVP compound IDs with InChI keys
nvp_inchi_map = read_tsv('ignore/nvp_inchi_map.tsv')
picks = read_tsv('../kd_moa_legacy/output/erics_picks_v2.tsv')
picks %>%
  left_join(nvp_inchi_map, by = c('compound' = 'nvp')) %>%
  mutate(compound = inchi_key) %>%
  select(-inchi_key, -status, -target) -> picks_sanitized
write_csv(picks_sanitized, 'data/4pt_picks.csv')
