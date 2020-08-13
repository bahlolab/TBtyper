# TBtyper

### **Notice**
`TBtyper` is still in the development stages and and should not be used in production. 
Documentation is currently poor and functionality is subject to change.

### Installation
`TBtyper` can be installed with `devtools` as follows:
```R 
devtools::install_github("bahlolab/TBtyper")
```
  
### Basic usage
1) Convert VCF file to GDS file using `SeqArray`
2) Extract allele counts from GDS
3) Fit phylotypes using TBtyper
4) Filter results
5) Save Results
```R
library(tidyverse)
library(SeqArray)
library(TBtyper)
library(future)

# set desired parallelism with future::plan
plan(multiprocess, workers = 4)

# 1) convert VCF file to GDS file
# Note: Variants should be called against the H37RV genome
vcf_fn <- '/path/to/my.vcf.gz'
gds_fn <- '/path/to/my.gds'
seqVCF2GDS(vcf_fn, gds_fn, storage.option = 'ZIP_RA', ignore.chr.prefix = NA_character_)
gds <- seqOpen(gds_fn, allow.duplicate = TRUE)

# 2) Extract allele counts from GDS
allele_counts <- get_allele_counts_gds(gds)

# 3) Fit phylotypes with TBtyper
results <- fit_phylotypes(allele_counts)

# 4) Filter for best match per sample
filtered_results <- 
  results %>% 
  filter(map_lgl(mix_prop, ~ min(.) > 0.01),
         abs_diff > 5,
         p_val_wsrst < 0.001,
         mix_n <= 3) %>% 
  group_by(sample_id) %>% 
  slice(which.max(mix_n)) %>% 
  ungroup() %>% 
  unnest_legacy() %>% 
  select(sample_id, mix_n, mix_prop, phylotype)
  
# 5) Save CSV
write_csv(filtered_results, 'TBtyper_results_filtered.csv')
```
