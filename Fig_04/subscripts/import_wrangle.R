## ------------------------------------------
## Wrangle genotyping data
## ------------------------------------------
## Filter fileset for APOE defining SNPs and convert to vcf using bash
# plink \
#   --snps rs429358, rs7412 \
#   --bed ../data/processed/GWAS/ukb_cal_chr19_v2.bed \
#   --bim ../data/processed/GWAS/ukb_snp_bim/ukb_snp_chr19_v2.bim \
#   --fam ../data/processed/GWAS/ukb62709_cal_chr19_v2_s488264.fam \
#   --recode vcf \
#   --output-missing-genotype '.'  \
#   --out ../data/processed/GWAS/01_rs429358_rs7412_only/rs429358_rs7412_only


if(!file.exists("../data/UKB/UKB_genotypes.rds")) {
  
  genotypes_raw <- fread("../data/UKB/rs429358_rs7412_only.vcf", 
                         sep = "\t", header = TRUE, skip = 6)
  genotypes <- data.table::transpose(genotypes_raw, keep.names = "eid")
  colnames(genotypes) <- c("eid", "rs429358", "rs7412")
  genotypes <- genotypes[-c(1:9), ]
  genotypes %>% 
    lazy_dt() %>% 
    group_by(rs7412, rs429358) %>% 
    count() %>% 
    as_tibble()
  
  genotypes_UKB <-
    genotypes %>%
    lazy_dt() %>%
    mutate(variant = case_when(rs429358 == "0/0" & rs7412 == ("0/0") ~ "E3;E3",
                               rs429358 == "0/0" & rs7412 == ("0/1") ~ "E2;E3",
                               rs429358 == "0/0" & rs7412 == ("1/1") ~ "E2;E2",
                               rs429358 == "0/1" & rs7412 == ("0/0") ~ "E3;E4",
                               rs429358 == "1/1" & rs7412 == ("0/0") ~ "E4;E4",
                               rs429358 == "0/1" & rs7412 == ("0/1") ~ "E2;E4"), 
           variant_carrier = case_when(variant == "E2;E2" ~ "E2 carrier",
                                       variant == "E3;E3" ~ "E3;E3",
                                       variant == "E4;E4" ~ "E4 carrier",
                                       variant == "E3;E4" ~ "E4 carrier",
                                       variant == "E2;E3" ~ "E2 carrier",
                                       variant == "E2;E4" ~ "E2;E4"), 
           eid = gsub("(\\d*)_\\d*", "\\1", eid)) %>%
    select(-rs429358, -rs7412) %>%
    as.data.table()
  
  saveRDS(genotypes_UKB, "../data/UKB/UKB_genotypes.rds")
}
genotypes_UKB <- readRDS("../data/UKB/UKB_genotypes.rds")


## ------------------------------------------
## Read in COVID-19 clinical data
## ------------------------------------------
# For details, see [here](http://biobank.ndph.ox.ac.uk/ukb/exinfo.cgi?src=COVID19_tests). 
# - specdate:	Date the specimen was taken	
# - spectype: Specimen type as recorded on the laboratory request form (e.g. nasal, nose and throat, sputum); encoding [here](http://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=1853)
# - laboratory: Laboratory that UKB the sample
# - origin: Whether the patient was an inpatient when the sample was taken. This is based on information provided on the specimen request form. (0, no evidence that inpatient; 1, evidence that inpatient)
# - result: Whether the sample was reported as positive or negative for SARS-CoV-2 (0, negative; 1, positive)

## ----------------------------------------
## Load Covid19 results
## ----------------------------------------
covid19_res <- 
  read_tsv("../data/UKB/covid19_result.txt") %>%
  as_tibble() %>%
  mutate(origin = case_when(origin == 0 ~ "outpatient", 
                            origin == 1 ~ "inpatient"), 
         result = case_when(result == 0 ~ "negative", 
                            result == 1 ~ "positive"), 
         eid = as.character(eid), 
         specdate = as.Date(specdate, format = "%d/%m/%Y")) 

## ----------------------------------------
## Get positive COVID19 cases and wrangle
## ----------------------------------------
covid19_pos_raw <- 
  covid19_res %>%
  filter(result == "positive") %>%
  
  ## Account for multiple tests per person
  group_by(eid) %>%
  
  ## For origin, use "inpatient" if any test taken as inpatient; for specdate take earliest positive date
  summarize(eid = eid[1], 
            result = result[1], 
            origin = ifelse("inpatient" %in% origin, "inpatient", "outpatient"), 
            specdate = min(specdate))

## ----------------------------------------
## Wrangle negative cases
## ----------------------------------------
covid19_neg <-
  covid19_res %>%
  filter(result == "negative") %>%
  group_by(eid) %>%
  summarize(eid = eid[1], 
            result = result[1], 
            origin = ifelse("inpatient" %in% origin, "inpatient", "outpatient"), 
            specdate = min(specdate))

## ----------------------------------------
## Bind pos and neg cases
## ----------------------------------------
## By binding pos cases first, `distinct` will preserve positive cases and eliminate negative cases if duplicated eid
cov19_bind <- 
  rbind(covid19_pos_raw, covid19_neg) %>%
  distinct(eid, .keep_all = TRUE) %>%
  full_join(., genotypes_UKB) %>%
  mutate(population = "UKB", 
         pos_vs_all = case_when(result == "positive" ~ "positive",
                                TRUE ~ "negative/untested"), 
         tested = case_when(is.na(result) ~ "no", 
                            TRUE ~ "yes")) %>%
  as.data.table()


## ------------------------------------------
## Add survival data
## ------------------------------------------
#First case recorded in UK: Jan 31, 2020

## Read death data; some records are duplicated with different death dates, as outlined in `./manuals/DeathLinkage.pdf`
death <- 
  read_tsv("../data/UKB/death.txt") %>%
  select(eid, date_of_death) %>%
  distinct(eid, .keep_all = TRUE)

## Table with death causes; multiple ICD10 codes can be present for one patient; 
## the main cause of death is assigned level = 1 (and will have `arr_index` = 0 in the table!)
death_cause <- 
  read_tsv("../data/UKB/death_cause.txt") %>%
  ## filter for only primary death causes
  filter(arr_index == 0) %>%
  mutate(OS.covid = case_when(cause_icd10 %in% c("U071", "U072") ~ 1, 
                              TRUE ~ 0), 
         OS.other = 1) %>%
  select(eid, cause_icd10, OS.covid, OS.other) %>%
  distinct(eid, .keep_all = TRUE)

survival <- 
  inner_join(death, death_cause, by = "eid") %>%
  mutate(eid = as.character(eid)) %>%
  as.data.table()

cov19 <- 
  merge(cov19_bind, survival, all.x = TRUE) %>%
  lazy_dt() %>%
  mutate(date_of_death = as.Date(date_of_death, format = "%d/%m/%Y")) %>%
  mutate(date_death_censor = case_when(is.na(date_of_death) ~ max(cov19_bind$specdate, na.rm = TRUE), 
                                       TRUE ~ date_of_death), 
         OS.covid = case_when(is.na(OS.covid) ~ 0, 
                              TRUE ~ OS.covid), 
         OS.other = case_when(is.na(OS.other) ~ 0, 
                              TRUE ~ OS.other)) %>% 
  as.data.table() %>%
  .[, OS.time := date_death_censor - specdate] %>%
  .[, OS.time.other := date_death_censor - as.Date("2019-01-01")] %>%
  .[OS.time < 0, OS.time := NA]  %>%
  .[OS.time.other < 0, OS.time.other := NA]


if (!file.exists("../data/UKB/UKB_filt.rds")){
  
  if (!file.exists("../data/UKB/UKB_main.rds")){
    UKB_main <- ukb_df("ukb42116")
    saveRDS(UKB_main, file = "../data/UKB/UKB_main.rds")
  }
  
  UKB_main <- readRDS("../data/UKB/UKB_main.rds")
  
  ## Filter UKB main dataframe (also see file `./select_categories.md`)
  categories_keep <- paste0(c("eid", "f31_0_", "f33_0_", "f34_0_", "f23_0_", "f50_0_", "f52_0_", "f53_0_", "f54_0_", "f55_0_", "f84_0_", "f87_0_", "f92_0_", "f93_0_", "f95_0_", "f134_0_", "f135_0_", "f136_0_", "f137_0_", "f189_0_", "f190_0_", "f1717_0_", "f1727_0_", "f2443_0_", "f2453_0_", "f2473_0_", "f2644_0_", "f2734_0_", "f2966_0_", "f3062_0_", "f3063_0_", "f4079_0_", "f4080_0_", "f6177_0_", "f12143_0_", "f12144_0_", "f20001_0_", "f20002_0_", "f20006_0_", "f20007_0_", "f20008_0_", "f20161_0_", "f21001_0_", "f21022_0_", "f22001_0_", "f22006_0_", "f22009_0_", "f23104_0_", "f30000_0_", "f30003_0_", "f30010_0_", "f30020_0_", "f30080_0_", "f30120_0_", "f30130_0_", "f30140_0_", "f30150_0_", "f30160_0_", "f30250_0_", "f30510_0_", "f30630_0_", "f30690_0_", "f30710_0_", "f30790_0_", "f30800_0_", "f30870_0_", "_f40000_", "f40005_0_", "f40006_0_", "f40007_0_", "f40008_0_", "f40009_0_", "f40011_0_", "f40012_0_", "f40013_0_", "f40018_0_", "f40019_0_", "f41202_0_", "f41270_0_", "f412800_", "f22009", "f22020"), collapse="|")
  
  categories_logical <- grepl(categories_keep, x = colnames(UKB_main))
  UKB_filt <- UKB_main[, categories_logical]
  UKB_filt <- as.data.table(UKB_filt)
  UKB_filt[, ':='(eid = as.character(eid),
                  sex = as.factor(sex_f31_0_0), 
                  DOB = as.Date(paste0("01", " ", month_of_birth_f52_0_0, " ", year_of_birth_f34_0_0), format = "%d %B %Y"), 
                  ethnicity = as.factor(genetic_ethnic_grouping_f22006_0_0))]
  
  saveRDS(UKB_filt, file = "../data/UKB/UKB_filt.rds", compress = FALSE)  
}

## ===========================================================
## Merge data and subset
## ===========================================================
UKB_filt <- readRDS("../data/UKB/UKB_filt.rds")
df_all <- merge(UKB_filt, cov19, all = TRUE, by = "eid")
df_all[, tested := fcase(is.na(result), "no",
                         default = "yes")]
df_all[, variant := as.character(variant)]
df_all[, variant_carrier := as.character(variant_carrier)]
df_all[, sex := as.character(sex)]
df_all[, result := factor(result)]
df_all[, origin := ordered(origin, levels =  c("outpatient", "inpatient"))]
df_all[, pos_vs_all := fcase(result == "positive", "positive",
                             default = "negativeUntested")]
df_all[, pos_vs_all := factor(pos_vs_all, ordered = FALSE)] 

df <- df_all[variant_carrier %chin% c("E2 carrier", "E3;E3", "E4 carrier")]
df[, variant := relevel(factor(variant, ordered = FALSE), ref = "E3;E3")]
df[, variant_carrier := relevel(factor(variant_carrier, ordered = FALSE), ref = "E3;E3")]

df_all[, variant := relevel(factor(variant, ordered = FALSE), ref = "E3;E3")]
df_all[, variant_carrier := relevel(factor(variant_carrier, ordered = FALSE), ref = "E3;E3")]


## ------------------------------------------
## Prepare data.table for positive COV19 patients
## ------------------------------------------
cov19_pos <- df[result == "positive"] 
cov19_pos[, age_at_inf :=  as.numeric(difftime(specdate, DOB, unit = "weeks"))/52.25]
cov19_pos[, age_bin := fcase(age_at_inf <= summary(cov19_pos$age_at_inf)[3], "young", 
                             age_at_inf > summary(cov19_pos$age_at_inf)[3], "old")]

cov19_pos_all <- df_all[result == "positive"]
cov19_pos_all[, age_at_inf :=  as.numeric(difftime(specdate, DOB, unit = "weeks"))/52.25]
cov19_pos_all[, age_bin := fcase(age_at_inf <= summary(cov19_pos$age_at_inf)[3], "young", 
                             age_at_inf > summary(cov19_pos$age_at_inf)[3], "old")]

## ===========================================================
## Create DTs for different populations
## ===========================================================
## Independent normal population (ARIC study; Blair et al., Neurology, 2005)
carrier_freq_ARIC <- c(rep("E2 carrier", 834), rep("E3;E3", 3648),
                       rep("E4 carrier", 1561), rep("E2;E4", 159))

variant_freq_ARIC <- c(rep("E2;E2", 41), rep("E2;E3", 793),rep("E3;E3", 3648), 
                       rep("E3;E4", 1434), rep("E4;E4", 127), rep("E2;E4", 159))

df_genotypes_ARIC <- data.table(variant_carrier = carrier_freq_ARIC,
                                variant = variant_freq_ARIC, 
                                tested = NA, 
                                result = NA, 
                                pos_vs_all = NA, 
                                population = "ARIC")

## Bind genotype freq DTs
df_genotypes_freqs <- 
  rbindlist(list(df_genotypes_ARIC, 
                 df[, .(variant_carrier, variant, tested, result, population, pos_vs_all)]), use.names = TRUE) %>%
  .[, variant := as.character(variant)] %>%
  
  ## Exclude E2;E4 patients
  .[!is.na(variant) & variant != "E2;E4"] %>%
  as_tibble()

cov19_unordered <-
  cov19_pos %>%
  .[, variant := factor(variant, ordered = FALSE, levels = APOE_variant_levels)]
