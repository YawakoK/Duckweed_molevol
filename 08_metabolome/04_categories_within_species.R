# Step 04: Functional categories, within-species HL > LL tests and the
#          category-level table for Fig. 6a-e (Supplementary Data 22).
#
# Input : output/SupplementaryData21_all_metabolite_log2FC_statistics.csv (step 03)
#         output/metabolite_log2_abundance_long.csv (step 03)
# Output: output/metabolite_category_annotation.csv (category of all 1,730 metabolites)
#         output/SupplementaryData22_category_metabolite_log2FC_statistics.csv
# Run from this directory: Rscript 04_categories_within_species.R
# Author: Natsu Katayama
#
# Categories are assigned from the metabolite names with the platform-specific
# regular expressions below. Fig. 6 panels:
#   a  LC-MS/MS  Aromatic / Phenylpropanoid, flavonoid subcategories
#   b  GC-MS/MS  Amino acids / Peptides
#   c  GC-MS/MS  Central carbon metabolism (glycolysis + TCA)
#   d  GC-MS/MS  Carbohydrates / Sugar metabolism + Polyols
#   e  CE-TOF MS Energy / Redox
# Within-species test: one-sided Welch t-test of log2y, HL > LL, for each
# metabolite and species; BH correction across all metabolites of the platform,
# separately for each species. A metabolite is "significantly increased under
# HL" in a species when log2FC > 0 and this FDR < 0.05.

library(dplyr)
library(tidyr)
library(readr)
library(stringr)

supp21 <- read_csv("output/SupplementaryData21_all_metabolite_log2FC_statistics.csv",
                   show_col_types = FALSE)
log2_long <- read_csv("output/metabolite_log2_abundance_long.csv", show_col_types = FALSE)

# ---- display labels ----------------------------------------------------------
make_label_name <- function(x) {
  x |>
    str_replace("^T-\\d+_?", "") |>                          # CE-TOF MS T-xxx_ prefix
    str_replace("^\\+_", "") |>                               # CE-TOF MS ion prefix
    str_replace_all("-[A-Za-z0-9]+TMS\\b", "") |>             # GC-MS/MS TMS derivative
    str_replace_all("-?meto\\b", "") |>                       # GC-MS/MS methoxime
    str_replace_all("\\([0-9]+\\)$", "") |>                   # peak/isomer index
    str_replace_all("\\[M[-+][^\\]]+\\]|\\[M\\][+-]?", "") |> # adduct notation
    str_replace_all("\\*+", "") |>
    str_replace_all("\\s*;\\s*$", "") |>
    str_squish() |>
    str_replace_all(";\\s*", "/")
}

clean_lcms_name <- function(x) {
  x |>
    str_replace_all("\\s+", " ") |>
    str_squish() |>
    str_remove(";\\s*$") |>
    word(1, sep = fixed(";"))
}

# ---- category definitions ----------------------------------------------------
gcms_photorespiration_regex <- "^(Glycine|Serine|Glyceric acid|Glycerate|Glycolic acid|Glycolate)(-\\d+TMS)?$"
gcms_central_carbon_regex <- "^(Citric acid|Isocitric acid|Aconitic acid|Malic acid|Fumaric acid|Succinic acid|Pyruvic acid|Lactic acid|2-oxoglutaric acid|alpha-ketoglutaric acid)(-\\d+TMS)?$"
flavonoid_regex <- "flavonoid|flavone|flavonol|flavan|flavanone|flavanol|isoflav|apigenin|luteolin|quercetin|kaempferol|anthocyan|anthocyanidin|cyanidin|delphinidin|pelargonidin|peonidin|petunidin|malvidin|riccionidin|salvianin|procyanidin|procyanidol|proanthocyanidin|catechin|epicatechin|gallocatechin|theasinensin|gambiriin|naringenin|hesperetin|eriodictyol|formononetin|chrysin|baicalein|myricetin|rutin"
anthocyanin_regex <- "anthocyan|anthocyanidin|cyanidin|delphinidin|pelargonidin|peonidin|petunidin|malvidin|riccionidin|salvianin"
proanthocyanidin_regex <- "procyanidin|procyanidol|proanthocyanidin|catechin|epicatechin|gallocatechin|theasinensin|gambiriin"
flavonol_flavone_regex <- "flavonol|flavone|apigenin|luteolin|quercetin|kaempferol|myricetin|nevadensin|formononetin|chrysin|baicalein|rutin"
flavonoid_subcategories <- c("Anthocyanins", "Flavonols / Flavones",
                             "Proanthocyanidins / Catechins", "Other flavonoids")

rx <- function(p) regex(p, ignore_case = TRUE)

add_category_cems <- function(df) {
  df %>%
    mutate(
      MetLabel = make_label_name(Metabolite),
      Category = case_when(
        str_detect(MetLabel, rx("^(ATP|ADP|AMP|NAD|NADP|NADH|NADPH|FAD|FMN)$")) ~ "Energy / Redox",
        str_detect(Metabolite, rx("Gly\\b|Ser\\b|Glycine|Serine")) ~ "Amino acids / Peptides",
        str_detect(Metabolite, rx("Citrate|Isocitrate|Malate|Succinate|Fumarate|Aconitate|Pyruvate|Lactate|DHAP|PEP|3PG|2PG")) ~ "Central carbon metabolism (glycolysis + TCA)",
        str_detect(Metabolite, rx("Putrescine|Spermidine|Agmatine|Cadaverine|Diaminopropane")) ~ "Polyamines",
        str_detect(Metabolite, rx("GSH|GSSG|Ascorb|Dehydroasc|Threonate|Threonic")) ~ "AsA/GSH",
        str_detect(Metabolite, rx("Sinapate|Benzoate|Shikimate|Quinate|Phenyl|Tyr\\b|Trp\\b|Phe\\b")) ~ "Aromatic / Phenylpropanoid",
        str_detect(Metabolite, rx("G6P|F6P|S7P|Rib|X5P|T6P|Suc6|UDP|CDP-choline|Choline|GPC")) ~ "Carbohydrates / Sugar metabolism",
        TRUE ~ "Other"
      ),
      SubCategory = NA_character_
    )
}

add_category_gcms <- function(df) {
  df %>%
    mutate(
      MetLabel = make_label_name(Metabolite),
      Category = case_when(
        str_detect(Metabolite, rx("Glycine|Serine|Proline|GABA|Alanine|Valine|Leucine|Isoleucine|Phenylalanine|Tyrosine|Tryptophan|Aspartic|Glutamic|Asparagine|Glutamine|Ornithine|Citrulline|Lysine|Histidine|Methionine")) ~ "Amino acids / Peptides",
        str_detect(MetLabel, rx(gcms_photorespiration_regex)) ~ "Photorespiration",
        str_detect(MetLabel, rx(gcms_central_carbon_regex)) ~ "Central carbon metabolism (glycolysis + TCA)",
        str_detect(Metabolite, rx("Glucose|Fructose|Sucrose|Trehalose|Raffinose")) ~ "Carbohydrates / Sugar metabolism",
        str_detect(Metabolite, rx("Mannitol|Sorbitol|Galactitol|Inositol|Xylitol")) ~ "Polyols",
        str_detect(Metabolite, rx("Ascorb|Threonate|Threonic")) ~ "AsA/GSH",
        TRUE ~ "Other"
      ),
      SubCategory = NA_character_
    )
}

add_category_lcms <- function(df) {
  df %>%
    mutate(
      Met_clean = clean_lcms_name(Metabolite),
      MetLabel = Met_clean,
      Category = case_when(
        str_detect(Met_clean, rx(paste(flavonoid_regex, "coumarin|scopoletin|chlorogen|caffe|ferul|sinap|cinnam|phenylprop|stilbene", sep = "|"))) ~ "Aromatic / Phenylpropanoid",
        str_detect(Met_clean, rx("\\bPC\\b|\\bPE\\b|Lyso|DG\\b|TG\\b|Cer|sphingo|phosphatidyl|galactolipid")) ~ "Lipids / Membrane",
        str_detect(Met_clean, rx("abscis|\\bABA\\b|jasmon|\\bJA\\b|salicy|\\bSA\\b|\\bIAA\\b|indole-3-acetic|gibberell|\\bGA\\b|auxin")) ~ "Hormones / Signals",
        str_detect(Met_clean, rx("diadinoxanthin|fucoxanthin|carotenone|apo-.*caroten")) ~ "Antioxidants",
        str_detect(Met_clean, rx("ascorb|dehydroasc|threonate|threonic|glutath|selenodiglutathione|ophthalmate")) ~ "AsA/GSH",
        str_detect(Met_clean, rx("violaxanthin|antheraxanthin|zeaxanthin|neoxanthin|lutein|beta-?carotene")) ~ "Antioxidants",
        str_detect(Met_clean, rx("glutam|aspart|glycine|serine|proline|gaba|valine|leucine|isoleucine|phenylalanine|tyrosine|tryptophan|ornithine|citrulline|lysine|histidine|methionine")) ~ "Amino acids / Peptides",
        TRUE ~ "Other"
      ),
      SubCategory = case_when(
        str_detect(Met_clean, rx(proanthocyanidin_regex)) ~ "Proanthocyanidins / Catechins",
        str_detect(Met_clean, rx(anthocyanin_regex)) ~ "Anthocyanins",
        str_detect(Met_clean, rx(flavonol_flavone_regex)) ~ "Flavonols / Flavones",
        str_detect(Met_clean, rx(flavonoid_regex)) ~ "Other flavonoids",
        str_detect(Met_clean, rx("coumarin|scopoletin")) ~ "Coumarins",
        str_detect(Met_clean, rx("chlorogen|caffe|ferul|sinap|cinnam|phenylprop|stilbene")) ~ "Phenylpropanoid (others)",
        str_detect(Met_clean, rx("diadinoxanthin|fucoxanthin")) ~ "Carotenoids (xanthophylls)",
        str_detect(Met_clean, rx("carotenone|apo-.*caroten")) ~ "Apocarotenoids",
        TRUE ~ NA_character_
      )
    ) %>%
    select(-Met_clean)
}

annotated <- bind_rows(
  add_category_cems(filter(supp21, Platform == "CE-TOF MS")),
  add_category_gcms(filter(supp21, Platform == "GC-MS/MS")),
  add_category_lcms(filter(supp21, Platform == "LC-MS/MS"))
)

write_csv(
  annotated %>% select(Platform, Metabolite, MetLabel, Category, SubCategory),
  "output/metabolite_category_annotation.csv"
)

# ---- within-species HL > LL tests ---------------------------------------------
within_species <- log2_long %>%
  group_by(Platform, Metabolite, Lineage) %>%
  summarise(
    p_HL_gt_LL = {
      hl <- log2y[Light == "HL"]
      ll <- log2y[Light == "LL"]
      tryCatch(t.test(hl, ll, alternative = "greater")$p.value, error = function(e) NA_real_)
    },
    .groups = "drop"
  ) %>%
  group_by(Platform, Lineage) %>%
  mutate(FDR_HL_gt_LL = p.adjust(p_HL_gt_LL, method = "BH")) %>%
  ungroup() %>%
  pivot_wider(names_from = Lineage, values_from = c(p_HL_gt_LL, FDR_HL_gt_LL))

annotated <- annotated %>%
  left_join(within_species, by = c("Platform", "Metabolite")) %>%
  mutate(
    Sp_sig_up = log2FC_Sp > 0 & FDR_HL_gt_LL_Sp < 0.05,
    Wa_sig_up = log2FC_Wa > 0 & FDR_HL_gt_LL_Wa < 0.05,
    hl_response_class = case_when(
      Sp_sig_up & Wa_sig_up ~ "Shared significant HL-up",
      Sp_sig_up & !Wa_sig_up ~ "Sp-specific significant HL-up",
      !Sp_sig_up & Wa_sig_up ~ "Wa-specific significant HL-up",
      TRUE ~ "Not significant / other"
    ),
    plot_class = case_when(
      sig_class == "Sp > Wa" ~ "Sp > Wa",
      sig_class == "Wa > Sp" ~ "Wa > Sp",
      TRUE ~ "Not significant / other"
    )
  )

# ---- Fig. 6 panels ---------------------------------------------------------------
panel_of <- function(Platform, Category, SubCategory) {
  case_when(
    Platform == "LC-MS/MS" & Category == "Aromatic / Phenylpropanoid" &
      SubCategory %in% flavonoid_subcategories ~ "Fig. 6A: Flavonoids",
    Platform == "GC-MS/MS" & Category == "Amino acids / Peptides" ~ "Fig. 6B: Amino acids / peptides",
    Platform == "GC-MS/MS" & Category == "Central carbon metabolism (glycolysis + TCA)" ~ "Fig. 6C: Central carbon metabolism",
    Platform == "GC-MS/MS" & Category %in% c("Carbohydrates / Sugar metabolism", "Polyols") ~ "Fig. 6D: Carbohydrates / sugar metabolism",
    Platform == "CE-TOF MS" & Category == "Energy / Redox" ~ "Fig. 6E: Energy / redox",
    TRUE ~ NA_character_
  )
}

supp22 <- annotated %>%
  mutate(panel_assignment = panel_of(Platform, Category, SubCategory)) %>%
  filter(!is.na(panel_assignment)) %>%
  transmute(
    panel_assignment, MS, Metabolite, MetLabel, Category, SubCategory,
    Sp_HL_response = log2FC_Sp, Wa_HL_response = log2FC_Wa,
    diff_Wa_minus_Sp = log2FC_Wa - log2FC_Sp,
    p_interaction, FDR, sig_class,
    p_HL_gt_LL_Sp, FDR_HL_gt_LL_Sp, p_HL_gt_LL_Wa, FDR_HL_gt_LL_Wa,
    hl_response_class, plot_class
  ) %>%
  arrange(panel_assignment, Metabolite)

write_csv(supp22, "output/SupplementaryData22_category_metabolite_log2FC_statistics.csv")

# ---- values cited in the main text ------------------------------------------------
cat("Metabolites per Fig. 6 panel (Supplementary Data 22: 115 in total):\n")
print(count(supp22, panel_assignment))

# -> Results "Lineage-specific metabolite accumulation": 61 flavonoid-related metabolites;
#    18 Sp > Wa and 12 Wa > Sp (Fig. 6a); amino acids 14 Wa > Sp vs 1 Sp > Wa (Fig. 6b);
#    central carbon 6 / 0, sugars and carbohydrates 10 / 0, energy/redox 4 / 0 (Fig. 6c-e)
print(supp22 %>% count(panel_assignment, plot_class) %>%
        pivot_wider(names_from = plot_class, values_from = n, values_fill = 0))

sp_up <- supp22$Sp_HL_response > 0 & supp22$FDR_HL_gt_LL_Sp < 0.05
wa_up <- supp22$Wa_HL_response > 0 & supp22$FDR_HL_gt_LL_Wa < 0.05
flav <- supp22$panel_assignment == "Fig. 6A: Flavonoids"
primary <- !flav

# -> Results: six flavonoid-related metabolites significantly increased under HL in both species
cat("Flavonoids increased under HL in both species:", sum(flav & sp_up & wa_up, na.rm = TRUE), "\n")
print(supp22 %>% filter(flav & sp_up & wa_up) %>% select(MetLabel, SubCategory))

# -> Results: 16 of 18 Sp-biased and 11 of 12 Wa-biased flavonoids also increased within species
cat("Flavonoids Sp > Wa also increased within Sp:",
    sum(flav & supp22$plot_class == "Sp > Wa" & sp_up, na.rm = TRUE), "of",
    sum(flav & supp22$plot_class == "Sp > Wa"), "\n")
cat("Flavonoids Wa > Sp also increased within Wa:",
    sum(flav & supp22$plot_class == "Wa > Sp" & wa_up, na.rm = TRUE), "of",
    sum(flav & supp22$plot_class == "Wa > Sp"), "\n")

# -> Results: 30 of 34 Wa-biased primary metabolites increased within Wa;
#    the single Sp-biased primary metabolite increased within Sp (Fig. 6b-e)
cat("Panels b-e, Wa > Sp also increased within Wa:",
    sum(primary & supp22$plot_class == "Wa > Sp" & wa_up, na.rm = TRUE), "of",
    sum(primary & supp22$plot_class == "Wa > Sp"), "\n")
cat("Panels b-e, Sp > Wa also increased within Sp:",
    sum(primary & supp22$plot_class == "Sp > Wa" & sp_up, na.rm = TRUE), "of",
    sum(primary & supp22$plot_class == "Sp > Wa"), "\n")
