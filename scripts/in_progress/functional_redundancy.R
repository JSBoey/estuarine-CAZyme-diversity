# Functional redundancy

# Main question: Is functional redundancy more prevalent in some conditions than
#                others?

# Environment ----
library(matrixStats)
library(tidyverse)
library(patchwork)
library(dendextend)
library(colorspace)

# Import data ----
load("results/transformed_counts.RData")
annotation <- read_tsv("results/curated_annotation_table.tsv")
metadata <- read_tsv("data/sample_metadata.txt")
dbcan <- read_tsv("https://bcb.unl.edu/dbCAN_sub/data/fam-substrate-mapping-08252022.tsv")

# Endo-lysis enzymes ----
# Write out dbcan mapping for those within the annotation
ann_dbcan <- unique(annotation$dbcan_label) %>% 
  str_split(";") %>% 
  unlist() %>% 
  str_trim() %>% 
  str_remove("_\\d+") %>% 
  unique() %>% 
  str_subset("(GH|PL)") %>% 
  .[order(.)]

ann_ec <- annotation %>% 
  filter(grepl("(GH|PL)", dbcan_label)) %>% 
  pull(cazyme_consensus_ec) %>% 
  unique() %>% 
  str_split(";") %>% 
  unlist() %>% 
  str_trim() %>% 
  unique()

# dbcan %>% 
#   filter(Family %in% ann_dbcan & EC_Number %in% ann_ec) %>% 
#   write_tsv(file = "manual_curation/endolytic_enzyme.tsv")

## Import curated endolytic CAZymes
endolytic <- readxl::read_xlsx("manual_curation/endolytic_enzyme_curated.xlsx") %>% 
  drop_na() %>% 
  select(Family, EC_Number, Endolytic) %>% 
  # Manually adding missing PL families
  bind_rows(
    tribble(
      ~Family, ~EC_Number, ~Endolytic,
      "PL24" , "4.2.2.-" , "Yes", # Ulvan lyase
      "PL25" , "4.2.2.-" , "Yes", # Ulvan lyase
      "PL28" , "4.2.2.-" , "Yes", # Ulvan lyase
      "PL35" , "4.2.2.-" , "Yes", # Chondroitin lyase
      "PL39" , "4.2.2.-" , "Yes", # Poly-MG-lyase (alginate)
      "PL39" , "4.2.2.11" , "Yes", # Poly-G-lyase (alginate)
      "PL39" , "4.2.2.3" , "Yes", # Poly-M-lyase (alginate)
      "PL40" , "4.2.2.-" , "Yes" # Ulvan lyase
    )
  )
  

# Collect annotations ----
# Node, CAZyme family, EC, and substrate information
sub_annotation_1 <- annotation %>% 
  select(node, bin, dbcan_label, contains("cazyme")) %>% 
  filter(grepl("(GH|PL)", dbcan_label)) %>%
  mutate(
    cazyme_consensus_ec = str_trim(str_remove(
      cazyme_consensus_ec, "3\\.2\\.1\\.201"
    )),
    cazyme_substrate_higher = str_replace_all(cazyme_substrate_higher, ", ", ";")
  ) %>%
  drop_na()

# Split nodes into CAZyme families, EC and substrates
sub_annotation_2 <- sub_annotation_1 %>% 
  mutate(
    across(
      contains(c("cazyme", "dbcan")), \(s) {
        str_split(s, ";") %>% 
          map(str_trim)
      }
    ),
    bin = str_remove(bin, "_pred"),
    family = map(dbcan_label, \(s) unique(str_remove_all(s, pattern = "_\\d+"))),
    endolytic = case_when(
      map_lgl(
        cazyme_consensus_ec, 
        ~ any(.x %in% unique(str_subset(endolytic$EC_Number, "-", T)))
      ) ~ "Endolytic",
      map2_lgl(
        cazyme_consensus_ec, family, \(ec, cz) {
          unchar_ec <- unique(str_subset(endolytic$EC_Number, "-"))
          unchar_cz <- filter(endolytic, grepl("-", EC_Number)) %>% 
            pull(Family) %>% 
            unique()
          any(ec %in% unchar_ec & cz %in% unchar_cz)
        }
      ) ~ "Endolytic"
    ))

## Enzyme-substrate relationship
## For non-ambiguous enzymes only
dbs <- dbcan %>% 
  mutate(
    Substrate_high_level = str_replace(
      Substrate_high_level, ", ", ";"
    )
  ) %>% 
  filter(EC_Number %in% ann_ec & Family %in% ann_dbcan) %>% 
  select(EC_Number, Substrate_high_level) %>% 
  filter(!grepl("-", EC_Number) & grepl("(3\\.2\\.1|4\\.2\\.2)", EC_Number)) %>% 
  group_by(EC_Number) %>% 
  summarise(
    Substrate_high_level = list(Substrate_high_level) %>% 
      map(\(s) unique(unlist(str_split(s, ";"))))
  )

# Sum/Binarise percentages ----
# Sum by per sample & per EC
sum_EC <- map(transform_count[grepl("TPM", names(transform_count))], \(m) {
  m <- as_tibble(as.data.frame(m), rownames = "node")
  d <- sub_annotation_2 %>% 
    select(node, bin, cazyme_consensus_ec) %>%
    unnest(cazyme_consensus_ec)
  md <- left_join(d, m, by = "node") %>% 
    group_by(bin, cazyme_consensus_ec) %>% 
    summarise(
      across(where(is.numeric), sum)
    ) %>% 
    filter(!grepl("-", cazyme_consensus_ec)) %>% 
    filter(grepl("3.2.1|4.2.2", cazyme_consensus_ec))
})

# Percentage data
percent_EC <- map(sum_EC, \(df) {
  mutate(df, across(where(is.numeric), \(x) x * 100 / sum(x)))
})

# Binary data
binary_EC <- map(sum_EC, \(df) {
  mutate(df, across(where(is.numeric), \(x) ifelse(x > 0, 1, 0)))
})

# MAG weighted binary data ----
wbinary_EC <- map(binary_EC, \(df) {
  group_by(df, cazyme_consensus_ec) %>% 
    summarise(across(where(is.numeric), sum))
})

# Convert to matrix
wbinary_EC_mt <- map(wbinary_EC, \(df) {
  mt <- as.matrix(select(df, where(is.numeric)))
  rownames(mt) <- df$cazyme_consensus_ec
  
  mt
})

# Plot MAG weighted binary data ----
## Prepare main data ----
plot_data <- map(wbinary_EC_mt, \(m) {
  # Split data by habitat type
  list(
    "water"    = m[, grepl("Fil", colnames(m))],
    "sediment" = m[, grepl("Sed", colnames(m))]
  )
}) %>% 
  # Flatten list
  list_flatten() %>% 
  # Change names
  set_names(., nm = str_remove(names(.), "_TPM"))

### Transform data
plot_data <- map(plot_data, \(m) log2(m + 1))

## Order data
order_data <- map(plot_data, \(m) {
  # Distance matrix
  d <- dist(m)
  # Hierarchical clustering
  h <- hclust(d, method = "ward.D2")
})

## Prepare annotation data ----
annotation_data <- tibble(
  "ec_number" = rownames(plot_data$wgs_sediment)
) %>% 
  left_join(dbs, by = c("ec_number" = "EC_Number")) %>% 
  rename_with(str_to_lower)

### Replace NULLs 
### 3.2.1.163 = 1,6-alpha-D-mannosidase (alpha-mannan)
### 3.2.1.197 = beta-1,2-mannosidase (beta-mannan)
### 3.2.1.216 = kojibiose hydrolase (alpha-glucan)

annotation_data <- annotation_data %>% 
  mutate(
    substrate_high_level = case_when(
      ec_number == "3.2.1.163" ~ list("alpha-mannan"),
      ec_number == "3.2.1.197" ~ list("beta-mannan"),
      ec_number == "3.2.1.216" ~ list("alpha-glucan"),
      TRUE ~ substrate_high_level
    )
  )

### Check substrates
unlist(annotation_data$substrate_high_level) %>% 
  unique() %>% 
  str_sort()

### Replace lower level substrate classifications with higher level ones

### cellulose -> beta-glucan
### {starch, glycogen, trehalose} -> alpha-glucan
### {chitosan, chitin} -> GlcNAc
### human milk polysaccharide -> host glycan
### {agarose, alginate, fucoidan, porphyran} -> algal glycan

annotation_data <- annotation_data %>% 
  transmute(
    ec_number,
    substrate = map(substrate_high_level, \(s) {
      a  <- c("starch", "glycogen", "trehalose")
      b  <- c("cellulose")
      al <- c("agarose", "alginate", "fucoidan", "porphyran")
      ch <- c("chitosan")
      hg <- c("human milk polysaccharide")
      
      case_when(
        s %in% a  ~ "alpha-glucan",
        s %in% b  ~ "beta-glucan",
        s %in% al ~ "algal glycan",
        s %in% ch ~ "chitin",
        s %in% hg ~ "host glycan",
        TRUE ~ s
      )
    }) %>% 
      map(\(s) str_to_sentence(str_sort(unique(s))))
  )

### Check substrates
unlist(annotation_data$substrate) %>% 
  unique() %>% 
  str_sort()

### Split annotation data
annotation_data <- annotation_data %>% 
  unnest_wider(substrate, names_sep = "_")

## Pivot data ----
### Main data
plot_data_long <- map2(plot_data, order_data, \(m, h) {
  o <- o <- order_data$wgs_water$labels[order_data$wgs_water$order]
  # Convert to tibble and extract row names
  as.data.frame(m) %>% 
    as_tibble(rownames = "ec_number") %>%
    # Pivot data
    pivot_longer(
      where(is.numeric), names_to = "sample", values_to = "nbin"
    )
})

### Annotations: row
annotation_data_long <- map(order_data, \(h) {
  o <- order_data$wgs_water$labels[order_data$wgs_water$order]
  # Pivot data
  annotation_data %>% 
    pivot_longer(-ec_number, names_to = "sub_h", values_to = "substrate")
})

## Plot code ----
## Rows are ordered according to hierarchical clustering result of the water 
##   column samples.
## Arrangement:
## Water + Sediment + row annotations
## col_annotations
## Legend
plot_list <- map(c("wgs" = "wgs", "wts" = "wts"), \(x) {
  
  # Subset data
  dat <- keep_at(plot_data_long, \(nm) grepl(x, nm))
  
  # Subset annotation
  ann <- keep_at(annotation_data_long, \(nm) grepl(x, nm))
  
  # Subset metadata
  env <- map(dat, \(df) {
    filter(metadata, sample %in% df$sample) %>% 
      select(sample, salinity)
  })
  
  # Subset hclust results from water samples
  h_nm <- paste0(x, "_water")
  h <- order_data[[h_nm]]
  o <- h$labels[h$order]
  
  # Factorise EC numbers
  f <- function(x) {
    mutate(x, ec_number = factor(ec_number, levels = o))
  }
  dat <- map(dat, f)
  ann <- map(ann, f)
  
  # Set global theme
  my_theme <- 
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text  = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank()
    )
  
  # Get maximums for unified scales
  f2 <- function(l, col) {
    map_dbl(l, \(df) max(df[[col]])) %>% 
      .[which.max(.)] %>% 
      round(.) + 1
  } 
  MAX_1 <- f2(dat, "nbin")
  MAX_2 <- f2(env, "salinity")
  
  # Rename plots for verbosity
  f3 <- function(nm) {
    str_remove(nm, "^.*_")
  }
  
  # Annotation (substrate) colours
  fill_1 <- c(
    "Algal glycan" = "#171738",
    "Alpha-galactan" = "#1E6BB3",
    "Alpha-glucan" = "#FFC20A",
    "Alpha-mannan" = "#EFD780",
    "Alpha-rhamnoside" = "#186318",
    "Arabinan" = "#D84797",
    "Arabinogalactan protein" = "#6A8EAE",
    "Beta-fucosides" = "#A0C1B9",
    "Beta-galactan" = "#475841",
    "Beta-glucan" = "#228B22",
    "Beta-glucuronan" = "#82A0C4",
    "Beta-mannan" = "#BEEF9E",
    "Chitin" = "#DBA159",
    "Exo-polysaccharide" = "#77867F",
    "Fructan" = "#FE5D26",
    "Glycolipid" = "#A02C18",
    "Host glycan" = "#355070",
    "Pectin" = "#FE9471",
    "Peptidoglycan" = "#473198",
    "Polyphenol" = "#6E0D25",
    "Raffinose" = "#8B687F",
    "Sialic acid" = "#9C89A9",
    "Sucrose" = "#7B435B",
    "Xylan" = "#D17B0F",
    "Xyloglucan" = "#EC9F05"
  )
  
  # Plot heatmaps
  p1 <- map(dat, \(x) {
    ggplot(x, aes(x = sample, y = ec_number)) +
      geom_tile(aes(fill = nbin)) +
      scale_fill_viridis_c(limits = c(0, MAX_1), option = 'rocket') +
      scale_x_discrete(expand = c(0, 0)) +
      scale_y_discrete(expand = c(0, 0)) +
      my_theme
  }) %>% 
    set_names(nm = f3)
  
  # Plot annotation 
  p2 <- filter(ann[[h_nm]],!is.na(substrate)) %>%
    ggplot(aes(x = sub_h, y = ec_number)) +
      geom_tile(aes(fill = substrate)) +
      # scale_fill_manual(values = fill_1) +
      scale_fill_viridis_d(option = "turbo") +
      scale_x_discrete(expand = c(0, 0)) +
      scale_y_discrete(expand = c(0, 0)) +
      my_theme
  
  # Plot salinity
  p3 <- map(env, \(x) {
    ggplot(x, aes(x = sample, y = 1)) +
      geom_tile(aes(fill = salinity)) +
      scale_fill_viridis_c(option = "D", limits = c(0, MAX_2)) +
      scale_x_discrete(expand = c(0, 0)) +
      scale_y_discrete(expand = c(0, 0)) +
      my_theme
  }) %>% 
    set_names(nm = f3)
  
  # Patch plots
  p1$water + p1$sediment + p2 + p3$water + p3$sediment + plot_spacer() +
    plot_layout(
      nrow = 2, heights = c(114, 4),
      ncol = 3, widths  = c(9, 27, 4),
      guides = 'collect'
    ) &
    theme(
      legend.position = 'bottom'
    )
  
})

plot_list

pdf("results/Figures/MAGs_with_ECs.heatmap.pdf", paper = "a4r")
plot_list
dev.off()

# Latent vs realised functional redundancy ----
bin_ratio <- c("water" = "water", "sediment" = "sediment") %>% 
  set_names(., nm = .) %>% 
  map(\(x) {
    # Get relevant matrices
    l <- keep_at(plot_data, \(nm) grepl(x, nm))
    # Subset columns to retain same columns
    col_nm <- map(l, colnames) %>% 
      reduce(dplyr::intersect)
    l <- map(l, \(m) m[, col_nm]) %>% 
      set_names(nm = \(s) str_remove(s, "_.*"))
    # Division
    r <- l[["wts"]] / l[["wgs"]]
    r
  })

bin_ratio_long <- map(bin_ratio, \(m) {
  df <- as_tibble(as.data.frame(m), rownames = "ec_number")
  pivot_longer(df, where(is.numeric), names_to = "sample", values_to = "bin_ratio")
}) %>% 
  bind_rows(.id = "type")

ggplot(bin_ratio_long, aes(x = sample, y = bin_ratio)) +
  geom_boxplot() +
  facet_wrap(~type, scales = 'free')


# Misc code ----
color_choices <- expand.grid(
  seq(0, 288, length.out = 5), sample(seq(0.9, 0.3, length.out = 5), 5, replace = T)
)
color_choices[, 3] <- rep(sample(seq(1, 0.6, length.out = 5), 5, replace = T), each = 5)
hsv_choices <- HSV(as.matrix(color_choices))
hex_choices <- hex(hsv_choices, fixup = TRUE)



for (r in 1:nrow(color_choices)) {
  print(color_choices[r, ])
}
