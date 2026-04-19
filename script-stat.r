# === # 1. INSTALAÇÃO DE PACOTES (Execute apenas na primeira vez) # === #

{install.packages("devtools")
  install.packages("BiocManager")
  install.packages("remotes")
  install.packages("tidyverse")
  install.packages("dplyr")
  install.packages("RColorBrewer")
  install.packages("patchwork")
  install.packages("ggplot2")
  install.packages("vegan")
  install.packages("permute")
  install.packages("tibble")
  install.packages("ggrepel")
  install.packages("forcats")
  install.packages(c("igraph", "ggraph", "Hmisc"))
  install.packages("ggVennDiagram")
  #
  BiocManager::install("phyloseq")
  BiocManager::install("ape")
  BiocManager::install("DESeq2")
  BiocManager::install("microbiome")
}

# === # 2. CARREGAR BIBLIOTECAS # === #

{library(phyloseq)
  library(tidyverse)
  library(ape)
  library(vegan)
  library(RColorBrewer)
  library(patchwork)
  library(ggplot2)
  library(permute)
  library(dplyr)
  library(tibble)
  library(ggrepel)
  library(forcats) 
  library(DESeq2)
  library(igraph)
  library(ggraph)
  library(Hmisc)}

# ==============================================================================
# === # 3. CARREGAMENTO DE DADOS # === #
# ==============================================================================

ARQUIVO_ABUNDANCIA <- "abundance.tsv"  # Nome do arquivo da tabela de abundância (TSV)
ARQUIVO_METADADOS  <- "metadados.tsv"         # Nome do arquivo de metadados (TSV)

# ==============================================================================
# === # 4. e 5. IMPORTAÇÃO, LIMPEZA E CRIAÇÃO DO PHYLOSEQ OBJECT # === #
# ==============================================================================

{
  print("--- Iniciando Importação ---")
  
  # 4.1. Carrega as tabelas usando read.delim para arquivos TSV (separados por tab)
  # check.names = FALSE impede que o R altere o nome das amostras
  tabela_bruta <- read.delim(ARQUIVO_ABUNDANCIA, header = TRUE, sep = "\t", check.names = FALSE) 
  
  # Lendo os metadados e já definindo a primeira coluna (SampleID) como o nome das linhas
  metadados_df <- read.delim(ARQUIVO_METADADOS, header = TRUE, sep = "\t", row.names = 1)
  
  # 4.2. Limpeza da Tabela
  # Remove a coluna "total" do final da tabela (se ela existir)
  tabela_limpa <- tabela_bruta %>% select(-total)
  
  # Separa a coluna "tax" nos 7 níveis taxonômicos (do Reino à Espécie) usando o ";"
  tax_separada <- tabela_limpa %>%
    select(tax) %>% 
    separate(tax, 
             into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), 
             sep = ";", 
             fill = "right",
             extra = "drop") %>% # "drop" ignora se por acaso houver algum nível extra
    as.matrix()
  
  # Limpar os espaços em branco que às vezes vêm colados nos nomes taxonômicos
  tax_separada <- trimws(tax_separada)
  
  # 4.3. Criação da Matriz Numérica de OTUs/ASVs
  # Removemos a coluna 'tax' para deixar apenas os números
  otu_matriz <- tabela_limpa %>%
    select(-tax) %>%
    as.matrix()
  
  # Força a matriz a ser numérica
  class(otu_matriz) <- "numeric"
  
  # Criar IDs padronizados (ASV_1, ASV_2...) para as linhas
  ids_otas <- paste0("ASV_", 1:nrow(otu_matriz))
  rownames(otu_matriz) <- ids_otas
  rownames(tax_separada) <- ids_otas
  
  # === # MONTAGEM DO OBJETO PHYLOSEQ # === #
  
  # Transforma as matrizes nos componentes exigidos pelo Phyloseq
  OTU <- otu_table(otu_matriz, taxa_are_rows = TRUE)
  TAX <- tax_table(tax_separada)
  SAMPLE <- sample_data(metadados_df)
  
  # Junta tudo no objeto final
  ps <- phyloseq(OTU, TAX, SAMPLE)
  
  print("Objeto Phyloseq criado com sucesso:")
  print(ps)
}

# ==============================================================================
# === # PAIRWISE COMPARISONS (WILCOXON) FOR UNBALANCED DESIGN # === #
# ==============================================================================

print("--- Starting Specific Pairwise Comparisons ---")

# Filter exact subsets for comparison
# 1. Only animals on Silage diet (to compare Angus vs Nelore)
data_silage <- subset(alpha_data, Diet == "Silage")

# 2. Only Nelore breed (to compare Silage vs Whole grain)
data_nelore <- subset(alpha_data, Race == "Nelore")

metrics_tested <- c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "Pielou")

# Create an empty dataframe to store results
pairwise_results <- data.frame()

for (metric in metrics_tested) {
  
  # --- COMPARISON 1: Breed Effect (Fixed on Silage) ---
  formula_breed <- as.formula(paste(metric, "~ Race"))
  test_breed <- wilcox.test(formula_breed, data = data_silage, exact = FALSE)
  
  row_breed <- data.frame(
    Metric = metric,
    Comparison = "Breed Effect (Angus vs Nelore on Silage)",
    Exact_P_Value = test_breed$p.value
  )
  
  # --- COMPARISON 2: Diet Effect (Fixed on Nelore) ---
  formula_diet <- as.formula(paste(metric, "~ Diet"))
  test_diet <- wilcox.test(formula_diet, data = data_nelore, exact = FALSE)
  
  row_diet <- data.frame(
    Metric = metric,
    Comparison = "Diet Effect (Silage vs Whole grain on Nelore)",
    Exact_P_Value = test_diet$p.value
  )
  
  # Bind to final dataframe
  pairwise_results <- rbind(pairwise_results, row_breed, row_diet)
}

# Format Exact P-Value to decimals (without scientific notation)
pairwise_results$Formatted_P_Value <- ifelse(
  pairwise_results$Exact_P_Value < 0.001,
  "< 0.001",
  sprintf("%.4f", pairwise_results$Exact_P_Value)
)

# Add significance column (p < 0.05)
pairwise_results$Significant <- ifelse(pairwise_results$Exact_P_Value < 0.05, "Yes", "No")

# Export clean table
write.table(pairwise_results, 
            file = "Table_Alpha_Pairwise_Wilcoxon.tsv", 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE, 
            dec = ".")

print("--- Pairwise Comparisons Completed! ---")
print("Table saved as 'Table_Alpha_Pairwise_Wilcoxon.tsv'")

# ==============================================================================
# === # 6.9. BOXPLOTS FOR PAIRWISE COMPARISONS (WILCOXON) # === #
# ==============================================================================

print("--- Generating Boxplots for Pairwise Comparisons ---")

# Custom function to create boxplots using the requested visual style
plot_wilcoxon_subset <- function(df, metric_col, title_text, var_x) {
  
  # Wilcoxon Test
  formula_stats <- as.formula(paste(metric_col, "~", var_x))
  test <- wilcox.test(formula_stats, data = df, exact = FALSE)
  p_val <- test$p.value
  
  # Format P-Value
  p_text <- ifelse(
    p_val < 0.001,
    "p < 0.001",
    paste("p =", sprintf("%.4f", p_val))
  )
  full_title <- paste0(title_text, " (Wilcoxon: ", p_text, ")")
  
  # Main Plot (Gray background, white grids, red jitter points)
  p <- ggplot(df, aes(x = .data[[var_x]], y = .data[[metric_col]], fill = .data[[var_x]])) +
    geom_boxplot(color = "black", outlier.shape = NA, alpha = 0.8) +
    # Jitter identical to the attached script
    geom_jitter(width = 0.2, size = 2.5, shape = 21, fill = "red", color = "black", stroke = 0.3, alpha = 0.8) +
    facet_wrap(~ Sample_Type) + 
    labs(x = NULL, y = NULL) +
    theme(
      text = element_text(family = "serif"),
      axis.text.x = element_text(size = 14, color = "black", face = "bold"),
      axis.text.y = element_text(size = 12, color = "black"),
      legend.position = "none",
      plot.margin = unit(c(0, 0.2, 0.2, 0.2), "cm"),
      panel.background = element_rect(fill = "gray90", color = "black"),
      panel.grid.major = element_line(color = "white", linewidth = 0.5),
      panel.grid.minor = element_line(color = "white", linewidth = 0.25),
      strip.background = element_rect(fill = "gray30", color = "black"), # Dark banner for Solid/Liquid
      strip.text = element_text(face = "bold", size = 12, color = "white")
    )
  
  # Top Title Banner (gray70) - identical to attached script
  lbl <- ggplot() + 
    annotate("text", x = 1, y = 1, label = full_title, size = 5.5, fontface = "bold", family = "serif") +
    theme_void() +
    theme(panel.background = element_rect(fill = "gray70", color = "black"))
  
  # Combine banner and plot
  return(lbl / p + plot_layout(heights = c(0.15, 1)))
}

# --- Function to create a Main Banner for the whole panel ---
create_main_banner <- function(title_text) {
  ggplot() +
    annotate("text", x = 1, y = 1, label = title_text, 
             size = 7, fontface = "bold", family = "serif") +
    theme_void() +
    theme(
      panel.background = element_rect(fill = "gray70", color = "black"),
      plot.margin = unit(c(0.2, 0, 0.2, 0), "cm")
    )
}

# ------------------------------------------------------------------------------
# === PANEL 1: BREED EFFECT (Only Animals on Silage Diet) ===
# ------------------------------------------------------------------------------

# Generate the 6 plots (X-axis = Race)
p_obs_r  <- plot_wilcoxon_subset(data_silage, "Observed", "Observed Richness", "Race") + scale_fill_manual(values = c("Angus" = "#D55E00", "Nelore" = "#0072B2"))
p_chao_r <- plot_wilcoxon_subset(data_silage, "Chao1", "Chao1", "Race") + scale_fill_manual(values = c("Angus" = "#D55E00", "Nelore" = "#0072B2"))
p_ace_r  <- plot_wilcoxon_subset(data_silage, "ACE", "ACE", "Race") + scale_fill_manual(values = c("Angus" = "#D55E00", "Nelore" = "#0072B2"))
p_sha_r  <- plot_wilcoxon_subset(data_silage, "Shannon", "Shannon Index", "Race") + scale_fill_manual(values = c("Angus" = "#D55E00", "Nelore" = "#0072B2"))
p_sim_r  <- plot_wilcoxon_subset(data_silage, "Simpson", "Simpson Index", "Race") + scale_fill_manual(values = c("Angus" = "#D55E00", "Nelore" = "#0072B2"))
p_pie_r  <- plot_wilcoxon_subset(data_silage, "Pielou", "Pielou's Evenness (J)", "Race") + scale_fill_manual(values = c("Angus" = "#D55E00", "Nelore" = "#0072B2"))

# Combine Breed panel with the top banner
banner_breed <- create_main_banner("Effect of Breed on Alpha Diversity (Silage Diet Only)")
panel_breed <- banner_breed / ((p_obs_r | p_chao_r | p_ace_r) / (p_sha_r | p_sim_r | p_pie_r)) + plot_layout(heights = c(0.1, 1))

print(panel_breed)

ggsave("Alpha_Diversity_Breed_Effect.tiff", plot = panel_breed, device = "tiff", 
       width = 18, height = 13, units = "in", dpi = 600, compression = "lzw", bg = "white")

# ------------------------------------------------------------------------------
# === PANEL 2: DIET EFFECT (Only Nelore Breed Animals) ===
# ------------------------------------------------------------------------------

# Define colors for Diets
colors_diet <- c("Silage" = "#009E73", "Whole grain" = "#E69F00")

# Generate the 6 plots (X-axis = Diet)
p_obs_d  <- plot_wilcoxon_subset(data_nelore, "Observed", "Observed Richness", "Diet") + scale_fill_manual(values = colors_diet)
p_chao_d <- plot_wilcoxon_subset(data_nelore, "Chao1", "Chao1", "Diet") + scale_fill_manual(values = colors_diet)
p_ace_d  <- plot_wilcoxon_subset(data_nelore, "ACE", "ACE", "Diet") + scale_fill_manual(values = colors_diet)
p_sha_d  <- plot_wilcoxon_subset(data_nelore, "Shannon", "Shannon Index", "Diet") + scale_fill_manual(values = colors_diet)
p_sim_d  <- plot_wilcoxon_subset(data_nelore, "Simpson", "Simpson Index", "Diet") + scale_fill_manual(values = colors_diet)
p_pie_d  <- plot_wilcoxon_subset(data_nelore, "Pielou", "Pielou's Evenness (J)", "Diet") + scale_fill_manual(values = colors_diet)

# Combine Diet panel with the top banner
banner_diet <- create_main_banner("Effect of Diet on Alpha Diversity (Nelore Breed Only)")
panel_diet <- banner_diet / ((p_obs_d | p_chao_d | p_ace_d) / (p_sha_d | p_sim_d | p_pie_d)) + plot_layout(heights = c(0.1, 1))

print(panel_diet)

ggsave("Alpha_Diversity_Diet_Effect.tiff", plot = panel_diet, device = "tiff", 
       width = 18, height = 13, units = "in", dpi = 600, compression = "lzw", bg = "white")

print("--- Images 'Alpha_Diversity_Breed_Effect.tiff' and 'Alpha_Diversity_Diet_Effect.tiff' saved successfully! ---")


# ==============================================================================
# === # 7. BETA DIVERSITY (PCoA & PERMANOVA) # === #
# ==============================================================================

print("--- Starting Beta Diversity Analysis ---")

# 7.1. Transform counts to Relative Abundance
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))

# 7.2. Subsets
ps_silage <- subset_samples(ps_rel, Diet == "Silage")
ps_nelore <- subset_samples(ps_rel, Race == "Nelore")

# 7.3. Function
plot_pcoa_custom <- function(physeq, var_test, title_text, color_mapping) {
  
  # Metadata
  meta <- as(sample_data(physeq), "data.frame")
  meta$Group <- as.factor(meta[[var_test]])
  sample_data(physeq) <- sample_data(meta)
  
  # Inicializa variáveis (🔥 evita erro)
  p_val <- NA
  r_sq <- NA
  bd_test <- NULL
  
  # Checagem de níveis
  if(length(unique(meta$Group)) < 2) {
    warning("Group has less than 2 levels.")
    
  } else {
    
    # Distância
    dist_bc <- phyloseq::distance(physeq, method = "bray")
    
    # PERMANOVA
    set.seed(123)
    perm <- vegan::adonis2(dist_bc ~ Group, data = meta, permutations = 999)
    perm_df <- as.data.frame(perm)
    
    # BETADISPER
    bd <- vegan::betadisper(dist_bc, meta$Group)
    bd_test <- anova(bd)
    
    print(bd_test)
    print(perm_df)
    
    # Extração robusta
    first_row <- perm_df[1, , drop = FALSE]
    
    if("Pr(>F)" %in% colnames(first_row)) {
      p_val <- first_row[["Pr(>F)"]]
    }
    
    if("R2" %in% colnames(first_row)) {
      r_sq <- first_row[["R2"]]
    }
  }
  
  # Texto
  p_text <- ifelse(!is.na(p_val) && p_val < 0.001, "p < 0.001", 
                   ifelse(!is.na(p_val), paste("p =", sprintf("%.4f", p_val)), "p = N/A"))
  
  r_text <- ifelse(!is.na(r_sq), paste("R² =", sprintf("%.3f", r_sq)), "R² = N/A")
  
  full_title <- paste0(title_text, "\n(PERMANOVA: ", p_text, " | ", r_text, ")")
  
  # PCoA
  dist_bc <- phyloseq::distance(physeq, method = "bray")
  ord <- ordinate(physeq, method = "PCoA", distance = dist_bc)
  
  eig_vals <- ord$values$Relative_eig * 100
  
  # Plot
  p <- plot_ordination(physeq, ord, color = "Group") +
    geom_point(size = 4, alpha = 0.9) +
    stat_ellipse(aes(group = Group), type = "t", linetype = 2, linewidth = 0.8) +
    facet_wrap(~ Sample_Type) + 
    scale_color_manual(values = color_mapping) +
    theme_bw() +
    labs(
      x = paste0("Axis 1 (", round(eig_vals[1], 1), "%)"),
      y = paste0("Axis 2 (", round(eig_vals[2], 1), "%)"),
      color = var_test
    )
  
  # Banner
  lbl <- ggplot() + 
    annotate("text", x = 1, y = 1, label = full_title, size = 5.5, 
             fontface = "bold", family = "serif", lineheight = 1.2) +
    theme_void() +
    theme(panel.background = element_rect(fill = "gray70", color = "black"),
          plot.margin = unit(c(0, 0, 0.1, 0), "cm"))
  
  # Stats (🔥 seguro contra erro)
  stats <- data.frame(
    Test = var_test,
    PERMANOVA_R2 = r_sq,
    PERMANOVA_p = p_val,
    BETADISPER_p = if(!is.null(bd_test)) bd_test$`Pr(>F)`[1] else NA
  )
  
  return(list(
    plot = (lbl / p + plot_layout(heights = c(0.12, 1))),
    stats = stats
  ))
}

# ------------------------------------------------------------------------------
# BREED
# ------------------------------------------------------------------------------
colors_breed <- c("Angus" = "#D55E00", "Nelore" = "#0072B2")

res_breed <- plot_pcoa_custom(ps_silage, "Race", 
                              "Beta Diversity: Breed Effect on Community Composition", 
                              colors_breed)

pcoa_breed <- res_breed$plot
stats_breed <- res_breed$stats

print(pcoa_breed)

ggsave("Beta_Diversity_Breed_PCoA.tiff", plot = pcoa_breed, device = "tiff", 
       width = 14, height = 8, units = "in", dpi = 600, compression = "lzw", bg = "white")

# ------------------------------------------------------------------------------
# DIET
# ------------------------------------------------------------------------------
colors_diet <- c("Silage" = "#009E73", "Whole grain" = "#E69F00")

res_diet <- plot_pcoa_custom(ps_nelore, "Diet", 
                             "Beta Diversity: Diet Effect on Community Composition", 
                             colors_diet)

pcoa_diet <- res_diet$plot
stats_diet <- res_diet$stats

print(pcoa_diet)

ggsave("Beta_Diversity_Diet_PCoA.tiff", plot = pcoa_diet, device = "tiff", 
       width = 14, height = 8, units = "in", dpi = 600, compression = "lzw", bg = "white")

# ------------------------------------------------------------------------------
# SAVE STATS
# ------------------------------------------------------------------------------
stats_all <- rbind(stats_breed, stats_diet)

write.table(stats_all,
            file = "Beta_Diversity_Stats.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

print("--- Beta Diversity Analysis Completed Successfully! ---")


# ==============================================================================
# === # 8. INDIVIDUAL SAMPLE ABUNDANCE (PHYLUM & GENUS) - WITH OTHERS # === #
# ==============================================================================

print("--- Generating Individual Sample Barplots with 'Others' group ---")

# 8.1. Preparação: Transformar para % e Ordenar Amostras
ps_perc <- transform_sample_counts(ps, function(x) 100 * x / sum(x))

# Extrair metadados para ordenar os níveis do SampleID (Angus primeiro, depois Nelore)
meta_ord <- as(sample_data(ps_perc), "data.frame") %>%
  arrange(Race, Diet) %>%
  mutate(SampleID = factor(rownames(.), levels = rownames(.)))

sample_data(ps_perc) <- sample_data(meta_ord)

# --- Função para Títulos (Banners) ---
create_taxa_banner <- function(title_text) {
  ggplot() +
    annotate("text", x = 1, y = 1, label = title_text, size = 6, fontface = "bold", family = "serif") +
    theme_void() +
    theme(panel.background = element_rect(fill = "gray70", color = "black"),
          plot.margin = unit(c(0.2, 0, 0.2, 0), "cm"))
}

# ------------------------------------------------------------------------------
# === 8.2. INDIVIDUAL PHYLUM PLOT (Top 10 + Others) ===
# ------------------------------------------------------------------------------

ps_phylum <- tax_glom(ps_perc, "Phylum")
df_phylum <- psmelt(ps_phylum)

# Identificar os 10 Filos mais abundantes globalmente
top10_phyla <- df_phylum %>%
  group_by(Phylum) %>%
  summarise(Sum = sum(Abundance)) %>%
  arrange(desc(Sum)) %>%
  head(10) %>%
  pull(Phylum)

# Agrupar o restante como "Others"
df_phylum$Phylum <- ifelse(df_phylum$Phylum %in% top10_phyla, as.character(df_phylum$Phylum), "Others")
df_phylum$Phylum <- factor(df_phylum$Phylum, levels = c(top10_phyla, "Others"))

p_phylum_ind <- ggplot(df_phylum, aes(x = SampleID, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", width = 0.95) +
  facet_grid(. ~ Race + Diet, scales = "free_x", space = "free_x", switch = "x") +
  # Cores: Paired + Cinza para "Others"
  scale_fill_manual(values = c(brewer.pal(min(length(top10_phyla), 12), "Paired"), "#D3D3D3")) +
  theme_bw() +
  labs(x = NULL, y = "Relative Abundance (%)") +
  theme(
    text = element_text(family = "serif"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8, color = "black"),
    axis.title.y = element_text(face = "bold", size = 12),
    strip.background = element_rect(fill = "gray30", color = "black"),
    strip.text = element_text(color = "white", face = "bold", size = 10),
    strip.placement = "outside",
    panel.spacing = unit(0.1, "lines"),
    panel.background = element_rect(fill = "gray90"),
    legend.position = "right"
  )

banner_ph <- create_taxa_banner("Relative abundance: Phylum Level")
panel_ph_ind <- banner_ph / p_phylum_ind + plot_layout(heights = c(0.08, 1))

ggsave("Relative_Abundance_Individual_Phylum.tiff", plot = panel_ph_ind, width = 18, height = 9, dpi = 600, compression = "lzw")

# ------------------------------------------------------------------------------
# === 8.3. INDIVIDUAL GENUS PLOT (Top 15 + Others) ===
# ------------------------------------------------------------------------------

ps_genus <- tax_glom(ps_perc, "Genus")
df_genus <- psmelt(ps_genus)

# Identificar os 15 Gêneros mais abundantes
top15_genus <- df_genus %>%
  group_by(Genus) %>%
  summarise(Sum = sum(Abundance)) %>%
  arrange(desc(Sum)) %>%
  head(15) %>%
  pull(Genus)

# Agrupar o restante como "Others"
df_genus$Genus <- ifelse(df_genus$Genus %in% top15_genus, as.character(df_genus$Genus), "Others")
df_genus$Genus <- factor(df_genus$Genus, levels = c(top15_genus, "Others"))

p_genus_ind <- ggplot(df_genus, aes(x = SampleID, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", width = 0.95) +
  facet_grid(. ~ Race + Diet, scales = "free_x", space = "free_x", switch = "x") +
  # Cores: Paired + Dark2 + Cinza para "Others"
  scale_fill_manual(values = c(brewer.pal(12, "Paired"), brewer.pal(3, "Dark2"), "#D3D3D3")) +
  theme_bw() +
  labs(x = NULL, y = "Relative Abundance (%)") +
  theme(
    text = element_text(family = "serif"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8, color = "black"),
    axis.title.y = element_text(face = "bold", size = 12),
    strip.background = element_rect(fill = "gray30", color = "black"),
    strip.text = element_text(color = "white", face = "bold", size = 10),
    strip.placement = "outside",
    panel.spacing = unit(0.1, "lines"),
    panel.background = element_rect(fill = "gray90"),
    legend.text = element_text(face = "italic")
  )

# Ajustar para que "Others" não fique em itálico na legenda (opcional, mas profissional)
p_genus_ind <- p_genus_ind + theme(legend.text = element_text(face = ifelse(levels(df_genus$Genus) == "Others", "plain", "italic")))

banner_ge <- create_taxa_banner("Relative abundance: Genus Level")
panel_ge_ind <- banner_ge / p_genus_ind + plot_layout(heights = c(0.08, 1))

ggsave("Relative_Abundance_Individual_Genus.tiff", plot = panel_ge_ind, width = 18, height = 10, dpi = 600, compression = "lzw")

print("--- Individual sample plots with 'Others' category saved successfully! ---")

# ==============================================================================
# === # 9. DIFFERENTIAL ABUNDANCE (DESeq2) - NELORE: SILAGE VS WHOLE GRAIN # === #
# ==============================================================================

print("--- Starting DESeq2 Analysis (Nelore: Silage vs Whole grain) ---")

# 9.1. Preparação: Filtrar Nelore e remover táxons com zero
ps_nelore_ds <- subset_samples(ps, Race == "Nelore")
ps_nelore_ds <- prune_taxa(taxa_sums(ps_nelore_ds) > 0, ps_nelore_ds)

# 9.2. Converter para objeto DESeq2
ds <- phyloseq_to_deseq2(ps_nelore_ds, ~ Diet)

# 9.3. Executar o teste
ds <- DESeq(ds, test = "Wald", fitType = "local")

# 9.4. Extrair resultados (Whole grain vs Silage)
res_ds <- results(ds, contrast = c("Diet", "Whole grain", "Silage"), cooksCutoff = FALSE)

# 9.5. Organizar e Filtrar (p-adj < 0.05 E |log2FC| >= 7)
res_df <- as.data.frame(res_ds)
res_df$ASV <- rownames(res_df)
tax_table_df <- as.data.frame(as.matrix(tax_table(ps_nelore_ds)))
tax_table_df$ASV <- rownames(tax_table_df)

res_final <- res_df %>%
  left_join(tax_table_df, by = "ASV") %>%
  filter(padj < 0.05) %>%
  filter(!is.na(Genus)) %>%
  # --- FILTRO DE LOG2 FOLD CHANGE >= x ---
  filter(abs(log2FoldChange) >= 7) %>% 
  group_by(Genus) %>%
  dplyr::summarize(log2FoldChange = mean(log2FoldChange), padj = min(padj)) %>%
  arrange(log2FoldChange)

# --- CÁLCULO DAS ESTATÍSTICAS PARA O BANNER ATUALIZADO ---
total_sig <- nrow(res_final)
up_grain  <- sum(res_final$log2FoldChange > 0)
up_silage <- sum(res_final$log2FoldChange < 0)

# Criando o texto de confiabilidade estatística
stats_text <- paste0("Total: ", total_sig, " genera (|Log2FC| >= 7, p-adj < 0.05)\n", 
                     "Enriched: ", up_grain, " in Whole grain, ", up_silage, " in Silage")
# 9.6. Gerar o Gráfico
res_final$Color <- ifelse(res_final$log2FoldChange > 0, "#E69F00", "#009E73") 

p_deseq <- ggplot(res_final, aes(x = reorder(Genus, log2FoldChange), y = log2FoldChange, fill = Color)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
  coord_flip() + 
  scale_fill_identity() +
  theme_bw() +
  labs(x = "Genus", y = "Log2 Fold Change (Whole grain / Silage)") +
  # Anotações dinâmicas
  annotate("text", x = 1, y = max(res_final$log2FoldChange) * 0.7, 
           label = "Enriched in Whole grain", fontface = "italic", family = "serif", color = "#E69F00") +
  annotate("text", x = nrow(res_final), y = min(res_final$log2FoldChange) * 0.7, 
           label = "Enriched in Silage", fontface = "italic", family = "serif", color = "#009E73") +
  theme(
    text = element_text(family = "serif"),
    axis.text.y = element_text(face = "italic", size = 10, color = "black"),
    axis.text.x = element_text(size = 11, color = "black"),
    panel.background = element_rect(fill = "gray90"),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_line(color = "white")
  )

# Banner do Título com Estatística e Teste
lbl_ds <- ggplot() + 
  annotate("text", x = 1, y = 1, 
           label = paste0("Differential Abundance: Nelore Diet Effect (Wald Test)\n", stats_text), 
           size = 5, fontface = "bold", family = "serif", lineheight = 1.2) +
  theme_void() + 
  theme(panel.background = element_rect(fill = "gray70", color = "black"),
        plot.margin = unit(c(0.2, 0, 0.2, 0), "cm"))

# Combinar e exibir
panel_deseq <- lbl_ds / p_deseq + plot_layout(heights = c(0.18, 1))
print(panel_deseq)

ggsave("DESeq2_Nelore_Diet_Filtered_LFC7.tiff", plot = panel_deseq, width = 12, height = 10, dpi = 600, compression = "lzw")

# Exportar Tabela
write.table(res_final, "Table_DESeq2_Nelore_Diet_LFC7.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

print("--- DESeq2 analysis with LFC >= 5 filter completed! ---")


# ==============================================================================
# === # 10. CORE MICROBIOME HEATMAP - NELORE (SILAGE vs WHOLE GRAIN) # === #
# ==============================================================================

# 10.1. Parâmetros de Definição do Core
CORE_RANK <- "Genus"         # Nível taxonômico
CORE_PREVALENCIA <- 0.8      # Presente em pelo menos 80% das amostras
CORE_ABUNDANCIA <- 0.0001    # Abundância mínima de 0.01%
VAR_AGRUPAMENTO <- "Diet"    # Variável para o facet e ordenação

print(paste("Analisando Core Microbiome do Nelore a nível de:", CORE_RANK))

# 10.2. Preparar dados: Filtrar Nelore e converter para Abundância Relativa
ps_nelore <- subset_samples(ps, Race == "Nelore")
ps_rel <- transform_sample_counts(ps_nelore, function(x) x / sum(x))

# 10.3. Filtrar Core (Lógica de Prevalência)
# Ao usar TRUE no final, o objeto 'ps_core' já nasce como um phyloseq filtrado
ps_core <- filter_taxa(ps_rel, function(x) sum(x > CORE_ABUNDANCIA) > (CORE_PREVALENCIA * nsamples(ps_rel)), TRUE)

print(paste("Número de táxons no Core do Nelore:", ntaxa(ps_core)))

# 10.4. Agrupamento Taxonômico e Preparação do Dataframe
# Se o número de táxons for > 0, seguimos para o glom
if(ntaxa(ps_core) > 0){
  ps_core_glom <- tax_glom(ps_core, taxrank = CORE_RANK)
  df_core <- psmelt(ps_core_glom)
  
  # Ordenar amostras por Dieta e criar fator para o eixo X
  df_core <- df_core %>% arrange(Diet)
  df_core$SampleID <- factor(df_core$Sample, levels = unique(df_core$Sample))
  
  # Tratamento de Nomes para o Eixo Y
  df_core$Taxon_Label <- as.character(df_core[[CORE_RANK]])
  df_core$Taxon_Label[is.na(df_core$Taxon_Label)] <- "Unclassified"
} else {
  stop("Nenhum táxon passou pelos critérios de Core. Tente diminuir a CORE_PREVALENCIA.")
}

# 10.5. Construção do Gráfico (Heatmap)
titulo_core <- paste0("Core Microbiome: Nelore Diet Effect (", CORE_RANK, ")")
stats_banner <- paste0("Taxa count: ", ntaxa(ps_core), " | Prev > ", CORE_PREVALENCIA*100, "% | Abund > ", CORE_ABUNDANCIA*100, "%")

P_HEATMAP <- ggplot(df_core, aes(x = SampleID, y = Taxon_Label, fill = log10(Abundance * 100 + 0.01))) + 
  geom_tile(color = "white", linewidth = 0.1) +
  # Escala Spectral (Azul = Pouco, Vermelho = Muito)
  scale_fill_distiller(palette = "Spectral", direction = -1, name = "Log10(%)") +
  # Divisão por Dieta na base
  facet_grid(~ Diet, scales = "free_x", space = "free_x", switch = "x") +
  labs(x = NULL, y = paste(CORE_RANK, "(Core)")) +
  theme(
    text = element_text(family = "serif"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8, color = "black"),
    axis.text.y = element_text(size = 10, face = "italic", color = "black"),
    legend.position = "right",
    # Estética das abas de Dieta
    strip.background = element_rect(fill = "gray30", color = "black"),
    strip.text = element_text(size = 10, face = "bold", color = "white"),
    strip.placement = "outside",
    panel.background = element_rect(fill = "gray90", color = NA),
    panel.grid = element_blank(),
    plot.margin = unit(c(0, 0.5, 0.5, 0.5), "cm")
  )

# Banner de Título com Estatísticas de Filtro
label_core <- ggplot() + 
  annotate("text", x = 1, y = 1.1, label = titulo_core, size = 6, fontface = "bold", family = "serif") +
  annotate("text", x = 1, y = 0.8, label = stats_banner, size = 4.5, fontface = "italic", family = "serif") +
  xlim(0.5, 1.5) + ylim(0.5, 1.5) + 
  theme_void() + 
  theme(panel.background = element_rect(fill = "gray70", color = NA))

# Combinar Banner e Heatmap
grafico_final_core <- label_core / P_HEATMAP + plot_layout(heights = c(0.15, 1))

# 10.6. Exibir e Salvar
print(grafico_final_core)

altura_dinamica <- max(6, length(unique(df_core$Taxon_Label)) * 0.35)
ggsave("Core_Heatmap_Nelore_Diet.tiff", plot = grafico_final_core, device = "tiff", 
       width = 14, height = altura_dinamica, units = "in", dpi = 600, compression = "lzw", bg = "white")

# 10.7. Exportar Tabela do Core
tabela_core_out <- df_core %>% 
  group_by(Taxon_Label, Phylum, Diet) %>% 
  summarise(Mean_Abundance_Pct = mean(Abundance * 100), .groups = "drop") %>% 
  arrange(Diet, desc(Mean_Abundance_Pct))

write.table(tabela_core_out, "Table_Core_Nelore_Heatmap.tsv", sep="\t", quote=FALSE, row.names=FALSE)

# ==============================================================================
# === # 6. ANÁLISE DE COBERTURA (GOOD'S COVERAGE) # === #
# ==============================================================================

print("--- Starting Good's Coverage Analysis ---")

# 6.1. Definições de Corte (Ajuste se necessário)
COV_CORTE <- 0.99  # Linha tracejada de cobertura ideal

# 6.2. Função de Cálculo
get_goods_coverage <- function(x) {
  n1 <- sum(x == 1)
  N <- sum(x)
  return(1 - (n1 / N))
}

# 6.3. Processamento de Dados
cobertura <- apply(otu_table(ps), 2, get_goods_coverage)
sample_data(ps)$Goods_Coverage <- cobertura

dados_cov <- data.frame(sample_data(ps)) %>% 
  tibble::rownames_to_column("SampleID") 

# 6.4. Gráfico de Cobertura
P_COV <- ggplot(dados_cov, aes(x = reorder(SampleID, Goods_Coverage), y = Goods_Coverage)) +
  geom_col(fill = "#4682B4", color = "black", width = 0.7) +
  geom_hline(yintercept = COV_CORTE, color = "red", linetype = "dashed", linewidth = 0.8) +
  coord_flip() +
  labs(x = NULL, y = "Good's Coverage Index") +
  theme_bw() +
  theme(
    text = element_text(family = "serif"),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.text.y = element_text(size = 9, color = "black"), 
    panel.background = element_rect(fill = "gray90"),
    panel.grid.major = element_line(color = "white", linewidth = 0.5),
    panel.grid.minor = element_line(color = "white", linewidth = 0.25)
  )

# Banner do Título
label_cov <- ggplot() + 
  annotate("text", x = 1, y = 1, label = "Sequencing Depth: Good's Coverage", size = 7,
           fontface = "bold", family = "serif", hjust = 0.5) +
  xlim(0, 2) + ylim(0.5, 1.5) +
  theme_void() +
  theme(panel.background = element_rect(fill = "gray70", color = NA))

grafico_final_cov <- label_cov / P_COV + plot_layout(heights = c(0.1, 1))

# Salvar
ggsave("QC_Goods_Coverage.tiff", plot = grafico_final_cov, device = "tiff",
       width = 12, height = 14, units = "in", dpi = 600, bg = "white", compression = "lzw")


# ==============================================================================
# === # 7. CURVA DE RAREFAÇÃO # === #
# ==============================================================================

print("--- Starting Rarefaction Curve Analysis ---")

# 7.1. Parâmetros de Rarefação
RARE_STEP <- 500  # Intervalo de leitura para a curva

ps_clean <- prune_samples(sample_sums(ps) > 0, ps)

if (taxa_are_rows(ps_clean)) {
  otu_mat <- t(as(otu_table(ps_clean), "matrix"))
} else {
  otu_mat <- as(otu_table(ps_clean), "matrix")
}
class(otu_mat) <- "matrix"

# 7.2. Cálculo da Curva (Vegan)
out_rare <- rarecurve(otu_mat, step = RARE_STEP, sample = min(rowSums(otu_mat)), label = FALSE)

# 7.3. Organização dos dados para ggplot
names(out_rare) <- rownames(otu_mat)
rare_df <- lapply(names(out_rare), function(x) {
  df <- data.frame(Reads = attr(out_rare[[x]], "Subsample"),
                   Richness = out_rare[[x]])
  df$SampleID <- x
  return(df)
}) %>% bind_rows()

# Adicionar metadados para colorir por Dieta ou Raça se desejar
meta_simple <- data.frame(sample_data(ps_clean)) %>% tibble::rownames_to_column("SampleID")
rare_df <- rare_df %>% left_join(meta_simple, by = "SampleID")

dados_labels <- rare_df %>%
  group_by(SampleID) %>%
  summarise(Max_Reads = max(Reads), Max_Richness = max(Richness))

# 7.4. Gráfico de Rarefação
P_RARE <- ggplot(rare_df, aes(x = Reads, y = Richness, group = SampleID, color = Diet)) +
  geom_line(linewidth = 0.6, alpha = 0.8) +
  geom_text(data = dados_labels, 
            aes(x = Max_Reads, y = Max_Richness, label = SampleID),
            color = "black", hjust = -0.1, size = 2.5, family = "serif", check_overlap = TRUE) + 
  geom_vline(xintercept = min(sample_sums(ps)), linetype = "dashed", color = "darkred", linewidth = 0.5) +
  scale_color_manual(values = c("Silage" = "#009E73", "Whole grain" = "#E69F00", "Angus" = "#D55E00")) + # Cores padronizadas
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.2))) +
  labs(x = "Sequencing Depth (Reads)", y = "Observed Richness (ASVs)") +
  theme_bw() +
  theme(
    text = element_text(family = "serif"),
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 11, face = "bold"),
    legend.position = "right",
    panel.background = element_rect(fill = "gray90"),
    panel.grid.major = element_line(color = "white", linewidth = 0.5),
    panel.grid.minor = element_line(color = "white", linewidth = 0.25)
  )

# Banner do Título
label_rare <- ggplot() + 
  annotate("text", x = 1, y = 1, label = "Rarefaction Curve: Species Saturation", size = 6,
           fontface = "bold", family = "serif", hjust = 0.5) +
  xlim(0, 2) + ylim(0.5, 1.5) +
  theme_void() +
  theme(panel.background = element_rect(fill = "gray70", color = NA))

grafico_final_rare <- label_rare / P_RARE + plot_layout(heights = c(0.1, 1))

# Salvar
ggsave("QC_Rarefaction_Curve.tiff", plot = grafico_final_rare, device = "tiff",
       width = 14, height = 8, units = "in", dpi = 600, bg = "white", compression = "lzw")

print("--- QC Analysis Completed Successfully! ---")

# ==============================================================================
# === # 12. EXPORTAÇÃO DE ABUNÂNCIAS RELATIVAS (CSV) # === #
# ==============================================================================

print("--- Exporting Relative Abundance Tables to CSV ---")

# 12.1. Preparação: Garantir que temos o objeto em porcentagem (0-100)
ps_perc <- transform_sample_counts(ps, function(x) 100 * x / sum(x))

# --- FUNÇÃO AUXILIAR PARA GERAR TABELA WIDE COM 'OTHERS' ---
export_taxa_csv <- function(physeq, rank, top_n, filename) {
  
  # Agrupar no nível desejado
  ps_glom <- tax_glom(physeq, rank)
  df <- psmelt(ps_glom)
  
  # Identificar os Top N
  top_taxa <- df %>%
    group_by(!!sym(rank)) %>%
    summarise(Sum = sum(Abundance)) %>%
    arrange(desc(Sum)) %>%
    head(top_n) %>%
    pull(!!sym(rank))
  
  # Agrupar como Others
  df$Taxa_Final <- ifelse(df[[rank]] %in% top_taxa, as.character(df[[rank]]), "Others")
  
  # Converter para formato largo (Amostras nas colunas)
  df_wide <- df %>%
    group_by(Sample, Taxa_Final) %>%
    summarise(Abundance = sum(Abundance), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = Sample, values_from = Abundance)
  
  # Salvar CSV
  write.csv(df_wide, filename, row.names = FALSE)
  print(paste("Arquivo salvo:", filename))
}

# ------------------------------------------------------------------------------
# 12.2. EXPORTANDO TABELAS DOS GRÁFICOS (Com 'Others')
# ------------------------------------------------------------------------------

# Filo (Top 10 + Others)
export_taxa_csv(ps_perc, "Phylum", 10, "Abundancia_Relativa_Filo_Top10.csv")

# Gênero (Top 15 + Others)
export_taxa_csv(ps_perc, "Genus", 15, "Abundancia_Relativa_Genero_Top15.csv")


# ------------------------------------------------------------------------------
# 12.3. EXPORTANDO TABELAS COMPLETAS (Sem filtros, todos os táxons)
# ------------------------------------------------------------------------------

# Função para exportar tudo sem o grupo "Others"
export_full_csv <- function(physeq, rank, filename) {
  ps_glom <- tax_glom(physeq, rank)
  
  # Extrair abundância
  otu <- as.data.frame(otu_table(ps_glom))
  if (taxa_are_rows(ps_glom)) {
    # Se táxons forem linhas, já está quase pronto
  } else {
    otu <- as.data.frame(t(otu))
  }
  
  # Extrair taxonomia
  tax <- as.data.frame(as.matrix(tax_table(ps_glom)))
  
  # Unir Taxonomia + Abundância
  full_tab <- cbind(tax, otu)
  
  write.csv(full_tab, filename, row.names = TRUE)
  print(paste("Tabela completa salva:", filename))
}

# Exportar tabelas completas para o material suplementar
export_full_csv(ps_perc, "Phylum", "Suplementar_Full_Phylum_Abundance.csv")
export_full_csv(ps_perc, "Genus", "Suplementar_Full_Genus_Abundance.csv")

print("--- All CSV files exported successfully! ---")







