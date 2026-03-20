# =======================================================================================
# File: src/wgd_pipeline.R
# Description: Módulo de processamento moldado a partir do pipeline RUN_DE_v2.r do Mateus
# ======================================================================================

# Carregamento de dependências
suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(DESeq2)
  library(dplyr)
  library(readr)
  library(tibble)
  library(readxl)
  library(apeglm)
})

run_wgd_analysis <- function(target_cancer, metadata_path = "data/raw/Samples_meta_WGD.xlsx") {
  
  cat("\n====================================================\n")
  cat("Iniciando processamento para o alvo:", target_cancer, "\n")
  cat("====================================================\n")
  
  # -------------------------------------------------------------------------
  # 1. IMPORTAÇÃO E FILTRAGEM DE METADADOS
  # -------------------------------------------------------------------------
  if (!file.exists(metadata_path)) {
    stop("Arquivo de metadados não encontrado em: ", metadata_path)
  }
  
  metadata <- read_excel(metadata_path)
  
  meta_clean <- metadata %>%
    dplyr::select(Sample, Type, Genome_doublings) %>%
    dplyr::filter(!is.na(Genome_doublings), !is.na(Type)) %>%
    dplyr::filter(Type == target_cancer) %>% 
    dplyr::mutate(WGD_status = ifelse(Genome_doublings > 0, "WGD_1_2", "WGD_0")) %>%
    dplyr::mutate(WGD_status = factor(WGD_status, levels = c("WGD_0", "WGD_1_2"))) %>%
    dplyr::mutate(Sample = substr(Sample, 1, 15)) # Padroniza barcode para 15 caracteres
  
  # Exportação para uso na metodologia
  write_tsv(meta_clean, paste0("data/processed/metadata_filtrada_", target_cancer, ".tsv"))
  
  # Verificação de segurança: n mínimo por grupo
  counts_per_group <- table(meta_clean$WGD_status)
  if (is.na(counts_per_group["WGD_0"]) || counts_per_group["WGD_0"] < 3 ||
      is.na(counts_per_group["WGD_1_2"]) || counts_per_group["WGD_1_2"] < 3) {
    stop(paste("Amostras insuficientes em um ou ambos os grupos para", target_cancer, 
               ". Contagens:", paste(names(counts_per_group), counts_per_group, sep="=", collapse=", ")))
  }
  
  # -------------------------------------------------------------------------
  # 2. DOWNLOAD & PREPARAÇÃO DOS DADOS BRUTOS (TCGA)
  # -------------------------------------------------------------------------
  tcga_project <- paste0("TCGA-", target_cancer)
  cat("Consultando GDC para contagens STAR do projeto", tcga_project, "...\n")
  
  query <- GDCquery(
    project = tcga_project,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts",
    experimental.strategy = "RNA-Seq"
  )
  
  # Reforço de rede: timeout extremo e motor libcurl
  options(timeout = 3600)
  options(download.file.method = "libcurl")
  
  # =========================================================================
  # O MÉTODO "VASSOURA" V2.0
  # =========================================================================
  root_dir <- getwd() 
  raw_dir <- file.path(root_dir, "data", "raw")
  gdc_target <- file.path(raw_dir, "GDCdata")
  
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(gdc_target, recursive = TRUE, showWarnings = FALSE)
  
  setwd(raw_dir) 
  
  tryCatch({
    # Fatiamento 
    cat("Baixando arquivos temporários na pasta raw...\n")
    GDCdownload(query, method = "api", files.per.chunk = 1)
    
    cat("Consolidando matrizes de contagem...\n")
    se <- GDCprepare(query)
    
    raw_counts_mat <- assay(se, "unstranded")
    gene_annotations <- as.data.frame(rowData(se)) %>%
      dplyr::select(gene_id, gene_name, gene_type)
    
  }, finally = {
    cat("Executando limpeza de diretórios (Método Vassoura V2.0)...\n")
    
    all_files <- list.files(raw_dir, full.names = TRUE)
    
    # RegEx expandido: Pega UUIDs, tar.gz, MANIFEST e logs
    tcga_junk <- all_files[
      grepl("^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$", basename(all_files)) | 
        grepl("\\.tar\\.gz$", basename(all_files)) |
        grepl("MANIFEST\\.txt$", basename(all_files)) |
        grepl("logs", basename(all_files))
    ]
    
    if(length(tcga_junk) > 0) {
      file.rename(from = tcga_junk, 
                  to = file.path(gdc_target, basename(tcga_junk)))
      cat("  ->", length(tcga_junk), "arquivos/pastas do TCGA movidos para GDCdata.\n")
    }
    
    setwd(root_dir) 
  })
  # =========================================================================
  
  colnames(raw_counts_mat) <- substr(colnames(raw_counts_mat), 1, 15)
  unique_samples <- !duplicated(colnames(raw_counts_mat))
  raw_counts_mat <- raw_counts_mat[, unique_samples]
  
  # -------------------------------------------------------------------------
  # 3. ALINHAMENTO ENTRE METADADOS E CONTAGENS
  # -------------------------------------------------------------------------
  shared_samples <- intersect(colnames(raw_counts_mat), meta_clean$Sample)
  if (length(shared_samples) < 6) {
    stop("Quantidade insuficiente de amostras correspondentes entre os dados do GDC e os metadados.")
  }
  
  meta_aligned <- meta_clean %>%
    dplyr::filter(Sample %in% shared_samples) %>%
    dplyr::arrange(match(Sample, shared_samples))
  
  counts_aligned <- raw_counts_mat[, meta_aligned$Sample]
  
  # Exportação das matrizes alinhadas para a metodologia
  write_tsv(meta_aligned, paste0("data/processed/meta_aligned_", target_cancer, ".tsv"))
  
  counts_export <- as.data.frame(counts_aligned) %>% tibble::rownames_to_column("gene_id")
  write_tsv(counts_export, paste0("data/processed/counts_aligned_", target_cancer, ".tsv"))
  
  # -------------------------------------------------------------------------
  # 4. EXECUÇÃO DO DESeq2 E SHRINKAGE
  # -------------------------------------------------------------------------
  cat("Executando o modelo DESeq2...\n")
  
  dds <- DESeqDataSetFromMatrix(
    countData = counts_aligned,
    colData   = data.frame(row.names = meta_aligned$Sample, WGD_status = meta_aligned$WGD_status),
    design    = ~ WGD_status 
  )
  
  # Pré-filtragem dinâmica
  min_samples <- round(0.05 * ncol(dds))
  keep <- rowSums(counts(dds) >= 10) >= min_samples
  dds <- dds[keep, ]
  
  dds <- DESeq(dds, quiet = TRUE)
  
  cat("Aplicando LFC Shrinkage (apeglm)...\n")
  coef_name <- "WGD_status_WGD_1_2_vs_WGD_0" 
  res_shrunken <- lfcShrink(dds, coef = coef_name, type = "apeglm", quiet = TRUE)
  
  res_df <- as.data.frame(res_shrunken) %>%
    tibble::rownames_to_column("gene_id") %>%
    dplyr::left_join(gene_annotations, by = "gene_id") %>%
    dplyr::relocate(gene_name, gene_type, .after = gene_id) %>%
    dplyr::arrange(padj)
  
  cat("Extraindo contagens normalizadas...\n")
  norm_counts_df <- as.data.frame(counts(dds, normalized = TRUE)) %>%
    tibble::rownames_to_column("gene_id") %>%
    dplyr::left_join(gene_annotations, by = "gene_id") %>%
    dplyr::relocate(gene_name, gene_type, .after = gene_id)
  
  # -------------------------------------------------------------------------
  # 5. EXPORTAÇÃO DE RESULTADOS
  # -------------------------------------------------------------------------
  output_file <- paste0("results/DE_WGD_1_2_vs_0_", target_cancer, ".tsv")
  norm_output_file <- paste0("results/NormCounts_WGD_1_2_vs_0_", target_cancer, ".tsv")
  
  write_tsv(res_df, output_file)
  write_tsv(norm_counts_df, norm_output_file)
  
  cat("Sucesso! Resultados exportados para a pasta results/.\n")
  
  # Limpeza forçada de memória para a próxima iteração
  rm(query, se, raw_counts_mat, gene_annotations, counts_aligned, dds, res_shrunken)
  gc()
  
  return(TRUE)
}