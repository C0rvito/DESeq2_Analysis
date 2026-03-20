# =========================================================================
# File: src/main.R
# Description: Script orquestrador principal. Inicializa o ambiente,
# define a lista de alvos e controla a execução paralela/sequencial do pipeline.
# =========================================================================

# 1. CRIAÇÃO DE DIRETÓRIOS (Executado apenas se não existirem)
dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("results", recursive = TRUE, showWarnings = FALSE)

# 2. CARREGAR A FUNÇÃO DE ANÁLISE DEseq PARA DUPLICAÇÃO DE GENOMA
source("src/wdg_pipeline.R")

# 3. DEFINIR A LISTA DE TUMORES ALVO
# Adicione ou remova as siglas do TCGA conforme a necessidade do projeto.
cancer_targets <- c("LUAD") 

# 4. LOOP DE EXECUÇÃO COM TRATAMENTO DE ERROS
cat("Iniciando orquestração do pipeline para", length(cancer_targets), "tipos de câncer.\n")

for (cancer in cancer_targets) {
  
  # O bloco tryCatch isola falhas. Se algum câncer buscado falhar a função só passará para o outro tipo de câncer.
  tryCatch({
    
    run_wgd_analysis(target_cancer = cancer, 
                     metadata_path = "data/raw/Samples_meta_WGD.xlsx")
    
  }, error = function(e) {
    cat("\n[FALHA] Erro ao processar o alvo", cancer, ":\n")
    cat("Motivo:", conditionMessage(e), "\n\n")
  })
}

cat("\nOrquestração finalizada.\n")