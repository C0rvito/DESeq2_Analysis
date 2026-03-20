# =========================================================================
# File: src/utils_drive.R
# Description: Módulo de integração com Google Drive e Sheets (Google One)
# =========================================================================

library(googledrive)
library(googlesheets4)
library(dplyr)

# 1. IDENTIFICAÇÃO DO APP
drive_auth_configure(path = ".secrets/client_tcc.json")
gs4_auth_configure(path = ".secrets/client_tcc.json")

cat("API configurada")

# 2. AUTENTICAÇÃO DO USUÁRIO
drive_auth()

gs4_auth(token = drive_token())

cat("Ativação concluida")
