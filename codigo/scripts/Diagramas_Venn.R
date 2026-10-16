library(dplyr)

script_arg <- grep("^--file=", commandArgs(), value = TRUE)

if (length(script_arg) > 0) {
  script_path <- sub("^--file=", "", script_arg[1])
  setwd(normalizePath(file.path(dirname(script_path), "..", "..")))
}

input_directory <- file.path("resultados", "DEGs", "originales")

if (!dir.exists(input_directory)) {
  stop("Ejecute este script desde la raíz del repositorio.")
}

output_directory <- file.path("resultados", "DEGs", "venn")

dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

# Extrae los nombres de los DEGs
leer_degs <- function(path_csv) {
  df <- read.csv(path_csv, row.names = 1)
  rownames(df)  # devuelve solo los IDs de genes
}

# ENFOQUE 1: Comparación directa por estación
# Cada archivo proviene de un contraste regional DESeq2:
# Aysén_estación vs Magallanes_estación

degs_win <- leer_degs(file.path(input_directory, "DEGs_Aysen_vs_Magallanes_WIN.csv"))
degs_spr <- leer_degs(file.path(input_directory, "DEGs_Aysen_vs_Magallanes_SPR.csv"))
degs_sum <- leer_degs(file.path(input_directory, "DEGs_Aysen_vs_Magallanes_SUM.csv"))

lista_enfoque1 <- list(
  Invierno  = degs_win,
  Primavera = degs_spr,
  Verano    = degs_sum
)

# DEGs completos por estación
write.csv(degs_win, file.path(output_directory, "DEGs_completos_WIN.csv"), row.names = FALSE)
write.csv(degs_spr, file.path(output_directory, "DEGs_completos_SPR.csv"), row.names = FALSE)
write.csv(degs_sum, file.path(output_directory, "DEGs_completos_SUM.csv"), row.names = FALSE)

# DEGs que se conservan en las tres estaciones (núcleo de divergencia regional)
nucleo_estable <- Reduce(intersect, lista_enfoque1)
write.csv(nucleo_estable, file.path(output_directory, "nucleo_divergencia_regional.csv"), row.names = FALSE)

# Exclusivos de cada estación
exclusivos_win <- setdiff(degs_win, union(degs_spr, degs_sum))
exclusivos_spr <- setdiff(degs_spr, union(degs_win, degs_sum))
exclusivos_sum <- setdiff(degs_sum, union(degs_win, degs_spr))

write.csv(exclusivos_win, file.path(output_directory, "DEGs_exclusivos_WIN.csv"), row.names = FALSE)
write.csv(exclusivos_spr, file.path(output_directory, "DEGs_exclusivos_SPR.csv"), row.names = FALSE)
write.csv(exclusivos_sum, file.path(output_directory, "DEGs_exclusivos_SUM.csv"), row.names = FALSE)

# ENFOQUE 2: Firma transcripcional por región
# Cada archivo proviene de un contraste estacional interno DESeq2:
# region_estacionA vs region_estacionB

# Contrastes internos de Aysén
aysen_win_spr <- leer_degs(file.path(input_directory, "DEGs_Aysen_WIN_vs_SPR.csv"))
aysen_win_sum <- leer_degs(file.path(input_directory, "DEGs_Aysen_WIN_vs_SUM.csv"))
aysen_spr_sum <- leer_degs(file.path(input_directory, "DEGs_Aysen_SPR_vs_SUM.csv"))

# DEGs completos por contraste interno de Aysén
write.csv(aysen_win_spr, file.path(output_directory, "DEGs_completos_Aysen_WIN_vs_SPR.csv"), row.names = FALSE)
write.csv(aysen_win_sum, file.path(output_directory, "DEGs_completos_Aysen_WIN_vs_SUM.csv"), row.names = FALSE)
write.csv(aysen_spr_sum, file.path(output_directory, "DEGs_completos_Aysen_SPR_vs_SUM.csv"), row.names = FALSE)

firma_aysen <- Reduce(union, list(aysen_win_spr, aysen_win_sum, aysen_spr_sum))

# Contrastes internos de Magallanes
mag_aut_win <- leer_degs(file.path(input_directory, "DEGs_Magallanes_AUT_vs_WIN.csv"))
mag_aut_spr <- leer_degs(file.path(input_directory, "DEGs_Magallanes_AUT_vs_SPR.csv"))
mag_aut_sum <- leer_degs(file.path(input_directory, "DEGs_Magallanes_AUT_vs_SUM.csv"))
mag_win_spr <- leer_degs(file.path(input_directory, "DEGs_Magallanes_WIN_vs_SPR.csv"))
mag_win_sum <- leer_degs(file.path(input_directory, "DEGs_Magallanes_WIN_vs_SUM.csv"))
mag_spr_sum <- leer_degs(file.path(input_directory, "DEGs_Magallanes_SPR_vs_SUM.csv"))

# DEGs completos por contraste interno de Magallanes
write.csv(mag_aut_win, file.path(output_directory, "DEGs_completos_Magallanes_AUT_vs_WIN.csv"), row.names = FALSE)
write.csv(mag_aut_spr, file.path(output_directory, "DEGs_completos_Magallanes_AUT_vs_SPR.csv"), row.names = FALSE)
write.csv(mag_aut_sum, file.path(output_directory, "DEGs_completos_Magallanes_AUT_vs_SUM.csv"), row.names = FALSE)
write.csv(mag_win_spr, file.path(output_directory, "DEGs_completos_Magallanes_WIN_vs_SPR.csv"), row.names = FALSE)
write.csv(mag_win_sum, file.path(output_directory, "DEGs_completos_Magallanes_WIN_vs_SUM.csv"), row.names = FALSE)
write.csv(mag_spr_sum, file.path(output_directory, "DEGs_completos_Magallanes_SPR_vs_SUM.csv"), row.names = FALSE)

firma_magallanes <- Reduce(union, list(
  mag_aut_win, mag_aut_spr, mag_aut_sum,
  mag_win_spr, mag_win_sum, mag_spr_sum
))

# Firmas regionales completas
write.csv(firma_aysen, file.path(output_directory, "firma_Aysen.csv"), row.names = FALSE)
write.csv(firma_magallanes, file.path(output_directory, "firma_Magallanes.csv"), row.names = FALSE)

# Genes exclusivos y compartidos entre regiones
exclusivos_aysen <- setdiff(firma_aysen, firma_magallanes)
exclusivos_magallanes <- setdiff(firma_magallanes, firma_aysen)
compartidos <- intersect(firma_aysen, firma_magallanes)

write.csv(exclusivos_aysen, file.path(output_directory, "DEGs_exclusivos_Aysen.csv"), row.names = FALSE)
write.csv(exclusivos_magallanes, file.path(output_directory, "DEGs_exclusivos_Magallanes.csv"), row.names = FALSE)
write.csv(compartidos, file.path(output_directory, "DEGs_compartidos_regiones.csv"), row.names = FALSE)

# Resumen
cat("Enfoque 1 - DEGs completos WIN:", length(degs_win), "genes\n")
cat("Enfoque 1 - DEGs completos SPR:", length(degs_spr), "genes\n")
cat("Enfoque 1 - DEGs completos SUM:", length(degs_sum), "genes\n")
cat("Enfoque 1 - Núcleo estable (3 estaciones):", length(nucleo_estable), "genes\n")
cat("Enfoque 2 - Firma completa Aysén:", length(firma_aysen), "genes\n")
cat("Enfoque 2 - Firma completa Magallanes:", length(firma_magallanes), "genes\n")
cat("Enfoque 2 - Exclusivos Aysén:", length(exclusivos_aysen), "genes\n")
cat("Enfoque 2 - Exclusivos Magallanes:", length(exclusivos_magallanes), "genes\n")
cat("Enfoque 2 - Compartidos:", length(compartidos), "genes\n")
