# MT_DEGsSalmoSalar

Proyecto de memoria de título de Tomás Cid.

## Alcance dentro del pipeline bioinformático

Este repositorio cubre las etapas **posteriores a la cuantificación** de la expresión génica. Parte de las matrices de conteos generadas previamente con FeatureCounts; no incluye el procesamiento de lecturas crudas, el control de calidad inicial, el recorte de lecturas, el alineamiento ni la cuantificación.

A partir de esas matrices, los scripts realizan el análisis de expresión diferencial con DESeq2 y comparan los conjuntos de genes expresados diferencialmente. Luego, los resultados se utilizan para la interpretación funcional con ClueGO y para generar tablas, gráficos y redes funcionales.

Los informes MultiQC de las etapas previas se conservan como referencia, pero este repositorio no ejecuta esas etapas del pipeline.

## Organización del repositorio

- `codigo/scripts/`: scripts de análisis en R y el notebook de exploración de sensores.
- `datos/featurecounts/`: matrices de conteos usadas en el análisis diferencial.
- `datos/metadatos/`: planillas de nombres y metadatos de archivos.
- `resultados/DEGs/originales/`: tablas originales de expresión diferencial.
- `resultados/DEGs/venn/`: resultados de intersecciones, exclusividad y núcleo regional.
- `resultados/ClueGO/<región>/`: tablas y matrices funcionales exportadas por ClueGO.
- `resultados/control_calidad/MultiQC/`: informes MultiQC, separados en `PreFastP` y `PostFastP`.
- `graficos/deseq2/`: gráficos PCA, volcano, heatmap e intersecciones del análisis DESeq2.
- `graficos/ClueGO/<región>/`: figuras PNG y SVG exportadas por ClueGO.
- `graficos/redes_funcionales/`: redes funcionales generales.
- `graficos/exploracion_sensores/`: figuras que genera el notebook de sensores.

Desde la raíz del repositorio, ejecuta los scripts con `Rscript codigo/scripts/Análisis_Transcriptomico.R` y `Rscript codigo/scripts/Diagramas_Venn.R`. Ambos localizan la raíz del proyecto a partir de la ruta del propio script. El análisis DESeq2 guarda sus tablas nuevas en `resultados/DEGs/deseq2/`; el script de diagramas/intersecciones lee las tablas de `resultados/DEGs/originales/` y guarda sus resultados en `resultados/DEGs/venn/`.
