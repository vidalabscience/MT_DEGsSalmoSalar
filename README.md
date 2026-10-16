# MT_DEGsSalmoSalar

Proyecto de memoria de título de Tomás Cid.

## Agradecimientos

Este trabajo se realizó con el apoyo de los siguientes proyectos:

- Fondecyt Iniciación 11261470.
- ANID – Programa Iniciativa Científica Milenio – NCN2024_047.

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

## Ejecución del código

Ejecuta los scripts R desde la raíz del repositorio:

```sh
Rscript codigo/scripts/Análisis_Transcriptomico.R
Rscript codigo/scripts/Diagramas_Venn.R
```

Antes de ejecutar `Análisis_Transcriptomico.R`, ajusta `select_regions` y `select_seasons` según el contraste: una región y dos estaciones para comparar estaciones, o dos regiones y una estación para comparar regiones. El orden indicado define el orden de la comparación. Las tablas de DESeq2 se guardan en `resultados/DEGs/deseq2/` y sus gráficos en `graficos/deseq2/`.

`Diagramas_Venn.R` procesa las tablas ubicadas en `resultados/DEGs/originales/` y guarda los resultados en `resultados/DEGs/venn/`.

Para `exploracion_sensores.ipynb`, abre el notebook en Jupyter o VS Code y ejecuta las celdas en orden. Requiere el archivo `dataset_sensores/Stats_Extraidos.xlsx`, que no está incluido en el repositorio.
