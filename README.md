# MT_DEGsSalmoSalar

Proyecto de memoria de título de Tomás Cid.

## Agradecimientos

Este trabajo se realizó con el apoyo de los siguientes proyectos:

- Fondecyt Iniciación 11261470.
- ANID – Programa Iniciativa Científica Milenio – NCN2024_047.

## Alcance dentro del pipeline bioinformático

Este repositorio incluye el código y los resultados de las etapas **posteriores a la cuantificación** de la expresión génica. Las matrices de conteos se deben obtener externamente; no se incluyen datos de entrada ni el procesamiento de lecturas crudas, el control de calidad inicial, el recorte, el alineamiento o la cuantificación.

Los análisis realizan expresión diferencial con DESeq2, comparación de conjuntos de genes e interpretación funcional con ClueGO. Se incluyen tablas, informes MultiQC de referencia, gráficos y redes funcionales; las etapas previas del pipeline no se ejecutan aquí.

## Organización del repositorio

- `codigo/scripts/`: scripts de análisis en R y notebook de exploración de sensores.
- `resultados/DEGs/originales/`: resultados DEGs de contrastes previos, reutilizados para calcular las intersecciones.
- `resultados/DEGs/venn/`: resultados de intersecciones, exclusividad y núcleo regional.
- `resultados/ClueGO/<región>/`: tablas y matrices funcionales exportadas por ClueGO.
- `resultados/control_calidad/MultiQC/`: informes MultiQC de referencia, separados en `PreFastP` y `PostFastP`.
- `resultados/graficos/`: gráficos DESeq2, figuras ClueGO, redes funcionales y gráficos de exploración de sensores.

## Ejecución del código

Ejecuta los scripts R desde la raíz del repositorio:

```sh
Rscript codigo/scripts/Análisis_Transcriptomico.R
Rscript codigo/scripts/Diagramas_Venn.R
```

Antes de ejecutar `Análisis_Transcriptomico.R`, ajusta `select_regions` y `select_seasons` según el contraste: una región y dos estaciones para comparar estaciones, o dos regiones y una estación para comparar regiones. El orden indicado define el orden de la comparación. Proporciona las matrices `counts_Aysen.txt` y `counts_Magallanes.txt` externamente en la carpeta `datos/featurecounts/` (excluida de Git) o define la variable `FEATURECOUNTS_DIR` con la ruta a esa carpeta. Las tablas de DESeq2 y sus gráficos se guardan en `resultados/DEGs/deseq2/` y `resultados/graficos/deseq2/`, respectivamente.

`Diagramas_Venn.R` procesa las tablas ubicadas en `resultados/DEGs/originales/` y guarda los resultados en `resultados/DEGs/venn/`.

Para `exploracion_sensores.ipynb`, abre el notebook en Jupyter o VS Code y ejecuta las celdas en orden. Requiere un archivo Excel de sensores externo, indicado en `SENSOR_DATASET_PATH` o en la ruta `dataset_sensores/Stats_Extraidos.xlsx`. Sus gráficos se guardan en `resultados/graficos/exploracion_sensores/`.
