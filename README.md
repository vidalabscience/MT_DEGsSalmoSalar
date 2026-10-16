# MT_DEGsSalmoSalar

Proyecto de memoria de título de Tomás Cid.

## Agradecimientos

Este trabajo se realizó con el apoyo de los siguientes proyectos:

- Fondecyt Iniciación 11261470.
- ANID – Programa Iniciativa Científica Milenio – NCN2024_047.

## Alcance y organización

Este repositorio contiene el código y los **gráficos** de las etapas posteriores a la cuantificación de expresión génica, además de los archivos de resultados exportados por ClueGO. No incluye matrices de conteos ni tablas CSV intermedias o finales que se utilizaron durante el proyecto.

- `codigo/scripts/`: scripts R y notebook de exploración de datos de sensores.
- `resultados/ClueGO/<región>/`: tablas, matrices y reportes exportados por ClueGO.
- `resultados/graficos/`: gráficos DESeq2, ClueGO, redes funcionales y exploración de sensores.

## Ejecución del código

Ejecuta los scripts R desde la raíz del repositorio:

```sh
Rscript codigo/scripts/Análisis_Transcriptomico.R
Rscript codigo/scripts/Diagramas_Venn.R
```

Antes de ejecutar `Análisis_Transcriptomico.R`, ajusta `select_regions` y `select_seasons` según el contraste y proporciona externamente las matrices de conteo requeridas. Los gráficos se guardan en `resultados/graficos/deseq2/`.

`Diagramas_Venn.R` requiere las tablas CSV de expresión diferencial generadas previamente y produce CSV con las intersecciones. Los CSV que genera y los CSV usados como entrada permanecen locales y no se incluyen en el repositorio. Los gráficos de Venn conservados están en `resultados/graficos/Venn/`.

Los resultados de ClueGO incluidos aquí se obtuvieron con el plugin ClueGO en Cytoscape, usando como entrada los CSV generados por `Diagramas_Venn.R`.

Para `exploracion_sensores.ipynb`, abre el notebook en Jupyter o VS Code y ejecuta las celdas en orden. Requiere un archivo Excel de sensores externo, indicado en `SENSOR_DATASET_PATH` o en la ruta `dataset_sensores/Stats_Extraidos.xlsx`. Sus gráficos se guardan en `resultados/graficos/exploracion_sensores/`.
