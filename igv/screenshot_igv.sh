#!/usr/bin/env bash

TEMPLATE_BAT="$1"
BAM1="$2"
BAM2="$3"
OUTDIR="${4:-screenshots}"
GENOME="${5:-hg38}"   # opcional, por defecto hg38

LOG="${OUTDIR%/}/igv.log"
mkdir -p "$OUTDIR"

# bat temporal

sed \
  -e "s|__BAM1__|$BAM1|g" \
  -e "s|__BAM2__|$BAM2|g" \
  -e "s|__OUTDIR__|$OUTDIR|g" \
  -e "s|__GENOME__|$GENOME|g" \
  "$TEMPLATE_BAT" > custom.bat

mamba env create -f igv.yaml -y

# 2. Ejecutar IGV directamente dentro del entorno sin "activarlo" manualmente
xvfb-run --auto-servernum mamba run -n igv_Retromet igv -b custom.bat

# 3. Mostrar log y borrar
echo "Log: $LOG"
mamba env remove -n igv_Retromet -y