
bedsecuencias=$1
reference=$2
bamselect=$3
sample=$4
outdir=$5

# Usa cores de SLURM si existen; si no, 1
CORES=${SLURM_CPUS_PER_TASK:-1}

# Directorios dentro de outdir
contigs_dir="$outdir/contigs"
allcontigs_dir="$contigs_dir/allcontigs"
customref_dir="$contigs_dir/custom_reference"
bamcustom_dir="$outdir/bamcustom"
plots_dir="$outdir/methylartist_plots"

mkdir -p "$contigs_dir" "$allcontigs_dir" "$customref_dir" "$bamcustom_dir" "$plots_dir"

# 1) Construcción de contigs + referencia custom
while IFS=$'\t' read -r -a columnas; do
    [[ ${#columnas[@]} -lt 4 ]] && continue

    id="${columnas[0]}"
    chr="${columnas[1]}"
    start="${columnas[2]}"
    secuencia="${columnas[3]}"

    Leftstart=$((start - 100000))
    Rightstart=$((start + 1))
    Rightend=$((Rightstart + 100000))

    # Evitar coordenadas negativas
    if (( Leftstart < 1 )); then Leftstart=1; fi

    samtools faidx "$reference" "$chr:$Leftstart-$start" > "$contigs_dir/${id}_left.fa"
    samtools faidx "$reference" "$chr:$Rightstart-$Rightend" > "$contigs_dir/${id}_right.fa"
    printf "%s\n" "$secuencia" > "$contigs_dir/${id}_center.fa"

    seqkit seq -w 0 "$contigs_dir/${id}_left.fa" \
      | sed "s/^>.*/>chr_${id}/" > "$contigs_dir/${id}_left2.fa"

    seqkit seq -w 0 "$contigs_dir/${id}_right.fa" \
      | grep -v '^>' > "$contigs_dir/${id}_right2.fa"

    cat "$contigs_dir/${id}_left2.fa" \
        "$contigs_dir/${id}_center.fa" \
        "$contigs_dir/${id}_right2.fa" > "$contigs_dir/${id}_merge.fa"

    seqkit seq -w 60 "$contigs_dir/${id}_merge.fa" > "$allcontigs_dir/contig_${id}.fa"

    rm -f \
      "$contigs_dir/${id}_merge.fa" \
      "$contigs_dir/${id}_left2.fa" \
      "$contigs_dir/${id}_center.fa" \
      "$contigs_dir/${id}_right2.fa" \
      "$contigs_dir/${id}_left.fa" \
      "$contigs_dir/${id}_right.fa"
done < "$bedsecuencias"

cat "$allcontigs_dir"/*.fa > "$customref_dir/custom.fa"
samtools faidx "$customref_dir/custom.fa"

# 2) Construcción de BAM custom (recortado + realineado)
samtools view -@ "$CORES" -h -L "$bamselect" "$sample" \
  | samtools fastq -@ "$CORES" -T '*' > "$bamcustom_dir/custom.fastq"

# minimap2 usa threads
minimap2 -k17 -t "$CORES" -ax map-ont -y "$customref_dir/custom.fa" "$bamcustom_dir/custom.fastq" \
  | samtools sort -@ "$CORES" -o "$bamcustom_dir/sorted_custom.bam"

samtools index "$bamcustom_dir/sorted_custom.bam"

# 3) Plots con methylartist (sobre el BAM custom y la custom ref)
while IFS=$'\t' read -r -a columnas; do
    [[ ${#columnas[@]} -lt 4 ]] && continue

    id="${columnas[0]}"
    secuencia="${columnas[3]}"

    # Aquí estás forzando start=100000 (lo mantengo igual)
    start=100000
    Leftstart=$((start - 1500))
    length=${#secuencia}
    insert_end=$((start + length))
    End=$((insert_end + 1500))

    methylartist locus \
        -b "$bamcustom_dir/sorted_custom.bam" \
        --ref "$customref_dir/custom.fa" \
        -i "chr_${id}:$Leftstart-$End" \
        -l "chr_${id}:$start-$insert_end" \
        -p 1,6,1,3,4 \
        --phased \
        --motif CG \
        -o "$plots_dir/${id}.jpeg"
done < "$bedsecuencias"
