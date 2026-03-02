#!/bin/bash
#
#SBATCH -p nadal-q
#SBATCH --chdir=/home/finiesta/prueba16febrero
#SBATCH --job-name=MeDuermo
#SBATCH --cpus-per-task=16  # Number of CPUs (max is 96)
#SBATCH --mem=100G 
#SBATCH --mail-type=ALL   # END/START/NONE
#SBATCH --mail-user=francisco.iniestam1@um.es
#SBATCH --error=salida_ejecucion/%j.err
#SBATCH --output=salida_ejecucion/%j.out

# Activate Anaconda work environment for OpenDrift
source ~/.bashrc
mamba activate samtools_env

bedsecuencias=$1
mkdir -p contigs
reference=$2
bamselect=$3
sample=$4
while IFS=$'\t' read -r -a columnas; do

	id="${columnas[0]}"
	chr="${columnas[1]}"
	start="${columnas[2]}"
    	secuencia="${columnas[3]}"
	Leftstart=$((start - 100000))
	Rightstart=$((start + 1))
	Rightend=$((Rightstart + 100000))

samtools faidx "$reference" "$chr:$Leftstart-$start" >"contigs/${id}_left.fa"
samtools faidx "$reference" "$chr:$Rightstart-$Rightend" >"contigs/${id}_right.fa"
echo "$secuencia" > "contigs/${id}_center.fa"
seqkit seq -w 0 "contigs/${id}_left.fa" | sed "s/^>.*/>chr_${id}/" > "contigs/${id}_left2.fa"
seqkit seq -w 0 "contigs/${id}_right.fa"|grep -v '^>'> "contigs/${id}_right2.fa" 
mkdir -p contigs/allcontigs
cat "contigs/${id}_left2.fa" "contigs/${id}_center.fa" "contigs/${id}_right2.fa">"contigs/${id}_merge.fa"
seqkit seq -w 60 "contigs/${id}_merge.fa"> "contigs/allcontigs/contig_${id}.fa"  

rm "contigs/${id}_merge.fa" "contigs/${id}_left2.fa" "contigs/${id}_center.fa" "contigs/${id}_right2.fa" "contigs/${id}_left.fa" "contigs/${id}_right.fa"
done < "$bedsecuencias"
mkdir -p contigs/custom_reference
cat contigs/allcontigs/*.fa>"contigs/custom_reference/custom.fa"
samtools faidx "contigs/custom_reference/custom.fa"

#Construccion de un bam recortado en base a fichero.bed
mkdir -p bamcustom

CORES=${SLURM_CPUS_PER_TASK:-1}
samtools view -@ "$CORES" -h -L "$bamselect" "$sample" | samtools fastq -@ "$CORES" -T '*' > bamcustom/custom.fastq
N_THREADS=12
minimap2 -k17 -t "$N_THREADS" -ax map-ont -y "contigs/custom_reference/custom.fa" "bamcustom/custom.fastq"|samtools sort -@ "$N_THREADS" -o "bamcustom/sorted_custom.bam"
samtools index "bamcustom/sorted_custom.bam"

#generamos plots con methylartist

mkdir -p methylartist_plots

while IFS=$'\t' read -r -a columnas; do

        id="${columnas[0]}"
        start=100000
        Leftstart=$((start - 1500))
        secuencia="${columnas[3]}"
        length=${#secuencia}
        insert_end=$((start + length))
        End=$((insert_end + 1500))

methylartist locus \
    -b "bamcustom/sorted_custom.bam" \
    --ref "contigs/custom_reference/custom.fa" \
    -i "chr_${id}:$Leftstart-$End" \
    -l "chr_${id}:$start-$insert_end" \
    -p 1,6,1,3,4 \
    --phased \
    --motif CG \
    -o "methylartist_plots/${id}.jpeg" 
done < "$bedsecuencias"
