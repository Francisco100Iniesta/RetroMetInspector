bedsecuencias=$1
reference=$2
sample=$3
gtf=$4
outdir=$5

mkdir -p methylartist_plots_ref

while IFS=$'\t' read -r -a columnas; do

        id="${columnas[0]}"
        chr="${columnas[1]}"
        start="${columnas[2]}"
        insert_end=$((start + 50))
        End=$((insert_end + 950))
	Begin=$((start - 1000))

methylartist locus \
    -b "$sample" \
    --ref "$reference" \
    -i "$chr:$Begin-$End" \
    -p 1,6,1,3,4 \
    --phased \
    --motif CG \
    -l "chr_${id}:$start-$insert_end" \
    -g "$gtf" \
    --labelgenes \
    -o "$outdir/${id}.jpeg"
done < "$bedsecuencias"

