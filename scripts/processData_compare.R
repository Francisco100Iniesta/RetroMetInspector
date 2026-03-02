library(data.table)
library(dplyr)
library(readr)
library(DT)
library(tidyr)
library(glue)
#Cargamos los datos de insertionTable generados en la carpeta de .rds por REtroinspector

met_window=snakemake@params[["met_window"]] #aqui tomaremos del configWindow la ventana de análisis lateral.
min_met_diff=snakemake@params[["min_met_diff"]]
min_reads=snakemake@params[["min_reads"]]
current_sample <- snakemake@wildcards[["sample"]]
samples<-snakemake@params[["samples"]]
data=readRDS(snakemake@input[["data"]])
soi=snakemake@params$soi
"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

                  PRIMERA FASE: PROCESADO DE INSERCIONES DE DATOS DE RETROINSOPECTOR"





"Hay demasiada info y ademas hay que filtrar por el SUPP_Trusty osea que este detectado en algun caller"

sample1<-soi[1]
sample2<-soi[2]
pos1<-match(sample1,samples)
pos2<-match(sample2,samples)

process_data_Sample1 <- data %>%
  filter(vapply(SUPP_mask_trusty, function(x) x[pos1], logical(1)),
         unlist(lapply(genotype, function(x) x[pos1] == 1)),
         SUPP_min3==1)

process_data_Sample2 <- data %>%
  filter(vapply(SUPP_mask_trusty, function(x) x[pos2], logical(1)),
         unlist(lapply(genotype, function(x) x[pos2] == 1)),
         SUPP_min3==1)

Sample1_df <- anti_join(process_data_Sample1, process_data_Sample2, by = "seqId") %>% 
  mutate(sample = sample1)


Sample2_df <- anti_join(process_data_Sample2, process_data_Sample1, by = "seqId") %>% 
  mutate(sample = sample2)

Unshared_insertions=bind_rows(Sample1_df,Sample2_df)


"Hasta aqui se formaria la lista de unshared"


process_data <- Unshared_insertions %>%
  filter(SUPP_trusty>0) %>%
  mutate(
    start=start,
    chr = seqnames,
    Svlen = SVLEN,
    longitud_secuencia = nchar(as.character(vcf_alt)) - 1,
    posicion = paste(seqnames, format(start, trim = TRUE), sep = ":"),
    id = seqId,
    strand=strand,
    clase=repeat.class,
    subclase=repeat.subclass,
    genotipo=genotype,
    sequence=vcf_alt,
    genotype=genotype
  ) %>%
  select(
    chr, 
    start,
    Svlen, 
    longitud_secuencia, 
    posicion, 
    id, 
    sequence,
    strand,
    genotipo,
    clase,
    subclase,
    genotipo,
    sample
  )


saveRDS(Unshared_insertions,file = snakemake@output[["unshared_insertion"]])

"Como con las pruebas que realicé prácticamente ventas de 2kbases enteras diluyen el efecto de la metilación es mejor acotar dos ventanas por cada flanco desde el 
punto de la inserción de esa forma contienes los valores reales de metilación por silenciamiento en el genoma que tienden a diluirse de 500 a 1000 bases.
Estoy probando flancos izq y derechos a 1000bases pero creo que tengo que ir probando"


leftBED<- process_data %>%
  mutate(
    end_L=start,
    start_L=(start-met_window)-1  ,#ajustamos al formato Bed que es 0 based
    chr = chr,
    id = id
  ) %>%
  select(
    chr, 
    start_L,
    end_L,
    id
  )

rightBED <- process_data %>%
  mutate(
    start_R = start-1,
    end_R   = (start + met_window)
  ) %>%
  select(chr, start_R, end_R, id)



write.table(leftBED, 
            file = snakemake@output[["left_bed"]], 
            sep = "\t",            # Define el tabulador como separador
            quote = FALSE,         # Evita que ponga comillas en los textos
            row.names = FALSE,     # No escribe los números de fila
            col.names = FALSE)     # No escribe los encabezados

write.table(rightBED, 
            file = snakemake@output[["right_bed"]], 
            sep = "\t",            # Define el tabulador como separador
            quote = FALSE,         # Evita que ponga comillas en los textos
            row.names = FALSE,     # No escribe los números de fila
            col.names = FALSE)     # No escribe los encabezados


#Construimos un fichero de coordenadas para IGV
genome <- "hg38"
bams <- c(glue("alns/{sample1}.bam"), glue("alns/{sample2}.bam"))
out_dir <- glue("igv/screenshots{sample1}vs{sample2}")

header <- c(
  "new",
  "genome __GENOME__",
  "load __BAM1__",
  "load __BAM2__",
  "snapshotDirectory __OUTDIR__",
  "group PHASE"
)

body<- process_data %>%mutate(
  snap=glue("snapshot {id}.png"),
  pos=glue("goto {chr}:{start - 1000}-{start +1000}"))%>%
  select(pos, snap) %>%
  t() %>%
  as.vector()

writeLines(c(header, body, "exit"), snakemake@output[["bat_file"]])
