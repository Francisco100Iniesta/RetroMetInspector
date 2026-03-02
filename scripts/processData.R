samples = snakemake@params[["samples"]]
library(data.table)
library(dplyr)
library(readr)
library(DT)
library(tidyr)
#Cargamos los datos de insertionTable generados en la carpeta de .rds por REtroinspector

met_window=snakemake@params[["met_window"]] #aqui tomaremos del configWindow la ventana de análisis lateral.
min_met_diff=snakemake@params[["min_met_diff"]]
min_reads=snakemake@params[["min_reads"]]
current_sample <- snakemake@wildcards[["sample"]]
samples<-snakemake@params[["samples"]]
data=readRDS(snakemake@input[["data"]])
"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

                  PRIMERA FASE: PROCESADO DE INSERCIONES DE DATOS DE RETROINSOPECTOR"





"Hay demasiada info y ademas hay que filtrar por el SUPP_Trusty osea que este detectado en algun caller"

pos<-match(current_sample,samples)

process_data <- data %>%
  filter(vapply(SUPP_mask_trusty, function(x) x[pos], logical(1))) %>%
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
    genotipo
  )

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
