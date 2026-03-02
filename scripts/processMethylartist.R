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


#En esta parte proceso los tsv de methylartist


LeftTSV<-read_tsv(snakemake@input[["left_tsv"]])
RightTSV<-read_tsv(snakemake@input[["right_tsv"]])

"Creo la columna de tasa de metilación diferencial por ventana y por flanco IZQ y DRCH"
col_ph1 <- grep("_ph1_m_methfrac$", names(RightTSV), value = TRUE)
col_ph2 <- grep("_ph2_m_methfrac$", names(RightTSV), value = TRUE)
col_ph1_read<- grep("_ph1_m_readcount$", names(RightTSV),value = TRUE)
col_ph2_read<- grep("_ph2_m_readcount$", names(RightTSV),value = TRUE)
RightTSV$diff <- abs(RightTSV[[col_ph1]] - RightTSV[[col_ph2]])
LeftTSV$diff <- abs(LeftTSV[[col_ph1]] - LeftTSV[[col_ph2]])

"Cambio la columna de seg_name por id para luego poder hacer un leftJoin con la tabla original"
RightTSV$id=RightTSV$seg_name 
LeftTSV$id=LeftTSV$seg_name 


Diferenciales_left<- LeftTSV %>%
  filter(diff>min_met_diff)%>%#aqui se puede tocar el umbral
  mutate(
    alelo1_Left=.data[[col_ph1]],
    alelo2_Left=.data[[col_ph2]],
    diff_L=diff,
    readsF1_L=.data[[col_ph1_read]],
    readsF2_L=.data[[col_ph2_read]],
    Ventana_L=seg_id,
    id=id
  )%>%
  select(
    id,
    Ventana_L,
    alelo1_Left,
    alelo2_Left,
    diff_L,
    readsF1_L,
    readsF2_L
  )


Diferenciales_right<- RightTSV %>%
  filter(diff>min_met_diff)%>%
  mutate(
    alelo1_Right=.data[[col_ph1]],
    alelo2_Right=.data[[col_ph2]],
    diff_R=diff,
    readsF1_R=.data[[col_ph1_read]],
    readsF2_R=.data[[col_ph2_read]],
    Ventana_R=seg_id,
    id=id
  )%>%
  select(
    id,
    Ventana_R,
    alelo1_Right,
    alelo2_Right,
    diff_R,
    readsF1_R,
    readsF2_R
  )


"Unimos ambos flancos en un dataframe final para luego vincular a toda la  info del retro"
resultadosALL<-full_join(Diferenciales_left, Diferenciales_right, by = "id")#Tabla con todas las diferencias detectadas(aqui habrá algo tipo tengo diferencias en lado izq pero no derecho y viceversa)
resultadosBOTH<-inner_join(Diferenciales_left, Diferenciales_right, by = "id")#tabla con efectos diferenciales a ambos lados

AllInfo<-resultadosALL%>%left_join(process_data, by = "id")#Uno esta tabla con toda la info del insertionTable

saveRDS(AllInfo, file = snakemake@output[["AllInfo"]])

#Añado lo que me dijo Javi de UNA SUPERTABLA



TOTAL_left<- LeftTSV %>%mutate(
    alelo1_Left=.data[[col_ph1]],
    alelo2_Left=.data[[col_ph2]],
    diff_L=diff,
    readsF1_L=.data[[col_ph1_read]],
    readsF2_L=.data[[col_ph2_read]],
    Ventana_L=seg_id,
    id=id
  )%>%
  select(
    id,
    Ventana_L,
    alelo1_Left,
    alelo2_Left,
    diff_L,
    readsF1_L,
    readsF2_L
  )


TOTAL_right<-RightTSV %>%mutate(
    alelo1_Right=.data[[col_ph1]],
    alelo2_Right=.data[[col_ph2]],
    diff_R=diff,
    readsF1_R=.data[[col_ph1_read]],
    readsF2_R=.data[[col_ph2_read]],
    Ventana_R=seg_id,
    id=id
  )%>%
  select(
    id,
    Ventana_R,
    alelo1_Right,
    alelo2_Right,
    diff_R,
    readsF1_R,
    readsF2_R
  )



resultados2<-full_join(TOTAL_left, TOTAL_right, by = "id")
resultadosTOTAL<-resultados2%>%left_join(process_data, by = "id")
saveRDS(resultadosTOTAL, file = snakemake@output[["TOTAL_INS"]])



#Generamos un bed-query para hacer consultas en bcftools

bcftools<-AllInfo%>%
  mutate(chr=chr,
         init=(start -100)-1,#este valor es muy dependiente de como este el mosaicismo de tu variable.
         end =(start +100),
  )%>%
  select(
    chr,
    init,
    end
  )

write.table(bcftools, 
            file = snakemake@output[["Query_bcftools"]], 
            sep = "\t",            # Define el tabulador como separador
            quote = FALSE,         # Evita que ponga comillas en los textos
            row.names = FALSE,     # No escribe los números de fila
            col.names = FALSE) 




#creamos consultas pasra TODAS


bcftools2<-resultadosTOTAL%>%
  mutate(chr=chr,
         init=(start -100)-1,#este valor es muy dependiente de como este el mosaicismo de tu variable.
         end =(start +100),
  )%>%
  select(
    chr,
    init,
    end
  )

write.table(bcftools2,
            file = snakemake@output[["Query_bcftools2"]],
            sep = "\t",            # Define el tabulador como separador
            quote = FALSE,         # Evita que ponga comillas en los textos
            row.names = FALSE,     # No escribe los números de fila
            col.names = FALSE)

























#Construimos un bed para plotear todos los resultados y otro con suppReads seleccionado para crear contigs con insertos.

"Asu vez aquí tengo las que voy a querer plotear con la secuencia del retro dentro es decir aquellas que se realizará un constructo artificial:
Se creará una referencia genómica donde cada chromosoma será un contig del hg38 + retrotransposon en el punto de insercción detectad.
Dichas inserciones tienen que pasar unos umbrales mínimos para plotear con methylartist dado que 3 reads es extramadamente bajo por defecto pondré 5 pero será customizable
"
candidatas<- AllInfo %>%
  filter((readsF1_L>=min_reads & readsF2_L>=min_reads)|(readsF1_R>=min_reads & readsF2_R>=min_reads))%>%
  mutate(
    chr=chr,
    InsertPoint=start,
    Class=clase,
    Subclass=subclase,
    Genotype=genotipo,
    MetLeftWindowDiff=diff_L,
    LeftHp1=alelo1_Left,
    LeftHp2=alelo2_Left,
    SuppReadsHp1L=readsF1_L,
    SuppReadsHp2L=readsF2_L,
    MetRightWindowDiff=diff_R,
    RightHp1=alelo1_Right,
    RightHp2=alelo2_Right,
    SuppReadsHp1R=readsF1_R,
    SuppReadsHp2R=readsF2_R,
    sequence=sequence,
    Id=id
  )%>%
  select(
    Id,
    chr,
    InsertPoint,
    Class,
    Subclass,
    Genotype,
    MetLeftWindowDiff,
    LeftHp1,
    LeftHp2,
    SuppReadsHp1L,
    SuppReadsHp2L,
    MetRightWindowDiff,
    RightHp1,
    RightHp2,
    SuppReadsHp1R,
    SuppReadsHp2R,
    sequence
  )

"Generamos un Bed con posiciones de interés para lluego extraer reads con: samtools view -h -L regiones.bed"
BedCandidatas<- candidatas%>%
  mutate(
    chr=chr,
    inicio=InsertPoint-100000,
    final=InsertPoint+100000,
    seqid=Id
  )%>%
  select(
    chr,
    inicio,
    final
  )

BedSecuencias<- candidatas%>% mutate(
  InserSeq=InsertPoint -1
)%>%select(
  Id,
  chr,
  InserSeq,
  sequence
)

write.table(BedCandidatas, 
            file = snakemake@output[["BedCandidatas"]], 
            sep = "\t",            # Define el tabulador como separador
            quote = FALSE,         # Evita que ponga comillas en los textos
            row.names = FALSE,     # No escribe los números de fila
            col.names = FALSE) 


write.table(BedSecuencias, 
            file = snakemake@output[["BedSecuencias"]], 
            sep = "\t",            # Define el tabulador como separador
            quote = FALSE,         # Evita que ponga comillas en los textos
            row.names = FALSE,     # No escribe los números de fila
            col.names = FALSE) 
