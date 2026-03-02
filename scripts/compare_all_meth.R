library(dplyr)
library(readr)
library(DT)
library(tidyr)
library(glue)
#Cargamos los datos de insertionTable generados en la carpeta de .rds por REtroinspector
MetWindow=snakemake@params[["met_window"]]
min_met_diff=snakemake@params[["min_met_diff"]]
min_reads=snakemake@params[["min_reads"]]
samples<-snakemake@params[["samples"]]
soi=snakemake@params$soi
sample1<-soi[1]
sample2<-soi[2]
TSVS2_Left=read_tsv(snakemake@input[["left_tsv2"]])
TSVS2_Right=read_tsv(snakemake@input[["right_tsv2"]])
TSVS1_Left=read_tsv(snakemake@input[["left_tsv1"]])
TSVS1_Right=read_tsv(snakemake@input[["right_tsv1"]])
unshared=readRDS(snakemake@input[["unshared_insertion"]])
"**********************************************************************"
#PARA EL S2

col_ph1_S2 <- grep("_ph1_m_methfrac$", names(TSVS2_Left), value = TRUE)
col_ph2_S2 <- grep("_ph2_m_methfrac$", names(TSVS2_Left), value = TRUE)
col_ph1_read_S2<- grep("_ph1_m_readcount$", names(TSVS2_Left),value = TRUE)
col_ph2_read_S2<- grep("_ph2_m_readcount$", names(TSVS2_Left),value = TRUE)
TSVS2_Left$diff <- abs(TSVS2_Left[[col_ph1_S2]] - TSVS2_Left[[col_ph2_S2]])
TSVS2_Right$diff <- abs(TSVS2_Right[[col_ph1_S2]] - TSVS2_Right[[col_ph2_S2]])


TSVS2_Left<-TSVS2_Left%>%mutate(
  
  InsertPoint=glue("{seg_chrom}:{seg_end}"),
  S2_L_DIFF=diff,
  Hp1_L_S2=.data[[col_ph1_S2]],
  Hp2_L_S2=.data[[col_ph2_S2]],
  readsF1_L_S2=.data[[col_ph1_read_S2]],
  readsF2_L_S2=.data[[col_ph2_read_S2]],
)%>%select(
  
  InsertPoint,
  Hp1_L_S2,
  Hp2_L_S2,
  S2_L_DIFF,
  readsF1_L_S2,
  readsF2_L_S2
)


TSVS2_Right<-TSVS2_Right%>%mutate(
  
  InsertPoint=glue("{seg_chrom}:{seg_end - MetWindow}"),
  S2_R_DIFF=diff,
  Hp1_R_S2=.data[[col_ph1_S2]],
  Hp2_R_S2=.data[[col_ph2_S2]],
  readsF1_R_S2=.data[[col_ph1_read_S2]],
  readsF2_R_S2=.data[[col_ph2_read_S2]],
)%>%select(
  
  InsertPoint,
  Hp1_R_S2,
  Hp2_R_S2,
  S2_R_DIFF,
  readsF1_R_S2,
  readsF2_R_S2
)





"**********************************************************************"
#PARA EL S1

col_ph1_S1 <- grep("_ph1_m_methfrac$", names(TSVS1_Left), value = TRUE)
col_ph2_S1<- grep("_ph2_m_methfrac$", names(TSVS1_Left), value = TRUE)
col_ph1_read_S1<- grep("_ph1_m_readcount$", names(TSVS1_Left),value = TRUE)
col_ph2_read_S1<- grep("_ph2_m_readcount$", names(TSVS1_Left),value = TRUE)
TSVS1_Left$diff <- abs(TSVS1_Left[[col_ph1_S1]] - TSVS1_Left[[col_ph2_S1]])
TSVS1_Right$diff <- abs(TSVS1_Right[[col_ph1_S1]] - TSVS1_Right[[col_ph2_S1]])



TSVS1_Left<-TSVS1_Left%>%mutate(
  
  InsertPoint=glue("{seg_chrom}:{seg_end}"),
  S1_L_DIFF=diff,
  Hp1_L_S1=.data[[col_ph1_S1]],
  Hp2_L_S1=.data[[col_ph2_S1]],
  readsF1_L_S1=.data[[col_ph1_read_S1]],
  readsF2_L_S1=.data[[col_ph2_read_S1]],
)%>%select(
  
  InsertPoint,
  Hp1_L_S1,
  Hp2_L_S1,
  S1_L_DIFF,
  readsF1_L_S1,
  readsF2_L_S1
)


TSVS1_Right<-TSVS1_Right%>%mutate(
  
  InsertPoint=glue("{seg_chrom}:{seg_end - MetWindow}"),
  S1_R_DIFF=diff,
  Hp1_R_S1=.data[[col_ph1_S1]],
  Hp2_R_S1=.data[[col_ph2_S1]],
  readsF1_R_S1=.data[[col_ph1_read_S1]],
  readsF2_R_S1=.data[[col_ph2_read_S1]],
)%>%select(
  
  InsertPoint,
  Hp1_R_S1,
  Hp2_R_S1,
  S1_R_DIFF,
  readsF1_R_S1,
  readsF2_R_S1
)

S1_meth=left_join(TSVS1_Right,TSVS1_Left,by = "InsertPoint")
S2_meth=left_join(TSVS2_Right,TSVS2_Left,by = "InsertPoint")





"HASTA AQUI TIENES LEFT AND RIGHT POR SAMPLE CON UN ID PARA AÑADIR AHORA A LA DEL S1 LE TIENES QUE DAR LA DE S2 QUE NO TENDRA ETC"
"SEGUNDA FASE: ANALISIS DE RESULTADOS DE METHYLARTIST + Bed bcftools"

unshared$Id=unshared$seqId
"CARGAMOS RUTAS DE LOS RDS CON TODOS LOS INSERTOS DETECTADOS"
FinalS1_all=readRDS(snakemake@input[["FinalTable2_S1"]])
FinalS2_all=readRDS(snakemake@input[["FinalTable2_S2"]])

SUBSET_S1_all=subset(unshared,sample==sample1,select = c(sample,Id))
SUBSET_S2_all=subset(unshared,sample==sample2,select = c(sample,Id))

S1_ALL=inner_join(FinalS1_all,SUBSET_S1_all,by="Id")
S2_ALL=inner_join(FinalS2_all,SUBSET_S2_all,by="Id")
#Unimos con el seqD de la insercion para fusionar informacion del que supuestamente no la tiene
S1_ALL=left_join(S1_ALL,S2_meth,by = "InsertPoint")
S2_ALL=left_join(S2_ALL,S1_meth,by = "InsertPoint")

saveRDS(S1_ALL,file = snakemake@output[["S1_ALL"]])
saveRDS(S2_ALL,file = snakemake@output[["S2_ALL"]])



#Lo mismo para los Criterios de Filtro estricto
FinalS1_Filter=readRDS(snakemake@input[["FinalTable_S1"]])
FinalS2_Filter=readRDS(snakemake@input[["FinalTable_S2"]])

SUBSET_S1=subset(unshared,sample==sample1,select = c(sample,Id))
SUBSET_S2=subset(unshared,sample==sample2,select = c(sample,Id))

S1_filtered=inner_join(FinalS1_Filter,SUBSET_S1,by="Id")
S2_filtered=inner_join(FinalS2_Filter,SUBSET_S2,by="Id")

S1_ALL_filtered=left_join(S1_filtered,S2_meth,by = "InsertPoint")
S2_ALL_filtered=left_join(S2_filtered,S1_meth,by = "InsertPoint")

saveRDS(S1_ALL_filtered,file = snakemake@output[["S1_filtered"]])
saveRDS(S2_ALL_filtered,file = snakemake@output[["S2_filtered"]])