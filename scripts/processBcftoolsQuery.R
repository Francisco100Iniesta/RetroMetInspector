library(dplyr)
library(readr)
library(DT)
library(tidyr)

#EN ESTE script se filtran los resultados de bcftools y se asocian genotipo y haplotipo
nombres<-c("chr","start","Type","Svlen","HP","GT","GQ","DR","DV")
queryBcftools<-read_tsv(snakemake@input[["Query"]], col_names=nombres)
filtradoBcftools<- queryBcftools %>%
  filter(Type == "INS")#Solo queremos Variantes estructurales catalogadas con Inserciones

filtradoBcftools$Svlen=as.numeric(filtradoBcftools$Svlen)


#lO MISMO PARA LA TABLA CON TODA LA INFORMACION
nombres<-c("chr","start","Type","Svlen","HP","GT","GQ","DR","DV")
queryBcftools2<-read_tsv(snakemake@input[["Query2"]], col_names=nombres)
filtradoBcftools2<- queryBcftools2 %>%
  filter(Type == "INS")#Solo queremos Variantes estructurales catalogadas con Inserciones

filtradoBcftools2$Svlen=as.numeric(filtradoBcftools2$Svlen)






window_bp <- 1000     # ventana +/- bp
use_svlen <- TRUE

# columnas clave
col_chr_small <- "chr"
col_pos_small <- "start"
col_len_small <- "Svlen"   # si use_svlen=TRUE

col_chr_big   <- "chr"
col_pos_big   <- "start"
col_len_big   <- "longitud_secuencia"

# qué campos quieres copiar del SMALL al BIG
fields_small_to_big <- c("HP", "GT") 

AllInfo<-readRDS(snakemake@input[["AllInfo"]])
df_big<-AllInfo
df_small<-filtradoBcftools
# ---------------- PREP: crear columnas en BIG ----------------


for (f in fields_small_to_big) {
  newcol <- paste0("small_", f)
  if (!newcol %in% names(df_big)) df_big[[newcol]] <- NA
}

# opcional: guardar también info del match
if (!"matched_from_small" %in% names(df_big)) df_big$matched_from_small <- NA_integer_
if (!"match_dist_pos" %in% names(df_big)) df_big$match_dist_pos <- NA_integer_
if (use_svlen && !"match_dist_len" %in% names(df_big)) df_big$match_dist_len <- NA_integer_

# ---------------- LOOP: small -> best big ----------------
for (i in seq_len(nrow(df_small))) {
  
  chr_i <- df_small[[col_chr_small]][i]
  pos_i <- df_small[[col_pos_small]][i]
  
  # candidatos en BIG por chr + ventana
  idx <- which(
    df_big[[col_chr_big]] == chr_i &
      abs(df_big[[col_pos_big]] - pos_i) <= window_bp
  )
  
  if (length(idx) == 0) next
  
  # scoring
  dist_pos <- abs(df_big[[col_pos_big]][idx] - pos_i)
  
  if (use_svlen) {
    len_i <- df_small[[col_len_small]][i]
    dist_len <- abs(df_big[[col_len_big]][idx] - len_i)
  } else {
    dist_len <- rep(0, length(idx))
  }
  
  # (opcional) si quieres usar DV/QUAL del BIG para desempatar:
  dv <- if ("DV" %in% names(df_big)) df_big$DV[idx] else rep(0, length(idx))
  dv[is.na(dv)] <- 0
  
  # elegir best: menor dist_pos, luego menor dist_len, luego mayor DV
  ord <- order(dist_pos, dist_len, -dv)
  best_j <- idx[ord[1]]
  
  # ---------- escribir SMALL -> BIG (en la fila best_j) ----------
  df_big$matched_from_small[best_j] <- i
  df_big$match_dist_pos[best_j] <- abs(df_big[[col_pos_big]][best_j] - pos_i)
  
  if (use_svlen) {
    df_big$match_dist_len[best_j] <- abs(df_big[[col_len_big]][best_j] - len_i)
  }
  
  for (f in fields_small_to_big) {
    df_big[[paste0("small_", f)]][best_j] <- df_small[[f]][i]
  }
}


df_big_process <- separate(df_big,small_HP, 
                      into = c("HP", "location", "readsSupp", "readsPass", "Filter1", "Filter2"), 
                      sep = ",")

df_big_process <- df_big_process %>%
  mutate(HP_final = case_when(
    grepl("NULL",Filter1)|grepl("NULL", Filter2) ~ "NULL", # <--- Aquí asigna el valor de otra casilla
    grepl("1/1",small_GT)| grepl("1\\|1",small_GT)~ "BOTH",
    TRUE              ~ HP
  ),
  
  GT_final = case_when(
    
    grepl("0/0",small_GT)| grepl("0\\|0",small_GT)~ "Hom-Ref",
    grepl("0/1",small_GT)| grepl("1\\/0",small_GT)|grepl("1\\|0",small_GT)|grepl("0\\|1",small_GT)~ "Het",
    grepl("1/1",small_GT)| grepl("1\\|1",small_GT)~ "Hom",
    grepl("./.",small_GT)~"No-Call",
    TRUE~small_GT
  )
  )

"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                                      PASO FINAL GENERAR EL REPORT FINAL"





ReportALL<-df_big_process%>%
  mutate(
    InsertPoint=posicion,
    Insert_HP=HP_final,
    Genotype=GT_final,
    Class=clase,
    Subclass=subclase,
    Genotype=GT_final,
    MetDiff_L=diff_L,
    LeftHp1=alelo1_Left,
    LeftHp2=alelo2_Left,
    ReadsHp1L=readsF1_L,
    ReadsHp2L=readsF2_L,
    MetDiff_R=diff_R,
    RightHp1=alelo1_Right,
    RightHp2=alelo2_Right,
    ReadsHp1R=readsF1_R,
    ReadsHp2R=readsF2_R,
    SvLen=Svlen,
    Id=id
  )%>%
  select(
    Id,
    InsertPoint,
    SvLen,
    Class,
    Subclass,
    Genotype,
    Insert_HP,
    MetDiff_L,
    LeftHp1,
    LeftHp2,
    ReadsHp1L,
    ReadsHp2L,
    MetDiff_R,
    RightHp1,
    RightHp2,
    ReadsHp1R,
    ReadsHp2R
  )


ReportALL$MetDiff_L<- round(ReportALL$MetDiff_L, digits = 2)
ReportALL$MetDiff_R<- round(ReportALL$MetDiff_R, digits = 2)
ReportALL$LeftHp1<- round(ReportALL$LeftHp1, digits = 2)
ReportALL$LeftHp2<- round(ReportALL$LeftHp2, digits = 2)
ReportALL$RightHp1<- round(ReportALL$RightHp1, digits = 2)
ReportALL$RightHp2<- round(ReportALL$RightHp2, digits = 2)


saveRDS(ReportALL, file = snakemake@output[["FinalTable"]])


#Lo mismo para la tabla completa

TOTAL_INS<-readRDS(snakemake@input[["TOTAL_INS"]])

df_big<-TOTAL_INS
df_small<-filtradoBcftools2
# ---------------- PREP: crear columnas en BIG ----------------


for (f in fields_small_to_big) {
  newcol <- paste0("small_", f)
  if (!newcol %in% names(df_big)) df_big[[newcol]] <- NA
}

# opcional: guardar también info del match
if (!"matched_from_small" %in% names(df_big)) df_big$matched_from_small <- NA_integer_
if (!"match_dist_pos" %in% names(df_big)) df_big$match_dist_pos <- NA_integer_
if (use_svlen && !"match_dist_len" %in% names(df_big)) df_big$match_dist_len <- NA_integer_

# ---------------- LOOP: small -> best big ----------------
for (i in seq_len(nrow(df_small))) {
  
  chr_i <- df_small[[col_chr_small]][i]
  pos_i <- df_small[[col_pos_small]][i]
  
  # candidatos en BIG por chr + ventana
  idx <- which(
    df_big[[col_chr_big]] == chr_i &
      abs(df_big[[col_pos_big]] - pos_i) <= window_bp
  )
  
  if (length(idx) == 0) next
  
  # scoring
  dist_pos <- abs(df_big[[col_pos_big]][idx] - pos_i)
  
  if (use_svlen) {
    len_i <- df_small[[col_len_small]][i]
    dist_len <- abs(df_big[[col_len_big]][idx] - len_i)
  } else {
    dist_len <- rep(0, length(idx))
  }
  
  # (opcional) si quieres usar DV/QUAL del BIG para desempatar:
  dv <- if ("DV" %in% names(df_big)) df_big$DV[idx] else rep(0, length(idx))
  dv[is.na(dv)] <- 0
  
  # elegir best: menor dist_pos, luego menor dist_len, luego mayor DV
  ord <- order(dist_pos, dist_len, -dv)
  best_j <- idx[ord[1]]
  
  # ---------- escribir SMALL -> BIG (en la fila best_j) ----------
  df_big$matched_from_small[best_j] <- i
  df_big$match_dist_pos[best_j] <- abs(df_big[[col_pos_big]][best_j] - pos_i)
  
  if (use_svlen) {
    df_big$match_dist_len[best_j] <- abs(df_big[[col_len_big]][best_j] - len_i)
  }
  
  for (f in fields_small_to_big) {
    df_big[[paste0("small_", f)]][best_j] <- df_small[[f]][i]
  }
}


df_big_process <- separate(df_big,small_HP, 
                      into = c("HP", "location", "readsSupp", "readsPass", "Filter1", "Filter2"), 
                      sep = ",")

df_big_process <- df_big_process %>%
  mutate(HP_final = case_when(
    grepl("NULL",Filter1)|grepl("NULL", Filter2) ~ "NULL", # <--- Aquí asigna el valor de otra casilla
    grepl("1/1",small_GT)| grepl("1\\|1",small_GT)~ "BOTH",
    TRUE              ~ HP
  ),
  
  GT_final = case_when(
    
    grepl("0/0",small_GT)| grepl("0\\|0",small_GT)~ "Hom-Ref",
    grepl("0/1",small_GT)| grepl("1\\/0",small_GT)|grepl("1\\|0",small_GT)|grepl("0\\|1",small_GT)~ "Het",
    grepl("1/1",small_GT)| grepl("1\\|1",small_GT)~ "Hom",
    grepl("./.",small_GT)~"No-Call",
    TRUE~small_GT
  )
  )

"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                                      PASO FINAL GENERAR EL REPORT FINAL"





ReportALL<-df_big_process%>%
  mutate(
    InsertPoint=posicion,
    Insert_HP=HP_final,
    Genotype=GT_final,
    Class=clase,
    Subclass=subclase,
    Genotype=GT_final,
    MetDiff_L=diff_L,
    LeftHp1=alelo1_Left,
    LeftHp2=alelo2_Left,
    ReadsHp1L=readsF1_L,
    ReadsHp2L=readsF2_L,
    MetDiff_R=diff_R,
    RightHp1=alelo1_Right,
    RightHp2=alelo2_Right,
    ReadsHp1R=readsF1_R,
    ReadsHp2R=readsF2_R,
    SvLen=Svlen,
    Id=id
  )%>%
  select(
    Id,
    InsertPoint,
    SvLen,
    Class,
    Subclass,
    Genotype,
    Insert_HP,
    MetDiff_L,
    LeftHp1,
    LeftHp2,
    ReadsHp1L,
    ReadsHp2L,
    MetDiff_R,
    RightHp1,
    RightHp2,
    ReadsHp1R,
    ReadsHp2R
  )


ReportALL$MetDiff_L<- round(ReportALL$MetDiff_L, digits = 2)
ReportALL$MetDiff_R<- round(ReportALL$MetDiff_R, digits = 2)
ReportALL$LeftHp1<- round(ReportALL$LeftHp1, digits = 2)
ReportALL$LeftHp2<- round(ReportALL$LeftHp2, digits = 2)
ReportALL$RightHp1<- round(ReportALL$RightHp1, digits = 2)
ReportALL$RightHp2<- round(ReportALL$RightHp2, digits = 2)

saveRDS(ReportALL, file = snakemake@output[["FinalTable2"]])
