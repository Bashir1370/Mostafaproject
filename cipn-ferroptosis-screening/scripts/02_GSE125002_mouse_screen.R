#!/usr/bin/env Rscript
# GSE125002 -- Mus musculus DRG -- Oxaliplatin vs Vehicle
# Starobova/Mueller/Vetter; L3-L5 DRG; C57BL/6; 3 biological replicates/group.
# Each biological replicate pooled DRG from four animals. Cisplatin samples are excluded.

source("cipn-ferroptosis-screening/scripts/lib_mouse_ferroptosis.R")
accession<-"GSE125002"; project_dir<-"cipn-ferroptosis-screening"; raw_dir<-file.path(project_dir,"data","raw",accession); result_dir<-file.path(project_dir,"results",accession);dir.create(raw_dir,recursive=TRUE,showWarnings=FALSE)
meta<-readr::read_csv(file.path(project_dir,"config","GSE125002_samples.csv"),show_col_types=FALSE)
url<-"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE125nnn/GSE125002/suppl/GSE125002_RAW.tar"; tarfile<-file.path(raw_dir,"GSE125002_RAW.tar")
if(!file.exists(tarfile)||file.info(tarfile)$size==0){options(timeout=max(1200,getOption("timeout")));download.file(url,tarfile,mode="wb")}
extdir<-file.path(raw_dir,"extracted");dir.create(extdir,showWarnings=FALSE);if(!length(list.files(extdir)))untar(tarfile,exdir=extdir)
read_htseq<-function(gsm){f<-list.files(extdir,pattern=gsm,full.names=TRUE,ignore.case=TRUE);if(length(f)!=1)stop("Expected one processed-count file for ",gsm,"; found ",length(f));x<-data.table::fread(f,header=FALSE,data.table=FALSE);if(ncol(x)<2)stop("Malformed HTSeq file: ",f);ids<-as.character(x[[1]]);vals<-suppressWarnings(as.numeric(x[[ncol(x)]]));keep<-!is.na(vals)&is.finite(vals)&vals>=0&!grepl("^__",ids);sym<-standardize_mouse_ids(ids[keep]);tibble(gene=sym,count=round(vals[keep])) |> filter(!is.na(gene),nzchar(gene)) |> group_by(gene) |> summarise(count=sum(count),.groups="drop") |> rename(!!gsm:=count)}
tabs<-lapply(meta$sample,read_htseq);ct<-Reduce(function(a,b)full_join(a,b,by="gene"),tabs) |> mutate(across(-gene,~replace_na(.x,0)));mat<-as.matrix(ct[,-1]);rownames(mat)<-ct$gene
run_mouse_ferroptosis_screen(mat,meta,accession,result_dir,dataset_note="C57BL/6 mouse L3-L5 whole DRG; n=3/group; each biological replicate pools four animals; oxaliplatin cumulative 40 ug intraplantar. Cisplatin samples excluded.")
