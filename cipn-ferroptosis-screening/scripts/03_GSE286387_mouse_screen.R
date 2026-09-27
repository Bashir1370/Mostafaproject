#!/usr/bin/env Rscript
# GSE286387 -- Mus musculus DRG -- chronic Oxaliplatin vs Vehicle
# C57BL/6 male mice; 10 mg/kg i.p. weekly x8 weeks; bulk bilateral lumbar DRG RNA-seq; n=5/group in RNA-seq.
# This loader deliberately accepts ONLY a genome-scale raw-count matrix. It refuses normalized-only or DEG-only tables.

source("cipn-ferroptosis-screening/scripts/lib_mouse_ferroptosis.R")
if(!requireNamespace("GEOquery",quietly=TRUE))stop("Install GEOquery before running GSE286387.")
accession<-"GSE286387";project_dir<-"cipn-ferroptosis-screening";raw_dir<-file.path(project_dir,"data","raw",accession);result_dir<-file.path(project_dir,"results",accession);dir.create(raw_dir,recursive=TRUE,showWarnings=FALSE)
message("[GSE286387] Fetching GEO metadata and supplementary files...")
gse<-GEOquery::getGEO(accession,GSEMatrix=TRUE,getGPL=FALSE);eset<-if(is.list(gse))gse[[1]] else gse;pd<-Biobase::pData(eset);write.csv(pd,file.path(raw_dir,"GEO_sample_metadata_audit.csv"),row.names=TRUE)
GEOquery::getGEOSuppFiles(accession,makeDirectory=TRUE,baseDir=raw_dir,fetch_files=TRUE)
files<-list.files(raw_dir,recursive=TRUE,full.names=TRUE);files<-files[grepl("\\.(txt|tsv|csv|gz)$",files,ignore.case=TRUE)]
read_candidate<-function(f){
  tryCatch({x<-data.table::fread(f,data.table=FALSE,check.names=FALSE);if(nrow(x)<5000||ncol(x)<11)return(NULL);num<-vapply(x,is.numeric,logical(1));if(sum(num)<10)return(NULL);gene_col<-which(!num)[1];if(is.na(gene_col))gene_col<-1;sm<-which(num);m<-as.matrix(x[,sm,drop=FALSE]);fin<-mean(is.finite(m));nonneg<-mean(m>=0,na.rm=TRUE);intfrac<-mean(abs(m-round(m))<1e-8,na.rm=TRUE);if(fin<.99||nonneg<.99||intfrac<.95)return(NULL);list(file=f,x=x,gene_col=gene_col,sample_cols=sm,score=nrow(x)*ncol(m))},error=function(e)NULL)
}
cands<-Filter(Negate(is.null),lapply(files,read_candidate));if(!length(cands))stop("No genome-scale raw-count matrix was found in GSE286387 supplementary files. The script refuses to run DESeq2/GSEA on normalized or DEG-only data. Inspect data/raw/GSE286387 and GEO_sample_metadata_audit.csv; if GEO exposes only FASTQ, quantify raw reads first.")
best<-cands[[which.max(vapply(cands,`[[`,numeric(1),"score"))]];x<-best$x;ids<-x[[best$gene_col]];m<-as.matrix(x[,best$sample_cols,drop=FALSE]);storage.mode(m)<-"numeric";rownames(m)<-standardize_mouse_ids(ids);colnames(m)<-names(x)[best$sample_cols]
message("Selected raw-count source: ",best$file," [",nrow(m)," genes x ",ncol(m)," numeric columns]")
# Match count columns to GEO samples using GSM IDs first, then titles.
gsms<-rownames(pd);titles<-as.character(pd$title);map_one<-function(cn){hit<-which(vapply(gsms,function(g)grepl(g,cn,fixed=TRUE),logical(1)));if(length(hit)==1)return(hit);hit<-which(vapply(titles,function(t)grepl(make.names(t),make.names(cn),fixed=TRUE)||grepl(make.names(cn),make.names(t),fixed=TRUE)),logical(1)));if(length(hit)==1)hit else NA_integer_}
idx<-vapply(colnames(m),map_one,integer(1));
# Derive treatment from all metadata text; keep only unambiguous control and 10 mg/kg oxaliplatin RNA-seq samples.
alltxt<-apply(pd,1,function(z)paste(z,collapse=" | "));grp<-rep(NA_character_,nrow(pd));grp[grepl("control|vehicle|dextrose|0 mg",alltxt,ignore.case=TRUE)]<-"Vehicle";grp[grepl("oxaliplatin|10 mg/kg|10mg/kg",alltxt,ignore.case=TRUE)]<-"Oxaliplatin"
if(sum(!is.na(idx))>=10){keep<-!is.na(idx);m<-m[,keep,drop=FALSE];idx<-idx[keep];meta<-tibble(sample=colnames(m),gsm=gsms[idx],geo_title=titles[idx],group=grp[idx])} else {
  # fallback for matrices labeled by simple sample names: require exactly 10 columns and infer group from column labels.
  if(ncol(m)!=10)stop("Could not map raw-count columns to GEO samples unambiguously.");cg<-rep(NA_character_,10);cg[grepl("control|vehicle|ctrl|con",colnames(m),ignore.case=TRUE)]<-"Vehicle";cg[grepl("oxa|oxaliplatin|10mg",colnames(m),ignore.case=TRUE)]<-"Oxaliplatin";meta<-tibble(sample=colnames(m),gsm=NA_character_,geo_title=colnames(m),group=cg)
}
if(sum(meta$group=="Vehicle",na.rm=TRUE)!=5||sum(meta$group=="Oxaliplatin",na.rm=TRUE)!=5)stop("GSE286387 design audit failed: expected 5 Vehicle + 5 Oxaliplatin RNA-seq samples. Review GEO_sample_metadata_audit.csv before proceeding.")
meta<-meta |> filter(!is.na(group)) |> mutate(replicate=ave(seq_along(group),group,FUN=seq_along));m<-m[,meta$sample,drop=FALSE];write_csv(meta,file.path(project_dir,"config","GSE286387_samples_runtime_audited.csv"))
run_mouse_ferroptosis_screen(m,meta,accession,result_dir,dataset_note="Chronic C57BL/6 male mouse OIPN; bilateral lumbar DRG; 10 mg/kg oxaliplatin i.p. weekly for 8 weeks vs vehicle; n=5/group RNA-seq. Raw-count source and GEO sample mapping are audited at runtime and exported.")
