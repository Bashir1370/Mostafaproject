#!/usr/bin/env Rscript
# Build and freeze the mouse gene-set layer used by BOTH GSE125002 and GSE286387.
# This is an explicit audit/update step. Main inferential scripts never call live MSigDB.

required <- c("msigdbr","babelgene","dplyr","readr","tibble")
miss <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]
if(length(miss)) stop("Missing packages: ", paste(miss, collapse=", "))
suppressPackageStartupMessages({library(dplyr);library(readr);library(tibble)})

project_dir <- "cipn-ferroptosis-screening"
config_dir <- file.path(project_dir,"config")
mouse_dir <- file.path(config_dir,"mouse")
dir.create(mouse_dir, recursive=TRUE, showWarnings=FALSE)

vinik_file <- file.path(config_dir,"vinik_2024_24_biomarkers.csv")
out_vinik <- file.path(mouse_dir,"vinik_2024_24_human_to_mouse_locked.csv")
out_vinik_audit <- file.path(mouse_dir,"vinik_2024_24_human_to_mouse_audit.csv")
out_msig <- file.path(mouse_dir,"msigdb_Hs_to_mouse_gene_sets_locked.csv")
out_manifest <- file.path(mouse_dir,"gene_set_manifest_mouse.csv")

if(!file.exists(out_vinik)) stop("Missing manually reviewed Vinik mouse lock: ", out_vinik)
if(any(file.exists(c(out_msig,out_manifest))) && Sys.getenv("FORCE_REBUILD_MOUSE_GENESETS") != "1") {
  stop("Mouse MSigDB locked files already exist. Refusing silent rebuild. Set FORCE_REBUILD_MOUSE_GENESETS=1 only for an explicit audited update.")
}

vinik <- read_csv(vinik_file, show_col_types=FALSE)
vinik_locked <- read_csv(out_vinik, show_col_types=FALSE)
if(nrow(vinik_locked)!=24 || !setequal(vinik$human_symbol,vinik_locked$human_symbol)) stop("Locked Vinik mouse file is not the same 24-gene source panel.")
ortho <- babelgene::orthologs(genes=unique(vinik$human_symbol), species="mouse", human=TRUE, min_support=2, top=TRUE) |> as_tibble()
human_col <- intersect(c("human_symbol","human_gene","human"),names(ortho))[1]; mouse_col <- intersect(c("symbol","mouse_symbol","ortholog_symbol"),names(ortho))[1]
if(is.na(human_col)||is.na(mouse_col)) stop("Unexpected babelgene columns: ",paste(names(ortho),collapse=", "))
live_map <- ortho |> transmute(human_symbol=.data[[human_col]], live_mouse_symbol=.data[[mouse_col]], support=if("support"%in%names(ortho)) support else NA_character_)
vinik_map <- vinik_locked |> left_join(live_map,by="human_symbol") |> mutate(matches_live_mapping=mouse_symbol==live_mouse_symbol)
write_csv(vinik_map,out_vinik_audit)
if(any(is.na(vinik_map$live_mouse_symbol))) stop("Live ortholog audit failed for one or more Vinik genes.")

# Keep source definition identical to the rat analysis: Human MSigDB -> target-species orthologs.
msig_mouse <- msigdbr::msigdbr(db_species="HS", species="Mus musculus")
required_sets <- c("GOBP_FERROPTOSIS","WP_FERROPTOSIS","HALLMARK_APOPTOSIS","HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY")
locked <- msig_mouse |> filter(gs_name %in% required_sets) |> transmute(signature=gs_name, mouse_gene_symbol=gene_symbol, gs_id, gs_source_species, db_version, num_ortholog_sources=if("num_ortholog_sources"%in%names(msig_mouse)) num_ortholog_sources else NA_integer_) |> distinct()
if(!all(required_sets %in% locked$signature)) stop("One or more required MSigDB sets are unavailable.")
if(any(table(locked$signature)<5)) stop("Unexpectedly small mouse gene set after mapping.")
write_csv(locked,out_msig)

human <- msigdbr::msigdbr(db_species="HS", species="Homo sapiens") |> filter(gs_name %in% required_sets)
manifest <- bind_rows(
  tibble(signature="VINIK_2024_24_FERROPTOSIS_BIOMARKERS",role="primary",source="Vinik et al. Advanced Science 2024",source_identifier="DOI 10.1002/advs.202307263; PMID 38441406",source_species="Homo sapiens",analysis_species="Mus musculus",source_member_count=24L,locked_mouse_member_count=n_distinct(vinik_locked$mouse_symbol),mapping="babelgene human->mouse; min_support=2; top=TRUE",version="Vinik 2024",interpretation_note="Human ferroptosis-vs-apoptosis biomarker panel; ortholog-mapped to mouse, not mouse-validated."),
  tibble(signature=required_sets,role=c("supporting","supporting","context_control","context_control"),source=c("MSigDB C5 GO Biological Process","MSigDB C2 CP WikiPathways","MSigDB Hallmark","MSigDB Hallmark"),source_identifier=c("GO:0097707","WP4313","HALLMARK_APOPTOSIS","HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY"),source_species="Homo sapiens",analysis_species="Mus musculus",source_member_count=vapply(required_sets,function(x)n_distinct(human$gene_symbol[human$gs_name==x]),integer(1)),locked_mouse_member_count=vapply(required_sets,function(x)n_distinct(locked$mouse_gene_symbol[locked$signature==x]),integer(1)),mapping="msigdbr Human MSigDB -> Mus musculus orthologs; frozen after audit",version=paste(unique(locked$db_version),collapse=";"),interpretation_note=c("Process-membership set; includes promoting and limiting genes; NES is not a direct activation score.","Mechanistic pathway membership; supporting evidence only.","Apoptosis context/specificity control; not ferroptosis evidence.","ROS context control; not ferroptosis evidence."))
)
write_csv(manifest,out_manifest)
message("Mouse gene-set audit locked successfully.")
print(manifest |> dplyr::select(signature, source_member_count, locked_mouse_member_count, version))
