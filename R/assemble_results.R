#!/usr/bin/env Rscript

#The overlap script is associated with "COMPREHENSIVE DISCOVERY OF CLINICALLY 
#RELEVANT FUSIONS IN PEDIATRIC CANCER" By S LaHaye, J Fitch, K Voytovich, et al.
#Developed at the Institute for Genomic Medicine at Nationwide Children's Hospital
#insert full citation

################################################################################ 
#Compatible with the following fusion detection algorithms: fusionMap, fusionCatcher, 
#jaffa, mapSplice, soapFuse, starFusion, tophatFusion, DRAGEN
################################################################################ 

srcdir <- dirname(sub("--file=", "", grep("--file=", commandArgs(), value = T)))
suppressWarnings({
  source(file.path(srcdir, "utils_import.R"))
  source(file.path(srcdir, "utils_equalize.R"))
  source(file.path(srcdir, "utils_equalize_helpers.R"))
  source(file.path(srcdir, "utils_add_ranks.R"))
  source(file.path(srcdir, "utils_results.R"))
  source(file.path(srcdir, "connect_to_db.R"))
  
  library(optparse, quietly = TRUE, warn.conflicts = FALSE)
  library(dplyr, quietly = TRUE, warn.conflicts = FALSE)
  library(readr, quietly = TRUE, warn.conflicts = FALSE)
  
  library(DBI, quietly = TRUE, warn.conflicts = FALSE)
  library(dbplyr, quietly = TRUE, warn.conflicts = FALSE)
  library(magrittr, quietly = TRUE, warn.conflicts = FALSE)
  library(purrr, quietly = TRUE, warn.conflicts = FALSE)
})

################################################################################ 

# Process Input

usage <- "assemble_results.R --sample s1"
description <- paste0("assemble_results.R reads a set of fusion detection result files and ", 
  "produces an aggregated report with the overlapping fusion events predicted ", 
  "in the input files. The program searches recursively from the current ", "directory to find known files. The input files can also be specified ", 
  "manually.")
option_list <- list(make_option("--sample", type = "character", help = "Name of the sample. (required)"), 
  make_option("--baseDir", type = "character", help = "Base directory to search for input files. [default = .]"), 
  make_option("--outReport", type = "character", help = "Location of the output report. [default = overlap_$sample.tsv]"), 
  make_option("--outReportfilt", type = "character", help = "Location of the output report. [default = overlap_filter_$sample.tsv]"), 
  make_option("--foutReport", type = "character", help = "Location of the output report. [default = filtered_overlap_whitelist_$sample.tsv]"), 
  make_option("--foutReport3", type = "character", help = "Location of the output report. [default = filtered_overlap_whitelist_3callers_$sample.tsv]"), 
  make_option("--outBreakpoints", type = "character", help = "Location of the breakpoint report. [default = breakpoints_$sample.tsv]"), 
  make_option("--dragen", type = "character", default = NULL, help = "Path to the dragen results file. [default = search for file named like 'DRAGEN.fusion_candidates.final']"),
  make_option("--fusionCatcher", type = "character", default = NULL, help = "Path to the fusionCatcher results file. [default = search for file named like 'final-list_candidate-fusion-genes.txt']"), 
  make_option("--fusionMap", type = "character", default = NULL, help = "Path to the fusionMap results file. [default = search for file named like 'FusionDetection.FusionReport.Table.txt']"), 
  make_option("--jaffa", type = "character", default = NULL, help = "Path to the jaffa results file. [default = search for file named like 'jaffa_results.csv']"), 
  make_option("--mapSplice", type = "character", default = NULL, help = "Path to the mapSplice results file. [default = search for file named like 'fusions_well_annotated.txt']"), 
  make_option("--soapFuse", type = "character", default = NULL, help = "Path to the soapFuse results file. [default = search for file named like '.final.Fusion.specific.for.genes']"), 
  make_option("--starFusion", type = "character", default = NULL, help = "Path to the starFusion results file. [default = search for file named like 'star-fusion.fusion_predictions.abridged.tsv']"), 
  make_option("--tophatFusion", type = "character", default = NULL, help = "Path to the tophatFusion results file. [default = search for file named like 'result.txt']. potential_fusion.txt must be in the same directory as result.txt for the import to work.")) 
opt <- parse_args(OptionParser(usage = usage, description = description, option_list = option_list))

stopifnot(is.character(opt$sample))
sample <- opt$sample

if (is.null(opt$baseDir)) {
  baseDir <- "."
} else {
  baseDir <- opt$baseDir
}
if (is.null(opt$outReport)) {
  outReport <- paste0("overlap_", sample, ".tsv")
} else {
  outReport <- opt$outReport
}
if (is.null(opt$outReportfilt)) {
  outReportfilt <- paste0("overlap_filtered__", sample, ".tsv")
}
if (is.null(opt$foutReport)) {
  foutReport <- paste0("filtered_overlap_whitelist_", sample, ".tsv")
}
if (is.null(opt$foutReport3)) {
  foutReport3 <- paste0("filtered_overlap_whitelist_3callers_", sample, ".tsv")
}
if (is.null(opt$outBreakpoints)) {
  outBreakpoints <- paste0("breakpoints_", sample, ".tsv")
}

# Assemble a list of files to merge

cat(paste0("Assembing a list of files to merge...", "\n"))

find_results <- function(base) {
  results_patterns <- list(starFusion = "star-fusion\\.fusion_predictions\\.abridged\\.tsv", 
    fusionMap = "FusionDetection\\.FusionReport\\.Table\\.txt", tophatFusion = "result\\.txt", 
    fusionCatcher = "final-list_candidate-fusion-genes\\.txt", jaffa = "jaffa_results\\.csv", 
    mapSplice = "fusions_well_annotated\\.txt", soapFuse = "\\.final\\.Fusion\\.specific\\.for\\.genes", 
    dragen = "DRAGEN\\.fusion_candidates\\.final")
  map(results_patterns, function(p) {
    list.files(base, pattern = p, recursive = TRUE, full.names = TRUE)
  })
}

list_results <- find_results(baseDir)

tools <- c("fusionCatcher", "fusionMap", "jaffa", "mapSplice", "soapFuse", 
  "starFusion", "tophatFusion", "dragen")
custom_tools <- match(tools, names(opt))
custom_tools <- custom_tools[!is.na(custom_tools)]
custom_input <- opt[custom_tools]

input_list <- list_results
for (tool in tools) {
  if (tool %in% names(custom_input)) {
    input_list[[tool]] <- custom_input[[tool]]
  }
  if (tool %in% names(input_list)) {
    if (length(input_list[[tool]]) > 0) {
      input_list[[tool]] <- input_list[[tool]][1]
    } else {
      input_list[[tool]] <- NULL
    }
  }
}

print(names(input_list))
cat(paste0("\n"))

# Put together all_data list

all_data <- vector("list", length(input_list))
size <- 0
for (tool in names(input_list)) {
  size <- size + 1
  all_data[[size]] <- list(tool = tool, file = input_list[[tool]])
}

if (size == 0) {
  stop("At least 1 tool must be provided as input.")
}


################################################################################ 

# Import datasets defined in all_data

calls <- vector("list", length(all_data))
for (i in length(all_data):1) {
  new_df <- import_results(all_data[[i]]$tool, all_data[[i]]$file)
  calls[[i]] <- list(tool = all_data[[i]]$tool, calls = nrow(new_df))
  if (all(dim(new_df) == c(0, 0))) {
    all_data[[i]] <- NULL
  } else {
    all_data[[i]]$df <- new_df
    all_data[[i]]$df <- equalize_results(all_data[[i]]$tool, all_data[[i]]$df)
    all_data[[i]]$df <- add_ranks(all_data[[i]]$tool, all_data[[i]]$df)
  }
}

if (length(all_data) == 0) {
  stop("At least 1 of the tools must have non-empty results.")
}

# Get the fusions

fusions <- list()
for (i in 1:length(all_data)) {
  tool <- all_data[[i]]$tool
  fusions[[tool]] <- unique(all_data[[i]]$df$UnorderedFusion)
}

################################################################################ 

# Get the 'population frequency' for all previously reported fusion pairs

con <- connect_to_db()

fusion_db <- tbl(con, "fusion")
analysis_db <- tbl(con, "analysis")

ordered_db <- left_join(fusion_db, analysis_db, "analysis_id") %>% distinct(sample_id, 
  gene1, gene2) %>% collect()

unordered_db <- ordered_db %>% mutate(unordered = get_unordered(gene1, gene2)) %>% 
  distinct(sample_id, unordered)

################################################################################ 

# Get the overlapping genes and results

overlapping <- get_overlapping_from_fusions(fusions)
report <- get_report(overlapping, all_data)
report <- order_report(report)
all_breakpoints <- get_report_all_breakpoints(overlapping, all_data)
all_breakpoints <- order_report_all_breakpoints(all_breakpoints, report)

################################################################################ 

# Add population frequencies to the report

pop_frequency <- function(unordered_fusions, unordered_db) {
  result <- purrr::map_chr(unordered_fusions, function(uf) {
    return(sum(unordered_db$unordered == uf)/length(unique(unordered_db$sample_id)))
  })
  return(result)
}

pop_count <- function(unordered_fusions, unordered_db) {
  result <- purrr::map_chr(unordered_fusions, function(uf) {
    return(sum(unordered_db$unordered == uf))
  })
}

report <- report %>% mutate(GenePairFrequency = pop_frequency(UnorderedFusion, unordered_db), 
  GenePairCount = pop_count(UnorderedFusion, unordered_db), SampleCount = length(unique(unordered_db$sample_id))) %>% 
  select(UnorderedFusion, NumTools, GenePairFrequency, GenePairCount, SampleCount, 
    everything())

################################################################################ 

# Create filtered report Create filtered report
whitelist <- c("ABI1+KMT2A", "ABL1+BCR", "ABL1+ETV6", "ABL1+NUP214", "ACE+KPNB1", 
  "ACSL3+ETV1", "ACTB+GLI1", "ACTN4+KMT2A", "ADGRG7+TFG", "AFDN+KMT2A", "AFF1+KMT2A", 
  "AFF3+KMT2A", "AFF4+KMT2A", "AGGF1+RAF1", "AGK+BRAF", "AKAP13+RET", "AKAP9+BRAF", 
  "AKT3+ARHGAP30", "AKT3+CEP170", "AKT3+PLD5", "AKT3+ZEB2", "ALK+ATIC", "ALK+CARS", 
  "ALK+CLTC", "ALK+DCTN1", "ALK+DDX6", "ALK+EML4", "ALK+GTF2IRD1", "ALK+HIP1", 
  "ALK+KIF5B", "ALK+KLC1", "ALK+MALAT1", "ALK+MSN", "ALK+NPM1", "ALK+PPFIBP1", 
  "ALK+RANBP2", "ALK+STRN", "ALK+TFG", "ALK+TPM3", "ALK+TPM4", "ALK+VCL", "AP3B1+BRAF", 
  "ARHGAP26+CLDN18", "ARHGAP26+KMT2A", "ARHGAP6+CLDN18", "ARID1A+MAST2", "ASPSCR1+TFE3", 
  "ATF1+EWSR1", "ATF1+FUS", "ATG7+BRAF", "ATP1B4+CHST11", "ATP5L+KMT2A", "BAIAP2L1+FGFR3", 
  "BCL2L11+BRAF", "BCOR+ZC3H7B", "BCR+JAK2", "BICC1+FGFR2", "BIRC6+LTBP1", "BRAF+CCNY", 
  "BRAF+CEP89", "BRAF+CLCN6", "BRAF+ERC1", "BRAF+FAM114A2", "BRAF+FAM131B", "BRAF+FCHSD1", 
  "BRAF+GATM", "BRAF+GNAI1", "BRAF+HERPUD1", "BRAF+KIAA1549", "BRAF+LSM14A", "BRAF+MACF1", 
  "BRAF+MKRN1", "BRAF+RNF130", "BRAF+RUNDC1", "BRAF+SLC45A3", "BRAF+SND1", "BRAF+SVOPL", 
  "BRAF+TAX1BP1", "BRAF+TRIM24", "BRAF+ZC3HAV1", "BRAF+ZSCAN30", "BRD3+NUTM1", 
  "BRD4+NUTM1", "C11orf1+SIK3", "C2CD2L+ZNF585B", "C8orf34+MET", "CADM2+MITF", 
  "CADM2+TFEB", "CANT1+ETV4", "CAPZA2+MET", "CBFA2T3+GLIS2", "CBFB+MYH11", "CCDC186+FGFR2", 
  "CCDC6+RET", "CCNB1IP1+HMGA2", "CD74+NRG1", "CD74+ROS1", "CDH11+USP6", "CDKN2D+WDFY2", 
  "CDX1+IRF2BP2", "CENPP+WNK2", "CEP170B+KMT2A", "CHCHD7+PLAG1", "CIC+DUX4L1", 
  "CIC+FOXO4", "CIPC+NGFR", "CLCN6+RAF1", "CLIP1+ROS1", "CLTC+ROS1", "CLTC+TFE3", 
  "CNTNAP2+GRM8", "COBL+SEPT14", "COL1A1+PDGFB", "COL1A1+USP6", "COL1A2+PLAG1", 
  "COL21A1+TFEB", "COX6C+HMGA2", "CPSF4L+ERBB4", "CREB1+EWSR1", "CREB3L1+FUS", 
  "CREB3L2+FUS", "CREBBP+KAT6A", "CREBBP+KMT2A", "CRTC1+MAML2", "CRTC3+MAML2", 
  "CTNNB1+PLAG1", "CXorf67+MBTD1", "DDIT3+EWSR1", "DDIT3+FUS", "DDX20+TBX15", "DHH+RHEBL1", 
  "DHX33+NLRP1", "DNAJB1+PRKACA", "DPM1+GRID1", "DVL2+TFE3", "EBF1+HMGA2", "EEFSEC+KMT2A", 
  "EGFR+SEC61G", "EGFR+SEPT14", "EIF3E+RSPO2", "ELAVL3+FGFR3", "ELK4+SLC45A3", 
  "ELL+KMT2A", "EP300+KMT2A", "EPHA3+LCLAT1", "EPHB1+MOBKL1B", "EPS15+KMT2A", "ERBB2+MTSS1", 
  "ERBB4+EZR", "ERC1+RET", "ERC1+ROS1", "ERG+EWSR1", "ERG+FUS", "ERG+NDRG1", "ERG+SLC45A3", 
  "ERG+TMPRSS2", "ESRP1+RAF1", "ETV1+EWSR1", "ETV1+HNRNPA2B1", "ETV1+KLK2", "ETV1+SLC45A3", 
  "ETV1+TMPRSS2", "ETV4+EWSR1", "ETV4+TMPRSS2", "ETV6+JAK2", "ETV6+MN1", "ETV6+NTRK3", 
  "ETV6+PDGFRB", "ETV6+RUNX1", "EWSR1+FEV", "EWSR1+FLI1", "EWSR1+MYB", "EWSR1+NFATC2", 
  "EWSR1+NR4A3", "EWSR1+PBX1", "EWSR1+POU5F1", "EWSR1+WT1", "EWSR1+YY1", "EWSR1+ZNF384", 
  "EWSR1+ZNF444", "EXOSC10+MTOR", "EZR+ROS1", "FGFR1+PLAG1", "FGFR1+TACC1", "FGFR2+FRK", 
  "FGFR2+OFD1", "FGFR2+SHTN1", "FGFR2+VCL", "FGFR3+TACC3", "FHIT+HMGA2", "FKBP15+RET", 
  "FLT3LG+RPS11", "FLT4+LBH", "FOXO1+PAX3", "FOXO1+PAX7", "FOXO3+KMT2A", "FOXP1+MITF", 
  "FRMD4B+MITF", "GABBR2+NOTCH1", "GAS7+KMT2A", "GOLGA5+RET", "GOPC+ROS1", "GOSR1+ZNF207", 
  "GPBP1L1+MAST2", "H2AFX+WDR18", "HACL1+RAF1", "HAS2+PLAG1", "HEY1+NCOA2", "HLA-A+ROS1", 
  "HMGA2+LHFPL6", "HMGA2+LPP", "HMGA2+NFIB", "HMGA2+PCBP2", "HMGA2+RAD51B", "HMGA2+SENP1", 
  "HMGA2+TSFM", "HMGA2+WIF1", "HNF1B+NOTCH1", "HOOK3+RET", "IGF2BP3+THADA", "INSL3+JAK3", 
  "IRF2BP2+NTRK1", "JAK2+PAX5", "JAK2+PCM1", "JAK2+SEC31A", "JAK2+SSBP2", "JAZF1+PHF1", 
  "JAZF1+SUZ12", "KDM5A+NUP98", "KIAA1797+p16INK4", "KIF5B+RET", "KMT2A+KNL1", 
  "KMT2A+MLLT1", "KMT2A+MLLT10", "KMT2A+MLLT11", "KMT2A+MLLT3", "KMT2A+MLLT4", 
  "KMT2A+MLLT6", "KMT2A+SEPT2", "KMT2A+SEPT5", "KMT2A+SEPT6", "KMT2A+SEPT9", "KMT2A+TET1", 
  "KSR1+TENM1", "KTN1+RET", "LANCL2+SEPT14", "LGR5+NUP107", "LIFR+PLAG1", "LMNA+NTRK1", 
  "LRIG3+ROS1", "LRRC37B+NF1", "LTK+UACA", "MAML3+TCF4", "MAML3+UBTF", "MAST1+NFIX", 
  "MAST1+NUP210", "MAST1+TADA2A", "MAST1+ZNF700", "MEAF6+PHF1", "MECOM+RPN1", "MECOM+RUNX1", 
  "MET+TFG", "MLLT10+PICALM", "MLLT10+PPP2R1B", "MYB+NFIB", "MYO5A+ROS1", "NAB2+STAT6", 
  "NACC2+NTRK2", "NBR1+WSB1", "NCOA1+PAX3", "NCOA2+PAX3", "NCOA4+RET", "NF1+RAB11FIP4", 
  "NONO+TFE3", "NOTCH1+SEC16A", "NR4A3+TAF15", "NR4A3+TFG", "NRG1+SLC3A2", "NSD1+NUP98", 
  "NTRK1+SQSTM1", "NTRK1+SSBP2", "NTRK1+TFG", "NTRK1+TP53", "NTRK1+TPM3", "NTRK1+TPR", 
  "NTRK2+QKI", "NTRK3+RBPMS", "NUP214+SET", "NUTM2A+YWHAE", "NUTM2B+YWHAE", "OLR1+XIAP", 
  "PAX8+PPARG", "PBX1+TCF3", "PCM1+RET", "PDCD6+TERT", "PLAG1+TCEA1", "PML+RARA", 
  "PPFIBP1+ROS1", "PRCC+TFE3", "PRKAR1A+RET", "PTPRK+RSPO3", "PWWP2A+ROS1", "RAF1+SRGAP3", 
  "RAF1+TRAK1", "RBM10+TFE3", "RELCH+RET", "RET+SPECC1L", "RET+TBL1XR1", "RET+TRIM24", 
  "RET+TRIM27", "RET+TRIM33", "RHOT1+TNKS", "ROS1+SDC4", "ROS1+SHTN1", "ROS1+SLC34A2", 
  "ROS1+TPM3", "ROS1+ZCCHC8", "RPS2P32+THADA", "RUNX1+RUNX1T1", "SFPQ+TFE3", "SLC16A14+TRIP12", 
  "SS18+SSX1", "SS18+SSX2", "SS18+SSX4", "STIL+TAL1", "TBL1XR1+TP63", "TCF7L2+VTI1A", 
  "TERT+TRIO", "TG+THADA")

report$KnownFusion <- ifelse(report$UnorderedFusion %in% whitelist, "yes", "no")  #create column identifying whitelist genes as 'known fusions'

freport <- filter(report, GenePairFrequency < 0.1 | !report$KnownFusion == "no")  # Filter based on database frequency for fusion < 10% and not in whitelist

freport$Distance = ifelse(freport$Chr1 == freport$Chr2 & freport$Strand1 == freport$Strand2, 
  abs((as.numeric(freport$Break1)) - (as.numeric(freport$Break2))), NA)

freport$ReadThrough = ifelse(freport$Distance > 2e+05 | is.na(freport$Distance), 
  "no", "yes")  #create readthrough column
freport <- filter(freport, freport$ReadThrough == "no" | !freport$KnownFusion == 
  "no")  # Filter out read through events unless in whitelist

freport_read <- filter(freport, freport$SupportingReads >= 4 | !freport$KnownFusion == 
  "no")  # At least 1 caller must have 4 reads of evidence or be in whitelist
fusions_with_enough_reads <- unique(freport_read$UnorderedFusion)
freport <- filter(freport, UnorderedFusion %in% fusions_with_enough_reads)  #filter fusions with out enough reads and not in whitelist


freport <- freport[c("UnorderedFusion", "NumTools", "KnownFusion", "GenePairFrequency", 
  "GenePairCount", "SampleCount", "Tool", "Rank", "Total", "SupportingReads", "Gene1", 
  "Chr1", "Break1", "Strand1", "Gene2", "Chr2", "Break2", "Strand2", "Distance", 
  "FMFusionJunctionSequence")]


comment_line <- paste0("# Sample: ", sample)
comment_line <- paste0(comment_line, "\n", "# NumToolsAggregated: ", length(calls))
for (i in 1:length(calls)) {
  comment_line <- paste0(comment_line, "\n", "# - ", calls[[i]]$tool, "Calls = ", 
    calls[[i]]$calls)
}

report_conn <- file(foutReport)
writeLines(comment_line, report_conn)
close(report_conn)
write_tsv(freport, foutReport, append = TRUE, col_names = TRUE)

freport3 <- filter(freport, freport$NumTools > 2 | !freport$KnownFusion == "no")  #filter out fusions with less than 3 callers unless in whitelist

report_conn <- file(foutReport3)
writeLines(comment_line, report_conn)
close(report_conn)
write_tsv(freport3, foutReport3, append = TRUE, col_names = TRUE)


# Create filtered report
Reportfilt <- filter(report, GenePairFrequency < 0.1)  # Database frequency for fusion < 10%.
Reportfilt$Distance = ifelse(Reportfilt$Chr1 == Reportfilt$Chr2 & Reportfilt$Strand1 == 
  Reportfilt$Strand2, abs((as.numeric(Reportfilt$Break1)) - (as.numeric(Reportfilt$Break2))), 
  NA)
Reportfilt <- filter(Reportfilt, Distance > 2e+05 | is.na(Distance))  # Filter out read through events
Reportfilt_read <- filter(Reportfilt, SupportingReads >= 4)  # At least 1 caller must have 4 reads of evidence.
fusions_with_enough_reads <- unique(Reportfilt_read$UnorderedFusion)
Reportfilt <- filter(Reportfilt, UnorderedFusion %in% fusions_with_enough_reads)
Reportfilt <- Reportfilt[c("UnorderedFusion", "NumTools", "GenePairFrequency", "GenePairCount", 
  "SampleCount", "Tool", "Rank", "Total", "SupportingReads", "Gene1", "Chr1", "Break1", 
  "Strand1", "Gene2", "Chr2", "Break2", "Strand2", "Distance", "FMFusionJunctionSequence")]
write_tsv(Reportfilt, outReportfilt)

################################################################################ 

# Write results

comment_line <- paste0("# Sample: ", sample)
comment_line <- paste0(comment_line, "\n", "# NumToolsAggregated: ", length(calls))
for (i in 1:length(calls)) {
  comment_line <- paste0(comment_line, "\n", "# - ", calls[[i]]$tool, "Calls = ", 
    calls[[i]]$calls)
}

report_conn <- file(outReport)
writeLines(comment_line, report_conn)
close(report_conn)
write_tsv(report, outReport, append = TRUE, col_names = TRUE)

breakpoints_conn <- file(outBreakpoints)
writeLines(comment_line, breakpoints_conn)
close(breakpoints_conn)
write_delim(all_breakpoints, outBreakpoints, append = TRUE, col_names = TRUE)
