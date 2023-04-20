#!/usr/bin/env Rscript
# Super script to remove single cell doublets and light chain only cells
#
# Author:  Edel Aron, Susanna Marquez, Hailong Meng
# Date:    2022.08.31
#
# Arguments:
#   -d           AIRR or Change-O formatted TSV (TAB) file(s). Comma separated.
#   -n           Sample name or run identifier which will be used as the output file prefix.
#                Defaults to a truncated version of the input filename.
#   -o           Output directory. Will be created if it does not exist.
#                Defaults to a directory matching the sample identifier in the current working directory.
#   -f           File format. One of 'airr' (default) or 'changeo'.
#   -h           Display help.

# Imports
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("alakazam"))
suppressPackageStartupMessages(library("airr"))
suppressPackageStartupMessages(library("dplyr"))

# Set defaults
OUTDIR <- "."
FORMAT <- "airr"

heavy_chains <- c("IGH", "TRB", "TRD")

# Define commmandline arguments
opt_list <- list(make_option(c("-d", "--db"), dest="DB",
                             help="Tabulated data files, in Change-O (TAB) or AIRR format (TSV)."),
                 make_option(c("-n", "--name"), dest="NAME",
                             help=paste("Sample name or run identifier which will be used as the output file prefix.",
                                        "\n\t\tDefaults to a truncated version of the first input filename.")),
                 make_option(c("-o", "--outdir"), dest="OUTDIR", default=OUTDIR,
                             help=paste("Output directory. Will be created if it does not exist.",
                                        "\n\t\tDefaults to the current working directory.")),
                 make_option(c("-f", "--format"), dest="FORMAT", default=FORMAT,
                             help=paste("File format. One of 'airr' (default) or 'changeo'."))
                 )

# Parse arguments
opt <- parse_args(OptionParser(option_list=opt_list))

cat("\n\n            START> singlecell-filter\n")

# Check input file
if (!("DB" %in% names(opt))) {
    stop("You must provide a database file with the -d option.")
}

# Create output directory
if (!(dir.exists(opt$OUTDIR))) {
    dir.create(opt$OUTDIR)
}

# Check write access
if (!(file.access(opt$OUTDIR, mode=2) == 0)) {
    stop("Output directory '", opt$OUTDIR, "' cannot be written to.")
}

# Reset parameters from opt (better for debugging)
DB <- opt$DB
OUTDIR <- opt$OUTDIR
FORMAT <- opt$FORMAT

db_files <- strsplit(DB,",")[[1]]

# Load data
if (FORMAT == "changeo") {
   db <- bind_rows(
      lapply(db_files, alakazam::readChangeoDb),
      .id="input_id"
   )
   v_call <- "V_CALL"
   ext <- "tab"
} else if (FORMAT == "airr") {
   db <- bind_rows(
      lapply(db_files, airr::read_rearrangement),
      .id="input_id"
   )
   v_call <- "v_call"
   ext <- "tsv"
}

input_size <- nrow(db)

if ("locus" %in% colnames(db) == F) {
   db[["locus"]] <- getLocus(db[[v_call]])
}

input_heavy <- db %>% filter(locus %in% heavy_chains) %>% nrow()
input_light <- db %>% filter(locus %in% heavy_chains == F) %>% nrow()

# Remove cells with multiple heavy chain
multi_heavy <- table(filter(db, locus %in% heavy_chains)[['cell_id']])
multi_heavy_cells <- names(multi_heavy)[multi_heavy > 1]

db <- filter(db, !cell_id %in% multi_heavy_cells)
#message("There are ", nrow(db), " rows in the data after filtering out cells with multiple heavy chains.")

# Split cells by heavy and light chains
heavy_cells <- filter(db, locus %in% heavy_chains)[['cell_id']]
light_cells <- filter(db, locus %in% heavy_chains == F )[['cell_id']]

# Identify and remove cells without heavy chain
no_heavy_cells <- light_cells[which(!light_cells %in% heavy_cells)]

db <- filter(db, !cell_id %in% no_heavy_cells)

output_size <- nrow(db)
output_heavy <- db %>% filter(locus %in% heavy_chains) %>% nrow()
output_light <- db %>% filter(locus %in% heavy_chains == F) %>% nrow()

#message("There are ", output_size, " rows in the data after filtering out cells without heavy chains.")


## Save files and print log
for (i in 1:length(db_files)) {
   cat("           INPUT",i,"> ", basename(db_files[i]), "\n", sep="")

   # Check and fill sample name
   if (!("NAME" %in% names(opt))) {
      n <- basename(db_files[i])
      this_name <- tools::file_path_sans_ext(n)
   } else {
      this_name <- paste(opt$NAME,i,sep="-")
   }
   out_file <- file.path(OUTDIR, paste0(this_name, "_sc-pass.",ext))
   cat("          OUTPUT",i,"> ", out_file, "\n", sep="")

   writeChangeoDb(db %>%
                     filter(input_id == i) %>%
                     select(-input_id),
                  out_file)
}


cat("          RECORDS> ", input_size, "\n", sep="")
cat("    RECORDS_HEAVY> ", input_heavy, "\n", sep="")
cat("    RECORDS_LIGHT> ", input_light, "\n", sep="")
cat("MULTI_HEAVY_CELLS> ", length(multi_heavy_cells), "\n", sep="")
cat("   NO_HEAVY_CELLS> ", length(no_heavy_cells), "\n", sep="")
cat("             PASS> ", output_size, "\n", sep="")
cat("       PASS_HEAVY> ", output_heavy, "\n", sep="")
cat("       PASS_LIGHT> ", output_light, "\n", sep="")
cat("              END> singlecell-filter\n\n")