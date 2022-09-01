#!/usr/bin/env Rscript
# Super script to run SHazaM distance to nearest tuning
#
# Author:  Jason Anthony Vander Heiden, Ruoyi Jiang
# Date:    2019.03.19
#
# Arguments:
#   -d           Change-O formatted TSV (TAB) file.
#   -m           Method to use for determining the optimal threshold.
#                Defaults to density.
#   -n           Sample name or run identifier which will be used as the output file prefix.
#                Defaults to a truncated version of the input filename.
#   -o           Output directory. Will be created if it does not exist.
#                Defaults to a directory matching the sample identifier in the current working directory.
#   -f           File format. One of 'airr' (default) or 'changeo'.
#   -p           Number of subprocesses for multiprocessing tools.
#                Defaults to the available processing units.
#   -h           Display help.

# Imports
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("alakazam"))
suppressPackageStartupMessages(library("airr"))
suppressPackageStartupMessages(library("scoper"))

# Set defaults
METHOD <- "nt"
OUTDIR <- "."
FORMAT <- "airr"
MODEL <- "gamma-gamma"
CUTOFF <- "optimal"
SPC <- 0.995
NPROC <- parallel::detectCores()

# Define commmandline arguments
opt_list <- list(make_option(c("-d", "--db"), dest="DB",
                             help="Tabulated data file, in Change-O (TAB) or AIRR format (TSV)."),
                 make_option(c("-t", "--threshold"), dest="THRESHOLD",
                             help=paste("distance threshold for clonal grouping.",
                                        "\n\t\t.One of 'nt' (nucleotide based clustering) or 'aa' (amino acid).")),                 
                 make_option(c("-m", "--method"), dest="METHOD", default=METHOD,
                             help=paste("Distance method for clonal assignment.",
                                        "\n\t\t.One of 'nt' (nucleotide based clustering) or 'aa' (amino acid).",
                                        "\n\t\tDefaults to 'nt'.")),                 
                 make_option(c("-n", "--name"), dest="NAME",
                             help=paste("Sample name or run identifier which will be used as the output file prefix.",
                                        "\n\t\tDefaults to a truncated version of the input filename.")),
                 make_option(c("-o", "--outdir"), dest="OUTDIR", default=OUTDIR,
                             help=paste("Output directory. Will be created if it does not exist.",
                                        "\n\t\tDefaults to the current working directory.")),
                 make_option(c("-f", "--format"), dest="FORMAT", default=FORMAT,
                             help=paste("File format. One of 'airr' (default) or 'changeo'.")),
                 make_option(c("-p", "--nproc"), dest="NPROC", default=NPROC,
                             help=paste("Number of subprocesses for multiprocessing tools.",
                                        "\n\t\tDefaults to the available processing units."))
)
                
# Parse arguments
opt <- parse_args(OptionParser(option_list=opt_list))

# Check input file
if (!("DB" %in% names(opt))) {
    stop("You must provide a database file with the -d option.")
}

# Check threshold
if (!("THRESHOLD" %in% names(opt))) {
   stop("You must provide a threshold value with the -t option.")
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
METHOD <- opt$METHOD
THRESHOLD <- opt$THRESHOLD
OUTDIR <- opt$OUTDIR
FORMAT <- opt$FORMAT
NPROC <- opt$NPROC
NAME <- opt$NAME

db_files <- strsplit(DB,",")[[1]]

# Load data
if (FORMAT == "changeo") {
   db <- bind_rows(
      lapply(db_files, alakazam::readChangeoDb),
      .id="input_id"
   )
   v_call <- "V_CALL"
   v_call <- "V_CALL"
   j_call <- "J_CALL"
   junction <- "JUNCTION"
   clone <- "CLONE"
   ext <- "tab"
} else if (FORMAT == "airr") {
   db <- bind_rows(
      lapply(db_files, airr::read_rearrangement),
      .id="input_id"
   )   
   
   v_call <- "v_call"
   j_call <- "j_call"
   junction <- "junction"
   clone <- "clone_id"
   ext <- "tsv"
}

input_size <- nrow(db)

if ("cell_id" %in% colnames(db)) {
   cell_id <- "cell_id"
} else {
   cell_id <- NULL
}

db <- suppressWarnings(suppressMessages(hierarchicalClones(db, 
                         threshold = THRESHOLD, 
                         method = METHOD, 
                         cell_id = cell_id, 
                         clone = clone,
                         only_heavy = FALSE, 
                         split_light = TRUE, 
                         summarize_clones = FALSE,
                         nproc=NPROC)))

# calculate and plot the rank-abundance curve
heavy_chains <- c("IGH", "TRB", "TRD")
db_h <- db %>%
   filter(getLocus(!!rlang::sym(v_call)) %in% heavy_chains)
ab <- estimateAbundance(db_h)
ab_plot <- plotAbundanceCurve(ab, colors = "steelblue", silent = T) + 
   labs(title="Clonal abundance distribution.")

## Save files and print log
for (i in 1:length(db_files)) {
   cat("INPUT",i,"> ", basename(db_files[i]), "\n", sep="")
   
   # Check and fill sample name
   if (!("NAME" %in% names(opt))) {
      n <- basename(db_files[i])
      this_name <- tools::file_path_sans_ext(n)
   } else {
      this_name <- paste(NAME,i,sep="-")
   }
   out_file <- file.path(OUTDIR, paste0(this_name, "_clone-pass.",ext))
   cat("OUTPUT",i,"> ", out_file, "\n", sep="")
   
   # Save db
   writeChangeoDb(db %>% 
                     filter(input_id == i) %>% 
                     select(-input_id), 
                  out_file)
   if (i == 1) {
      #Save abundance
      save(ab, file=file.path(OUTDIR, paste0(this_name, "_clone-pass_abundance.RData")))
      ggsave(plot=ab_plot,
             device="pdf",
             width=6,
             height=3,
             filename=file.path(OUTDIR, paste0(this_name, "_clone-pass_abundance.pdf")))
   }
}
cat("INPUT_SIZE> ", input_size, "\n", sep="")
cat("OUTPUT_SIZE> ", nrow(db), "\n", sep="")
cat("NUM_CLONES> ", length(unique(db[[clone]])), "\n", sep="")