#!/usr/bin/env Rscript
# Super script to run SHazaM 0.1.7 distance to nearest tuning
#
# Author:  Jason Anthony Vander Heiden
# Date:    2017.05.26
#
# Arguments:
#   -d  Change-O formatted TSV (TAB) file.
#   -m  Method.
#       Defaults to gmm.
#   -n  Sample name or run identifier which will be used as the output file prefix.
#       Defaults to a truncated version of the input filename.
#   -o  Output directory.
#       Defaults to current directory.
#   -p  Number of subprocesses for multiprocessing tools.
#       Defaults to the available processing units.
#   -h  Display help.

# Imports
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("methods"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("alakazam"))
suppressPackageStartupMessages(library("shazam"))

# Set defaults
NPROC <- shazam::getnproc()

# Define commmandline arguments
opt_list <- list(make_option(c("-d", "--db"), dest="DB",
                             help="Change-O formatted TSV (TAB) file."),
                 make_option(c("-m", "--method"), dest="METHOD", default="gmm",
                             help=paste("Threshold inferrence to use. One of gmm or dens.", 
                                        "\n\t\tDefaults to gmm.")),
                 make_option(c("-n", "--name"), dest="NAME",
                             help=paste("Sample name or run identifier which will be used as the output file prefix.",
                                        "\n\t\tDefaults to a truncated version of the input filename.")),
                 make_option(c("-o", "--outdir"), dest="OUTDIR", default=".",
                             help=paste("Output directory.", "Defaults to the sample name.")),
                 make_option(c("-p", "--nproc"), dest="NPROC", default=NPROC,
                             help=paste("Number of subprocesses for multiprocessing tools.",
                                        "\n\t\tDefaults to the available processing units.")))
# Parse arguments
opt <- parse_args(OptionParser(option_list=opt_list))

# Check input file
if (!("DB" %in% names(opt))) {
    stop("You must provide a Change-O database file with the -d option.")
}

# Check and fill sample name
if (!("NAME" %in% names(opt))) {
    n <- basename(opt$DB)
    opt$NAME <- tools::file_path_sans_ext(basename(opt$DB))
}

# Create output directory
if (!(dir.exists(opt$OUTDIR))) {
    dir.create(opt$OUTDIR)
}

# Load data
db <- as.data.frame(readChangeoDb(opt$DB))

# Calculate distance to nearest and threshold
db <- distToNearest(db, model="ham", first=FALSE, normalize="len", nproc=opt$NPROC)
threshold <- findThreshold(db$DIST_NEAREST, method=opt$METHOD)

# Extract relevant slots into data_frame
slots <- slotNames(threshold)
slots <- slots[!(slots %in% c("x", "xdens", "ydens"))]
.extract <- function(x) {
    data_frame(PARAMETER=x, VALUE=as.character(slot(threshold, x)))
}
thresh_df <- bind_rows(lapply(slots, .extract))

# Print and save threshold table
cat("THRESHOLD> ", threshold@threshold, "\n", sep="")
write_tsv(thresh_df, file.path(opt$OUTDIR, paste0(opt$NAME, "_threshold-values.tab")))

# Plot
p1 <- plot(threshold, binwidth=0.02, silent=TRUE)
ggsave(file.path(opt$OUTDIR, paste0(opt$NAME, "_threshold-plot.pdf")), plot=p1, 
       device="pdf", width=6, height=4, useDingbats=FALSE)
