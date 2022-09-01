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
#   --model      Model when "-m gmm" is specified.
#                Defaults to "gamma-gamma".
#   --cutoff     Method to use for threshold selection.
#                Defaults to "optimal".
#   --spc        Specificity required for threshold selection. Applies only when
#                method='gmm' and cutoff='user'.
#                Defaults to 0.995.
#   --subsample  Number of distances to downsample to before threshold calculation.
#                By default, subsampling is not performed.
#   --repeats    Number of times to repeat the threshold calculation (with plotting).
#   -h           Display help.

# Imports
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("methods"))
suppressPackageStartupMessages(library("tibble"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("tidyr"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("alakazam"))
suppressPackageStartupMessages(library("shazam"))
suppressPackageStartupMessages(library("airr"))

# Set defaults
METHOD <- "density"
OUTDIR <- "."
FORMAT <- "airr"
MODEL <- "gamma-gamma"
CUTOFF <- "optimal"
SPC <- 0.995
NPROC <- parallel::detectCores()
SUBSAMPLE <- NULL
REPEATS <- 1

# Define commmandline arguments
opt_list <- list(make_option(c("-d", "--db"), dest="DB",
                             help="Tabulated data file, in Change-O (TAB) or AIRR format (TSV)."),
                 make_option(c("-m", "--method"), dest="METHOD", default=METHOD,
                             help=paste("Threshold inferrence to use. One of gmm, density, or none.",
                                        "\n\t\tIf none, the distance-to-nearest distribution is plotted without threshold detection.",
                                        "\n\t\tDefaults to density.")),
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
                                        "\n\t\tDefaults to the available processing units.")),
                 make_option(c("--model"), dest="MODEL", default=MODEL,
                             help=paste("Model to use for the gmm model.",
                                        "\n\t\tOne of gamma-gamma, gamma-norm, norm-norm or norm-gamma.",
                                        "\n\t\tDefaults to gamma-gamma.")),
                 make_option(c("--cutoff"), dest="CUTOFF", default=CUTOFF,
                             help=paste("Method to use for threshold selection.",
                                        "\n\t\tOne of optimal, intersect or user.",
                                        "\n\t\tDefaults to optimal.")),
                 make_option(c("--spc"), dest="SPC", default=SPC,
                             help=paste("Specificity required for threshold selection.",
                                        "\n\t\tApplies only when method='gmm' and cutoff='user'.",
                                        "\n\t\tDefaults to 0.995.")),
                 make_option(c("--subsample"), dest="SUBSAMPLE", default=SUBSAMPLE,
                             help=paste("Number of distances to downsample the data to before threshold calculation.",
                                        "\n\t\tBy default, subsampling is not performed.")),
                 make_option(c("--repeats"), dest="REPEATS", default=REPEATS,
                             help=paste("Number of times to recalculate.",
                                        "\n\t\tDefaults to 1.")))

# Parse arguments
opt <- parse_args(OptionParser(option_list=opt_list))

# Check input file
if (!("DB" %in% names(opt))) {
    stop("You must provide a database file with the -d option.")
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

# Check write access
if (!(file.access(opt$OUTDIR, mode=2) == 0)) {
    stop("Output directory '", opt$OUTDIR, "' cannot be written to.")
}

# Reset parameters from opt (better for debugging)
DB <- opt$DB
METHOD <- opt$METHOD
OUTDIR <- opt$OUTDIR
FORMAT <- opt$FORMAT
MODEL <- opt$MODEL
CUTOFF <- opt$CUTOFF
SPC <- opt$SPC
NPROC <- opt$NPROC
NAME <- opt$NAME
SUBSAMPLE <- opt$SUBSAMPLE
REPEATS <- opt$REPEATS

# Load data
if (FORMAT == "changeo") {
    db <- as.data.frame(alakazam::readChangeoDb(DB))
    v_call <- "V_CALL"
    j_call <- "J_CALL"
    junction <- "JUNCTION"
} else if (FORMAT == "airr") {
    db <- airr::read_rearrangement(DB)
    v_call <- "v_call"
    j_call <- "j_call"
    junction <- "junction"
}

if ("cell_id" %in% colnames(db)) {
   cell_id <- "cell_id"
} else {
   cell_id <- NULL
}

# Check alakazam version, to determine column names
if (numeric_version(packageVersion("alakazam")) > numeric_version('0.3.0')) {
    # lower case
    dist_nearest <- "dist_nearest"
} else {
    # upper case
    dist_nearest <- "DIST_NEAREST"
}

# Calculate distance-to-nearest
db <- suppressWarnings(distToNearest(db,
                                    sequenceColumn=junction,
                                    vCallColumn=v_call,
                                    jCallColumn=j_call,
                                    cellIdColumn=cell_id,
                                    model="ham", first=FALSE,
                                    normalize="len", nproc=NPROC))

# Simply plot and exit for method="none"
if (METHOD == "none") {
    # Plot distToNearest distribution
    p1 <- ggplot(filter(db, !is.na(dist_nearest), dist_nearest > 0), aes(x=dist_nearest)) +
        baseTheme() +
        xlab("Distance") +
        ylab("Density") +
        geom_histogram(aes(y=..density..), binwidth=0.02, fill="gray40", color="white")
    f <- file.path(OUTDIR, paste0(NAME, "_threshold-plot.pdf"))
    suppressWarnings(ggsave(f, plot=p1, width=6, height=4))

    quit()
}

# Open plot device
f <- file.path(OUTDIR, paste0(NAME, "_threshold-plot.pdf"))
pdf(f, width=6, height=4, useDingbats=FALSE)

# Repeat threshold calculations and plot
threshold_list <- list()
for(i in 1:REPEATS) {
    # Subsample distances
    if (is.null(SUBSAMPLE) || length(db[[dist_nearest]]) < SUBSAMPLE){
        sampling <- db[[dist_nearest]]
    } else {
        sampling <- sample(db[[dist_nearest]], SUBSAMPLE, replace=FALSE)
    }

    # Calculate threshold
    threshold <- findThreshold(sampling,
                                method=METHOD,
                                model=MODEL,
                                cutoff=CUTOFF,
                                spc=SPC)

    # Build results data.frame
    slots <- slotNames(threshold)
    slots <- slots[!(slots %in% c("x", "xdens", "ydens"))]
    .extract <- function(x) {
        return(tibble(PARAMETER=x, VALUE=as.character(slot(threshold, x))))
    }
    threshold_list[[as.character(i)]] <- bind_rows(lapply(slots, .extract))
    # Plot histogram
    suppressWarnings(plot(threshold, binwidth=0.02, silent=FALSE))
}
# Close plot
dev.off()

# Build data.frame of replicates
thresh_df <- bind_rows(threshold_list, .id = "REPEAT") %>%
    spread(PARAMETER, VALUE) %>%
    select(-REPEAT)

# Print and save threshold table
cat("THRESHOLD_AVG> ", mean(as.numeric(thresh_df$threshold), na.rm = TRUE), "\n", sep="")
write_tsv(thresh_df, file.path(OUTDIR, paste0(NAME, "_threshold-values.tab")))