#!/usr/bin/env Rscript
# Super script to run SHazaM distance to nearest tuning
#
# Author:  Jason Anthony Vander Heiden, Ruoyi Jiang
# Date:    2019.03.19
#
# Arguments:
#   -d           Change-O formatted TSV (TAB) file.
#   -m           Method.
#                Defaults to density.
#   -n           Sample name or run identifier which will be used as the output file prefix.
#                Defaults to a truncated version of the input filename.
#   -o           Output directory. Will be created if it does not exist.
#                Defaults to a directory matching the sample identifier in the current working directory.
#   -f           File format. One of 'changeo' (default) or 'airr'.
#   -p           Number of subprocesses for multiprocessing tools.
#                Defaults to the available processing units.
#   --model      Model when "-m gmm" is specified.
#                Defaults to "gamma-gamma".
#   --subsample  Number of distances to downsample to before threshold calculation.
#                By default, subsampling is not performed.
#   --repeats    Number of times to repeat the threshold calculation (with plotting).
#   -h           Display help.

# Imports
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("methods"))
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
FORMAT <- "changeo"
MODEL <- "gamma-gamma"
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
                             help=paste("File format. One of 'changeo' (default) or 'airr'.")),
                 make_option(c("-p", "--nproc"), dest="NPROC", default=NPROC,
                             help=paste("Number of subprocesses for multiprocessing tools.",
                                        "\n\t\tDefaults to the available processing units.")),
                 make_option(c("--model"), dest="MODEL", default=MODEL,
                             help=paste("Model to use for the gmm model.",
                                        "\n\t\tOne of gamma-gamma, gamma-norm, norm-norm or norm-gamma.",
                                        "\n\t\tDefaults to gamma-gamma.")),
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

# Calculate distance-to-nearest
db <- distToNearest(db, sequenceColumn=junction, vCallColumn=v_call, jCallColumn=j_call,
                    model="ham", first=FALSE, normalize="len", nproc=NPROC)

# Simply plot and exit for method="none"
if (METHOD == "none") {
    # Plot distToNearest distribution
    p1 <- ggplot(db, aes_string(x="DIST_NEAREST")) +
        baseTheme() +
        xlab("Distance") +
        ylab("Density") +
        geom_histogram(aes_string(y="..density.."), fill="gray40", color="white")
    ggsave(file.path(OUTDIR, paste0(NAME, "_threshold-plot.pdf")), plot=p1, width=6, height=4)

    quit()
}

# Open plot device
plot_file <- file.path(OUTDIR, paste0(NAME, "_threshold-plot.pdf"))
pdf(plot_file, width=6, height=4, useDingbats=FALSE)

# Repeat threshold calculations and plot
threshold_list <- list()
for(i in 1:REPEATS){
    # Subsample distances
    if (is.null(SUBSAMPLE) || length(db$DIST_NEAREST) < SUBSAMPLE){
        sampling <- db$DIST_NEAREST
    } else {
        sampling <- sample(db$DIST_NEAREST, SUBSAMPLE, replace=FALSE)
    }

    # Calculate threshold
    threshold <- findThreshold(sampling, method=METHOD, model=MODEL)

    # Build results data.frame
    slots <- slotNames(threshold)
    slots <- slots[!(slots %in% c("x", "xdens", "ydens"))]
    .extract <- function(x) {
        return(data_frame(PARAMETER=x, VALUE=as.character(slot(threshold, x))))
    }
    threshold_list[[as.character(i)]] <- bind_rows(lapply(slots, .extract))
    # Plot histogram
    plot(threshold, binwidth=0.02, silent=FALSE)
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