#!/usr/bin/env Rscript
# Super script to run TIgGER polymorphism detection and genotyping
#
# Author:  Jason Anthony Vander Heiden
# Date:    2018.10.05
#
# Arguments:
#   -d  Tabulated data, in Change-O (TAB) or AIRR (TSV) format.
#   -r  FASTA file containing IMGT-gapped V segment reference germlines.
#       Defaults to /usr/local/share/germline   s/imgt/human/vdj/imgt_human_IGHV.fasta.
#   -v  Name of the output field containing genotyped V assignments.
#       Defaults to V_CALL_GENOTYPED.
#   -x  Minimum number of sequences in the mutation/coordinate range. Samples with 
#       insufficient sequences will be excluded. Defaults to 50.
#   -y  Minimum number of sequences required to analyze a germline allele. Defaults to 200.
#   -n  Sample name or run identifier which will be used as the output file prefix.
#       Defaults to a truncated version of the input filename.
#   -u  Whether to use '-r' to find which samples are unmutated. Default is True.
#   -o  Output directory. Will be created if it does not exist.
#       Defaults to the current working directory.
#   -f  File format. One of 'airr' (default) or 'changeo'.
#   -p  Number of subprocesses for multiprocessing tools.
#       Defaults to the available processing units.
#   -h  Display help.

# Imports
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("alakazam"))
suppressPackageStartupMessages(library("shazam"))
suppressPackageStartupMessages(library("tigger"))

# Set defaults
NPROC <- parallel::detectCores()
FORMAT <- "airr"
MIN_SEQS <- 50
GERMLINE_MIN <- 200
FIND_UNMUTATED <- TRUE

# Define commmandline arguments
opt_list <- list(make_option(c("-d", "--db"), dest="DB",
                             help="Change-O formatted TSV (TAB) file."),
                 make_option(c("-r", "--ref"), dest="REF",
                             default="/usr/local/share/germlines/imgt/human/vdj/imgt_human_IGHV.fasta",
                             help=paste("FASTA file containing IMGT-gapped V segment reference germlines.",
                                        "\n\t\tDefaults to /usr/local/share/germlines/imgt/human/vdj/imgt_human_IGHV.fasta.")),
                 make_option(c("-v", "--vfield"), dest="VFIELD",
                             default="V_CALL_GENOTYPED",
                             help=paste("Name of the output field containing genotyped V assignments.",
                                        "\n\t\tDefaults to V_CALL_GENOTYPED.")),
                 make_option(c("-x", "--minseq"), dest="MIN_SEQS",
                             default=MIN_SEQS,
                             help=paste("Minimum number of sequences in the mutation/coordinate range.",
                                        "\n\t\tSamples with insufficient sequences will be excluded.",
                                        "\n\t\tDefaults to 50.")),
                 make_option(c("-y", "--mingerm"), dest="GERMLINE_MIN",
                             default=GERMLINE_MIN,
                             help=paste("Minimum number of sequences required to analyze a germline allele.",
                                        "\n\t\tDefaults to 200.")),
                 make_option(c("-n", "--name"), dest="NAME",
                             help=paste("Sample name or run identifier which will be used as the output file prefix.",
                                        "\n\t\tDefaults to a truncated version of the input filename.")),
                 make_option(c("-u", "--find-unmutated"), dest="FIND_UNMUTATED",
                             default=FIND_UNMUTATED,
                             help=paste("Whether to use '-r' to find which samples are unmutated.",
                                        "\n\t\tDefaults to TRUE.")),                 
                 make_option(c("-o", "--outdir"), dest="OUTDIR", default=".",
                             help=paste("Output directory. Will be created if it does not exist.",
                                        "\n\t\tDefaults to the current working directory.")),
                 make_option(c("-f", "--format"), dest="FORMAT", default=FORMAT,
                             help=paste("File format. One of 'airr' (default) or 'changeo'.")),                 
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

# Check unmutated
opt$FIND_UNMUTATED <- tolower(as.character(opt$FIND_UNMUTATED))
if (opt$FIND_UNMUTATED %in% c("true", "t", "1")) {
    opt$FIND_UNMUTATED <- TRUE
} else if  (opt$FIND_UNMUTATED %in% c("false", "f", "0"))  {
    opt$FIND_UNMUTATED <- FALSE
} else {
    stop("Unexpected find unmutated.")
}
 
# Create output directory
if (!(dir.exists(opt$OUTDIR))) {
    dir.create(opt$OUTDIR)
}

# Check write access
if (!(file.access(opt$OUTDIR, mode=2) == 0)) {
    stop("Output directory '", opt$OUTDIR, "' cannot be written to.")
}

# Check alakazam version, to determine column names
if (numeric_version(packageVersion("alakazam")) > numeric_version('0.3.0')) {
    # lower case
    v_call_genotyped <- "v_call_genotyped"
} else {
    # upper case
    v_call_genotyped <- "V_CALL_GENOTYPED"
}

# Load data
if (opt$FORMAT == "changeo") {
    db <- as.data.frame(alakazam::readChangeoDb(opt$DB))
    v_call <- "V_CALL"
    j_call <- "J_CALL"
    junction <- "JUNCTION"
    junction_length <- "JUNCTION_LENGTH"
    sequence_alignment <- "SEQUENCE_IMGT"
    ext <- "tab"
} else if (opt$FORMAT == "airr") {
    db <- airr::read_rearrangement(opt$DB)
    v_call <- "v_call"
    j_call <- "j_call"
    junction <- "junction"
    junction_length <- "junction_length"
    sequence_alignment <- "sequence_alignment"
    ext <- "tsv"
}

igv <- readIgFasta(opt$REF)

# Identify polymorphisms and genotype
nv <- findNovelAlleles(db, germline_db=igv, v_call=v_call, j_call=j_call,
                       seq=sequence_alignment, junction=junction,
                       junction_length=junction_length,
                       min_seqs=opt$MIN_SEQS, germline_min=opt$GERMLINE_MIN,
                       nproc=opt$NPROC)
gt <- inferGenotype(db, germline_db=igv, novel=nv,
                    v_call=v_call, seq=sequence_alignment, 
                    find_unmutated = opt$FIND_UNMUTATED)

# Write genotype FASTA file
gt_seq <- genotypeFasta(gt, germline_db=igv, novel=nv)
writeFasta(gt_seq, file.path(opt$OUTDIR, paste0(opt$NAME, "_genotype.fasta")))

# Modify allele calls
db <- reassignAlleles(db, gt_seq, v_call=v_call, seq=sequence_alignment)

# Rename genotyped V call column if necessary
if (opt$VFIELD != v_call_genotyped) {
    db[[opt$VFIELD]] <- db[[v_call_genotyped]]
    db <- dplyr::select(db, -v_call_genotyped)
}

# Write genotyped data
writeChangeoDb(db, file.path(opt$OUTDIR, paste0(opt$NAME, "_genotyped.",ext)))

# Plot genotype
plot_file <- file.path(opt$OUTDIR, paste0(opt$NAME, "_genotype.pdf"))
pdf(plot_file, width=7, height=10, useDingbats=FALSE)
plotGenotype(gt, silent=FALSE)
dev.off()
