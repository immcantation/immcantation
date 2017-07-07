#!/usr/bin/env Rscript
# Super script to run TIgGER 0.2.10 genotyping
#
# Author:  Jason Anthony Vander Heiden
# Date:    2017.07.07
#
# Arguments:
#   -d  Change-O formatted TSV (TAB) file.
#   -r  FASTA file containing IMGT-gapped V segment reference germlines.
#       Defaults to /usr/local/share/germlines/imgt/human/vdj/imgt_human_IGHV.fasta.
#   -n  Sample name or run identifier which will be used as the output file prefix.
#       Defaults to a truncated version of the input filename.
#   -o  Output directory.
#       Defaults to current directory.
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
NPROC <- shazam::getnproc()

# Define commmandline arguments
opt_list <- list(make_option(c("-d", "--db"), dest="DB",
                             help="Change-O formatted TSV (TAB) file."),
                 make_option(c("-r", "--ref"), dest="REF",
                             default="/usr/local/share/germlines/imgt/human/vdj/imgt_human_IGHV.fasta",
                             help=paste("FASTA file containing IMGT-gapped V segment reference germlines.",
                                        "\n\t\tDefaults to /usr/local/share/germlines/imgt/human/vdj/imgt_human_IGHV.fasta.")),
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
igv <- readIgFasta(opt$REF)

# Identify polymorphisms and genotype
nv <- findNovelAlleles(db, germline_db=igv, nproc=opt$NPROC)
gt <- inferGenotype(db, germline_db=igv, novel_df=nv)

# Write genotype FASTA file
gt_seq <- genotypeFasta(gt, germline_db=igv, novel_df=nv)
writeFasta(gt_seq, file.path(opt$OUTDIR, paste0(opt$NAME, "_genotype.fasta")))

# Modify allele calls and write db
db <- cbind(db, reassignAlleles(db, gt_seq))
writeChangeoDb(db, file.path(opt$OUTDIR, paste0(opt$NAME, "_genotyped.tab")))

# Plot genotype
ggsave(file.path(opt$OUTDIR, paste0(opt$NAME, "_genotype.pdf")),
       plotGenotype(gt, silent=TRUE))

