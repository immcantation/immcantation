#!/usr/bin/env Rscript
# Build an R package
#
# Author:  Jason Anthony Vander Heiden
# Date:    2018.12.08
#
# Arguments:
#   -p    Package source directory.
#   -h    Display help.

# Imports
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("devtools"))
suppressPackageStartupMessages(library("pkgbuild"))

# Set defaults
PKG_DIR <- "."
UPGRADE <- "default"

# Define commmandline arguments
opt_list <- list(make_option(c("-p", "--package"), dest="PKG_DIR", default=PKG_DIR,
                             help="Package source directory. Defaults to current directory."),
                make_option(c("-u", "--upgrade"), dest="UPGRADE", default=UPGRADE,
                             help="Whether package dependencies should be upgraded.")
                             )

# Parse arguments
opt <- parse_args(OptionParser(option_list=opt_list))

# Build
setwd(opt$PKG_DIR)
install_deps(dependencies=TRUE, upgrade=opt$UPGRADE)
compile_dll()
document()
install(build_vignettes=TRUE, upgrade=opt$UPGRADE)
