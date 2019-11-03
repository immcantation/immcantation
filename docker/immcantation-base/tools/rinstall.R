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

# Define commmandline arguments
opt_list <- list(make_option(c("-p", "--package"), dest="PKG_DIR", default=PKG_DIR,
                             help="Package source directory. Defaults to current directory."))

# Parse arguments
opt <- parse_args(OptionParser(option_list=opt_list))

# Build
setwd(opt$PKG_DIR)
# Added devtools.ellipsis_action: https://github.com/r-lib/devtools/issues/2109
options(devtools.ellipsis_action = rlang::signal)
install_deps(dependencies=TRUE, upgrade=TRUE, clean=TRUE)
compile_dll()
document()
install(build_vignettes=TRUE, clean=TRUE)
