#!/usr/bin/env Rscript
# Build an R package
#
# Author:  Jason Anthony Vander Heiden
# Date:    2018.12.08
#
# Arguments:
#   -p    Package source directory.
#   -u    Upgrade argument to pass to devtools install
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
                             help="Whether package dependencies should be upgraded."))

# Parse arguments
opt <- parse_args(OptionParser(option_list=opt_list))

# Build
setwd(opt$PKG_DIR)
install_deps(dependencies=TRUE, upgrade=opt$UPGRADE)
compile_dll()
document()
# Since devtools 2.5.0, install() (unlike install_deps()) requires upgrade to
# be a single TRUE/FALSE/NA, not a "never"/"always" string.
upgrade_flag <- switch(opt$UPGRADE,
                       "always"=TRUE,
                       "never"=FALSE,
                       "TRUE"=TRUE,
                       "FALSE"=FALSE,
                       FALSE)
install(build_vignettes=TRUE, upgrade=upgrade_flag)
