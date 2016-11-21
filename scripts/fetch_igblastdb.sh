#!/usr/bin/env bash
# Download IgBLAST database files
#
# Author:  Mohamed Uduman
# Date:    2016.11.21
#
# Arguments:
#   $1 = folder to download database to

# Parameters
DB_DIR=$1

# Fetch internal_data
wget -q --show-progress -nd -r ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/internal_data \
    -P ${DB_DIR}/internal_data

# Fetch database
wget -q --show-progress -nd -r ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/database \
    -P ${DB_DIR}/database

# Fetch optional_file
wget -q --show-progress -nd -r ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/optional_file  \
    -P ${DB_DIR}/optional_file