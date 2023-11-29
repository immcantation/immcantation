#!/usr/bin/env bash
set -o nounset
# Use user provided passord
#echo "rstudio:${PASSWORD}" | chpasswd
echo "${RS_USER}:${PASSWORD}" | chpasswd
echo "Starting server"
exec /usr/bin/rserver --server-working-dir=/data --server-daemonize=false --www-port=8787 >> /var/log/rstudio-server/rstudio-server.log 2>&1  
