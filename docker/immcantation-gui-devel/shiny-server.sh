#!/bin/sh

# Create log dir
mkdir -p /var/log/shiny-server
chown shiny.shiny /var/log/shiny-server

mkdir -p /var/lib/shiny-server/bookmarks
chown shiny.shiny /var/lib/shiny-server/bookmarks

exec shiny-server >> /var/log/shiny-server.log 2>&1 
