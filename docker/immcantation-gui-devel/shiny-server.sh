#!/usr/bin/env bash

exec shiny-server >> /var/log/shiny-server/shiny-server.log 2>&1 
