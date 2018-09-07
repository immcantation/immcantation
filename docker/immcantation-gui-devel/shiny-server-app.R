library(shiny)
dir <- system.file('shiny-app', package = 'immcantation')
setwd(dir)
shiny::shinyAppDir('.')
