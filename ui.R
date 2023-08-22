library(shiny)
library(shinyWidgets)
library(survival)
library(tidyverse)
library(ggsurvfit)
library(plotly)
library(shinyjs)

shinyUI(fluidPage(
    withMathJax(),
    useShinyjs(),
    # Application title
    titlePanel("Survival models"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(verticalLayout(
          radioButtons("radioBtnMode", label = "Choose a mode:", choices = c("Descriptive analysis", "Model analysis")),
          pickerInput(inputId =  "dataset", label = "Data sets", choices = c("All data", "Progression-free survival", "Overall survival")),
          hidden(pickerInput(inputId = "dist", label = "Models",
                      choices = c("All distributions", "Weibull", "Exponential", "Log-logistic", "Rayleigh"),
          ))
        )),

        # Main panel
        mainPanel(
          uiOutput("mainPanelContent")
          )
        )
    )
)
