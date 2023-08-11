#v0.1 - 20210416 - first working version

################
### DEFINE UI ###
#################

#make UI
ui <- shiny::navbarPage(
  title = "Root fitting",
  position = "fixed-top",
  collapsible = TRUE,

  # Biomech: data tab
  shiny::tabPanel(
    "Biomechanics: data",
    tags$style(type="text/css", "body {padding-top: 70px;}"),
    fluidRow(
      fileInput(
        "upload",
        "Upload a file",
        accept = c(".csv", "text/csv", "text/comma-separated-values"),
      ),
      numericInput(
        "coldr",
        "Diameter column",
        1,
        min = 1,
        max = 1000,
        step = 1
      ),
      numericInput(
        "colepsru",
        "Strain column",
        2,
        min = 1,
        max = 1000,
        step = 1
      ),
      numericInput(
        "coltru",
        "Strength column",
        3,
        min = 1,
        max = 1000,
        step = 1
      )
    ),
    rHandsontableOutput("datatable")
  ),

  # Fit tab
  shiny::tabPanel(
    "Biomechanics: fit",
    tags$style(type="text/css", "body {padding-top: 70px;}"),
    radioButtons(
      "weightpowerlaw",
      "Power-law fit weighting",
      choiceNames = c("Equal", "Root cross-sectional area"),
      choiceValues = c("equal", "area"),
      selected = "area"
    ),
    radioButtons(
      "weightcopula",
      "Copula weighting",
      choiceNames = c("Equal", "Root cross-sectional area"),
      choiceValues = c("equal", "area"),
      selected = "area"
    ),
    rHandsontableOutput("fittable"),
    plotOutput("plot_drepsru"),
    plotOutput("plot_drtru"),
    plotOutput("plot_epsrutru"),
  )
)
