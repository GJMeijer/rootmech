ui <- shiny::navbarPage(
  title = "Power law fitting of root strength data",
  position = "fixed-top",
  collapsible = TRUE,

  # FIRST TAB
  shiny::tabPanel(
    "Analysis",
    tags$style(type="text/css", "body {padding-top: 70px;}"),

    sidebarLayout(
      sidebarPanel(
        # settings - fitting method
        selectInput(
          "fittype",
          "Power law fit type",
          choices = df_fittype$name,
          selected = df_fittype$name[df_fittype$label == "weibull"],
          multiple = FALSE,
          selectize = FALSE
        ),
        # settings - weighting exponent
        sliderInput(
          "weighting_exponent",
          "Power law weighting exponent",
          min = -3,
          max = 3,
          value = 0,
          step = 0.1
        ),
        # settings - confidence level (slider)
        sliderInput(
          "confidence",
          "Confidence level of fit [%]",
          min = 0,
          max = 100,
          value = 95,
          step = 1
        ),
        # settings - prediction interval (slider)
        sliderInput(
          "prediction",
          "Prediction interval level [%]",
          min = 0,
          max = 100,
          value = 95,
          step = 1
        ),
        # data select
        selectInput(
          "dataset",
          "Data",
          choices = dataset_options,
          selected = dataset_options[10],
          multiple = FALSE,
          selectize = FALSE
        ),
        # DOI of data
        htmlOutput("label_doi"),
        # editable dable
        rHandsontableOutput("table")
      ),
      mainPanel(
        # plot
        plotlyOutput("plot_fit"),
        tableOutput("aa"),
        #tableOutput("bb"),
        plotlyOutput("plot_ks"),
        # show fitting parameters
        br(),
        htmlOutput("multiplier"),
        htmlOutput("power"),
        htmlOutput("fitpar1"),
        htmlOutput("fitpar2"),
        htmlOutput("pAD")
      )
    )
  ),

  shiny::tabPanel(
    "Documentation",
    tags$style(type="text/css", "body {padding-top: 70px;}"),
    shiny::withMathJax(),
    shiny::includeMarkdown('documentation.rmd')
  ),

  footer = fluidRow(
    hr(),
    column(
      4,
      "G. J. Meijer",
      br(),
      "University of Bath",
      br(),
      tags$a(href = "mailto:gjm36@bath.ac.uk", "gjm36@bath.ac.uk")
    ),
    column(
      4,
      "January 2024"
    ),
    column(
      4,
      "Application developed as part of publication",
      tags$a(href = "test", "Improving power law fitting"),
      ". ",
      "Please consider citing this source when using this tool"
    )
  )

  # "January 2024 - G.J. Meijer - University of Bath - gjm36@bath.ac.uk - Application built as part of publication x")
)
