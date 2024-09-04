ui <- shiny::navbarPage(
  title = "Power law fitting of root tensile strength data",
  position = "fixed-top",
  collapsible = TRUE,

  # TAB: FIT ####
  shiny::tabPanel(
    "Fit custom data",
    shiny::tags$style(type="text/css", "body {padding-top: 70px;}"),

    shiny::sidebarLayout(
      shiny::sidebarPanel(
        # settings - fitting method
        shiny::selectInput(
          "fittype",
          "Power law fit model",
          choices = df_fittype$name,
          selected = df_fittype$name[df_fittype$label == "gamma"],
          multiple = FALSE,
          selectize = FALSE
        ),
        # settings - weighting exponent
        shiny::sliderInput(
          "weighting_exponent",
          "Power law weighting exponent",
          min = -4,
          max = 4,
          value = 0,
          step = 0.1
        ),
        # settings - confidence level (slider)
        shiny::sliderInput(
          "confidence",
          "Confidence level of fit [%]",
          min = 0,
          max = 100,
          value = 95,
          step = 1
        ),
        # settings - prediction interval (slider)
        shiny::sliderInput(
          "prediction",
          "Prediction interval level [%]",
          min = 0,
          max = 100,
          value = 95,
          step = 1
        ),
        # editable dable
        rhandsontable::rHandsontableOutput("table")
      ),
      shiny::mainPanel(
        # render results
        shiny::uiOutput("label_multiplier"),
        shiny::uiOutput("label_exponent"),
        shiny::uiOutput("label_intradiameter1"),
        shiny::uiOutput("label_intradiameter2"),
        shiny::textOutput("label_loglikelihood"),
        shiny::textOutput("label_ks"),
        # results table
         # shiny::tableOutput("table_results"),
        # plot power law
        plotly::plotlyOutput("plot_fit"),
        # plot KS
        plotly::plotlyOutput("plot_ks")
      )
    )
  ),

  # TAB: POWER LAW DATA ####
  shiny::tabPanel(
    "Fit results",
    shiny::sidebarLayout(
      shiny::sidebarPanel(
        # settings - fitting method
        shiny::selectInput(
          "fittype2",
          "Power law fit model",
          choices = df_fittype$name,
          selected = df_fittype$name[df_fittype$label == "gamma"],
          multiple = FALSE,
          selectize = FALSE
        ),
        # functional group
        shiny::selectInput(
          "species_group",
          "Plant group",
          choices = unique(fit_opts$functional_group2),
          selected = "Conifers",
          multiple = TRUE,
          selectize = TRUE
        ),
        # Family - selectize (default - all)
        shiny::selectInput(
          "species_family",
          "Plant family",
          choices = unique(fit_opts$family),
          selected = unique(fit_opts$family)[1],
          multiple = TRUE,
          selectize = TRUE
        ),
        # Species
        shiny::selectInput(
          "species_species",
          "Species",
          choices = unique(fit_opts$species),
          selected = unique(fit_opts$species)[1],
          multiple = TRUE,
          selectize = TRUE
        ),
        # Log plot
        shiny::strong("Plot settings"),
        shiny::checkboxInput(
          "logplot",
          "Logarithmic axes",
          value = TRUE
        )
      ),
      shiny::mainPanel(
        # plot
        plotly::plotlyOutput("plot_powerlaw"),
        # table
        shiny::br(),
        DT::DTOutput("table_powerlaw") #shiny::tableOutput
      )
    )
  ),


  # TAB: DOCUMENTATION ####
  shiny::tabPanel(
    "Documentation",
    shiny::tags$style(type="text/css", "body {padding-top: 70px;}"),
    shiny::withMathJax(),
    shiny::includeMarkdown("./www/documentation.rmd")
  ),

  footer = shiny::fluidRow(
    shiny::hr(),
    shiny::column(
      4,
      "G. J. Meijer",
      br(),
      "University of Bath",
      br(),
      tags$a(href = "mailto:gjm36@bath.ac.uk", "gjm36@bath.ac.uk")
    ),
    shiny::column(
      4,
      "September 2024"
    ),
    shiny::column(
      4,
      "Application developed as part of publication",
      tags$a(href = "test", "Improving power law fitting of root tensile strength–diameter relationships"),
      ". ",
      "Please consider citing this source when using this tool"
    )
  )

  # "January 2024 - G.J. Meijer - University of Bath - gjm36@bath.ac.uk - Application built as part of publication x")
)
