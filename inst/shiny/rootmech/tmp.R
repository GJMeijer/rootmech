server <- function(input, output, session) {
  # DATA SELECTION AND SETTINGS ####
  # Initiate reactive values
  val <- reactiveValues(
    df = data.frame(dr = 0, tru = 0),
    custom = custom_default,
    table = custom_default
  )
  # get fit type
  fittype <- reactive({
    df_fittype$label[df_fittype$name == input$fittype]
  })
  # fractions for confidence and prediction intervals
  confidence <- reactive({
    input$confidence/100
  })
  prediction <- reactive({
    input$prediction/100
  })

  # REFERENCE DATA SOURCE ####
  # DOI
  output$label_doi <- renderUI({
    if (input$dataset == user_data_name) {
      NULL
    } else {
      link <- df$link[df$label == input$dataset]
      tags$a(href = url, link, target = "_blank")
    }
  })

  # REACTIVE DATA ####
  # update data in table when different species selected
  observeEvent(input$dataset, {
    if (input$dataset == user_data_name) {
      val$table <- val$custom
    } else {
      val$table <- dm %>%
        filter(dataset_id == df$dataset_id[df$label == input$dataset])
      select(dr, tru)
    }
    val$df <- val$table %>% drop_na()
  })
  # update data to fit when custom table is altered + store custom data
  observeEvent(input$table, {
    if (!is.null(input$table) & input$dataset == user_data_name) {
      val$table <- as.data.frame(hot_to_r(input$table)) %>%
        rename(dr = "Diameter [mm]", tru = "Strength [MPa]")
      val$custom <- val$table
      val$df <- val$table %>% drop_na()
    }
  })

  # render datatable
  output$table <- renderRHandsontable({
    if (input$dataset == user_data_name) {
      rhandsontable(
        val$table %>% rename(
          "Diameter [mm]" = dr,
          "Strength [MPa]" = tru
        )
      )
    } else {
      rhandsontable(
        val$table %>% rename(
          "Diameter [mm]" = dr,
          "Strength [MPa]" = tru
        )
      ) %>%
        hot_col("Diameter [mm]", readOnly = TRUE) %>%
        hot_col("Strength [MPa]", readOnly = TRUE)
    }
  })

  # FITS ####
  # check if perfect fit - co-linearity of data in log-space
  perfect <- reactive({
    check_perfectfit(val$df$dr, val$df$tru)
  })
  # generate fits
  ft <- reactive({
    # check input
    validate(
      need(
        length(unique(val$df$dr)) >= 2,
        "Need at least two unique diameter values to generate fit"
      ),
      need(
        all(val$df$dr > 0),
        "All diameter values must be >0"
      ),
      need(
        all(val$df$tru > 0),
        "All strength values must be >0"
      )
    )
    # generate fits
    if (perfect()$perfect == TRUE) {
      perfect()$par %>%
        mutate(
          sd_multiplier = 0,
          sd_power = ifelse(fittype() == "normal_force", -2, 0),
          sdlog = 0,
          shape = Inf
        )
    } else if (fittype() == "weibull") {
      power_weibull_fit(val$df$dr, val$df$tru)
    } else if (fittype() == "normal_strength") {
      power_normal_fit(val$df$dr, val$df$tru, sd_power = 0)
    } else if (fittype() == "normal_force") {
      Ar <- pi/4*val$df$dr^2
      power_normal_fit(val$df$dr, Ar*val$df$tru, sd_power = 0) %>%
        mutate(
          multiplier = multiplier/(pi/4),
          power = power - 2,
          sd_multiplier = sd_multiplier/(pi/4),
          sd_power = sd_power - 2
        )
    } else if (fittype() == "normal_scaled") {
      ft <- power_normal_fit(val$df$dr, val$df$tru, sd_power = "scaled")
    } else if (fittype() == "lognormal_uncorrected") {
      ft <- power_lognormal_fit(val$df$dr, val$df$tru, correction = FALSE)
    } else if (fittype() == "lognormal") {
      ft <- power_lognormal_fit(val$df$dr, val$df$tru, correction = TRUE)
    }
  })
  # output: fit parameters
  output$multiplier <- renderUI(
    HTML(paste0(
      "Average strength power law, multiplier: ",
      "t", tags$sub("r,u,0"), " = ",
      numeric2character(ft()$multiplier, digits),
      " MPa"
    ))
  )
  output$power <- renderUI(
    HTML(paste0(
      "Average strength power law, power coefficient: ",
      "\u03b2", tags$sub("t"), " = ",
      numeric2character(ft()$power, digits)
    ))
  )
  output$fitpar1 <- renderUI(
    if (fittype() == "weibull") {
      HTML(paste0(
        "Weibull shape parameter: ",
        "\u03ba", tags$sub("t"), " = ",
        numeric2character(ft()$shape, digits)
      ))
    } else if (fittype() %in% c("normal_strength", "normal_scaled", "normal_force")) {
      HTML(paste0(
        "Standard deviation power law, multiplier: ",
        "\u03c3", tags$sub("t,0"), " = ",
        numeric2character(ft()$sd_multiplier, digits),
        " MPa"
      ))
    } else if (fittype() %in% c("lognormal_uncorrected", "lognormal")) {
      HTML(paste0(
        "Log-transformed standard deviation: ",
        "\u03c3", tags$sub("L"), " = ",
        numeric2character(ft()$sdlog, digits)
      ))
    }
  )
  output$fitpar2 <- renderUI(
    if (fittype() %in% c("normal_strength", "normal_scaled", "normal_force")) {
      HTML(paste0(
        "Standard deviation power law, power coefficient: ",
        "\u03b2", tags$sub("\u03c3,t"), " = ",
        numeric2character(ft()$sd_power, digits)
      ))
    } else {
      ""
    }
  )

  # ANDERSON-DARLING ####
  # predictions
  tru_pred <- reactive({
    ft()$multiplier*val$df$dr^ft()$power
  })
  # calculate Anderson-Darling confidence level
  pAD <- reactive({
    if (perfect()$perfect == TRUE) {
      NA
    } else if (fittype() == "weibull") {
      anderson_darling_weibull(
        val$df$tru/tru_pred(),
        shape = ft()$shape
      )$p
    } else if (fittype() %in% c("normal_strength", "normal_scaled", "normal_force")) {
      anderson_darling_normal(
        (val$df$tru - tru_pred())/(ft()$sd_multiplier*val$df$dr^ft()$sd_power),
        mu = 0,
        sd = 1
      )$p
    } else if (fittype() == "lognormal") {
      anderson_darling_normal(
        (log(val$df$tru) - log(ft()$multiplier) - power*log(val$df$dr) + 0.5*ft()$sdlog^2)/ft()$sdlog,
        mu = 0,
        sd = 1
      )$p
    } else if (fittype() == "lognormal_uncorrected") {
      anderson_darling_normal(
        (log(val$df$tru) - log(ft()$multiplier) - power*log(val$df$dr))/ft()$sdlog,
        mu = 0,
        sd = 1
      )$p
    }
  })
  # create label
  output$pAD <- renderText({
    paste0(
      "Anderson-Darling confidence level: ",
      "p", tags$sub("AD"), " = ",
      numeric2character(pAD(), digits)
    )
  })

  # CONFIDENCE INTERVALS  ####
  # calculate covariance matrix
  Sigma <- reactive({
    if (perfect()$perfect == TRUE) {
      matrix(Inf, nrow = 2, ncol = 2)
    } else if (fittype() == "weibull") {
      power_weibull_covariancematrix(
        val$df$dr, val$df$tru,
        ft()$multiplier, ft()$power, ft()$shape
      )
    } else if (fittype() %in% c("normal_strength", "normal_scaled", "normal_force")) {
      power_normal_covariancematrix(
        val$df$dr, val$df$tru,
        ft()$multiplier, ft()$power, ft()$sd_multiplier, ft()$sd_power
      )
    } else if (fittype() == "lognormal_uncorrected") {
      power_lognormal_covariancematrix(
        val$df$dr, val$df$tru,
        ft()$multiplier, ft()$power, ft()$sdlog,
        correction = FALSE
      )
    } else if (fittype() == "lognormal") {
      power_lognormal_covariancematrix(
        val$df$dr, val$df$tru,
        ft()$multiplier, ft()$power, ft()$sdlog,
        correction = TRUE
      )
    }
  })
  # generate fit lines and confidence intervals
  df_confidence <- reactive({
    xp = seq(min(val$df$dr), max(val$df$dr), l = 251)
    if (perfect()$perfect == TRUE) {
      tibble(
        x = xp,
        y = ft()$multiplier*x^ft()$power,
        ymin = y,
        ymax = y
      )
    } else {
      power_confidence(
        xp,
        Sigma(),
        ft()$multiplier, ft()$power,
        level = confidence()
      )
    }
  })

  # PREDICTION INTERVALS ####
  # generate prediction interval
  df_prediction <- reactive({
    x <- seq(min(val$df$dr), max(val$df$dr), l = 251)
    if (perfect()$perfect == TRUE) {
      y <- ft()$multiplier*x^ft()$power
      data.frame(x = x, y = y, ymin = y, ymax = y)
    } else if (fittype() == "weibull") {
      power_weibull_predictioninterval(
        x,
        ft()$multiplier, ft()$power, ft()$shape,
        level = prediction()
      )
    } else if (fittype() %in% c("normal_strength", "normal_scaled", "normal_force")) {
      power_normal_predictioninterval(
        x,
        ft()$multiplier, ft()$power, ft()$sd_multiplier, ft()$sd_power,
        level = prediction()
      )
    } else if (fittype() == "lognormal_uncorrected") {
      power_lognormal_predictioninterval(
        x, ft()$multiplier, ft()$power, ft()$sdlog,
        level = prediction(),
        correction = FALSE
      )
    } else if (fittype() == "lognormal") {
      power_lognormal_predictioninterval(
        x,
        ft()$multiplier, ft()$power, ft()$sdlog,
        level = prediction(),
        correction = TRUE
      )
    }
  })



  # PLOT ####
  # plot limits
  # xlim <- reactive({
  #   round_range(1.1*val$df$dr, lim = c(0, NA), ticks = 5)$lim
  # })
  # ylim <- reactive({
  #   round_range(1.1*val$df$tru, lim = c(0, NA), ticks = 5)$lim
  # })
  # generate plotly object with fit
  output$plot <- plotly::renderPlotly({
    validate(
      need(nrow(ft() > 0), "Plot requires valid fit")
    )
    plotly::plot_ly() %>%
      plotly::add_ribbons(    # prediction interval ribbon
        data = df_prediction(),
        name = paste0(input$prediction, "% Prediction interval"),
        x = ~x,
        ymin = ~ymin,
        ymax = ~ymax,
        line = list(width = 0)
      ) %>%
      plotly::add_ribbons(    # confidence interval ribbon
        data = df_confidence(),
        name = paste0(input$confidence, "% Confidence interval of fit"),
        x = ~x,
        ymin = ~ymin,
        ymax = ~ymax,
        line = list(width = 0)
      ) %>%
      plotly::add_trace(        # add power law fit line
        data = df_confidence(),
        name = "Power law fit",
        x = ~x,
        y = ~y,
        type = "scatter",
        mode = "lines"
      ) %>%
      plotly::add_trace(      # add points
        data = val$df,
        name = "Data points",
        x = ~dr,
        y = ~tru,
        type = "scatter",
        mode = "markers"
      ) %>%
      plotly::layout(      # add layout
        xaxis = list(
          range = round_range(1.1*val$df$dr, lim = c(0, NA), ticks = 5)$lim,
          title = "Root diameter [mm]",
          rangemode = "tozero"
        ),
        yaxis = list(
          range = round_range(1.1*val$df$tru, lim = c(0, NA), ticks = 5)$lim,
          title = "Tensile strength [MPa]",
          rangemode = "tozero"
        )
      )
  })
}
