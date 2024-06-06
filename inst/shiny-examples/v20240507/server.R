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
      tags$a(
        href = df$link[df$label == input$dataset],
        df$link_label[df$label == input$dataset],
        target = "_blank"
      )
    }
  })

  # REACTIVE DATA ####
  # update data in table when different species selected
  observeEvent(input$dataset, {
    if (input$dataset == user_data_name) {
      val$table <- val$custom
    } else {
      d_id <- df$dataset_id[df$label == input$dataset]
      val$table <- dm %>%
        filter(dataset_id == d_id) %>%
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
  # weights
  w <- reactive({
    n <- length(val$df$dr)
    weights <- val$df$dr^input$weighting_exponent
    weights*n/sum(weights)
  })
  # generate fits
  ft <- reactive({
    # check input
    validate(
      need(
        length(unique(val$df$dr)) >= 3,
        "Need at least three unique diameter values"
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
    # check if perfect fit - co-linearity of data in log-space
    perfect <- check_perfectfit(val$df$dr, val$df$tru)
    # generate fits
    if (perfect$colinear == TRUE) {
      ft <- perfect
    } else {
      ft <- powerlaw_fit_iterativeweighting(
        val$df$dr, val$df$tru, fittype(), weights = w()
      )
      ft$colinear <- FALSE
      ft
    }
    # add diameter range
    ft$drmin <- min(val$df$dr, na.rm = TRUE)
    ft$drmax <- max(val$df$dr, na.rm = TRUE)
    # return
    ft
  })

  # PREDICTION INTERVAL ####
  df_pred <- reactive({
    xp <- seq(ft()$drmin, ft()$drmax, l = 101)
    if (ft()$colinear == TRUE) {
      yp <- ft()$multiplier*xp^ft()$exponent
      data.frame(
        x = xp,
        y = yp,
        ymin = yp,
        ymax = yp
      )
    } else {
      powerlaw_predictioninterval(
        xp,
        fittype(),
        ft()$multiplier,
        ft()$exponent,
        sd_multiplier = if ("sd_multiplier" %in% names(ft())) ft()$sd_multiplier else NULL,
        sd_exponent = if ("sd_exponent" %in% names(ft())) ft()$sd_exponent else NULL,
        sdlog = if ("sdlog" %in% names(ft())) ft()$sdlog else NULL,
        shape = if ("shape" %in% names(ft())) ft()$shape else NULL,
        scale = if ("scale" %in% names(ft())) ft()$scale else NULL,
        width = if ("width" %in% names(ft())) ft()$width else NULL,
        level = prediction()
      )
    }
  })

  # COVARIANCE MATRIX ####
  covariance <- reactive({
    if (ft()$colinear == TRUE) {
      matrix(
        0,
        nrow = 2, ncol = 2,
        dimnames = list(c("multiplier", "exponent"), c("multiplier", "exponent"))
      )
    } else {
      n <- length(val$df$dr)
      weights <- val$df$dr^input$weighting_exponent
      w <- weights*n/sum(weights)
      powerlaw_covariancematrix(
        val$df$dr,
        val$df$tru,
        fittype(),
        method = "fisher",
        multiplier = ft()$multiplier,
        exponent = ft()$exponent,
        sd_multiplier = if ("sd_multiplier" %in% names(ft())) ft()$sd_multiplier else NULL,
        sd_exponent = if ("sd_exponent" %in% names(ft())) ft()$sd_exponent else NULL,
        sdlog = if ("sdlog" %in% names(ft())) ft()$sdlog else NULL,
        shape = if ("shape" %in% names(ft())) ft()$shape else NULL,
        scale = if ("scale" %in% names(ft())) ft()$scale else NULL,
        width = if ("width" %in% names(ft())) ft()$width else NULL,
        n = 100,
        weights = w()
      )
    }
  })

  # CONFIDENCE INTERVAL ####
  df_conf <- reactive({
    xp <- seq(ft()$drmin, ft()$drmax, l = 101)
    if (ft()$colinear == TRUE) {
      yp <- ft()$multiplier*xp^ft()$exponent
      data.frame(
        x = xp,
        y = yp,
        ymin = yp,
        ymax = yp
      )
    } else {
      powerlaw_confidenceinterval(
        xp,
        ft()$multiplier,
        ft()$exponent,
        covariance(),
        level = confidence()
      )
    }
  })

  # GENERATE PLOTLY FIT PLOT ####
  ## plotly plot with datapoints ####
  plot1a <- reactive({
    plt <- plot_ly()
    plt <- add_markers(
      plt,
      type = "scatter",
      mode = "markers",
      x = val$df$dr,
      y = val$df$tru,
      name = "Data"
    )
    layout(
      plt,
      xaxis = list(
        title = "Root diameter [mm]",
        range = round_range(1.05*val$df$dr, lim = c(0, NA))$lim
      ),
      yaxis = list(
        title = "Tensile strength [MPa]",
        range = round_range(1.05*val$df$tru, lim = c(0, NA))$lim
      ),
      legend = list(
        title = list(
          text = NA
        )
      )
    )
  })
  ## add best fit ####
  plot1b <- reactive({
    xp <- seq(min(val$df$dr, na.rm = TRUE), max(val$df$dr, na.rm = TRUE), l = 101)
    yp <- ft()$multiplier*xp^ft()$exponent
    add_trace(
      plot1a(),
      type = "scatter",
      mode = "lines",
      x = xp,
      y = yp,
      name = "Fit"
    )
  })
  ## add prediction interval ####
  plot1c <- reactive({
    add_ribbons(
      plot1b(),
      x = df_pred()$x,
      ymin = df_pred()$ymin,
      ymax = df_pred()$ymax,
      name = paste0(input$prediction, "% Prediction interval"),
      line = list(
        width = 0
      ),
      opacity = 0.5
    )
  })
  ## add confidence interval ####
  plot1d <- reactive({
    add_ribbons(
      plot1c(),
      x = df_conf()$x,
      ymin = df_conf()$ymin,
      ymax = df_conf()$ymax,
      name = paste0(input$confidence, "% Confidence interval"),
      line = list(
        width = 0
      ),
      opacity = 0.5
    )
  })
  ## render final plot ####
  output$plot_fit <- renderPlotly({
    plot1d()
  })

  # GENERATE KS PLOT ####
  # Get all KS information
  ks <- reactive({
    validate(
      need(
        ft()$colinear == FALSE,
        "Zero intra-diameter variation, KS = 0"
      )
    )
    powerlaw_ks(
      val$df$dr,
      val$df$tru,
      fittype(),
      ft()$multiplier,
      ft()$exponent,
      sd_multiplier = if ("sd_multiplier" %in% names(ft())) ft()$sd_multiplier else NULL,
      sd_exponent = if ("sd_exponent" %in% names(ft())) ft()$sd_exponent else NULL,
      sdlog = if ("sdlog" %in% names(ft())) ft()$sdlog else NULL,
      shape = if ("shape" %in% names(ft())) ft()$shape else NULL,
      scale = if ("scale" %in% names(ft())) ft()$scale else NULL,
      width = if ("width" %in% names(ft())) ft()$width else NULL,
      weights = w()
    )
  })
  output$aa <- renderTable({
    as.data.frame(ft())
  })
  output$bb <- renderTable({
    length(ks())
  })
  # generate a plot
  output$plot_ks <- renderPlotly({
    plt <- plot_ly()
    plt <- add_trace(
      plt,
      type = "scatter",
      mode = "lines",
      x = ks()$data$x,
      y = ks()$data$y,
      name = "Data"
    )
    plt <- add_trace(
      plt,
      type = "scatter",
      mode = "lines",
      x = ks()$fit$x,
      y = ks()$fit$y,
      name = "Fit"
    )
    plt <- add_trace(
      plt,
      type = "scatter",
      mode = "lines+markers",
      x = ks()$difference$x,
      y = ks()$difference$y,
      name = "KS distance"
    )
    layout(
      plt,
      xaxis = list(
        title = "Observations"
      ),
      yaxis = list(
        title = "Cumulative density",
        range = c(0, 1)
      ),
      legend = list(
        title = list(
          text = NA
        )
      )
    )
  })
}
