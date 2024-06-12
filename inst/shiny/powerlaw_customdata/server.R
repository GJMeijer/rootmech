server <- function(input, output, session) {
  # DATA SELECTION AND SETTINGS ####
  # Initiate reactive values
  val <- shiny::reactiveValues(
    df = df_maize
  )
  # get fit type
  fittype <- shiny::reactive({
    df_fittype$label[df_fittype$name == input$fittype]
  })
  # fractions for confidence and prediction intervals
  confidence <- shiny::reactive({
    input$confidence/100
  })
  prediction <- shiny::reactive({
    input$prediction/100
  })

  # REACTIVE DATA ####
  # update data to fit when custom table is altered + store custom data
  shiny::observeEvent(input$table, {
    tmp <- as.data.frame(rhandsontable::hot_to_r(input$table))
    colnames(tmp) <- c("dr", "tru", "weight")
    val$df <- tmp
  })
  # change weights if weight slider used
  shiny::observeEvent(input$weighting_exponent, {
    w <- val$df$dr^input$weighting_exponent
    val$df$weight <- w/sum(w)*nrow(val$df)
  })
  # render datatable
  output$table <- shiny::debounce(
    rhandsontable::renderRHandsontable({
      tmp <- val$df
      colnames(tmp) <- c("Diameter [mm]", "Strength [MPa]", "Weight")
      rhandsontable::rhandsontable(
        tmp,
        rowHeaders = NULL,
        useTypes = TRUE
      )
    }),
    500
  )

  # FITS ####
  # data to use
  df_use <- shiny::reactive({
    # check input
    shiny::validate(
      shiny::need(
        length(unique(val$df$dr)) >= 3,
        "Need at least three unique diameter values"
      ),
      shiny::need(
        all(val$df$dr > 0),
        "All diameter values must be >0"
      ),
      shiny::need(
        all(val$df$tru > 0),
        "All strength values must be >0"
      ),
      shiny::need(
        all(val$df$weight >= 0),
        "All weights values must be >=0"
      )
    )
    # filter data - only use data
    val$df[(val$df$dr > 0) & (val$df$tru > 0) & (val$df$weight > 0), ]
  })
  # generate fits
  ft <- shiny::reactive({
    # check if perfect fit - co-linearity of data in log-space
    perfect <- check_perfectfit(df_use()$dr, df_use()$tru)
    # generate fits
    if (perfect$colinear == TRUE) {
      ft <- perfect
    } else {
      ft <- powerlaw_fit_iterativeweighting(
        df_use()$dr, df_use()$tru, fittype(), weights = df_use()$weight
      )
      ft$colinear <- FALSE
    }
    # add diameter range
    ft$drmin <- min(df_use()$dr, na.rm = TRUE)
    ft$drmax <- max(df_use()$dr, na.rm = TRUE)
    # return
    ft
  })

  # PREDICTION INTERVAL ####
  df_pred <- shiny::reactive({
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
  covariance <- shiny::reactive({
    if (ft()$colinear == TRUE) {
      matrix(
        0,
        nrow = 2, ncol = 2,
        dimnames = list(c("multiplier", "exponent"), c("multiplier", "exponent"))
      )
    } else {
      n <- length(df_use()$dr)
      powerlaw_covariancematrix(
        df_use()$dr,
        df_use()$tru,
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
        weights = df_use()$weight
      )
    }
  })

  # CONFIDENCE INTERVAL ####
  df_conf <- shiny::reactive({
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
  plot1a <- shiny::reactive({
    plt <- plotly::plot_ly()
    plt <- plotly::add_markers(
      plt,
      type = "scatter",
      mode = "markers",
      x = df_use()$dr,
      y = df_use()$tru,
      name = "Data"
    )
    plotly::layout(
      plt,
      xaxis = list(
        title = "Root diameter [mm]",
        range = round_range(1.05*df_use()$dr, lim = c(0, NA))$lim
      ),
      yaxis = list(
        title = "Tensile strength [MPa]",
        range = round_range(1.05*df_use()$tru, lim = c(0, NA))$lim
      ),
      legend = list(
        title = list(
          text = NA
        )
      )
    )
  })
  ## add best fit ####
  plot1b <- shiny::reactive({
    xp <- seq(min(df_use()$dr, na.rm = TRUE), max(df_use()$dr, na.rm = TRUE), l = 101)
    yp <- ft()$multiplier*xp^ft()$exponent
    plotly::add_trace(
      plot1a(),
      type = "scatter",
      mode = "lines",
      x = xp,
      y = yp,
      name = "Fit"
    )
  })
  ## add prediction interval ####
  plot1c <- shiny::reactive({
    plotly::add_ribbons(
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
  plot1d <- shiny::reactive({
    plotly::add_ribbons(
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
  output$plot_fit <- plotly::renderPlotly({
    plot1d()
  })

  # GENERATE KS PLOT ####
  # Get all KS information
  ks <- shiny::reactive({
    shiny::validate(
      shiny::need(
        ft()$colinear == FALSE,
        "Zero intra-diameter variation, KS = 0"
      )
    )
    powerlaw_ks(
      df_use()$dr,
      df_use()$tru,
      fittype(),
      ft()$multiplier,
      ft()$exponent,
      sd_multiplier = if ("sd_multiplier" %in% names(ft())) ft()$sd_multiplier else NULL,
      sd_exponent = if ("sd_exponent" %in% names(ft())) ft()$sd_exponent else NULL,
      sdlog = if ("sdlog" %in% names(ft())) ft()$sdlog else NULL,
      shape = if ("shape" %in% names(ft())) ft()$shape else NULL,
      scale = if ("scale" %in% names(ft())) ft()$scale else NULL,
      width = if ("width" %in% names(ft())) ft()$width else NULL,
      weights = df_use()$weight
    )
  })
  output$table_results <- shiny::renderTable({
    dfout <- as.data.frame(ft())
    dfout <- dfout[, !(names(dfout) %in% c("colinear", "drmin", "drmax"))]
    dfout
  })
  # generate a plot
  output$plot_ks <- plotly::renderPlotly({
    plt <- plotly::plot_ly()
    plt <- plotly::add_trace(
      plt,
      type = "scatter",
      mode = "lines",
      x = ks()$data$x,
      y = ks()$data$y,
      name = "Data"
    )
    plt <- plotly::add_trace(
      plt,
      type = "scatter",
      mode = "lines",
      x = ks()$fit$x,
      y = ks()$fit$y,
      name = "Fit"
    )
    plt <- plotly::add_trace(
      plt,
      type = "scatter",
      mode = "lines+markers",
      x = ks()$difference$x,
      y = ks()$difference$y,
      name = "KS distance"
    )
    plotly::layout(
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

  # GENERATE FIT RESULTS LABELS ####
  # KS distance
  output$label_ks <- renderText(paste0(
    "Kolmogorov-Smirnov distance, KS = ",
    numeric2character(ft()$ks_distance, digits = digits)
  ))
  # Loglikelihood
  output$label_loglikelihood <- renderText(paste0(
    "Loglikelihood = ",
    numeric2character(ft()$loglikelihood, digits = digits)
  ))
  # power law multiplier
  output$label_multiplier <- renderUI(shiny::HTML(paste0(
    "Power law multiplier, t<sub>r,u,0</sub> = ",
    numeric2character(ft()$multiplier, digits = digits),
    " MPa"
  )))
  # power law exponent
  output$label_exponent <- renderUI(shiny::HTML(paste0(
    "Power law exponent, &#946;<sub>t</sub> = ",
    numeric2character(ft()$exponent, digits = digits)
  )))
  # variation parameter 1
  output$label_intradiameter1 <- renderUI({
    if (fittype() == "gumbel") {
      shiny::HTML(paste0(
        "Scale parameter, &#952;<sub>0</sub> = ",
        numeric2character(ft()$scale, digits = digits),
        " MPa"
      ))
    } else if (fittype() == "gamma") {
      shiny::HTML(paste0(
        "Shape parameter, k<sub>t</sub> = ",
        numeric2character(ft()$shape, digits = digits)
      ))
    } else if (fittype() == "logistic") {
      shiny::HTML(paste0(
        "Scale parameter, s<sub>0</sub> = ",
        numeric2character(ft()$scale, digits = digits),
        " MPa"
      ))
    } else if (fittype() %in% c("lognormal", "lognormal_uncorrected")) {
      shiny::HTML(paste0(
        "Log-standard deviation, &#963;<sub>L</sub> = ",
        numeric2character(ft()$sdlog, digits = digits)
      ))
    } else if (fittype() %in% c("normal_strength", "normal_force", "normal_scaled")) {
      shiny::HTML(paste0(
        "Standard deviation power law multiplier, &#963;<sub>0</sub> = ",
        numeric2character(ft()$sd_multiplier, digits = digits),
        " MPa"
      ))
    } else if (fittype() == "uniform") {
      shiny::HTML(paste0(
        "Width parameter, w<sub>t</sub> = ",
        numeric2character(ft()$width, digits = digits),
        " MPa"
      ))
    } else if (fittype() == "weibull") {
      shiny::HTML(paste0(
        "Shape parameter, &#954;<sub>t</sub> = ",
        numeric2character(ft()$shape, digits = digits)
      ))
    } else {
      shiny::HTML(NULL)
    }
  })
  # variation paramter 2
  output$label_intradiameter2 <- renderUI({
    if (fittype() %in% c("normal_strength", "normal_force", "normal_scaled")) {
      shiny::HTML(paste0(
        "Standard deviation power law exponent, &#946;<sub>t</sub> = ",
        numeric2character(ft()$sd_exponent, digits = digits)
      ))
    } else {
      shiny::HTML(NULL)
    }
  })
}
