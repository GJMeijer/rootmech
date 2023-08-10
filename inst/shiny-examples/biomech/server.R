server <- function(input, output) {

  # reactive values
  values <- reactiveValues(data = data_initial)

  observeEvent(input$upload, {
    tmp <- read.csv(input$upload$datapath)
    values$data <- tmp[, c(input$coldr, input$colepsru, input$coltru)]
  })

  # update data in reactive values if table changes
  observe({
    if(!is.null(input$datatable)){
      values$data <- as.data.frame(hot_to_r(input$datatable))
      output$datatable <- renderRHandsontable({
        rhandsontable(
          values$data,
          colHeaders = c("Root diameter [mm]", "Tensile strain to peak [mm/mm]", "Tensile strength [MPa]")
        )
      })
    }
  })

  # generate rhandsontable object
  output$datatable <- renderRHandsontable({
    rhandsontable(
      values$data,
      colHeaders = c("Root diameter [mm]", "Tensile strain to peak [mm/mm]", "Tensile strength [MPa]")
    )
  })

  # generate biomechanical fitting parameters
  bmfit <- reactive({
    biomech_fit(
      values$data$dr,
      values$data$epsru,
      values$data$tru,
      weights_power = if (input$weightpowerlaw == "equal") {
        rep(1, length(values$data$dr))
      } else {
        values$data$dr^2
      },
      weights_copula = if (input$weightcopula == "equal") {
        rep(1, length(values$data$dr))
      } else {
        values$data$dr^2
      }
    )
  })

  # generate plot - dr-epsru
  lab_drepsru <- reactive({
    create_labels(
      c(bmfit()$epsru0, bmfit()$betaeps, bmfit()$kappaeps),
      prefix = c(
        "epsilon[r*','*u*','*0]==",
        "beta[epsilon]==",
        "kappa[epsilon]=="
      ),
      suffix = c("~'mm/mm'", "", ""),
      parse = TRUE
    )
  })
  output$plot_drepsru <- renderPlot({
    power_weibull_plot(
      values$data$dr,
      values$data$epsru,
      multiplier = bmfit()$epsru0,
      power = bmfit()$betaeps,
      shape = bmfit()$kappaeps,
      ylab = expression("Ultimate tensile strain"~epsilon[r*','*u]~"[mm/mm]"),
      labels = lab_drepsru()
    )
  })

  # generate plot - dr-tru
  lab_drtru <- reactive({
    create_labels(
      c(bmfit()$tru0, bmfit()$betat, bmfit()$kappat),
      prefix = c(
        "t[r*','*u*','*0]==",
        "beta[t]==",
        "kappa[t]=="
      ),
      suffix = c("~'MPa'", "", ""),
      parse = TRUE
    )
  })
  output$plot_drtru <- renderPlot({
    power_weibull_plot(
      values$data$dr,
      values$data$tru,
      multiplier = bmfit()$tru0,
      power = bmfit()$betat,
      shape = bmfit()$kappat,
      labels = lab_drtru()
    )
  })

  # normalised data
  df_copula <- reactive({
    data.frame(
      x = values$data$epsru/(bmfit()$epsru0*values$data$dr^bmfit()$betaeps),
      y = values$data$tru/(bmfit()$tru0*values$data$dr^bmfit()$betat)
    )
  })
  # generate plot - copula
  lab_epsrutru <- reactive({
    create_labels(
      bmfit()$rho,
      prefix = "rho==",
      suffix = "",
      parse = TRUE
    )
  })
  output$plot_epsrutru <- renderPlot({
    copula_gaussian_plot(
      df_copula()$x,
      df_copula()$y,
      bmfit()$kappaeps,
      bmfit()$kappat,
      rho = bmfit()$rho,
      labels = lab_epsrutru()
    )
  })

  # data table with fitting results
  output$fittable <- renderRHandsontable({
    rhandsontable(
      bmfit(),
      colHeaders = c("epsru0 [mm/mm]", "betaeps", "kappaeps",
                     "tru0 [MPa]", "betat", "kappat",
                     "rho")
    )
  })
}
