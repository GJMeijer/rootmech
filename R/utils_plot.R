# Functions for generating nice plots for data+fits
# G. J. Meijer, 08/08/2023


#' Add a label in the corner of an existing ggplot
#'
#' @description
#' Add a text label an existing ggplot object. For example. use to annotate
#' a subplot indicator to an exsiting plot
#'
#' @param plt ggplot object
#' @param label vector with character strings to add
#' @param location two-part vector indicating the relative location in the
#'   plot for placing labels. Follows same logic as legend placement in ggplot.
#' @param justification two-part vector indicating the horizontal and
#'   vertical justification of the labels. Follwos the same logic as legend
#'   justification in ggplot
#' @param label_spacing if `label` consists of multiple elements, the relative
#'   `y`-spacing between each element
#' @param size plot size of text
#' @param color plot color of text
#' @param parse if `TRUE`, labels are parsed when plotted
#' @return ggplot object with added annotation
#'
annotate_label <- function(
    plt,
    label,
    location = c(0.02, 0.98),
    justification = c(0, 1),
    label_spacing = 0.10,
    size = 3,
    color = "black",
    parse = FALSE
) {
  # get current ggplot limits
  xlim <- plt$coordinates$limits$x
  ylim <- plt$coordinates$limits$y
  # plus or minus
  pm <- 1 - 2*justification
  # get coordinates for plotting
  x <- xlim[1] + location[1]*diff(xlim)
  y <- ylim[1] + (location[2] + pm[2]*(label_spacing*(seq(length(label)) - 1)))*diff(ylim)
  # add annotation and return
  plt + ggplot2::annotate(
    "text",
    x = x,
    y = y,
    label = label,
    hjust = justification[1],
    vjust = justification[2],
    color = color,
    size = size,
    parse = parse
  )
}


#' Determine rounded plot axis limits
#'
#' @description
#' Determine plot axis limits based on input values, preferred tick step sizes
#' and a maximum number of ticks
#'
#' @md
#' @param x values on axis that should be within the plot range
#' @param lim two-parameter vector with fixed axis limits. If `NA`, an axis
#'   limit can be freely chosen
#' @param step array with possible step sizes. These values are relative to
#'   the order of magnitude of the difference in input `x`, i.e. relative to
#'   `10^floor(log10(max(x) - min(x)))`
#' @param ticks the maximum number of ticks that is acceptable
#' @return a list with field:
#'
#'   * `lim`: two-parameter vector with min and max limit of the plot range
#'   * `breaks`: vector with plot breakpoints
#'
#' @examples
#' # data
#' x <- 10*rnorm(100)
#'
#' # basic range
#' z <- round_range(x)
#' plot(x, seq(length(x)), xlim = z$lim, xaxs = "i", xaxt = "n")
#' axis(side = 1, at = z$breaks)
#'
#' # fix the lower bound
#' z <- round_range(x, lim = c(-100, NA))
#' plot(x, seq(length(x)), xlim = z$lim, xaxs = "i", xaxt = "n")
#' axis(side = 1, at = z$breaks)
#'
round_range <- function(
    x,
    lim = c(NA, NA),
    step = c(0.1, 0.2, 0.5, 1, 2, 5),
    ticks = 7
) {
  # get all limits
  lims <- c(
    min(min(x), lim[1], na.rm = TRUE),
    max(max(x), lim[2], na.rm = TRUE)
  )
  # order of magnitude of difference
  oom <- floor(log10(diff(lims)))
  # round scaled upper and lower values using potential step sizes (step)
  z0 <- floor((lims[1]/10^oom)/step)*step
  z1 <- ceiling((lims[2]/10^oom)/step)*step
  # calculate number of ticks with distance <step> required for each option
  ntick <- ceiling((z1 - z0)/step)
  ntick[ntick > ticks] <- 0
  # pick option that has the most ticks but number still smaller than <n>
  i <- which.max(ntick)
  # set range
  range <- c(z0[i], z1[i])*10^oom
  # substitute fixed values
  range[!is.na(lim)] <- lim[!is.na(lim)]
  # return range and ticks
  list(
    lim = range,
    breaks = seq(range[1], range[2], step[i]*10^oom)
  )
}


#' Round a double to a number of significant digits and convert to character
#'
#' @description
#' Round a number to a number of significant digits and convert to character
#' string. Always returns `digits` number of significant digits, also when
#' there are trailing zeros
#'
#' @param x numeric value to convert to character string
#' @param digits number of significant digits
#' @return character string
#' @examples
#' numeric2character(pi, digits = 4)
#' numeric2character(sqrt(5)/100000, digits = 3)
#'
numeric2character <- function(x, digits = 3) {
  gsub(
    "\\.$", "",
    formatC(
      signif(x, digits = digits),
      digits = digits,
      format = "fg",
      flag = "#"
    )
  )
}


#' Create annotation labels for plotting
#'
#' @description
#' Creates character strings for annotating a plot. Combines a prefix, a
#' numeric value (will be converted to a string using the `numeric2character()`
#' function) and a character suffix.
#'
#' @param values array with values to be rounded
#' @param prefix,suffix arrays with character prefix and suffixes
#' @param parse if `TRUE`, the combined character string is assumed as an
#'   parseable expression
#' @param sep seperator to put between prefix, value and suffixes
#' @param digits number of significant digits to use when rounding values
#' @param x,y relative location of label in the plot (0 <= x,y <= 1)
#' @param hjust,vjust horizontal and vertical justification of labels in the
#'   plot
#' @return dataframe with fields for labels (`label`), relative position in
#'   plot (`x`, `y`), text justification (`hjust`, `vjust`) and whether the
#'   labels should be parsed or not (`parse`)
#' @examples
#' create_annotations(
#'   c(0.23, -0.1),
#'   prefix = c("a=", ""),
#'   suffix = c(" mm", " N"),
#'   parse = FALSE,
#'   digits = 4
#' )
#'
create_annotations <- function(
    values,
    prefix = NULL,
    suffix = NULL,
    parse = TRUE,
    sep = "",
    digits = 3,
    x = 0.98,
    y = 0.98 - seq(0, length(values) - 1)*0.10,
    hjust = 1,
    vjust = 1
) {
  # values to string
  val_chr <- numeric2character(values, digits = digits)
  # generate
  if (parse == TRUE) {
    labels <- paste(prefix, "'", val_chr, "'", suffix, sep = sep)
  } else {
    labels <- paste(prefix, val_chr, suffix, sep = sep)
  }
  # return dataframe
  data.frame(
    label = labels,
    x = x,
    y = y,
    hjust = hjust,
    vjust = vjust,
    parse = parse
  )
}


#' List with general plot settings
#'
#' @description
#' Generate a list with 'standard' plot settings.
#'
#' @param color_meas color for measurements
#' @param fill_meas fill color for measurements polygons
#' @param linewidth_meas linewidth for measurements lines
#' @param linetype_meas linetype for measurements lines
#' @param size_meas point size for measurements
#' @param shape_meas shape for measurement points
#' @param alpha_meas transparancy values for measurements polygons
#' @param color_fit color for fitting lines
#' @param fill_fit fill color for fitting polygons
#' @param linewidth_fit linewidth for fitting lines
#' @param linetype_fit linetype for fitting lines
#' @param shape_fit shape for fitting points
#' @param alpha_fit transparancy values for fitting polygons
#' @param color_ann color for text annotations
#' @param size_ann font size for text annotations
#' @return list with fields (input parameter names) and values
#' @examples
#' plot_settings()
#'
plot_settings <- function(
    color_meas = "#377EB8",
    fill_meas = color_meas,
    size_meas = 1,
    linewidth_meas = 0.5,
    linetype_meas = 1,
    shape_meas = 20,
    alpha_meas = 0.25,
    color_fit = "#E41A1C",
    fill_fit = color_fit,
    linewidth_fit = 0.5,
    linetype_fit = 1,
    shape_fit = 20,
    alpha_fit = 0.25,
    color_ann = color_fit,
    size_ann = 3
) {
  list(
    color_meas = color_meas,
    fill_meas = fill_meas,
    size_meas = size_meas,
    linewidth_meas = linewidth_meas,
    linetype_meas = linetype_meas,
    shape_meas = shape_meas,
    alpha_meas = alpha_meas,
    color_fit = color_fit,
    fill_fit = fill_fit,
    linewidth_fit = linewidth_fit,
    linetype_fit = linetype_fit,
    shape_fit = shape_fit,
    alpha_fit = alpha_fit,
    color_ann = color_ann,
    size_ann = size_ann
  )
}

