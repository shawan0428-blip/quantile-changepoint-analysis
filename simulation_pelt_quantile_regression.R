# PELT change-point detection and segment-wise quantile regression.
# Put the definition of PELT.quantloss() in PELT_quantloss.R, or load it first.

if (!requireNamespace("quantreg", quietly = TRUE)) {
  stop('Install quantreg first: install.packages("quantreg")')
}

if (!exists("PELT.quantloss", mode = "function") &&
    file.exists("PELT_quantloss.R")) {
  source("PELT_quantloss.R")
}
if (!exists("PELT.quantloss", mode = "function")) {
  stop("Load PELT.quantloss() first, or run from the folder containing PELT_quantloss.R.")
}

set.seed(123)

true_change <- 50L
taus <- c(0.25, 0.50, 0.75)
min_seg_len <- 20L

# Both the AR coefficient and the deterministic trend change after time 50.
series <- c(
  arima.sim(n = 50, model = list(ar = 0.7)) + 0.01 * (1:50),
  arima.sim(n = 50, model = list(ar = -0.7)) + 0.03 * (51:100)
)
n <- length(series)

changes <- PELT.quantloss(
  series,
  taus = taus,
  pen = 12 * log(n),
  minseglen = min_seg_len
)

# Assume each change point is the LAST observation of its segment.
# Ignore 0 and n if the detector includes them as boundary markers.
if (is.null(changes)) changes <- numeric(0)
if (!is.numeric(changes) || !is.null(dim(changes)) ||
    any(!is.finite(changes)) || any(changes != floor(changes)) ||
    any(changes < 0 | changes > n)) {
  stop("Expected PELT.quantloss() to return a numeric vector of segment-end indices.")
}
changes <- sort(unique(as.integer(changes[changes > 0 & changes < n])))
boundaries <- c(0L, changes, n)

if (any(diff(boundaries) < min_seg_len)) {
  stop("A segment is shorter than minseglen; check the detector's output convention.")
}

message("Detected change points: ",
        if (length(changes)) paste(changes, collapse = ", ") else "none")

segments <- lapply(seq_len(length(boundaries) - 1L), function(k) {
  seq.int(boundaries[k] + 1L, boundaries[k + 1L])
})

fit_segment <- function(index, series, taus) {
  time <- index
  values <- series[index]
  trend_fit <- lm(values ~ time)
  trend <- fitted(trend_fit)
  
  # Keep the original two-stage fit: detrend the response, then regress
  # on the raw lagged series. No lagged pair crosses a segment boundary.
  qr_data <- data.frame(
    detrended = values[-1L] - trend[-1L],
    lag_y = head(values, -1L)
  )
  
  quantile_fits <- vapply(taus, function(tau) {
    fit <- quantreg::rq(detrended ~ lag_y, tau = tau, data = qr_data)
    trend[-1L] + as.numeric(fitted(fit))
  }, numeric(length(values) - 1L))
  colnames(quantile_fits) <- paste0("tau_", taus)
  
  list(time = time, trend = trend, quantiles = quantile_fits)
}

segment_fits <- lapply(segments, fit_segment, series = series, taus = taus)

plot_analysis <- function(series, changes, fits, taus, true_change) {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mfrow = c(ceiling((length(taus) + 1L) / 2L), 2L),
      mar = c(4, 4, 3, 1))
  
  # Use the same limits for every panel, including all fitted values.
  y_limits <- range(series, unlist(lapply(fits, function(fit) {
    c(fit$trend, fit$quantiles)
  })))
  change_labels <- c(
    paste0("True change (", true_change, ")"),
    if (length(changes)) "Detected change" else "Detected change: none"
  )
  