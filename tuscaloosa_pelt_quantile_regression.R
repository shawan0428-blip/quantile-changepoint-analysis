tusc <- data.frame(
  year = as.integer(year_text),
  tmp = suppressWarnings(as.numeric(raw_data[[2L]]))
)
if (any(!is.finite(tusc$tmp))) {
  stop("Temperatures must be numeric and non-missing. Check the input data.")
}
tusc <- tusc[order(tusc$year), ]
rownames(tusc) <- NULL
n <- nrow(tusc)
if (n < min_seg_len) stop("There are fewer observations than minseglen.")

# A lagged observation should represent the previous year.
if (any(diff(tusc$year) != 1L)) {
  stop("Expected consecutive annual observations; check for duplicate or missing years.")
}

changes <- PELT.quantloss(
  tusc$tmp,
  taus = taus,
  pen = 18 * log(n),
  minseglen = min_seg_len
)

# Treat a change point as the last observation of the preceding segment.
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

if (length(changes)) {
  message("Detected change points (last observation of each preceding segment):")
  print(data.frame(index = changes, year = tusc$year[changes]))
} else {
  message("No change points detected; fitting the full series as one segment.")
}

segments <- lapply(seq_len(length(boundaries) - 1L), function(k) {
  seq.int(boundaries[k] + 1L, boundaries[k + 1L])
})

fit_segment <- function(index, data, taus) {
  segment <- data[index, ]
  trend_fit <- lm(tmp ~ year, data = segment)
  trend <- fitted(trend_fit)
  
  # Retain the original two-stage model: detrended response on raw lagged
  # temperature, with lagged pairs kept within each segment.
  qr_data <- data.frame(
    detrended = segment$tmp[-1L] - trend[-1L],
    lag_tmp = head(segment$tmp, -1L)
  )
  quantile_fits <- vapply(taus, function(tau) {
    fit <- quantreg::rq(detrended ~ lag_tmp, tau = tau, data = qr_data)
    trend[-1L] + as.numeric(fitted(fit))
  }, numeric(nrow(segment) - 1L))
  colnames(quantile_fits) <- paste0("tau_", taus)
  
  list(year = segment$year, tmp = segment$tmp,
       trend = trend, quantiles = quantile_fits)
}

segment_fits <- lapply(segments, fit_segment, data = tusc, taus = plot_taus)

plot_analysis <- function(data, changes, fits, taus, reference_years) {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  n_panels <- length(fits) + 1L
  par(mfrow = c(ceiling(n_panels / 2L), 2L), mar = c(4, 4, 3, 1))
  
  y_limits <- range(data$tmp, unlist(lapply(fits, function(fit) {
    c(fit$trend, fit$quantiles)
  })))
  quantile_colors <- colorRampPalette(c("blue", "red"))(length(taus))
  reference_years <- reference_years[
    reference_years >= min(data$year) & reference_years <= max(data$year)
  ]
  
  plot(data$year, data$tmp, type = "o", pch = 16, cex = 0.3,
       col = "grey55", xlab = "Year", ylab = "Temperature",
       main = "Annual temperature: Tuscaloosa, AL", ylim = y_limits)
  for (fit in fits) {
    lines(fit$year, fit$trend, col = "darkgreen", lwd = 1.5)
  }
  if (length(changes)) {
    abline(v = data$year[changes], col = "black", lty = 3)
  }
  if (length(reference_years)) {
    abline(v = reference_years, col = "blue", lty = 2)
  }
  legend("topleft",
         legend = c("Observed temperature", "OLS trend",
                    if (length(changes)) "Detected changes" else "Detected changes: none",
                    if (length(reference_years)) "Reference years" else "Reference years: outside range"),
         col = c("grey55", "darkgreen", "black", "blue"),
         lty = c(1, 1, 3, 2), bty = "n", cex = 0.7)
  
  for (k in seq_along(fits)) {
    fit <- fits[[k]]
    title <- sprintf("Segment %d (%d-%d)", k, min(fit$year), max(fit$year))
    plot(fit$year, fit$tmp, type = "o", pch = 16, cex = 0.3,
         col = "grey55", xlab = "Year", ylab = "Temperature",
         main = title, ylim = y_limits)
    lines(fit$year, fit$trend, col = "darkgreen", lwd = 1.5)
    
    for (i in seq_along(taus)) {
      lines(fit$year[-1L], fit$quantiles[, i],
            col = quantile_colors[i], lwd = 1.5)
    }
    legend("topleft", legend = c("OLS trend", paste("Quantile", taus)),
           col = c("darkgreen", quantile_colors), lty = 1,
           lwd = 1.5, bty = "n", cex = 0.7)
  }
}

plot_analysis(tusc, changes, segment_fits, plot_taus, reference_years)
