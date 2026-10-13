# SPDX-License-Identifier: GPL-2.0-or-later
# Run from the package root: Rscript --vanilla test/reference/generate.R
# R is needed only to regenerate these frozen reference values, never by FPM.
for (file in list.files('upstream/R', pattern = '\\.R$', full.names = TRUE)) source(file)
quantiles <- qnorm(ppoints(20, a = 3/8))
samples <- list(
  c(-1.6, -1.1, -.9, -.8, -.6, -.5, -.4, -.3, -.2, -.1, 0, .1, .2, .3, .4, .5, .7, .9, 1.2, 1.8),
  sin((1:8) * 1.7) + cos((1:8) * .31),
  quantiles,
  quantiles + .15 * quantiles^2,
  quantiles + .25 * quantiles^2,
  quantiles + .30 * quantiles^2,
  1:100,
  1:101,
  c(rep(0, 199), 1),
  c(rep(-1, 100), rep(1, 100)),
  c(rep(0, 1999), 1),
  sqrt(-2 * log((((1:5000) * 37) %% 997 + .5)/997)) *
    cos(2 * pi * (((1:5000) * 67) %% 991 + .5)/991),
  c(-2, -.7, .1, .8, 1.7)
)
statistic <- matrix(0, 8, length(samples))
p_value <- statistic
n_classes <- matrix(0L, 8, length(samples))
df <- n_classes
statuses <- n_classes
for (j in seq_along(samples)) {
  x <- samples[[j]]
  fs <- list(ad.test, cvm.test, lillie.test, pearson.test,
             function(x) pearson.test(x, adjust = FALSE),
             function(x) pearson.test(x, n.classes = 7),
             function(x) pearson.test(x, n.classes = 7, adjust = FALSE), sf.test)
  for (i in seq_along(fs)) {
    result <- tryCatch(suppressWarnings(fs[[i]](x)), error = function(e) NULL)
    if (is.null(result)) {
      statuses[i, j] <- 1L
    } else {
      statistic[i, j] <- result$statistic
      p_value[i, j] <- result$p.value
      if (!is.null(result$n.classes)) n_classes[i, j] <- result$n.classes
      if (!is.null(result$df)) df[i, j] <- result$df
    }
  }
}
output <- file('test/reference_values.f90', 'w')
writeLines(c('! SPDX-License-Identifier: GPL-2.0-or-later',
             paste('! Generated with', R.version.string, 'and nortest 1.0-4.'),
             '! Regenerate using test/reference/generate.R, from the package root.',
             'module nortest_reference',
             '   use nortest, only: dp',
             '   implicit none',
             '   integer, parameter :: reference_count = 13'), output)
emit <- function(name, x, is_integer = FALSE, shape = NULL) {
  if (is_integer) {
    values <- as.character(as.vector(x))
    declaration <- 'integer'
  } else {
    values <- paste0(formatC(as.vector(x), digits = 17, format = 'e'), '_dp')
    declaration <- 'real(dp)'
  }
  dims <- if (is.null(shape)) length(x) else paste(shape, collapse = ', ')
  intro <- paste0('   ', declaration, ', parameter :: ', name, '(', dims, ') = ')
  if (!is.null(shape)) intro <- paste0(intro, 'reshape(')
  writeLines(paste0(intro, '[ &'), output)
  for (i in seq(1L, length(values), by = 2L)) {
    last <- min(i + 1L, length(values))
    line <- paste(values[i:last], collapse = ', ')
    if (last < length(values)) line <- paste0(line, ', &') else line <- paste0(line, ' &')
    writeLines(paste0('      ', line), output)
  }
  ending <- if (is.null(shape)) '   ]' else paste0('   ], [', dims, '])')
  writeLines(ending, output)
}
emit('normal_quantiles', quantiles)
emit('reference_n', lengths(samples), TRUE)
emit('reference_statistic', statistic, shape = c(8, length(samples)))
emit('reference_p_value', p_value, shape = c(8, length(samples)))
emit('reference_classes', n_classes, TRUE, c(8, length(samples)))
emit('reference_df', df, TRUE, c(8, length(samples)))
emit('reference_status', statuses, TRUE, c(8, length(samples)))
writeLines('end module nortest_reference', output)
close(output)
cat('Generated', length(samples) * 8L, 'R reference evaluations.\n')
