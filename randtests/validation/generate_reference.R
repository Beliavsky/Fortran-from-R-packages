# SPDX-License-Identifier: GPL-2.0-or-later
# Run from the randtests package root: Rscript validation/generate_reference.R
# Uses only the pinned upstream snapshot. No R package installation is required.
for (path in sort(list.files('upstream/R', full.names = TRUE))) source(path)
inputs <- list(c(4, 1, 8, 2, 6, 3, 7, 5), c(1, 1, 3, 2, 2, 5, 4, 4), 8:1)
alternatives <- c('two.sided', 'left.sided', 'right.sided')
functions <- list(difference.sign.test, rank.test, turning.point.test, cox.stuart.test,
                  function(x, alternative) runs.test(x, alternative, pvalue = 'normal'),
                  function(x, alternative) runs.test(x, alternative, pvalue = 'exact'),
                  function(x, alternative) bartels.rank.test(x, alternative, pvalue = 'normal'),
                  function(x, alternative) bartels.rank.test(x, alternative, pvalue = 'beta'))
reference <- array(0, c(3, 8, 3, 3))
rows <- list()
for (d in seq_along(inputs)) for (a in seq_along(alternatives)) for (f in seq_along(functions)) {
  result <- functions[[f]](inputs[[d]], alternatives[a])
  sample_size <- if (f %in% 5:6) result$parameter[['n']] else as.numeric(result$parameter[1])
  reference[, f, a, d] <- c(as.numeric(result$statistic), min(1, result$p.value), sample_size)
  rows[[length(rows) + 1]] <- data.frame(dataset = d, alternative = alternatives[a],
    method = f, statistic = as.numeric(result$statistic), p_value = result$p.value, n = sample_size)
}
write.csv(do.call(rbind, rows), 'validation/r_reference.csv', row.names = FALSE)
lines <- c('! SPDX-License-Identifier: GPL-2.0-or-later',
 '! Generated from pinned randtests 1.0.2 by validation/generate_reference.R.',
 'module r_reference', '   use randtests, only: dp', '   implicit none',
 '   real(dp), parameter :: expected(3, 8, 3, 3) = reshape([ &')
numbers <- sprintf('%.17e_dp', as.numeric(reference))
for (i in seq(1, length(numbers), by = 2)) {
  group <- numbers[i:min(i + 1, length(numbers))]
  suffix <- if (i + 1 >= length(numbers)) ' & ' else ', &'
  lines <- c(lines, paste0('      ', paste(group, collapse = ', '), suffix))
}
lines <- c(lines, '      ], [3, 8, 3, 3])', 'end module r_reference')
writeLines(lines, 'test/r_reference.f90')
cat('Generated', length(rows), 'R parity cases with three checked fields each.\n')
