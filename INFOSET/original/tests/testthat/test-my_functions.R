test_that("infoset function works as expected", {
  skip_on_cran()  # Skip this test on CRAN

  # Your test code here
  gross.ret<-as.data.frame(lapply(sample.data, g_ret))
  infoset(gross.ret$ETF_1, plot_cp = "T")
  LR <- LR_cp(sample.data, FT= 1290, ov = 125)
  df <- as.data.frame(matrix(unlist(LR), nrow = length(LR), ncol = ncol(sample.data)))
  colnames(df) <- c(paste("tw", rep(1:16)))
  plot(df[,1], pch=19, col=asset.label$label, ylab="LR_cp", xlab="ETFs")

  })
