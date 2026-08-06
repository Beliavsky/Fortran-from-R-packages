test_that("Valid Input", {
  data(indipsa)
  ipsa<-diff(log(indipsa))
  ipsa<-100*ipsa
  expect_error(tvGarchKalmanFit(ipsa, c = c(0.05,0.05), alpha = c(0.05,0.05), beta = c(0.05,0.05),
                                type = c("trigonometrico","trigonometric","trigonometric"), trig ="cos", arg = "3*(1-log(u))"),"invalid type for c", fixed=T)
  expect_error(tvGarchKalmanFit(ipsa, c = c(0.05,0.05), alpha = c(0.05,0.05), beta = c(0.05,0.05),
                                type = c("trigonometric","trigonometrico","trigonometric"), trig ="cos", arg = "3*(1-log(u))"),"invalid type for a", fixed=T)
  expect_error(tvGarchKalmanFit(ipsa, c = c(0.05,0.05), alpha = c(0.05,0.05), beta = c(0.05,0.05),
                                type = c("trigonometric","trigonometric","trigonometrico"), trig ="cos", arg = "3*(1-log(u))"),"invalid type for b", fixed=T)
})

test_that("Data is loaded correctly", {
  expect_true(exists("indipsa"))
  expect_equal(ncol(indipsa), 1)  # verificar el número de columnas
  expect_equal(nrow(indipsa), 3186) # Verificar el número de filas
})

test_that("Example with real data",{
  data(indipsa)
  ipsa<-diff(log(indipsa))
  ipsa<-100*ipsa
  fit <- tvGarchKalmanFit(ipsa, c = c(0.05,0.05), alpha = c(0.05,0.05), beta = c(0.05,0.05),
                          type = c("trigonometric","trigonometric","trigonometric"), trig ="cos", arg = "3*(1-log(u))")
  result = c(0.03001, 0.01237, 0.14990, 0.02919, 0.82541, -0.04145)
  expect_equal(fit,result)
  model<-tvGarchKalmanPrint(fit, ipsa, c = c(0.05,0.05), alpha = c(0.05,0.05), beta = c(0.05,0.05), type = c("trigonometric","trigonometric","trigonometric"),
                            trig ="cos", arg = "3*(1-log(u))", predict = 10)
  expect_equal(model$loglik, 2375.377)
})
