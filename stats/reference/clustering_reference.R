# Generate deterministic reference values for test/test_clustering.f90.
x <- matrix(c(
  0.0, 0.0,
  0.8, 0.2,
  1.1, 1.0,
  5.0, 4.8,
  5.9, 5.1,
  6.2, 6.0
), ncol = 2, byrow = TRUE)

print_vector <- function(name, value) {
  cat(name, "=", paste(sprintf("%.17g", value), collapse = ","), "\n")
}
print_integer_vector <- function(name, value) {
  cat(name, "=", paste(value, collapse = ","), "\n")
}
print_matrix <- function(name, value, integer = FALSE) {
  cat(name, "\n")
  for (i in seq_len(nrow(value))) {
    formatted <- if (integer) value[i, ] else sprintf("%.17g", value[i, ])
    cat(paste(formatted, collapse = ","), "\n")
  }
}

print_vector("dist_euclidean", as.vector(dist(x, method = "euclidean")))
print_vector("dist_manhattan", as.vector(dist(x, method = "manhattan")))

initial <- matrix(c(0.0, 0.0, 6.2, 6.0), ncol = 2, byrow = TRUE)
fit <- kmeans(x, centers = initial, algorithm = "Lloyd", iter.max = 20)
print_matrix("kmeans_centers", fit$centers)
print_integer_vector("kmeans_cluster", fit$cluster)
print_integer_vector("kmeans_size", fit$size)
print_vector("kmeans_withinss", fit$withinss)
print_vector("kmeans_sums", c(fit$totss, fit$tot.withinss, fit$betweenss, fit$iter))

tree <- hclust(dist(x), method = "complete")
print_matrix("hclust_merge", tree$merge, integer = TRUE)
print_vector("hclust_height", tree$height)
print_integer_vector("hclust_order", tree$order)
print_integer_vector("cutree_k2", cutree(tree, k = 2))
