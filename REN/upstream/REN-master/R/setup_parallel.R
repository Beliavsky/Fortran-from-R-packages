#' Setup Parallel Processing for Portfolio Analysis
#'
#' This function sets up parallel processing by loading necessary libraries, allowing the user to specify the
#' number of cores to use, and creating a parallel backend for faster computation.
#'
#' @param num_cores The default number of cores to use for parallel processing. Default is 7.
#'
#' @return A parallel cluster object that can be used with functions that support parallel computation.
#'
#' @details This function allows the user to specify the number of cores for parallel processing either through the
#' argument \code{num_cores} or via interactive user input. The function also loads a set of libraries required for
#' portfolio analysis.
#'
#' @examples
#' # Set up parallel processing with a specified number of cores
#' cl <- setup_parallel(num_cores = 2)  # Use 2 cores for the example
#' print(cl)  # Print the cluster information
#' parallel::stopCluster(cl)  # Stop the cluster after use to clean up
#' @export
setup_parallel <- function(num_cores = 7) {
  # Load necessary libraries
  Packages <- c("glmnet", "quadprog", "doParallel", "lubridate", "Matrix", "tictoc", "corpcor", "ggplot2", "reshape2", "foreach", "stats", "parallel")
  lapply(Packages, library, character.only = TRUE)

  # Allow the user to override the number of cores
  user_input <- as.numeric(readline(prompt = sprintf("Enter the number of cores to use (default is %d): ", num_cores)))
  if (!is.na(user_input) && user_input > 0) {
    num_cores <- user_input
  }

  # Set up parallel processing
  cl <- parallel::makeCluster(num_cores)  # Create a cluster
  doParallel::registerDoParallel(cl)

  return(cl)
}
