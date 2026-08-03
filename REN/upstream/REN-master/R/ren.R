#' Main Function for Portfolio Analysis
#'
#' This function integrates data preparation, parallel setup, and portfolio analysis.
#' It takes raw data as input, prepares it using `prepare_data`, sets up parallel processing
#' using `setup_parallel`, and performs the analysis using `perform_analysis`.
#'
#' @param dat A data frame or matrix where the first column is the date and the remaining columns are the data.
#' @param date_column_index The index of the date column in the input data. Default is 1.
#' @param start_date A character string specifying the start date for filtering the data in 'YYYYMMDD' format.
#' Default is '19990101'.
#' @param end_date A character string specifying the end date for filtering the data in 'YYYYMMDD' format.
#' Default is '20231231'.
#' @param num_cores The number of cores to use for parallel processing. Default is 2.
#'
#' @return The results from `perform_analysis`, including plots and performance metrics.
#'
#' @examples
#' \dontrun{
#' # load the sample data
#' dat(FF25)
#' # Run the main function
#' result <- ren(FF25)
#'
#' # Display results
#' print(result$cumulative_return_plot)
#' print(result$cumulative_turnover_plot)
#' }
#'
#' @export
ren <- function(dat, date_column_index = 1, start_date = '19990101', end_date = '20231231', num_cores = 2) {
  # Prepare data
  prepared_data <- prepare_data(dat, date_column_index, start_date, end_date)

  # Setup parallel processing
  cl <- setup_parallel(num_cores)

  # Perform analysis
  result <- perform_analysis(prepared_data$x, prepared_data$mon, prepared_data$count, prepared_data$Date, num_cores)

  # Clean up parallel cluster
  parallel::stopCluster(cl)

  # Return the result
  return(result)
}
