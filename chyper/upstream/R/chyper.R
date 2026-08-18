#' @importFrom stats runif dhyper
#' @title Probability mass function for conditional hypergeometric distributions
#'
#' @description Calculates the PMF of a conditional hypergeometric distribution: the distribution of how many items are in the overlap of all samples when samples of arbitrary size are each taken without replacement from populations of arbitrary size.
#'
#' @param k an integer or vector of integers representing the overlap size
#' @param s an integer representing the size of the intersecting population
#' @param n a vector of integers representing the sizes of each non-intersecting population
#' @param m a vector of integers representing the sample sizes
#' @param verbose T/F should intermediate messages be printed?
#'
#' @return The probability of sampling \code{k} of the same items in all samples
#'
#' @examples dchyper(c(3,5), 10, c(12,13,14), c(7,8,9))
#'
#' @export dchyper

# The following implements the conditional hypergeometric PMF as a function
# parameterized by the number or numbers of people sampled in all samples
# (the value of interest), the size of the overlapping population
# (an integer), the size of the unique populations (a vector of integers at
# least 2 long for 2 populations), and the number of people sampled from
# each population (a vector of integers at least 2 long for 2 samples).

dchyper <- function(k, s, n, m, verbose=T) {
  # Check that the inputs are valid

  # Ensure at least two unique populations are given
  if (length(n) < 2) {stop("Need two regions (n too short)")}

  # Check that only one overlapping population is given
  if (length(s) > 1) {stop("Only 1 intersection region (s too long)")}

  # Check that the number of unique populations and the number of samples are the same
  if (length(m) != length(n)) {stop("n and m must match in length")}

  # Check that all population and sample sizes are non-negative
  if (s < 0 | any(n < 0) | any(m < 0)) {stop("Only positive values allowed")}

  # Check that sample sizes are not bigger than the populations they are from
  if (any(m > (n + s))) {stop("Sample size cannot be bigger than population")}

  # Check that all inputs are integers
  if (!all(c(k, s, n, m) == floor(c(k, s, n, m)))) {stop("Inputs must be integers")}

  # Get the max value of interest in the possible overlap range
  top <- max(k[k<=min(s,m) & k >= 0], 0)

  # Store the rest of the inputs
  store <- k

  # Let k be a vector of all integers up to the maximum value of interest
  k <- 0:top

  # Reorder samples and populations for optimal speed
  # Evaluate the smallest population first to get the most 0s in the rest of the data
  df <- data.frame(n, m)
  df <- df[order(df$m, decreasing = T),]
  n <- df$n
  m <- df$m

  # Calculate the probabilities of each number of overlaps
  # for one level (i.e. P(Xi=xi) for 0 to mi) from 0 to the largest
  # possible overlap size
  calculateConsecutive <- function(s, n, m) {
    # If there is only one population, apply the base case of
    # a HGeom(s, n1, m1)
    if (length(n) == 1) {
      PO1 <- dhyper(0:min(m, s),
                    s,
                    n,
                    m)
      if (verbose) {message("Finished 1 level...")}
      return(PO1)
    }
    else {
      # Otherwise get all the values from 0 to the max overlap size
      x1 <- 0:min(m[2], s)

      # For use in determining whether there is rounding-to-zero error
      wasPos <- rep(F,length(x1))

      # Get the previous level
      PO1 <- calculateConsecutive(s, n[2:length(n)], m[2:length(m)])

      # Create the current level
      PO2 <- vector(length=length(0:min(m[1], s)))

      # Find P(Xi=0)
      save <- dhyper(0,
                     x1,
                     n[1] + s - x1,
                     m[1])
      PO2[1] <- sum(save  * PO1)

      # Count how many numbers were positive on the prior iteration
      numWasPosPrior <- sum(wasPos)

      # Update to show if the value at an index has ever been positive
      wasPos <- wasPos | (save!=0)

      # Check how many of the current terms are now positive
      numWasPos <- sum(wasPos)

      # Update the index
      i <- 1

      # While the number of positive terms is increasing each round, continue to
      # recalculate with dhyper to avoid the round-to-zero issue explained in the
      # vignette
      while ((numWasPosPrior < numWasPos | (PO2[i+1] == 0)) & i <= min(m[1], s)) {
        # Conditional probability of current level
        save <- dhyper(i, x1[-i:-1], n[1] + s - x1[-i:-1], m[1])

        # Update how many were positive
        numWasPosPrior <- numWasPos
        wasPos <- wasPos[-1] | (save!=0)
        numWasPos <- sum(wasPos)

        # Store the new value
        PO2[i + 1] <- sum(save * PO1[-i:-1])

        # Increment index
        i <- i + 1
      }
      # Quick variable swap
      j <- i

      # Ensure that the whole level hasn't been found yet
      if (j <= min(m[1], s)) {
        # Update more efficiently with the indexing shifting and coefficient
        # multiplication as described in the vignette
        for(i in j:min(m[1], s)) {
          save <- save[-1] * (x1[-i:-1] - (i-1)) / ((i - 1) + 1) *
            (m[1] - (i - 1)) / (n[1] + s - x1[-i:-1] - (m[1] - ((i - 1) + 1)))
          PO2[i + 1] <- sum(save * PO1[-i:-1])
        }
      }
    }

    # Give status update
    if (verbose) {message("Finished another level...")}

    # Return one level (i.e. P(Xi=xi) for 0 to mi)
    return(PO2)
  }

  # Do the same as calculateConsecutive except stop at the
  # maximum desired value
  outerConsecutive <- function(k, s, n, m) {

    # Get all the values from 0 to the max overlap size
    x1 <- 0:min(m[2], s)

    # For use in determining whether there is rounding-to-zero error
    wasPos <- rep(F,length(x1))

    # Get the previous level
    PO1 <- calculateConsecutive(s, n[2:length(n)], m[2:length(m)])

    # Create the current level
    PO2 <- vector(length=length(0:min(m[1], s, k)))

    # Find P(Xi=0)
    save <- dhyper(0,
                   x1,
                   n[1] + s - x1,
                   m[1])
    PO2[1] <- sum(save  * PO1)

    # Count how many numbers were positive on the prior iteration
    numWasPosPrior <- sum(wasPos)

    # Update to show if the value at an index has ever been positive
    wasPos <- wasPos | (save!=0)

    # Check how many of the current terms are now positive
    numWasPos <- sum(wasPos)

    # Update the index
    i <- 1

    # While the number of positive terms is increasing each round, continue to
    # recalculate with dhyper to avoid the round-to-zero issue explained in the
    # vignette
    while ((numWasPosPrior < numWasPos | (PO2[i+1] == 0)) & i <= min(m[1], s, k)) {
      # Conditional probability of current level
      save <- dhyper(i, x1[-i:-1], n[1] + s - x1[-i:-1], m[1])

      # Update how many were positive
      numWasPosPrior <- numWasPos
      wasPos <- wasPos[-1] | (save!=0)
      numWasPos <- sum(wasPos)

      # Store the new value
      PO2[i + 1] <- sum(save * PO1[-i:-1])

      # Increment index
      i <- i + 1
    }
    # Quick variable swap
    j <- i

    # Ensure that the whole level hasn't been found yet
    if (j <= min(m[1], s, k)) {

      # Update more efficiently with the indexing shifting and coefficient
      # multiplication as described in the vignette up to the maximum value
      # desired
      for(i in j:min(m[1], s, k)) {
        save <- save[-1] * (x1[-i:-1] - (i-1)) / ((i - 1) + 1) *
          (m[1] - (i - 1)) / (n[1] + s - x1[-i:-1] - m[1] + (i - 1) + 1)
        PO2[i + 1] <- sum(save * PO1[-i:-1])
      }
    }

    # Give status update
    if (verbose) {message("Finished another level...")}

    # Return one level (i.e. P(Xi=xi) for 0 to mi)
    return(PO2)
  }

  # Use calculateConsecutive if only one top level probability is requested
  useConsecutive <- function(k, s, n, m) {
    # Get the previous level
    PO1 <- calculateConsecutive(s, n[-1], m[-1])

    # Calculate the single top-level value
    x1 <- 0:min(m[2], s)
    save <- dhyper(k,
                   x1,
                   n[1] + s - x1,
                   m[1])
    PO2 <- sum(save * PO1)

    # Return one value from the top level
    return(PO2)
  }

  # If the length of the desired values is only one, use the faster single
  # value subfunction
  if (length(store) == 1) {
    return(useConsecutive(store, s, n, m))
  }

  # Otherwise use the calculate-all-with-updating
  else {
    allValues <- outerConsecutive(top, s, n, m)

    # Return 0 if the input is outside the possible range; otherwise return the calculated value
    return(ifelse(store > min(s,m), 0, ifelse(store < 0, 0, allValues[pmax(store,0) + 1])))
  }
}

#' @title Cumulative density function for conditional hypergeometric distributions
#'
#' @description Calculates the CDF of a conditional hypergeometric distribution: the distribution of how many items are in the overlap of all samples when samples of arbitrary size are each taken without replacement from populations of arbitrary size.
#'
#' @param k an integer or vector of integers representing the overlap size
#' @param s an integer representing the size of the intersecting population
#' @param n a vector of integers representing the sizes of each non-intersecting population
#' @param m a vector of integers representing the sample sizes
#' @param verbose T/F should intermediate messages be printed?
#'
#' @return The probability of sampling \code{k} or less of the same items in all samples
#'
#' @examples pchyper(c(3,5), 10, c(12,13,14), c(7,8,9))
#'
#' @export pchyper

# The following implements the conditional hypergeometric CDF as a function
# parameterized by the number or numbers of people sampled in all samples
# (the value of interest), the size of the overlapping population
# (an integer), the size of the unique populations (a vector of integers at
# least 2 long for 2 populations), and the number of people sampled from
# each population (a vector of integers at least 2 long for 2 samples).

pchyper <- function(k, s, n, m, verbose=T) {
  # Check that the inputs are valid

  # Ensure at least two unique populations are given
  if (length(n) < 2) {stop("Need two regions (n too short)")}

  # Check that only one overlapping population is given
  if (length(s) > 1) {stop("Only 1 intersection region (s too long)")}

  # Check that the number of unique populations and the number of samples are the same
  if (length(m) != length(n)) {stop("n and m must match in length")}

  # Check that all population and sample sizes are non-negative
  if (s < 0 | any(n < 0) | any(m < 0)) {stop("Only positive values allowed")}

  # Check that sample sizes are not bigger than the populations they are from
  if (any(m > (n + s))) {stop("Sample size cannot be bigger than population")}

  # Check that all inputs are integers
  if (!all(c(k, s, n, m) == floor(c(k, s, n, m)))) {stop("Inputs must be integers")}

  # Reorder samples and populations for optimal speed
  # Evaluate the smallest population first to get the most 0s in the rest of the data
  df <- data.frame(n, m)
  df <- df[order(df$m, decreasing = T),]
  n <- df$n
  m <- df$m

  # Get the max value of interest in the possible overlap range
  top <- max(k[k<=min(s,m) & k >= 0], 0)

  # Store the rest of the inputs
  store <- k

  # Let k be a vector of all integers up to the maximum value of interest
  k <- 0:top


  # Calculate the probabilities of each number of overlaps
  # for one level (i.e. P(Xi=xi) for 0 to mi) from 0 to the largest
  # possible overlap size
  calculateConsecutive <- function(s, n, m) {
    # If there is only one population, apply the base case of
    # a HGeom(s, n1, m1)
    if (length(n) == 1) {
      PO1 <- dhyper(0:min(m, s),
                    s,
                    n,
                    m)
      if (verbose) {message("Finished 1 level...")}
      return(PO1)
    }
    else {
      # Otherwise get all the values from 0 to the max overlap size
      x1 <- 0:min(m[2], s)

      # For use in determining whether there is rounding-to-zero error
      wasPos <- rep(F,length(x1))

      # Get the previous level
      PO1 <- calculateConsecutive(s, n[2:length(n)], m[2:length(m)])

      # Create the current level
      PO2 <- vector(length=length(0:min(m[1], s)))

      # Find P(Xi=0)
      save <- dhyper(0,
                     x1,
                     n[1] + s - x1,
                     m[1])
      PO2[1] <- sum(save  * PO1)

      # Count how many numbers were positive on the prior iteration
      numWasPosPrior <- sum(wasPos)

      # Update to show if the value at an index has ever been positive
      wasPos <- wasPos | (save!=0)

      # Check how many of the current terms are now positive
      numWasPos <- sum(wasPos)

      # Update the index
      i <- 1

      # While the number of positive terms is increasing each round, continue to
      # recalculate with dhyper to avoid the round-to-zero issue explained in the
      # vignette
      while ((numWasPosPrior < numWasPos | (PO2[i+1] == 0)) & i <= min(m[1], s)) {
        # Conditional probability of current level
        save <- dhyper(i, x1[-i:-1], n[1] + s - x1[-i:-1], m[1])

        # Update how many were positive
        numWasPosPrior <- numWasPos
        wasPos <- wasPos[-1] | (save!=0)
        numWasPos <- sum(wasPos)

        # Store the new value
        PO2[i + 1] <- sum(save * PO1[-i:-1])

        # Increment index
        i <- i + 1
      }
      # Quick variable swap
      j <- i

      # Ensure that the whole level hasn't been found yet
      if (j <= min(m[1], s)) {
        # Update more efficiently with the indexing shifting and coefficient
        # multiplication as described in the vignette
        for(i in j:min(m[1], s)) {
          save <- save[-1] * (x1[-i:-1] - (i-1)) / ((i - 1) + 1) *
            (m[1] - (i - 1)) / (n[1] + s - x1[-i:-1] - (m[1] - ((i - 1) + 1)))
          PO2[i + 1] <- sum(save * PO1[-i:-1])
        }
      }
    }

    # Give status update
    if (verbose) {message("Finished another level...")}

    # Return one level (i.e. P(Xi=xi) for 0 to mi)
    return(PO2)
  }

  # Calculate the cumulative probability of each number of overlaps
  # for the top level up to the maximum desired number
  outerConsecutive <- function(k, s, n, m) {

    # Get all the values from 0 to the max overlap size
    x1 <- 0:min(m[2], s)

    # For use in determining whether there is rounding-to-zero error
    wasPos <- rep(F,length(x1))

    # Get the previous level
    PO1 <- calculateConsecutive(s, n[2:length(n)], m[2:length(m)])

    # Create the current level
    PO2 <- vector(length=length(0:min(m[1], s, k)))

    # Find P(Xi=0)
    save <- dhyper(0,
                   x1,
                   n[1] + s - x1,
                   m[1])
    PO2[1] <- sum(save  * PO1)

    # Count how many numbers were positive on the prior iteration
    numWasPosPrior <- sum(wasPos)

    # Update to show if the value at an index has ever been positive
    wasPos <- wasPos | (save!=0)

    # Check how many of the current terms are now positive
    numWasPos <- sum(wasPos)

    # Update the index
    i <- 1

    # While the number of positive terms is increasing each round, continue to
    # recalculate with dhyper to avoid the round-to-zero issue explained in the
    # vignette
    while ((numWasPosPrior < numWasPos | (PO2[i+1] == 0)) & i <= min(m[1], s, k)) {
      # Conditional probability of current level
      save <- dhyper(i, x1[-i:-1], n[1] + s - x1[-i:-1], m[1])

      # Update how many were positive
      numWasPosPrior <- numWasPos
      wasPos <- wasPos[-1] | (save!=0)
      numWasPos <- sum(wasPos)

      # Store the new value and add to the previous for the cumulative value
      PO2[i + 1] <- PO2[i] + sum(save * PO1[-i:-1])

      # If the cumulative probability reaches 1, set the rest to 1 also
      if (PO2[i + 1] == 1) {
        PO2[(i + 1):length(PO2)] <- 1
        i <- min(m[1], s, k) + 1
        break
      }

      # Update the index
      i <- i + 1
    }

    # Quick variable swap
    j <- i

    # Ensure that the whole level hasn't been found yet
    if (j <= min(m[1], s, k)) {

      # Update more efficiently with the indexing shifting and coefficient
      # multiplication as described in the vignette up to the maximum value
      # desired
      for(i in j:min(m[1], s, k)) {
        save <- save[-1] * ((x1[-i:-1] - (i-1)) / ((i - 1) + 1) *
                              (m[1] - (i - 1)) /
                              (n[1] + s - x1[-i:-1] - m[1] + (i - 1) + 1))
        PO2[i + 1] <- PO2[i] + sum(save * PO1[-i:-1])

        # If the cumulative probability reaches 1, set the rest to 1 also
        if (PO2[i + 1] == 1) {
          PO2[(i + 1):length(PO2)] <- 1
          break
        }
      }
    }

    # Give status update
    if (verbose) {message("Finished another level...")}

    # Return the top level
    return(PO2)
  }

  # Get the cumulative probabilities up to the max value requested
  temp <- outerConsecutive(top, s, n, m)

  # Return 0 if the input is less than the possible range; 1 if it is more than the range,
  # and otherwise return the calculated value
  output <- ifelse(store < 0, 0, ifelse(store > min(s,m), 1, temp[pmax(store,0) + 1]))
  return (output)
}

#' @title Quantile function for conditional hypergeometric distributions
#'
#' @description Calculates the quantile function of a conditional hypergeometric distribution: the distribution of how many items are in the overlap of all samples when samples of arbitrary size are each taken without replacement from populations of arbitrary size.
#'
#' @param p the desired quantile or quantiles
#' @param s an integer representing the size of the intersecting population
#' @param n a vector of integers representing the sizes of each non-intersecting population
#' @param m a vector of integers representing the sample sizes
#' @param verbose T/F should intermediate messages be printed?
#'
#' @return The minimum integer (or integers for a vector input) such that the input probability is less than or equal to the probability of sampling that many of the same items in all samples.
#'
#' @examples qchyper(c(0,0.9,1), 10, c(12,13,14), c(7,8,9))
#'
#' @export qchyper

# The following implements the conditional hypergeometric quantile function
# parameterized by the proportion(s) of interest, the size of the
# overlapping population (an integer), the size of the unique populations (a
# vector of integers at least 2 long for 2 populations), and the number of
# people sampled from each population (a vector of integers at least 2 long
# for 2 samples).

qchyper <- function(p, s, n, m, verbose=T) {
  # Check that the inputs are valid

  # Ensure at least two unique populations are given
  if (length(n) < 2) {stop("Need two regions (n too short)")}

  # Check that only one overlapping population is given
  if (length(s) > 1) {stop("Only 1 intersection region (s too long)")}

  # Check that the number of unique populations and the number of samples are the same
  if (length(m) != length(n)) {stop("n and m must match in length")}

  # Check that all population and sample sizes are non-negative
  if (s < 0 | any(n < 0) | any(m < 0)) {stop("Only positive values allowed")}

  # Check that sample sizes are not bigger than the populations they are from
  if (any(m > (n + s))) {stop("Sample size cannot be bigger than population")}

  # Check that all inputs are integers
  if (!all(c(s, n, m) == floor(c(s, n, m)))) {stop("Inputs must be integers")}

  # Return a single quantile
  getEach <- function(p, cumsum, s, m) {

    # If the probability is more than 1, return the max overlap size
    if(p >= 1) {
      return (min(s, m))

    # If the probability is less than 0, return 0
    } else if (p <= 0) {
      return(0)
    }

    # Otherwise return the first value with a cumulative probability
    # greater than the input
    return(which(cumsum > p)[1] - 1)
  }

  # Store the largest input not greater than or equal to 1
  if(length(p[p<1]) >= 1) {
    biggest <- max(p[p<1])
  }
  # If all are greater than 1, set the max at 1
  else {
    biggest <- 1
  }

  # If any inputs are less than 1, do the actual calculations
  if (any(p < 1)) {
    # Reorder samples and populations for optimal speed
    # Evaluate the smallest population first to get the most 0s in the rest of the data
    df <- data.frame(n, m)
    df <- df[order(df$m, decreasing = T),]
    n <- df$n
    m <- df$m

    # Calculate the probabilities of each number of overlaps
    # for one level (i.e. P(Xi=xi) for 0 to mi) from 0 to the largest
    # possible overlap size
    calculateConsecutive <- function(s, n, m) {
      # If there is only one population, apply the base case of
      # a HGeom(s, n1, m1)
      if (length(n) == 1) {
        PO1 <- dhyper(0:min(m, s),
                      s,
                      n,
                      m)
        if (verbose) {message("Finished 1 level...")}
        return(PO1)
      }
      else {
        # Otherwise get all the values from 0 to the max overlap size
        x1 <- 0:min(m[2], s)

        # For use in determining whether there is rounding-to-zero error
        wasPos <- rep(F,length(x1))

        # Get the previous level
        PO1 <- calculateConsecutive(s, n[2:length(n)], m[2:length(m)])

        # Create the current level
        PO2 <- vector(length=length(0:min(m[1], s)))

        # Find P(Xi=0)
        save <- dhyper(0,
                       x1,
                       n[1] + s - x1,
                       m[1])
        PO2[1] <- sum(save  * PO1)

        # Count how many numbers were positive on the prior iteration
        numWasPosPrior <- sum(wasPos)

        # Update to show if the value at an index has ever been positive
        wasPos <- wasPos | (save!=0)

        # Check how many of the current terms are now positive
        numWasPos <- sum(wasPos)

        # Update the index
        i <- 1

        # While the number of positive terms is increasing each round, continue to
        # recalculate with dhyper to avoid the round-to-zero issue explained in the
        # vignette
        while ((numWasPosPrior < numWasPos | (PO2[i+1] == 0)) & i <= min(m[1], s)) {
          # Conditional probability of current level
          save <- dhyper(i, x1[-i:-1], n[1] + s - x1[-i:-1], m[1])

          # Update how many were positive
          numWasPosPrior <- numWasPos
          wasPos <- wasPos[-1] | (save!=0)
          numWasPos <- sum(wasPos)

          # Store the new value
          PO2[i + 1] <- sum(save * PO1[-i:-1])

          # Increment index
          i <- i + 1
        }
        # Quick variable swap
        j <- i

        # Ensure that the whole level hasn't been found yet
        if (j <= min(m[1], s)) {
          # Update more efficiently with the indexing shifting and coefficient
          # multiplication as described in the vignette
          for(i in j:min(m[1], s)) {
            save <- save[-1] * (x1[-i:-1] - (i-1)) / ((i - 1) + 1) *
              (m[1] - (i - 1)) / (n[1] + s - x1[-i:-1] - (m[1] - ((i - 1) + 1)))
            PO2[i + 1] <- sum(save * PO1[-i:-1])
          }
        }
      }

      # Give status update
      if (verbose) {message("Finished another level...")}

      # Return one level (i.e. P(Xi=xi) for 0 to mi)
      return(PO2)
    }

    # Calculate the cumulative probability of each number of overlaps
    # for the top level up to the maximum desired number
    outerConsecutive <- function(biggest, s, n, m) {

      # Get all the values from 0 to the max overlap size
      x1 <- 0:min(m[2], s)

      # For use in determining whether there is rounding-to-zero error
      wasPos <- rep(F,length(x1))

      # Get the previous level
      PO1 <- calculateConsecutive(s, n[2:length(n)], m[2:length(m)])

      # Create the current level
      PO2 <- vector(length=length(0:min(m[1], s)))

      # Find P(Xi=0)
      save <- dhyper(0,
                     x1,
                     n[1] + s - x1,
                     m[1])
      PO2[1] <- sum(save  * PO1)

      # Count how many numbers were positive on the prior iteration
      numWasPosPrior <- sum(wasPos)

      # Update to show if the value at an index has ever been positive
      wasPos <- wasPos | (save!=0)

      # Check how many of the current terms are now positive
      numWasPos <- sum(wasPos)

      # Update the index
      i <- 1

      # While the number of positive terms is increasing each round, continue to
      # recalculate with dhyper to avoid the round-to-zero issue explained in the
      # vignette
      while ((numWasPosPrior < numWasPos | (PO2[i+1] == 0)) &
             i <= min(m[1], s) & PO2[i] <= biggest) {

        # Conditional probability of current level
        save <- dhyper(i, x1[-i:-1], n[1] + s - x1[-i:-1], m[1])

        # Update how many were positive
        numWasPosPrior <- numWasPos
        wasPos <- wasPos[-1] | (save!=0)
        numWasPos <- sum(wasPos)

        # Store the new value and add to the previous for the cumulative value
        PO2[i + 1] <- PO2[i] + sum(save * PO1[-i:-1])

        # If the cumulative probability reaches 1, set the rest to 1 also
        if (PO2[i + 1] == 1) {
          PO2[(i + 1):length(PO2)] <- 1
          i <- min(m[1], s) + 1
          break
        }

        # Update index
        i <- i + 1
      }

      # Ensure that the whole level hasn't been found yet and we aren't at the end
      while (i <= min(m[1], s) & PO2[i] <= biggest) {

        # Update more efficiently with the indexing shifting and coefficient
        # multiplication as described in the vignette up to the maximum value
        # desired
        save <- save[-1] * ((x1[-i:-1] - (i-1)) / ((i - 1) + 1) *
                              (m[1] - (i - 1)) /
                              (n[1] + s - x1[-i:-1] - m[1] + (i - 1) + 1))
        PO2[i + 1] <- PO2[i] + sum(save * PO1[-i:-1])

        # If the cumulative probability reaches 1, set the rest to 1 also
        if (PO2[i + 1] == 1) {
          PO2[(i + 1):length(PO2)] <- 1
          break
        }
        i <- i+1
      }

      # Give status update
      if (verbose) {message("Finished another level...")}

      # Return the cumulative sums
      return(PO2)
    }

    # Get all the cumulative sums up to the biggest value needed
    cumsum <- outerConsecutive(biggest, s, n, m)

    # Get the CDF value for each input
    return(sapply(p, getEach, cumsum, s, m))
  }
  # If no inputs are less than 1, return the maximum possible overlap size
  # for all inputs
  else {
    return(rep(min(s, m), length(p)))
  }
}

#' @title Random number generator for conditional hypergeometric distributions
#'
#' @description Generates random numbers from a conditional hypergeometric distribution: the distribution of how many items are in the overlap of all samples when samples of arbitrary size are each taken without replacement from populations of arbitrary size.
#'
#' @param size the number of random numbers to generate
#' @param s an integer representing the size of the intersecting population
#' @param n a vector of integers representing the sizes of each non-intersecting population
#' @param m a vector of integers representing the sample sizes
#' @param verbose T/F should intermediate messages be printed?
#'
#' @return A vector of random numbers generated from the PMF of the conditional hypergeometric distribution specified by the parameters
#'
#' @examples rchyper(100, 10, c(12,13,14), c(7,8,9))
#'
#' @export rchyper


# The following implements a random number generator from a conditional
# hypergeometric distribution parameterized by the number of draws, the
# size of the overlapping population (an integer), the size of the
# unique populations (a vector of integers at least 2 long for 2 populations),
# and the number of people sampled from each population (a vector of integers
# at least 2 long for 2 samples).

rchyper <- function(size, s, n, m, verbose=T) {
  # Check that the input size is a positive integer
  if (size != floor(size) | !(size > 0)) {stop("Size must be a positive integer")}

  # Generate uniform random variables from 0 to 1 and send those to the quantile function
  inputs <- runif(size, 0, 1)
  return(qchyper(inputs, s, n, m, verbose))
}

#' @title Maximum likelihood estimator for overlap size in conditional hypergeometric distributions
#'
#' @description Calculates the MLE of the overlap size in a conditional hypergeometric distribution: the distribution of how many items are in the overlap of all samples when samples of arbitrary size are each taken without replacement from populations of arbitrary size.
#'
#' @param k the observed overlaps
#' @param n a vector of integers representing the sizes of each non-intersecting population
#' @param m a vector of integers representing the sample sizes
#' @param verbose T/F should intermediate messages be printed?
#'
#' @return The maximum likelihood estimator of the intersecting population size
#'
#' @examples mleS(c(0,0,1,1,0,2,0), c(12,13,14), c(7,8,9))
#'
#' @export mleS

# This function finds the maximum likelihood estimator of the overlapping
# population size given an input vector of observed sample overlaps
# (integers), the unique population sizes (a vector of integers), and the
# sample sizes from each population (a vector of integers).

mleS <- function(k, n, m, verbose=T) {
  # Check that inputs are valid

  # Check that no observed samples are larger than the sample sizes
  if (any(k > min(m))) {
    stop("Bad input: input larger than sample size")
  }
  # Check that observed samples are non-negative
  if (any(k < 0)) {
    stop("Bad input: input must be non-negative")
  }

  # Show it has started to run
  if (verbose) {message("Computing...")}

  # This function gets the log likelihood given observed values and a set
  # of parameters by multiplying all the PMFs together
  likelihood <- function(k, s, n, m, verbose) {
    sum(log(dchyper(k, s, n, m, verbose)))
  }

  # Placeholder to ensure a first value is calculated
  probCurrent <- -Inf

  # Get the smallest possible overlap population size
  index <- ifelse(all((m-n) < 0), 0, min((m - n)[(m - n) >= 0]))

  # Get the likelihood of this first overlap population size
  probNext <- likelihood(k, index, n, m, verbose)

  # The likelihood function will be monotonically decreasing or unimodal, so
  # increase the overlap size until the likelihood starts to decrease at
  # which point the maximum value has been reached the index before
  while (probCurrent < probNext | probCurrent == -Inf) {
    probCurrent <- probNext
    index <- index + 1
    probNext <- likelihood(k, index, n, m, verbose)
  }

  # Correct and return the index
  return(index - 1)
}

#' @title Maximum likelihood estimator for a unique population size in conditional hypergeometric distributions
#'
#' @description Calculates the MLE of a unique population size in a conditional hypergeometric distribution: the distribution of how many items are in the overlap of all samples when samples of arbitrary size are each taken without replacement from populations of arbitrary size.
#'
#' @param population the index of the unique population to estimate
#' @param k the observed overlaps
#' @param s an integer representing the size of the intersecting population
#' @param n a vector of integers representing the sizes of each non-intersecting population where the value of the unknown population size should be any integer as a placeholder
#' @param m a vector of integers representing the sample sizes
#' @param verbose T/F should intermediate messages be printed?
#'
#' @return The maximum likelihood estimator of the unknown unique population size
#'
#' @examples mleN(1, c(0,0,1,1,0,2,0), 8, c(0,13,14), c(7,8,9))
#'
#' @export mleN

# This function finds the maximum likelihood estimator of an unknown unique
# population size given the index of the unknown population size, an input
# vector of sample overlaps (integers), the overlapping population size, a
# vector of unique population sizes (a vector of integers where the unknown
# value should be any placeholder integer), and the sample sizes from each
# population (a vector of integers).

mleN <- function(population, k, s, n, m, verbose=T) {
  # Check that inputs are valid

  # Make sure the observed overlaps are not larger than the sample size
  if (any(k > min(m))) {stop("Bad input: input larger than sample size")}

  # Make sure the observed overlaps are not larger than the overlap size
  if (any(k > s)) {stop("Bad input: input larger than overlap size")}

  # Make sure the observed overlaps are non-negative
  if (any(k < 0)) {stop("Bad input: input must be non-negative")}

  # Make sure the unknown population index is valid
  if (population > length(n) | population <= 0) {stop("Invalid population selected")}

  # Make sure only one unknown population is selected
  if (length(population) != 1) {stop("Population length must be 1")}

  # If all inputs are 0, the most likely sample size is infinite
  if(all(k==0)) {return(Inf)}

  # Give status update
  if (verbose) {message("Computing...")}

  # This function gets the log likelihood given observed values and a set
  # of parameters by multiplying all the PMFs together
  likelihood <- function(population, n0, k, s, n, m, verbose) {
    n[population] <- n0
    sum(log(dchyper(k, s, n, m, verbose)))
  }

  # Placeholder to ensure a first value is calculated
  probCurrent <- -Inf

  # Get the smallest possible unique population size
  index <- max(m[population] - s, 0)

  # Get the likelihood of this first unique population size
  probNext <- likelihood(population, index, k, s, n, m, verbose)

  # The likelihood function will be monotonically decreasing or unimodal, so
  # increase the unique population size until the likelihood starts to decrease
  # at which point the maximum value has been reached the index before
  while (probCurrent < probNext | probCurrent == -Inf) {
    probCurrent <- probNext
    index <- index + 1
    probNext <- likelihood(population, index, k, s, n, m, verbose)
  }

  # Correct and return the index
  return(index - 1)
}

#' @title Maximum likelihood estimator for sample size in conditional hypergeometric distributions
#'
#' @description Calculates the MLE of a sample size in a conditional hypergeometric distribution: the distribution of how many items are in the overlap of all samples when samples of arbitrary size are each taken without replacement from populations of arbitrary size.
#'
#' @param population the index of the unknown sample size
#' @param k the observed overlaps
#' @param s an integer representing the size of the intersecting population
#' @param n a vector of integers representing the sizes of each non-intersecting population
#' @param m a vector of integers representing the sample sizes where the value of the unknown sample size should be any integer as a placeholder
#' @param verbose T/F should intermediate messages be printed?
#'
#' @return The maximum likelihood estimator of the unknown sample size
#'
#' @examples mleM(1, c(0,0,1,1,0,2,0), 8, c(12,13,14), c(0,8,9))
#'
#' @export mleM

# This function finds the maximum likelihood estimator of an unknown sample
# size given the index of the unknown sample size, an input vector of sample
# overlaps (integers), the overlapping population size, a vector of unique
# population sizes (integers), and the sample sizes from each population (a
# vector of integers where the unknown value should be any placeholder integer).

mleM <- function(population, k, s, n, m, verbose=T) {
  # Check that inputs are valid

  # Make sure the observed overlaps are not larger than the overlap size
  if (any(k > s)) {stop("Bad input: input larger than overlap size")}

  # Make sure the observed overlaps are non-negative
  if (any(k < 0)) {stop("Bad input: input must be non-negative")}

  # Make sure only one unknown population is selected
  if (length(population) != 1) {stop("Population length must be 1")}

  # Make sure the unknown sample size index is valid
  if (population > length(m) | population <= 0) {stop("Invalid population selected")}

  # Give status update
  if (verbose) {message("Computing...")}

  # This function gets the log likelihood given observed values and a set
  # of parameters by multiplying all the PMFs together
  likelihood <- function(population, m0, k, s, n, m, verbose) {
    m[population] <- m0
    sum(log(dchyper(k, s, n, m, verbose)))
  }

  # Placeholder to ensure a first value is calculated
  probCurrent <- -Inf

  # Get the smallest possible sample size
  index <- max(k)

  # Get the likelihood of this first sample size
  probNext <- likelihood(population, index, k, s, n, m, verbose)

  # The likelihood function will usually be monotonically decreasing or unimodal,
  # so increase the sample size until the likelihood starts to decrease
  # at which point the maximum value has been reached the index before.
  # However, sometimes the likelihood will monotonically increase until it reaches
  # the max value allowed by the population sizes, so break if that value is hit
  # Either way, the likelihood will never decrease unless a maximum likelihood
  # preceeds it.
  while (probCurrent < probNext | probCurrent == -Inf) {
    probCurrent <- probNext
    index <- index + 1
    if (index > s + n[population]) {
      break
    }
    probNext <- likelihood(population, index, k, s, n, m, verbose)
  }

  # Correct and return the index
  return(index - 1)
}


#' @title P-values from a conditional hypergeometric distribution
#'
#' @description Calculates p-values from a conditional hypergeometric distribution: the distribution of how many items are in the overlap of all samples when samples of arbitrary size are each taken without replacement from populations of arbitrary size.
#'
#' @param k an integer or vector of integers representing the overlap size
#' @param s an integer representing the size of the intersecting population
#' @param n a vector of integers representing the sizes of each non-intersecting population
#' @param m a vector of integers representing the sample sizes
#' @param tail whether the p-value should be from the upper or lower tail (options: "upper", "lower")
#' @param verbose T/F should intermediate messages be printed?
#'
#' @return The probability of getting the k or more (or less if tail="lower") overlaps by chance from the conditional hypergeometric distribution specified by the parameters
#'
#' @examples pvalchyper(c(1,2), 8, c(12,13,14), c(7,8,9), "upper")
#'
#' @export pvalchyper

# This function finds the probability of getting a number of overlaps
# equal to or greater than (or equal to or less then for tail="lower")
# the observed value (or values) k.  It is parameterized by a vector of
# observed values, the size of the overlapping population, the unique
# population sizes (a vector of integers), the sample sizes (a vector of
# integers), and which tail the p-value should be taken from ("lower" or
# "upper").

pvalchyper <- function(k, s, n, m, tail="upper", verbose=T) {
  # Check the inputs

  # Check that the inputs are valid

  # Ensure at least two unique populations are given
  if (length(n) < 2) {stop("Need two regions (n too short)")}

  # Check that only one overlapping population is given
  if (length(s) > 1) {stop("Only 1 intersection region (s too long)")}

  # Check that the number of unique populations and the number of samples are the same
  if (length(m) != length(n)) {stop("n and m must match in length")}

  # Check that all population and sample sizes are non-negative
  if (s < 0 | any(n < 0) | any(m < 0)) {stop("Only positive values allowed")}

  # Check that sample sizes are not bigger than the populations they are from
  if (any(m > (n + s))) {stop("Sample size cannot be bigger than population")}

  # Check that all inputs are integers
  if (!all(c(k, s, n, m) == floor(c(k, s, n, m)))) {stop("Inputs must be integers")}

  # Check that the tail is a valid input
  if (!(tail=="upper" | tail=="lower")) {stop("Tail must be upper or lower")}

  # If requesting the upper tail, get the CDF of the value one below the requested
  if (tail=="upper") {
    temp <- k - 1

    # If the requested value is 0 the upper tail density is 1, otherwise it is
    # the CDF
    output <- ifelse(k==0, 1, 1 - pchyper(temp, s, n, m, verbose))

    # Sometimes the tail density falls slightly below 0 because of rounding errors,
    # so return the smallest value R can store as the p-value in that case
    return(ifelse(output < 0, 10^(-16), output))

  # If requesting the lower tail, just return the CDF at that index and make sure
  # it is not negative.
  } else {
    output <- pchyper(k, s, n, m)
    return(ifelse(output < 0, 10^(-16), output))
  }

}
