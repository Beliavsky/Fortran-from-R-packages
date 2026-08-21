pso <-
  function(func,
           S = 350,
           lim_inf,
           lim_sup,
           e = 0.0001,
           data = NULL,
           N = 500, prop = 0.2) {
    b_lo = min(lim_inf)
    b_up = max(lim_sup)
    integer_max = .Machine$integer.max
    
    if (length(lim_sup) != length(lim_inf)) {
      stop("The vectors lim_inf and lim_sup must have the same dimension.")
    }
    dimension = length(lim_sup)
    swarm_xi = swarm_pi = swarm_vi = matrix(NA, nrow = S, ncol = dimension)
    
    # The best position of the particles.
    g = runif(n = dimension, min = lim_inf, max = lim_sup)
    
    # Objective function calculated in g.
    f_g = func(par = as.vector(g), x = as.vector(data))
    
    if (NaN %in% f_g == TRUE || Inf %in% abs(f_g) == TRUE) {
      while (NaN %in% f_g == TRUE || Inf %in% abs(f_g) == TRUE) {
        g = runif(n = dimension, min = lim_inf, max = lim_sup)
        f_g = func(par = g, x = as.vector(data))
      }
    }
    
    # Here begins initialization of the algorithm.
    x_i = mapply(runif,
                 n = S,
                 min = lim_inf,
                 max = lim_sup)
    
    # Initializing the best position of particularities i to initial position.
    
    swarm_pi = swarm_xi = x_i
    f_pi = apply(
      X = x_i,
      MARGIN = 1,
      FUN = func,
      x = as.vector(data)
    )
    
    is.integer0 <- function(x) {
      is.integer(x) && length(x) == 0L
    }
    
    if (NaN %in% f_pi == TRUE || Inf %in% abs(f_pi)) {
      while (NaN %in% f_pi == TRUE || Inf %in% abs(f_pi)) {
        id_inf_fpi = which(abs(f_pi) == Inf)
        if (is.integer0(id_inf_fpi) != TRUE) {
          f_pi[id_inf_fpi] = integer_max
        }
        id_nan_fpi = which(f_pi == NaN)
        if (is.integer0(id_nan_fpi) != TRUE) {
          x_i[id_nan_fpi,] = mapply(runif,
                                    n = length(id_nan_fpi),
                                    min = lim_inf,
                                    max = lim_sup)
          swarm_pi = swarm_xi = x_i
          f_pi = apply(
            X = x_i,
            MARGIN = 1,
            FUN = func,
            x = as.vector(data)
          )
        }
      }
    }
     
    minimo_fpi = min(f_pi)
    
    if (minimo_fpi < f_g)
      g = x_i[which.min(f_pi),]
    
    # Initializing the speeds of the particles.
    
    swarm_vi = mapply(runif,
                      n = S,
                      min = -abs(rep(abs(b_up - b_lo), dimension)),
                      max = abs(rep(abs(b_up - b_lo), dimension)))
    
    # Here ends the initialization of the algorithm
    
    omega = 0.5
    phi_p = 0.5
    phi_g = 0.5
    
    m = 1
    vector_f_g <- vector()
    #vector_par <- vector()
    while (1)  {
      # r_p and r_g are randomized numbers in (0.1).
      r_p = runif(n = dimension, min = 0, max = 1)
      r_g = runif(n = dimension, min = 0, max = 1)
      
      # Updating the vector speed.
      swarm_vi = omega * swarm_vi + phi_p * r_p * (swarm_pi - swarm_xi) +
        phi_g * r_g * (g - swarm_xi)
      
      # Updating the position of each particle.
      swarm_xi = swarm_xi + swarm_vi
      
      myoptim = function(...)
        tryCatch(
          optim(...),
          error = function(e)
            NA
        )
      
      f_xi = apply(
        X = as.matrix(swarm_xi),
        MARGIN = 1,
        FUN = func,
        x = as.vector(data)
      )
      f_pi = apply(
        X = as.matrix(swarm_pi),
        MARGIN = 1,
        FUN = func,
        x = as.vector(data)
      )
      f_g =  func(par = g, x = as.vector(data))
      #vector_par[m] = g 
      if (NaN %in% f_xi == TRUE || NaN %in% f_pi == TRUE) {
        while (NaN %in% f_xi == TRUE) {
          id_comb = c(which(is.na(f_xi) == TRUE), which(is.na(f_pi) == TRUE))
          if (is.integer0(id_comb) != TRUE) {
            new_xi = mapply(runif,
                            n = length(id_comb),
                            min = lim_inf,
                            max = lim_sup)
            swarm_pi[id_comb,] = swarm_xi[id_comb,] = new_xi
            if (length(id_comb) > 1) {
              f_xi[id_comb] = apply(
                X = as.matrix(swarm_xi[id_comb,]),
                MARGIN = 1,
                FUN = func,
                x = as.vector(data)
              )
              f_pi[id_comb] = apply(
                X = as.matrix(swarm_pi[id_comb,]),
                MARGIN = 1,
                FUN = func,
                x = as.vector(data)
              )
            } else{
              f_xi[id_comb] = func(par = new_xi, x = as.vector(data))
            }
          }
        }
      }
      
      if (Inf %in% abs(f_xi) == TRUE) {
        f_xi[which(is.infinite(f_xi))] = integer_max
      }
      if (Inf %in% abs(f_pi) == TRUE) {
        f_pi[which(is.infinite(f_pi))] = integer_max
      }
      
      # There are values below the lower limit of restrictions?
      id_test_inf =
        which(apply(swarm_xi < t(matrix(
          rep(lim_inf, S), dimension, S
        )), 1, sum) >= 1)
      id_test_sup =
        which(apply(swarm_xi > t(matrix(
          rep(lim_sup, S), dimension, S
        )), 1, sum) >= 1)
      
      if (is.integer0(id_test_inf) != TRUE) {
        swarm_pi[id_test_inf,] = swarm_xi[id_test_inf,] =
          mapply(runif,
                 n = length(id_test_inf),
                 min = lim_inf,
                 max = lim_sup)
      }
      
      if (is.integer0(id_test_sup) != TRUE) {
        swarm_pi[id_test_sup,] = swarm_xi[id_test_sup,] =
          mapply(runif,
                 n = length(id_test_sup),
                 min = lim_inf,
                 max = lim_sup)
      }
      
      if (is.integer0(which((f_xi <= f_pi) == TRUE))) {
        swarm_pi[which((f_xi <= f_pi)),] = swarm_pi[which((f_xi <= f_pi)),]
      }
      
      if (f_xi[which.min(f_xi)] <= f_pi[which.min(f_pi)]) {
        swarm_pi[which.min(f_pi),] = swarm_xi[which.min(f_xi),]
        if (f_pi[which.min(f_pi)] < f_g)
          g = swarm_pi[which.min(f_pi),]
      } # Here ends the block if.
      
      vector_f_g[m] = f_g
      m = m + 1

      if(length(vector_f_g)>=N){
        n_var = ceiling(length(vector_f_g)*prop)
        id_var = seq(length(vector_f_g),length(vector_f_g)-n_var,-1)
        if(var(vector_f_g[id_var])<=e){
          break
        }
      }
      
    } # Here ends the block while.
    
    f_x = apply(
      X = swarm_xi,
      MARGIN = 1,
      FUN = func,
      x = as.vector(data)
    )
    list(par = g, f = vector_f_g)
    
  } # Here ends the function.