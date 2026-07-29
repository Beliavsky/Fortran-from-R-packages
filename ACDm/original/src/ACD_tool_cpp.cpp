#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export(name = ".group_running_less_than")]]
NumericVector group_running_less_than(NumericVector x, 
                             double c,
                             bool group_zeores = false){
  // x is vector of durations
  
  int N = x.size();
  NumericVector return_vector(N);
  int group_counter = 1;
  
  if(group_zeores == false){
    if(x(0) < c){
      return_vector(0) = group_counter;
    } else{
      return_vector(0) = 0;
    }
    
    for(int i = 1; i < N ; i++){
      
      if(x(i) < c){
        return_vector(i) = group_counter;
        if(return_vector(i - 1) == 0) return_vector(i - 1) = group_counter;
      } else{
        return_vector(i) = 0;
        if(return_vector(i - 1) != 0) group_counter++;
      }
      
      if (i % 10000 == 0) Rcpp::checkUserInterrupt();
    }
  } else if(group_zeores == true){
    
    return_vector(0) = 1;
    
    for(int i = 1; i < N ; i++){
      if(x(i) >= c) group_counter++;
      return_vector(i) = group_counter;
    }
  }
  
  return(return_vector);
}

// [[Rcpp::export(name = ".time_rolling_ntrans_left")]]
SEXP time_rolling_ntrans_left(NumericVector x,
                            double window_size){
  
  // calculates a rolling mean, 
  // where the window is a time interval instead of number of close observations
  // this function uses a window to the left of each observation
  //
  // x must be  a vector of durations
  
  int N = x.size();
  NumericVector time(N); //cumsum of x
  NumericVector return_vector(N);
  
  NumericVector first_id(N); 
  
  // calculates the 'time' variable as the cumulative sum of the durations (less the first duration): 
  time(0) = 0;
  for(int i = 1; i < N ; i++){
    time(i) = time(i - 1) + x(i);
  }
  
  // the index of the observation most left in the window for the last observation:
  double last_incl_id = 0; 
  
  int last_na = 0;
  return_vector(0) = NA_REAL;
  
  // sets the first observations, with a time less than the windowsize, to NA:
  for(int i = 1; i < N; i++){
    if(time(i) > window_size){
      last_na = i - 1;
      break;
    }
    return_vector(i) = NA_REAL;
  }
  
  for(int i = last_na + 1; i < N; i++){
    
    double time_i = time(i);
    
    for(int i2 = last_incl_id; i2 <= i; i2++){
      if(time_i - time(i2) < window_size){
        last_incl_id = i2;
        break;
      } 
    }
    
    return_vector(i) = i - last_incl_id + 1;
    
  }
  
  return(return_vector);
}

// [[Rcpp::export(name = ".time_rolling_ntrans_right")]]
SEXP time_rolling_ntrans_right(NumericVector x,
                            double window_size){
  
  // calculates a rolling mean, 
  // where the window is a time interval instead of number of close observations
  // this function uses a window to the right of each observation
  //
  // x must be  a vector of durations
  
  int N = x.size();
  NumericVector time(N); //cumsum of x
  NumericVector return_vector(N);
  
  NumericVector first_id(N); 
  
  // calculates the 'time' variable as the cumulative sum of the durations (less the first duration): 
  time(0) = 0;
  for(int i = 1; i < N ; i++){
    time(i) = time(i - 1) + x(i);
  }
  
  // the index of the observation most right in the window for the last observation:
  double last_incl_id = 0; 
  
  int first_na = 0;
  
  // sets the first observations, with a time less than the windowsize, to NA:
  for(int i = N - 1; i > 0; i--){
    if(time(N - 1) - time(i) > window_size){
      first_na = i + 1;
      break;
    }
    return_vector(i) = NA_REAL;
  }

  for(int i = 0; i < first_na; i++){
    
    double time_i = time(i);
    
    for(int i2 = last_incl_id; i2 < N; i2++){
      if(time(i2) - time_i > window_size){
        last_incl_id = i2 - 1;
        break;
      } 
    }
    
    return_vector(i) = last_incl_id - i + 1;
    
  }
  
  return(return_vector);
}

// [[Rcpp::export(name = ".roll_time_mean")]]
NumericVector roll_time_mean(NumericVector x,
                             double width,
                             StringVector align){
  
  // calculates a rolling mean, 
  // where the window is a time interval instead of number of close observations
  // this function uses a window to the right of each observation
  //
  // x must be  a vector of durations
  
  
  int N = x.size();
  NumericVector time(N); //cumsum of x
  NumericVector return_vector(N);
  
  NumericVector first_id(N); 
  
  // calculates the 'time' variable as the cumulative sum of the durations (less the first duration): 
  time(0) = 0;
  for(int i = 1; i < N ; i++){
    time(i) = time(i - 1) + x(i);
  }
  
  
  if(align(0) == "right"){
    
    // the index of the observation most right in the window for the last observation:
    double max_incl_id = 0; 
    
    int first_na = 0;
    
    // sets the first observations, with a time less than the windowsize, to NA:
    for(int i = N - 1; i > 0; i--){
      if(time(N - 1) - time(i) > width){
        first_na = i + 1;
        break;
      }
      return_vector(i) = NA_REAL;
    }
    
    for(int i = 0; i < first_na; i++){
      
      double time_i = time(i);
      
      for(int i2 = max_incl_id; i2 < N; i2++){
        if(time(i2) - time_i > width){
          max_incl_id = i2 - 1;
          break;
        } 
      }
      
      // return_vector(i) = max_incl_id - i + 1;
      
      return_vector(i) = 0;
      for(int i2 = i; i2 <= max_incl_id; i2++){
        return_vector(i) += x(i2);
      }
      return_vector(i) /= max_incl_id - i + 1;
      
    }
  } else if(align(0) == "left"){
    // the index of the observation most left in the window for the last observation:
    double min_incl_id = 0; 
    
    int last_na = 0;
    return_vector(0) = NA_REAL;
    
    // sets the first observations, with a time less than the windowsize, to NA:
    for(int i = 1; i < N; i++){
      if(time(i) > width){
        last_na = i - 1;
        break;
      }
      return_vector(i) = NA_REAL;
    }
    
    for(int i = last_na + 1; i < N; i++){
      
      double time_i = time(i);
      
      for(int i2 = min_incl_id; i2 <= i; i2++){
        if(time_i - time(i2) < width){
          min_incl_id = i2;
          break;
        } 
      }
      
      // return_vector(i) = i - min_incl_id + 1;
      
      return_vector(i) = 0;
      for(int i2 = min_incl_id; i2 <= i; i2++){
        return_vector(i) += x(i2);
      }
      return_vector(i) /= i - min_incl_id + 1;
    }
  } else if(align(0) == "center"){
    
    double half_width = width / 2;
    
    int first_na = 0;
    int last_na = 0;
    return_vector(0) = NA_REAL;
    // sets the first observations, with a time less than the windowsize, to NA:
    for(int i = 1; i < N; i++){
      if(time(i) > half_width){
        last_na = i - 1;
        break;
      }
      return_vector(i) = NA_REAL;
    }
    // sets the first observations, with a time less than the windowsize, to NA:
    for(int i = N - 1; i > 0; i--){
      if(time(N - 1) - time(i) > half_width){
        first_na = i + 1;
        break;
      }
      return_vector(i) = NA_REAL;
    }
    
    double min_incl_id = 0; 
    double max_incl_id = 0; 
    
    for(int i = last_na + 1; i < first_na; i++){
      
      double time_i = time(i);
      
      // get the smallest id that is within the window:
      for(int i2 = min_incl_id; i2 <= i; i2++){
        if(time_i - time(i2) < half_width){
          min_incl_id = i2;
          break;
        } 
      }
      // get the margest id that is within the window:
      for(int i2 = max_incl_id; i2 < N; i2++){
        if(time(i2) - time_i > half_width){
          max_incl_id = i2 - 1;
          break;
        } 
      }
      
      // return_vector(i) = i - last_incl_id + 1;
      
      return_vector(i) = 0;
      for(int i2 = min_incl_id; i2 <= max_incl_id; i2++){
        return_vector(i) += x(i2);
      }
      return_vector(i) /= max_incl_id - min_incl_id + 1;
    }
    
  }
  
  return(return_vector);
}

// [[Rcpp::export(name = ".roll_time_count")]]
NumericVector roll_time_count(NumericVector x,
                             double width,
                             StringVector align){
  
  // calculates a rolling mean, 
  // where the window is a time interval instead of number of close observations
  // this function uses a window to the right of each observation
  //
  // x must be  a vector of durations
  
  
  int N = x.size();
  NumericVector time(N); //cumsum of x
  NumericVector return_vector(N);
  
  NumericVector first_id(N); 
  
  // calculates the 'time' variable as the cumulative sum of the durations (less the first duration): 
  time(0) = 0;
  for(int i = 1; i < N ; i++){
    time(i) = time(i - 1) + x(i);
  }
  
  
  if(align(0) == "right"){
    
    // the index of the observation most right in the window for the last observation:
    double max_incl_id = 0; 
    
    int first_na = 0;
    
    // sets the first observations, with a time less than the windowsize, to NA:
    for(int i = N - 1; i > 0; i--){
      if(time(N - 1) - time(i) > width){
        first_na = i + 1;
        break;
      }
      return_vector(i) = NA_REAL;
    }
    
    for(int i = 0; i < first_na; i++){
      
      double time_i = time(i);
      
      for(int i2 = max_incl_id; i2 < N; i2++){
        if(time(i2) - time_i > width){
          max_incl_id = i2 - 1;
          break;
        } 
      }
      
      return_vector(i) = max_incl_id - i + 1;
    }
  } else if(align(0) == "left"){
    // the index of the observation most left in the window for the last observation:
    double min_incl_id = 0; 
    
    int last_na = 0;
    return_vector(0) = NA_REAL;
    
    // sets the first observations, with a time less than the windowsize, to NA:
    for(int i = 1; i < N; i++){
      if(time(i) > width){
        last_na = i - 1;
        break;
      }
      return_vector(i) = NA_REAL;
    }
    
    for(int i = last_na + 1; i < N; i++){
      
      double time_i = time(i);
      
      for(int i2 = min_incl_id; i2 <= i; i2++){
        if(time_i - time(i2) < width){
          min_incl_id = i2;
          break;
        } 
      }
      
      return_vector(i) = i - min_incl_id + 1;
      
    }
  } else if(align(0) == "center"){
    
    double half_width = width / 2;
    
    int first_na = 0;
    int last_na = 0;
    return_vector(0) = NA_REAL;
    // sets the first observations, with a time less than the windowsize, to NA:
    for(int i = 1; i < N; i++){
      if(time(i) > half_width){
        last_na = i - 1;
        break;
      }
      return_vector(i) = NA_REAL;
    }
    // sets the first observations, with a time less than the windowsize, to NA:
    for(int i = N - 1; i > 0; i--){
      if(time(N - 1) - time(i) > half_width){
        first_na = i + 1;
        break;
      }
      return_vector(i) = NA_REAL;
    }
    
    double min_incl_id = 0; 
    double max_incl_id = 0; 
    
    for(int i = last_na + 1; i < first_na; i++){
      
      double time_i = time(i);
      
      // get the smallest id that is within the window:
      for(int i2 = min_incl_id; i2 <= i; i2++){
        if(time_i - time(i2) < half_width){
          min_incl_id = i2;
          break;
        } 
      }
      // get the margest id that is within the window:
      for(int i2 = max_incl_id; i2 < N; i2++){
        if(time(i2) - time_i > half_width){
          max_incl_id = i2 - 1;
          break;
        } 
      }
      
      return_vector(i) = max_incl_id - min_incl_id + 1;
      
    }
    
  }
  
  return(return_vector);
}

