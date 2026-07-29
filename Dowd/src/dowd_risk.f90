! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module dowd_risk
  use dowd_kinds, only: dp, pi
  use dowd_math, only: mean_value, sample_sd, sorted_copy, quantile_linear, &
       normal_pdf, normal_cdf, normal_quantile, student_t_pdf, student_t_quantile, &
       random_normal, sort_in_place
  implicit none
  private

  integer, parameter, public :: kernel_gaussian = 1
  integer, parameter, public :: kernel_box = 2
  integer, parameter, public :: kernel_triangular = 3
  integer, parameter, public :: kernel_epanechnikov = 4

  public :: normal_var, normal_es, student_t_var, student_t_es
  public :: lognormal_var, lognormal_es, log_student_t_var, log_student_t_es
  public :: historical_var, historical_es, historical_var_es
  public :: cornish_fisher_var, cornish_fisher_es
  public :: gumbel_var, gumbel_es, frechet_var, frechet_es
  public :: generalized_pareto_var, generalized_pareto_es
  public :: hill_estimator, hill_quantile_estimator, pickands_estimator
  public :: kernel_var, kernel_es, kde_bandwidth_nrd
  public :: bootstrap_var_es, bootstrap_confidence_interval
  public :: boxcox_lambda, boxcox_var, boxcox_es
  public :: spectral_risk_normal, quantile_standard_error_normal
  public :: quantile_standard_error_t

contains

  pure subroutine validate_probability(cl)
    real(dp), intent(in) :: cl
    if (cl <= 0.0_dp .or. cl >= 1.0_dp) error stop "confidence level must be in (0,1)"
  end subroutine validate_probability

  pure real(dp) function normal_var(mu, sigma, cl, hp) result(value)
    real(dp), intent(in) :: mu, sigma, cl, hp
    call validate_probability(cl)
    if (sigma < 0.0_dp .or. hp <= 0.0_dp) error stop "normal_var: invalid sigma or holding period"
    value = -mu*hp - sigma*sqrt(hp)*normal_quantile(1.0_dp-cl)
  end function normal_var

  pure real(dp) function normal_es(mu, sigma, cl, hp) result(value)
    real(dp), intent(in) :: mu, sigma, cl, hp
    real(dp) :: z
    call validate_probability(cl)
    if (sigma < 0.0_dp .or. hp <= 0.0_dp) error stop "normal_es: invalid sigma or holding period"
    z = normal_quantile(cl)
    value = sigma*sqrt(hp)*normal_pdf(z)/(1.0_dp-cl)-mu*hp
  end function normal_es

  real(dp) function student_t_var(mu, sigma, nu, cl, hp) result(value)
    real(dp), intent(in) :: mu, sigma, nu, cl, hp
    real(dp) :: scale
    call validate_probability(cl)
    if (sigma < 0.0_dp .or. nu <= 2.0_dp .or. hp <= 0.0_dp) &
      error stop "student_t_var: invalid parameters"
    scale = sigma*sqrt((nu-2.0_dp)/nu)
    value = -mu*hp - scale*sqrt(hp)*student_t_quantile(1.0_dp-cl,nu)
  end function student_t_var

  real(dp) function student_t_es(mu, sigma, nu, cl, hp) result(value)
    real(dp), intent(in) :: mu, sigma, nu, cl, hp
    real(dp) :: q, scale, tail_mean
    call validate_probability(cl)
    if (sigma < 0.0_dp .or. nu <= 2.0_dp .or. hp <= 0.0_dp) &
      error stop "student_t_es: invalid parameters"
    q = student_t_quantile(cl,nu)
    scale = sigma*sqrt((nu-2.0_dp)/nu)
    tail_mean = (nu+q*q)/(nu-1.0_dp)*student_t_pdf(q,nu)/(1.0_dp-cl)
    value = scale*sqrt(hp)*tail_mean-mu*hp
  end function student_t_es

  pure real(dp) function lognormal_var(investment, mu, sigma, cl, hp) result(value)
    real(dp), intent(in) :: investment, mu, sigma, cl, hp
    real(dp) :: lower_return
    call validate_probability(cl)
    if (investment < 0.0_dp .or. sigma < 0.0_dp .or. hp <= 0.0_dp) &
      error stop "lognormal_var: invalid parameters"
    lower_return = exp(mu*hp+sigma*sqrt(hp)*normal_quantile(1.0_dp-cl))
    value = investment*(1.0_dp-lower_return)
  end function lognormal_var

  pure real(dp) function lognormal_es(investment, mu, sigma, cl, hp) result(value)
    real(dp), intent(in) :: investment, mu, sigma, cl, hp
    real(dp) :: alpha, z, conditional_wealth
    call validate_probability(cl)
    if (investment < 0.0_dp .or. sigma < 0.0_dp .or. hp <= 0.0_dp) &
      error stop "lognormal_es: invalid parameters"
    alpha = 1.0_dp-cl
    z = normal_quantile(alpha)
    conditional_wealth = investment*exp(mu*hp+0.5_dp*sigma*sigma*hp) * &
      normal_cdf(z-sigma*sqrt(hp))/alpha
    value = investment-conditional_wealth
  end function lognormal_es

  real(dp) function log_student_t_var(investment, mu, sigma, nu, cl, hp) result(value)
    real(dp), intent(in) :: investment, mu, sigma, nu, cl, hp
    real(dp) :: scale, lr
    call validate_probability(cl)
    if (nu <= 2.0_dp) error stop "log_student_t_var: nu must exceed 2"
    scale = sigma*sqrt((nu-2.0_dp)/nu)
    lr = mu*hp+scale*sqrt(hp)*student_t_quantile(1.0_dp-cl,nu)
    value = investment*(1.0_dp-exp(lr))
  end function log_student_t_var

  real(dp) function log_student_t_es(investment, mu, sigma, nu, cl, hp, n_slices) result(value)
    real(dp), intent(in) :: investment, mu, sigma, nu, cl, hp
    integer, intent(in), optional :: n_slices
    integer :: i, n
    real(dp) :: p, alpha, scale, wealth_sum
    call validate_probability(cl)
    n = 2000
    if (present(n_slices)) n = max(100,n_slices)
    alpha = 1.0_dp-cl
    scale = sigma*sqrt((nu-2.0_dp)/nu)
    wealth_sum = 0.0_dp
    do i = 1, n
      p = alpha*(real(i,dp)-0.5_dp)/real(n,dp)
      wealth_sum = wealth_sum+investment*exp(mu*hp+scale*sqrt(hp)*student_t_quantile(p,nu))
    end do
    value = investment-wealth_sum/real(n,dp)
  end function log_student_t_es

  real(dp) function historical_var(profit_loss, cl) result(value)
    real(dp), intent(in) :: profit_loss(:), cl
    call validate_probability(cl)
    value = quantile_linear(-profit_loss,cl)
  end function historical_var

  real(dp) function historical_es(profit_loss, cl, n_slices) result(value)
    real(dp), intent(in) :: profit_loss(:), cl
    integer, intent(in), optional :: n_slices
    real(dp), allocatable :: losses(:)
    real(dp) :: p, total
    integer :: i, n
    call validate_probability(cl)
    n = max(1000,size(profit_loss)*4)
    if (present(n_slices)) n = max(100,n_slices)
    losses = -profit_loss
    total = 0.0_dp
    do i = 1, n
      p = cl+(1.0_dp-cl)*(real(i,dp)-0.5_dp)/real(n,dp)
      total = total+quantile_linear(losses,p)
    end do
    value = total/real(n,dp)
  end function historical_es

  subroutine historical_var_es(profit_loss, cl, var_value, es_value)
    real(dp), intent(in) :: profit_loss(:), cl
    real(dp), intent(out) :: var_value, es_value
    var_value = historical_var(profit_loss,cl)
    es_value = historical_es(profit_loss,cl)
  end subroutine historical_var_es

  pure real(dp) function cornish_fisher_adjustment(z, skewness, kurtosis) result(a)
    real(dp), intent(in) :: z, skewness, kurtosis
    a = (z*z-1.0_dp)*skewness/6.0_dp + &
        (z**3-3.0_dp*z)*(kurtosis-3.0_dp)/24.0_dp - &
        (2.0_dp*z**3-5.0_dp*z)*skewness*skewness/36.0_dp
  end function cornish_fisher_adjustment

  pure real(dp) function cornish_fisher_var(mu, sigma, skewness, kurtosis, cl) result(value)
    real(dp), intent(in) :: mu, sigma, skewness, kurtosis, cl
    real(dp) :: z
    call validate_probability(cl)
    z = normal_quantile(1.0_dp-cl)
    value = -sigma*(z+cornish_fisher_adjustment(z,skewness,kurtosis))-mu
  end function cornish_fisher_var

  real(dp) function cornish_fisher_es(mu, sigma, skewness, kurtosis, cl, n_slices) result(value)
    real(dp), intent(in) :: mu, sigma, skewness, kurtosis, cl
    integer, intent(in), optional :: n_slices
    integer :: i, n
    real(dp) :: p, total
    call validate_probability(cl)
    n = 1000
    if (present(n_slices)) n = max(100,n_slices)
    total = 0.0_dp
    do i = 1, n
      p = cl+(1.0_dp-cl)*(real(i,dp)-0.5_dp)/real(n,dp)
      total = total+cornish_fisher_var(mu,sigma,skewness,kurtosis,p)
    end do
    value = total/real(n,dp)
  end function cornish_fisher_es

  pure real(dp) function gumbel_var(mu, sigma, block_size, cl) result(value)
    real(dp), intent(in) :: mu, sigma, block_size, cl
    call validate_probability(cl)
    if (sigma < 0.0_dp .or. block_size <= 0.0_dp) error stop "gumbel_var: invalid parameters"
    value = mu-sigma*log(-block_size*log(cl))
  end function gumbel_var

  real(dp) function gumbel_es(mu, sigma, block_size, cl, n_slices) result(value)
    real(dp), intent(in) :: mu, sigma, block_size, cl
    integer, intent(in), optional :: n_slices
    integer :: i, n
    real(dp) :: p
    n = 1000
    if (present(n_slices)) n = max(100,n_slices)
    value = 0.0_dp
    do i = 1, n
      p = cl+(1.0_dp-cl)*(real(i,dp)-0.5_dp)/real(n,dp)
      value = value+gumbel_var(mu,sigma,block_size,p)
    end do
    value = value/real(n,dp)
  end function gumbel_es

  pure real(dp) function frechet_var(mu, sigma, tail_index, block_size, cl) result(value)
    real(dp), intent(in) :: mu, sigma, tail_index, block_size, cl
    call validate_probability(cl)
    if (sigma < 0.0_dp .or. tail_index <= 0.0_dp .or. block_size <= 0.0_dp) &
      error stop "frechet_var: invalid parameters"
    value = mu-(sigma/tail_index)*(1.0_dp-(-block_size*log(cl))**(-tail_index))
  end function frechet_var

  real(dp) function frechet_es(mu, sigma, tail_index, block_size, cl, n_slices) result(value)
    real(dp), intent(in) :: mu, sigma, tail_index, block_size, cl
    integer, intent(in), optional :: n_slices
    integer :: i, n
    real(dp) :: p
    n = 1000
    if (present(n_slices)) n = max(100,n_slices)
    value = 0.0_dp
    do i = 1, n
      p = cl+(1.0_dp-cl)*(real(i,dp)-0.5_dp)/real(n,dp)
      value = value+frechet_var(mu,sigma,tail_index,block_size,p)
    end do
    value = value/real(n,dp)
  end function frechet_es

  real(dp) function generalized_pareto_var(data, beta, xi, threshold_probability, cl) result(value)
    real(dp), intent(in) :: data(:), beta, xi, threshold_probability, cl
    real(dp), allocatable :: x(:)
    real(dp) :: u
    integer :: n, nu, idx
    call validate_probability(cl)
    if (beta <= 0.0_dp .or. threshold_probability <= 0.0_dp .or. threshold_probability >= 1.0_dp) &
      error stop "generalized_pareto_var: invalid parameters"
    x = sorted_copy(data)
    n = size(x)
    nu = max(1,min(n-1,int(floor(threshold_probability*real(n,dp)))))
    idx = max(1,n-nu)
    u = x(idx)
    if (abs(xi) < sqrt(epsilon(1.0_dp))) then
      value = u-beta*log((1.0_dp-cl)/threshold_probability)
    else
      value = u+(beta/xi)*(((1.0_dp-cl)/threshold_probability)**(-xi)-1.0_dp)
    end if
  end function generalized_pareto_var

  real(dp) function generalized_pareto_es(data, beta, xi, threshold_probability, cl) result(value)
    real(dp), intent(in) :: data(:), beta, xi, threshold_probability, cl
    real(dp) :: var_value
    if (xi >= 1.0_dp) error stop "generalized_pareto_es: ES is infinite for xi >= 1"
    var_value = generalized_pareto_var(data,beta,xi,threshold_probability,cl)
    value = (var_value+beta-xi*generalized_pareto_threshold(data,threshold_probability))/(1.0_dp-xi)
  end function generalized_pareto_es

  real(dp) function generalized_pareto_threshold(data, threshold_probability) result(u)
    real(dp), intent(in) :: data(:), threshold_probability
    real(dp), allocatable :: x(:)
    integer :: n, nu
    x = sorted_copy(data)
    n = size(x)
    nu = max(1,min(n-1,int(floor(threshold_probability*real(n,dp)))))
    u = x(max(1,n-nu))
  end function generalized_pareto_threshold

  real(dp) function hill_estimator(data, tail_size) result(value)
    real(dp), intent(in) :: data(:)
    integer, intent(in) :: tail_size
    real(dp), allocatable :: x(:)
    integer :: n, k, i
    x = pack(data,data>0.0_dp)
    if (size(x) < 2) error stop "hill_estimator: at least two positive observations required"
    call sort_in_place(x)
    n = size(x)
    k = min(max(1,tail_size),n-1)
    value = 0.0_dp
    do i = n-k+1, n
      value = value+log(x(i)/x(n-k))
    end do
    value = value/real(k,dp)
  end function hill_estimator

  real(dp) function hill_quantile_estimator(data, tail_index, in_sample_probability, cl) result(value)
    real(dp), intent(in) :: data(:), tail_index, in_sample_probability, cl
    real(dp), allocatable :: x(:)
    integer :: n, k
    call validate_probability(cl)
    x = sorted_copy(data)
    n = size(x)
    k = max(1,min(n-1,int(floor(in_sample_probability*real(n,dp)))))
    value = x(n-k)*(real(n,dp)*(1.0_dp-cl)/real(k,dp))**(-tail_index)
  end function hill_quantile_estimator

  real(dp) function pickands_estimator(data, tail_size) result(value)
    real(dp), intent(in) :: data(:)
    integer, intent(in) :: tail_size
    real(dp), allocatable :: x(:)
    real(dp) :: numerator, denominator
    integer :: n, k
    x = sorted_copy(data)
    n = size(x)
    k = tail_size
    if (k < 1 .or. 4*k >= n) error stop "pickands_estimator: tail_size must satisfy 4*k < n"
    numerator = x(n-k)-x(n-2*k)
    denominator = x(n-2*k)-x(n-4*k)
    if (numerator <= 0.0_dp .or. denominator <= 0.0_dp) error stop "pickands_estimator: degenerate spacings"
    value = log(numerator/denominator)/log(2.0_dp)
  end function pickands_estimator

  real(dp) function kde_bandwidth_nrd(x) result(h)
    real(dp), intent(in) :: x(:)
    real(dp) :: s, iqr, scale
    integer :: n
    n = size(x)
    if (n < 2) error stop "kde_bandwidth_nrd: insufficient data"
    s = sample_sd(x)
    iqr = quantile_linear(x,0.75_dp)-quantile_linear(x,0.25_dp)
    scale = min(s,iqr/1.34_dp)
    if (scale <= 0.0_dp) scale = max(s,1.0_dp)
    h = 0.9_dp*scale*real(n,dp)**(-0.2_dp)
    h = max(h,sqrt(epsilon(1.0_dp))*max(1.0_dp,abs(mean_value(x))))
  end function kde_bandwidth_nrd

  elemental real(dp) function kernel_cdf_scalar(u, kernel) result(v)
    real(dp), intent(in) :: u
    integer, intent(in) :: kernel
    select case(kernel)
    case(kernel_gaussian)
      v = normal_cdf(u)
    case(kernel_box)
      if (u <= -1.0_dp) then
        v = 0.0_dp
      else if (u >= 1.0_dp) then
        v = 1.0_dp
      else
        v = 0.5_dp*(u+1.0_dp)
      end if
    case(kernel_triangular)
      if (u <= -1.0_dp) then
        v = 0.0_dp
      else if (u < 0.0_dp) then
        v = 0.5_dp*(u+1.0_dp)**2
      else if (u < 1.0_dp) then
        v = 1.0_dp-0.5_dp*(1.0_dp-u)**2
      else
        v = 1.0_dp
      end if
    case(kernel_epanechnikov)
      if (u <= -1.0_dp) then
        v = 0.0_dp
      else if (u >= 1.0_dp) then
        v = 1.0_dp
      else
        v = 0.5_dp+0.75_dp*(u-u**3/3.0_dp)
      end if
    case default
      v = normal_cdf(u)
    end select
  end function kernel_cdf_scalar

  real(dp) function kde_cdf(x0, data, h, kernel) result(p)
    real(dp), intent(in) :: x0, data(:), h
    integer, intent(in) :: kernel
    p = sum(kernel_cdf_scalar((x0-data)/h,kernel))/real(size(data),dp)
  end function kde_cdf

  real(dp) function kernel_var(data, cl, kernel, bandwidth) result(value)
    real(dp), intent(in) :: data(:), cl
    integer, intent(in), optional :: kernel
    real(dp), intent(in), optional :: bandwidth
    integer :: k, iter
    real(dp) :: h, lo, hi, mid, s
    call validate_probability(cl)
    k = kernel_gaussian
    if (present(kernel)) k = kernel
    h = kde_bandwidth_nrd(data)
    if (present(bandwidth)) h = bandwidth
    s = max(sample_sd(data),h)
    lo = minval(data)-10.0_dp*s
    hi = maxval(data)+10.0_dp*s
    do iter = 1, 100
      mid = 0.5_dp*(lo+hi)
      if (kde_cdf(mid,data,h,k) < cl) then
        lo = mid
      else
        hi = mid
      end if
    end do
    value = 0.5_dp*(lo+hi)
  end function kernel_var

  real(dp) function kernel_es(data, cl, kernel, bandwidth, n_slices) result(value)
    real(dp), intent(in) :: data(:), cl
    integer, intent(in), optional :: kernel, n_slices
    real(dp), intent(in), optional :: bandwidth
    integer :: i, n, k
    real(dp) :: p, h
    call validate_probability(cl)
    k = kernel_gaussian
    if (present(kernel)) k = kernel
    h = kde_bandwidth_nrd(data)
    if (present(bandwidth)) h = bandwidth
    n = 500
    if (present(n_slices)) n = max(100,n_slices)
    value = 0.0_dp
    do i = 1, n
      p = cl+(1.0_dp-cl)*(real(i,dp)-0.5_dp)/real(n,dp)
      value = value+kernel_var(data,p,k,h)
    end do
    value = value/real(n,dp)
  end function kernel_es

  subroutine bootstrap_var_es(profit_loss, number_resamples, cl, mean_var, mean_es)
    real(dp), intent(in) :: profit_loss(:), cl
    integer, intent(in) :: number_resamples
    real(dp), intent(out) :: mean_var, mean_es
    real(dp), allocatable :: sample(:)
    real(dp) :: u
    integer :: b, i, n, idx
    n = size(profit_loss)
    if (number_resamples <= 0 .or. n == 0) error stop "bootstrap_var_es: invalid size"
    allocate(sample(n))
    mean_var = 0.0_dp
    mean_es = 0.0_dp
    do b = 1, number_resamples
      do i = 1, n
        call random_number(u)
        idx = min(n,1+int(u*real(n,dp)))
        sample(i) = profit_loss(idx)
      end do
      mean_var = mean_var+historical_var(sample,cl)
      mean_es = mean_es+historical_es(sample,cl,max(200,n))
    end do
    mean_var = mean_var/real(number_resamples,dp)
    mean_es = mean_es/real(number_resamples,dp)
  end subroutine bootstrap_var_es

  subroutine bootstrap_confidence_interval(profit_loss, number_resamples, cl, &
      lower_probability, upper_probability, var_interval, es_interval)
    real(dp), intent(in) :: profit_loss(:), cl, lower_probability, upper_probability
    integer, intent(in) :: number_resamples
    real(dp), intent(out) :: var_interval(2), es_interval(2)
    real(dp), allocatable :: sample(:), vars(:), ess(:)
    real(dp) :: u
    integer :: b, i, n, idx
    n = size(profit_loss)
    allocate(sample(n),vars(number_resamples),ess(number_resamples))
    do b = 1, number_resamples
      do i = 1, n
        call random_number(u)
        idx = min(n,1+int(u*real(n,dp)))
        sample(i) = profit_loss(idx)
      end do
      vars(b) = historical_var(sample,cl)
      ess(b) = historical_es(sample,cl,max(200,n))
    end do
    var_interval = [quantile_linear(vars,lower_probability),quantile_linear(vars,upper_probability)]
    es_interval = [quantile_linear(ess,lower_probability),quantile_linear(ess,upper_probability)]
  end subroutine bootstrap_confidence_interval

  real(dp) function boxcox_lambda(positive_data) result(best_lambda)
    real(dp), intent(in) :: positive_data(:)
    real(dp), allocatable :: y(:)
    real(dp) :: lambda, objective, best_objective, variance
    integer :: j, n
    if (minval(positive_data) <= 0.0_dp) error stop "boxcox_lambda: data must be positive"
    n = size(positive_data)
    allocate(y(n))
    best_objective = -huge(1.0_dp)
    best_lambda = 1.0_dp
    do j = 0, 400
      lambda = -2.0_dp+0.01_dp*real(j,dp)
      if (abs(lambda) < 1.0e-10_dp) then
        y = log(positive_data)
      else
        y = (positive_data**lambda-1.0_dp)/lambda
      end if
      variance = sum((y-mean_value(y))**2)/real(n,dp)
      objective = -0.5_dp*real(n,dp)*log(max(variance,tiny(1.0_dp))) + &
                  (lambda-1.0_dp)*sum(log(positive_data))
      if (objective > best_objective) then
        best_objective = objective
        best_lambda = lambda
      end if
    end do
  end function boxcox_lambda

  real(dp) function boxcox_var(profit_loss, cl, lambda_out) result(value)
    real(dp), intent(in) :: profit_loss(:), cl
    real(dp), intent(out), optional :: lambda_out
    real(dp), allocatable :: losses(:), positive(:), transformed(:)
    real(dp) :: shift, lambda, mu, sigma, tq
    call validate_probability(cl)
    losses = -profit_loss
    shift = 1.0_dp-minval(losses)
    positive = losses+shift
    lambda = boxcox_lambda(positive)
    if (present(lambda_out)) lambda_out = lambda
    if (abs(lambda) < 1.0e-10_dp) then
      transformed = log(positive)
    else
      transformed = (positive**lambda-1.0_dp)/lambda
    end if
    mu = mean_value(transformed)
    sigma = sample_sd(transformed)
    tq = mu+sigma*normal_quantile(cl)
    if (abs(lambda) < 1.0e-10_dp) then
      value = exp(tq)-shift
    else
      value = max(0.0_dp,1.0_dp+lambda*tq)**(1.0_dp/lambda)-shift
    end if
  end function boxcox_var

  real(dp) function boxcox_es(profit_loss, cl) result(value)
    real(dp), intent(in) :: profit_loss(:), cl
    real(dp), allocatable :: losses(:)
    real(dp) :: threshold
    integer :: count_tail
    losses = -profit_loss
    threshold = boxcox_var(profit_loss,cl)
    count_tail = count(losses >= threshold)
    if (count_tail == 0) then
      value = threshold
    else
      value = sum(losses,mask=losses>=threshold)/real(count_tail,dp)
    end if
  end function boxcox_es

  real(dp) function spectral_risk_normal(mu, sigma, gamma, number_slices) result(value)
    real(dp), intent(in) :: mu, sigma, gamma
    integer, intent(in) :: number_slices
    real(dp) :: p, weight, total
    integer :: i, n
    if (sigma < 0.0_dp .or. gamma <= 0.0_dp) error stop "spectral_risk_normal: invalid parameters"
    n = max(100,number_slices)
    total = 0.0_dp
    do i = 1, n
      p = (real(i,dp)-0.5_dp)/real(n,dp)
      weight = exp(-(1.0_dp-p)/gamma)/(gamma*(1.0_dp-exp(-1.0_dp/gamma)))
      total = total+(mu+sigma*normal_quantile(p))*weight
    end do
    value = total/real(n,dp)
  end function spectral_risk_normal

  pure real(dp) function quantile_standard_error_normal(prob, n, mu, sigma, bin_size) result(value)
    real(dp), intent(in) :: prob, mu, sigma, bin_size
    integer, intent(in) :: n
    real(dp) :: x, freq
    if (prob <= 0.0_dp .or. prob >= 1.0_dp .or. n <= 0 .or. sigma <= 0.0_dp .or. bin_size <= 0.0_dp) &
      error stop "quantile_standard_error_normal: invalid parameters"
    x = mu+sigma*normal_quantile(prob)
    freq = normal_cdf((x+0.5_dp*bin_size-mu)/sigma)-normal_cdf((x-0.5_dp*bin_size-mu)/sigma)
    value = sqrt(prob*(1.0_dp-prob)/(real(n,dp)*freq*freq))
  end function quantile_standard_error_normal

  real(dp) function quantile_standard_error_t(prob, n, mu, sigma, nu, bin_size) result(value)
    real(dp), intent(in) :: prob, mu, sigma, nu, bin_size
    integer, intent(in) :: n
    real(dp) :: x, lo, hi, freq
    if (prob <= 0.0_dp .or. prob >= 1.0_dp .or. n <= 0 .or. sigma <= 0.0_dp .or. &
        nu <= 0.0_dp .or. bin_size <= 0.0_dp) error stop "quantile_standard_error_t: invalid parameters"
    x = mu+sigma*student_t_quantile(prob,nu)
    lo = (x-0.5_dp*bin_size-mu)/sigma
    hi = (x+0.5_dp*bin_size-mu)/sigma
    freq = t_cdf_local(hi,nu)-t_cdf_local(lo,nu)
    value = sqrt(prob*(1.0_dp-prob)/(real(n,dp)*freq*freq))
  contains
    real(dp) function t_cdf_local(z,df) result(p)
      use dowd_math, only: student_t_cdf
      real(dp), intent(in) :: z,df
      p = student_t_cdf(z,df)
    end function t_cdf_local
  end function quantile_standard_error_t

end module dowd_risk
