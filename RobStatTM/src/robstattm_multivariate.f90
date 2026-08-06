! SPDX-License-Identifier: GPL-3.0-or-later
module robstattm_multivariate
  use robstattm_kinds, only : dp
  use robstattm_types, only : location_scale_result, covariance_result, projection_result, &
    robstattm_success, robstattm_invalid_argument, robstattm_singular, robstattm_no_convergence
  use robstattm_psi, only : rho_weight, rho_prime, rho_second, scale_m, rho_shr, weight_shr, &
    tuning_for_efficiency
  use robstattm_utils, only : covariance_to_correlation, weighted_center_covariance, &
    normalize_determinant, squared_mahalanobis, matrix_inverse, mean_value
  use rrcov_types, only : rrcov_covariance_result => covariance_result
  use rrcov_robust, only : cov_classic_rrcov => cov_classic, cov_mve
  use rrcov_stats, only : median, mad_scale, covariance_matrix, chi_square_quantile
  use rrcov_sort, only : order_smallest
  use rrcov_random, only : random_unit_vector, seed_random
  implicit none
  private
  public :: loc_scale_m, fast_mve, init_pp, kurt_sd_new
  public :: cov_classic, cov_rob, cov_rob_mm, cov_rob_rocke
  public :: rocke_weight, rocke_rho, rocke_scale, shr_scale
contains
  subroutine loc_scale_m(x, result, family, efficiency, max_iter, tol)
    real(dp), intent(in) :: x(:)
    type(location_scale_result), intent(out) :: result
    character(len=*), intent(in), optional :: family
    real(dp), intent(in), optional :: efficiency, tol
    integer, intent(in), optional :: max_iter
    character(len=16) :: fam
    real(dp) :: eff, tolerance, tuning, old_location, new_location, initial_scale
    real(dp), allocatable :: residuals(:), weights(:), psi(:), psip(:)
    integer :: i, mi, iter
    fam = 'mopt'
    if (present(family)) fam = family
    eff = 0.95_dp
    if (present(efficiency)) eff = efficiency
    mi = 50
    if (present(max_iter)) mi = max_iter
    tolerance = 1.0e-4_dp
    if (present(tol)) tolerance = tol
    if (size(x) == 0) then
      result%status = robstattm_invalid_argument
      return
    end if
    tuning = tuning_for_efficiency(eff, fam)
    old_location = median(x)
    initial_scale = mad_scale(x)
    if (initial_scale <= 1.0e-10_dp) then
      result%location = old_location
      result%scale = 0.0_dp
      result%standard_error = 0.0_dp
      result%converged = .true.
      return
    end if
    allocate(residuals(size(x)), weights(size(x)), psi(size(x)), psip(size(x)))
    do iter = 1, mi
      result%iterations = iter
      residuals = (x - old_location) / initial_scale
      do i = 1, size(x)
        weights(i) = rho_weight(residuals(i), fam, tuning)
      end do
      if (sum(weights) <= tiny(1.0_dp)) exit
      new_location = sum(weights * x) / sum(weights)
      if (abs(new_location - old_location) / initial_scale <= tolerance) then
        old_location = new_location
        result%converged = .true.
        exit
      end if
      old_location = new_location
    end do
    residuals = (x - old_location) / initial_scale
    do i = 1, size(x)
      psi(i) = rho_prime(residuals(i), fam, tuning)
      psip(i) = rho_second(residuals(i), fam, tuning)
    end do
    result%location = old_location
    result%standard_error = initial_scale * sqrt(mean_value(psi * psi) / &
      real(size(x),dp)) / max(abs(mean_value(psip)), sqrt(tiny(1.0_dp)))
    result%scale = scale_m(x - old_location, 0.5_dp, fam, tuning)
    result%status = merge(robstattm_success, robstattm_no_convergence, result%converged)
  end subroutine loc_scale_m

  subroutine fast_mve(x, result, nsamp, seed)
    real(dp), intent(in) :: x(:, :)
    type(covariance_result), intent(out) :: result
    integer, intent(in), optional :: nsamp, seed
    type(rrcov_covariance_result) :: fit
    integer :: trials, random_seed
    trials = 500
    if (present(nsamp)) trials = nsamp
    random_seed = 12345
    if (present(seed)) random_seed = seed
    call cov_mve(x, fit, alpha=0.5_dp, nsamp=trials, seed=random_seed, reweight=.false.)
    call copy_rrcov_covariance(fit, result)
    result%method = 'fast MVE'
    result%scale = median(fit%distances)
    result%converged = fit%status == 0
  end subroutine fast_mve

  subroutine init_pp(x, result, random_multiplier, fixed_multiplier, minimum_directions, seed)
    real(dp), intent(in) :: x(:, :)
    type(projection_result), intent(out) :: result
    integer, intent(in), optional :: random_multiplier, fixed_multiplier, minimum_directions, seed
    real(dp), allocatable :: center0(:), scales(:), z(:, :), direction(:), projections(:)
    real(dp), allocatable :: location(:), covariance(:, :), selected(:, :)
    integer, allocatable :: subset(:)
    integer :: n, p, nr, nf, nd, i, j, h, status, random_seed
    real(dp) :: loc, sc, cutoff
    n = size(x,1)
    p = size(x,2)
    if (n <= p .or. p < 1) then
      result%status = robstattm_invalid_argument
      return
    end if
    nr = 20
    if (present(random_multiplier)) nr = max(1, random_multiplier)
    nf = 10
    if (present(fixed_multiplier)) nf = max(1, fixed_multiplier)
    nd = 1000
    if (present(minimum_directions)) nd = max(1, minimum_directions)
    nd = max(nd, nr * p)
    random_seed = 24681357
    if (present(seed)) random_seed = seed
    call seed_random(random_seed)
    allocate(center0(p), scales(p), z(n,p), direction(p), projections(n), result%outlyingness(n))
    result%outlyingness = 0.0_dp
    do j = 1, p
      center0(j) = median(x(:,j))
      scales(j) = max(mad_scale(x(:,j)), 1.0e-12_dp)
      z(:,j) = (x(:,j)-center0(j))/scales(j)
    end do
    do i = 1, min(n, 2*nf*p)
      direction = z(i,:)
      if (sqrt(sum(direction*direction)) <= tiny(1.0_dp)) cycle
      direction = direction / sqrt(sum(direction*direction))
      projections = matmul(z, direction)
      loc = median(projections)
      sc = max(mad_scale(projections), 1.0e-12_dp)
      result%outlyingness = max(result%outlyingness, abs(projections-loc)/sc)
    end do
    do i = 1, nd
      call random_unit_vector(direction)
      projections = matmul(z, direction)
      loc = median(projections)
      sc = max(mad_scale(projections), 1.0e-12_dp)
      result%outlyingness = max(result%outlyingness, abs(projections-loc)/sc)
    end do
    h = max(p+1, min(n, (n+p+1)/2))
    allocate(subset(h), selected(h,p))
    call order_smallest(result%outlyingness, h, subset)
    selected = x(subset,:)
    allocate(location(p))
    covariance = covariance_matrix(selected, unbiased=.true., center=location, status=status)
    if (status /= 0) then
      result%status = robstattm_singular
      return
    end if
    call squared_mahalanobis(x, location, covariance, result%distances, status)
    cutoff = chi_square_quantile(0.975_dp, real(p,dp))
    allocate(result%outliers(n), result%center(p), result%covariance(p,p))
    result%outliers = result%distances > cutoff
    result%center = location
    result%covariance = covariance
    result%directions = nd + min(n,2*nf*p)
    result%status = merge(robstattm_success, robstattm_singular, status == 0)
  end subroutine init_pp

  subroutine kurt_sd_new(x, result, random_multiplier, fixed_multiplier, minimum_directions, seed)
    real(dp), intent(in) :: x(:, :)
    type(projection_result), intent(out) :: result
    integer, intent(in), optional :: random_multiplier, fixed_multiplier, minimum_directions, seed
    call init_pp(x, result, random_multiplier, fixed_multiplier, minimum_directions, seed)
  end subroutine kurt_sd_new

  subroutine cov_classic(x, result, unbiased, correlation)
    real(dp), intent(in) :: x(:, :)
    type(covariance_result), intent(out) :: result
    logical, intent(in), optional :: unbiased, correlation
    type(rrcov_covariance_result) :: fit
    logical :: ub, corr
    ub = .true.
    if (present(unbiased)) ub = unbiased
    corr = .false.
    if (present(correlation)) corr = correlation
    call cov_classic_rrcov(x, fit, unbiased=ub)
    call copy_rrcov_covariance(fit, result)
    if (corr) call covariance_to_correlation(result%covariance, result%correlation)
    result%method = 'classical'
    result%converged = fit%status == 0
  end subroutine cov_classic

  subroutine cov_rob(x, result, method, max_iter, tol, correlation)
    real(dp), intent(in) :: x(:, :)
    type(covariance_result), intent(out) :: result
    character(len=*), intent(in), optional :: method
    integer, intent(in), optional :: max_iter
    real(dp), intent(in), optional :: tol
    logical, intent(in), optional :: correlation
    character(len=16) :: name
    name = 'auto'
    if (present(method)) name = method
    if (trim(name) == 'auto') then
      if (size(x,2) >= 10) then
        name = 'Rocke'
      else
        name = 'MM'
      end if
    end if
    if (trim(name) == 'Rocke' .or. trim(name) == 'rocke') then
      call cov_rob_rocke(x, result, max_iter=max_iter, tol=tol, correlation=correlation)
    else
      call cov_rob_mm(x, result, max_iter=max_iter, tol=tol, correlation=correlation)
    end if
  end subroutine cov_rob

  subroutine cov_rob_mm(x, result, max_iter, tol, correlation, seed)
    real(dp), intent(in) :: x(:, :)
    type(covariance_result), intent(out) :: result
    integer, intent(in), optional :: max_iter, seed
    real(dp), intent(in), optional :: tol
    logical, intent(in), optional :: correlation
    type(projection_result) :: initial
    real(dp), allocatable :: center(:), old_center(:), covariance(:, :), old_covariance(:, :)
    real(dp), allocatable :: distances(:), weights(:), inverse_old(:, :), difference(:), ratio(:, :)
    real(dp) :: delta, constant, sigma, tolerance, difference1, difference2, parameter_change, factor
    integer :: n, p, mi, i, status, random_seed, iter
    logical :: corr
    n = size(x,1)
    p = size(x,2)
    if (n <= p .or. p < 1) then
      result%status = robstattm_invalid_argument
      return
    end if
    mi = 50
    if (present(max_iter)) mi = max_iter
    tolerance = 1.0e-4_dp
    if (present(tol)) tolerance = tol
    corr = .false.
    if (present(correlation)) corr = correlation
    random_seed = 12345
    if (present(seed)) random_seed = seed
    call init_pp(x, initial, minimum_directions=max(200,20*p), seed=random_seed)
    if (initial%status /= 0) then
      result%status = initial%status
      return
    end if
    center = initial%center
    covariance = initial%covariance
    call normalize_determinant(covariance, status)
    call squared_mahalanobis(x, center, covariance, distances, status)
    delta = 0.5_dp * (1.0_dp - real(p,dp)/real(n,dp))
    constant = 4.5041_dp/real(p,dp) - 1.1117_dp*real(p,dp)/real(n,dp) + 0.61161_dp
    sigma = shr_scale(distances, delta) * max(constant, 0.1_dp)
    old_center = 0.0_dp * center
    allocate(old_covariance(p,p))
    old_covariance = 0.0_dp
    do i = 1, p
      old_covariance(i,i) = 1.0_dp
    end do
    allocate(weights(n), difference(p))
    parameter_change = huge(1.0_dp)
    do iter = 1, mi
      result%iterations = iter
      weights = weight_shr(distances/max(sigma,tiny(1.0_dp)))
      call weighted_center_covariance(x, weights, center, covariance)
      call normalize_determinant(covariance, status)
      if (status /= 0) exit
      call squared_mahalanobis(x, center, covariance, distances, status)
      difference = center-old_center
      inverse_old = matrix_inverse(old_covariance, status)
      if (status /= 0) exit
      difference1 = dot_product(difference,matmul(inverse_old,difference))
      ratio = matmul(inverse_old,covariance)
      do i = 1, p
        ratio(i,i) = ratio(i,i)-1.0_dp
      end do
      difference2 = maxval(abs(ratio))
      parameter_change = max(difference1,difference2)
      old_center = center
      old_covariance = covariance
      if (parameter_change <= tolerance) exit
    end do
    call consistency_scale(distances, p, factor)
    covariance = covariance * factor
    distances = distances / factor
    result%center = center
    result%covariance = covariance
    result%distances = distances
    result%weights = weights
    result%scale = sigma
    result%objective = shr_objective(distances/max(sigma,tiny(1.0_dp)))
    result%converged = parameter_change <= tolerance
    result%status = merge(robstattm_success, robstattm_no_convergence, result%converged)
    result%method = 'MM-SHR'
    if (corr) call covariance_to_correlation(covariance,result%correlation)
  end subroutine cov_rob_mm

  subroutine cov_rob_rocke(x, result, initial_method, max_steps, proportion_minimum, q, max_iter, tol, correlation, seed)
    real(dp), intent(in) :: x(:, :)
    type(covariance_result), intent(out) :: result
    character(len=*), intent(in), optional :: initial_method
    integer, intent(in), optional :: max_steps, proportion_minimum, q, max_iter, seed
    real(dp), intent(in), optional :: tol
    logical, intent(in), optional :: correlation
    type(projection_result) :: projection
    type(covariance_result) :: mve
    real(dp), allocatable :: center(:), old_center(:), covariance(:, :), old_covariance(:, :)
    real(dp), allocatable :: distances(:), weights(:), difference(:), inverse_old(:, :), ratio(:, :), sorted_difference(:)
    real(dp) :: gamma0, gamma, alpha, delta, sigma, old_sigma, tolerance
    real(dp) :: parameter_change, scale_change, difference1, difference2, fraction, factor
    integer :: n, p, mi, steps, propmin, exponent, status, i, line_step, random_seed, target, iter
    logical :: corr
    character(len=8) :: start
    n = size(x,1)
    p = size(x,2)
    if (n <= p .or. p < 1) then
      result%status = robstattm_invalid_argument
      return
    end if
    start = 'K'
    if (present(initial_method)) start = initial_method
    steps = 5
    if (present(max_steps)) steps = max_steps
    propmin = 2
    if (present(proportion_minimum)) propmin = proportion_minimum
    exponent = 2
    if (present(q)) exponent = q
    mi = 50
    if (present(max_iter)) mi = max_iter
    tolerance = 1.0e-4_dp
    if (present(tol)) tolerance = tol
    corr = .false.
    if (present(correlation)) corr = correlation
    random_seed = 12345
    if (present(seed)) random_seed = seed
    if (trim(start) == 'mve' .or. trim(start) == 'MVE') then
      call fast_mve(x,mve,nsamp=500,seed=random_seed)
      center = mve%center
      covariance = mve%covariance
    else
      call init_pp(x,projection,minimum_directions=max(500,20*p),seed=random_seed)
      center = projection%center
      covariance = projection%covariance
    end if
    call normalize_determinant(covariance,status)
    call squared_mahalanobis(x,center,covariance,distances,status)
    call rocke_constant(p,n,start,gamma0,alpha)
    delta = 0.5_dp*(1.0_dp-real(p,dp)/real(n,dp))
    sigma = rocke_scale(distances,gamma0,exponent,delta)
    allocate(sorted_difference(n))
    sorted_difference = abs(distances/max(sigma,tiny(1.0_dp))-1.0_dp)
    call sort_in_place(sorted_difference)
    target = min(n,max(1,propmin*p))
    gamma = max(gamma0,sorted_difference(target))
    old_sigma = rocke_scale(distances,gamma,exponent,delta)
    old_center = center
    old_covariance = covariance
    allocate(weights(n),difference(p))
    parameter_change = huge(1.0_dp)
    scale_change = huge(1.0_dp)
    do iter = 1, mi
      result%iterations = iter
      weights = rocke_weight(distances/max(old_sigma,tiny(1.0_dp)),gamma,exponent)
      call weighted_center_covariance(x,weights,center,covariance,real(n,dp))
      call normalize_determinant(covariance,status)
      if(status/=0)exit
      call squared_mahalanobis(x,center,covariance,distances,status)
      sigma = rocke_scale(distances,gamma,exponent,delta)
      line_step=0
      fraction=1.0_dp
      do while(sigma>old_sigma .and. line_step<steps)
        fraction=0.5_dp*fraction
        center=fraction*center+(1.0_dp-fraction)*old_center
        covariance=fraction*covariance+(1.0_dp-fraction)*old_covariance
        call normalize_determinant(covariance,status)
        call squared_mahalanobis(x,center,covariance,distances,status)
        sigma=rocke_scale(distances,gamma,exponent,delta)
        line_step=line_step+1
      end do
      difference=center-old_center
      inverse_old=matrix_inverse(old_covariance,status)
      if(status/=0)exit
      difference1=dot_product(difference,matmul(inverse_old,difference))/real(p,dp)
      ratio=matmul(inverse_old,covariance)
      do i=1,p
        ratio(i,i)=ratio(i,i)-1.0_dp
      end do
      difference2=maxval(abs(ratio))
      parameter_change=max(difference1,difference2)
      scale_change=1.0_dp-sigma/max(old_sigma,tiny(1.0_dp))
      old_center=center
      old_covariance=covariance
      old_sigma=sigma
      if(parameter_change<=tolerance .and. abs(scale_change)<=tolerance)exit
      if(scale_change<0.0_dp)exit
    end do
    call consistency_scale(distances,p,factor)
    covariance=covariance*factor
    distances=distances/factor
    result%center=center
    result%covariance=covariance
    result%distances=distances
    result%weights=weights
    result%scale=sigma
    result%gamma=gamma
    result%objective=sigma
    result%converged=parameter_change<=tolerance .and. abs(scale_change)<=tolerance
    result%status=merge(robstattm_success,robstattm_no_convergence,result%converged)
    result%method='Rocke S'
    if(corr)call covariance_to_correlation(covariance,result%correlation)
  end subroutine cov_rob_rocke

  elemental function rocke_weight(t, gamma, q) result(value)
    real(dp), intent(in) :: t, gamma
    integer, intent(in) :: q
    real(dp) :: value, s
    s=(t-1.0_dp)/max(gamma,tiny(1.0_dp))
    if(abs(s)>1.0_dp)then
      value=0.0_dp
    else
      value=max(0.0_dp,1.0_dp-s**q)
    end if
  end function rocke_weight

  elemental function rocke_rho(t,gamma,q) result(value)
    real(dp),intent(in)::t,gamma
    integer,intent(in)::q
    real(dp)::value,u
    u=(t-1.0_dp)/max(gamma,tiny(1.0_dp))
    if(u>=1.0_dp)then
      value=1.0_dp
    else if(u< -1.0_dp)then
      value=0.0_dp
    else
      value=u/(2.0_dp*real(q,dp))*(real(q+1,dp)-u**q)+0.5_dp
    end if
  end function rocke_rho

  function rocke_scale(x,gamma,q,delta,tol) result(scale)
    real(dp),intent(in)::x(:),gamma,delta
    integer,intent(in)::q
    real(dp),intent(in),optional::tol
    real(dp)::scale,lo,hi,mid,value,tolerance
    integer::i
    if(real(count(abs(x)<=1.0e-16_dp),dp)>real(size(x),dp)*(1.0_dp-delta))then
      scale=0.0_dp
      return
    end if
    tolerance=1.0e-6_dp
    if(present(tol))tolerance=tol
    lo=max(minval(pack(abs(x),abs(x)>0.0_dp))/(1.0_dp+gamma),tiny(1.0_dp))
    hi=maxval(abs(x))*10.0_dp+1.0_dp
    do i=1,120
      mid=0.5_dp*(lo+hi)
      value=mean_rocke_rho(x/mid,gamma,q)-delta
      if(value>0.0_dp)then
        lo=mid
      else
        hi=mid
      end if
      if(abs(hi-lo)<=tolerance*(1.0_dp+mid))exit
    end do
    scale=0.5_dp*(lo+hi)
  end function rocke_scale

  function shr_scale(distances,delta,initial,max_iter,tol) result(scale)
    real(dp),intent(in)::distances(:),delta
    real(dp),intent(in),optional::initial,tol
    integer,intent(in),optional::max_iter
    real(dp)::scale,s1,s2,s3,y1,y2,den,tolerance,median_distance
    integer::i,mi
    if(real(count(distances<1.0e-16_dp),dp)/real(size(distances),dp)>=0.5_dp)then
      scale=0.0_dp
      return
    end if
    median_distance=median(distances)
    s1=median_distance
    if(present(initial))s1=initial
    tolerance=1.0e-4_dp
    if(present(tol))tolerance=tol
    mi=50
    if(present(max_iter))mi=max_iter
    y1=shr_fixed_point(s1,distances,delta)
    do while(y1>s1 .and. s1<1000.0_dp*max(median_distance,tiny(1.0_dp)))
      s1=2.0_dp*s1
      y1=shr_fixed_point(s1,distances,delta)
    end do
    s2=y1
    do i=1,mi
      y2=shr_fixed_point(s2,distances,delta)
      den=s2-y2+y1-s1
      if(abs(den)<tolerance*max(s1,tiny(1.0_dp)))then
        s3=s2
        exit
      end if
      s3=(y1*s2-s1*y2)/den
      s1=s2
      s2=max(s3,tiny(1.0_dp))
      y1=y2
    end do
    scale=s2
  end function shr_scale

  function shr_fixed_point(scale,distances,delta) result(value)
    real(dp),intent(in)::scale,distances(:),delta
    real(dp)::value
    integer::i
    value=0.0_dp
    do i=1,size(distances)
      value=value+rho_shr(distances(i)/max(scale,tiny(1.0_dp)))
    end do
    value=scale*value/(real(size(distances),dp)*delta)
  end function shr_fixed_point

  function shr_objective(d) result(value)
    real(dp),intent(in)::d(:)
    real(dp)::value
    integer::i
    value=0.0_dp
    do i=1,size(d)
      value=value+rho_shr(d(i))
    end do
  end function shr_objective

  function mean_rocke_rho(x,gamma,q) result(value)
    real(dp),intent(in)::x(:),gamma
    integer,intent(in)::q
    real(dp)::value
    integer::i
    value=0.0_dp
    do i=1,size(x)
      value=value+rocke_rho(x(i),gamma,q)
    end do
    value=value/real(size(x),dp)
  end function mean_rocke_rho

  subroutine rocke_constant(p,n,initial,gamma,alpha)
    integer,intent(in)::p,n
    character(len=*),intent(in)::initial
    real(dp),intent(out)::gamma,alpha
    real(dp)::beta(3)
    if(trim(initial)=='mve'.or.trim(initial)=='MVE')then
      beta=[-5.4358_dp,-0.50303_dp,0.4214_dp]
    else
      beta=[-6.1357_dp,-1.0078_dp,0.81564_dp]
    end if
    if(p>=15)then
      alpha=exp(beta(1)+beta(2)*log(real(p,dp))+beta(3)*log(real(n,dp)))
      gamma=min(chi_square_quantile(1.0_dp-alpha,real(p,dp))/real(p,dp)-1.0_dp,1.0_dp)
    else
      gamma=1.0_dp
      alpha=1.0e-6_dp
    end if
  end subroutine rocke_constant

  subroutine consistency_scale(distances,p,factor)
    real(dp),intent(in)::distances(:)
    integer,intent(in)::p
    real(dp),intent(out)::factor
    factor=median(distances)/max(chi_square_quantile(0.5_dp,real(p,dp)),tiny(1.0_dp))
    factor=max(factor,tiny(1.0_dp))
  end subroutine consistency_scale

  subroutine copy_rrcov_covariance(source,target)
    type(rrcov_covariance_result),intent(in)::source
    type(covariance_result),intent(out)::target
    target%center=source%center
    target%covariance=source%covariance
    target%distances=source%distances
    target%weights=source%weights
    target%subset=source%subset
    target%objective=source%objective
    target%iterations=source%iterations
    target%status=source%status
    target%converged=source%status==0
    target%method=source%method
  end subroutine copy_rrcov_covariance

  subroutine sort_in_place(x)
    real(dp),intent(inout)::x(:)
    real(dp)::temp
    integer::i,j
    do i=2,size(x)
      temp=x(i)
      j=i-1
      do while(j>=1)
        if(x(j)<=temp)exit
        x(j+1)=x(j)
        j=j-1
      end do
      x(j+1)=temp
    end do
  end subroutine sort_in_place
end module robstattm_multivariate
