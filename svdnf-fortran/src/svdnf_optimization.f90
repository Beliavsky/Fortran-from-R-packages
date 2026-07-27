! SPDX-License-Identifier: GPL-3.0-only
module svdnf_optimization
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use svdnf_kinds, only : dp
  use svdnf_types, only : svm_dynamics, grid_type, filter_result, optimization_result, model_custom, &
    model_heston, model_bates, model_duffie_pan_singleton, model_taylor, &
    model_taylor_leverage, model_pitt_malik_doucet, model_capm_sv
  use svdnf_models, only : parameter_vector, parameter_count, set_parameter_vector, validate_dynamics
  use svdnf_filter, only : dnf_filter
  use svdnf_stats, only : mean_value, standard_deviation
  implicit none
  private
  public :: dnf_optimize, dnf_optimize_custom, dnf_optim, initial_guess, numerical_hessian
  public :: custom_parameter_setter

  abstract interface
    function objective_function(x) result(value)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp) :: value
    end function objective_function

    subroutine custom_parameter_setter(dynamics, parameters, ok)
      import dp, svm_dynamics
      type(svm_dynamics), intent(inout) :: dynamics
      real(dp), intent(in) :: parameters(:)
      logical, intent(out) :: ok
    end subroutine custom_parameter_setter
  end interface

  type(svm_dynamics), save :: context_dynamics
  type(grid_type), save :: context_grid
  real(dp), allocatable, save :: context_data(:), context_factors(:,:)
  logical, save :: context_has_grid = .false., context_has_factors = .false.
  integer, save :: context_n = 50, context_k = 20, context_r = 1
  logical, save :: context_is_custom = .false.
  procedure(custom_parameter_setter), pointer, save :: context_setter => null()

contains

  function dnf_optimize(dynamics, data, initial_parameters, factors, grids, n, k, r, &
      max_iterations, tolerance, calculate_hessian) result(output)
    type(svm_dynamics), intent(in) :: dynamics
    real(dp), intent(in) :: data(:)
    real(dp), intent(in), optional :: initial_parameters(:)
    real(dp), intent(in), optional :: factors(:,:)
    type(grid_type), intent(in), optional :: grids
    integer, intent(in), optional :: n, k, r, max_iterations
    real(dp), intent(in), optional :: tolerance
    logical, intent(in), optional :: calculate_hessian
    type(optimization_result) :: output
    real(dp), allocatable :: start(:), optimum(:), inverse(:,:)
    real(dp) :: best_value, tol_value
    integer :: max_iter
    logical :: hessian_requested, inverse_ok, parameter_ok
    type(svm_dynamics) :: fitted_dynamics

    if (dynamics%model_id == model_custom) then
      output%message = 'Automatic parameter packing is not defined for custom callback models.'
      return
    end if
    if (present(initial_parameters)) then
      allocate(start(size(initial_parameters)))
      start = initial_parameters
    else
      allocate(start(parameter_count(dynamics)))
      start = initial_guess(dynamics,data,factors)
    end if
    if (size(start) == 0) then
      output%message = 'No optimizable parameters were found.'
      return
    end if

    call set_objective_context(dynamics,data,factors,grids,n,k,r)
    context_is_custom = .false.
    nullify(context_setter)
    max_iter = 250
    tol_value = 1.0e-6_dp
    if (present(max_iterations)) max_iter = max_iterations
    if (present(tolerance)) tol_value = tolerance
    call nelder_mead(dnf_objective,start,optimum,best_value,output%iterations,output%evaluations, &
      output%converged,max_iter,tol_value)
    output%parameters = optimum
    output%log_likelihood = -best_value

    hessian_requested = .true.
    if (present(calculate_hessian)) hessian_requested = calculate_hessian
    if (hessian_requested) then
      output%hessian = numerical_hessian(dnf_objective,optimum)
      call invert_matrix(output%hessian,inverse,inverse_ok)
      if (inverse_ok) then
        allocate(output%standard_errors(size(optimum)))
        output%standard_errors = sqrt(max(diagonal(inverse),0.0_dp))
      end if
    end if

    fitted_dynamics = dynamics
    call set_parameter_vector(fitted_dynamics,optimum,parameter_ok)
    if (.not. parameter_ok) then
      output%message = 'The optimized parameter vector could not be applied.'
      return
    end if
    if (present(factors) .and. present(grids)) then
      output%filter = dnf_filter(fitted_dynamics,data,factors=factors,grids=grids,n=n,k=k,r=r)
    else if (present(factors)) then
      output%filter = dnf_filter(fitted_dynamics,data,factors=factors,n=n,k=k,r=r)
    else if (present(grids)) then
      output%filter = dnf_filter(fitted_dynamics,data,grids=grids,n=n,k=k,r=r)
    else
      output%filter = dnf_filter(fitted_dynamics,data,n=n,k=k,r=r)
    end if
    output%ok = output%filter%ok
    if (output%ok) then
      if (output%converged) then
        output%message = 'Optimization converged.'
      else
        output%message = 'Optimization stopped at the iteration limit.'
      end if
    else
      output%message = output%filter%message
    end if
  end function dnf_optimize

  function dnf_optimize_custom(dynamics, data, initial_parameters, setter, factors, grids, n, k, r, &
      max_iterations, tolerance, calculate_hessian) result(output)
    type(svm_dynamics), intent(in) :: dynamics
    real(dp), intent(in) :: data(:), initial_parameters(:)
    procedure(custom_parameter_setter) :: setter
    real(dp), intent(in), optional :: factors(:,:)
    type(grid_type), intent(in), optional :: grids
    integer, intent(in), optional :: n, k, r, max_iterations
    real(dp), intent(in), optional :: tolerance
    logical, intent(in), optional :: calculate_hessian
    type(optimization_result) :: output
    real(dp), allocatable :: optimum(:), inverse(:,:)
    real(dp) :: best_value, tol_value
    integer :: max_iter
    logical :: hessian_requested, inverse_ok, parameter_ok
    type(svm_dynamics) :: fitted_dynamics

    if (dynamics%model_id /= model_custom) then
      output%message = 'dnf_optimize_custom requires a custom dynamics object.'
      return
    end if
    call set_objective_context(dynamics,data,factors,grids,n,k,r)
    context_is_custom = .true.
    context_setter => setter
    max_iter=250
    tol_value=1.0e-6_dp
    if (present(max_iterations)) max_iter=max_iterations
    if (present(tolerance)) tol_value=tolerance
    call nelder_mead(dnf_objective,initial_parameters,optimum,best_value,output%iterations, &
      output%evaluations,output%converged,max_iter,tol_value)
    output%parameters=optimum
    output%log_likelihood=-best_value
    hessian_requested=.true.
    if (present(calculate_hessian)) hessian_requested=calculate_hessian
    if (hessian_requested) then
      output%hessian=numerical_hessian(dnf_objective,optimum)
      call invert_matrix(output%hessian,inverse,inverse_ok)
      if (inverse_ok) then
        allocate(output%standard_errors(size(optimum)))
        output%standard_errors=sqrt(max(diagonal(inverse),0.0_dp))
      end if
    end if
    fitted_dynamics=dynamics
    call setter(fitted_dynamics,optimum,parameter_ok)
    if (.not. parameter_ok) then
      output%message='The custom parameter setter rejected the optimized vector.'
      return
    end if
    output%filter=filter_from_context(fitted_dynamics)
    output%ok=output%filter%ok
    if (output%ok) then
      if (output%converged) then
        output%message='Custom-model optimization converged.'
      else
        output%message='Custom-model optimization stopped at the iteration limit.'
      end if
    else
      output%message=output%filter%message
    end if
  end function dnf_optimize_custom

  function dnf_optim(dynamics, data, par, factors, grids, n, k, r, max_iterations, tolerance) result(output)
    type(svm_dynamics), intent(in) :: dynamics
    real(dp), intent(in) :: data(:)
    real(dp), intent(in), optional :: par(:), factors(:,:)
    type(grid_type), intent(in), optional :: grids
    integer, intent(in), optional :: n, k, r, max_iterations
    real(dp), intent(in), optional :: tolerance
    type(optimization_result) :: output
    if (present(par) .and. present(factors) .and. present(grids)) then
      output=dnf_optimize(dynamics,data,par,factors,grids,n,k,r,max_iterations,tolerance)
    else if (present(par) .and. present(factors)) then
      output=dnf_optimize(dynamics,data,initial_parameters=par,factors=factors,n=n,k=k,r=r, &
        max_iterations=max_iterations,tolerance=tolerance)
    else if (present(par) .and. present(grids)) then
      output=dnf_optimize(dynamics,data,initial_parameters=par,grids=grids,n=n,k=k,r=r, &
        max_iterations=max_iterations,tolerance=tolerance)
    else if (present(factors) .and. present(grids)) then
      output=dnf_optimize(dynamics,data,factors=factors,grids=grids,n=n,k=k,r=r, &
        max_iterations=max_iterations,tolerance=tolerance)
    else if (present(par)) then
      output=dnf_optimize(dynamics,data,initial_parameters=par,n=n,k=k,r=r, &
        max_iterations=max_iterations,tolerance=tolerance)
    else if (present(factors)) then
      output=dnf_optimize(dynamics,data,factors=factors,n=n,k=k,r=r, &
        max_iterations=max_iterations,tolerance=tolerance)
    else if (present(grids)) then
      output=dnf_optimize(dynamics,data,grids=grids,n=n,k=k,r=r, &
        max_iterations=max_iterations,tolerance=tolerance)
    else
      output=dnf_optimize(dynamics,data,n=n,k=k,r=r,max_iterations=max_iterations,tolerance=tolerance)
    end if
  end function dnf_optim

  subroutine set_objective_context(dynamics,data,factors,grids,n,k,r)
    type(svm_dynamics), intent(in) :: dynamics
    real(dp), intent(in) :: data(:)
    real(dp), intent(in), optional :: factors(:,:)
    type(grid_type), intent(in), optional :: grids
    integer, intent(in), optional :: n,k,r
    context_dynamics = dynamics
    context_data = data
    context_has_factors = present(factors)
    context_has_grid = present(grids)
    if (allocated(context_factors)) deallocate(context_factors)
    if (present(factors)) context_factors = factors
    if (present(grids)) context_grid = grids
    context_n=50; context_k=20; context_r=1
    if (present(n)) context_n=n
    if (present(k)) context_k=k
    if (present(r)) context_r=r
  end subroutine set_objective_context

  function dnf_objective(parameters) result(value)
    real(dp), intent(in) :: parameters(:)
    real(dp) :: value
    type(svm_dynamics) :: trial
    type(filter_result) :: filtered
    logical :: valid_parameters, valid_dynamics
    character(len=160) :: message
    trial=context_dynamics
    if (context_is_custom) then
      if (.not. associated(context_setter)) then
        value=1.0e150_dp
        return
      end if
      call context_setter(trial,parameters,valid_parameters)
    else
      call set_parameter_vector(trial,parameters,valid_parameters)
    end if
    if (.not. valid_parameters) then
      value=1.0e150_dp
      return
    end if
    call validate_dynamics(trial,valid_dynamics,message)
    if (.not. valid_dynamics) then
      value=1.0e150_dp
      return
    end if
    filtered=filter_from_context(trial)
    if (filtered%ok .and. ieee_is_finite(filtered%log_likelihood)) then
      value=-filtered%log_likelihood
    else
      value=1.0e150_dp
    end if
  end function dnf_objective

  function filter_from_context(dynamics) result(filtered)
    type(svm_dynamics), intent(in) :: dynamics
    type(filter_result) :: filtered
    if (context_has_factors .and. context_has_grid) then
      filtered=dnf_filter(dynamics,context_data,factors=context_factors,grids=context_grid)
    else if (context_has_factors) then
      filtered=dnf_filter(dynamics,context_data,factors=context_factors,n=context_n,k=context_k,r=context_r)
    else if (context_has_grid) then
      filtered=dnf_filter(dynamics,context_data,grids=context_grid)
    else
      filtered=dnf_filter(dynamics,context_data,n=context_n,k=context_k,r=context_r)
    end if
  end function filter_from_context

  function initial_guess(dynamics, data, factors) result(parameters)
    type(svm_dynamics), intent(in) :: dynamics
    real(dp), intent(in) :: data(:)
    real(dp), intent(in), optional :: factors(:,:)
    real(dp), allocatable :: parameters(:)
    type(svm_dynamics) :: guess
    real(dp), allocatable :: adjusted(:)
    real(dp) :: variance_value, mean_return, log_variance, correlation, xmean, ymean
    integer :: k
    guess = dynamics
    adjusted = data
    if (present(factors) .and. allocated(guess%coefs)) then
      if (size(factors,1)==size(data) .and. size(factors,2)==size(guess%coefs)) then
        adjusted = adjusted - matmul(factors,guess%coefs)
      end if
    end if
    mean_return = mean_value(adjusted)
    variance_value = max(standard_deviation(adjusted)**2,1.0e-6_dp)
    log_variance = log(variance_value)
    select case (guess%model_id)
    case (model_heston,model_bates,model_duffie_pan_singleton)
      guess%mu = mean_return/max(guess%h,1.0e-6_dp)
      guess%theta = max(variance_value/max(guess%h,1.0e-6_dp),1.0e-4_dp)
      guess%kappa = max(guess%kappa,1.0_dp)
      guess%sigma = max(guess%sigma,0.1_dp)
      guess%rho = min(max(guess%rho,-0.9_dp),0.9_dp)
      if (guess%model_id /= model_heston) then
        guess%omega=max(guess%omega,0.01_dp)
        guess%delta=max(guess%delta,0.005_dp)
      end if
      if (guess%model_id == model_duffie_pan_singleton) guess%nu=max(guess%nu,0.001_dp)
    case (model_taylor,model_taylor_leverage,model_pitt_malik_doucet,model_capm_sv)
      guess%theta = log_variance
      guess%phi = min(max(guess%phi,0.5_dp),0.98_dp)
      guess%sigma=max(guess%sigma,0.1_dp)
      if (size(adjusted)>2) then
        xmean=mean_value(adjusted(1:size(adjusted)-1)**2)
        ymean=mean_value(adjusted(2:size(adjusted))**2)
        correlation=sum((adjusted(1:size(adjusted)-1)**2-xmean)* &
          (adjusted(2:size(adjusted))**2-ymean))
        correlation=correlation/max(sqrt(sum((adjusted(1:size(adjusted)-1)**2-xmean)**2)* &
          sum((adjusted(2:size(adjusted))**2-ymean)**2)),tiny(1.0_dp))
        guess%phi=min(max(correlation,0.2_dp),0.98_dp)
      end if
      if (guess%model_id == model_pitt_malik_doucet) then
        guess%p=min(max(guess%p,0.001_dp),0.2_dp)
        guess%delta=max(guess%delta,0.01_dp)
      end if
      if (guess%model_id == model_capm_sv .and. present(factors)) then
        k=size(factors,2)
        if (.not. allocated(guess%coefs)) allocate(guess%coefs(k),source=0.0_dp)
      end if
    end select
    parameters = parameter_vector(guess)
  end function initial_guess

  function numerical_hessian(function_value, x) result(hessian)
    procedure(objective_function) :: function_value
    real(dp), intent(in) :: x(:)
    real(dp), allocatable :: hessian(:,:)
    real(dp), allocatable :: xp(:), xm(:), xpp(:), xpm(:), xmp(:), xmm(:), step(:)
    real(dp) :: f0
    integer :: i, j, n
    n=size(x)
    allocate(hessian(n,n),xp(n),xm(n),xpp(n),xpm(n),xmp(n),xmm(n),step(n))
    step = max(1.0e-4_dp*max(abs(x),1.0_dp),1.0e-6_dp)
    f0=function_value(x)
    do i=1,n
      xp=x; xm=x
      xp(i)=xp(i)+step(i); xm(i)=xm(i)-step(i)
      hessian(i,i)=(function_value(xp)-2.0_dp*f0+function_value(xm))/step(i)**2
      do j=i+1,n
        xpp=x; xpm=x; xmp=x; xmm=x
        xpp(i)=xpp(i)+step(i); xpp(j)=xpp(j)+step(j)
        xpm(i)=xpm(i)+step(i); xpm(j)=xpm(j)-step(j)
        xmp(i)=xmp(i)-step(i); xmp(j)=xmp(j)+step(j)
        xmm(i)=xmm(i)-step(i); xmm(j)=xmm(j)-step(j)
        hessian(i,j)=(function_value(xpp)-function_value(xpm)-function_value(xmp)+ &
          function_value(xmm))/(4.0_dp*step(i)*step(j))
        hessian(j,i)=hessian(i,j)
      end do
    end do
  end function numerical_hessian

  subroutine nelder_mead(function_value, start, optimum, best_value, iterations, evaluations, &
      converged, max_iterations, tolerance)
    procedure(objective_function) :: function_value
    real(dp), intent(in) :: start(:)
    real(dp), allocatable, intent(out) :: optimum(:)
    real(dp), intent(out) :: best_value
    integer, intent(out) :: iterations, evaluations
    logical, intent(out) :: converged
    integer, intent(in) :: max_iterations
    real(dp), intent(in) :: tolerance
    real(dp), allocatable :: simplex(:,:), values(:), centroid(:), reflected(:), expanded(:), contracted(:)
    real(dp) :: step, fr, fe, fc
    integer :: n, i, worst, second_worst
    n=size(start)
    allocate(simplex(n,n+1),values(n+1),centroid(n),reflected(n),expanded(n),contracted(n))
    simplex(:,1)=start
    do i=1,n
      simplex(:,i+1)=start
      step=0.05_dp*max(abs(start(i)),1.0_dp)
      simplex(i,i+1)=simplex(i,i+1)+step
    end do
    do i=1,n+1
      values(i)=function_value(simplex(:,i))
    end do
    evaluations=n+1
    converged=.false.
    do iterations=1,max_iterations
      call sort_simplex(simplex,values)
      if (maxval(abs(simplex(:,2:n+1)-spread(simplex(:,1),2,n)))/ &
          max(1.0_dp,maxval(abs(simplex(:,1)))) < tolerance .and. &
          maxval(abs(values(2:n+1)-values(1))) < tolerance) then
        converged=.true.
        exit
      end if
      worst=n+1; second_worst=n
      centroid=sum(simplex(:,1:n),dim=2)/real(n,dp)
      reflected=centroid+(centroid-simplex(:,worst))
      fr=function_value(reflected); evaluations=evaluations+1
      if (fr < values(1)) then
        expanded=centroid+2.0_dp*(reflected-centroid)
        fe=function_value(expanded); evaluations=evaluations+1
        if (fe < fr) then
          simplex(:,worst)=expanded; values(worst)=fe
        else
          simplex(:,worst)=reflected; values(worst)=fr
        end if
      else if (fr < values(second_worst)) then
        simplex(:,worst)=reflected; values(worst)=fr
      else
        if (fr < values(worst)) then
          contracted=centroid+0.5_dp*(reflected-centroid)
        else
          contracted=centroid+0.5_dp*(simplex(:,worst)-centroid)
        end if
        fc=function_value(contracted); evaluations=evaluations+1
        if (fc < min(fr,values(worst))) then
          simplex(:,worst)=contracted; values(worst)=fc
        else
          do i=2,n+1
            simplex(:,i)=simplex(:,1)+0.5_dp*(simplex(:,i)-simplex(:,1))
            values(i)=function_value(simplex(:,i))
          end do
          evaluations=evaluations+n
        end if
      end if
    end do
    call sort_simplex(simplex,values)
    optimum=simplex(:,1)
    best_value=values(1)
  end subroutine nelder_mead

  subroutine sort_simplex(simplex,values)
    real(dp), intent(inout) :: simplex(:,:), values(:)
    integer :: i,j
    real(dp) :: value
    real(dp), allocatable :: column(:)
    allocate(column(size(simplex,1)))
    do i=2,size(values)
      value=values(i); column=simplex(:,i); j=i-1
      do while (j>=1)
        if (values(j)<=value) exit
        values(j+1)=values(j); simplex(:,j+1)=simplex(:,j); j=j-1
      end do
      values(j+1)=value; simplex(:,j+1)=column
    end do
  end subroutine sort_simplex

  subroutine invert_matrix(a,inverse,ok)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: inverse(:,:)
    logical, intent(out) :: ok
    real(dp), allocatable :: aug(:,:)
    real(dp) :: pivot_value, factor
    integer :: n,i,j,pivot
    n=size(a,1)
    if (size(a,2)/=n) then
      allocate(inverse(0,0)); ok=.false.; return
    end if
    allocate(aug(n,2*n),inverse(n,n))
    aug(:,1:n)=a; aug(:,n+1:)=0.0_dp
    do i=1,n; aug(i,n+i)=1.0_dp; end do
    ok=.true.
    do i=1,n
      pivot=i
      do j=i+1,n
        if (abs(aug(j,i))>abs(aug(pivot,i))) pivot=j
      end do
      if (abs(aug(pivot,i))<=1.0e-12_dp) then
        ok=.false.; inverse=0.0_dp; return
      end if
      if (pivot/=i) call swap_rows(aug,i,pivot)
      pivot_value=aug(i,i); aug(i,:)=aug(i,:)/pivot_value
      do j=1,n
        if (j==i) cycle
        factor=aug(j,i); aug(j,:)=aug(j,:)-factor*aug(i,:)
      end do
    end do
    inverse=aug(:,n+1:)
  end subroutine invert_matrix

  subroutine swap_rows(a,i,j)
    real(dp), intent(inout) :: a(:,:)
    integer, intent(in) :: i,j
    real(dp), allocatable :: temp(:)
    allocate(temp(size(a,2)))
    temp=a(i,:); a(i,:)=a(j,:); a(j,:)=temp
  end subroutine swap_rows

  pure function diagonal(a) result(values)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable :: values(:)
    integer :: i,n
    n=min(size(a,1),size(a,2)); allocate(values(n))
    do i=1,n; values(i)=a(i,i); end do
  end function diagonal

end module svdnf_optimization
