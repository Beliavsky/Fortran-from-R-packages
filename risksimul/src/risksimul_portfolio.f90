! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
! Derived from riskSimul 0.1.2 by Wolfgang Hormann and Ismail Basoglu.
module risksimul_portfolio
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use ghyp_kinds, only : dp
   use ghyp_linalg, only : cholesky_lower
   use ghyp_special, only : student_cdf
   use ghyp_model, only : ghyp_ad, ghyp_moments, moments_result
   use ghyp_distribution, only : dghyp, qghyp
   use risksimul_types, only : portfolio_model, importance_parameters, &
      marginal_t, marginal_gh
   use risksimul_math, only : student_quantile, inverse_table_quantile
   implicit none
   private

   public :: new_portfolio, new_portfobj
   public :: return_copula, portfolio_return_one
   public :: tail_loss_response, excess_response
   public :: algorithm_2, algorithm_3, touch_value

   interface return_copula
      module procedure return_copula_matrix
      module procedure portfolio_return_one
   end interface return_copula

   interface new_portfobj
      module procedure new_portfolio
   end interface new_portfobj

contains

   function lowercase(text) result(out)
      character(len=*), intent(in) :: text
      character(len=len(text)) :: out
      integer :: i, code
      do i = 1, len(text)
         code = iachar(text(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            out(i:i) = achar(code+32)
         else
            out(i:i) = text(i:i)
         end if
      end do
   end function lowercase

   function new_portfolio(nu, correlation, marginal_type, parameters, scale, &
      weight, gh_grid_size) result(model)
      real(dp), intent(in) :: nu
      real(dp), intent(in) :: correlation(:,:), parameters(:,:)
      character(len=*), intent(in) :: marginal_type
      real(dp), intent(in), optional :: scale(:), weight(:)
      integer, intent(in), optional :: gh_grid_size
      type(portfolio_model) :: model
      integer :: d, i, grid_size
      logical :: chol_ok
      character(len=:), allocatable :: family

      d = size(correlation,1)
      if (d < 1 .or. size(correlation,2) /= d) then
         model%message = 'correlation matrix must be square and nonempty'
         return
      end if
      if (nu <= 2.0_dp .or. .not. ieee_is_finite(nu)) then
         model%message = 'copula degrees of freedom must exceed two'
         return
      end if
      if (.not. all(ieee_is_finite(correlation)) .or. &
          maxval(abs(correlation-transpose(correlation))) > 1.0e-10_dp) then
         model%message = 'correlation matrix must be finite and symmetric'
         return
      end if
      if (maxval(abs([(correlation(i,i)-1.0_dp,i=1,d)])) > 1.0e-8_dp) then
         model%message = 'correlation matrix diagonal must equal one'
         return
      end if
      call cholesky_lower(correlation,model%cholesky,chol_ok)
      if (.not. chol_ok) then
         model%message = 'correlation matrix must be positive definite'
         return
      end if

      allocate(model%correlation(d,d),model%marginal_parameters(size(parameters,1), &
         size(parameters,2)),model%scale(d),model%weight(d))
      model%correlation = correlation
      model%marginal_parameters = parameters
      model%copula_df = nu
      model%scale = 1.0_dp
      if (present(scale)) then
         if (size(scale) /= d .or. .not. all(ieee_is_finite(scale))) then
            model%message = 'scale vector has incompatible size or nonfinite values'
            return
         end if
         model%scale = scale
      end if
      if (present(weight)) then
         if (size(weight) /= d .or. .not. all(ieee_is_finite(weight))) then
            model%message = 'weight vector has incompatible size or nonfinite values'
            return
         end if
         model%weight = weight
      else
         if (abs(sum(model%scale)) <= tiny(1.0_dp)) then
            model%message = 'default weights require a nonzero sum of scale values'
            return
         end if
         model%weight = model%scale/sum(model%scale)
      end if
      if (abs(sum(model%weight)) <= tiny(1.0_dp)) then
         model%message = 'portfolio weights must have a nonzero sum'
         return
      end if

      family = trim(lowercase(adjustl(marginal_type)))
      select case(family)
      case('t','student','student-t')
         model%marginal_family = marginal_t
         if (size(parameters,1) /= d .or. size(parameters,2) /= 3) then
            model%message = 't marginal parameters must be D by 3: mu, sigma, df'
            return
         end if
         if (any(parameters(:,2) <= 0.0_dp) .or. any(parameters(:,3) <= 0.0_dp) .or. &
             .not. all(ieee_is_finite(parameters))) then
            model%message = 't marginal scales and degrees of freedom must be positive'
            return
         end if
      case('gh','generalized hyperbolic','generalized-hyperbolic')
         model%marginal_family = marginal_gh
         if (size(parameters,1) /= d .or. size(parameters,2) /= 5) then
            model%message = 'GH parameters must be D by 5: lambda, alpha, beta, delta, mu'
            return
         end if
         allocate(model%gh_models(d),model%gh_tables(d))
         grid_size = 2049
         if (present(gh_grid_size)) grid_size = max(129,gh_grid_size)
         do i = 1, d
            model%gh_models(i) = ghyp_ad(parameters(i,1),parameters(i,2), &
               parameters(i,4),[parameters(i,3)],[parameters(i,5)], &
               reshape([1.0_dp],[1,1]))
            if (.not. model%gh_models(i)%ok) then
               model%message = 'invalid generalized-hyperbolic marginal parameters'
               return
            end if
            call build_inverse_table(model%gh_models(i),grid_size,model%gh_tables(i)%x, &
               model%gh_tables(i)%cdf,model%gh_tables(i)%ready)
            if (.not. model%gh_tables(i)%ready) then
               model%message = 'could not construct GH inverse table'
               return
            end if
         end do
      case default
         model%message = 'marginal type must be t or GH'
         return
      end select
      model%ok = .true.
   end function new_portfolio

   subroutine build_inverse_table(gh_model, grid_size, x, cdf, ok)
      use ghyp_model, only : ghyp_model_type
      type(ghyp_model_type), intent(in) :: gh_model
      integer, intent(in) :: grid_size
      real(dp), allocatable, intent(out) :: x(:), cdf(:)
      logical, intent(out) :: ok
      type(moments_result) :: moments
      real(dp) :: lo, hi, dx, area
      integer :: i

      moments = ghyp_moments(gh_model)
      if (.not. moments%ok .or. .not. ieee_is_finite(moments%mean(1)) .or. &
          moments%covariance(1,1) <= 0.0_dp) then
         ok = .false.
         return
      end if
      lo = qghyp(1.0e-8_dp,gh_model)
      hi = qghyp(1.0_dp-1.0e-8_dp,gh_model)
      if (.not. ieee_is_finite(lo) .or. .not. ieee_is_finite(hi) .or. hi <= lo) then
         lo = moments%mean(1)-16.0_dp*sqrt(moments%covariance(1,1))
         hi = moments%mean(1)+16.0_dp*sqrt(moments%covariance(1,1))
      end if
      allocate(x(grid_size),cdf(grid_size))
      dx = (hi-lo)/real(grid_size-1,dp)
      do i = 1, grid_size
         x(i) = lo+real(i-1,dp)*dx
      end do
      cdf(1) = 0.0_dp
      do i = 2, grid_size
         area = 0.5_dp*dx*(dghyp(x(i-1),gh_model)+dghyp(x(i),gh_model))
         cdf(i) = cdf(i-1)+max(area,0.0_dp)
      end do
      if (cdf(grid_size) <= tiny(1.0_dp)) then
         ok = .false.
         return
      end if
      cdf = cdf/cdf(grid_size)
      cdf(1) = 0.0_dp
      cdf(grid_size) = 1.0_dp
      ok = .true.
   end subroutine build_inverse_table

   function marginal_quantile(model, asset, probability) result(value)
      type(portfolio_model), intent(in) :: model
      integer, intent(in) :: asset
      real(dp), intent(in) :: probability
      real(dp) :: value, p
      p = min(1.0_dp-epsilon(1.0_dp),max(epsilon(1.0_dp),probability))
      if (model%marginal_family == marginal_t) then
         value = model%marginal_parameters(asset,1)+ &
            model%marginal_parameters(asset,2)*student_quantile(p, &
            model%marginal_parameters(asset,3))
      else
         value = inverse_table_quantile(p,model%gh_tables(asset)%x, &
            model%gh_tables(asset)%cdf)
      end if
   end function marginal_quantile

   function portfolio_return_one(z, y, model) result(value)
      real(dp), intent(in) :: z(:), y
      type(portfolio_model), intent(in) :: model
      real(dp) :: value
      real(dp), allocatable :: transformed(:)
      real(dp) :: t_value, probability, marginal
      integer :: i, d

      d = model%dimension()
      if (.not. model%ok .or. size(z) /= d .or. y <= 0.0_dp) then
         value = 0.0_dp
         return
      end if
      allocate(transformed(d))
      transformed = matmul(model%cholesky,z)/sqrt(y/model%copula_df)
      value = 0.0_dp
      do i = 1, d
         t_value = transformed(i)
         probability = student_cdf(t_value,model%copula_df)
         marginal = marginal_quantile(model,i,probability)
         value = value+model%weight(i)*exp(model%scale(i)*marginal)
      end do
   end function portfolio_return_one

   function return_copula_matrix(z, y, model) result(values)
      real(dp), intent(in) :: z(:,:), y(:)
      type(portfolio_model), intent(in) :: model
      real(dp), allocatable :: values(:)
      integer :: n, j

      n = size(z,2)
      allocate(values(n))
      if (size(z,1) /= model%dimension() .or. size(y) /= n) then
         values = 0.0_dp
         return
      end if
      do j = 1, n
         values(j) = portfolio_return_one(z(:,j),y(j),model)
      end do
   end function return_copula_matrix

   pure function tail_loss_response(portfolio_return, threshold) result(value)
      real(dp), intent(in) :: portfolio_return, threshold
      real(dp) :: value
      if (portfolio_return < threshold) then
         value = 1.0_dp
      else
         value = 0.0_dp
      end if
   end function tail_loss_response

   pure function excess_response(portfolio_return, threshold) result(value)
      real(dp), intent(in) :: portfolio_return, threshold
      real(dp) :: value
      if (portfolio_return < threshold) then
         value = 1.0_dp-portfolio_return
      else
         value = 0.0_dp
      end if
   end function excess_response

   function touch_value(r, direction, model, threshold) result(value)
      real(dp), intent(in) :: r, direction(:), threshold
      type(portfolio_model), intent(in) :: model
      real(dp) :: value
      value = portfolio_return_one(r*direction,model%copula_df,model)-threshold+1.0e-5_dp
   end function touch_value

   function algorithm_2(direction, model, threshold) result(result)
      real(dp), intent(in) :: direction(:), threshold
      type(portfolio_model), intent(in) :: model
      type(importance_parameters) :: result
      real(dp), allocatable :: unit_direction(:)
      real(dp) :: norm_value, lo, hi, mid, flo, fhi, fmid, r0, y0
      integer :: iter

      norm_value = sqrt(dot_product(direction,direction))
      if (.not. model%ok .or. size(direction) /= model%dimension() .or. &
          norm_value <= tiny(1.0_dp)) then
         result%message = 'invalid importance-sampling direction'
         return
      end if
      unit_direction = direction/norm_value
      lo = -1000.0_dp
      hi = 1.0e-5_dp
      flo = touch_value(lo,unit_direction,model,threshold)
      fhi = touch_value(hi,unit_direction,model,threshold)
      do iter = 1, 30
         if (flo*fhi <= 0.0_dp) exit
         lo = 2.0_dp*lo
         flo = touch_value(lo,unit_direction,model,threshold)
      end do
      if (flo*fhi > 0.0_dp) then
         result%message = 'could not bracket the rare-event boundary'
         return
      end if
      do iter = 1, 160
         mid = 0.5_dp*(lo+hi)
         fmid = touch_value(mid,unit_direction,model,threshold)
         if (flo*fmid <= 0.0_dp) then
            hi = mid
            fhi = fmid
         else
            lo = mid
            flo = fmid
         end if
         if (abs(hi-lo) <= 1.0e-10_dp*max(1.0_dp,abs(mid))) exit
      end do
      r0 = 0.5_dp*(lo+hi)
      y0 = (model%copula_df-2.0_dp)/(1.0_dp+r0*r0/model%copula_df)
      if (y0 <= 0.0_dp) then
         result%message = 'direction optimization produced an invalid gamma mean'
         return
      end if
      allocate(result%shift(model%dimension()))
      result%shift = r0*sqrt(y0/model%copula_df)*unit_direction
      result%gamma_mean = y0
      result%objective = (0.5_dp*model%copula_df-1.0_dp)*(log(y0)-1.0_dp)
      result%ok = .true.
   end function algorithm_2

   function algorithm_3(model, threshold, max_iterations) result(result)
      type(portfolio_model), intent(in) :: model
      real(dp), intent(in) :: threshold
      integer, intent(in), optional :: max_iterations
      type(importance_parameters) :: result
      type(importance_parameters) :: trial_result, best_result
      real(dp), allocatable :: direction(:), trial(:)
      real(dp) :: norm_value, step, best_value, trial_value
      integer :: d, i, iter, n_iter, sign_index
      logical :: improved

      if (.not. model%ok) then
         result%message = 'invalid portfolio model'
         return
      end if
      d = model%dimension()
      allocate(direction(d),trial(d))
      direction = matmul(model%cholesky,model%scale*model%weight)
      direction = abs(direction)
      norm_value = sqrt(dot_product(direction,direction))
      if (norm_value <= tiny(1.0_dp)) direction = 1.0_dp
      direction = direction/sqrt(dot_product(direction,direction))
      best_result = algorithm_2(direction,model,threshold)
      if (.not. best_result%ok) then
         result = best_result
         return
      end if
      best_value = -best_result%objective
      step = 0.25_dp
      n_iter = 160
      if (present(max_iterations)) n_iter = max(1,max_iterations)

      do iter = 1, n_iter
         improved = .false.
         do i = 1, d
            do sign_index = -1, 1, 2
               trial = direction
               trial(i) = max(0.0_dp,trial(i)+real(sign_index,dp)*step)
               norm_value = sqrt(dot_product(trial,trial))
               if (norm_value <= tiny(1.0_dp)) cycle
               trial = trial/norm_value
               trial_result = algorithm_2(trial,model,threshold)
               if (.not. trial_result%ok) cycle
               trial_value = -trial_result%objective
               if (trial_value < best_value) then
                  direction = trial
                  best_result = trial_result
                  best_value = trial_value
                  improved = .true.
               end if
            end do
         end do
         if (.not. improved) step = 0.5_dp*step
         if (step < 1.0e-5_dp) exit
      end do
      result = best_result
   end function algorithm_3

end module risksimul_portfolio
