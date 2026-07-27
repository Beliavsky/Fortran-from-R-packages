! SPDX-License-Identifier: GPL-3.0-only
module svdnf_filter
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use svdnf_kinds, only : dp
  use svdnf_types, only : svm_dynamics, grid_type, filter_result, percentile_result
  use svdnf_models, only : validate_dynamics, evaluate_mu_y, evaluate_sigma_y, &
    evaluate_mu_x, evaluate_sigma_x, jump_probability
  use svdnf_grids, only : grid_maker, validate_grid
  use svdnf_stats, only : normal_pdf, normal_cdf, gamma_cdf
  implicit none
  private
  public :: dnf_filter, dnf, probability_components
  public :: transition_matrix, extract_vol_percentile, extract_vol_perc

contains

  function dnf_filter(dynamics, data, factors, grids, n, k, r) result(output)
    type(svm_dynamics), intent(in) :: dynamics
    real(dp), intent(in) :: data(:)
    real(dp), intent(in), optional :: factors(:,:)
    type(grid_type), intent(in), optional :: grids
    integer, intent(in), optional :: n, k, r
    type(filter_result) :: output
    type(grid_type) :: work_grid
    real(dp), allocatable :: adjusted(:), transition(:,:,:,:), mean_y(:,:,:,:), sd_y(:,:,:,:)
    real(dp), allocatable :: numerator(:)
    integer :: nstate, njump, ncount, tmax, t, i, j, m, q
    real(dp) :: likelihood, density
    logical :: ok
    character(len=160) :: message

    call validate_dynamics(dynamics,ok,message)
    if (.not. ok) then
      output%message = message
      return
    end if
    if (size(data) < 1) then
      output%message = 'At least one return observation is required.'
      return
    end if
    if (present(grids)) then
      work_grid = grids
    else
      work_grid = grid_maker(dynamics,n,k,r)
    end if
    call validate_grid(work_grid,ok,message)
    if (.not. ok) then
      output%message = message
      return
    end if

    adjusted = data
    if (present(factors)) then
      if (.not. allocated(dynamics%coefs)) then
        output%message = 'Factors were supplied but the dynamics have no coefficients.'
        return
      end if
      if (size(factors,1) /= size(data) .or. size(factors,2) /= size(dynamics%coefs)) then
        output%message = 'The factor matrix must be observations by coefficients.'
        return
      end if
      adjusted = adjusted - matmul(factors,dynamics%coefs)
    else if (allocated(dynamics%coefs)) then
      if (size(dynamics%coefs) > 0) then
        output%message = 'Factor coefficients are present but no factors were supplied.'
        return
      end if
    end if

    nstate = size(work_grid%var_mid_points)
    njump = size(work_grid%jump_mid_points)
    ncount = size(work_grid%jump_counts)
    tmax = size(data)
    allocate(transition(nstate,nstate,njump,ncount))
    allocate(mean_y(nstate,nstate,njump,ncount),sd_y(nstate,nstate,njump,ncount))
    call probability_components(dynamics,work_grid,transition,mean_y,sd_y)

    allocate(output%likelihoods(tmax),output%filter_grid(nstate,tmax+1),numerator(nstate))
    output%filter_grid(:,1) = 1.0_dp/real(nstate,dp)
    output%log_likelihood = 0.0_dp
    do t = 1, tmax
      numerator = 0.0_dp
      do q = 1, ncount
        do m = 1, njump
          do j = 1, nstate
            do i = 1, nstate
              density = normal_pdf(adjusted(t),mean_y(i,j,m,q),sd_y(i,j,m,q))
              numerator(i) = numerator(i) + output%filter_grid(j,t)*transition(i,j,m,q)*density
            end do
          end do
        end do
      end do
      likelihood = sum(numerator)
      output%likelihoods(t) = likelihood
      if (.not. ieee_is_finite(likelihood)) then
        output%message = 'A non-finite likelihood contribution was produced.'
        output%log_likelihood = -huge(1.0_dp)
        output%grids = work_grid
        output%dynamics = dynamics
        output%data = data
        return
      end if
      if (likelihood <= tiny(1.0_dp)) then
        output%filter_grid(:,t+1) = output%filter_grid(:,t)
        output%log_likelihood = -huge(1.0_dp)
        output%message = 'A likelihood contribution underflowed to zero.'
        output%grids = work_grid
        output%dynamics = dynamics
        output%data = data
        return
      end if
      output%filter_grid(:,t+1) = numerator/likelihood
      output%log_likelihood = output%log_likelihood + log(likelihood)
    end do
    output%grids = work_grid
    output%dynamics = dynamics
    output%data = data
    output%ok = .true.
  end function dnf_filter

  function dnf(dynamics, data, factors, grids, n, k, r) result(output)
    type(svm_dynamics), intent(in) :: dynamics
    real(dp), intent(in) :: data(:)
    real(dp), intent(in), optional :: factors(:,:)
    type(grid_type), intent(in), optional :: grids
    integer, intent(in), optional :: n, k, r
    type(filter_result) :: output
    if (present(factors) .and. present(grids)) then
      output = dnf_filter(dynamics,data,factors=factors,grids=grids,n=n,k=k,r=r)
    else if (present(factors)) then
      output = dnf_filter(dynamics,data,factors=factors,n=n,k=k,r=r)
    else if (present(grids)) then
      output = dnf_filter(dynamics,data,grids=grids,n=n,k=k,r=r)
    else
      output = dnf_filter(dynamics,data,n=n,k=k,r=r)
    end if
  end function dnf

  subroutine probability_components(dynamics, grids, transition, mean_y, sd_y)
    type(svm_dynamics), intent(in) :: dynamics
    type(grid_type), intent(in) :: grids
    real(dp), intent(out) :: transition(:,:,:,:), mean_y(:,:,:,:), sd_y(:,:,:,:)
    real(dp), allocatable :: var_bounds(:), jump_bounds(:)
    real(dp) :: x, xprev, jump_size, mux, sigx, muy, sigy, eps, qv, qj, pn
    integer :: i, j, m, q, count, nstate, njump, ncount
    nstate = size(grids%var_mid_points)
    njump = size(grids%jump_mid_points)
    ncount = size(grids%jump_counts)
    call make_state_bounds(grids%var_mid_points,var_bounds)
    call make_jump_bounds(grids%jump_mid_points,jump_bounds)
    do q = 1, ncount
      count = grids%jump_counts(q)
      pn = jump_probability(dynamics,count)
      do m = 1, njump
        jump_size = grids%jump_mid_points(m)
        if (njump == 1 .and. abs(jump_size) <= epsilon(1.0_dp)) then
          qj = 1.0_dp
        else if (count == 0) then
          if (jump_bounds(m) <= 0.0_dp .and. jump_bounds(m+1) >= 0.0_dp) then
            qj = 1.0_dp
          else
            qj = 0.0_dp
          end if
        else if (dynamics%nu > 0.0_dp) then
          qj = gamma_cdf(jump_bounds(m+1),real(count,dp),dynamics%nu) - &
            gamma_cdf(max(jump_bounds(m),0.0_dp),real(count,dp),dynamics%nu)
        else
          qj = 0.0_dp
        end if
        do j = 1, nstate
          xprev = grids%var_mid_points(j)
          mux = evaluate_mu_x(dynamics,xprev)
          sigx = evaluate_sigma_x(dynamics,xprev)
          muy = evaluate_mu_y(dynamics,xprev)
          sigy = evaluate_sigma_y(dynamics,xprev)
          do i = 1, nstate
            x = grids%var_mid_points(i)
            qv = normal_cdf(var_bounds(i+1),mux+jump_size,sigx) - &
              normal_cdf(var_bounds(i),mux+jump_size,sigx)
            transition(i,j,m,q) = max(qv,0.0_dp)*max(pn,0.0_dp)*max(qj,0.0_dp)
            if (sigx > tiny(1.0_dp)) then
              eps = (x-mux-jump_size)/sigx
            else
              eps = 0.0_dp
            end if
            mean_y(i,j,m,q) = muy + dynamics%rho*sigy*eps + &
              dynamics%alpha*real(count,dp) + dynamics%rho_z*jump_size
            sd_y(i,j,m,q) = sqrt(max((1.0_dp-dynamics%rho**2)*sigy**2 + &
              real(count,dp)*dynamics%delta**2,tiny(1.0_dp)))
          end do
        end do
      end do
    end do
  end subroutine probability_components

  function transition_matrix(dynamics, grids) result(matrix_value)
    type(svm_dynamics), intent(in) :: dynamics
    type(grid_type), intent(in) :: grids
    real(dp), allocatable :: matrix_value(:,:)
    real(dp), allocatable :: transition(:,:,:,:), mean_y(:,:,:,:), sd_y(:,:,:,:)
    integer :: nstate, njump, ncount, j, m, q
    nstate = size(grids%var_mid_points)
    njump = size(grids%jump_mid_points)
    ncount = size(grids%jump_counts)
    allocate(transition(nstate,nstate,njump,ncount))
    allocate(mean_y(nstate,nstate,njump,ncount),sd_y(nstate,nstate,njump,ncount))
    call probability_components(dynamics,grids,transition,mean_y,sd_y)
    allocate(matrix_value(nstate,nstate))
    matrix_value = 0.0_dp
    do q = 1, ncount
      do m = 1, njump
        matrix_value = matrix_value + transition(:,:,m,q)
      end do
    end do
    do j = 1, nstate
      if (sum(matrix_value(:,j)) > 0.0_dp) matrix_value(:,j)=matrix_value(:,j)/sum(matrix_value(:,j))
    end do
  end function transition_matrix

  function extract_vol_percentile(filtered, p, prediction) result(output)
    type(filter_result), intent(in) :: filtered
    real(dp), intent(in), optional :: p
    logical, intent(in), optional :: prediction
    type(percentile_result) :: output
    real(dp) :: probability, cumulative, total
    real(dp), allocatable :: trans(:,:)
    integer :: t, i, nt, nstate
    logical :: use_prediction
    if (.not. filtered%ok) then
      output%message = 'The supplied filter result is not valid.'
      return
    end if
    probability = 0.5_dp
    if (present(p)) probability = p
    if (probability < 0.0_dp .or. probability > 1.0_dp) then
      output%message = 'The percentile probability must be between zero and one.'
      return
    end if
    use_prediction = .false.
    if (present(prediction)) use_prediction = prediction
    nstate = size(filtered%filter_grid,1)
    if (use_prediction) then
      nt = size(filtered%filter_grid,2)-1
      trans = transition_matrix(filtered%dynamics,filtered%grids)
      allocate(output%distributions(nstate,nt))
      do t = 1, nt
        output%distributions(:,t) = matmul(trans,filtered%filter_grid(:,t))
      end do
    else
      nt = size(filtered%filter_grid,2)
      allocate(output%distributions(nstate,nt))
      output%distributions = filtered%filter_grid
    end if
    allocate(output%values(nt))
    do t = 1, nt
      total = sum(output%distributions(:,t))
      if (total <= 0.0_dp) then
        output%values(t) = filtered%grids%var_mid_points(1)
      else
        output%distributions(:,t) = output%distributions(:,t)/total
        cumulative = 0.0_dp
        output%values(t) = filtered%grids%var_mid_points(nstate)
        do i = 1, nstate
          cumulative = cumulative + output%distributions(i,t)
          if (cumulative >= probability) then
            output%values(t) = filtered%grids%var_mid_points(i)
            exit
          end if
        end do
      end if
    end do
    output%prediction = use_prediction
    output%ok = .true.
  end function extract_vol_percentile

  function extract_vol_perc(filtered, p, pred) result(output)
    type(filter_result), intent(in) :: filtered
    real(dp), intent(in), optional :: p
    logical, intent(in), optional :: pred
    type(percentile_result) :: output
    output = extract_vol_percentile(filtered,p,pred)
  end function extract_vol_perc

  subroutine make_state_bounds(points,bounds)
    real(dp), intent(in) :: points(:)
    real(dp), allocatable, intent(out) :: bounds(:)
    integer :: i, n
    n = size(points)
    allocate(bounds(n+1))
    bounds(1) = points(1)-0.5_dp*(points(2)-points(1))
    do i = 2, n
      bounds(i)=0.5_dp*(points(i-1)+points(i))
    end do
    bounds(n+1)=huge(1.0_dp)
  end subroutine make_state_bounds

  subroutine make_jump_bounds(points,bounds)
    real(dp), intent(in) :: points(:)
    real(dp), allocatable, intent(out) :: bounds(:)
    integer :: i, n
    n=size(points)
    allocate(bounds(n+1))
    if (n == 1) then
      bounds(1)=-huge(1.0_dp)
      bounds(2)=huge(1.0_dp)
      return
    end if
    bounds(1)=0.5_dp*(-points(1)+points(1))
    do i=2,n
      bounds(i)=0.5_dp*(points(i-1)+points(i))
    end do
    bounds(n+1)=huge(1.0_dp)
  end subroutine make_jump_bounds

end module svdnf_filter
