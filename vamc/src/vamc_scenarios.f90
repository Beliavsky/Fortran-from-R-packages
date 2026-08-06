module vamc_scenarios
  use vamc_kinds, only: dp
  use vamc_status, only: status_type, vamc_invalid_argument, vamc_dimension_error
  use vamc_math, only: cholesky_upper, seed_rng, normal_random, uniform_random, sample_integer, random_permutation
  implicit none
  private
  public :: gen_index_scen, gen_fund_scen, random_fund_map
  interface gen_fund_scen
    module procedure gen_fund_scen_2d
    module procedure gen_fund_scen_3d
  end interface gen_fund_scen
contains
  subroutine gen_index_scen(covariance, num_scenarios, num_steps, dt, forward_curve, scenarios, seed, &
                            source_compatible_drift, status)
    real(dp), intent(in) :: covariance(:,:), dt, forward_curve(:)
    integer, intent(in) :: num_scenarios, num_steps
    real(dp), allocatable, intent(out) :: scenarios(:,:,:)
    integer, intent(in), optional :: seed
    logical, intent(in), optional :: source_compatible_drift
    type(status_type), intent(inout), optional :: status
    real(dp), allocatable :: chol(:,:), z(:,:), correlated(:,:), mu(:,:), variance_term(:)
    logical :: source_mode
    integer :: nindex, i, j, k
    type(status_type) :: local_status
    if (present(status)) call status%clear()
    if (num_scenarios < 1 .or. num_steps < 1 .or. dt <= 0.0_dp) then
      allocate(scenarios(0,0,0))
      if (present(status)) call status%set(vamc_invalid_argument, 'num_scenarios, num_steps, and dt must be positive.')
      return
    end if
    nindex = size(covariance,1)
    if (nindex < 1 .or. size(covariance,2) /= nindex .or. size(forward_curve) < num_steps) then
      allocate(scenarios(0,0,0))
      if (present(status)) call status%set(vamc_dimension_error, 'Scenario input dimensions do not conform.')
      return
    end if
    if (present(seed)) call seed_rng(seed)
    call cholesky_upper(covariance, chol, local_status)
    if (.not. local_status%ok()) then
      allocate(scenarios(0,0,0))
      if (present(status)) call status%set(local_status%code, local_status%message)
      return
    end if
    source_mode = .true.
    if (present(source_compatible_drift)) source_mode = source_compatible_drift
    allocate(variance_term(nindex), mu(num_steps,nindex), z(num_steps,nindex), correlated(num_steps,nindex))
    if (source_mode) then
      do j = 1, nindex
        variance_term(j) = 0.5_dp * sum(chol(j,:)**2)
      end do
    else
      do j = 1, nindex
        variance_term(j) = 0.5_dp * covariance(j,j)
      end do
    end if
    do j = 1, nindex
      mu(:,j) = (forward_curve(1:num_steps) - variance_term(j)) * dt
    end do
    allocate(scenarios(num_scenarios,num_steps,nindex))
    do i = 1, num_scenarios
      do j = 1, nindex
        do k = 1, num_steps
          z(k,j) = normal_random()
        end do
      end do
      correlated = matmul(z, chol)
      scenarios(i,:,:) = exp(mu + sqrt(dt) * correlated)
    end do
  end subroutine gen_index_scen

  subroutine gen_fund_scen_3d(fund_map, index_scenarios, fund_scenarios, status)
    real(dp), intent(in) :: fund_map(:,:), index_scenarios(:,:,:)
    real(dp), allocatable, intent(out) :: fund_scenarios(:,:,:)
    type(status_type), intent(inout), optional :: status
    integer :: nscen, nstep, nindex, nfund, i
    if (present(status)) call status%clear()
    nscen = size(index_scenarios,1)
    nstep = size(index_scenarios,2)
    nindex = size(index_scenarios,3)
    nfund = size(fund_map,1)
    if (size(fund_map,2) /= nindex) then
      allocate(fund_scenarios(0,0,0))
      if (present(status)) call status%set(vamc_dimension_error, 'Funds from index scenarios must align with fund map.')
      return
    end if
    allocate(fund_scenarios(nscen,nstep,nfund))
    do i = 1, nscen
      fund_scenarios(i,:,:) = matmul(index_scenarios(i,:,:), transpose(fund_map))
    end do
  end subroutine gen_fund_scen_3d

  subroutine gen_fund_scen_2d(fund_map, index_scenarios, fund_scenarios, status)
    real(dp), intent(in) :: fund_map(:,:), index_scenarios(:,:)
    real(dp), allocatable, intent(out) :: fund_scenarios(:,:)
    type(status_type), intent(inout), optional :: status
    if (present(status)) call status%clear()
    if (size(fund_map,2) /= size(index_scenarios,2)) then
      allocate(fund_scenarios(0,0))
      if (present(status)) call status%set(vamc_dimension_error, 'Funds from index scenarios must align with fund map.')
      return
    end if
    allocate(fund_scenarios(size(index_scenarios,1),size(fund_map,1)))
    fund_scenarios = matmul(index_scenarios, transpose(fund_map))
  end subroutine gen_fund_scen_2d

  subroutine random_fund_map(num_indices, num_random_funds, fund_map, seed, status)
    integer, intent(in) :: num_indices, num_random_funds
    real(dp), allocatable, intent(out) :: fund_map(:,:)
    integer, intent(in), optional :: seed
    type(status_type), intent(inout), optional :: status
    integer, allocatable :: order(:)
    integer :: i, j, nselect
    real(dp) :: total
    if (present(status)) call status%clear()
    if (num_indices < 1 .or. num_random_funds < 0) then
      allocate(fund_map(0,0))
      if (present(status)) call status%set(vamc_invalid_argument, 'Invalid fund-map dimensions.')
      return
    end if
    if (present(seed)) call seed_rng(seed)
    allocate(fund_map(num_indices+num_random_funds,num_indices), source=0.0_dp)
    do i = 1, num_indices
      fund_map(i,i) = 1.0_dp
    end do
    allocate(order(num_indices))
    do i = 1, num_random_funds
      order = [(j,j=1,num_indices)]
      call random_permutation(order)
      nselect = max(1, min(num_indices, ceiling(uniform_random()*real(num_indices,dp))))
      do j = 1, nselect
        fund_map(num_indices+i,order(j)) = real(sample_integer(1,10),dp)
      end do
      total = sum(fund_map(num_indices+i,:))
      fund_map(num_indices+i,:) = fund_map(num_indices+i,:) / total
    end do
  end subroutine random_fund_map
end module vamc_scenarios
