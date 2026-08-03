! SPDX-License-Identifier: GPL-3.0-or-later
module qbc_calibration
   use qbc_kinds, only : dp
   use qbc_status, only : qbc_success, qbc_invalid_argument, qbc_size_mismatch, qbc_no_convergence
   use qbc_types, only : qbc_bond, qbc_swap, qbc_curve, qbc_calibration_result, qbc_rate_continuous
   use qbc_bonds, only : valuation_bonds
   use qbc_swaps, only : valuation_swaps
   use qbc_curves, only : curve_from_market_yields
   use qbc_optimization, only : qbc_optimizer_result, bounded_nelder_mead
   implicit none
   private
   public :: curve_calibration, curve_calculation, bootstrap_curve, basis_curve

   type :: bond_calibration_data
      type(qbc_bond), allocatable :: bonds(:)
      real(dp), allocatable :: market_prices(:)
      real(dp), allocatable :: weights(:)
      real(dp), allocatable :: knot_terms(:)
      integer :: approximation = 2
   end type bond_calibration_data

   type :: basis_calibration_data
      type(qbc_swap), allocatable :: swaps(:)
      type(qbc_curve) :: local_curve
      type(qbc_curve) :: foreign_curve
      real(dp), allocatable :: market_values(:)
      real(dp), allocatable :: weights(:)
      real(dp), allocatable :: knot_terms(:)
      integer :: approximation = 2
   end type basis_calibration_data

contains

   subroutine curve_calibration(market_terms, market_yields, nodes, result, approximation, &
                                rate_type, frequency, forward_output, status)
      real(dp), intent(in) :: market_terms(:), market_yields(:), nodes(:)
      type(qbc_calibration_result), intent(out) :: result
      integer, intent(in), optional :: approximation, rate_type, frequency
      logical, intent(in), optional :: forward_output
      integer, intent(out), optional :: status
      integer :: st
      call curve_from_market_yields(market_terms, market_yields, nodes, approximation, rate_type, &
                                    frequency, result%curve, forward_output, st)
      result%objective = 0.0_dp
      result%convergence = merge(0, 1, st == qbc_success)
      result%iterations = 0
      if (st == qbc_success) then
         result%message = 'market yields interpolated'
      else
         result%message = 'invalid market-yield curve'
      end if
      if (present(status)) status = st
   end subroutine curve_calibration

   subroutine curve_calculation(market_terms, yield_history, nodes, curves, approximation, &
                                rate_type, frequency, forward_output, status)
      real(dp), intent(in) :: market_terms(:), yield_history(:, :), nodes(:)
      type(qbc_curve), allocatable, intent(out) :: curves(:)
      integer, intent(in), optional :: approximation, rate_type, frequency
      logical, intent(in), optional :: forward_output
      integer, intent(out), optional :: status
      type(qbc_calibration_result) :: fit
      integer :: i, st
      if (size(yield_history, 1) /= size(market_terms)) then
         allocate(curves(0)); st = qbc_size_mismatch
      else
         allocate(curves(size(yield_history, 2)))
         st = qbc_success
         do i = 1, size(yield_history, 2)
            call curve_calibration(market_terms, yield_history(:, i), nodes, fit, approximation, &
                                   rate_type, frequency, forward_output, st)
            curves(i) = fit%curve
            if (st /= qbc_success) exit
         end do
      end if
      if (present(status)) status = st
   end subroutine curve_calculation

   subroutine bootstrap_curve(bonds, market_prices, knot_terms, initial_rates, result, weights, &
                              approximation, lower, upper, status)
      type(qbc_bond), intent(in) :: bonds(:)
      real(dp), intent(in) :: market_prices(:), knot_terms(:), initial_rates(:)
      type(qbc_calibration_result), intent(out) :: result
      real(dp), intent(in), optional :: weights(:), lower(:), upper(:)
      integer, intent(in), optional :: approximation
      integer, intent(out), optional :: status
      type(bond_calibration_data) :: data
      type(qbc_optimizer_result) :: solution
      real(dp), allocatable :: lo(:), hi(:)
      integer :: n, st
      n = size(knot_terms)
      if (size(bonds) /= size(market_prices) .or. size(initial_rates) /= n .or. n == 0) then
         st = qbc_size_mismatch
         result%convergence = 1
         result%message = 'size mismatch'
         if (present(status)) status = st
         return
      end if
      allocate(data%bonds(size(bonds)), data%market_prices(size(market_prices)), &
               data%weights(size(market_prices)), data%knot_terms(n))
      data%bonds = bonds
      data%market_prices = market_prices
      data%weights = 1.0_dp
      if (present(weights)) then
         if (size(weights) /= size(market_prices)) then
            st = qbc_size_mismatch
            if (present(status)) status = st
            return
         end if
         data%weights = weights
      end if
      data%knot_terms = knot_terms
      data%approximation = 2
      if (present(approximation)) data%approximation = approximation

      allocate(lo(n),hi(n))
      lo=-0.25_dp;hi=2.0_dp
      if (present(lower)) then
         if (size(lower)==n) lo=lower
      end if
      if (present(upper)) then
         if (size(upper)==n) hi=upper
      end if
      call bounded_nelder_mead(bond_objective,initial_rates,lo,hi,solution,data,1.0e-9_dp,4000)
      allocate(result%curve%terms(n), result%curve%rates(n))
      result%curve%terms = knot_terms
      if (allocated(solution%x)) then
         result%curve%rates = solution%x
      else
         result%curve%rates = initial_rates
      end if
      result%curve%approximation = data%approximation
      result%curve%rate_type = qbc_rate_continuous
      result%curve%frequency = 1
      result%objective = solution%value
      result%iterations = solution%iterations
      result%convergence = solution%status
      result%message = solution%message
      st = solution%status
      if (present(status)) status = st
   end subroutine bootstrap_curve

   subroutine bond_objective(x, value, raw_data)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: value
      class(*), intent(in), optional :: raw_data
      type(qbc_curve) :: curve
      real(dp) :: model_price, residual
      integer :: i
      value = huge(1.0_dp)
      if (.not. present(raw_data)) return
      select type (data => raw_data)
      type is (bond_calibration_data)
         allocate(curve%terms(size(data%knot_terms)), curve%rates(size(x)))
         curve%terms = data%knot_terms
         curve%rates = x
         curve%approximation = data%approximation
         curve%rate_type = qbc_rate_continuous
         curve%frequency = 1
         value = 0.0_dp
         do i = 1, size(data%bonds)
            model_price = valuation_bonds(data%bonds(i), [data%bonds(i)%coupon_rate], curve)
            residual = model_price - data%market_prices(i)
            value = value + data%weights(i) * residual * residual
         end do
      class default
         value = huge(1.0_dp)
      end select
   end subroutine bond_objective

   subroutine basis_curve(swaps, local_curve, foreign_curve, market_values, knot_terms, &
                          initial_rates, result, weights, approximation, status)
      type(qbc_swap), intent(in) :: swaps(:)
      type(qbc_curve), intent(in) :: local_curve, foreign_curve
      real(dp), intent(in) :: market_values(:), knot_terms(:), initial_rates(:)
      type(qbc_calibration_result), intent(out) :: result
      real(dp), intent(in), optional :: weights(:)
      integer, intent(in), optional :: approximation
      integer, intent(out), optional :: status
      type(basis_calibration_data) :: data
      type(qbc_optimizer_result) :: solution
      real(dp), allocatable :: lo(:), hi(:)
      integer :: n, st
      n = size(knot_terms)
      if (size(swaps) /= size(market_values) .or. size(initial_rates) /= n .or. n == 0) then
         st = qbc_size_mismatch
         result%message = 'size mismatch'
         if (present(status)) status = st
         return
      end if
      allocate(data%swaps(size(swaps)), data%market_values(size(market_values)), &
               data%weights(size(market_values)), data%knot_terms(n))
      data%swaps = swaps
      data%local_curve = local_curve
      data%foreign_curve = foreign_curve
      data%market_values = market_values
      data%weights = 1.0_dp
      if (present(weights)) then
         if (size(weights) == size(market_values)) data%weights = weights
      end if
      data%knot_terms = knot_terms
      data%approximation = 2
      if (present(approximation)) data%approximation = approximation

      allocate(lo(n),hi(n))
      lo=-1.0_dp;hi=2.0_dp
      call bounded_nelder_mead(basis_objective,initial_rates,lo,hi,solution,data,1.0e-9_dp,4000)
      allocate(result%curve%terms(n), result%curve%rates(n))
      result%curve%terms = knot_terms
      if (allocated(solution%x)) then
         result%curve%rates = solution%x
      else
         result%curve%rates = initial_rates
      end if
      result%curve%approximation = data%approximation
      result%curve%rate_type = qbc_rate_continuous
      result%curve%frequency = 1
      result%objective = solution%value
      result%iterations = solution%iterations
      result%convergence = solution%status
      result%message = solution%message
      st = solution%status
      if (present(status)) status = st
   end subroutine basis_curve

   subroutine basis_objective(x, value, raw_data)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: value
      class(*), intent(in), optional :: raw_data
      type(qbc_curve) :: basis
      real(dp) :: model_value, residual
      integer :: i
      value = huge(1.0_dp)
      if (.not. present(raw_data)) return
      select type (data => raw_data)
      type is (basis_calibration_data)
         allocate(basis%terms(size(data%knot_terms)), basis%rates(size(x)))
         basis%terms = data%knot_terms
         basis%rates = x
         basis%approximation = data%approximation
         basis%rate_type = qbc_rate_continuous
         basis%frequency = 1
         value = 0.0_dp
         do i = 1, size(data%swaps)
            model_value = valuation_swaps(data%swaps(i), data%local_curve, data%foreign_curve, basis)
            residual = model_value - data%market_values(i)
            value = value + data%weights(i) * residual * residual
         end do
      class default
         value = huge(1.0_dp)
      end select
   end subroutine basis_objective

end module qbc_calibration
