! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2012 Marius Hofert, Ivan Kojadinovic, Martin Maechler and Jun Yan
module copula_simulation
  use copula_kinds, only : dp, i8
  use copula_types
  use copula_special, only : normal_cdf, student_cdf, log_one_plus, exp_minus_one
  use copula_linalg, only : cholesky_lower
  use copula_random, only : seed_random, random_uniform, random_normal, random_gamma, &
    random_exponential, random_positive_stable
  use copula_families, only : conditional_cdf
  implicit none
  private
  public :: random_copula
contains
  subroutine random_copula(n, model, samples, ok, seed)
    integer, intent(in) :: n
    type(copula_model), intent(in) :: model
    real(dp), allocatable, intent(out) :: samples(:,:)
    logical, intent(out) :: ok
    integer(i8), intent(in), optional :: seed
    type(copula_model) :: base_model
    real(dp), allocatable :: l(:,:), z(:), x(:)
    real(dp) :: frailty, alpha, chi_square
    logical :: good
    integer :: i, j, d, count
    d = model%dimension
    allocate(samples(max(0,n),max(0,d)))
    samples = 0.0_dp
    ok = n >= 1 .and. model%valid()
    if (.not. ok) return
    if (present(seed)) call seed_random(seed)
    base_model = model
    base_model%rotation = rotation_none
    select case (model%family)
    case (family_independence)
      do i = 1, n
        do j = 1, d
          samples(i,j) = random_uniform()
        end do
      end do
    case (family_gaussian, family_student)
      call cholesky_lower(model%correlation,l,good)
      if (.not. good) then
        ok = .false.
        return
      end if
      allocate(z(d),x(d))
      do i = 1, n
        do j = 1, d
          z(j) = random_normal()
        end do
        x = matmul(l,z)
        if (model%family == family_gaussian) then
          samples(i,:) = normal_cdf(x)
        else
          chi_square = random_gamma(0.5_dp*model%df,2.0_dp)
          samples(i,:) = student_cdf(x/sqrt(chi_square/model%df),model%df)
        end if
      end do
    case (family_clayton)
      do i = 1, n
        frailty = random_gamma(1.0_dp/model%theta,1.0_dp)
        do j = 1, d
          samples(i,j) = (1.0_dp+random_exponential()/frailty)**(-1.0_dp/model%theta)
        end do
      end do
    case (family_gumbel)
      alpha = 1.0_dp/model%theta
      do i = 1, n
        frailty = random_positive_stable(alpha)
        do j = 1, d
          samples(i,j) = exp(-(random_exponential()/frailty)**alpha)
        end do
      end do
    case (family_frank)
      if (d > 2 .and. model%theta > 0.0_dp) then
        do i = 1, n
          count = random_logseries(1.0_dp-exp(-model%theta))
          do j = 1, d
            samples(i,j) = -log_one_plus(exp_minus_one(-model%theta)*exp(-random_exponential()/real(count,dp)))/model%theta
          end do
        end do
      else
        call random_bivariate_inversion(n,base_model,samples)
      end if
    case default
      if (d /= 2) then
        ok = .false.
        return
      end if
      call random_bivariate_inversion(n,base_model,samples)
    end select
    if (model%rotation /= rotation_none) then
      select case (model%rotation)
      case (rotation_90)
        samples(:,1) = 1.0_dp-samples(:,1)
      case (rotation_180)
        samples = 1.0_dp-samples
      case (rotation_270)
        samples(:,2) = 1.0_dp-samples(:,2)
      end select
    end if
  contains
    subroutine random_bivariate_inversion(number, cop, output)
      integer, intent(in) :: number
      type(copula_model), intent(in) :: cop
      real(dp), intent(out) :: output(:,:)
      integer :: row, it
      real(dp) :: lower, upper, middle, wanted, first
      do row = 1, number
        first = random_uniform()
        wanted = random_uniform()
        lower = epsilon(1.0_dp)
        upper = 1.0_dp-epsilon(1.0_dp)
        do it = 1, 70
          middle = 0.5_dp*(lower+upper)
          if (conditional_cdf(middle,first,cop) < wanted) then
            lower = middle
          else
            upper = middle
          end if
        end do
        output(row,1) = first
        output(row,2) = 0.5_dp*(lower+upper)
      end do
    end subroutine random_bivariate_inversion
  end subroutine random_copula

  integer function random_logseries(probability) result(k)
    real(dp), intent(in) :: probability
    real(dp) :: target, cumulative, mass
    k = 1
    target = random_uniform()
    mass = -probability/log(1.0_dp-probability)
    cumulative = mass
    do while (target > cumulative .and. k < 100000)
      k = k+1
      mass = mass*probability*real(k-1,dp)/real(k,dp)
      cumulative = cumulative+mass
    end do
  end function random_logseries
end module copula_simulation
