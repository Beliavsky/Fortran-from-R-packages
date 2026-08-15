!-------------------------------------------------------------------------------
!
! MODULE: Delaporte
!
! AUTHOR: Avraham Adler <Avraham.Adler@gmail.com>
!
! DESCRIPTION: Probability mass, cumulative distribution, and quantile
!              functions for the Delaporte distribution. Random number
!              generation and method of moments functions as well.
!
! HISTORY:
!          Version 1.0: 2016-11-20
!                       Porting from C++ code in Delaporte package for R.
!          Version 1.1: 2017-03-01
!                       Various tweaks.
!          Version 1.2: 2017-08-13
!                       Corrected MoMdelap code.
!          Version 1.3: 2017-11-20
!                       Updates.
!          Version 1.4: 2018-06-18
!                       Added skew bias correction option to MoMdelap.
!          Version 1.5: 2018-11-20
!                       Zapping absolute values <= EPS to 0.
!          Version 1.6: 2018-12-10
!                       Replaced zapping with setting min to 0 and max to 1
!                       as appropriate. Less monkeying with values this way.
!          Version 2.0: 2021-01-03
!                       Setting limits as < 0 to be more consistent with R
!                       defaults for d/p/q/r functions. Returning NaN for NaN
!                       inputs per R defaults. Trapping for INFTY more
!                       consistently with base R. Should use ieee_arithmetic
!                       once current oldrelease gets deprecated and min GCC
!                       version is > 5. Using iso_fortran_env for INT64 to allow
!                       wider domain for d/pdelap.
!          Version 3.0: 2023-01-29
!                       Updated to rely on Fortran 2008 intrinsics and use
!                       ieee_artithmetic.
!          Version 4.0: 2023-08-08
!                       Added OpenMP thread control functionality.
!          Version 4.1: 2024-04-04
!                       Use "source" when allocating arrays in qdelap_f.
!          Version 5.0: 2024-06-17
!                       OpenMP is still significantly faster than extending
!                       parameter vectors and applying the elemental singleton
!                       function. Use imk helper function to calculate position
!                       for vector recycling. Turn floating point error cleaner
!                       into a function. Move lt and lg for p/delap into
!                       existing loops. While elemental functions can be even
!                       faster than OpenMP loops, they are still another loop.
!                       Saving the overhead by calling the conditionals inside
!                       the necessary loops is still faster. This is not useful
!                       for qdelap as we need to know every "real" percentile
!                       first. So run those before anything else. Other minor
!                       tweaks to code for efficiency. SIMD instructions do not
!                       save that much time and prevent compilation with current
!                       versions of flang (through 19) as it does not have full
!                       OpenMP 4.5 implementation yet.
!          Version 5.1: 2025-07-17
!                       Made checks of the lt and lg variables passed from C to
!                       be against _c_int variables, which they should be.
!          Version 5.2: 2025-12-31
!                       Declared intent of threads variable in rdelap_f.
!
! LICENSE:
!   Copyright (c) 2016, Avraham Adler
!   All rights reserved.
!
!   Redistribution and use in source and binary forms, with or without
!   modification, are permitted provided that the following conditions are met:
!       1. Redistributions of source code must retain the above copyright
!          notice, this list of conditions and the following disclaimer.
!       2. Redistributions in binary form must reproduce the above copyright
!          notice, this list of conditions and the following disclaimer in the
!          documentation and/or other materials provided with the distribution.
!
!   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
!   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
!   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
!   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
!   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
!   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!   POSSIBILITY OF SUCH DAMAGE.
!-------------------------------------------------------------------------------

module delaporte
    use, intrinsic :: iso_c_binding
    use, intrinsic :: iso_fortran_env
    use, intrinsic :: ieee_arithmetic
    !$ use omp_lib
    use utils

    implicit none
    private
    public :: ddelap_f, pdelap_f, qdelap_f, rdelap_f, momdelap_f

contains

!-------------------------------------------------------------------------------
! FUNCTION: ddelap_f_s
!
! DESCRIPTION: Calculate the Delaporte probability mass function for a single
!              observation and return the value or its log. Calculated through
!              explicit summation. Follows R convention that real observations
!              are errors and have 0 probability, so calls floor to build to
!              the last integer. Implements hard floor of 0 and hard ceiling of
!              1 to prevent spurious floating point errors. The single vector
!              trick from pdelap actually slows ddelap down in almost all case
!              unless the passed vectors are very close to one another and small
!              in magnitude, so it is not worth programming for now.
!-------------------------------------------------------------------------------

    pure elemental function ddelap_f_s(x, alpha, beta, lambda) result(pmf)

    real(kind = c_double), intent(in)   :: x, alpha, beta, lambda
    real(kind = c_double)               :: pmf, ii, kk
    integer(INT64)                      :: i, k                   

        if (alpha <= ZERO .or. beta <= ZERO .or. lambda <= ZERO .or. x < ZERO &
            .or. ieee_is_nan(alpha + beta + lambda + x)) then
            pmf = ieee_value(x, ieee_quiet_nan)
        else
            pmf = ZERO
            k = floor(x, INT64)
            kk = real(k, c_double)
            if (x < MAXD .and. x == kk) then
                do i = 0_INT64, k
                    ii = real(i, c_double)
                    pmf = pmf + exp(log_gamma(alpha + ii) + ii * log(beta) &
                    + (kk - ii) * log(lambda) - lambda - log_gamma(alpha) &
                    - log_gamma(ii + ONE) - (alpha + ii) * log1p(beta) &
                    - log_gamma(kk - ii + ONE))
                end do
            pmf = cFPe(pmf)                       ! Clear floating point errors
            end if
        end if

    end function ddelap_f_s

!-------------------------------------------------------------------------------
! ROUTINE: ddelap_f
!
! DESCRIPTION: Vector-based PMF allowing parameter vector recycling and called 
!              from C. As Fortran starts its indices at 1, for the mod function
!              to properly recycle the vectors, the index needs to be reduced by
!              one, mod applied, and then increased by one again. Follows R
!              convention that real observations are errors and have 0
!              probability, so returns 0 for non-integer without calling
!              summation loop. 
!-------------------------------------------------------------------------------

    subroutine ddelap_f(x, nx, a, na, b, nb, l, nl, lg, threads, pmfv) &
               bind(C, name="ddelap_f_")
                        
    integer(kind = c_int), intent(in), value     :: nx, na, nb, nl
    real(kind = c_double), intent(in)            :: x(nx), a(na), b(nb), l(nl)
    integer(kind = c_int), intent(in)            :: lg, threads
    real(kind = c_double), intent(out)           :: pmfv(nx)
    integer                                      :: i
    
        !$omp parallel do num_threads(threads) default(shared) private(i) &
        !$omp schedule(static)
        do i = 1, nx
            pmfv(i) = ddelap_f_s(x(i), a(imk(i, na)), b(imk(i, nb)), &
            l(imk(i, nl)))
            if (lg == 1_c_int) pmfv(i) = log(pmfv(i))
        end do
        !$omp end parallel do
        
        if (any(ieee_is_nan(pmfv))) call rwarn("NaNs produced")

    end subroutine ddelap_f
    
!-------------------------------------------------------------------------------
! FUNCTION: pdelap_f_s
!
! DESCRIPTION: Calculate the Delaporte cumulative distribution function for a
!              single observation and return the value or its log. Calculated
!              through explicit summation. Follows R convention that real
!              observations are errors and have 0 probability, so calls floor to
!              build to last integer. Implements hard floor of 0 and hard
!              ceiling of 1 to prevent spurious floating point errors.
!-------------------------------------------------------------------------------

    pure elemental function pdelap_f_s(q, alpha, beta, lambda) result(cdf)

    real(kind = c_double), intent(in)   :: q, alpha, beta, lambda
    real(kind = c_double)               :: cdf
    integer(INT64)                      :: i, k

        if (alpha <= ZERO .or. beta <= ZERO .or. lambda <= ZERO .or. q < ZERO &
            .or. ieee_is_nan(alpha + beta + lambda + q)) then
            cdf = ieee_value(q, ieee_quiet_nan)
        else if (.not. ieee_is_finite(q)) then
            cdf = ONE
        else
            k = floor(q, INT64)
            cdf = exp(-lambda) / ((beta + ONE) ** alpha)
            do i = 1_INT64, k
                cdf = cdf + ddelap_f_s(real(i, c_double), alpha, beta, lambda)
            end do
            cdf = cFPe(cdf)                       ! Clear floating point errors
        end if

    end function pdelap_f_s

!-------------------------------------------------------------------------------
! ROUTINE: pdelap_f
!
! DESCRIPTION: Vector-based CDF allowing parameter vector recycling and called
!              from C. If parameters are all singletons (not vectors) then the
!              idea is to find the largest value in the vector and build the PDF
!              up to that point. Building the vector has each succesive value
!              piggyback off of the prior instead of calling p_delap_f_s each
!              time which increases the speed dramatically. Once created,
!              remaining values are simple lookups off of the svec vector.
!              Otherwise, each entry will need to build its own pmf value by
!              calling p_delap_f_s on each entry. Implements hard floor of 0 and
!              hard ceiling of 1 to prevent spurious floating point errors.
!-------------------------------------------------------------------------------

    subroutine pdelap_f(q, nq, a, na, b, nb, l, nl, lt, lg, threads, pmfv) &
               bind(C, name="pdelap_f_")
                        
    integer(kind = c_int), intent(in), value    :: nq, na, nb, nl
    real(kind = c_double), intent(in)           :: q(nq), a(na), b(nb), l(nl)
    integer(kind = c_int), intent(in)           :: lg, lt, threads
    real(kind = c_double), intent(out)          :: pmfv(nq)
    real(kind = c_double), allocatable          :: svec(:)
    integer                                     :: i, k

! If there are any complications at all, don't use the fast version. pdelap_f_s
! and ddelap_f_s are more robust to improper entries.

        if (na > 1 .or. nb > 1 .or. nl > 1 .or. minval(q) < ZERO .or. &
            maxval(q) > REAL(MAXVECSIZE, c_double) &
            .or. any(ieee_is_nan(q))) then
            !$omp parallel do num_threads(threads) default(shared) private(i) &
            !$omp schedule(static)
                do i = 1, nq
                    pmfv(i) = pdelap_f_s(q(i), a(imk(i, na)), b(imk(i, nb)), &
                    l(imk(i, nl)))
                    if (lt == 0_c_int) pmfv(i) = HALF - pmfv(i) + HALF!See dpq.h
                    if (lg == 1_c_int) pmfv(i) = log(pmfv(i))
                end do
            !$omp end parallel do
        else if (a(1) <= ZERO .or. b(1) <= ZERO .or. l(1) <= ZERO .or. &
                 ieee_is_nan(a(1) + b(1) + l(1))) then
            pmfv = ieee_value(q, ieee_quiet_nan)
        else
            k = floor(maxval(q))
            allocate (svec(k + 1))
            svec(1) = cFPe(exp(-l(1)) / ((b(1) + ONE) ** a(1)))
            do i = 2, k + 1
                svec(i) = cFPe(svec(i - 1) + &
                ddelap_f_s(real(i - 1, c_double), a(1), b(1), l(1)))
            end do
            do i = 1, nq
                pmfv(i) = svec(floor(q(i)) + 1)
                if (lt == 0_c_int) pmfv(i) = HALF - pmfv(i) + HALF  ! See dpq.h
                if (lg == 1_c_int) pmfv(i) = log(pmfv(i))
            end do
            deallocate(svec)
        end if
        
        if (any(ieee_is_nan(pmfv))) call rwarn("NaNs produced")
        
    end subroutine pdelap_f

!-------------------------------------------------------------------------------
! FUNCTION: qdelap_f_s
!
! DESCRIPTION: Calculate the Delaporte quantile function for a single 
!              observation and return the value. Calculated through explicit
!              summation. Returns NaN and Inf where appropriate.
!-------------------------------------------------------------------------------

    pure elemental function qdelap_f_s(p, alpha, beta, lambda) result(value)

    real(kind = c_double), intent(in)   :: p, alpha, beta, lambda
    real(kind = c_double)               :: testcdf, value

        if (alpha <= ZERO .or. beta <= ZERO .or. lambda <= ZERO .or. p < ZERO &
          .or. ieee_is_nan(alpha + beta + lambda + p)) then
            value = ieee_value(p, ieee_quiet_nan)
        else if (p >= ONE) then
            value = ieee_value(p, ieee_positive_inf)
        else
            value = ZERO
            testcdf = exp(-lambda) / ((beta + ONE) ** alpha)
            do while (p > testcdf)
                value = value + ONE
                testcdf = testcdf + ddelap_f_s(value, alpha, beta, lambda)
            end do
        end if

    end function qdelap_f_s

!-------------------------------------------------------------------------------
! ROUTINE: qdelap_f
!
! DESCRIPTION: Vector-based quantile function with parameter vector recycling.
!              If parameters are all singletons (not vectors) then the idea is
!              to find the largest value in the vector and build the PDF up to
!              that point. Building the vector has each succesive value
!              piggyback off of the prior instead of calling p_delap_f_s each
!              time which increases the speed dramatically. Once created,
!              remaining values are lookups off of the svec vector.
!              Otherwise, each entry will need to build its own pmf value by
!              calling q_delap_f_s on each entry.
!-------------------------------------------------------------------------------

    subroutine qdelap_f(p, np, a, na, b, nb, l, nl, lt, lg, threads, obsv) &
               bind(C, name="qdelap_f_")

    integer(kind = c_int), intent(in), value        :: np, na, nb, nl
    real(kind = c_double), intent(in)               :: a(na), b(nb), l(nl)
    integer(kind = c_int), intent(in)               :: lg, lt, threads
    real(kind = c_double), intent(inout)            :: p(np) 
    real(kind = c_double), intent(out)              :: obsv(np)
    real(kind = c_double), allocatable              :: svec(:), tvec(:)
    real(kind = c_double)                           :: x
    integer                                         :: i

        if (lg == 1_c_int) p = exp(p)

        if (lt == 0_c_int) p = HALF - p + HALF  ! See See dpq.h in R source code

        if(na == 1 .and. nb == na .and. nl == nb) then
            if (a(1) <= ZERO .or. b(1) <= ZERO .or. l(1) <= ZERO) then
                obsv = ieee_value(p, ieee_quiet_nan)
            else
                x = maxval(p, 1, p < 1)
                allocate(svec(1), source = exp(-l(1)) / ((b(1) + ONE) ** a(1)))
                i = 1
                do
                    if (svec(i) >= x) exit
                    i = i + 1
                    allocate(tvec(1:i), source = ZERO)
                    tvec(1:i-1) = svec
                    call move_alloc(tvec, svec)
                    svec(i) = svec(i - 1) + &
                        ddelap_f_s(real(i - 1, c_double), a(1), b(1), l(1))
                end do
                do i = 1, np
                    if (p(i) < ZERO .or. ieee_is_nan(p(i))) then
                        obsv(i) = ieee_value(p(i), ieee_quiet_nan)
                    else if (p(i) >= ONE) then
                        obsv(i) = ieee_value(p(i), ieee_positive_inf)
                    else
                        obsv(i) = real(minloc(svec, dim = 1, &
                            mask = svec >= p(i)) - 1)
                    end if
                end do
                deallocate(svec)
            end if
        else
            !$omp parallel do num_threads(threads) default(shared) private(i) &
            !$omp schedule(static)
            do i = 1, np
                obsv(i) = qdelap_f_s(p(i), a(imk(i, na)), b(imk(i, nb)), &
                          l(imk(i, nl)))
            end do
            !$omp end parallel do
        end if
        
    end subroutine qdelap_f

!-------------------------------------------------------------------------------
! ROUTINE: rdelap_f
!
! DESCRIPTION: Vector-based random number generator with parameter vector
!              recycling. It calls a C procedure to generate uniform random
!              variates which jibe with R's own internals and then calls
!              qdelap_f on the uniforms. This allows qdelap's singleton mode to
!              activate if appropriate. This is the single routine slower in
!              this Fortran implementation than the prior C++ implementation, as
!              the vector creation and pushback is more efficient in C++ STL
!              than the ballet between allocate and move_alloc in Fortran. On
!              vector-valued parameters Fortran is faster than C++. Technically
!              this is a slowdown in qdelap, not rdelap, but the C++ version of
!              qdelap did not use the vector lookup trick; it was only
!              programmed in rdelap, wheras now the Fortran version of qdelap
!              uses the trick for a net speedup. Only rdelap suffers slightly.
!-------------------------------------------------------------------------------

    subroutine rdelap_f(n, a, na, b, nb, l, nl, threads, vars) &
               bind(C, name="rdelap_f_")

    external unifrnd

    integer(kind = c_int), intent(in), value           :: n, na, nb, nl
    real(kind = c_double), intent(in)                  :: a(na), b(nb), l(nl)
    real(kind = c_double), intent(out)                 :: vars(n)
    integer(kind = c_int), intent(in)                  :: threads
    real(kind = c_double)                              :: p(n)

        call unifrnd(n, p)
        call qdelap_f(p, n, a, na, b, nb, l, nl, 1, 0, threads, vars)

    end subroutine rdelap_f

!-------------------------------------------------------------------------------
! ROUTINE: momdelap_f
!
! DESCRIPTION: Calculates method of moments estimates of parameters for a 
!              Delaporte distribution based on supplied vector. Based on
!              algorithms of Welford, Knuth, and Cook.
!              https://www.johndcook.com/blog/skewness_kurtosis/
!-------------------------------------------------------------------------------

    pure subroutine momdelap_f(obs, n, tp, params) bind(C, name="momdelap_f_")

    integer(kind = c_int), intent(in), value :: n
    integer(kind = c_int), intent(in)        :: tp
    real(kind = c_double), intent(in)        :: obs(n)
    real(kind = c_double), intent(out)       :: params(3)
    real(kind = c_double)                    :: nnm1, P, Mu_D, M2, M3, T1
    real(kind = c_double)                    :: delta, delta_i, nn, Var_D
    real(kind = c_double)                    :: Skew_D, VmM_D, ii
    integer                                  :: i

        nn = real(n, c_double)
        nnm1 = nn - ONE
        select case (tp)
            case (1_c_int)
                P = ONE
            case (2_c_int)
                P = sqrt(nn * nnm1) / (nn - TWO)
            case (3_c_int)
                P = (nnm1 / nn) ** THREEHALFS
            case default
                P = sqrt(nn * nnm1) / (nn - TWO)
        end select
        Mu_D = ZERO
        M2 = ZERO
        M3 = ZERO
        do i = 1, n
            ii = real(i, c_double)
            delta = obs(i) - Mu_D
            delta_i = delta / ii
            T1 = delta * delta_i * (ii - ONE)
            Mu_D = Mu_D + delta_i
            M3 = M3 + (T1 * delta_i * (ii - TWO) - THREE * delta_i * M2)
            M2 = M2 + T1
        end do
        Var_D = M2 / nnm1
        Skew_D = P * sqrt(nn) * M3 / (M2 ** THREEHALFS)
        VmM_D = Var_D - Mu_D
        params(2) = HALF * (Skew_D * (Var_D ** THREEHALFS) - Mu_D - THREE &
                            * VmM_D) / VmM_D
        params(1) = VmM_D / (params(2) ** 2)
        params(3) = Mu_D - params(1) * params(2)
 
    end subroutine momdelap_f

end module delaporte
