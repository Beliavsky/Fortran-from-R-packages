!-------------------------------------------------------------------------------
!
! MODULE: Utils
!
! AUTHOR: Avraham Adler <Avraham.Adler@gmail.com>
!
! DESCRIPTION: Utility functions and definitions for Delaporte package
!
! HISTORY:
!          Version 1.0: 2016-11-20
!          Version 1.1: 2017-03-01
!          Version 1.2: 2017-11-20
!                       Reformatted to 80 columns
!          Version 1.3: 2018-11-20
!                       Added cleanzeros function to handle EPS issues for right
!                       tail. See Issue #1
!          Version 1.4: 2018-12-10
!                       Cleanzeros removed in favor of floor & ceiling of 0 & 1
!          Version 1.5: 2021-01-03
!                       Added parameter for pdelap max vector size
!          Version 2.0: 2023-01-29
!                       Updated to rely on Fortran 2008 intrinsics
!          Version 2.1: 2023-08-08
!                       Added OpenMP control functions
!          Version 3.0: 2023-09-28
!                       Converted log1p to one based on its Taylor expansion.
!                       This is also the degree 2 polynomial minimax
!                       approximation.
!          Version 3.1: 2024-05-21
!                       Added pure header to log1p.
!          Version 4.0: 2024-06-17
!                       Added imk helper function. A smidgen faster—I'm not sure
!                       why, perhaps due to pre-compilation in module—and easier
!                       to read. Turn FP error cleaning into a function.
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

module utils
    use, intrinsic :: iso_c_binding
    use, intrinsic :: iso_fortran_env
    !$ use omp_lib
    implicit none

    real(kind = c_double), parameter :: ZERO = 0._c_double
    real(kind = c_double), parameter :: HALF = 0.5_c_double
    real(kind = c_double), parameter :: ONE = 1._c_double
    real(kind = c_double), parameter :: THREEHALFS = 1.5_c_double
    real(kind = c_double), parameter :: TWO = 2._c_double
    real(kind = c_double), parameter :: THREE = 3._c_double
    real(kind = c_double), parameter :: EPS = 2.2204460492503131e-16_c_double
    real(kind = c_double), parameter :: MAXD = REAL(HUGE(1_INT64), c_double)
    integer, parameter               :: MAXVECSIZE = 16384

contains

!-------------------------------------------------------------------------------
! FUNCTION: log1p
!
! DESCRIPTION: Fortran 2008 does not have log1p as an intrinsic. This uses the
!              Taylor expansion for small x to reduce relative error.
!-------------------------------------------------------------------------------

    pure elemental function log1p(x) result(y)

        real(kind = c_double), intent(in) :: x
        real(kind = c_double)             :: y

        if (abs(x) <= 1.e-4_c_double) then
            y = (-x * HALF + ONE) * x
        else
            y = log(x + ONE)
        end if
            
    end function log1p
    
!-------------------------------------------------------------------------------
! FUNCTION: gOMPT
!
! DESCRIPTION: Gets current OMP Threads
!-------------------------------------------------------------------------------

    subroutine gOMPT_f(n) bind(C, name="gOMPT_f_")
    
    integer(kind = c_int), intent(out) :: n
    
        n = 1_c_int
        !$ n = omp_get_max_threads()
    
    end subroutine gOMPT_f
    
    
!-------------------------------------------------------------------------------
! FUNCTION: sOMPT
!
! DESCRIPTION: Sets current OMP Threads
!-------------------------------------------------------------------------------

    subroutine sOMPT_f(n) bind(C, name="sOMPT_f_")
    
    integer(kind = c_int), intent(in) :: n
    
        !$ call omp_set_num_threads(n)
    
    end subroutine sOMPT_f    

!-------------------------------------------------------------------------------
! FUNCTION: imk (i mod k)
!
! DESCRIPTION: Calculates mod(i - 1, k) + 1 for vector recyling.
!-------------------------------------------------------------------------------

    pure elemental function imk(i, k) result(j)

    integer(kind = c_int), intent(in) :: i, k
    integer(kind = c_int)             :: j
    
        j = mod(i - 1, k) + 1
    
    end function imk
    
!-------------------------------------------------------------------------------
! FUNCTION: cFPe (clearFPerrors)
!
! DESCRIPTION: Restricts solutions to [0, 1] and eliminates spurious FP errors.
!-------------------------------------------------------------------------------

    pure elemental function cFPe(x) result(y)

    real(kind = c_double), intent(in) :: x
    real(kind = c_double)             :: y
    
        y = max(min(x, ONE), ZERO)
    
    end function cFPe
    
end module utils ! # nocov covr doesn't always pick up the end module
