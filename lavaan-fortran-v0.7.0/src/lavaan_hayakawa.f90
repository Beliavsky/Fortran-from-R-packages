module lavaan_hayakawa
   use lavaan_kinds, only : dp
   use lavaan_linalg, only : vech
   implicit none
   private

   type, public :: hayakawa_trace_result
      real(dp) :: trace_ugamma = 0.0_dp
      real(dp) :: trace_ugamma2 = 0.0_dp
      real(dp) :: mv_df = 0.0_dp
      real(dp) :: mv_scaling = 1.0_dp
      real(dp) :: chisq_mv = huge(1.0_dp)
      real(dp) :: ss_df = 0.0_dp
      real(dp) :: ss_scaling = 1.0_dp
      real(dp) :: ss_shift = 0.0_dp
      real(dp) :: chisq_ss = huge(1.0_dp)
      integer :: nobs = 0
      integer :: nstat = 0
      integer :: status = 0
   end type hayakawa_trace_result

   public :: hayakawa_trace_corrected, hayakawa_adjusted_tests

contains

   subroutine hayakawa_trace_corrected(u, data, result, meanstructure)
      real(dp), intent(in) :: u(:, :), data(:, :)
      type(hayakawa_trace_result), intent(out) :: result
      logical, intent(in), optional :: meanstructure
      real(dp), allocatable :: yc(:, :), z(:, :), zc(:, :), gram(:, :), outer(:, :), vv(:)
      real(dp), allocatable :: mean(:), zmean(:)
      real(dp) :: tr_h, tr_h2, tr_d2, rn
      integer :: n, p, pstar, q, i
      logical :: with_mean

      n = size(data, 1)
      p = size(data, 2)
      with_mean = .false.
      if (present(meanstructure)) with_mean = meanstructure
      pstar = p * (p + 1) / 2
      q = pstar + merge(p, 0, with_mean)

      result%nobs = n
      result%nstat = q
      if (n < 4 .or. p < 1) then
         result%status = -1
         return
      end if
      if (size(u, 1) /= q .or. size(u, 2) /= q) then
         result%status = -2
         return
      end if

      allocate(mean(p), yc(n, p), z(n, q))
      mean = sum(data, dim=1) / real(n, dp)
      do i = 1, n
         yc(i, :) = data(i, :) - mean
         outer = spread(yc(i, :), 2, p) * spread(yc(i, :), 1, p)
         vv = vech(outer)
         if (with_mean) then
            z(i, 1:p) = data(i, :)
            z(i, p+1:q) = vv
         else
            z(i, :) = vv
         end if
      end do

      allocate(zmean(q), zc(n, q))
      zmean = sum(z, dim=1) / real(n, dp)
      do i = 1, n
         zc(i, :) = z(i, :) - zmean
      end do

      ! Gram matrix M = Zc U Zc'.  This is the form used by the current
      ! lavaan implementation and avoids constructing U**(1/2).
      gram = matmul(matmul(zc, u), transpose(zc))
      tr_h = 0.0_dp
      tr_h2 = sum(gram * gram)
      tr_d2 = 0.0_dp
      do i = 1, n
         tr_h = tr_h + gram(i, i)
         tr_d2 = tr_d2 + gram(i, i) * gram(i, i)
      end do

      rn = real(n, dp)
      result%trace_ugamma = tr_h / real(n - 1, dp)
      result%trace_ugamma2 = &
         (real((n - 2) * (n - 1), dp) * tr_h2 - &
          real(n * (n - 1), dp) * tr_d2 + tr_h * tr_h) / &
         (rn * real(n - 1, dp) * real(n - 2, dp) * real(n - 3, dp))
      result%status = 0
   end subroutine hayakawa_trace_corrected

   subroutine hayakawa_adjusted_tests(u,data,chisq,df,result,meanstructure)
      real(dp),intent(in)::u(:,:),data(:,:),chisq,df
      type(hayakawa_trace_result),intent(out)::result
      logical,intent(in),optional::meanstructure
      real(dp)::a
      if(present(meanstructure)) then
         call hayakawa_trace_corrected(u,data,result,meanstructure)
      else
         call hayakawa_trace_corrected(u,data,result)
      end if
      if(result%status/=0) return
      if(df<=0.0_dp .or. result%trace_ugamma<=0.0_dp .or. result%trace_ugamma2<=0.0_dp) then
         result%status=-3
         return
      end if
      result%mv_df=result%trace_ugamma*result%trace_ugamma/result%trace_ugamma2
      result%mv_scaling=result%trace_ugamma2/result%trace_ugamma
      result%chisq_mv=chisq/result%mv_scaling
      result%ss_df=df
      a=sqrt(df/result%trace_ugamma2)
      result%ss_scaling=1.0_dp/a
      result%ss_shift=df-a*result%trace_ugamma
      result%chisq_ss=chisq*a+result%ss_shift
   end subroutine hayakawa_adjusted_tests

end module lavaan_hayakawa
