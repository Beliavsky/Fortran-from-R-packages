module misc_tools_stats
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, ieee_quiet_nan
   use misc_tools_kinds, only : dp
   use misc_tools_special, only : student_t_two_sided_p
   implicit none
   private

   type, public :: coef_table_result
      real(dp), allocatable :: table(:,:)
   end type coef_table_result

   public :: ddnorm, median_value, col_medians, row_medians
   public :: col_medians_3d, col_medians_4d
   public :: r_squared, std_er, coef_table
   public :: n_obs_matrix, n_param_vector

contains

   pure elemental real(dp) function normal_pdf(x,mean,sd) result(p)
      real(dp), intent(in) :: x,mean,sd
      real(dp), parameter :: pi = acos(-1.0_dp)
      real(dp) :: z
      if (sd <= 0.0_dp) then
         p = 0.0_dp
      else
         z = (x-mean)/sd
         p = exp(-0.5_dp*z*z)/(sd*sqrt(2.0_dp*pi))
      end if
   end function normal_pdf

   pure elemental real(dp) function ddnorm(x,mean,sd) result(deriv)
      real(dp), intent(in) :: x
      real(dp), intent(in), optional :: mean,sd
      real(dp) :: mu,sigma
      mu = 0.0_dp
      sigma = 1.0_dp
      if (present(mean)) mu = mean
      if (present(sd)) sigma = sd
      if (sigma <= 0.0_dp) then
         deriv = 0.0_dp
      else
         deriv = -normal_pdf(x,mu,sigma)*(x-mu)/(sigma*sigma)
      end if
   end function ddnorm

   subroutine sort_real(x)
      real(dp), intent(inout) :: x(:)
      real(dp) :: v
      integer :: i,j
      do i = 2, size(x)
         v = x(i)
         j = i-1
         do while (j >= 1)
            if (x(j) <= v) exit
            x(j+1) = x(j)
            j = j-1
         end do
         x(j+1) = v
      end do
   end subroutine sort_real

   real(dp) function median_value(x,na_rm,status) result(med)
      real(dp), intent(in) :: x(:)
      logical, intent(in), optional :: na_rm
      integer, intent(out), optional :: status
      logical :: rm
      real(dp), allocatable :: y(:)
      integer :: n

      if (present(status)) status = 0
      rm = .false.
      if (present(na_rm)) rm = na_rm

      if (rm) then
         y = pack(x,.not. ieee_is_nan(x))
      else
         if (any(ieee_is_nan(x))) then
            med = ieee_value(0.0_dp,ieee_quiet_nan)
            if (present(status)) status = 2
            return
         end if
         y = x
      end if
      n = size(y)
      if (n == 0) then
         med = ieee_value(0.0_dp,ieee_quiet_nan)
         if (present(status)) status = 1
         return
      end if
      call sort_real(y)
      if (mod(n,2) == 1) then
         med = y((n+1)/2)
      else
         med = 0.5_dp*(y(n/2)+y(n/2+1))
      end if
   end function median_value

   subroutine col_medians(x,result,na_rm)
      real(dp), intent(in) :: x(:,:)
      real(dp), allocatable, intent(out) :: result(:)
      logical, intent(in), optional :: na_rm
      integer :: j
      allocate(result(size(x,2)))
      do j = 1, size(x,2)
         result(j) = median_value(x(:,j),na_rm)
      end do
   end subroutine col_medians

   subroutine row_medians(x,result,na_rm)
      real(dp), intent(in) :: x(:,:)
      real(dp), allocatable, intent(out) :: result(:)
      logical, intent(in), optional :: na_rm
      integer :: i
      allocate(result(size(x,1)))
      do i = 1, size(x,1)
         result(i) = median_value(x(i,:),na_rm)
      end do
   end subroutine row_medians

   subroutine col_medians_3d(x,result,na_rm)
      real(dp), intent(in) :: x(:,:,:)
      real(dp), allocatable, intent(out) :: result(:,:)
      logical, intent(in), optional :: na_rm
      integer :: j,k
      allocate(result(size(x,2),size(x,3)))
      do k = 1, size(x,3)
         do j = 1, size(x,2)
            result(j,k) = median_value(x(:,j,k),na_rm)
         end do
      end do
   end subroutine col_medians_3d

   subroutine col_medians_4d(x,result,na_rm)
      real(dp), intent(in) :: x(:,:,:,:)
      real(dp), allocatable, intent(out) :: result(:,:,:)
      logical, intent(in), optional :: na_rm
      integer :: j,k,l
      allocate(result(size(x,2),size(x,3),size(x,4)))
      do l = 1, size(x,4)
         do k = 1, size(x,3)
            do j = 1, size(x,2)
               result(j,k,l) = median_value(x(:,j,k,l),na_rm)
            end do
         end do
      end do
   end subroutine col_medians_4d

   real(dp) function r_squared(y,resid,status) result(r2)
      real(dp), intent(in) :: y(:),resid(:)
      integer, intent(out), optional :: status
      real(dp) :: mean_y,tss,rss
      if (present(status)) status = 0
      if (size(y) /= size(resid) .or. size(y) == 0) then
         r2 = 0.0_dp
         if (present(status)) status = 1
         return
      end if
      mean_y = sum(y)/real(size(y),dp)
      tss = sum((y-mean_y)**2)
      rss = sum(resid*resid)
      if (tss <= tiny(1.0_dp)) then
         r2 = 0.0_dp
         if (present(status)) status = 2
      else
         r2 = 1.0_dp-rss/tss
      end if
   end function r_squared

   subroutine std_er(cov,se,status)
      real(dp), intent(in) :: cov(:,:)
      real(dp), allocatable, intent(out) :: se(:)
      integer, intent(out), optional :: status
      integer :: i,n
      if (present(status)) status = 0
      if (size(cov,1) /= size(cov,2)) then
         allocate(se(0))
         if (present(status)) status = 1
         return
      end if
      n = size(cov,1)
      allocate(se(n))
      do i = 1, n
         if (cov(i,i) < 0.0_dp) then
            se(i) = 0.0_dp
            if (present(status)) status = 2
         else
            se(i) = sqrt(cov(i,i))
         end if
      end do
   end subroutine std_er

   function coef_table(coef,std_err,df) result(res)
      real(dp), intent(in) :: coef(:),std_err(:)
      real(dp), intent(in), optional :: df
      type(coef_table_result) :: res
      integer :: i,n
      real(dp) :: t

      if (size(coef) /= size(std_err)) error stop "coef_table: length mismatch"
      n = size(coef)
      allocate(res%table(n,4))
      res%table(:,1) = coef
      res%table(:,2) = std_err
      do i = 1, n
         if (std_err(i) <= 0.0_dp) then
            res%table(i,3) = huge(1.0_dp)
         else
            res%table(i,3) = coef(i)/std_err(i)
         end if
         if (present(df)) then
            t = abs(res%table(i,3))
            res%table(i,4) = student_t_two_sided_p(t,df)
         else
            res%table(i,4) = -1.0_dp
         end if
      end do
   end function coef_table

   pure integer function n_obs_matrix(x) result(n)
      real(dp), intent(in) :: x(:,:)
      n = size(x,1)
   end function n_obs_matrix

   pure integer function n_param_vector(coef) result(n)
      real(dp), intent(in) :: coef(:)
      n = size(coef)
   end function n_param_vector

end module misc_tools_stats
