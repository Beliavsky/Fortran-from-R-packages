! SPDX-License-Identifier: GPL-3.0-or-later
module cec2005benchmark
   use, intrinsic :: iso_fortran_env, only : real64
   implicit none
   private

   integer, parameter, public :: dp = real64
   integer, parameter :: nfunc = 10
   real(dp), parameter :: pi = 3.1415926535897932384626433832795029_dp
   real(dp), parameter :: euler_e = 2.7182818284590452353602874713526625_dp
   real(dp), parameter :: comp_scale = 2000.0_dp

   type, public :: cec2005_context
      integer :: function_id = 0
      integer :: n = 0
      logical :: initialized = .false.
      logical :: noise_enabled = .true.
      character(len=:), allocatable :: data_dir
      real(dp) :: global_bias = 0.0_dp
      real(dp), allocatable :: o(:,:), g(:,:), l(:,:,:)
      real(dp) :: sigma(nfunc) = 1.0_dp
      real(dp) :: lambda(nfunc) = 1.0_dp
      real(dp) :: bias(nfunc) = 0.0_dp
      real(dp) :: norm_f(nfunc) = 0.0_dp
      real(dp), allocatable :: a5(:,:), b5(:)
      real(dp), allocatable :: a12(:,:), b12(:,:), alpha12(:)
   contains
      procedure :: init => cec_init
      procedure :: evaluate => cec_evaluate
      procedure :: evaluate_batch => cec_evaluate_batch
      procedure :: set_noise => cec_set_noise
   end type cec2005_context

   public :: cec2005_eval, cec2005_eval_batch, cec2005_seed

contains

   subroutine cec2005_seed(seed)
      integer, intent(in) :: seed
      integer :: nseed, i
      integer, allocatable :: put(:)
      call random_seed(size=nseed)
      allocate(put(nseed))
      do i = 1, nseed
         put(i) = modulo(seed + 104729*i, huge(1)-1)
         if (put(i) <= 0) put(i) = i
      end do
      call random_seed(put=put)
   end subroutine cec2005_seed

   subroutine cec_set_noise(self, enabled)
      class(cec2005_context), intent(inout) :: self
      logical, intent(in) :: enabled
      self%noise_enabled = enabled
      if (self%initialized .and. self%function_id >= 15) call compute_norms(self)
   end subroutine cec_set_noise

   subroutine cec_init(self, function_id, dimension, data_dir, noise_enabled, ierr)
      class(cec2005_context), intent(inout) :: self
      integer, intent(in) :: function_id, dimension
      character(len=*), intent(in), optional :: data_dir
      logical, intent(in), optional :: noise_enabled
      integer, intent(out), optional :: ierr
      integer :: i, stat

      stat = 0
      if (function_id < 1 .or. function_id > 25) stat = 1
      if (.not. supported_dimension(dimension)) stat = 2
      if (present(ierr)) ierr = stat
      if (stat /= 0) return

      self%function_id = function_id
      self%n = dimension
      self%initialized = .false.
      if (present(data_dir)) then
         self%data_dir = trim(data_dir)
      else
         self%data_dir = 'data'
      end if
      self%noise_enabled = .true.
      if (present(noise_enabled)) self%noise_enabled = noise_enabled

      if (allocated(self%o)) deallocate(self%o, self%g, self%l)
      if (allocated(self%a5)) deallocate(self%a5, self%b5)
      if (allocated(self%a12)) deallocate(self%a12, self%b12, self%alpha12)
      allocate(self%o(nfunc,dimension), self%g(dimension,dimension), &
               self%l(nfunc,dimension,dimension))
      self%o = 0.0_dp
      self%g = 0.0_dp
      self%l = 0.0_dp
      do i = 1, dimension
         self%g(i,i) = 1.0_dp
         self%l(:,i,i) = 1.0_dp
      end do
      self%sigma = 1.0_dp
      self%lambda = 1.0_dp
      do i = 1, nfunc
         self%bias(i) = 100.0_dp * real(i-1,dp)
      end do
      self%norm_f = 0.0_dp
      self%global_bias = 0.0_dp

      call initialize_function(self, stat)
      if (present(ierr)) ierr = stat
      if (stat /= 0) return
      if (function_id >= 15) call compute_norms(self)
      self%initialized = .true.
   end subroutine cec_init

   logical function supported_dimension(n) result(ok)
      integer, intent(in) :: n
      ok = n == 2 .or. n == 10 .or. n == 30 .or. n == 50
   end function supported_dimension

   function path_join(dir, name) result(path)
      character(len=*), intent(in) :: dir, name
      character(len=:), allocatable :: path
      integer :: m
      m = len_trim(dir)
      if (m > 0 .and. (dir(m:m) == '/' .or. dir(m:m) == achar(92))) then
         path = trim(dir)//trim(name)
      else
         path = trim(dir)//'/'//trim(name)
      end if
   end function path_join

   subroutine open_data(self, name, unit, ierr)
      class(cec2005_context), intent(in) :: self
      character(len=*), intent(in) :: name
      integer, intent(out) :: unit
      integer, intent(out) :: ierr
      character(len=:), allocatable :: path
      path = path_join(self%data_dir, name)
      open(newunit=unit, file=path, status='old', action='read', iostat=ierr)
   end subroutine open_data

   function dim_file(stem, n) result(name)
      character(len=*), intent(in) :: stem
      integer, intent(in) :: n
      character(len=:), allocatable :: name
      character(len=16) :: sn
      write(sn,'(i0)') n
      name = trim(stem)//trim(sn)//'.txt'
   end function dim_file

   subroutine read_o(self, name, ierr)
      class(cec2005_context), intent(inout) :: self
      character(len=*), intent(in) :: name
      integer, intent(out) :: ierr
      integer :: u
      call open_data(self, name, u, ierr)
      if (ierr /= 0) return
      read(u,*,iostat=ierr) self%o(1,:)
      close(u)
   end subroutine read_o

   subroutine read_o_rows(self, name, ierr)
      class(cec2005_context), intent(inout) :: self
      character(len=*), intent(in) :: name
      integer, intent(out) :: ierr
      integer :: u, i
      real(dp) :: row(100)
      call open_data(self, name, u, ierr)
      if (ierr /= 0) return
      do i = 1, nfunc
         read(u,*,iostat=ierr) row
         if (ierr /= 0) exit
         self%o(i,:) = row(1:self%n)
      end do
      close(u)
   end subroutine read_o_rows

   subroutine read_g(self, name, ierr)
      class(cec2005_context), intent(inout) :: self
      character(len=*), intent(in) :: name
      integer, intent(out) :: ierr
      integer :: u, i
      call open_data(self, name, u, ierr)
      if (ierr /= 0) return
      do i = 1, self%n
         read(u,*,iostat=ierr) self%g(i,:)
         if (ierr /= 0) exit
      end do
      close(u)
   end subroutine read_g

   subroutine read_l(self, name, ierr)
      class(cec2005_context), intent(inout) :: self
      character(len=*), intent(in) :: name
      integer, intent(out) :: ierr
      integer :: u, i, j
      call open_data(self, name, u, ierr)
      if (ierr /= 0) return
      do i = 1, nfunc
         do j = 1, self%n
            read(u,*,iostat=ierr) self%l(i,j,:)
            if (ierr /= 0) exit
         end do
         if (ierr /= 0) exit
      end do
      close(u)
   end subroutine read_l

   subroutine initialize_function(self, ierr)
      class(cec2005_context), intent(inout) :: self
      integer, intent(out) :: ierr
      integer :: u, i, idx
      real(dp), allocatable :: tmpa(:,:), tmpb(:,:), tmpalpha(:)
      character(len=:), allocatable :: fn
      ierr = 0

      select case(self%function_id)
      case(1)
         call read_o(self, 'sphere_func_data.txt', ierr); self%bias(1) = -450.0_dp
      case(2,4)
         call read_o(self, 'schwefel_102_data.txt', ierr); self%bias(1) = -450.0_dp
      case(3)
         fn = dim_file('elliptic_M_D', self%n); call read_g(self, fn, ierr)
         if (ierr == 0) call read_o(self, 'high_cond_elliptic_rot_data.txt', ierr)
         self%bias(1) = -450.0_dp
      case(5)
         allocate(self%a5(self%n,self%n), self%b5(self%n))
         call open_data(self, 'schwefel_206_data.txt', u, ierr)
         if (ierr /= 0) return
         block
            real(dp) :: row(100)
            do i = 1, nfunc
               read(u,*,iostat=ierr) row
               if (ierr /= 0) exit
               self%o(i,:) = row(1:self%n)
            end do
            if (ierr == 0) then
               do i = 1, self%n
                  read(u,*,iostat=ierr) row
                  if (ierr /= 0) exit
                  self%a5(i,:) = row(1:self%n)
               end do
            end if
         end block
         close(u)
         if (ierr /= 0) return
         idx = (self%n + 3)/4
         self%o(1,1:idx) = -100.0_dp
         idx = (3*self%n)/4
         self%o(1,max(idx,1):self%n) = 100.0_dp
         self%b5 = matmul(self%a5, self%o(1,:))
         self%bias(1) = -310.0_dp
      case(6)
         call read_o(self, 'rosenbrock_func_data.txt', ierr)
         if (ierr == 0) self%o = self%o - 1.0_dp
         self%bias(1) = 390.0_dp
      case(7)
         fn = dim_file('griewank_M_D', self%n); call read_g(self, fn, ierr)
         if (ierr == 0) call read_o(self, 'griewank_func_data.txt', ierr)
         self%bias(1) = -180.0_dp
      case(8)
         fn = dim_file('ackley_M_D', self%n); call read_g(self, fn, ierr)
         if (ierr == 0) call read_o(self, 'ackley_func_data.txt', ierr)
         if (ierr == 0) then
            do i = 1, self%n/2
               self%o(1,2*i-1) = -32.0_dp
            end do
         end if
         self%bias(1) = -140.0_dp
      case(9)
         call read_o(self, 'rastrigin_func_data.txt', ierr); self%bias(1) = -330.0_dp
      case(10)
         fn = dim_file('rastrigin_M_D', self%n); call read_g(self, fn, ierr)
         if (ierr == 0) call read_o(self, 'rastrigin_func_data.txt', ierr)
         self%bias(1) = -330.0_dp
      case(11)
         fn = dim_file('weierstrass_M_D', self%n); call read_g(self, fn, ierr)
         if (ierr == 0) call read_o(self, 'weierstrass_data.txt', ierr)
         self%bias(1) = 90.0_dp
      case(12)
         allocate(self%a12(self%n,self%n), self%b12(self%n,self%n), self%alpha12(self%n))
         allocate(tmpa(100,100), tmpb(100,100), tmpalpha(100))
         call open_data(self, 'schwefel_213_data.txt', u, ierr)
         if (ierr /= 0) return
         do i = 1, 100
            read(u,*,iostat=ierr) tmpa(i,:)
            if (ierr /= 0) exit
         end do
         if (ierr == 0) then
            do i = 1, 100
               read(u,*,iostat=ierr) tmpb(i,:)
               if (ierr /= 0) exit
            end do
         end if
         if (ierr == 0) read(u,*,iostat=ierr) tmpalpha
         close(u)
         if (ierr /= 0) return
         self%a12 = tmpa(1:self%n,1:self%n)
         self%b12 = tmpb(1:self%n,1:self%n)
         self%alpha12 = tmpalpha(1:self%n)
         self%bias(1) = -460.0_dp
      case(13)
         call read_o(self, 'EF8F2_func_data.txt', ierr)
         if (ierr == 0) self%o = self%o - 1.0_dp
         self%bias(1) = -130.0_dp
      case(14)
         fn = dim_file('E_ScafferF6_M_D', self%n); call read_g(self, fn, ierr)
         if (ierr == 0) call read_o(self, 'E_ScafferF6_func_data.txt', ierr)
         self%bias(1) = -300.0_dp
      case(15)
         call read_o_rows(self, 'hybrid_func1_data.txt', ierr)
         self%lambda = [1.0_dp,1.0_dp,10.0_dp,10.0_dp,1.0_dp/12.0_dp,1.0_dp/12.0_dp, &
                        5.0_dp/32.0_dp,5.0_dp/32.0_dp,1.0_dp/20.0_dp,1.0_dp/20.0_dp]
         self%global_bias = 120.0_dp
      case(16,17)
         call read_o_rows(self, 'hybrid_func1_data.txt', ierr)
         if (ierr == 0) then
            fn = dim_file('hybrid_func1_M_D',self%n); call read_l(self,fn,ierr)
         end if
         self%lambda = [1.0_dp,1.0_dp,10.0_dp,10.0_dp,1.0_dp/12.0_dp,1.0_dp/12.0_dp, &
                        5.0_dp/32.0_dp,5.0_dp/32.0_dp,1.0_dp/20.0_dp,1.0_dp/20.0_dp]
         self%global_bias = 120.0_dp
      case(18,19,20)
         call read_o_rows(self, 'hybrid_func2_data.txt', ierr)
         if (ierr == 0 .and. self%function_id == 20) then
            do i = 1, self%n/2
               self%o(1,2*i) = 5.0_dp
            end do
         end if
         if (ierr == 0) then
            fn = dim_file('hybrid_func2_M_D',self%n); call read_l(self,fn,ierr)
         end if
         if (ierr == 0) self%o(nfunc,:) = 0.0_dp
         self%sigma = [1.0_dp,2.0_dp,1.5_dp,1.5_dp,1.0_dp,1.0_dp,1.5_dp,1.5_dp,2.0_dp,2.0_dp]
         self%lambda = [5.0_dp/16.0_dp,5.0_dp/32.0_dp,2.0_dp,1.0_dp,1.0_dp/10.0_dp, &
                        1.0_dp/20.0_dp,20.0_dp,10.0_dp,1.0_dp/6.0_dp,1.0_dp/12.0_dp]
         if (self%function_id == 19) then
            self%sigma(1) = 0.1_dp
            self%lambda(1) = 0.5_dp/32.0_dp
         end if
         self%global_bias = 10.0_dp
      case(21,22,23)
         call read_o_rows(self, 'hybrid_func3_data.txt', ierr)
         if (ierr == 0) then
            if (self%function_id == 22) then
               fn = dim_file('hybrid_func3_HM_D',self%n)
            else
               fn = dim_file('hybrid_func3_M_D',self%n)
            end if
            call read_l(self,fn,ierr)
         end if
         self%sigma = [1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,2.0_dp,2.0_dp,2.0_dp,2.0_dp,2.0_dp]
         self%lambda = [1.0_dp/4.0_dp,1.0_dp/20.0_dp,5.0_dp,1.0_dp,5.0_dp,1.0_dp, &
                        50.0_dp,10.0_dp,1.0_dp/8.0_dp,1.0_dp/40.0_dp]
         self%global_bias = 360.0_dp
      case(24,25)
         call read_o_rows(self, 'hybrid_func4_data.txt', ierr)
         if (ierr == 0) then
            fn = dim_file('hybrid_func4_M_D',self%n); call read_l(self,fn,ierr)
         end if
         self%sigma = 2.0_dp
         self%lambda = [10.0_dp,1.0_dp/4.0_dp,1.0_dp,5.0_dp/32.0_dp,1.0_dp, &
                        1.0_dp/20.0_dp,1.0_dp/10.0_dp,1.0_dp,1.0_dp/20.0_dp,1.0_dp/20.0_dp]
         self%global_bias = 260.0_dp
      end select
   end subroutine initialize_function

   subroutine transform(self, x, count, z)
      class(cec2005_context), intent(in) :: self
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: count
      real(dp), intent(out) :: z(:)
      real(dp) :: t2(self%n), t3(self%n)
      integer :: i, j
      t2 = (x - self%o(count,:))/self%lambda(count)
      do j = 1, self%n
         t3(j) = 0.0_dp
         do i = 1, self%n
            t3(j) = t3(j) + self%g(i,j)*t2(i)
         end do
      end do
      do j = 1, self%n
         z(j) = 0.0_dp
         do i = 1, self%n
            z(j) = z(j) + self%l(count,i,j)*t3(i)
         end do
      end do
   end subroutine transform

   subroutine transform_norm(self, count, z)
      class(cec2005_context), intent(in) :: self
      integer, intent(in) :: count
      real(dp), intent(out) :: z(:)
      real(dp) :: t2(self%n), t3(self%n)
      integer :: i, j
      t2 = 5.0_dp/self%lambda(count)
      do j = 1, self%n
         t3(j) = 0.0_dp
         do i = 1, self%n
            t3(j) = t3(j) + self%g(i,j)*t2(i)
         end do
      end do
      do j = 1, self%n
         z(j) = 0.0_dp
         do i = 1, self%n
            z(j) = z(j) + self%l(count,i,j)*t3(i)
         end do
      end do
   end subroutine transform_norm

   subroutine calc_weights(self, x, w)
      class(cec2005_context), intent(in) :: self
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: w(nfunc)
      real(dp) :: s, wmax
      integer :: i
      do i = 1, nfunc
         s = sum((x-self%o(i,:))**2)
         w(i) = exp(-s/(2.0_dp*real(self%n,dp)*self%sigma(i)**2))
      end do
      wmax = maxval(w)
      s = 0.0_dp
      do i = 1, nfunc
         if (abs(w(i)-wmax) > tiny(1.0_dp)) w(i) = w(i)*(1.0_dp-wmax**10)
         s = s + w(i)
      end do
      if (abs(s) <= tiny(1.0_dp)) then
         w = 1.0_dp/real(nfunc,dp)
      else
         w = w/s
      end if
   end subroutine calc_weights

   real(dp) function normal_deviate() result(z)
      real(dp) :: u1, u2
      call random_number(u1); call random_number(u2)
      u1 = max(u1, tiny(1.0_dp))
      z = sqrt(-2.0_dp*log(u1))*cos(2.0_dp*pi*u2)
   end function normal_deviate

   real(dp) function sphere(x) result(f)
      real(dp), intent(in) :: x(:)
      f = sum(x*x)
   end function sphere

   real(dp) function schwefel_102(x) result(f)
      real(dp), intent(in) :: x(:)
      real(dp) :: s
      integer :: i
      f = 0.0_dp; s = 0.0_dp
      do i = 1, size(x)
         s = s + x(i); f = f + s*s
      end do
   end function schwefel_102

   real(dp) function rosenbrock(x) result(f)
      real(dp), intent(in) :: x(:)
      integer :: i
      f = 0.0_dp
      do i = 1, size(x)-1
         f = f + 100.0_dp*(x(i)*x(i)-x(i+1))**2 + (x(i)-1.0_dp)**2
      end do
   end function rosenbrock

   real(dp) function griewank(x) result(f)
      real(dp), intent(in) :: x(:)
      real(dp) :: p
      integer :: i
      p = 1.0_dp
      do i = 1, size(x)
         p = p*cos(x(i)/sqrt(real(i,dp)))
      end do
      f = 1.0_dp + sum(x*x)/4000.0_dp - p
   end function griewank

   real(dp) function ackley(x) result(f)
      real(dp), intent(in) :: x(:)
      real(dp) :: s1, s2, rn
      rn = real(size(x),dp)
      s1 = -0.2_dp*sqrt(sum(x*x)/rn)
      s2 = sum(cos(2.0_dp*pi*x))/rn
      f = 20.0_dp + euler_e - 20.0_dp*exp(s1) - exp(s2)
   end function ackley

   real(dp) function rastrigin(x) result(f)
      real(dp), intent(in) :: x(:)
      f = sum(x*x - 10.0_dp*cos(2.0_dp*pi*x) + 10.0_dp)
   end function rastrigin

   real(dp) function weierstrass(x) result(f)
      real(dp), intent(in) :: x(:)
      integer :: i, k
      f = 0.0_dp
      do i = 1, size(x)
         do k = 0, 20
            f = f + 0.5_dp**k*cos(2.0_dp*pi*3.0_dp**k*(x(i)+0.5_dp))
         end do
      end do
   end function weierstrass

   real(dp) function elliptic(x) result(f)
      real(dp), intent(in) :: x(:)
      integer :: i, n
      n = size(x); f = 0.0_dp
      do i = 1, n
         f = f + x(i)*x(i)*1.0e6_dp**(real(i-1,dp)/real(n-1,dp))
      end do
   end function elliptic

   real(dp) function ef8f2(x) result(f)
      real(dp), intent(in) :: x(:)
      integer :: i, j, n
      real(dp) :: t
      n = size(x); f = 0.0_dp
      do i = 1, n
         j = merge(i+1,1,i<n)
         t = 100.0_dp*(x(i)*x(i)-x(j))**2 + (x(i)-1.0_dp)**2
         f = f + t*t/4000.0_dp - cos(t) + 1.0_dp
      end do
   end function ef8f2

   real(dp) function escaffer(x) result(f)
      real(dp), intent(in) :: x(:)
      integer :: i, j, n
      real(dp) :: r2, t1, t2
      n = size(x); f = 0.0_dp
      do i = 1, n
         j = merge(i+1,1,i<n)
         r2 = x(i)*x(i)+x(j)*x(j)
         t1 = sin(sqrt(r2))**2
         t2 = 1.0_dp + 0.001_dp*r2
         f = f + 0.5_dp + (t1-0.5_dp)/(t2*t2)
      end do
   end function escaffer

   real(dp) function nc_schaffer_pair(x, y) result(f)
      real(dp), intent(in) :: x, y
      real(dp) :: t(2), r(2), q, b, r2
      integer :: i, a
      t = [x,y]
      do i = 1,2
         if (abs(t(i)) >= 0.5_dp) then
            q = 2.0_dp*t(i); a = int(q); b = abs(q-real(a,dp))
            if (b < 0.5_dp) then
               r(i) = real(a,dp)/2.0_dp
            else if (q <= 0.0_dp) then
               r(i) = real(a-1,dp)/2.0_dp
            else
               r(i) = real(a+1,dp)/2.0_dp
            end if
         else
            r(i) = t(i)
         end if
      end do
      r2 = r(1)*r(1)+r(2)*r(2)
      f = 0.5_dp + (sin(sqrt(r2))**2-0.5_dp)/(1.0_dp+0.001_dp*r2)**2
   end function nc_schaffer_pair

   real(dp) function nc_rastrigin_fn(x) result(f)
      real(dp), intent(in) :: x(:)
      real(dp) :: y(size(x)), q, b
      integer :: i, a
      do i = 1, size(x)
         if (abs(x(i)) >= 0.5_dp) then
            q = 2.0_dp*x(i); a = int(q); b = abs(q-real(a,dp))
            if (b < 0.5_dp) then
               y(i) = real(a,dp)/2.0_dp
            else if (q <= 0.0_dp) then
               y(i) = real(a-1,dp)/2.0_dp
            else
               y(i) = real(a+1,dp)/2.0_dp
            end if
         else
            y(i) = x(i)
         end if
      end do
      f = rastrigin(y)
   end function nc_rastrigin_fn

   subroutine compute_norms(self)
      class(cec2005_context), intent(inout) :: self
      real(dp) :: z(self%n), zero(self%n)
      integer :: i
      zero = 0.0_dp
      select case(self%function_id)
      case(15:17)
         do i=1,2; call transform_norm(self,i,z); self%norm_f(i)=rastrigin(z); end do
         do i=3,4; call transform_norm(self,i,z); self%norm_f(i)=weierstrass(z)-weierstrass(zero); end do
         do i=5,6; call transform_norm(self,i,z); self%norm_f(i)=griewank(z); end do
         do i=7,8; call transform_norm(self,i,z); self%norm_f(i)=ackley(z); end do
         do i=9,10; call transform_norm(self,i,z); self%norm_f(i)=sphere(z); end do
      case(18:20)
         do i=1,2; call transform_norm(self,i,z); self%norm_f(i)=ackley(z); end do
         do i=3,4; call transform_norm(self,i,z); self%norm_f(i)=rastrigin(z); end do
         do i=5,6; call transform_norm(self,i,z); self%norm_f(i)=sphere(z); end do
         do i=7,8; call transform_norm(self,i,z); self%norm_f(i)=weierstrass(z)-weierstrass(zero); end do
         do i=9,10; call transform_norm(self,i,z); self%norm_f(i)=griewank(z); end do
      case(21:23)
         do i=1,2; call transform_norm(self,i,z); self%norm_f(i)=escaffer(z); end do
         do i=3,4; call transform_norm(self,i,z); self%norm_f(i)=rastrigin(z); end do
         do i=5,6; call transform_norm(self,i,z); self%norm_f(i)=ef8f2(z); end do
         do i=7,8; call transform_norm(self,i,z); self%norm_f(i)=weierstrass(z)-weierstrass(zero); end do
         do i=9,10; call transform_norm(self,i,z); self%norm_f(i)=griewank(z); end do
      case(24,25)
         call transform_norm(self,1,z); self%norm_f(1)=weierstrass(z)-weierstrass(zero)
         call transform_norm(self,2,z); self%norm_f(2)=escaffer(z)
         call transform_norm(self,3,z); self%norm_f(3)=ef8f2(z)
         call transform_norm(self,4,z); self%norm_f(4)=ackley(z)
         call transform_norm(self,5,z); self%norm_f(5)=rastrigin(z)
         call transform_norm(self,6,z); self%norm_f(6)=griewank(z)
         call transform_norm(self,7,z); self%norm_f(7)=nc_escaffer(z)
         call transform_norm(self,8,z); self%norm_f(8)=nc_rastrigin_fn(z)
         call transform_norm(self,9,z); self%norm_f(9)=elliptic(z)
         call transform_norm(self,10,z); self%norm_f(10)=sphere(z)
         if (self%noise_enabled) self%norm_f(10)=self%norm_f(10)*(1.0_dp+0.1_dp*abs(normal_deviate()))
      end select
   end subroutine compute_norms

   real(dp) function nc_escaffer(x) result(f)
      real(dp), intent(in) :: x(:)
      integer :: i, j, n
      n=size(x); f=0.0_dp
      do i=1,n
         j=merge(i+1,1,i<n)
         f=f+nc_schaffer_pair(x(i),x(j))
      end do
   end function nc_escaffer

   function rounded_for_f23(self, x) result(y)
      class(cec2005_context), intent(in) :: self
      real(dp), intent(in) :: x(:)
      real(dp) :: y(size(x)), q, b
      integer :: i, a
      do i=1,size(x)
         if (abs(x(i)-self%o(1,i)) >= 0.5_dp) then
            q=2.0_dp*x(i); a=int(q); b=abs(q-real(a,dp))
            if (b < 0.5_dp) then
               y(i)=real(a,dp)/2.0_dp
            else if (q <= 0.0_dp) then
               y(i)=real(a-1,dp)/2.0_dp
            else
               y(i)=real(a+1,dp)/2.0_dp
            end if
         else
            y(i)=x(i)
         end if
      end do
   end function rounded_for_f23

   real(dp) function composition_value(self, x, kind) result(f)
      class(cec2005_context), intent(in) :: self
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: kind
      real(dp) :: z(self%n), zero(self%n), b(nfunc), w(nfunc), xx(self%n)
      integer :: i
      zero=0.0_dp; b=0.0_dp; xx=x
      if (kind == 3 .and. self%function_id == 23) xx=rounded_for_f23(self,x)
      select case(kind)
      case(1)
         call transform(self,xx,1,z); b(1)=rastrigin(z)
         call transform(self,xx,2,z); b(2)=rastrigin(z)
         do i=3,4; call transform(self,xx,i,z); b(i)=weierstrass(z)-weierstrass(zero); end do
         do i=5,6; call transform(self,xx,i,z); b(i)=griewank(z); end do
         do i=7,8; call transform(self,xx,i,z); b(i)=ackley(z); end do
         do i=9,10; call transform(self,xx,i,z); b(i)=sphere(z); end do
      case(2)
         do i=1,2; call transform(self,xx,i,z); b(i)=ackley(z); end do
         do i=3,4; call transform(self,xx,i,z); b(i)=rastrigin(z); end do
         do i=5,6; call transform(self,xx,i,z); b(i)=sphere(z); end do
         do i=7,8; call transform(self,xx,i,z); b(i)=weierstrass(z)-weierstrass(zero); end do
         do i=9,10; call transform(self,xx,i,z); b(i)=griewank(z); end do
      case(3)
         do i=1,2; call transform(self,xx,i,z); b(i)=escaffer(z); end do
         do i=3,4; call transform(self,xx,i,z); b(i)=rastrigin(z); end do
         do i=5,6; call transform(self,xx,i,z); b(i)=ef8f2(z); end do
         do i=7,8; call transform(self,xx,i,z); b(i)=weierstrass(z)-weierstrass(zero); end do
         do i=9,10; call transform(self,xx,i,z); b(i)=griewank(z); end do
      case(4)
         call transform(self,xx,1,z); b(1)=weierstrass(z)-weierstrass(zero)
         call transform(self,xx,2,z); b(2)=escaffer(z)
         call transform(self,xx,3,z); b(3)=ef8f2(z)
         call transform(self,xx,4,z); b(4)=ackley(z)
         call transform(self,xx,5,z); b(5)=rastrigin(z)
         call transform(self,xx,6,z); b(6)=griewank(z)
         call transform(self,xx,7,z); b(7)=nc_escaffer(z)
         call transform(self,xx,8,z); b(8)=nc_rastrigin_fn(z)
         call transform(self,xx,9,z); b(9)=elliptic(z)
         call transform(self,xx,10,z); b(10)=sphere(z)
         if (self%noise_enabled) b(10)=b(10)*(1.0_dp+0.1_dp*abs(normal_deviate()))
      end select
      b = b*comp_scale/self%norm_f
      call calc_weights(self,xx,w)
      f = self%global_bias + sum(w*(b+self%bias))
      if (self%function_id == 17 .and. self%noise_enabled) then
         f = (f-self%global_bias)*(1.0_dp+0.2_dp*abs(normal_deviate())) + self%global_bias
      end if
   end function composition_value

   real(dp) function cec_evaluate(self, x, ierr) result(f)
      class(cec2005_context), intent(inout) :: self
      real(dp), intent(in) :: x(:)
      integer, intent(out), optional :: ierr
      real(dp) :: z(self%n), s1, s2
      integer :: i, j, stat
      stat=0; f=huge(1.0_dp)
      if (.not.self%initialized .or. size(x)/=self%n) stat=1
      if (present(ierr)) ierr=stat
      if (stat/=0) return
      select case(self%function_id)
      case(1)
         call transform(self,x,1,z); f=sphere(z)+self%bias(1)
      case(2)
         call transform(self,x,1,z); f=schwefel_102(z)+self%bias(1)
      case(3)
         call transform(self,x,1,z); f=elliptic(z)+self%bias(1)
      case(4)
         call transform(self,x,1,z); f=schwefel_102(z)
         if (self%noise_enabled) f=f*(1.0_dp+0.4_dp*abs(normal_deviate()))
         f=f+self%bias(1)
      case(5)
         f=-huge(1.0_dp)
         do i=1,self%n
            f=max(f,abs(dot_product(self%a5(i,:),x)-self%b5(i)))
         end do
         f=f+self%bias(1)
      case(6)
         call transform(self,x,1,z); f=rosenbrock(z)+self%bias(1)
      case(7)
         call transform(self,x,1,z); f=griewank(z)+self%bias(1)
      case(8)
         call transform(self,x,1,z); f=ackley(z)+self%bias(1)
      case(9,10)
         call transform(self,x,1,z); f=rastrigin(z)+self%bias(1)
      case(11)
         call transform(self,x,1,z); f=weierstrass(z)-weierstrass(0.0_dp*z)+self%bias(1)
      case(12)
         f=0.0_dp
         do i=1,self%n
            s1=0.0_dp; s2=0.0_dp
            do j=1,self%n
               s1=s1+self%a12(i,j)*sin(self%alpha12(j))+self%b12(i,j)*cos(self%alpha12(j))
               s2=s2+self%a12(i,j)*sin(x(j))+self%b12(i,j)*cos(x(j))
            end do
            f=f+(s1-s2)**2
         end do
         f=f+self%bias(1)
      case(13)
         call transform(self,x,1,z); f=ef8f2(z)+self%bias(1)
      case(14)
         call transform(self,x,1,z); f=escaffer(z)+self%bias(1)
      case(15:17)
         f=composition_value(self,x,1)
      case(18:20)
         f=composition_value(self,x,2)
      case(21:23)
         f=composition_value(self,x,3)
      case(24:25)
         f=composition_value(self,x,4)
      end select
   end function cec_evaluate

   subroutine cec_evaluate_batch(self, x, f, ierr)
      class(cec2005_context), intent(inout) :: self
      real(dp), intent(in) :: x(:,:)
      real(dp), intent(out) :: f(:)
      integer, intent(out), optional :: ierr
      integer :: i, stat
      stat=0
      if (size(x,2)/=self%n .or. size(f)/=size(x,1)) stat=1
      if (present(ierr)) ierr=stat
      if (stat/=0) return
      do i=1,size(x,1)
         f(i)=self%evaluate(x(i,:))
      end do
   end subroutine cec_evaluate_batch

   real(dp) function cec2005_eval(function_id, x, data_dir, noise_enabled, ierr) result(f)
      integer, intent(in) :: function_id
      real(dp), intent(in) :: x(:)
      character(len=*), intent(in), optional :: data_dir
      logical, intent(in), optional :: noise_enabled
      integer, intent(out), optional :: ierr
      type(cec2005_context) :: ctx
      integer :: stat
      character(len=:), allocatable :: d
      logical :: noise
      d='data'; if(present(data_dir)) d=data_dir
      noise=.true.; if(present(noise_enabled)) noise=noise_enabled
      call ctx%init(function_id,size(x),d,noise,stat)
      if(stat==0) then
         f=ctx%evaluate(x,stat)
      else
         f=huge(1.0_dp)
      end if
      if(present(ierr)) ierr=stat
   end function cec2005_eval

   subroutine cec2005_eval_batch(function_id, x, f, data_dir, noise_enabled, ierr)
      integer, intent(in) :: function_id
      real(dp), intent(in) :: x(:,:)
      real(dp), intent(out) :: f(:)
      character(len=*), intent(in), optional :: data_dir
      logical, intent(in), optional :: noise_enabled
      integer, intent(out), optional :: ierr
      type(cec2005_context) :: ctx
      integer :: stat
      character(len=:), allocatable :: d
      logical :: noise
      d='data'; if(present(data_dir)) d=data_dir
      noise=.true.; if(present(noise_enabled)) noise=noise_enabled
      call ctx%init(function_id,size(x,2),d,noise,stat)
      if(stat==0) call ctx%evaluate_batch(x,f,stat)
      if(present(ierr)) ierr=stat
   end subroutine cec2005_eval_batch

end module cec2005benchmark
