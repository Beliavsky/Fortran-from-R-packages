module rvmf
   use, intrinsic :: iso_fortran_env, only : int64
   implicit none
   private

   integer, parameter, public :: dp = kind(1.0d0)
   real(dp), parameter :: pi = acos(-1.0_dp)

   public :: rvmf_sample, rvmf_angle_sample, dvmf_angle, dvmf_angle_vec, log_chf

contains

   pure real(dp) function logaddexp(a, b) result(c)
      real(dp), intent(in) :: a, b
      real(dp) :: m
      m = max(a, b)
      c = m + log(exp(a-m) + exp(b-m))
   end function logaddexp

   real(dp) function log_chf(kappa, d1) result(logm)
      ! log M(d1/2, d1, 2*kappa), evaluated by a positive-term
      ! log-sum-exp Taylor series. This is the same confluent
      ! hypergeometric normalization used by upstream rvMF.
      real(dp), intent(in) :: kappa, d1
      real(dp) :: a, b, z, logterm, newsum
      integer :: j
      real(dp), parameter :: reltol = 1.0e-14_dp
      integer, parameter :: maxit = 1000000

      if (kappa < 0.0_dp .or. d1 <= 0.0_dp) error stop "log_chf: invalid parameters"
      if (kappa <= 0.0_dp) then
         logm = 0.0_dp
         return
      end if

      a = 0.5_dp*d1
      b = d1
      z = 2.0_dp*kappa
      logterm = 0.0_dp
      logm = 0.0_dp
      do j = 0, maxit-1
         logterm = logterm + log(a + real(j,dp)) + log(z) - &
              log(b + real(j,dp)) - log(real(j+1,dp))
         newsum = logaddexp(logm, logterm)
         if (logterm - newsum < log(reltol) .and. j > int(z)) then
            logm = newsum
            return
         end if
         logm = newsum
      end do
      error stop "log_chf: series did not converge"
   end function log_chf

   subroutine random_normal(z)
      real(dp), intent(out) :: z
      real(dp) :: u1, u2
      call random_number(u1)
      call random_number(u2)
      u1 = max(u1, tiny(1.0_dp))
      z = sqrt(-2.0_dp*log(u1))*cos(2.0_dp*pi*u2)
   end subroutine random_normal

   recursive real(dp) function gamma_rng(shape) result(x)
      real(dp), intent(in) :: shape
      real(dp) :: d, c, z, u, v, g
      if (shape <= 0.0_dp) error stop "gamma_rng: shape must be positive"
      if (shape < 1.0_dp) then
         call random_number(u)
         x = gamma_rng(shape + 1.0_dp)*max(u,tiny(1.0_dp))**(1.0_dp/shape)
         return
      end if
      d = shape - 1.0_dp/3.0_dp
      c = 1.0_dp/sqrt(9.0_dp*d)
      do
         call random_normal(z)
         v = 1.0_dp + c*z
         if (v <= 0.0_dp) cycle
         v = v*v*v
         call random_number(u)
         if (u < 1.0_dp - 0.0331_dp*z**4) exit
         if (log(max(u,tiny(1.0_dp))) < 0.5_dp*z*z + d*(1.0_dp-v+log(v))) exit
      end do
      g = d*v
      x = g
   end function gamma_rng

   real(dp) function beta_rng(a, b) result(x)
      real(dp), intent(in) :: a, b
      real(dp) :: ga, gb
      if (a <= 0.0_dp .or. b <= 0.0_dp) error stop "beta_rng: shapes must be positive"
      ga = gamma_rng(a)
      gb = gamma_rng(b)
      x = ga/(ga+gb)
   end function beta_rng

   integer function digit64(m, k) result(d)
      integer(int64), intent(in) :: m
      integer, intent(in) :: k
      d = int(ibits(m, 30-6*k, 6))
   end function digit64

   subroutine build_probabilities(kappa, m, log_conf, ptab, offset)
      real(dp), intent(in) :: kappa, log_conf
      integer, intent(in) :: m
      integer(int64), allocatable, intent(out) :: ptab(:)
      integer, intent(out) :: offset
      integer :: i, imax, last, first
      real(dp) :: p, t, logp
      real(dp), parameter :: two30 = 1073741824.0_dp
      real(dp), parameter :: two31 = 2147483648.0_dp

      if (kappa <= 0.0_dp) error stop "build_probabilities: kappa must be positive"
      if (m < 1) error stop "build_probabilities: invalid dimension"

      if (log_conf <= log(two30)) then
         p = exp(-log_conf)
         t = p
         i = 1
         do while (t*two31 > 1.0_dp)
            t = t*(2.0_dp*kappa*(0.5_dp*real(m,dp)+real(i-1,dp)) / &
                 ((real(m+i-1,dp))*real(i,dp)))
            i = i + 1
         end do
         last = i-2
         offset = 0
         allocate(ptab(last+1))
         ptab(1) = int(exp(-log_conf)*two30 + 0.5_dp, int64)
         p = exp(-log_conf)
         do i = 1, last
            p = p*(2.0_dp*kappa*(0.5_dp*real(m,dp)+real(i-1,dp)) / &
                 (real(m+i-1,dp)*real(i,dp)))
            ptab(i+1) = int(p*two30 + 0.5_dp, int64)
         end do
      else
         imax = ceiling((-real(m,dp)-1.0_dp+2.0_dp*kappa + &
              sqrt((real(m,dp)+1.0_dp-2.0_dp*kappa)**2 - &
              4.0_dp*real(m,dp)*(1.0_dp-kappa)))/2.0_dp)
         imax = max(0, imax)
         logp = -log_conf + log_gamma(real(m,dp)) - log_gamma(0.5_dp*real(m,dp)) + &
              real(imax,dp)*log(2.0_dp*kappa) + log_gamma(0.5_dp*real(m,dp)+real(imax,dp)) - &
              log_gamma(real(imax+1,dp)) - log_gamma(real(m+imax,dp))
         p = exp(logp)
         t = p
         i = imax+1
         do while (t*two31 > 1.0_dp)
            t = t*(2.0_dp*kappa*(0.5_dp*real(m,dp)+real(i-1,dp)) / &
                 (real(m+i-1,dp)*real(i,dp)))
            i = i + 1
         end do
         last = i-2
         t = p
         first = 0
         do i = imax-1, 0, -1
            t = t*(real(i+1,dp)*real(m+i,dp) / &
                 (2.0_dp*kappa*(0.5_dp*real(m,dp)+real(i,dp))))
            if (t*two31 < 1.0_dp) then
               first = i+1
               exit
            end if
         end do
         offset = first
         allocate(ptab(last-offset+1))
         ptab = 0_int64
         ptab(imax-offset+1) = int(p*two30 + 0.5_dp, int64)
         t = p
         do i = imax+1, last
            t = t*(2.0_dp*kappa*(0.5_dp*real(m,dp)+real(i-1,dp)) / &
                 (real(m+i-1,dp)*real(i,dp)))
            ptab(i-offset+1) = int(t*two30 + 0.5_dp, int64)
         end do
         t = p
         do i = imax-1, offset, -1
            t = t*(real(i+1,dp)*real(m+i,dp) / &
                 (2.0_dp*kappa*(0.5_dp*real(m,dp)+real(i,dp))))
            ptab(i-offset+1) = int(t*two30 + 0.5_dp, int64)
         end do
      end if
   end subroutine build_probabilities

   subroutine build_tables(ptab, aa, bb, cc, dd, ee, t1, t2, t3, t4, t5)
      integer(int64), intent(in) :: ptab(:)
      integer, allocatable, intent(out) :: aa(:), bb(:), cc(:), dd(:), ee(:)
      integer(int64), intent(out) :: t1, t2, t3, t4, t5
      integer :: i, j, d1, d2, d3, d4, d5
      integer :: na, nb, nc, nd, ne, ia, ib, ic, id, ie

      na=0; nb=0; nc=0; nd=0; ne=0
      do i=1,size(ptab)
         na=na+digit64(ptab(i),1); nb=nb+digit64(ptab(i),2)
         nc=nc+digit64(ptab(i),3); nd=nd+digit64(ptab(i),4)
         ne=ne+digit64(ptab(i),5)
      end do
      allocate(aa(na),bb(nb),cc(nc),dd(nd),ee(ne))
      t1 = shiftl(int(na,int64),24)
      t2 = t1 + shiftl(int(nb,int64),18)
      t3 = t2 + shiftl(int(nc,int64),12)
      t4 = t3 + shiftl(int(nd,int64),6)
      t5 = t4 + int(ne,int64)
      ia=1; ib=1; ic=1; id=1; ie=1
      do i=1,size(ptab)
         d1=digit64(ptab(i),1); d2=digit64(ptab(i),2); d3=digit64(ptab(i),3)
         d4=digit64(ptab(i),4); d5=digit64(ptab(i),5)
         do j=1,d1; aa(ia)=i-1; ia=ia+1; end do
         do j=1,d2; bb(ib)=i-1; ib=ib+1; end do
         do j=1,d3; cc(ic)=i-1; ic=ic+1; end do
         do j=1,d4; dd(id)=i-1; id=id+1; end do
         do j=1,d5; ee(ie)=i-1; ie=ie+1; end do
      end do
   end subroutine build_tables

   subroutine rvmf_angle_sample(x, p, kappa)
      real(dp), intent(out) :: x(:)
      integer, intent(in) :: p
      real(dp), intent(in) :: kappa
      integer(int64), allocatable :: ptab(:)
      integer, allocatable :: aa(:),bb(:),cc(:),dd(:),ee(:)
      integer(int64) :: t1,t2,t3,t4,t5,j64
      integer :: i, u, offset
      real(dp) :: log_conf, ur, a

      if (p < 2) error stop "rvmf_angle_sample: p must be at least 2"
      if (kappa <= 0.0_dp) error stop "rvmf_angle_sample: kappa must be positive"
      log_conf = log_chf(kappa, real(p-1,dp))
      call build_probabilities(kappa,p-1,log_conf,ptab,offset)
      call build_tables(ptab,aa,bb,cc,dd,ee,t1,t2,t3,t4,t5)
      if (t5 <= 0_int64) error stop "rvmf_angle_sample: empty lookup table"
      do i=1,size(x)
         call random_number(ur)
         j64 = min(int(ur*real(t5,dp),int64),t5-1_int64)
         if (j64 < t1) then
            u = aa(int(shiftr(j64,24))+1)
         else if (j64 < t2) then
            u = bb(int(shiftr(j64-t1,18))+1)
         else if (j64 < t3) then
            u = cc(int(shiftr(j64-t2,12))+1)
         else if (j64 < t4) then
            u = dd(int(shiftr(j64-t3,6))+1)
         else
            u = ee(int(j64-t4)+1)
         end if
         a = 0.5_dp*real(p-1,dp) + real(u+offset,dp)
         x(i) = 2.0_dp*beta_rng(a,0.5_dp*real(p-1,dp))-1.0_dp
      end do
   end subroutine rvmf_angle_sample

   real(dp) function dvmf_angle(r, p, kappa) result(f)
      real(dp), intent(in) :: r, kappa
      integer, intent(in) :: p
      real(dp) :: a, logint, lr
      if (p < 2 .or. kappa < 0.0_dp .or. abs(r) > 1.0_dp) then
         f = 0.0_dp
         return
      end if
      a = 0.5_dp*real(p-3,dp)
      if (abs(r) >= 1.0_dp) then
         if (a > 0.0_dp) then
            f = 0.0_dp
            return
         else if (a < 0.0_dp) then
            f = huge(1.0_dp)
            return
         end if
      end if
      logint = -kappa + real(p-2,dp)*log(2.0_dp) + &
           2.0_dp*log_gamma(0.5_dp*real(p-1,dp)) - log_gamma(real(p-1,dp)) + &
           log_chf(kappa,real(p-1,dp))
      if (abs(r) < 1.0_dp) then
         lr = kappa*r + a*log(1.0_dp-r*r) - logint
      else
         lr = kappa*r - logint
      end if
      f = exp(lr)
   end function dvmf_angle

   subroutine dvmf_angle_vec(f, r, p, kappa)
      real(dp), intent(out) :: f(:)
      real(dp), intent(in) :: r(:), kappa
      integer, intent(in) :: p
      integer :: i
      if (size(f) /= size(r)) error stop "dvmf_angle_vec: shape mismatch"
      do i=1,size(r)
         f(i)=dvmf_angle(r(i),p,kappa)
      end do
   end subroutine dvmf_angle_vec

   subroutine rotation_from_pole(mu, rot)
      real(dp), intent(in) :: mu(:)
      real(dp), intent(out) :: rot(:,:)
      integer :: p, i, j
      real(dp), allocatable :: b(:), ca(:)
      real(dp) :: ab, nca, theta
      p=size(mu)
      if (size(rot,1)/=p .or. size(rot,2)/=p) error stop "rotation_from_pole: shape mismatch"
      allocate(b(p),ca(p)); b=0.0_dp; b(p)=1.0_dp
      ab = mu(p)
      ca = b - mu*ab
      nca = sqrt(sum(ca*ca))
      rot=0.0_dp
      do i=1,p; rot(i,i)=1.0_dp; end do
      if (nca < 32.0_dp*epsilon(1.0_dp)) then
         if (ab < 0.0_dp) rot = -rot
         return
      end if
      ca=ca/nca
      theta=acos(max(-1.0_dp,min(1.0_dp,ab)))
      do i=1,p
         do j=1,p
            rot(i,j)=rot(i,j)+sin(theta)*(mu(i)*ca(j)-ca(i)*mu(j)) + &
                 (cos(theta)-1.0_dp)*(mu(i)*mu(j)+ca(i)*ca(j))
         end do
      end do
   end subroutine rotation_from_pole

   subroutine rvmf_sample(x, mu, kappa)
      ! Each row x(i,:) is vMF(mu,kappa). mu need not be normalized.
      real(dp), intent(out) :: x(:,:)
      real(dp), intent(in) :: mu(:), kappa
      integer :: n,p,i,j
      real(dp), allocatable :: mun(:), w(:), v(:), s(:), rot(:,:)
      real(dp) :: normmu, nv, z

      n=size(x,1); p=size(x,2)
      if (p /= size(mu) .or. p < 2) error stop "rvmf_sample: shape/dimension mismatch"
      if (kappa < 0.0_dp) error stop "rvmf_sample: kappa must be nonnegative"
      normmu=sqrt(sum(mu*mu))
      if (normmu <= 0.0_dp) error stop "rvmf_sample: mu must be nonzero"
      allocate(mun(p),v(p-1),s(p),rot(p,p))
      mun=mu/normmu

      if (kappa <= 0.0_dp) then
         do i=1,n
            do j=1,p
               call random_normal(z); x(i,j)=z
            end do
            x(i,:)=x(i,:)/sqrt(sum(x(i,:)*x(i,:)))
         end do
         return
      end if

      allocate(w(n)); call rvmf_angle_sample(w,p,kappa)
      call rotation_from_pole(mun,rot)
      do i=1,n
         do j=1,p-1
            call random_normal(v(j))
         end do
         nv=sqrt(sum(v*v)); v=v/nv
         s(1:p-1)=sqrt(max(0.0_dp,1.0_dp-w(i)*w(i)))*v
         s(p)=w(i)
         x(i,:)=matmul(s,transpose(rot))
         x(i,:)=x(i,:)/sqrt(sum(x(i,:)*x(i,:)))
      end do
   end subroutine rvmf_sample

end module rvmf
