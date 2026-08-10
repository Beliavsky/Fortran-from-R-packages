! SPDX-License-Identifier: GPL-2.0-only
module nlsic_linalg
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   use nlsic_kinds, only : dp
   implicit none
   private
   public :: dense_solve, symmetric_eigen, pseudoinverse, matrix_rank
   public :: null_space, least_norm_ls, vecnorm, identity_matrix
   public :: matrix_inverse, student_t_quantile
contains

   pure real(dp) function vecnorm(x) result(v)
      real(dp), intent(in) :: x(:)
      v = sqrt(max(0.0_dp, dot_product(x, x)))
   end function vecnorm

   function identity_matrix(n) result(a)
      integer, intent(in) :: n
      real(dp), allocatable :: a(:,:)
      integer :: i
      allocate(a(n,n)); a = 0.0_dp
      do i = 1, n
         a(i,i) = 1.0_dp
      end do
   end function identity_matrix

   subroutine dense_solve(a, b, x, info)
      real(dp), intent(in) :: a(:,:), b(:)
      real(dp), intent(out) :: x(:)
      integer, intent(out) :: info
      real(dp), allocatable :: m(:,:), rhs(:), row(:)
      real(dp) :: piv, fac, tmp, scale
      integer :: n, i, j, k, p
      n = size(a,1); info = 0; x = 0.0_dp
      if (size(a,2) /= n .or. size(b) /= n .or. size(x) /= n) then
         info = -1; return
      end if
      if (n == 0) return
      allocate(m(n,n), rhs(n), row(n)); m = a; rhs = b
      scale = max(1.0_dp, maxval(abs(m)))
      do k = 1, n-1
         p = k; piv = abs(m(k,k))
         do i = k+1, n
            if (abs(m(i,k)) > piv) then
               piv = abs(m(i,k)); p = i
            end if
         end do
         if (piv <= 100.0_dp*epsilon(1.0_dp)*scale) then
            info = k; return
         end if
         if (p /= k) then
            row = m(k,:); m(k,:) = m(p,:); m(p,:) = row
            tmp = rhs(k); rhs(k) = rhs(p); rhs(p) = tmp
         end if
         do i = k+1, n
            fac = m(i,k)/m(k,k); m(i,k) = 0.0_dp
            do j = k+1, n
               m(i,j) = m(i,j) - fac*m(k,j)
            end do
            rhs(i) = rhs(i) - fac*rhs(k)
         end do
      end do
      if (abs(m(n,n)) <= 100.0_dp*epsilon(1.0_dp)*scale) then
         info = n; return
      end if
      x(n) = rhs(n)/m(n,n)
      do i = n-1, 1, -1
         x(i) = (rhs(i)-dot_product(m(i,i+1:n),x(i+1:n)))/m(i,i)
      end do
   end subroutine dense_solve

   subroutine matrix_inverse(a, ainv, info)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(out) :: ainv(:,:)
      integer, intent(out) :: info
      real(dp), allocatable :: rhs(:), col(:)
      integer :: n, j, st
      n = size(a,1); info = 0; ainv = 0.0_dp
      if (size(a,2) /= n .or. size(ainv,1) /= n .or. size(ainv,2) /= n) then
         info = -1; return
      end if
      allocate(rhs(n), col(n))
      do j = 1, n
         rhs = 0.0_dp; rhs(j) = 1.0_dp
         call dense_solve(a, rhs, col, st)
         if (st /= 0) then
            info = st; return
         end if
         ainv(:,j) = col
      end do
   end subroutine matrix_inverse

   subroutine symmetric_eigen(a, eval, evec, info)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(out) :: eval(:), evec(:,:)
      integer, intent(out) :: info
      real(dp), allocatable :: d(:,:)
      real(dp) :: app, aqq, apq, tau, t, c, s, dip, diq, vip, viq
      real(dp) :: off, eps, tmp, scale
      integer :: n, p, q, i, sweep, imax
      n = size(a,1); info = 0
      if (size(a,2) /= n .or. size(eval) /= n .or. size(evec,1) /= n .or. size(evec,2) /= n) then
         info = -1; return
      end if
      if (n == 0) return
      allocate(d(n,n)); d = 0.5_dp*(a+transpose(a)); evec = 0.0_dp
      do i = 1, n
         evec(i,i) = 1.0_dp
      end do
      eps = 50.0_dp*epsilon(1.0_dp)
      do sweep = 1, max(80, 50*n*n)
         off = 0.0_dp
         do p = 1, n-1
            do q = p+1, n
               off = max(off, abs(d(p,q)))
            end do
         end do
         scale = max(1.0_dp, maxval(abs(d)))
         if (off <= eps*scale) exit
         do p = 1, n-1
            do q = p+1, n
               apq = d(p,q)
               if (abs(apq) <= eps*max(1.0_dp,abs(d(p,p))+abs(d(q,q)))) cycle
               app = d(p,p); aqq = d(q,q); tau = (aqq-app)/(2.0_dp*apq)
               if (tau >= 0.0_dp) then
                  t = 1.0_dp/(tau+sqrt(1.0_dp+tau*tau))
               else
                  t = -1.0_dp/(-tau+sqrt(1.0_dp+tau*tau))
               end if
               c = 1.0_dp/sqrt(1.0_dp+t*t); s = t*c
               do i = 1, n
                  if (i == p .or. i == q) cycle
                  dip = d(i,p); diq = d(i,q)
                  d(i,p) = c*dip-s*diq; d(p,i) = d(i,p)
                  d(i,q) = s*dip+c*diq; d(q,i) = d(i,q)
               end do
               d(p,p) = c*c*app - 2.0_dp*s*c*apq + s*s*aqq
               d(q,q) = s*s*app + 2.0_dp*s*c*apq + c*c*aqq
               d(p,q) = 0.0_dp; d(q,p) = 0.0_dp
               do i = 1, n
                  vip = evec(i,p); viq = evec(i,q)
                  evec(i,p) = c*vip-s*viq
                  evec(i,q) = s*vip+c*viq
               end do
            end do
         end do
      end do
      do i = 1, n
         eval(i) = d(i,i)
      end do
      do i = 1, n-1
         imax = i
         do p = i+1, n
            if (eval(p) > eval(imax)) imax = p
         end do
         if (imax /= i) then
            tmp = eval(i); eval(i) = eval(imax); eval(imax) = tmp
            do q = 1, n
               tmp = evec(q,i); evec(q,i) = evec(q,imax); evec(q,imax) = tmp
            end do
         end if
      end do
   end subroutine symmetric_eigen

   subroutine pseudoinverse(a, ap, rank, rcond, info)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(out) :: ap(:,:)
      integer, intent(out) :: rank
      real(dp), intent(in), optional :: rcond
      integer, intent(out), optional :: info
      real(dp), allocatable :: ata(:,:), ev(:), v(:,:), av(:)
      real(dp) :: tol, smax, cutoff
      integer :: m, n, i, st
      m = size(a,1); n = size(a,2); ap = 0.0_dp; rank = 0; st = 0
      if (size(ap,1) /= n .or. size(ap,2) /= m) then
         st = -1; if (present(info)) info = st; return
      end if
      if (m == 0 .or. n == 0) then
         if (present(info)) info = 0; return
      end if
      allocate(ata(n,n), ev(n), v(n,n), av(m)); ata = matmul(transpose(a),a)
      call symmetric_eigen(ata, ev, v, st)
      tol = 1.0e-10_dp
      if (present(rcond)) then
         if (rcond > 1.0_dp) tol = 1.0_dp/rcond
      end if
      smax = sqrt(max(0.0_dp, ev(1)))
      cutoff=max((tol*max(1.0_dp,smax))**2, &
         100.0_dp*epsilon(1.0_dp)*max(1.0_dp,abs(ev(1))))
      do i = 1, n
         if (ev(i) > cutoff) then
            av = matmul(a,v(:,i))
            ap = ap + spread(v(:,i)/ev(i),2,m)*spread(av,1,n)
            rank = rank+1
         end if
      end do
      if (present(info)) info = st
   end subroutine pseudoinverse

   integer function matrix_rank(a, rcond) result(rank)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(in), optional :: rcond
      real(dp), allocatable :: ap(:,:)
      integer :: st
      allocate(ap(size(a,2),size(a,1)))
      if (present(rcond)) then
         call pseudoinverse(a,ap,rank,rcond,st)
      else
         call pseudoinverse(a,ap,rank,info=st)
      end if
   end function matrix_rank

   subroutine null_space(a, z, rank, rcond, info)
      real(dp), intent(in) :: a(:,:)
      real(dp), allocatable, intent(out) :: z(:,:)
      integer, intent(out) :: rank
      real(dp), intent(in), optional :: rcond
      integer, intent(out), optional :: info
      real(dp), allocatable :: ata(:,:), ev(:), v(:,:)
      real(dp) :: tol, smax, cutoff
      integer :: n, st
      n = size(a,2)
      if (n == 0) then
         rank = 0; allocate(z(0,0)); if (present(info)) info = 0; return
      end if
      allocate(ata(n,n),ev(n),v(n,n)); ata = matmul(transpose(a),a)
      call symmetric_eigen(ata,ev,v,st)
      tol = 1.0e-10_dp
      if (present(rcond)) then
         if (rcond > 1.0_dp) tol = 1.0_dp/rcond
      end if
      smax=sqrt(max(0.0_dp,ev(1)))
      cutoff=max((tol*max(1.0_dp,smax))**2, &
         100.0_dp*epsilon(1.0_dp)*max(1.0_dp,abs(ev(1))))
      rank=count(ev>cutoff)
      allocate(z(n,n-rank))
      if (n-rank > 0) z = v(:,rank+1:n)
      if (present(info)) info = st
   end subroutine null_space

   subroutine least_norm_ls(a, b, x, rank, rcond, mnorm, x0, info)
      real(dp), intent(in) :: a(:,:), b(:)
      real(dp), intent(out) :: x(:)
      integer, intent(out) :: rank
      real(dp), intent(in), optional :: rcond
      real(dp), intent(in), optional :: mnorm(:,:), x0(:)
      integer, intent(out), optional :: info
      real(dp), allocatable :: ap(:,:), xp(:), z(:,:), mb(:,:), rhs(:), pinv_mb(:,:), zcoef(:), base(:)
      integer :: m, n, st, rankz, rankm
      m = size(a,1); n = size(a,2); x = 0.0_dp; st = 0
      if (size(b) /= m .or. size(x) /= n) then
         st = -1; rank = 0; if (present(info)) info = st; return
      end if
      allocate(ap(n,m),xp(n))
      if (present(rcond)) then
         call pseudoinverse(a,ap,rank,rcond,st)
      else
         call pseudoinverse(a,ap,rank,info=st)
      end if
      xp = matmul(ap,b)
      if (rank == n) then
         x = xp; if (present(info)) info = st; return
      end if
      if (present(rcond)) then
         call null_space(a,z,rankz,rcond,st)
      else
         call null_space(a,z,rankz,info=st)
      end if
      if (size(z,2) == 0) then
         x = xp; if (present(info)) info = st; return
      end if
      allocate(base(n)); base = xp
      if (present(x0)) base = xp-x0
      if (present(mnorm)) then
         allocate(mb(size(mnorm,1),size(z,2)),rhs(size(mnorm,1)))
         mb = matmul(mnorm,z); rhs = -matmul(mnorm,base)
      else
         allocate(mb(n,size(z,2)),rhs(n)); mb = z; rhs = -base
      end if
      allocate(pinv_mb(size(z,2),size(mb,1)),zcoef(size(z,2)))
      if (present(rcond)) then
         call pseudoinverse(mb,pinv_mb,rankm,rcond,st)
      else
         call pseudoinverse(mb,pinv_mb,rankm,info=st)
      end if
      zcoef=matmul(pinv_mb,rhs)
      x=xp+matmul(z,zcoef)
      if (present(info)) info=st
   end subroutine least_norm_ls

   real(dp) function beta_cont_frac(a,b,x) result(cf)
      real(dp), intent(in) :: a,b,x
      integer, parameter :: maxit=300
      real(dp), parameter :: fpmin=1.0e-300_dp, eps=3.0e-14_dp
      real(dp) :: qab,qap,qam,c,d,h,aa,del
      integer :: m,m2
      qab=a+b; qap=a+1.0_dp; qam=a-1.0_dp
      c=1.0_dp; d=1.0_dp-qab*x/qap
      if(abs(d)<fpmin) d=fpmin
      d=1.0_dp/d; h=d
      do m=1,maxit
         m2=2*m
         aa=m*(b-m)*x/((qam+m2)*(a+m2))
         d=1.0_dp+aa*d; if(abs(d)<fpmin) d=fpmin
         c=1.0_dp+aa/c; if(abs(c)<fpmin) c=fpmin
         d=1.0_dp/d; h=h*d*c
         aa=-(a+m)*(qab+m)*x/((a+m2)*(qap+m2))
         d=1.0_dp+aa*d; if(abs(d)<fpmin) d=fpmin
         c=1.0_dp+aa/c; if(abs(c)<fpmin) c=fpmin
         d=1.0_dp/d; del=d*c; h=h*del
         if(abs(del-1.0_dp)<eps) exit
      end do
      cf=h
   end function beta_cont_frac

   real(dp) function regularized_beta(x,a,b) result(v)
      real(dp), intent(in) :: x,a,b
      real(dp) :: bt
      if(x<=0.0_dp) then; v=0.0_dp; return; end if
      if(x>=1.0_dp) then; v=1.0_dp; return; end if
      bt=exp(log_gamma(a+b)-log_gamma(a)-log_gamma(b)+a*log(x)+b*log(1.0_dp-x))
      if(x < (a+1.0_dp)/(a+b+2.0_dp)) then
         v=bt*beta_cont_frac(a,b,x)/a
      else
         v=1.0_dp-bt*beta_cont_frac(b,a,1.0_dp-x)/b
      end if
      v=max(0.0_dp,min(1.0_dp,v))
   end function regularized_beta

   real(dp) function student_t_cdf(t,nu) result(p)
      real(dp), intent(in) :: t,nu
      real(dp) :: x,ib
      if(nu<=0.0_dp) then; p=ieee_value(0.0_dp,ieee_quiet_nan); return; end if
      if(abs(t)<=tiny(1.0_dp)) then; p=0.5_dp; return; end if
      x=nu/(nu+t*t); ib=regularized_beta(x,0.5_dp*nu,0.5_dp)
      if(t>0.0_dp) then; p=1.0_dp-0.5_dp*ib; else; p=0.5_dp*ib; end if
   end function student_t_cdf

   real(dp) function student_t_quantile(p,nu) result(q)
      real(dp), intent(in) :: p,nu
      real(dp) :: lo,hi,mid,pm
      integer :: i
      if(p<=0.0_dp) then; q=-huge(1.0_dp); return; end if
      if(p>=1.0_dp) then; q=huge(1.0_dp); return; end if
      if(abs(p-0.5_dp)<=tiny(1.0_dp)) then; q=0.0_dp; return; end if
      lo=-1.0_dp; hi=1.0_dp
      do while(student_t_cdf(lo,nu)>p); lo=2.0_dp*lo; end do
      do while(student_t_cdf(hi,nu)<p); hi=2.0_dp*hi; end do
      do i=1,160
         mid=0.5_dp*(lo+hi); pm=student_t_cdf(mid,nu)
         if(pm<p) then; lo=mid; else; hi=mid; end if
      end do
      q=0.5_dp*(lo+hi)
   end function student_t_quantile
end module nlsic_linalg
