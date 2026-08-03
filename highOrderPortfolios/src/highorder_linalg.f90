! SPDX-License-Identifier: GPL-3.0-only
module highorder_linalg
   use fitheavytail_kinds, only: dp
   use fitheavytail_linalg, only: symmetric_eigen
   use fitheavytail_status, only: ht_success
   implicit none
   private
   public :: project_simplex, solve_simplex_qp, psd_projection
   public :: spectral_radius, vector_norm2, cholesky_upper, identity

contains

   pure function vector_norm2(x) result(v)
      real(dp), intent(in) :: x(:)
      real(dp) :: v
      v = sqrt(max(0.0_dp, dot_product(x,x)))
   end function vector_norm2

   subroutine identity(a)
      real(dp), intent(out) :: a(:,:)
      integer :: i
      a = 0.0_dp
      do i = 1, min(size(a,1),size(a,2))
         a(i,i) = 1.0_dp
      end do
   end subroutine identity

   subroutine project_simplex(v,w)
      real(dp), intent(in) :: v(:)
      real(dp), intent(out) :: w(:)
      real(dp), allocatable :: u(:)
      real(dp) :: cssv, theta
      integer :: i, j, n, rho
      n = size(v)
      if (size(w) /= n) error stop 'project_simplex: size mismatch'
      allocate(u(n))
      u = v
      do i = 2, n
         theta = u(i)
         j = i - 1
         do while (j >= 1)
            if (u(j) >= theta) exit
            u(j+1) = u(j)
            j = j - 1
         end do
         u(j+1) = theta
      end do
      cssv = 0.0_dp
      rho = 1
      do i = 1, n
         cssv = cssv + u(i)
         theta = (cssv - 1.0_dp)/real(i,dp)
         if (u(i) > theta) rho = i
      end do
      theta = (sum(u(1:rho))-1.0_dp)/real(rho,dp)
      w = max(v-theta,0.0_dp)
      if (sum(w) > 0.0_dp) then
         w = w/sum(w)
      else
         w = 1.0_dp/real(n,dp)
      end if
   end subroutine project_simplex

   function spectral_radius(a) result(rho)
      real(dp), intent(in) :: a(:,:)
      real(dp) :: rho
      real(dp), allocatable :: eval(:), evec(:,:)
      integer :: status, n
      n = size(a,1)
      if (size(a,2) /= n) then
         rho = huge(1.0_dp)
         return
      end if
      allocate(eval(n),evec(n,n))
      call symmetric_eigen(0.5_dp*(a+transpose(a)),eval,evec,status)
      if (status /= ht_success) then
         rho = maxval(sum(abs(a),dim=2))
      else
         rho = maxval(abs(eval))
      end if
   end function spectral_radius

   subroutine psd_projection(a,ap)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(out) :: ap(:,:)
      real(dp), allocatable :: eval(:), evec(:,:)
      real(dp) :: shift
      integer :: status, i, n
      n = size(a,1)
      if (size(a,2) /= n .or. any(shape(ap) /= shape(a))) then
         ap = 0.0_dp
         return
      end if
      allocate(eval(n),evec(n,n))
      call symmetric_eigen(0.5_dp*(a+transpose(a)),eval,evec,status)
      if (status /= ht_success) then
         ap = 0.5_dp*(a+transpose(a))
         shift=0.0_dp
         do i=1,n
            shift=min(shift,ap(i,i))
         end do
         shift=max(0.0_dp,-shift)+1.0e-12_dp
         do i=1,n
            ap(i,i)=ap(i,i)+shift
         end do
         return
      end if
      eval = max(eval,0.0_dp)
      ap = matmul(evec,matmul(diag_matrix(eval),transpose(evec)))
      ap = 0.5_dp*(ap+transpose(ap))
   end subroutine psd_projection

   function diag_matrix(d) result(a)
      real(dp), intent(in) :: d(:)
      real(dp) :: a(size(d),size(d))
      integer :: i
      a = 0.0_dp
      do i=1,size(d)
         a(i,i)=d(i)
      end do
   end function diag_matrix

   subroutine solve_simplex_qp(qmat,qvec,x0,x,max_iter,tol)
      real(dp), intent(in) :: qmat(:,:), qvec(:), x0(:)
      real(dp), intent(out) :: x(:)
      integer, intent(in), optional :: max_iter
      real(dp), intent(in), optional :: tol
      real(dp), allocatable :: y(:), xnew(:), grad(:)
      real(dp) :: lips, step, t, tnew, eps
      integer :: k, n, niter
      n = size(x0)
      niter = 1000
      if (present(max_iter)) niter = max_iter
      eps = 1.0e-10_dp
      if (present(tol)) eps = tol
      allocate(y(n),xnew(n),grad(n))
      call project_simplex(x0,x)
      y = x
      t = 1.0_dp
      lips = max(spectral_radius(qmat),1.0e-10_dp)
      step = 1.0_dp/lips
      do k=1,niter
         grad = matmul(qmat,y)-qvec
         call project_simplex(y-step*grad,xnew)
         if (vector_norm2(xnew-x) <= eps*max(1.0_dp,vector_norm2(x))) then
            x = xnew
            exit
         end if
         tnew = 0.5_dp*(1.0_dp+sqrt(1.0_dp+4.0_dp*t*t))
         y = xnew + ((t-1.0_dp)/tnew)*(xnew-x)
         x = xnew
         t = tnew
      end do
   end subroutine solve_simplex_qp

   subroutine cholesky_upper(a,r,status)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(out) :: r(:,:)
      integer, intent(out) :: status
      real(dp) :: s
      integer :: i,j,k,n
      n=size(a,1)
      r=0.0_dp
      status=0
      if(size(a,2)/=n .or. any(shape(r)/=shape(a))) then
         status=1
         return
      end if
      do j=1,n
         s=a(j,j)
         do k=1,j-1
            s=s-r(k,j)*r(k,j)
         end do
         if(s<=tiny(1.0_dp)) then
            status=2
            return
         end if
         r(j,j)=sqrt(s)
         do i=j+1,n
            s=a(j,i)
            do k=1,j-1
               s=s-r(k,j)*r(k,i)
            end do
            r(j,i)=s/r(j,j)
         end do
      end do
   end subroutine cholesky_upper

end module highorder_linalg
