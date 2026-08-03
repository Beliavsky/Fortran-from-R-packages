! SPDX-License-Identifier: GPL-3.0-only
module sgt_updates
   use sgt_kinds, only : dp
   use sgt_status, only : sgt_ok, sgt_invalid_input, sgt_size_mismatch
   use sgt_linalg, only : symmetric_eigen_jacobi, inverse_symmetric_positive_definite, &
      diagonal_of_triple, isotonic_increasing, outer_product
   use sgt_operators, only : L, A, Lstar, Astar, Linv, vecLmat, vec
   implicit none
   private
   public :: initialize_weights_naive, initialize_weights_qp
   public :: laplacian_w_update, joint_w_update, bipartite_w_update
   public :: laplacian_u_update, bipartite_v_update
   public :: laplacian_lambda_update, bipartite_psi_update
   public :: spectral_reconstruction
contains
   pure function spectral_reconstruction(u,values) result(m)
      real(dp), intent(in) :: u(:,:), values(:)
      real(dp), allocatable :: m(:,:)
      real(dp), allocatable :: scaled(:,:)
      integer :: j
      if (size(u,2)/=size(values)) then
         allocate(m(0,0)); return
      end if
      allocate(scaled(size(u,1),size(u,2)))
      scaled=u
      do j=1,size(values)
         scaled(:,j)=scaled(:,j)*values(j)
      end do
      allocate(m(size(u,1),size(u,1)))
      m=matmul(scaled,transpose(u))
   end function spectral_reconstruction

   subroutine initialize_weights_naive(sinv,w,status,add_floor)
      real(dp), intent(in) :: sinv(:,:)
      real(dp), allocatable, intent(out) :: w(:)
      integer, intent(out), optional :: status
      real(dp), intent(in), optional :: add_floor
      real(dp) :: floor_value
      if (size(sinv,1)<2 .or. size(sinv,2)/=size(sinv,1)) then
         allocate(w(0)); if (present(status)) status=sgt_size_mismatch; return
      end if
      w=max(Linv(sinv),0.0_dp)
      floor_value=0.0_dp; if (present(add_floor)) floor_value=max(0.0_dp,add_floor)
      if (floor_value>0.0_dp) w=w+floor_value
      if (maxval(w)<=tiny(1.0_dp)) w=1.0_dp/real(size(w),dp)
      if (present(status)) status=sgt_ok
   end subroutine initialize_weights_naive

   subroutine initialize_weights_qp(sinv,w,status,max_iterations,tolerance)
      real(dp), intent(in) :: sinv(:,:)
      real(dp), allocatable, intent(out) :: w(:)
      integer, intent(out), optional :: status
      integer, intent(in), optional :: max_iterations
      real(dp), intent(in), optional :: tolerance
      real(dp), allocatable :: r(:,:),h(:,:),b(:)
      real(dp) :: lips,step,tol,den
      integer :: nedge,k,maxiter,local_status
      if (size(sinv,1)<2 .or. size(sinv,2)/=size(sinv,1)) then
         allocate(w(0)); if (present(status)) status=sgt_size_mismatch; return
      end if
      r=vecLmat(size(sinv,1))
      h=matmul(transpose(r),r)
      b=matmul(transpose(r),vec(sinv))
      nedge=size(b); allocate(w(nedge)); w=max(Linv(sinv),0.0_dp)
      if (maxval(w)<=tiny(1.0_dp)) w=1.0_dp/real(nedge,dp)
      lips=maxval(sum(abs(h),dim=2)); step=1.0_dp/max(lips,tiny(1.0_dp))
      maxiter=10000; if (present(max_iterations)) maxiter=max(1,max_iterations)
      tol=1e-10_dp; if (present(tolerance)) tol=tolerance
      local_status=4
      block
         real(dp) :: wnew(nedge),grad(nedge)
         do k=1,maxiter
            grad=matmul(h,w)-b
            wnew=max(w-step*grad,0.0_dp)
            den=max(1.0_dp,sqrt(dot_product(w,w)))
            if (sqrt(dot_product(wnew-w,wnew-w))/den<tol) then
               w=wnew; local_status=sgt_ok; exit
            end if
            w=wnew
         end do
      end block
      if (present(status)) status=local_status
   end subroutine initialize_weights_qp

   subroutine laplacian_w_update(w,lw,u,beta,lambda,kmat,wnew,status)
      real(dp), intent(in) :: w(:),lw(:,:),u(:,:),beta,lambda(:),kmat(:,:)
      real(dp), allocatable, intent(out) :: wnew(:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: c(:),grad(:),mgrad(:),target(:,:)
      real(dp) :: numerator,denominator,step
      if (beta<=0.0_dp .or. size(u,2)/=size(lambda) .or. size(w)/=size(Lstar(lw))) then
         allocate(wnew(0)); if (present(status)) status=sgt_invalid_input; return
      end if
      target=spectral_reconstruction(u,lambda)-kmat/beta
      c=Lstar(target)
      grad=Lstar(lw)-c
      mgrad=Lstar(L(grad))
      denominator=dot_product(grad,mgrad)
      allocate(wnew(size(w)))
      if (denominator<=100.0_dp*epsilon(1.0_dp)*max(1.0_dp,dot_product(grad,grad))) then
         wnew=w
      else
         numerator=dot_product(w,mgrad)-dot_product(c,grad)
         step=numerator/denominator
         if (.not.(step>0.0_dp)) step=1.0_dp/max(1.0_dp,maxval(abs(mgrad)))
         wnew=max(w-step*grad,0.0_dp)
      end if
      if (present(status)) status=sgt_ok
   end subroutine laplacian_w_update

   subroutine joint_w_update(w,lw,aw,u,v,lambda,psi,beta,nu,kmat,wnew,status)
      real(dp), intent(in) :: w(:),lw(:,:),aw(:,:),u(:,:),v(:,:),lambda(:),psi(:)
      real(dp), intent(in) :: beta,nu,kmat(:,:)
      real(dp), allocatable, intent(out) :: wnew(:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: c1(:),c2(:),mw(:),pw(:),g1(:),g2(:),grad(:),mg1(:),pg2(:)
      real(dp), allocatable :: ulut(:,:),vpvt(:,:)
      real(dp) :: numerator,denominator,step
      if (beta<0.0_dp .or. nu<0.0_dp .or. beta+nu<=0.0_dp .or. &
         size(aw,1)/=size(lw,1) .or. size(aw,2)/=size(lw,2)) then
         allocate(wnew(0)); if (present(status)) status=sgt_invalid_input; return
      end if
      ulut=spectral_reconstruction(u,lambda)
      vpvt=spectral_reconstruction(v,psi)
      c1=Lstar(beta*ulut-kmat)
      c2=nu*Astar(vpvt)
      mw=Lstar(lw); pw=2.0_dp*w
      g1=beta*mw-c1; mg1=Lstar(L(g1))
      g2=nu*pw-c2; pg2=2.0_dp*g2
      grad=g1+g2
      denominator=dot_product(grad,beta*mg1+nu*pg2)
      allocate(wnew(size(w)))
      if (denominator<=100.0_dp*epsilon(1.0_dp)*max(1.0_dp,dot_product(grad,grad))) then
         wnew=w
      else
         numerator=dot_product(beta*mw+nu*pw-(c1+c2),grad)
         step=numerator/denominator
         if (.not.(step>0.0_dp)) step=1.0_dp/max(1.0_dp,maxval(abs(beta*mg1+nu*pg2)))
         wnew=max(w-step*grad,0.0_dp)
      end if
      if (present(status)) status=sgt_ok
   end subroutine joint_w_update

   subroutine bipartite_w_update(w,aw,v,nu,psi,kmat,jmat,lips,wnew,status)
      real(dp), intent(in) :: w(:),aw(:,:),v(:,:),nu,psi(:),kmat(:,:),jmat(:,:),lips
      real(dp), allocatable, intent(out) :: wnew(:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: grad_h(:),inv_lj(:,:),grad(:),target(:,:)
      integer :: local_status
      if (nu<=0.0_dp .or. lips<0.0_dp .or. size(aw,1)/=size(v,1) .or. &
         size(aw,2)/=size(v,1)) then
         allocate(wnew(0)); if (present(status)) status=sgt_invalid_input; return
      end if
      target=spectral_reconstruction(v,psi)
      grad_h=2.0_dp*w-Astar(target)
      call inverse_symmetric_positive_definite(L(w)+jmat,inv_lj,local_status)
      allocate(wnew(size(w)))
      if (local_status/=sgt_ok) then
         wnew=w; if (present(status)) status=local_status; return
      end if
      grad=Lstar(inv_lj+kmat)+nu*grad_h
      wnew=max(w-grad/(2.0_dp*nu+lips),0.0_dp)
      if (present(status)) status=sgt_ok
   end subroutine bipartite_w_update

   subroutine laplacian_u_update(lw,k,u,status)
      real(dp), intent(in) :: lw(:,:)
      integer, intent(in) :: k
      real(dp), allocatable, intent(out) :: u(:,:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: values(:),vectors(:,:)
      integer :: n,local_status
      n=size(lw,1)
      if (k<0 .or. k>=n) then
         allocate(u(0,0)); if (present(status)) status=sgt_invalid_input; return
      end if
      call symmetric_eigen_jacobi(lw,values,vectors,local_status)
      allocate(u(n,n-k)); u=vectors(:,k+1:n)
      if (present(status)) status=local_status
   end subroutine laplacian_u_update

   subroutine bipartite_v_update(aw,z,v,status)
      real(dp), intent(in) :: aw(:,:)
      integer, intent(in) :: z
      real(dp), allocatable, intent(out) :: v(:,:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: values(:),vectors(:,:)
      integer :: n,q,h,local_status
      n=size(aw,1); q=n-z
      if (z<0 .or. z>=n .or. mod(q,2)/=0) then
         allocate(v(0,0)); if (present(status)) status=sgt_invalid_input; return
      end if
      call symmetric_eigen_jacobi(aw,values,vectors,local_status)
      h=q/2
      allocate(v(n,q))
      v(:,1:h)=vectors(:,1:h)
      v(:,h+1:q)=vectors(:,n-h+1:n)
      if (present(status)) status=local_status
   end subroutine bipartite_v_update

   subroutine laplacian_lambda_update(lb,ub,beta,u,lw,k,lambda,status)
      real(dp), intent(in) :: lb,ub,beta,u(:,:),lw(:,:)
      integer, intent(in) :: k
      real(dp), allocatable, intent(out) :: lambda(:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: d(:)
      integer :: q
      q=size(lw,1)-k
      if (q<1 .or. size(u,2)/=q .or. beta<=0.0_dp .or. ub<lb) then
         allocate(lambda(0)); if (present(status)) status=sgt_invalid_input; return
      end if
      allocate(d(q))
      d=diagonal_of_triple(u,lw)
      allocate(lambda(q))
      lambda=0.5_dp*(d+sqrt(d*d+4.0_dp/beta))
      lambda=min(max(lambda,lb),ub)
      ! Numerical eigensolvers already order d; enforce monotonicity harmlessly.
      block
         real(dp) :: fitted(q)
         call isotonic_increasing(lambda,fitted)
         lambda=min(max(fitted,lb),ub)
      end block
      if (present(status)) status=sgt_ok
   end subroutine laplacian_lambda_update

   subroutine bipartite_psi_update(v,aw,psi,status,lb,ub)
      real(dp), intent(in) :: v(:,:),aw(:,:)
      real(dp), allocatable, intent(out) :: psi(:)
      integer, intent(out), optional :: status
      real(dp), intent(in), optional :: lb,ub
      real(dp), allocatable :: c(:),ctilde(:),reversed(:),fit(:)
      real(dp) :: lower,upper
      integer :: n,h,i
      n=size(v,2)
      if (n<2 .or. mod(n,2)/=0) then
         allocate(psi(0)); if (present(status)) status=sgt_invalid_input; return
      end if
      lower=-huge(1.0_dp); upper=huge(1.0_dp)
      if (present(lb)) lower=lb
      if (present(ub)) upper=ub
      allocate(c(n))
      c=diagonal_of_triple(v,aw); h=n/2
      allocate(ctilde(h),reversed(h),fit(h),psi(n))
      do i=1,h
         ctilde(i)=0.5_dp*(c(n-i+1)-c(i))
         reversed(i)=ctilde(h-i+1)
      end do
      call isotonic_increasing(reversed,fit)
      do i=1,h
         psi(i)=-fit(h-i+1)
         psi(h+i)=fit(i)
      end do
      psi=min(max(psi,lower),upper)
      if (present(status)) status=sgt_ok
   end subroutine bipartite_psi_update
end module sgt_updates
