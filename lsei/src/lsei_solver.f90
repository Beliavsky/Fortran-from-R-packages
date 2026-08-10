! SPDX-License-Identifier: GPL-2.0-or-later
module lsei_solver
   use lsei_kinds, only : dp
   use lsei_types
   use lsei_linalg
   use lsei_nnls, only : nnls_solve, pnnls_solve, ldp_solve
   implicit none
   private
   public :: lsi_solve, lsei_solve, qp_solve, pnnqp_solve, hfti_solve
contains
   subroutine lsi_solve(a,b,e,f,res,lower,upper,tol)
      real(dp), intent(in) :: a(:,:),b(:)
      real(dp), intent(in), optional :: e(:,:),f(:),lower(:),upper(:)
      type(ls_result), intent(out) :: res
      real(dp), intent(in), optional :: tol
      real(dp), allocatable :: c0(:,:),d0(:),e0(:,:),f0(:)
      integer :: n
      n=size(a,2); allocate(c0(0,n),d0(0))
      if(present(e)) then
         if(.not.present(f)) then; call invalid_result(res,n,size(a,1)); return; end if
         allocate(e0(size(e,1),n),f0(size(f))); e0=e; f0=f
      else
         allocate(e0(0,n),f0(0))
      end if
      call constrained_ls(a,b,c0,d0,e0,f0,res,lower,upper,tol)
   end subroutine lsi_solve

   subroutine lsei_solve(a,b,c,d,e,f,res,lower,upper,tol)
      real(dp), intent(in) :: a(:,:),b(:)
      real(dp), intent(in), optional :: c(:,:),d(:),e(:,:),f(:),lower(:),upper(:)
      type(ls_result), intent(out) :: res
      real(dp), intent(in), optional :: tol
      real(dp), allocatable :: c0(:,:),d0(:),e0(:,:),f0(:)
      integer :: n
      n=size(a,2)
      if(present(c)) then
         if(.not.present(d)) then; call invalid_result(res,n,size(a,1)); return; end if
         allocate(c0(size(c,1),n),d0(size(d))); c0=c; d0=d
      else
         allocate(c0(0,n),d0(0))
      end if
      if(present(e)) then
         if(.not.present(f)) then; call invalid_result(res,n,size(a,1)); return; end if
         allocate(e0(size(e,1),n),f0(size(f))); e0=e; f0=f
      else
         allocate(e0(0,n),f0(0))
      end if
      call constrained_ls(a,b,c0,d0,e0,f0,res,lower,upper,tol)
   end subroutine lsei_solve

   subroutine constrained_ls(a,b,c,d,e,f,res,lower,upper,tol)
      real(dp), intent(in) :: a(:,:),b(:),c(:,:),d(:),e(:,:),f(:)
      type(ls_result), intent(out) :: res
      real(dp), intent(in), optional :: lower(:),upper(:),tol
      real(dp), allocatable :: cp(:,:),x0(:),z(:,:),gg(:,:),hh(:),gr(:,:),hr(:)
      real(dp), allocatable :: ar(:,:),br(:),hz(:,:),dv(:),zr(:),x(:)
      type(ls_result) :: feas
      real(dp) :: eps
      integer :: m,n,rankc,st,q
      m=size(a,1); n=size(a,2); eps=1.0e-12_dp; if(present(tol)) eps=tol
      if(size(b)/=m .or. size(c,2)/=n .or. size(d)/=size(c,1) .or. &
         size(e,2)/=n .or. size(f)/=size(e,1)) then
         call invalid_result(res,n,m); return
      end if
      allocate(x0(n)); x0=0.0_dp
      if(size(c,1)>0) then
         allocate(cp(n,size(c,1))); call pseudoinverse(c,cp,rankc,1.0e-14_dp,st); x0=matmul(cp,d)
         if(maxval(abs(matmul(c,x0)-d))>eps*max(1.0_dp,maxval(abs(d)))) then
            call invalid_result(res,n,m); res%mode=LSEI_INFEASIBLE; res%x=x0; return
         end if
         call null_space(c,z,rankc,1.0e-14_dp,st)
      else
         rankc=0; z=identity_matrix(n)
      end if
      call augment_bounds(e,f,n,gg,hh,lower,upper)
      q=size(z,2)
      if(q==0) then
         allocate(res%x(n),res%residuals(m),res%dual(n),res%index(n)); res%x=x0
         res%residuals=b-matmul(a,res%x); res%rnorm=norm2(res%residuals); res%objective=res%rnorm**2
         res%dual=matmul(transpose(a),res%residuals); res%index=[(st,st=1,n)]; res%rank=rankc
         if(size(gg,1)>0 .and. minval(matmul(gg,res%x)-hh)<-eps) then
            res%mode=LSEI_INFEASIBLE
         else
            res%mode=LSEI_SUCCESS
         end if
         return
      end if
      allocate(gr(size(gg,1),q),hr(size(gg,1))); gr=matmul(gg,z); hr=hh-matmul(gg,x0)
      allocate(zr(q)); zr=0.0_dp
      if(size(gr,1)>0 .and. minval(-hr)<-eps) then
         call ldp_solve(gr,hr,feas,eps)
         if(.not.feas%succeeded()) then; call invalid_result(res,n,m); res%mode=LSEI_INFEASIBLE; return; end if
         zr=feas%x
      end if
      allocate(ar(m,q),br(m),hz(q,q),dv(q)); ar=matmul(a,z); br=b-matmul(a,x0)
      hz=matmul(transpose(ar),ar); dv=matmul(transpose(ar),br)
      call reduced_qp_active(hz,dv,gr,hr,zr,eps,st)
      allocate(x(n)); x=x0+matmul(z,zr)
      allocate(res%x(n),res%residuals(m),res%dual(n),res%index(n)); res%x=x
      where(abs(res%x)<100.0_dp*epsilon(1.0_dp)) res%x=0.0_dp
      res%residuals=b-matmul(a,res%x); res%rnorm=norm2(res%residuals); res%objective=res%rnorm**2
      res%dual=matmul(transpose(a),res%residuals); res%index=[(rankc,rankc=1,n)]; res%rank=matrix_rank(a,1.0e-14_dp)
      if(st==0) then; res%mode=LSEI_SUCCESS; else; res%mode=LSEI_NUMERICAL; end if
   end subroutine constrained_ls

   subroutine reduced_qp_active(h,d,g,b,z,tol,info)
      real(dp), intent(in) :: h(:,:),d(:),g(:,:),b(:),tol
      real(dp), intent(inout) :: z(:)
      integer, intent(out) :: info
      logical, allocatable :: active(:)
      real(dp), allocatable :: grad(:),p(:),kkt(:,:),rhs(:),sol(:),lambda(:)
      real(dp) :: alpha,slack,den,pnorm,scale,ridge
      integer :: n,m,na,i,j,k,block,drop,st,iter
      n=size(z); m=size(g,1); allocate(active(m),grad(n),p(n)); active=.false.; info=0
      if(m>0) then
         do i=1,m; if(abs(dot_product(g(i,:),z)-b(i))<=100.0_dp*tol) active(i)=.true.; end do
      end if
      scale=max(1.0_dp,maxval(abs(h))); ridge=100.0_dp*epsilon(1.0_dp)*scale
      do iter=1,max(100,20*(n+m+1))
         grad=matmul(h,z)-d; na=count(active)
         allocate(kkt(n+na,n+na),rhs(n+na),sol(n+na),lambda(na)); kkt=0.0_dp; rhs=0.0_dp
         kkt(1:n,1:n)=h
         do i=1,n; kkt(i,i)=kkt(i,i)+ridge; end do
         rhs(1:n)=-grad; k=0
         do i=1,m
            if(active(i)) then
               k=k+1; kkt(1:n,n+k)=-g(i,:); kkt(n+k,1:n)=g(i,:)
            end if
         end do
         call dense_solve(kkt,rhs,sol,st)
         if(st/=0) then; info=st; return; end if
         p=sol(1:n); if(na>0) lambda=sol(n+1:n+na)
         deallocate(kkt,rhs,sol)
         pnorm=norm2(p)
         if(pnorm<=tol*(1.0_dp+norm2(z))) then
            drop=0; k=0
            do i=1,m
               if(active(i)) then
                  k=k+1
                  if(lambda(k)<-100.0_dp*tol) then
                     if(drop==0) then; drop=i
                     else
                        j=count(active(1:i-1))
                        if(lambda(k)<lambda(j)) drop=i
                     end if
                  end if
               end if
            end do
            deallocate(lambda)
            if(drop==0) return
            active(drop)=.false.; cycle
         end if
         deallocate(lambda)
         alpha=1.0_dp; block=0
         do i=1,m
            if(active(i)) cycle
            den=dot_product(g(i,:),p)
            if(den < -tol) then
               slack=dot_product(g(i,:),z)-b(i)
               if(slack/(-den)<alpha) then; alpha=max(0.0_dp,slack/(-den)); block=i; end if
            end if
         end do
         z=z+alpha*p
         if(block>0 .and. alpha<1.0_dp-10.0_dp*epsilon(1.0_dp)) active(block)=.true.
      end do
      info=1
   end subroutine reduced_qp_active

   subroutine qp_solve(q,p,res,c,d,e,f,lower,upper,tol)
      real(dp), intent(in) :: q(:,:),p(:)
      type(ls_result), intent(out) :: res
      real(dp), intent(in), optional :: c(:,:),d(:),e(:,:),f(:),lower(:),upper(:),tol
      real(dp), allocatable :: ev(:),v(:,:),a(:,:),bb(:)
      real(dp) :: eps
      integer :: n,st,kr,i
      n=size(q,1); eps=1.0e-15_dp; if(present(tol)) eps=tol
      if(size(q,2)/=n .or. size(p)/=n) then; call invalid_result(res,n,0); return; end if
      allocate(ev(n),v(n,n)); call symmetric_eigen(q,ev,v,st)
      if(ev(1)<=0.0_dp) then; call invalid_result(res,n,0); res%mode=LSEI_NUMERICAL; return; end if
      kr=count(ev>=ev(1)*eps); allocate(a(kr,n),bb(kr)); a=0.0_dp; bb=0.0_dp
      do i=1,kr
         a(i,:)=sqrt(max(ev(i),0.0_dp))*v(:,i)
         bb(i)=-dot_product(p,v(:,i))/sqrt(max(ev(i),tiny(1.0_dp)))
      end do
      call lsei_solve(a,bb,c,d,e,f,res,lower,upper)
      if(allocated(res%x)) res%objective=0.5_dp*dot_product(res%x,matmul(q,res%x))+dot_product(p,res%x)
   end subroutine qp_solve

   subroutine pnnqp_solve(q,p,kfree,res,sum_value,tol)
      real(dp), intent(in) :: q(:,:),p(:)
      integer, intent(in) :: kfree
      type(ls_result), intent(out) :: res
      real(dp), intent(in), optional :: sum_value,tol
      real(dp), allocatable :: ev(:),v(:,:),a(:,:),bb(:)
      real(dp) :: eps
      integer :: n,st,kr,i
      n=size(q,1); eps=1.0e-20_dp; if(present(tol)) eps=tol
      allocate(ev(n),v(n,n)); call symmetric_eigen(q,ev,v,st); kr=count(ev>=ev(1)*eps)
      allocate(a(kr,n),bb(kr))
      do i=1,kr
         a(i,:)=sqrt(max(ev(i),0.0_dp))*v(:,i)
         bb(i)=-dot_product(p,v(:,i))/sqrt(max(ev(i),tiny(1.0_dp)))
      end do
      if(present(sum_value)) then; call pnnls_solve(a,bb,kfree,res,sum_value); else; call pnnls_solve(a,bb,kfree,res); end if
      if(allocated(res%x)) res%objective=0.5_dp*dot_product(res%x,matmul(q,res%x))+dot_product(p,res%x)
   end subroutine pnnqp_solve

   subroutine hfti_solve(a,b,res,tol)
      real(dp), intent(in) :: a(:,:),b(:,:)
      type(hfti_result), intent(out) :: res
      real(dp), intent(in), optional :: tol
      real(dp), allocatable :: x(:,:),rn(:)
      integer, allocatable :: piv(:)
      integer :: n
      n=size(a,2); allocate(x(n,size(b,2)),rn(size(b,2)),piv(n))
      if(size(b,1)/=size(a,1)) then; res%mode=LSEI_BAD_DIMENSIONS; return; end if
      if (present(tol)) then
         call least_squares_pivoted(a, b, x, res%krank, rn, piv, tol)
      else
         call least_squares_pivoted(a, b, x, res%krank, rn, piv)
      end if
      allocate(res%x(n,size(b,2)),res%transformed_b(max(size(a,1),n),size(b,2)),res%pivot(n),res%rnorm(size(b,2)))
      res%x=x; res%pivot=piv; res%rnorm=rn; res%transformed_b=0.0_dp
      res%transformed_b(1:n,:)=x; res%mode=LSEI_SUCCESS
   end subroutine hfti_solve

   subroutine augment_bounds(e,f,n,gg,hh,lower,upper)
      real(dp), intent(in) :: e(:,:),f(:)
      integer, intent(in) :: n
      real(dp), allocatable, intent(out) :: gg(:,:),hh(:)
      real(dp), intent(in), optional :: lower(:),upper(:)
      integer :: nl,nu,m,i,k
      m=size(e,1); nl=0; nu=0
      if(present(lower)) then
         if(size(lower)==1) then; if(abs(lower(1))<huge(1.0_dp)/10.0_dp) nl=n
         else; nl=count(abs(lower)<huge(1.0_dp)/10.0_dp); end if
      end if
      if(present(upper)) then
         if(size(upper)==1) then; if(abs(upper(1))<huge(1.0_dp)/10.0_dp) nu=n
         else; nu=count(abs(upper)<huge(1.0_dp)/10.0_dp); end if
      end if
      allocate(gg(m+nl+nu,n),hh(m+nl+nu)); gg=0.0_dp; hh=0.0_dp
      if(m>0) then; gg(1:m,:)=e; hh(1:m)=f; end if; k=m
      if(present(lower)) then
         do i=1,n
            if(size(lower)==1) then
               if(abs(lower(1))>=huge(1.0_dp)/10.0_dp) cycle; k=k+1; gg(k,i)=1.0_dp; hh(k)=lower(1)
            else
               if(abs(lower(i))>=huge(1.0_dp)/10.0_dp) cycle; k=k+1; gg(k,i)=1.0_dp; hh(k)=lower(i)
            end if
         end do
      end if
      if(present(upper)) then
         do i=1,n
            if(size(upper)==1) then
               if(abs(upper(1))>=huge(1.0_dp)/10.0_dp) cycle; k=k+1; gg(k,i)=-1.0_dp; hh(k)=-upper(1)
            else
               if(abs(upper(i))>=huge(1.0_dp)/10.0_dp) cycle; k=k+1; gg(k,i)=-1.0_dp; hh(k)=-upper(i)
            end if
         end do
      end if
   end subroutine augment_bounds

   subroutine invalid_result(res,n,m)
      type(ls_result), intent(out) :: res
      integer, intent(in) :: n,m
      allocate(res%x(max(0,n)),res%residuals(max(0,m)),res%dual(max(0,n)),res%index(max(0,n)))
      res%x=0.0_dp; res%residuals=0.0_dp; res%dual=0.0_dp; res%index=0; res%mode=LSEI_BAD_DIMENSIONS
   end subroutine invalid_result
end module lsei_solver
