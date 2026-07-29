! SPDX-License-Identifier: GPL-3.0-or-later
module cccp_cones
   use cccp_kinds, only : dp
   use cccp_types, only : cone_constraint, cone_nnoc, cone_socc, cone_psdc
   use cccp_linalg, only : spd_inverse_logdet, symmetric_eigenvalues, vector_norm2
   implicit none
   private

   public :: nnoc, socc, psdc, nlfc
   public :: cone_barrier, cones_interior, cone_slacks, cone_duals, cone_dual_gradient
   public :: cone_identity, cone_jordan_product, cone_jordan_inverse
   public :: cone_inner_product, cone_norm, cone_max_step
   public :: minimum_cone_slack, phase_shift_needed, barrier_parameter

contains

   function nnoc(g, h) result(c)
      real(dp), intent(in) :: g(:,:), h(:)
      type(cone_constraint) :: c
      c%kind = cone_nnoc
      c%dim = size(h)
      allocate(c%g(size(g,1),size(g,2)), c%h(size(h)))
      c%g = g
      c%h = h
   end function nnoc

   function nlfc(g, h) result(c)
      real(dp), intent(in) :: g(:,:), h(:)
      type(cone_constraint) :: c
      c = nnoc(g, h)
   end function nlfc

   function socc(fmat, g, d, f) result(c)
      real(dp), intent(in) :: fmat(:,:), g(:), d(:), f
      type(cone_constraint) :: c
      integer :: m, n
      m = size(fmat,1)
      n = size(fmat,2)
      c%kind = cone_socc
      c%dim = m + 1
      allocate(c%g(m+1,n), c%h(m+1))
      c%g(1,:) = -d
      c%g(2:,:) = -fmat
      c%h(1) = f
      c%h(2:) = g
   end function socc

   function psdc(flist, f0) result(c)
      real(dp), intent(in) :: flist(:,:,:), f0(:,:)
      type(cone_constraint) :: c
      integer :: m, n, j
      m = size(f0,1)
      n = size(flist,3)
      c%kind = cone_psdc
      c%dim = m
      allocate(c%g(m*m,n), c%h(m*m))
      c%h = reshape(f0, [m*m])
      do j = 1, n
         c%g(:,j) = reshape(flist(:,:,j), [m*m])
      end do
   end function psdc

   pure integer function barrier_parameter(cones) result(nu)
      type(cone_constraint), intent(in) :: cones(:)
      integer :: k
      nu = 0
      do k = 1, size(cones)
         select case (cones(k)%kind)
         case (cone_nnoc)
            nu = nu + cones(k)%dim
         case (cone_socc)
            nu = nu + 2
         case (cone_psdc)
            nu = nu + cones(k)%dim
         end select
      end do
   end function barrier_parameter

   subroutine block_residual(c, x, s, phase, tau)
      type(cone_constraint), intent(in) :: c
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: s(:)
      logical, intent(in), optional :: phase
      real(dp), intent(in), optional :: tau
      logical :: use_phase
      real(dp) :: tv
      integer :: i, m

      s = c%h - matmul(c%g, x)
      use_phase = .false.
      if (present(phase)) use_phase = phase
      tv = 0.0_dp
      if (present(tau)) tv = tau
      if (.not. use_phase) return
      select case (c%kind)
      case (cone_nnoc)
         s = s + tv
      case (cone_socc)
         s(1) = s(1) + tv
      case (cone_psdc)
         m = c%dim
         do i = 1, m
            s(i + (i-1)*m) = s(i + (i-1)*m) + tv
         end do
      end select
   end subroutine block_residual

   logical function cones_interior(cones, x, phase, tau) result(ok)
      type(cone_constraint), intent(in) :: cones(:)
      real(dp), intent(in) :: x(:)
      logical, intent(in), optional :: phase
      real(dp), intent(in), optional :: tau
      real(dp), allocatable :: s(:), w(:)
      integer :: k, info, m
      logical :: ph

      ok = .true.
      ph = .false.
      if (present(phase)) ph = phase
      do k = 1, size(cones)
         allocate(s(size(cones(k)%h)))
         if (present(tau)) then
            call block_residual(cones(k), x, s, ph, tau)
         else
            call block_residual(cones(k), x, s, ph)
         end if
         select case (cones(k)%kind)
         case (cone_nnoc)
            if (minval(s) <= 0.0_dp) ok = .false.
         case (cone_socc)
            if (s(1) <= vector_norm2(s(2:))) ok = .false.
         case (cone_psdc)
            m = cones(k)%dim
            allocate(w(m))
            call symmetric_eigenvalues(reshape(s,[m,m]), w, info)
            if (info /= 0 .or. minval(w) <= 0.0_dp) ok = .false.
            deallocate(w)
         case default
            ok = .false.
         end select
         deallocate(s)
         if (.not. ok) exit
      end do
   end function cones_interior

   subroutine cone_barrier(cones, x, phi, grad, hess, ok, phase, tau)
      type(cone_constraint), intent(in) :: cones(:)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: phi
      real(dp), intent(out) :: grad(:)
      real(dp), intent(out) :: hess(:,:)
      logical, intent(out) :: ok
      logical, intent(in), optional :: phase
      real(dp), intent(in), optional :: tau
      real(dp), allocatable :: s(:), gs(:), hs(:,:), jmat(:,:), sinv(:,:)
      real(dp) :: delta, logdet, tv
      integer :: k, m, n, i, j, info
      logical :: ph

      n = size(x)
      phi = 0.0_dp
      grad = 0.0_dp
      hess = 0.0_dp
      ok = .true.
      ph = .false.
      if (present(phase)) ph = phase
      tv = 0.0_dp
      if (present(tau)) tv = tau

      do k = 1, size(cones)
         m = size(cones(k)%h)
         allocate(s(m), gs(m), hs(m,m), jmat(m,n))
         call block_residual(cones(k), x(1:size(cones(k)%g,2)), s, ph, tv)
         jmat = 0.0_dp
         jmat(:,1:size(cones(k)%g,2)) = -cones(k)%g
         if (ph .and. n > size(cones(k)%g,2)) then
            select case (cones(k)%kind)
            case (cone_nnoc)
               jmat(:,n) = 1.0_dp
            case (cone_socc)
               jmat(1,n) = 1.0_dp
            case (cone_psdc)
               do i = 1, cones(k)%dim
                  jmat(i+(i-1)*cones(k)%dim,n) = 1.0_dp
               end do
            end select
         end if

         select case (cones(k)%kind)
         case (cone_nnoc)
            if (minval(s) <= 0.0_dp) then
               ok = .false.
            else
               phi = phi - sum(log(s))
               gs = -1.0_dp / s
               hs = 0.0_dp
               do i = 1, m
                  hs(i,i) = 1.0_dp / (s(i)*s(i))
               end do
            end if
         case (cone_socc)
            delta = s(1)*s(1) - dot_product(s(2:),s(2:))
            if (delta <= 0.0_dp .or. s(1) <= 0.0_dp) then
               ok = .false.
            else
               phi = phi - log(delta)
               gs(1) = -2.0_dp*s(1)/delta
               gs(2:) = 2.0_dp*s(2:)/delta
               hs = 0.0_dp
               hs(1,1) = -2.0_dp/delta
               do i = 2, m
                  hs(i,i) = 2.0_dp/delta
               end do
               gs = -gs * delta
               hs = hs + matmul(reshape(gs,[m,1]),reshape(gs,[1,m]))/(delta*delta)
               gs = -gs / delta
            end if
         case (cone_psdc)
            m = cones(k)%dim
            allocate(sinv(m,m))
            call spd_inverse_logdet(reshape(s,[m,m]), sinv, logdet, info)
            if (info /= 0) then
               ok = .false.
            else
               phi = phi - logdet
               gs = -reshape(sinv,[m*m])
               hs = 0.0_dp
               do i = 1, m*m
                  do j = 1, m*m
                     hs(i,j) = sinv(1 + mod((i-1)/m, m), 1 + mod(j-1, m)) * &
                        sinv(1 + mod((j-1)/m, m), 1 + mod(i-1, m))
                  end do
               end do
            end if
            deallocate(sinv)
         case default
            ok = .false.
         end select

         if (.not. ok) then
            deallocate(s,gs,hs,jmat)
            return
         end if
         grad = grad + matmul(transpose(jmat), gs)
         hess = hess + matmul(transpose(jmat), matmul(hs, jmat))
         deallocate(s,gs,hs,jmat)
      end do
   end subroutine cone_barrier

   subroutine cone_slacks(cones, x, s, offsets)
      type(cone_constraint), intent(in) :: cones(:)
      real(dp), intent(in) :: x(:)
      real(dp), allocatable, intent(out) :: s(:)
      integer, allocatable, intent(out) :: offsets(:,:)
      integer :: k, pos, m, total
      total = 0
      do k = 1, size(cones)
         total = total + size(cones(k)%h)
      end do
      allocate(s(total), offsets(size(cones),2))
      pos = 1
      do k = 1, size(cones)
         m = size(cones(k)%h)
         offsets(k,:) = [pos, pos+m-1]
         call block_residual(cones(k), x, s(pos:pos+m-1))
         pos = pos + m
      end do
   end subroutine cone_slacks

   subroutine cone_duals(cones, x, tbar, z)
      type(cone_constraint), intent(in) :: cones(:)
      real(dp), intent(in) :: x(:), tbar
      real(dp), allocatable, intent(out) :: z(:)
      real(dp), allocatable :: s(:), sinv(:,:)
      real(dp) :: delta, logdet
      integer :: total, k, pos, m, info
      total = 0
      do k = 1, size(cones)
         total = total + size(cones(k)%h)
      end do
      allocate(z(total))
      pos = 1
      do k = 1, size(cones)
         m = size(cones(k)%h)
         allocate(s(m))
         call block_residual(cones(k),x,s)
         select case(cones(k)%kind)
         case(cone_nnoc)
            z(pos:pos+m-1) = 1.0_dp/(tbar*s)
         case(cone_socc)
            delta=s(1)*s(1)-dot_product(s(2:),s(2:))
            z(pos)=2.0_dp*s(1)/(tbar*delta)
            z(pos+1:pos+m-1)=-2.0_dp*s(2:)/(tbar*delta)
         case(cone_psdc)
            allocate(sinv(cones(k)%dim,cones(k)%dim))
            call spd_inverse_logdet(reshape(s,[cones(k)%dim,cones(k)%dim]),sinv,logdet,info)
            if(info==0) then
               z(pos:pos+m-1)=reshape(sinv,[m])/tbar
            else
               z(pos:pos+m-1)=0.0_dp
            end if
            deallocate(sinv)
         end select
         pos=pos+m
         deallocate(s)
      end do
   end subroutine cone_duals


   subroutine cone_dual_gradient(cones,z,grad)
      type(cone_constraint),intent(in)::cones(:)
      real(dp),intent(in)::z(:)
      real(dp),intent(out)::grad(:)
      integer::k,pos,m
      grad=0.0_dp;pos=1
      do k=1,size(cones)
         m=size(cones(k)%h)
         grad=grad+matmul(transpose(cones(k)%g),z(pos:pos+m-1))
         pos=pos+m
      end do
   end subroutine cone_dual_gradient

   real(dp) function minimum_cone_slack(cones,x) result(ans)
      type(cone_constraint), intent(in) :: cones(:)
      real(dp), intent(in) :: x(:)
      real(dp), allocatable :: s(:), w(:)
      integer :: k, m, info
      ans=huge(1.0_dp)
      do k=1,size(cones)
         allocate(s(size(cones(k)%h)))
         call block_residual(cones(k),x,s)
         select case(cones(k)%kind)
         case(cone_nnoc)
            ans=min(ans,minval(s))
         case(cone_socc)
            ans=min(ans,s(1)-vector_norm2(s(2:)))
         case(cone_psdc)
            m=cones(k)%dim; allocate(w(m))
            call symmetric_eigenvalues(reshape(s,[m,m]),w,info)
            if(info==0) ans=min(ans,minval(w))
            deallocate(w)
         end select
         deallocate(s)
      end do
      if(size(cones)==0) ans=huge(1.0_dp)
   end function minimum_cone_slack

   real(dp) function phase_shift_needed(cones,x) result(tau)
      type(cone_constraint), intent(in) :: cones(:)
      real(dp), intent(in) :: x(:)
      real(dp), allocatable :: s(:), w(:)
      integer :: k,m,info
      tau=1.0_dp
      do k=1,size(cones)
         allocate(s(size(cones(k)%h)))
         call block_residual(cones(k),x,s)
         select case(cones(k)%kind)
         case(cone_nnoc)
            tau=max(tau,1.0_dp-minval(s))
         case(cone_socc)
            tau=max(tau,1.0_dp+vector_norm2(s(2:))-s(1))
         case(cone_psdc)
            m=cones(k)%dim; allocate(w(m))
            call symmetric_eigenvalues(reshape(s,[m,m]),w,info)
            if(info==0) tau=max(tau,1.0_dp-minval(w))
            deallocate(w)
         end select
         deallocate(s)
      end do
   end function phase_shift_needed

   pure real(dp) function cone_inner_product(kind,s,z,dim) result(ans)
      integer,intent(in)::kind,dim
      real(dp),intent(in)::s(:),z(:)
      real(dp)::sm(dim,dim),zm(dim,dim)
      select case(kind)
      case(cone_nnoc,cone_socc)
         ans=dot_product(s,z)
      case(cone_psdc)
         sm=reshape(s,[dim,dim]); zm=reshape(z,[dim,dim])
         ans=sum(sm*zm)
      case default
         ans=0.0_dp
      end select
   end function cone_inner_product

   pure real(dp) function cone_norm(kind,s,dim) result(ans)
      integer,intent(in)::kind,dim
      real(dp),intent(in)::s(:)
      ans=sqrt(max(0.0_dp,cone_inner_product(kind,s,s,dim)))
   end function cone_norm

   pure subroutine cone_identity(kind,dim,e)
      integer,intent(in)::kind,dim
      real(dp),intent(out)::e(:)
      integer::i
      e=0.0_dp
      select case(kind)
      case(cone_nnoc)
         e=1.0_dp
      case(cone_socc)
         e(1)=1.0_dp
      case(cone_psdc)
         do i=1,dim
            e(i+(i-1)*dim)=1.0_dp
         end do
      end select
   end subroutine cone_identity

   pure subroutine cone_jordan_product(kind,dim,s,z,p)
      integer,intent(in)::kind,dim
      real(dp),intent(in)::s(:),z(:)
      real(dp),intent(out)::p(:)
      real(dp)::sm(dim,dim),zm(dim,dim),pm(dim,dim)
      select case(kind)
      case(cone_nnoc)
         p=s*z
      case(cone_socc)
         p(1)=dot_product(s,z)
         p(2:)=s(1)*z(2:)+z(1)*s(2:)
      case(cone_psdc)
         sm=reshape(s,[dim,dim]); zm=reshape(z,[dim,dim])
         pm=0.5_dp*(matmul(sm,zm)+matmul(zm,sm))
         p=reshape(pm,[dim*dim])
      end select
   end subroutine cone_jordan_product

   subroutine cone_jordan_inverse(kind,dim,s,z,p,info)
      integer,intent(in)::kind,dim
      real(dp),intent(in)::s(:),z(:)
      real(dp),intent(out)::p(:)
      integer,intent(out)::info
      real(dp)::aa,dd,sm(dim,dim),zm(dim,dim),ainv(dim,dim),logdet
      integer::i,j
      info=0
      select case(kind)
      case(cone_nnoc)
         if(any(abs(z)<=tiny(1.0_dp))) then; info=1; return; end if
         p=s/z
      case(cone_socc)
         aa=z(1)*z(1)-dot_product(z(2:),z(2:))
         if(aa<=0.0_dp) then; info=1; return; end if
         dd=dot_product(s(2:),z(2:))
         p(1)=(s(1)*z(1)-dd)/aa
         p(2:)=(aa/z(1)*s(2:)+(dd/z(1)-s(1))*z(2:))/aa
      case(cone_psdc)
         sm=reshape(s,[dim,dim]); zm=reshape(z,[dim,dim])
         call spd_inverse_logdet(zm,ainv,logdet,info)
         if(info/=0)return
         do i=1,dim
            do j=1,dim
               p(i+(j-1)*dim)=2.0_dp*sm(i,j)/(zm(i,i)+zm(j,j))
            end do
         end do
      case default
         info=1
      end select
   end subroutine cone_jordan_inverse

   real(dp) function cone_max_step(kind,dim,s) result(a)
      integer,intent(in)::kind,dim
      real(dp),intent(in)::s(:)
      real(dp),allocatable::w(:)
      integer::info
      select case(kind)
      case(cone_nnoc)
         a=-minval(s)
      case(cone_socc)
         a=vector_norm2(s(2:))-s(1)
      case(cone_psdc)
         allocate(w(dim)); call symmetric_eigenvalues(reshape(s,[dim,dim]),w,info)
         if(info==0) then; a=-minval(w); else; a=huge(1.0_dp); end if
         deallocate(w)
      case default
         a=huge(1.0_dp)
      end select
   end function cone_max_step

end module cccp_cones
