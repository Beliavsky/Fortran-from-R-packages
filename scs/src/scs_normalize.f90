! SPDX-License-Identifier: GPL-3.0-only
! Ruiz/L2 scaling translated from upstream SCS.
module scs_normalize
   use scs_kinds, only : dp
   use scs_types, only : scs_matrix, scs_cone, scs_scaling, scs_solution
   use scs_linalg, only : norm_inf, norm_sq, safe_div_pos
   use scs_cones, only : enforce_cone_boundaries
   implicit none
   private
   public :: normalize_a_p, normalize_b_c, normalize_sol, unnormalize_sol
   public :: unnormalize_primal, unnormalize_dual
   real(dp), parameter :: min_norm_factor=1.0e-4_dp, max_norm_factor=1.0e4_dp
contains
   pure real(dp) function apply_limit(x) result(v)
      real(dp),intent(in)::x
      v=x
      if(v<min_norm_factor)v=1.0_dp
      if(v>max_norm_factor)v=max_norm_factor
   end function apply_limit

   subroutine normalize_a_p(P,has_p,A,k,scal)
      type(scs_matrix),intent(inout)::P,A
      logical,intent(in)::has_p
      type(scs_cone),intent(in)::k
      type(scs_scaling),intent(out)::scal
      real(dp),allocatable::dt(:),et(:)
      integer::pass
      allocate(scal%D(A%m),scal%E(A%n),dt(A%m),et(A%n))
      scal%D=1.0_dp;scal%E=1.0_dp;scal%primal_scale=1.0_dp;scal%dual_scale=1.0_dp
      do pass=1,25
         call compute_ruiz(P,has_p,A,k,dt,et);call rescale(P,has_p,A,dt,et,scal)
      end do
      call compute_l2(P,has_p,A,k,dt,et);call rescale(P,has_p,A,dt,et,scal)
   end subroutine normalize_a_p

   subroutine compute_ruiz(P,has_p,A,k,dt,et)
      type(scs_matrix),intent(in)::P,A
      logical,intent(in)::has_p
      type(scs_cone),intent(in)::k
      real(dp),intent(out)::dt(:),et(:)
      integer::c,j,r
      real(dp)::w,nm
      dt=0.0_dp
      do c=1,A%n
         do j=A%p(c),A%p(c+1)-1
            r=A%i(j);dt(r)=max(dt(r),abs(A%x(j)))
         end do
      end do
      call enforce_cone_boundaries(k,dt,.false.)
      dt=1.0_dp/sqrt([(apply_limit(dt(j)),j=1,size(dt))])
      et=0.0_dp
      if(has_p)then
         do c=1,P%n
            do j=P%p(c),P%p(c+1)-1
               r=P%i(j);w=abs(P%x(j));et(c)=max(et(c),w);if(r/=c)et(r)=max(et(r),w)
            end do
         end do
      end if
      do c=1,A%n
         if(A%p(c+1)>A%p(c))then;nm=maxval(abs(A%x(A%p(c):A%p(c+1)-1)));else;nm=0.0_dp;end if
         et(c)=1.0_dp/sqrt(apply_limit(max(et(c),nm)))
      end do
   end subroutine compute_ruiz

   subroutine compute_l2(P,has_p,A,k,dt,et)
      type(scs_matrix),intent(in)::P,A
      logical,intent(in)::has_p
      type(scs_cone),intent(in)::k
      real(dp),intent(out)::dt(:),et(:)
      integer::c,j,r
      real(dp)::w
      dt=0.0_dp
      do c=1,A%n
         do j=A%p(c),A%p(c+1)-1
            r=A%i(j);dt(r)=dt(r)+A%x(j)*A%x(j)
         end do
      end do
      dt=sqrt(dt);call enforce_cone_boundaries(k,dt,.true.)
      do j=1,size(dt);dt(j)=1.0_dp/sqrt(apply_limit(dt(j)));end do
      et=0.0_dp
      if(has_p)then
         do c=1,P%n
            do j=P%p(c),P%p(c+1)-1
               r=P%i(j);w=P%x(j)*P%x(j);et(c)=et(c)+w;if(r/=c)et(r)=et(r)+w
            end do
         end do
      end if
      do c=1,A%n
         if(A%p(c+1)>A%p(c))et(c)=et(c)+norm_sq(A%x(A%p(c):A%p(c+1)-1))
         et(c)=1.0_dp/sqrt(apply_limit(sqrt(et(c))))
      end do
   end subroutine compute_l2

   subroutine rescale(P,has_p,A,dt,et,scal)
      type(scs_matrix),intent(inout)::P,A
      logical,intent(in)::has_p
      real(dp),intent(in)::dt(:),et(:)
      type(scs_scaling),intent(inout)::scal
      integer::c,j,r
      do c=1,A%n
         do j=A%p(c),A%p(c+1)-1;r=A%i(j);A%x(j)=A%x(j)*dt(r)*et(c);end do
      end do
      if(has_p)then
         do c=1,P%n
            do j=P%p(c),P%p(c+1)-1;r=P%i(j);P%x(j)=P%x(j)*et(r)*et(c);end do
         end do
      end if
      scal%D=scal%D*dt;scal%E=scal%E*et
   end subroutine rescale

   subroutine normalize_b_c(scal,b,c)
      type(scs_scaling),intent(inout)::scal
      real(dp),intent(inout)::b(:),c(:)
      real(dp)::sigma,nmb,nmc
      c=c*scal%E;b=b*scal%D;nmc=norm_inf(c);nmb=norm_inf(b);sigma=max(nmc,nmb)
      if(sigma<min_norm_factor)sigma=1.0_dp
      if(sigma>max_norm_factor)sigma=max_norm_factor
      sigma=safe_div_pos(1.0_dp,sigma);c=c*sigma;b=b*sigma;scal%primal_scale=sigma;scal%dual_scale=sigma
   end subroutine normalize_b_c

   subroutine normalize_sol(scal,sol)
      type(scs_scaling),intent(in)::scal;type(scs_solution),intent(inout)::sol
      sol%x=sol%x/(scal%E/scal%dual_scale);sol%y=sol%y/(scal%D/scal%primal_scale);sol%s=sol%s*(scal%D*scal%dual_scale)
   end subroutine normalize_sol
   subroutine unnormalize_sol(scal,sol)
      type(scs_scaling),intent(in)::scal;type(scs_solution),intent(inout)::sol
      sol%x=sol%x*(scal%E/scal%dual_scale);sol%y=sol%y*(scal%D/scal%primal_scale);sol%s=sol%s/(scal%D*scal%dual_scale)
   end subroutine unnormalize_sol
   subroutine unnormalize_primal(scal,r)
      type(scs_scaling),intent(in)::scal;real(dp),intent(inout)::r(:);r=r/(scal%D*scal%dual_scale)
   end subroutine unnormalize_primal
   subroutine unnormalize_dual(scal,r)
      type(scs_scaling),intent(in)::scal;real(dp),intent(inout)::r(:);r=r/(scal%E*scal%primal_scale)
   end subroutine unnormalize_dual
end module scs_normalize
