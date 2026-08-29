! SPDX-License-Identifier: GPL-3.0-only
module matrixdist_structure
   use r_compat, only: dp, runif1
   use matrixdist_types, only: ph_type, bivph_type
   implicit none
   private
   public :: random_structure, random_structure_bivph
contains
   function random_structure(p,structure,scale_factor) result(x)
      integer,intent(in)::p
      character(len=*),intent(in),optional::structure
      real(dp),intent(in),optional::scale_factor
      type(ph_type)::x
      character(len=32)::st
      real(dp)::sc,sm,u
      logical,allocatable::alegal(:),slegal(:,:)
      integer::i,j
      st='general'
      if(present(structure))st=trim(structure)
      sc=1.0_dp
      if(present(scale_factor))sc=scale_factor
      allocate(alegal(p),slegal(p,p),x%alpha(p),x%s(p,p))
      alegal=.false.
      slegal=.false.
      x%alpha=0.0_dp
      x%s=0.0_dp
      select case(trim(st))
      case('general')
      alegal=.true.
      slegal=.true.
      case('hyperexponential')
      alegal=.true.
      do i=1,p
      slegal(i,i)=.true.
      end do
      case('gerlang')
         alegal(1)=.true.
         do i=1,p-1
         slegal(i,i+1)=.true.
         end do
         slegal(p,p)=.true.
      case('coxian')
         alegal(1)=.true.
         do i=1,p-1
         slegal(i,i)=.true.
         slegal(i,i+1)=.true.
         end do
         slegal(p,p)=.true.
      case('gcoxian')
         alegal=.true.
         do i=1,p-1
         slegal(i,i)=.true.
         slegal(i,i+1)=.true.
         end do
         slegal(p,p)=.true.
      case default;error stop "random_structure: unknown structure"
      end select
      sm=0.0_dp
      do i=1,p
      if(alegal(i))then
      u=runif1()
      x%alpha(i)=u
      sm=sm+u
      end if
      end do
      if(sm>0.0_dp)x%alpha=x%alpha/sm
      do i=1,p
         do j=1,p
            if(i/=j .and. slegal(i,j))then
            u=runif1()
            x%s(i,j)=u
            x%s(i,i)=x%s(i,i)-u
            end if
         end do
      end do
      do i=1,p
      if(slegal(i,i))x%s(i,i)=x%s(i,i)-runif1()
      end do
      x%s=x%s*sc
   end function random_structure

   function random_structure_bivph(p1,p2,scale_factor) result(x)
      integer,intent(in)::p1,p2
      real(dp),intent(in),optional::scale_factor
      type(bivph_type)::x
      real(dp)::sc,sm,u
      integer::i,j
      sc=1.0_dp
      if(present(scale_factor))sc=scale_factor
      allocate(x%alpha(p1),x%s11(p1,p1),x%s12(p1,p2),x%s22(p2,p2))
      x%s11=0.0_dp
      x%s12=0.0_dp
      x%s22=0.0_dp
      sm=0.0_dp
      do i=1,p1
      u=runif1()
      x%alpha(i)=u
      sm=sm+u
      end do
      x%alpha=x%alpha/sm
      do i=1,p1
         do j=1,p1
         if(i/=j)then
         u=runif1()
         x%s11(i,j)=u
         x%s11(i,i)=x%s11(i,i)-u
         end if
         end do
         do j=1,p2
         u=runif1()
         x%s12(i,j)=u
         x%s11(i,i)=x%s11(i,i)-u
         end do
      end do
      do i=1,p2
         do j=1,p2
         if(i/=j)then
         u=runif1()
         x%s22(i,j)=u
         x%s22(i,i)=x%s22(i,i)-u
         end if
         end do
         x%s22(i,i)=x%s22(i,i)-runif1()
      end do
      x%s11=x%s11*sc
      x%s12=x%s12*sc
      x%s22=x%s22*sc
   end function random_structure_bivph
end module matrixdist_structure
