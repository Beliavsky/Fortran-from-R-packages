! Derived from locfit src/lf_fitfun.c, GPL-2-or-later.
module locfit_basis
  use locfit_kinds, only : dp
  use locfit_constants, only : ksph,kprod,klm,kzeon,stangl,mxdeg
  implicit none
  private
  public :: basis_size, polynomial_basis, derivative_basis

contains

  pure integer function basis_size(dim, degree, kernel_type) result(p)
    integer,intent(in)::dim,degree,kernel_type
    integer::i
    select case(kernel_type)
    case(ksph)
      p=1
      do i=1,degree
        p=p*(dim+i)/i
      end do
    case(kprod)
      p=dim*degree+1
    case(klm)
      p=dim
    case(kzeon)
      p=1
    case default
      p=0
    end select
  end function basis_size

  pure subroutine angular_terms(dx,scale,degree,nder,ff)
    real(dp),intent(in)::dx,scale
    integer,intent(in)::degree,nder
    real(dp),intent(out)::ff(0:mxdeg)
    ff=0.0_dp
    select case(nder)
    case(0)
      ff(0)=1.0_dp
      if(degree>=1)ff(1)=sin(dx/scale)*scale
      if(degree>=2)ff(2)=(1.0_dp-cos(dx/scale))*scale*scale
    case(1)
      if(degree>=1)ff(1)=cos(dx/scale)
      if(degree>=2)ff(2)=sin(dx/scale)*scale
    case(2)
      if(degree>=1)ff(1)=-sin(dx/scale)/scale
      if(degree>=2)ff(2)=cos(dx/scale)
    end select
  end subroutine angular_terms

  pure subroutine derivative_basis(x,target,degree,kernel_type,style,scale,deriv,f)
    real(dp),intent(in)::x(:),target(:),scale(:)
    integer,intent(in)::degree,kernel_type
    integer,intent(in)::style(:)
    integer,intent(in)::deriv(:) ! 1-based variable indices, repeated for higher derivatives
    real(dp),intent(out)::f(:)
    integer::d,nd,m,i,j,k
    integer::ct(size(x))
    real(dp)::ff(0:mxdeg,size(x)),dx(size(x))
    d=size(x); nd=size(deriv); m=0; f=0.0_dp
    if(kernel_type==kzeon)then; f(1)=1.0_dp; return; end if
    if(kernel_type==klm)then; f(1:d)=x; return; end if
    m=m+1; f(m)=merge(1.0_dp,0.0_dp,nd==0)
    if(degree==0)return
    ct=0; dx=x-target
    do i=1,nd
      if(deriv(i)>=1 .and. deriv(i)<=d)ct(deriv(i))=ct(deriv(i))+1
    end do
    ff=0.0_dp
    do i=1,d
      if(style(i)==stangl)then
        call angular_terms(dx(i),scale(i),degree,ct(i),ff(:,i))
      else
        if(ct(i)<=degree)then
          ff(ct(i),i)=1.0_dp
          do j=ct(i)+1,degree
            ff(j,i)=ff(j-1,i)*dx(i)/real(j-ct(i),dp)
          end do
        end if
      end if
    end do
    if(d==1 .or. kernel_type==kprod)then
      do j=1,degree
        do i=1,d
          m=m+1
          if(ct(i)==nd)f(m)=ff(j,i)
        end do
      end do
      return
    end if
    ! Full spherical polynomial basis, matching locfit ordering through degree 3.
    do i=1,d
      m=m+1; if(ct(i)==nd)f(m)=ff(1,i)
    end do
    if(degree==1)return
    do i=1,d
      m=m+1; if(ct(i)==nd)f(m)=ff(2,i)
      do j=i+1,d
        m=m+1; if(ct(i)+ct(j)==nd)f(m)=ff(1,i)*ff(1,j)
      end do
    end do
    if(degree==2)return
    do i=1,d
      m=m+1; if(ct(i)==nd)f(m)=ff(3,i)
      do k=i+1,d
        m=m+1; if(ct(i)+ct(k)==nd)f(m)=ff(2,i)*ff(1,k)
      end do
      do j=i+1,d
        m=m+1; if(ct(i)+ct(j)==nd)f(m)=ff(1,i)*ff(2,j)
        do k=j+1,d
          m=m+1
          if(ct(i)+ct(j)+ct(k)==nd)f(m)=ff(1,i)*ff(1,j)*ff(1,k)
        end do
      end do
    end do
  end subroutine derivative_basis

  pure subroutine polynomial_basis(x,target,degree,kernel_type,style,scale,f)
    real(dp),intent(in)::x(:),target(:),scale(:)
    integer,intent(in)::degree,kernel_type,style(:)
    real(dp),intent(out)::f(:)
    integer,allocatable::empty(:)
    allocate(empty(0))
    call derivative_basis(x,target,degree,kernel_type,style,scale,empty,f)
  end subroutine polynomial_basis

end module locfit_basis
