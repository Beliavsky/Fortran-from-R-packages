! SPDX-License-Identifier: GPL-2.0-or-later
module rootsolve_sparse
  use rootsolve_kinds, only : dp
  use rootsolve_types, only : steady_rhs
  use rootsolve_derivatives, only : perturb_value, jacobian_full
  use rootsolve_linalg, only : solve_linear, vector_norm2
  implicit none
  private
  public :: build_grid_pattern, discover_pattern, sparse_jacobian, solve_csr
contains

  subroutine build_grid_pattern(nspec, dims, ndims, cyclic, rowptr, colind)
    integer, intent(in) :: nspec, dims(3), ndims
    logical, intent(in) :: cyclic(3)
    integer, allocatable, intent(out) :: rowptr(:), colind(:)
    integer, allocatable :: tmp(:)
    integer :: nx, ny, nz, nc, n, ix, iy, iz, s, q, row, cell, c2, k, cap

    if (nspec < 1 .or. ndims < 1 .or. ndims > 3) error stop 'build_grid_pattern: invalid dimensions'
    nx = dims(1)
    ny = 1
    nz = 1
    if (ndims >= 2) ny = dims(2)
    if (ndims >= 3) nz = dims(3)
    if (nx < 1 .or. ny < 1 .or. nz < 1) error stop 'build_grid_pattern: zero grid dimension'
    nc = nx*ny*nz
    n = nspec*nc
    cap = n*(nspec+2*ndims)
    allocate(rowptr(n+1), tmp(max(1,cap)))
    k = 1
    rowptr(1) = 1
    do iz = 1, nz
      do iy = 1, ny
        do ix = 1, nx
          cell = cell_index(ix,iy,iz,nx,ny)
          do s = 1, nspec
            row = var_index(cell,s,nspec)
            do q = 1, nspec
              call push_unique(tmp,k,var_index(cell,q,nspec),rowptr(row))
            end do
            call add_neighbor(ix-1,iy,iz,1)
            call add_neighbor(ix+1,iy,iz,1)
            if (ndims >= 2) then
              call add_neighbor(ix,iy-1,iz,2)
              call add_neighbor(ix,iy+1,iz,2)
            end if
            if (ndims >= 3) then
              call add_neighbor(ix,iy,iz-1,3)
              call add_neighbor(ix,iy,iz+1,3)
            end if
            call sort_slice(tmp, rowptr(row), k-1)
            rowptr(row+1) = k
          end do
        end do
      end do
    end do
    allocate(colind(k-1))
    colind = tmp(:k-1)

  contains
    integer function cell_index(i,j,l,nx0,ny0) result(v)
      integer, intent(in) :: i,j,l,nx0,ny0
      v = i + nx0*((j-1)+ny0*(l-1))
    end function cell_index

    integer function var_index(c,s0,ns0) result(v)
      integer, intent(in) :: c,s0,ns0
      v = (c-1)*ns0+s0
    end function var_index

    subroutine add_neighbor(i,j,l,axis)
      integer, intent(in) :: i,j,l,axis
      integer :: ii,jj,ll
      ii=i; jj=j; ll=l
      if (axis == 1) then
        if (ii < 1) then
          if (.not. cyclic(1)) return
          ii=nx
        else if (ii > nx) then
          if (.not. cyclic(1)) return
          ii=1
        end if
      else if (axis == 2) then
        if (jj < 1) then
          if (.not. cyclic(2)) return
          jj=ny
        else if (jj > ny) then
          if (.not. cyclic(2)) return
          jj=1
        end if
      else
        if (ll < 1) then
          if (.not. cyclic(3)) return
          ll=nz
        else if (ll > nz) then
          if (.not. cyclic(3)) return
          ll=1
        end if
      end if
      c2=cell_index(ii,jj,ll,nx,ny)
      call push_unique(tmp,k,var_index(c2,s,nspec),rowptr(row))
    end subroutine add_neighbor
  end subroutine build_grid_pattern

  subroutine discover_pattern(func, t, y, pert, drop_tol, rowptr, colind)
    procedure(steady_rhs) :: func
    real(dp), intent(in) :: t, y(:), pert, drop_tol
    integer, allocatable, intent(out) :: rowptr(:), colind(:)
    real(dp), allocatable :: jac(:,:)
    integer, allocatable :: tmp(:)
    real(dp) :: scale
    integer :: n, i, j, k

    n = size(y)
    allocate(jac(n,n), rowptr(n+1), tmp(n*n))
    call jacobian_full(y, func, jac, time=t, pert=pert)
    k = 1
    rowptr(1) = 1
    do i = 1, n
      scale = max(1.0_dp, maxval(abs(jac(i,:))))
      do j = 1, n
        if (abs(jac(i,j)) > drop_tol*scale .or. i == j) then
          tmp(k)=j
          k=k+1
        end if
      end do
      rowptr(i+1)=k
    end do
    allocate(colind(k-1))
    colind=tmp(:k-1)
  end subroutine discover_pattern

  subroutine sparse_jacobian(func, t, y, f0, pert, rowptr, colind, values)
    procedure(steady_rhs) :: func
    real(dp), intent(in) :: t, y(:), f0(:), pert
    integer, intent(in) :: rowptr(:), colind(:)
    real(dp), intent(out) :: values(:)
    integer, allocatable :: color(:), cols(:)
    real(dp), allocatable :: yp(:), fp(:), delta(:)
    integer :: n, ncolor, c, j, i, p, nc

    n=size(y)
    if (size(values) /= size(colind)) error stop 'sparse_jacobian: shape mismatch'
    allocate(color(n), yp(n), fp(n), delta(n), cols(n))
    do j=1,n
      delta(j)=perturb_value(y(j),pert)
    end do
    call color_columns(rowptr,colind,n,color,ncolor)
    values=0.0_dp
    do c=1,ncolor
      yp=y
      nc=0
      do j=1,n
        if(color(j)==c)then
          yp(j)=yp(j)+delta(j)
          nc=nc+1
          cols(nc)=j
        end if
      end do
      call func(t,yp,fp)
      do i=1,n
        do p=rowptr(i),rowptr(i+1)-1
          j=colind(p)
          if(color(j)==c) values(p)=(fp(i)-f0(i))/delta(j)
        end do
      end do
    end do
  end subroutine sparse_jacobian

  subroutine solve_csr(rowptr, colind, values, b, x, info)
    integer, intent(in) :: rowptr(:), colind(:)
    real(dp), intent(in) :: values(:), b(:)
    real(dp), intent(out) :: x(:)
    integer, intent(out) :: info
    real(dp), allocatable :: dense(:,:)
    integer :: n, i, p

    n=size(b)
    if (n <= 400) then
      allocate(dense(n,n))
      dense=0.0_dp
      do i=1,n
        do p=rowptr(i),rowptr(i+1)-1
          dense(i,colind(p))=values(p)
        end do
      end do
      call solve_linear(dense,b,x,info)
    else
      call bicgstab(rowptr,colind,values,b,x,info)
    end if
  end subroutine solve_csr

  subroutine bicgstab(rowptr,colind,a,b,x,info)
    integer,intent(in)::rowptr(:),colind(:)
    real(dp),intent(in)::a(:),b(:)
    real(dp),intent(out)::x(:)
    integer,intent(out)::info
    real(dp),allocatable::r(:),r0(:),p(:),v(:),s(:),t(:),ph(:),sh(:),diag(:)
    real(dp)::rho,rho_old,alpha,omega,beta,tol,den
    integer::n,i,k,pos,maxit
    n=size(b)
    allocate(r(n),r0(n),p(n),v(n),s(n),t(n),ph(n),sh(n),diag(n))
    x=0.0_dp
    call csr_matvec(rowptr,colind,a,x,v)
    r=b-v
    r0=r
    p=0.0_dp
    v=0.0_dp
    rho_old=1.0_dp
    alpha=1.0_dp
    omega=1.0_dp
    diag=1.0_dp
    do i=1,n
      do pos=rowptr(i),rowptr(i+1)-1
        if(colind(pos)==i)then
          if(abs(a(pos))>sqrt(tiny(1.0_dp))) diag(i)=a(pos)
          exit
        end if
      end do
    end do
    tol=1.0e-10_dp*max(1.0_dp,vector_norm2(b))
    maxit=max(200,4*n)
    info=1
    do k=1,maxit
      rho=dot_product(r0,r)
      if(abs(rho)<=tiny(1.0_dp))exit
      if(k==1)then
        p=r
      else
        beta=(rho/rho_old)*(alpha/omega)
        p=r+beta*(p-omega*v)
      end if
      ph=p/diag
      call csr_matvec(rowptr,colind,a,ph,v)
      den=dot_product(r0,v)
      if(abs(den)<=tiny(1.0_dp))exit
      alpha=rho/den
      s=r-alpha*v
      if(vector_norm2(s)<=tol)then
        x=x+alpha*ph
        info=0
        return
      end if
      sh=s/diag
      call csr_matvec(rowptr,colind,a,sh,t)
      den=dot_product(t,t)
      if(den<=tiny(1.0_dp))exit
      omega=dot_product(t,s)/den
      x=x+alpha*ph+omega*sh
      r=s-omega*t
      if(vector_norm2(r)<=tol)then
        info=0
        return
      end if
      if(abs(omega)<=tiny(1.0_dp))exit
      rho_old=rho
    end do
  end subroutine bicgstab

  subroutine csr_matvec(rowptr,colind,a,x,y)
    integer,intent(in)::rowptr(:),colind(:)
    real(dp),intent(in)::a(:),x(:)
    real(dp),intent(out)::y(:)
    integer::i,p
    y=0.0_dp
    do i=1,size(y)
      do p=rowptr(i),rowptr(i+1)-1
        y(i)=y(i)+a(p)*x(colind(p))
      end do
    end do
  end subroutine csr_matvec

  subroutine color_columns(rowptr,colind,n,color,ncolor)
    integer,intent(in)::rowptr(:),colind(:),n
    integer,intent(out)::color(:),ncolor
    logical,allocatable::forbid(:)
    integer::j,k,c
    color=0
    ncolor=0
    allocate(forbid(n))
    do j=1,n
      forbid=.false.
      do k=1,j-1
        if(columns_conflict(j,k,rowptr,colind))then
          if(color(k)>0)forbid(color(k))=.true.
        end if
      end do
      c=1
      do while(c<=ncolor)
        if(.not.forbid(c))exit
        c=c+1
      end do
      if(c>ncolor)ncolor=c
      color(j)=c
    end do
  end subroutine color_columns

  logical function columns_conflict(c1,c2,rowptr,colind) result(conflict)
    integer,intent(in)::c1,c2,rowptr(:),colind(:)
    integer::i,p
    logical::a,b
    conflict=.false.
    do i=1,size(rowptr)-1
      a=.false.; b=.false.
      do p=rowptr(i),rowptr(i+1)-1
        if(colind(p)==c1)a=.true.
        if(colind(p)==c2)b=.true.
      end do
      if(a.and.b)then
        conflict=.true.
        return
      end if
    end do
  end function columns_conflict

  subroutine push_unique(a,k,v,start)
    integer,intent(inout)::a(:),k
    integer,intent(in)::v,start
    integer::i
    do i=start,k-1
      if(a(i)==v)return
    end do
    if(k>size(a))error stop 'build_grid_pattern: internal capacity exceeded'
    a(k)=v
    k=k+1
  end subroutine push_unique

  subroutine sort_slice(a,lo,hi)
    integer,intent(inout)::a(:)
    integer,intent(in)::lo,hi
    integer::i,j,v
    do i=lo+1,hi
      v=a(i);j=i-1
      do while(j>=lo)
        if(a(j)<=v)exit
        a(j+1)=a(j);j=j-1
      end do
      a(j+1)=v
    end do
  end subroutine sort_slice
end module rootsolve_sparse
