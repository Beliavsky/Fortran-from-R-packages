! SPDX-License-Identifier: MIT
module bekks_matrix
  use bekks_kinds, only: dp
  use bekks_linalg, only: general_inverse
  implicit none
  private
  public :: elimination_mat, commutation_mat, duplication_mat
  public :: diag_selection_mat, cut_mat_symmetric, cut_mat_asymmetric
  public :: vech_lower, unvech_lower, y_lag_cr, extract_csd

contains
  function vech_lower(a) result(v)
    real(dp), intent(in) :: a(:,:)
    real(dp) :: v(size(a,1)*(size(a,1)+1)/2)
    integer :: i,j,k,n
    n=size(a,1); k=0
    do j=1,n
      do i=j,n
        k=k+1; v(k)=a(i,j)
      end do
    end do
  end function vech_lower

  function unvech_lower(v,n) result(a)
    real(dp), intent(in) :: v(:)
    integer, intent(in) :: n
    real(dp) :: a(n,n)
    integer :: i,j,k
    a=0.0_dp; k=0
    do j=1,n
      do i=j,n
        k=k+1; a(i,j)=v(k)
      end do
    end do
  end function unvech_lower

  function elimination_mat(n) result(l)
    integer, intent(in) :: n
    real(dp) :: l(n*(n+1)/2,n*n)
    integer :: i,j,k,col
    l=0.0_dp; k=0
    do j=1,n
      do i=j,n
        k=k+1; col=i+(j-1)*n; l(k,col)=1.0_dp
      end do
    end do
  end function elimination_mat

  function commutation_mat(n) result(km)
    integer, intent(in) :: n
    real(dp) :: km(n*n,n*n)
    integer :: i,j
    km=0.0_dp
    do i=1,n
      do j=1,n
        km(i+n*(j-1),j+n*(i-1))=1.0_dp
      end do
    end do
  end function commutation_mat

  function duplication_mat(n) result(d)
    integer, intent(in) :: n
    real(dp) :: d(n*n,n*(n+1)/2)
    integer :: i,j,k
    d=0.0_dp; k=0
    do j=1,n
      do i=j,n
        k=k+1
        d(i+(j-1)*n,k)=1.0_dp
        if(i/=j)d(j+(i-1)*n,k)=1.0_dp
      end do
    end do
  end function duplication_mat

  function diag_selection_mat(n) result(s)
    integer, intent(in) :: n
    real(dp) :: s(n*n,n)
    integer :: i
    s=0.0_dp
    do i=1,n; s(i+(i-1)*n,i)=1.0_dp; end do
  end function diag_selection_mat

  function cut_mat_symmetric(n) result(c)
    integer, intent(in) :: n
    real(dp) :: c(n*(n+1)/2+2*n*n,n*(n+1)/2+2*n)
    integer :: nc,i,off_full,off_diag
    nc=n*(n+1)/2; c=0.0_dp
    do i=1,nc;c(i,i)=1.0_dp;end do
    off_full=nc; off_diag=nc
    do i=1,n
      c(off_full+i+(i-1)*n,off_diag+i)=1.0_dp
      c(off_full+n*n+i+(i-1)*n,off_diag+n+i)=1.0_dp
    end do
  end function cut_mat_symmetric

  function cut_mat_asymmetric(n) result(c)
    integer, intent(in) :: n
    real(dp) :: c(n*(n+1)/2+3*n*n,n*(n+1)/2+3*n)
    integer :: nc,i
    nc=n*(n+1)/2; c=0.0_dp
    do i=1,nc;c(i,i)=1.0_dp;end do
    do i=1,n
      c(nc+i+(i-1)*n,nc+i)=1.0_dp
      c(nc+n*n+i+(i-1)*n,nc+n+i)=1.0_dp
      c(nc+2*n*n+i+(i-1)*n,nc+2*n+i)=1.0_dp
    end do
  end function cut_mat_asymmetric

  function y_lag_cr(y,p) result(out)
    real(dp), intent(in) :: y(:,:)
    integer, intent(in) :: p
    real(dp), allocatable :: out(:,:)
    integer :: t,n,i,k,row
    t=size(y,1); n=size(y,2)
    if(p<1 .or. t<=p)then
      allocate(out(0,0)); return
    end if
    allocate(out(t-p,n*p+1)); out(:,1)=1.0_dp
    do row=1,t-p
      do k=1,p
        do i=1,n
          out(row,1+(k-1)*n+i)=y(row+p-k,i)
        end do
      end do
    end do
  end function y_lag_cr

  function extract_csd(h) result(csd)
    real(dp), intent(in) :: h(:,:,:)
    real(dp), allocatable :: csd(:,:)
    integer :: t,n,i
    t=size(h,3); n=size(h,1); allocate(csd(t,n))
    do i=1,n; csd(:,i)=sqrt(max(h(i,i,:),0.0_dp)); end do
  end function extract_csd
end module bekks_matrix
