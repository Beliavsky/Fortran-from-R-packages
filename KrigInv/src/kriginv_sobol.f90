module kriginv_sobol
  use, intrinsic :: iso_fortran_env, only : int64
  use kriginv_kinds, only : dp
  use kriginv_sobol_data, only : sobol_initial, sobol_degree, sobol_a
  implicit none
  private
  integer(int64), parameter :: mask32=int(z'FFFFFFFF',int64)
  real(dp), parameter :: two32=4294967296.0_dp
  public :: sobol_points
contains
  function sobol_points(n,dim,start) result(x)
    integer, intent(in) :: n,dim
    integer(int64), intent(in), optional :: start
    real(dp), allocatable :: x(:,:)
    integer(int64), allocatable :: v(:,:)
    integer(int64) :: s0,idx,g,state,t
    integer :: i,j,k,maxbit,s,a
    if(n<0 .or. dim<1 .or. dim>1111) error stop 'kriginv: invalid Sobol dimensions'
    s0=1_int64; if(present(start)) s0=start
    if(n==0) then; allocate(x(0,dim)); return; end if
    idx=s0+int(n-1,int64); maxbit=1
    do while(shiftr(idx,maxbit)>0_int64 .and. maxbit<32)
      maxbit=maxbit+1
    end do
    allocate(v(maxbit,dim)); v=0_int64
    do i=1,maxbit; v(i,1)=shiftl(1_int64,32-i); end do
    do j=2,dim
      s=sobol_degree(j); a=sobol_a(j)
      do i=1,min(s,maxbit)
        v(i,j)=iand(shiftl(sobol_initial(i,j),32-i),mask32)
      end do
      do i=s+1,maxbit
        t=ieor(v(i-s,j),shiftr(v(i-s,j),s))
        do k=1,s-1
          if(btest(a,s-1-k)) t=ieor(t,v(i-k,j))
        end do
        v(i,j)=iand(t,mask32)
      end do
    end do
    allocate(x(n,dim))
    do i=1,n
      idx=s0+int(i-1,int64); g=ieor(idx,shiftr(idx,1))
      do j=1,dim
        state=0_int64
        do k=1,maxbit
          if(btest(g,k-1)) state=ieor(state,v(k,j))
        end do
        x(i,j)=real(iand(state,mask32),dp)/two32
      end do
    end do
  end function sobol_points
end module kriginv_sobol
