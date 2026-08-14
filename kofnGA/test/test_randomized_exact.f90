program test_randomized_exact
  use kofnga, only : dp, i64, kofnga_control, kofnga_result, kofn_ga
  implicit none
  integer, parameter :: n=12,k=4,ncase=20
  real(dp) :: values(n)
  type(kofnga_control) :: ctl
  type(kofnga_result) :: res
  integer :: c,i,j,t,exact(k),ord(n)

  ctl%popsize=100; ctl%keepbest=10; ctl%ngen=100; ctl%tourneysize=10; ctl%mutprob=0.04_dp
  do c=1,ncase
    do i=1,n
      values(i)=real(mod(7919*i+104729*c+37*i*c,10007),dp)/10007.0_dp + 1.0e-5_dp*real(i,dp)
    end do
    ord=[(i,i=1,n)]
    do i=2,n
      t=ord(i); j=i-1
      do while(j>=1)
        if(values(ord(j))<=values(t)) exit
        ord(j+1)=ord(j); j=j-1
      end do
      ord(j+1)=t
    end do
    exact=ord(1:k)
    call sort_indices(exact)
    ctl%seed=1000_i64+int(c,i64)
    call kofn_ga(n,k,obj,res,ctl)
    if(any(res%bestsol/=exact)) then
      write(*,*) 'case',c,'expected',exact,'got',res%bestsol
      error stop 'randomized exact optimum mismatch'
    end if
  end do
  print *, 'test_randomized_exact: PASS'
contains
  function obj(s) result(v)
    integer,intent(in)::s(:)
    real(dp)::v
    v=sum(values(s))
  end function obj
  subroutine sort_indices(x)
    integer,intent(inout)::x(:)
    integer::a,b,t
    do a=2,size(x)
      t=x(a); b=a-1
      do while(b>=1)
        if(x(b)<=t) exit
        x(b+1)=x(b); b=b-1
      end do
      x(b+1)=t
    end do
  end subroutine sort_indices
end program test_randomized_exact
