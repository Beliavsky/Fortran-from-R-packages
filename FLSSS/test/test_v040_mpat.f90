program test_v040_mpat
  use flsss_mod
  use flsss_mpat, only : search_md_i8_mpat
  use flsss_fast_search, only : search_md_i8_packed
  use flsss_search, only : search_md_i8
  implicit none

  call test_complete_solution_set()
  call test_public_mpat()
  call test_auto_central_target()
  call test_imposed_bounds()
  call test_randomized_equivalence()
  print *, 'test_v040_mpat: PASS'

contains

  subroutine test_complete_solution_set()
    integer(i8) :: v(15,6), target(6), me(6)
    integer :: i,j,k,idx(5)
    type(subset_solutions) :: ref,fast
    allocate(ref%sol(0),fast%sol(0))
    do j=1,6
      do i=1,15
        v(i,j)=int(i*(j+2)+i*i/(17+j)+3*j,i8)
      end do
    end do
    idx=[2,5,8,11,14]
    target=0_i8
    do k=1,5
      target=target+v(idx(k),:)
    end do
    me=1_i8
    call search_md_i8(v,5,target,me,ref,100000,10.0_dp)
    call search_md_i8_mpat(v,5,target,me,fast,100000,10.0_dp)
    call require_same_sets(ref,fast,'mPAT complete set')
    if(trim(fast%engine)/='mpat-packed') error stop 'mPAT engine tag'
    if(fast%tri_entries<=0_i8) error stop 'mPAT triangular cache'
    if(fast%bound_updates<=0_i8) error stop 'mPAT bound updates'
    if(fast%mpat_splits<=0_i8) error stop 'mPAT split count'
  end subroutine test_complete_solution_set

  subroutine test_public_mpat()
    real(dp) :: v(18,7),target(7),me(7)
    integer :: i,j,k,idx(5)
    type(integerized_search_result) :: a,b
    do j=1,7
      do i=1,18
        v(i,j)=real(i*(j+1)+i*i/(23+j),dp)/4.0_dp-3.0_dp
      end do
    end do
    idx=[2,6,9,13,17]
    target=0.0_dp
    do k=1,5
      target=target+v(idx(k),:)
    end do
    me=0.03_dp
    a=mflsss_par_integerized(5,v,target,me,solution_need=100000,engine='dfs')
    b=mflsss_par_integerized(5,v,target,me,solution_need=100000,engine='mpat')
    call require_same_sets(a%solution,b%solution,'public mPAT')
    if(trim(b%solution%engine)/='mpat-packed') error stop 'public mPAT engine'
  end subroutine test_public_mpat

  subroutine test_auto_central_target()
    real(dp) :: v(40,5),target(5),me(5)
    integer :: i,j,k,idx(6)
    type(integerized_search_result) :: r
    do j=1,5
      do i=1,40
        v(i,j)=real(i*(j+2)+i*i/(30+j),dp)
      end do
    end do
    idx=[8,13,18,23,28,33]
    target=0.0_dp
    do k=1,6
      target=target+v(idx(k),:)
    end do
    me=0.01_dp
    r=mflsss_par_integerized(6,v,target,me,solution_need=1,engine='auto')
    if(r%solution%size()/=1) error stop 'auto central solution'
    if(trim(r%solution%engine)/='mpat-packed') error stop 'auto did not select mPAT for central target'
  end subroutine test_auto_central_target

  subroutine test_imposed_bounds()
    real(dp) :: v(14,4),target(4),me(4)
    integer :: i,j,k,idx(4),lb(4),ub(4)
    type(integerized_search_result) :: a,b
    do j=1,4
      do i=1,14
        v(i,j)=real(i*(j+2)+i*i/(20+j),dp)
      end do
    end do
    idx=[3,6,9,12]
    target=0.0_dp
    do k=1,4
      target=target+v(idx(k),:)
    end do
    me=0.01_dp
    lb=[1,4,7,10]
    ub=[5,8,11,14]
    a=mflsss_par_impose_bounds_integerized(4,v,target,me,lb,ub,solution_need=10000,engine='dfs')
    b=mflsss_par_impose_bounds_integerized(4,v,target,me,lb,ub,solution_need=10000,engine='mpat')
    call require_same_sets(a%solution,b%solution,'bounded mPAT')
  end subroutine test_imposed_bounds

  subroutine test_randomized_equivalence()
    integer(i8) :: v(13,4),target(4),me(4)
    integer :: c,i,j,k,idx(4)
    type(subset_solutions) :: a,b
    allocate(a%sol(0),b%sol(0))
    do c=1,12
      do j=1,4
        do i=1,13
          v(i,j)=int((j+1)*i+(1+mod(c+j,3))*i*i/19+c*j,i8)
        end do
      end do
      idx=[1+mod(c,3),5+mod(c,2),9,12]
      target=0_i8
      do k=1,4
        target=target+v(idx(k),:)
      end do
      me=int(mod(c,2),i8)
      call search_md_i8(v,4,target,me,a,100000,10.0_dp)
      call search_md_i8_mpat(v,4,target,me,b,100000,10.0_dp)
      call require_same_sets(a,b,'randomized mPAT')
    end do
  end subroutine test_randomized_equivalence

  subroutine require_same_sets(a,b,label)
    type(subset_solutions),intent(in)::a,b
    character(len=*),intent(in)::label
    integer::i
    if(a%size()/=b%size())then
      print *,trim(label),a%size(),b%size()
      error stop 'solution-set size mismatch'
    end if
    do i=1,a%size()
      if(.not.has_solution(b,a%sol(i)%idx))then
        print *,trim(label),' missing ',a%sol(i)%idx
        error stop 'solution-set mismatch'
      end if
    end do
  end subroutine require_same_sets

  logical function has_solution(r,idx) result(found)
    type(subset_solutions),intent(in)::r
    integer,intent(in)::idx(:)
    integer::i
    found=.false.
    do i=1,r%size()
      if(size(r%sol(i)%idx)/=size(idx))cycle
      if(all(r%sol(i)%idx==idx))then
        found=.true.
        return
      end if
    end do
  end function has_solution

end program test_v040_mpat
