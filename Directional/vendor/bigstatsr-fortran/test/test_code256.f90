program test_code256
    use bigstatsr
    implicit none
    type(fbm_code256) :: x
    integer :: c1(4),c2(4),c3(4),map(0:255)
    integer, allocatable :: cr(:,:),cc(:,:)
    map=0
    map(0)=1; map(1)=2; map(2)=3
    x=create_code256('test_code256.bk',4,3)
    c1=[0,1,2,1]; c2=[1,1,0,2]; c3=[2,1,2,0]
    call x%write_col(1,c1); call x%write_col(2,c2); call x%write_col(3,c3)
    cr=big_counts_rows(x,map)
    cc=big_counts_cols(x,map)
    call check(all(cr(:,1)==[1,1,1]),'row counts')
    call check(all(cc(:,1)==[1,2,1]),'column counts')
    call execute_command_line('rm -f test_code256.bk')
    print *, 'test_code256: PASS'
contains
    subroutine check(ok,msg)
        logical,intent(in)::ok
        character(len=*),intent(in)::msg
        if(.not.ok) then
            print *, 'FAIL: ',trim(msg)
            error stop 1
        end if
    end subroutine check
end program test_code256
