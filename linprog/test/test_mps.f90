program test_mps
    use linprog
    implicit none
    real(dp) :: c(3), b(3), a(3,3)
    character(len=8) :: vn(3), rn(3)
    type(mps_model) :: m
    integer :: u

    c = [1800.0_dp,600.0_dp,600.0_dp]
    b = [40.0_dp,90.0_dp,2500.0_dp]
    a = reshape([0.7_dp,1.5_dp,50.0_dp, 0.35_dp,1.0_dp,12.5_dp, &
                 0.0_dp,3.0_dp,20.0_dp], [3,3])
    vn = ['Cows    ','Bulls   ','Pigs    ']
    rn = ['Land    ','Stable  ','Labor   ']
    call writeMps('test_roundtrip.mps', c, b, a, 'Steinhauser', vn, rn)
    call readMps('test_roundtrip.mps', m, .true., .true.)
    if (trim(m%name) /= 'Steinhauser') error stop 'MPS name'
    if (maxval(abs(m%cvec-c)) > 1.0e-12_dp) error stop 'MPS cvec'
    if (maxval(abs(m%bvec-b)) > 1.0e-12_dp) error stop 'MPS bvec'
    if (maxval(abs(m%amat-a)) > 1.0e-12_dp) error stop 'MPS Amat'
    if (.not. m%has_result .or. abs(m%result%opt-93600.0_dp) > 1.0e-7_dp) then
        error stop 'MPS solve'
    end if

    open(newunit=u,file='test_bounds.mps',status='replace',action='write')
    write(u,'(a)') 'NAME          BoundTest'
    write(u,'(a)') 'ROWS'
    write(u,'(a)') ' N  obj'
    write(u,'(a)') ' G  minx'
    write(u,'(a)') 'COLUMNS'
    write(u,'(a)') '    x obj 1'
    write(u,'(a)') '    x minx 1'
    write(u,'(a)') 'RHS'
    write(u,'(a)') '    RHS minx 2'
    write(u,'(a)') 'BOUNDS'
    write(u,'(a)') ' UP BND x 3'
    write(u,'(a)') 'ENDATA'
    close(u)
    call readMps('test_bounds.mps', m, .false., .false.)
    if (size(m%bvec) /= 2) error stop 'MPS UP row count'
    if (maxval(abs(m%bvec-[-2.0_dp,3.0_dp])) > 1.0e-12_dp) error stop 'MPS G/UP rhs'
    if (maxval(abs(m%amat(:,1)-[-1.0_dp,1.0_dp])) > 1.0e-12_dp) error stop 'MPS G/UP A'
    print *, 'test_mps: PASS'
end program test_mps
