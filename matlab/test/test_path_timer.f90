! SPDX-License-Identifier: Artistic-2.0
! Derived from the R matlab package; see COPYRIGHTS and upstream/.
program test_path_timer
    use matlab, only : dp, fileparts, fileparts_result, fullfile, strcmp, tic, toc
    use test_support
    implicit none

    type(fileparts_result) :: p
    real(dp) :: elapsed
    integer :: i
    real(dp) :: sink

    p = fileparts('/home/luser/myfile.ext.gz')
    call assert_true(p%pathstr == '/home/luser', 'fileparts path')
    call assert_true(p%name == 'myfile.ext', 'fileparts name')
    call assert_true(p%ext == '.gz', 'fileparts extension')
    p = fileparts('.profile')
    call assert_true(p%name == '' .and. p%ext == '.profile', 'hidden file')
    call assert_true(fullfile('a', 'b', 'c') == 'a/b/c', 'fullfile')
    call assert_true(strcmp('foo', 'foo') .and. .not. strcmp('foo', 'bar'), 'strcmp')

    call tic()
    sink = 0.0_dp
    do i = 1, 100000
        sink = sink + sqrt(real(i, dp))
    end do
    elapsed = toc(.false.)
    call assert_true(elapsed >= 0.0_dp .and. sink > 0.0_dp, 'timer')

    write(*, '(a)') 'test_path_timer: PASS'
end program test_path_timer
