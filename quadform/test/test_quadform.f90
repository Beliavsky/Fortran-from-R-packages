program test_quadform
    use quadform, only : dp, ht, cprod, tcprod, quad_form, quad_form_inv, quad_form_chol, &
        quad3_form, quad3_form_ab, quad3_form_bc, quad3_form_inv, quad3_tform, quad3_tform_ab, &
        quad3_tform_bc, quad_tform, quad_tform_inv, quad_diag, quad_tdiag, quad3_diag, quad3_tdiag, &
        quad_trace, quad_ttrace, qf, q3, qt, qd, qtr
    implicit none

    call test_real()
    call test_complex()
    call test_inverse_errors()
    print '(a)', 'test_quadform: PASS'

contains

    subroutine test_real()
        real(dp) :: m(2,2), x(2,3), x1(2,3), y(3,2), y1(3,2)
        real(dp) :: l(2,2), spd(2,2)
        real(dp), allocatable :: a(:,:), b(:,:)
        real(dp) :: d(3), td(3)
        integer :: info

        m = reshape([2.0_dp, -0.5_dp, 1.25_dp, 3.0_dp], [2,2])
        x = reshape([1.0_dp,2.0_dp, -1.0_dp,0.5_dp, 2.0_dp,-3.0_dp], [2,3])
        x1 = reshape([0.5_dp,-1.0_dp, 2.0_dp,1.5_dp, -0.5_dp,4.0_dp], [2,3])
        y = transpose(x)
        y1 = transpose(x1)

        call assert_mat(ht(x), transpose(x), 1.0e-12_dp)
        call assert_mat(cprod(x), matmul(transpose(x),x), 1.0e-12_dp)
        call assert_mat(cprod(x,x1), matmul(transpose(x),x1), 1.0e-12_dp)
        call assert_mat(tcprod(y,y1), matmul(y,transpose(y1)), 1.0e-12_dp)
        call assert_mat(quad_form(m,x), matmul(transpose(x),matmul(m,x)), 1.0e-12_dp)
        call assert_mat(qf(m,x), quad_form(m,x), 1.0e-12_dp)
        call assert_mat(quad3_form(m,x,x1), matmul(transpose(x),matmul(m,x1)), 1.0e-12_dp)
        call assert_mat(quad3_form_ab(m,x,x1), quad3_form(m,x,x1), 1.0e-12_dp)
        call assert_mat(quad3_form_bc(m,x,x1), quad3_form(m,x,x1), 1.0e-12_dp)
        call assert_mat(q3(m,x,x1), quad3_form(m,x,x1), 1.0e-12_dp)
        call assert_mat(quad3_tform(m,y,y1), matmul(y,matmul(m,transpose(y1))), 1.0e-12_dp)
        call assert_mat(quad3_tform_ab(m,y,y1), quad3_tform(m,y,y1), 1.0e-12_dp)
        call assert_mat(quad3_tform_bc(m,y,y1), quad3_tform(m,y,y1), 1.0e-12_dp)
        call assert_mat(quad_tform(m,y), matmul(y,matmul(m,transpose(y))), 1.0e-12_dp)
        call assert_mat(qt(m,y), quad_tform(m,y), 1.0e-12_dp)

        d = quad_diag(m,x)
        call assert_vec(d, diag_real(quad_form(m,x)), 1.0e-12_dp)
        call assert_vec(qd(m,x), d, 1.0e-12_dp)
        td = quad_tdiag(m,y)
        call assert_vec(td, diag_real(quad_tform(m,y)), 1.0e-12_dp)
        call assert_vec(quad3_diag(m,x,x1), diag_real(quad3_form(m,x,x1)), 1.0e-12_dp)
        call assert_vec(quad3_tdiag(m,y,y1), diag_real(quad3_tform(m,y,y1)), 1.0e-12_dp)
        call assert_scalar(quad_trace(m,x), sum(d), 1.0e-12_dp)
        call assert_scalar(qtr(m,x), quad_trace(m,x), 1.0e-12_dp)
        call assert_scalar(quad_ttrace(m,y), sum(td), 1.0e-12_dp)

        spd = reshape([4.0_dp,1.0_dp,1.0_dp,3.0_dp],[2,2])
        l = 0.0_dp
        l(1,1) = 2.0_dp
        l(2,1) = 0.5_dp
        l(2,2) = sqrt(2.75_dp)
        call assert_mat(quad_form_chol(l,x), quad_form(spd,x), 1.0e-12_dp)

        a = quad_form_inv(spd,x,info)
        if (info /= 0) error stop 'real quad_form_inv returned failure'
        b = matmul(transpose(x), solve2_real(spd,x))
        call assert_mat(a,b,1.0e-11_dp)
        a = quad3_form_inv(spd,x,x1,info)
        if (info /= 0) error stop 'real quad3_form_inv returned failure'
        b = matmul(transpose(x), solve2_real(spd,x1))
        call assert_mat(a,b,1.0e-11_dp)
        call assert_mat(quad_tform_inv(spd,y,info), matmul(y,solve2_real(spd,transpose(y))), 1.0e-11_dp)
        if (info /= 0) error stop 'real quad_tform_inv returned failure'
    end subroutine test_real

    subroutine test_complex()
        complex(dp) :: m(2,2), x(2,3), x1(2,3), y(3,2), y1(3,2)
        complex(dp), allocatable :: a(:,:), b(:,:)
        integer :: info

        m(1,1)=cmplx(2.0_dp,0.5_dp,dp); m(2,1)=cmplx(-0.5_dp,0.2_dp,dp)
        m(1,2)=cmplx(1.25_dp,-0.3_dp,dp); m(2,2)=cmplx(3.0_dp,0.1_dp,dp)
        x(1,:)=[cmplx(1.0_dp,0.5_dp,dp),cmplx(-1.0_dp,0.2_dp,dp),cmplx(2.0_dp,-0.5_dp,dp)]
        x(2,:)=[cmplx(2.0_dp,-0.3_dp,dp),cmplx(0.5_dp,1.0_dp,dp),cmplx(-3.0_dp,0.7_dp,dp)]
        x1(1,:)=[cmplx(0.5_dp,0.1_dp,dp),cmplx(2.0_dp,-0.2_dp,dp),cmplx(-0.5_dp,0.6_dp,dp)]
        x1(2,:)=[cmplx(-1.0_dp,0.4_dp,dp),cmplx(1.5_dp,0.3_dp,dp),cmplx(4.0_dp,-0.8_dp,dp)]
        y = transpose(x)
        y1 = transpose(x1)

        call assert_cmat(ht(x), transpose(conjg(x)), 1.0e-12_dp)
        call assert_cmat(cprod(x), matmul(transpose(conjg(x)),x), 1.0e-12_dp)
        call assert_cmat(cprod(x,x1), matmul(transpose(conjg(x)),x1), 1.0e-12_dp)
        call assert_cmat(tcprod(y,y1), matmul(y,transpose(conjg(y1))), 1.0e-12_dp)
        call assert_cmat(quad_form(m,x), matmul(transpose(conjg(x)),matmul(m,x)), 1.0e-12_dp)
        call assert_cmat(quad3_form(m,x,x1), matmul(transpose(conjg(x)),matmul(m,x1)), 1.0e-12_dp)
        call assert_cmat(quad3_form_ab(m,x,x1), quad3_form(m,x,x1), 1.0e-12_dp)
        call assert_cmat(quad3_form_bc(m,x,x1), quad3_form(m,x,x1), 1.0e-12_dp)
        call assert_cmat(quad3_tform(m,y,y1), matmul(y,matmul(m,transpose(conjg(y1)))), 1.0e-12_dp)
        call assert_cmat(quad_tform(m,y), matmul(y,matmul(m,transpose(conjg(y)))), 1.0e-12_dp)
        call assert_cvec(quad_diag(m,x), diag_complex(quad_form(m,x)), 1.0e-12_dp)
        call assert_cvec(quad_tdiag(m,y), diag_complex(quad_tform(m,y)), 1.0e-12_dp)
        call assert_cvec(quad3_diag(m,x,x1), diag_complex(quad3_form(m,x,x1)), 1.0e-12_dp)
        call assert_cvec(quad3_tdiag(m,y,y1), diag_complex(quad3_tform(m,y,y1)), 1.0e-12_dp)
        call assert_cscalar(quad_trace(m,x), sum(quad_diag(m,x)), 1.0e-12_dp)
        call assert_cscalar(quad_ttrace(m,y), sum(quad_tdiag(m,y)), 1.0e-12_dp)

        a = quad_form_inv(m,x,info)
        if (info /= 0) error stop 'complex quad_form_inv returned failure'
        b = matmul(transpose(conjg(x)),solve2_complex(m,x))
        call assert_cmat(a,b,1.0e-10_dp)
        a = quad3_form_inv(m,x,x1,info)
        if (info /= 0) error stop 'complex quad3_form_inv returned failure'
        b = matmul(transpose(conjg(x)),solve2_complex(m,x1))
        call assert_cmat(a,b,1.0e-10_dp)
    end subroutine test_complex

    subroutine test_inverse_errors()
        real(dp) :: m(2,2), x(2,1)
        real(dp), allocatable :: z(:,:)
        integer :: info
        m = reshape([1.0_dp,2.0_dp,2.0_dp,4.0_dp],[2,2])
        x(:,1) = [1.0_dp,2.0_dp]
        z = quad_form_inv(m,x,info)
        if (info == 0) error stop 'singular matrix not detected'
        if (size(z) /= 0) error stop 'failure result should be empty'
    end subroutine test_inverse_errors

    function diag_real(a) result(d)
        real(dp), intent(in) :: a(:,:)
        real(dp) :: d(min(size(a,1),size(a,2)))
        integer :: i
        do i=1,size(d); d(i)=a(i,i); end do
    end function diag_real

    function diag_complex(a) result(d)
        complex(dp), intent(in) :: a(:,:)
        complex(dp) :: d(min(size(a,1),size(a,2)))
        integer :: i
        do i=1,size(d); d(i)=a(i,i); end do
    end function diag_complex

    function solve2_real(a,b) result(xout)
        real(dp), intent(in) :: a(2,2),b(:,:)
        real(dp) :: xout(2,size(b,2)),det
        det=a(1,1)*a(2,2)-a(1,2)*a(2,1)
        xout(1,:)=(a(2,2)*b(1,:)-a(1,2)*b(2,:))/det
        xout(2,:)=(-a(2,1)*b(1,:)+a(1,1)*b(2,:))/det
    end function solve2_real

    function solve2_complex(a,b) result(xout)
        complex(dp), intent(in) :: a(2,2),b(:,:)
        complex(dp) :: xout(2,size(b,2)),det
        det=a(1,1)*a(2,2)-a(1,2)*a(2,1)
        xout(1,:)=(a(2,2)*b(1,:)-a(1,2)*b(2,:))/det
        xout(2,:)=(-a(2,1)*b(1,:)+a(1,1)*b(2,:))/det
    end function solve2_complex


    subroutine assert_mat(a,b,tol)
        real(dp), intent(in) :: a(:,:),b(:,:),tol
        if (any(shape(a)/=shape(b)) .or. maxval(abs(a-b))>tol) error stop 'real matrix assertion failed'
    end subroutine assert_mat
    subroutine assert_vec(a,b,tol)
        real(dp), intent(in) :: a(:),b(:),tol
        if (size(a)/=size(b) .or. maxval(abs(a-b))>tol) error stop 'real vector assertion failed'
    end subroutine assert_vec
    subroutine assert_scalar(a,b,tol)
        real(dp), intent(in) :: a,b,tol
        if (abs(a-b)>tol) error stop 'real scalar assertion failed'
    end subroutine assert_scalar
    subroutine assert_cmat(a,b,tol)
        complex(dp), intent(in) :: a(:,:),b(:,:)
        real(dp), intent(in) :: tol
        if (any(shape(a)/=shape(b)) .or. maxval(abs(a-b))>tol) error stop 'complex matrix assertion failed'
    end subroutine assert_cmat
    subroutine assert_cvec(a,b,tol)
        complex(dp), intent(in) :: a(:),b(:)
        real(dp), intent(in) :: tol
        if (size(a)/=size(b) .or. maxval(abs(a-b))>tol) error stop 'complex vector assertion failed'
    end subroutine assert_cvec
    subroutine assert_cscalar(a,b,tol)
        complex(dp), intent(in) :: a,b
        real(dp), intent(in) :: tol
        if (abs(a-b)>tol) error stop 'complex scalar assertion failed'
    end subroutine assert_cscalar

end program test_quadform
