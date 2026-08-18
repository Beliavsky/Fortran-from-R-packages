! SPDX-License-Identifier: GPL-2.0-or-later
program test_linalg
    use tensora
    implicit none
    type(tensor_t) :: a, ai, ident, b, x, s, ch, pwr
    type(tensor_svd_t) :: sv
    real(dp) :: amat(3,3), bmat(3,2), spd(3,3)
    complex(dp), allocatable :: u(:,:), v(:,:), rec(:,:)
    real(dp), allocatable :: sval(:)
    integer :: k

    amat = reshape([4.0_dp,1.0_dp,2.0_dp, 1.0_dp,3.0_dp,0.5_dp, 2.0_dp,0.5_dp,5.0_dp],[3,3])
    a = tensor(reshape(amat,[9]),[3,3],['row','col'])
    ai = inv_tensor(a,['row'])
    ident = mul_tensor(a,['col'],ai,['col'])
    call assert_identity(reshape(ident%data,[3,3]),1.0e-10_dp,'inverse')

    bmat = reshape([1.0_dp,2.0_dp,3.0_dp,4.0_dp,0.0_dp,2.0_dp],[3,2])
    b = tensor(reshape(bmat,[6]),[3,2],['row','rhs'])
    x = solve_tensor(a,b,['row'],['row'])
    ident = mul_tensor(a,['col'],x,['col'])
    call assert_close(ident%data,b%data,1.0e-10_dp,'solve')

    sv = svd_tensor(a,['row'],['col'])
    allocate(u(3,3),v(3,3),sval(3),rec(3,3))
    u=reshape(sv%u%data,[3,3])
    v=reshape(sv%v%data,[3,3])
    sval=real(sv%d%data,dp)
    rec=(0.0_dp,0.0_dp)
    do k=1,3
        rec=rec+sval(k)*spread(u(:,k),2,3)*spread(conjg(v(:,k)),1,3)
    end do
    call assert_close(reshape(rec,[9]),cmplx(reshape(amat,[9]),0.0_dp,dp),1.0e-10_dp,'svd')

    spd = matmul(transpose(amat),amat) + 0.5_dp*identity3()
    s = tensor(reshape(spd,[9]),[3,3],['i','j'])
    ch = chol_tensor(s,['i'],['j'])
    ident = mul_tensor(ch,['lambda'],mark_tensor(ch,mark="'"),["lambda'"])
    call assert_close(ident%data,s%data,1.0e-9_dp,'cholesky')

    pwr = power_tensor(s,['i'],['j'],0.5_dp)
    ident = mul_tensor(pwr,['j'],pwr,['i'])
    call assert_close(ident%data,s%data,2.0e-8_dp,'power')

    print '(a)', 'test_linalg: PASS'
contains
    function identity3() result(q)
        real(dp)::q(3,3)
        integer::i
        q=0.0_dp
        do i=1,3
            q(i,i)=1.0_dp
        end do
    end function
    subroutine assert_identity(m,tol,msg)
        complex(dp),intent(in)::m(:,:)
        real(dp),intent(in)::tol
        character(len=*),intent(in)::msg
        complex(dp)::q(size(m,1),size(m,2))
        integer::i
        q=(0.0_dp,0.0_dp)
        do i=1,min(size(q,1),size(q,2))
            q(i,i)=(1.0_dp,0.0_dp)
        end do
        if(maxval(abs(m-q))>tol) then
            print *,msg,maxval(abs(m-q))
            error stop 'identity mismatch'
        end if
    end subroutine
    subroutine assert_close(x,y,tol,msg)
        complex(dp),intent(in)::x(:),y(:)
        real(dp),intent(in)::tol
        character(len=*),intent(in)::msg
        if(maxval(abs(x-y))>tol) then
            print *,msg,maxval(abs(x-y))
            error stop 'close mismatch'
        end if
    end subroutine
end program test_linalg
