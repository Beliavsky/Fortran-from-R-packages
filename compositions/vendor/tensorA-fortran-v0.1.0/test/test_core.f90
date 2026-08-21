! SPDX-License-Identifier: GPL-2.0-or-later
program test_core
    use tensora
    implicit none
    type(tensor_t) :: a, b, c, d, tr, del, dg, rp, sl, bd, un
    real(dp) :: amat(2,3), bmat(3,2), expect(2,2)
    integer :: i

    amat = reshape([1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp,6.0_dp],[2,3])
    bmat = reshape([2.0_dp,0.0_dp,1.0_dp,3.0_dp,4.0_dp,5.0_dp],[3,2])
    a = tensor(reshape(amat,[6]), [2,3], ['a','b'])
    b = tensor(reshape(bmat,[6]), [3,2], ['b','c'])
    c = mul_tensor(a, ['b'], b, ['b'])
    expect = matmul(amat,bmat)
    call assert_shape(c,[2,2])
    call assert_close(c%data,cmplx(reshape(expect,[4]),0.0_dp,dp),1.0e-12_dp,'matrix contraction')

    d = reorder_tensor(c,['c'])
    call assert_axis(d,['c','a'])
    call assert_close(d%data,cmplx(reshape(transpose(expect),[4]),0.0_dp,dp),1.0e-12_dp,'reorder')

    tr = trace_tensor(delta_tensor([2,3],['a','b']), [1,2], [3,4])
    if (abs(real(tr%data(1),dp)-6.0_dp) > 1.0e-12_dp) error stop 'trace(delta)'

    del = delta_tensor([2,2],['a','b'])
    if (count(abs(del%data-cmplx(1.0_dp,0.0_dp,dp)) < 1.0e-14_dp) /= 4) error stop 'delta count'
    dg = diag_tensor(tensor([1.0_dp,2.0_dp,3.0_dp,4.0_dp],[2,2],['a','b']))
    if (abs(sum(real(dg%data,dp))-10.0_dp) > 1.0e-12_dp) error stop 'diag sum'

    rp = repeat_tensor(a,3,pos=2,name='k')
    call assert_shape(rp,[2,3,3])
    sl = slice_tensor(rp,2,[2],drop=.true.)
    call assert_shape(sl,[2,3])
    call assert_close(sl%data,a%data,1.0e-12_dp,'repeat/slice')

    bd = bind_tensor(a,1,a,1)
    call assert_shape(bd,[4,3])
    do i=1,3
        call assert_close(bd%data((i-1)*4+1:i*4), &
            [a%data((i-1)*2+1:i*2),a%data((i-1)*2+1:i*2)],1.0e-12_dp,'bind')
    end do

    un = untensor_tensor(tensor([(real(i,dp),i=1,24)],[2,3,4],['a','b','c']),[1,2],'ab',1)
    call assert_shape(un,[6,4])
    call assert_close(un%data,cmplx([(real(i,dp),i=1,24)],0.0_dp,dp),1.0e-12_dp,'untensor')

    print '(a)', 'test_core: PASS'
contains
    subroutine assert_shape(x,s)
        type(tensor_t),intent(in)::x
        integer,intent(in)::s(:)
        if(size(x%shape)/=size(s) .or. any(x%shape/=s)) error stop 'shape mismatch'
    end subroutine
    subroutine assert_axis(x,s)
        type(tensor_t),intent(in)::x
        character(len=*),intent(in)::s(:)
        integer::k
        if(size(x%axis)/=size(s)) error stop 'axis size mismatch'
        do k=1,size(s)
            if(trim(x%axis(k))/=trim(s(k))) error stop 'axis mismatch'
        end do
    end subroutine
    subroutine assert_close(x,y,tol,msg)
        complex(dp),intent(in)::x(:),y(:)
        real(dp),intent(in)::tol
        character(len=*),intent(in)::msg
        if(size(x)/=size(y)) error stop 'assert_close size'
        if(maxval(abs(x-y))>tol) then
            print *, msg, maxval(abs(x-y))
            error stop 'assert_close'
        end if
    end subroutine
end program test_core
