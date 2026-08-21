! SPDX-License-Identifier: GPL-2.0-or-later
program test_named_algebra
    use tensora
    implicit none
    type(tensor_t) :: a,b,c,x,y,z,r,br,dg,e
    real(dp) :: am(2,3), bm(3,4), expect(2,4)
    complex(dp) :: xc(2,2), yc(2,2), zc(2,2)
    integer :: i

    am=reshape([(real(i,dp),i=1,6)],[2,3])
    bm=reshape([(real(i,dp)/3.0_dp,i=1,12)],[3,4])
    a=tensor(reshape(am,[6]),[2,3],['i','j'])
    b=tensor(reshape(bm,[12]),[3,4],['j','k'])
    c=einstein_pair(a,b)
    expect=matmul(am,bm)
    call close(c%data,cmplx(reshape(expect,[8]),0.0_dp,dp),1.0e-12_dp,'einstein')

    x=tensor([1.0_dp,2.0_dp],[2],['i'])
    y=tensor([10.0_dp,20.0_dp,30.0_dp],[3],['j'])
    r=x+y
    if(any(r%shape/=[2,3])) error stop 'broadcast shape'
    call close(r%data,cmplx([11.0_dp,12.0_dp,21.0_dp,22.0_dp,31.0_dp,32.0_dp],0.0_dp,dp), &
        1.0e-12_dp,'broadcast')

    dg=tensor([2.0_dp,3.0_dp],[2],['i'])
    e=diagmul_tensor(a,['i'],dg,['i'])
    call close(e%data,cmplx([2.0_dp,6.0_dp,6.0_dp,12.0_dp,10.0_dp,18.0_dp],0.0_dp,dp),1.0e-12_dp,'diagmul')

    xc=reshape([cmplx(1.0_dp,1.0_dp,dp),cmplx(2.0_dp,-1.0_dp,dp), &
        cmplx(0.5_dp,0.2_dp,dp),cmplx(3.0_dp,0.0_dp,dp)],[2,2])
    yc=reshape([cmplx(2.0_dp,0.0_dp,dp),cmplx(1.0_dp,0.5_dp,dp), &
        cmplx(-1.0_dp,0.3_dp,dp),cmplx(0.5_dp,-0.5_dp,dp)],[2,2])
    x=tensor(reshape(xc,[4]),[2,2],['i','j'])
    y=tensor(reshape(yc,[4]),[2,2],[character(len=2) :: '^j','k'])
    r=riemann_pair(x,y)
    zc=matmul(xc,yc)
    call close(r%data,reshape(zc,[4]),1.0e-12_dp,'riemann complex')

    br=mark_tensor(a,'*',['j'])
    if(trim(br%axis(2))/='j*') error stop 'mark'
    if(trim(contraname('^alpha'))/='alpha') error stop 'contraname 1'
    if(trim(contraname('alpha'))/='^alpha') error stop 'contraname 2'

    z = 2.0_dp*x
    call close(z%data, 2.0_dp*x%data, 1.0e-12_dp, 'real scalar multiply')
    z = z/2.0_dp
    call close(z%data, x%data, 1.0e-12_dp, 'real scalar divide')

    print '(a)', 'test_named_algebra: PASS'
contains
    subroutine close(v,w,tol,msg)
      complex(dp),intent(in)::v(:),w(:)
      real(dp),intent(in)::tol
      character(len=*),intent(in)::msg
      if(maxval(abs(v-w))>tol) then
        print *,msg,maxval(abs(v-w))
        error stop 'mismatch'
      end if
    end subroutine
end program test_named_algebra
