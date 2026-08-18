! SPDX-License-Identifier: GPL-2.0-or-later
program test_batch_stats
    use tensora
    implicit none
    type(tensor_t) :: x,y,z,mu,v,a,ai,b,sol,chk,pinv,rec
    real(dp) :: xd(12), yd(12), aa(8), bb(4)
    integer :: k

    xd=[1.0_dp,2.0_dp,3.0_dp, 2.0_dp,4.0_dp,6.0_dp, &
        3.0_dp,6.0_dp,9.0_dp, 4.0_dp,8.0_dp,12.0_dp]
    x=tensor(xd,[3,2,2],[character(len=5) :: 'obs','var','batch'])
    mu=mean_tensor(x,['obs'])
    call close(mu%data,cmplx([2.0_dp,4.0_dp,6.0_dp,8.0_dp],0.0_dp,dp),1.0e-12_dp,'mean')
    v=var_tensor(x,['obs'],['batch'])
    if(any(v%shape/=[2,2,2])) error stop 'var shape'
    if(abs(real(v%data(1),dp)-1.0_dp)>1.0e-12_dp) error stop 'var batch1'
    if(abs(real(v%data(5),dp)-9.0_dp)>1.0e-12_dp) error stop 'var batch2'

    yd=[1.0_dp,0.0_dp, 0.0_dp,1.0_dp, 2.0_dp,1.0_dp, 1.0_dp,3.0_dp, &
        2.0_dp,0.0_dp, 0.0_dp,2.0_dp]
    y=tensor(yd,[2,2,3],['i','j','k'])
    z=mul_tensor(y,['j'],y,['i'],['k'])
    do k=1,3
      call matmul_slice_check(y,z,k)
    end do

    aa=[2.0_dp,0.0_dp,0.0_dp,3.0_dp, 4.0_dp,1.0_dp,1.0_dp,2.0_dp]
    a=tensor(aa,[2,2,2],['i','j','k'])
    ai=inv_tensor(a,['i'],['k'])
    z=mul_tensor(a,['j'],ai,['j'],['k'])
    do k=1,2
      call identity_slice(z,k)
    end do

    bb=[1.0_dp,2.0_dp,3.0_dp,4.0_dp]
    b=tensor(bb,[2,2],[character(len=3) :: 'i','rhs'])
    sol=solve_tensor(a,b,['i'],['i'],['k'])
    chk=mul_tensor(a,['j'],sol,['j'],['k'])
    do k=1,2
      call close(chk%data((k-1)*4+1:k*4),b%data,1.0e-10_dp,'batch solve')
    end do

    a=tensor([1.0_dp,2.0_dp,3.0_dp, 2.0_dp,4.0_dp,6.0_dp],[3,2],['i','j'])
    pinv=inv_tensor(a,['i'],allow_singular=.true.)
    rec=mul_tensor(a,['j'],pinv,['j'])
    rec=mul_tensor_pos(rec,[2],a,[1])
    call close(rec%data,a%data,1.0e-8_dp,'pseudoinverse identity')

    print '(a)', 'test_batch_stats: PASS'
contains
    subroutine matmul_slice_check(t,r,k)
      type(tensor_t),intent(in)::t,r
      integer,intent(in)::k
      complex(dp)::m(2,2),e(2,2),got(2,2)
      m=reshape(t%data((k-1)*4+1:k*4),[2,2])
      e=matmul(m,m)
      got=reshape(r%data((k-1)*4+1:k*4),[2,2])
      if(maxval(abs(e-got))>1.0e-12_dp) error stop 'batched mul'
    end subroutine
    subroutine identity_slice(t,k)
      type(tensor_t),intent(in)::t
      integer,intent(in)::k
      complex(dp)::m(2,2),q(2,2)
      q=(0.0_dp,0.0_dp)
      q(1,1)=1
      q(2,2)=1
      m=reshape(t%data((k-1)*4+1:k*4),[2,2])
      if(maxval(abs(m-q))>1.0e-10_dp) error stop 'batch inverse'
    end subroutine
    subroutine close(v,w,tol,msg)
      complex(dp),intent(in)::v(:),w(:)
      real(dp),intent(in)::tol
      character(len=*),intent(in)::msg
      if(maxval(abs(v-w))>tol) then
        print *,msg,maxval(abs(v-w))
        error stop 'mismatch'
      end if
    end subroutine
end program test_batch_stats
