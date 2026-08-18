! SPDX-License-Identifier: GPL-2.0-or-later
program test_geometry
    use tensora
    implicit none
    type(tensor_t) :: x,g,y,z,d
    integer,allocatable :: pos(:,:)

    x=tensor([2.0_dp,6.0_dp],[2],['i'])
    g=tensor([2.0_dp,0.0_dp,0.0_dp,3.0_dp],[2,2],['i','j'])
    y=drag_tensor(x,g,['i'])
    if(trim(y%axis(1))/='^i') error stop 'drag axis'
    call close(y%data,cmplx([1.0_dp,2.0_dp],0.0_dp,dp),1.0e-10_dp,'drag cov to contra')
    z=drag_tensor(y,g,['^i'])
    if(trim(z%axis(1))/='i') error stop 'drag back axis'
    call close(z%data,x%data,1.0e-10_dp,'drag back')

    d=delta_tensor([2],['i'])
    z=riemann_pair(d,contraname_tensor(d))
    if(z%nelem()/=1 .or. abs(real(z%data(1),dp)-2.0_dp)>1.0e-12_dp) error stop 'riemann delta'

    pos=pos_tensor([2,3])
    if(any(pos(1,:)/=[1,1]) .or. any(pos(2,:)/=[2,1]) .or. any(pos(3,:)/=[1,2])) then
        error stop 'pos_tensor order'
    end if

    print '(a)', 'test_geometry: PASS'
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
end program test_geometry
