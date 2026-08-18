! SPDX-License-Identifier: GPL-2.0-only
program test_numeric
    use elliptic, only : dp, residue, integrate_segments
    implicit none
    integer::fails
    complex(dp)::r,v,pts(4)
    fails=0
    r=residue(fpole,(0.2_dp,0.1_dp),0.5_dp,nsub=2000)
    if(abs(r-exp((0.2_dp,0.1_dp)))>2e-9_dp)then;fails=fails+1;print *,'FAIL residue',r;end if
    pts=[(0.0_dp,0.0_dp),(1.0_dp,0.0_dp),(1.0_dp,1.0_dp),(0.0_dp,1.0_dp)]
    v=integrate_segments(fid,pts,.true.,400)
    if(abs(v)>2e-10_dp)then;fails=fails+1;print *,'FAIL contour',v;end if
    if(fails>0)error stop 1
    print *, 'test_numeric: PASS'
contains
    function fpole(z) result(v)
        complex(dp),intent(in)::z;complex(dp)::v
        v=exp(z)
    end function
    function fid(z) result(v)
        complex(dp),intent(in)::z;complex(dp)::v
        v=z*z
    end function
end program
