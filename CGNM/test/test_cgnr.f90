program test_cgnr
   use cgnm, only : dp,cgnr_ata_atb,cgnr_ata_atb_reg
   implicit none
   real(dp) :: a(4,2),b(4),x(2),xr(2)
   a=reshape([1._dp,0._dp,1._dp,2._dp, 0._dp,1._dp,1._dp,-1._dp],[4,2])
   b=matmul(a,[2._dp,-1._dp])
   call cgnr_ata_atb(a,b,x)
   if(maxval(abs(x-[2._dp,-1._dp]))>1.e-10_dp) error stop 'cgnr'
   call cgnr_ata_atb_reg(a,b,1._dp,xr)
   if(any(abs(xr)>abs(x)+1.e-12_dp)) error stop 'regularization'
end program test_cgnr
