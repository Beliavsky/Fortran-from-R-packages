program test_special
   use rebayes_kinds, only : dp
   use rebayes_math, only : noncentral_t_pdf
   implicit none
   real(dp)::f
   f=noncentral_t_pdf(0.18964679350340852_dp,28.30332725682117_dp,1.7321530778591772_dp)
   if(abs(f-0.12070484951976014_dp)>5.0e-7_dp)then
      print *,f
      error stop "noncentral t reference"
   end if
   print *,"test_special: PASS"
end program test_special
