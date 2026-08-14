program test_utils
   use rebayes_kinds, only : dp
   use rebayes_math, only : huber_eps, huber_pdf
   use rebayes_utils
   implicit none
   real(dp)::x(3),f(3),q(3),out(3),smooth(3),eps
   x=[0.0_dp,1.0_dp,2.0_dp];f=[0.2_dp,0.5_dp,0.3_dp];q=[0.1_dp,0.6_dp,0.95_dp]
   call kw_quantiles(x,f,q,out)
   if(any(abs(out-[0.0_dp,1.0_dp,2.0_dp])>1.0e-12_dp))error stop "quantiles"
   call kw_smooth(x,f,0.5_dp,1,smooth)
   if(any(smooth<0.0_dp))error stop "smooth"
   eps=huber_eps(1.345_dp)
   if(eps<=0.0_dp.or.eps>=1.0_dp)error stop "huber eps"
   if(huber_pdf(0.0_dp,1.0_dp,1.345_dp,eps)<=0.0_dp)error stop "huber pdf"
   print *,"test_utils: PASS"
end program test_utils
