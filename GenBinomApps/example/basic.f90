program basic
   use genbinomapps
   implicit none
   integer(i64) :: sizes(3),x(3),nreq
   real(dp) :: probs(3),d(3)
   type(confidence_interval) :: ci

   sizes = [100_i64,100_i64,200_i64]
   probs = [0.001_dp,0.005_dp,0.01_dp]
   x = [0_i64,1_i64,2_i64]
   call dgbinom_vec(x,sizes,probs,d)

   print '(a,3f12.8)', "P(X=0:2): ", d
   print '(a,f12.8)', "P(X>2):   ", pgbinom(2.0_dp,sizes,probs,.false.)

   ci = clopper_pearson_ci(5_i64,100000_i64,0.05_dp,"upper")
   print '(a,es14.6)', "Upper CP limit: ", ci%upper

   nreq = n_clopper_pearson(8_i64,0.0002_dp)
   print '(a,i0)', "Required n: ", nreq
end program basic
