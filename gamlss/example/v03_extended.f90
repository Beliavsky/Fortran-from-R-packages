program v03_extended
   use gamlss
   implicit none
   integer,parameter::ng=8,m=6,n=ng*m,nt=80
   real(dp)::y(n),xf(n,2),zr(n,2),xv,truth(ng,2)
   real(dp)::lo(nt),up(nt),entry(nt),xi(nt,1),xs(nt,1),u
   integer::grp(n),cens(nt),i,j,g
   type(random_effects_result_t)::re
   type(gamlss_fit_result_t)::truncfit

   i=0
   do g=1,ng
      truth(g,1)=0.45_dp*real(g-(ng+1)/2,dp)/real(ng,dp)
      truth(g,2)=0.30_dp*truth(g,1)+0.08_dp*sin(real(g,dp))
      do j=1,m
         i=i+1;grp(i)=100+g;xv=-1.0_dp+2.0_dp*real(j-1,dp)/real(m-1,dp)
         xf(i,:)=[1.0_dp,xv];zr(i,:)=[1.0_dp,xv]
         y(i)=1.0_dp+0.75_dp*xv+truth(g,1)+truth(g,2)*xv+0.04_dp*sin(real(i,dp))
      end do
   end do
   call fit_gamlss_random_effects(y,xf,zr,grp,GAMLSS_NO,re,max_outer=10)
   print '(a,2f10.4)', 'Random-slope fixed coefficients: ',re%model%mu%coefficients(1:2)
   print '(a,3f10.5)', 'Random covariance (v11,v12,v22): ', &
      re%covariance(1,1),re%covariance(1,2),re%covariance(2,2)

   xi=1.0_dp;xs=1.0_dp;entry=0.0_dp
   do i=1,nt
      u=(real(i,dp)-0.5_dp)/real(nt,dp)
      lo(i)=qNO(pNO(0.0_dp,0.4_dp,1.1_dp)+(1.0_dp-pNO(0.0_dp,0.4_dp,1.1_dp))*u,0.4_dp,1.1_dp)
      up(i)=lo(i);cens(i)=CENS_EXACT
   end do
   call fit_gamlss_censored(lo,up,cens,xi,GAMLSS_NO,truncfit,x_sigma=xs,entry=entry)
   print '(a,2f10.4)', 'Delayed-entry normal (mu,sigma): ',truncfit%fitted_mu(1),truncfit%fitted_sigma(1)
end program v03_extended
