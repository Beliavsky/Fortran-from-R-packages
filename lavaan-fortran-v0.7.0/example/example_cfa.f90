program example_cfa
   use lavaan
   implicit none
   type(ram_model) :: model
   type(ram_free_map) :: map
   type(sem_fit_result) :: fit
   real(dp) :: lambda(4,1), beta(1,1), psi(1,1), theta(4,4), cov(4,4), mu(4)

   lambda(:,1)=[1.0_dp,0.6_dp,0.6_dp,0.6_dp]
   beta=0.0_dp
   psi=1.0_dp
   theta=0.0_dp
   theta(1,1)=0.8_dp
   theta(2,2)=0.8_dp
   theta(3,3)=0.8_dp
   theta(4,4)=0.8_dp
   model=ram_from_lisrel(lambda,beta,psi,theta)
   allocate(map%matrix_id(8),map%row(8),map%col(8))
   map%matrix_id=[ram_a,ram_a,ram_a,ram_s,ram_s,ram_s,ram_s,ram_s]
   map%row=[2,3,4,1,2,3,4,5]
   map%col=[5,5,5,1,2,3,4,5]

   cov=reshape([ &
      1.7_dp,0.96_dp,1.08_dp,0.84_dp, &
      0.96_dp,1.368_dp,0.864_dp,0.672_dp, &
      1.08_dp,0.864_dp,1.372_dp,0.756_dp, &
      0.84_dp,0.672_dp,0.756_dp,1.288_dp], [4,4])
   mu=0.0_dp
   call fit_ram_cov(model,map,cov,mu,500,fit,'ML')
   write(*,'(a,3f10.4)') 'free loadings: ',fit%par(1:3)
   write(*,'(a,f10.4)') 'chi-square: ',fit%chisq
   write(*,'(a,f10.4)') 'rmsea:      ',fit%rmsea
end program example_cfa
