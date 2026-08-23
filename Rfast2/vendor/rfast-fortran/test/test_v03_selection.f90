program test_v03_selection
   use rfast
   use zigg, only : ziggurat_rng
   implicit none
   integer, parameter :: n=240,p=4
   real(dp)::x(n,p),yb(n),yp(n),yq(n),yn(n),yg(n),yw(n),ymv(n,2),u,mu,pr,e
   integer::ym(n),i,j,failures
   type(ziggurat_rng)::rng
   type(selection_result)::s
   failures=0
   call rng%set_seed(987654321)
   do i=1,n
      do j=1,p;x(i,j)=rng%rnorm();end do
      pr=1.0_dp/(1.0_dp+exp(-1.5_dp*x(i,1)))
      yb(i)=merge(1.0_dp,0.0_dp,rng%runi()<pr)
      mu=exp(0.2_dp+0.65_dp*x(i,1));e=0.35_dp*rng%rnorm();yp(i)=real(max(0,nint(mu+e)),dp);yq(i)=yp(i)
      yn(i)=exp(0.25_dp+0.8_dp*x(i,1))+0.06_dp*rng%rnorm()+0.08_dp
      yg(i)=exp(0.15_dp+0.75_dp*x(i,1))*(0.85_dp+0.15_dp*abs(rng%rnorm()))
      u=max(1.0e-8_dp,min(1.0_dp-1.0e-8_dp,rng%runi()))
      yw(i)=exp(0.1_dp+0.7_dp*x(i,1))*(-log(u))**(1.0_dp/1.8_dp)
      ymv(i,1)=1.2_dp*x(i,1)+0.4_dp*rng%rnorm();ymv(i,2)=-0.9_dp*x(i,1)+0.4_dp*rng%rnorm()
      pr=1.0_dp/(1.0_dp+exp(-1.0_dp*x(i,1)));u=rng%runi()
      if(u<0.30_dp*(1.0_dp-pr))then;ym(i)=1
      else if(u<0.30_dp*(1.0_dp-pr)+0.55_dp*pr)then;ym(i)=2
      else;ym(i)=3;end if
   end do
   s=omp_glm(yb,x,OMP_LOGISTIC);call first_one(s,'omp logistic')
   s=omp_glm(yp,x,OMP_POISSON);call first_one(s,'omp poisson')
   s=omp_glm(yq,x,OMP_QUASIPOISSON);call first_one(s,'omp quasipoisson')
   s=omp_glm(yb,x,OMP_QUASIBINOMIAL);call first_one(s,'omp quasibinomial')
   s=omp_glm(yn,x,OMP_NORMLOG);call first_one(s,'omp normlog')
   s=omp_glm(yg,x,OMP_GAMMA);call first_one(s,'omp gamma')
   s=omp_glm(yw,x,OMP_WEIBULL);call first_one(s,'omp weibull')
   s=omp_multivariate(ymv,x);call first_one(s,'omp multivariate')
   s=omp_multinomial(ym,x);call first_one(s,'omp multinomial')
   if(failures==0)then
      print *,'test_v03_selection: PASS'
   else
      print *,'test_v03_selection: FAIL',failures
      error stop 1
   end if
contains
   subroutine first_one(s,name)
      type(selection_result),intent(in)::s
      character(*),intent(in)::name
      if(s%status/=0.or..not.allocated(s%selected).or.size(s%selected)<1)then
         print *,'FAIL: ',trim(name),' status/allocation',s%status;failures=failures+1
      else if(s%selected(1)/=1)then
         print *,'FAIL: ',trim(name),' selected',s%selected(1);failures=failures+1
      end if
   end subroutine first_one
end program test_v03_selection
