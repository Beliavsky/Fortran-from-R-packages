program compare_rugarch
 use rugarch, only: dp, distribution_pdf, distribution_cdf, distribution_quantile, &
  dist_norm,dist_std,dist_ged,dist_snorm,dist_sstd,dist_sged,dist_jsu
 implicit none
 integer,parameter::n=1001,nd=7
 integer,parameter::kinds(nd)=[dist_norm,dist_std,dist_ged,dist_snorm,dist_sstd,dist_sged,dist_jsu]
 character(5),parameter::names(nd)=[character(5)::'norm','std','ged','snorm','sstd','sged','jsu']
 real(dp),parameter::shapes(nd)=[5._dp,7._dp,1.6_dp,5._dp,7._dp,1.6_dp,1.8_dp]
 real(dp),parameter::skews(nd)=[1._dp,1._dp,1._dp,1.3_dp,1.3_dp,1.3_dp,.5_dp]
 real(dp)::x(n),p(n),y(n),t0,t1,value; integer::i,j,k,u,reps; character(512)::out
 call get_command_argument(1,out); if(len_trim(out)==0)out='fortran_results.csv'
 open(newunit=u,file=trim(out),status='replace'); write(u,'(a)')'case,value,seconds,abs_tol,rel_tol'
 do i=1,n; x(i)=-4._dp+8._dp*real(i-1,dp)/real(n-1,dp); p(i)=.001_dp+.998_dp*real(i-1,dp)/real(n-1,dp); end do
 do k=1,nd
  reps=1000
  call cpu_time(t0); do j=1,reps; do i=1,n; y(i)=distribution_pdf(x(i),kinds(k),shapes(k),skews(k)); end do; value=wcheck(y); end do
  call cpu_time(t1); call emit(trim(names(k))//'_density',value,t1-t0)
  call cpu_time(t0); do j=1,reps; do i=1,n; y(i)=distribution_cdf(x(i),kinds(k),shapes(k),skews(k)); end do; value=wcheck(y); end do
  call cpu_time(t1); call emit(trim(names(k))//'_cdf',value,t1-t0)
  reps=50; if(k==1 .or. k==4)reps=2000
  call cpu_time(t0); do j=1,reps; do i=1,n; y(i)=distribution_quantile(p(i),kinds(k),shapes(k),skews(k)); end do; value=wcheck(y); end do
  call cpu_time(t1); call emit(trim(names(k))//'_quantile',value,t1-t0)
 end do
 close(u)
contains
 real(dp) function wcheck(q); real(dp),intent(in)::q(:); integer::ii; wcheck=sum(q*[(real(ii,dp),ii=1,size(q))]); end function
 subroutine emit(name,v,s); character(*),intent(in)::name; real(dp),intent(in)::v,s
  write(u,'(a,",",es24.16,",",es16.8,",",es12.4,",",es12.4)')trim(name),v,s,5e-5_dp,2e-7_dp
 end subroutine
end program
