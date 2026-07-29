! SPDX-License-Identifier: GPL-3.0-or-later
program test_distributions
  use acdm
  implicit none
  integer :: failures
  failures=0
  call check_dist(DIST_EXPONENTIAL,[real(dp)::],0.4493289641172216_dp,0.5506710358827784_dp,1.203972804325936_dp)
  call check_dist(DIST_WEIBULL,[1.5_dp],0.6229327842207739_dp,0.4586762441723429_dp,1.253658410942133_dp)
  call check_dist(DIST_BURR,[2.5_dp,0.3_dp],0.8769351827830919_dp,0.3909687031173490_dp,1.192389832701695_dp)
  call check_dist(DIST_GENGAMMA,[2.0_dp,1.2_dp],0.7456299360340816_dp,0.4325386721068839_dp,1.219091834791879_dp)
  call check_dist(DIST_GENF,[2.0_dp,3.0_dp,1.5_dp],0.7734747661348950_dp,0.4856912248642452_dp,1.141748375235077_dp)
  call check_dist(DIST_QWEIBULL,[1.5_dp,1.2_dp],0.5992047298407931_dp,0.5069420193881261_dp,1.190873049857338_dp)
  call check_dist(DIST_MIXQWE,[0.8_dp,1.5_dp,1.2_dp,0.5_dp],0.5396187604361979_dp,0.5209050661727738_dp,1.200624141566930_dp)
  call check_dist(DIST_MIXQWW,[0.8_dp,1.5_dp,1.2_dp,1.0_dp,1.5_dp],0.6070810815889165_dp,0.4985884250726392_dp,1.199690954117487_dp)
  call check_dist(DIST_BIRNBAUM_SAUNDERS,[0.8_dp],0.6221335915188049_dp,0.5271544108875760_dp,1.148976561904265_dp)
  call check_inversion(DIST_MIXINVGAUSS,[0.7_dp,2.0_dp,0.5_dp])
  call check_component_identities
  if(failures>0) error stop 'test_distributions failed'
  print '(a)','test_distributions: PASS'
contains
  subroutine check_dist(code,para,pdf_ref,cdf_ref,q_ref)
    integer,intent(in)::code
    real(dp),intent(in)::para(:),pdf_ref,cdf_ref,q_ref
    real(dp)::f,c,q
    f=distribution_pdf(0.8_dp,code,para,.true.)
    c=distribution_cdf(0.8_dp,code,para,.true.)
    q=distribution_quantile(0.7_dp,code,para,.true.)
    call assert_close(f,pdf_ref,2.0e-9_dp,'pdf')
    call assert_close(c,cdf_ref,2.0e-9_dp,'cdf')
    call assert_close(q,q_ref,5.0e-8_dp,'quantile')
    call assert_close(distribution_cdf(q,code,para,.true.),0.7_dp,2.0e-8_dp,'inverse')
  end subroutine
  subroutine check_inversion(code,para)
    integer,intent(in)::code
    real(dp),intent(in)::para(:)
    real(dp)::p,q
    integer::i
    do i=1,9
      p=real(i,dp)/10.0_dp
      q=distribution_quantile(p,code,para,.true.)
      call assert_close(distribution_cdf(q,code,para,.true.),p,2.0e-7_dp,'mixture inverse')
      if(distribution_pdf(q,code,para,.true.)<=0.0_dp) then
        failures=failures+1;print *,'nonpositive pdf',code,p,q
      end if
    end do
  end subroutine
  subroutine check_component_identities
    real(dp)::x,p
    x=qburr(0.37_dp,0.8_dp,2.5_dp,0.3_dp);p=pburr(x,0.8_dp,2.5_dp,0.3_dp)
    call assert_close(p,0.37_dp,2.0e-12_dp,'burr direct inverse')
    x=qgengamma(0.63_dp,1.2_dp,2.0_dp,0.7_dp)
    call assert_close(pgengamma(x,1.2_dp,2.0_dp,0.7_dp),0.63_dp,2.0e-10_dp,'gengamma inverse')
    x=qgenf(0.42_dp,2.0_dp,3.0_dp,1.5_dp,0.8_dp)
    call assert_close(pgenf(x,2.0_dp,3.0_dp,1.5_dp,0.8_dp),0.42_dp,2.0e-9_dp,'genf inverse')
    x=qqweibull(0.58_dp,1.4_dp,0.8_dp,1.1_dp)
    call assert_close(pqweibull(x,1.4_dp,0.8_dp,1.1_dp),0.58_dp,2.0e-12_dp,'qweib q<1')
  end subroutine
  subroutine assert_close(x,y,tol,label)
    real(dp),intent(in)::x,y,tol
    character(*),intent(in)::label
    if(abs(x-y)>tol*max(1.0_dp,abs(y))) then
      failures=failures+1
      print '(a,2es24.14)','FAIL '//trim(label)//': ',x,y
    end if
  end subroutine
end program
