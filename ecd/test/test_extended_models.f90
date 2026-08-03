! SPDX-License-Identifier: Artistic-2.0
program test_extended_models
  use ecd_api
  implicit none
  type(ecld_model) :: d1,d3,d4
  real(dp) :: a,b,r
  integer :: st

  d1=ecld_new(lambda=1.0_dp,sigma=1.0_dp)
  d3=ecld_new(lambda=3.0_dp,sigma=0.4_dp)
  d4=ecld_new(lambda=4.0_dp,sigma=1.0_dp)

  call check_close('erfq',erfq(5.0_dp,1),0.19621886146307763_dp,2.0e-12_dp)
  call check_close('erfq sum',erfq_sum(5.0_dp,1),0.19621886146111664_dp,2.0e-12_dp)

  a=ecld_ogf_star(d1,1.3_dp)
  b=ecld_ogf_star_analytic(d1,1.3_dp)
  call check_close('lambda1 star ogf',a,b,5.0e-11_dp)
  a=ecld_ogf_star(d4,1.3_dp)
  b=ecld_ogf_star_analytic(d4,1.3_dp)
  call check_close('lambda4 star ogf',a,b,5.0e-11_dp)

  call check_close('incomplete moment order zero', &
    ecld_incomplete_moment(d3,0.5_dp,0,'c'),0.36933574624430943_dp,2.0e-11_dp)

  d4=ecld_new(lambda=4.0_dp,sigma=0.18_dp,mu=0.01_dp)
  a=ecld_ogf_quartic(d4,0.03_dp,'c',.false.)
  b=ecld_ogf_quartic(d4,0.03_dp,'p',.false.)
  call check_close('quartic put-call parity',a-b, &
    ecld_mgf_quartic(d4)-exp(0.03_dp),2.0e-11_dp)
  r=ecld_quartic_q(d4,0.0_dp,'p',st)
  call check_status('quartic Q status',st)
  call check_close('quartic Q',r,4.3740334509868175_dp,2.0e-10_dp)

  r=ecld_quartic_sn0_atm_ki(st)
  call check_status('quartic atm status',st)
  call check_close('quartic atm strike',r,-11.480999814198681_dp,2.0e-8_dp)
  call check_close('quartic rho/sd',ecld_quartic_sn0_rho_stdev(), &
    1.0480670968248731_dp,2.0e-12_dp)
  call check_close('quartic skew',ecld_quartic_sn0_skew(), &
    -0.44952825163743793_dp,2.0e-12_dp)
  r=ecld_quartic_sn0_max_rnv(st)
  call check_status('quartic max status',st)
  call check_close('quartic max ratio',r,0.29619969695404613_dp,2.0e-8_dp)

  print '(a)','test_extended_models: PASS'
contains
  subroutine check_close(name,x,y,tol)
    character(len=*),intent(in)::name
    real(dp),intent(in)::x,y,tol
    if(abs(x-y)>tol*(1.0_dp+abs(y)))then
      write(*,'(a,2es24.15)')trim(name)//' failed: ',x,y
      error stop 1
    end if
  end subroutine check_close
  subroutine check_status(name,s)
    character(len=*),intent(in)::name
    integer,intent(in)::s
    if(s/=ecd_ok)then
      write(*,'(a,i0)')trim(name)//' failed: ',s
      error stop 1
    end if
  end subroutine check_status
end program test_extended_models
