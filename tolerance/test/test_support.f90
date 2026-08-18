program test_support
  use tolerance
  implicit none
  real(dp)::kt(1,1,1),p,a
  real(dp)::x(8)=[-1.2_dp,-0.7_dp,-0.2_dp,0.1_dp,0.4_dp,0.8_dp,1.1_dp,1.5_dp]
  type(tolerance_interval)::bt
  integer::fail
  fail=0
  call k_table([20],[0.05_dp],[0.95_dp],1,'HE',50,kt)
  if(abs(kt(1,1,1)-k_factor(20,0.05_dp,0.95_dp,1,'HE'))>1e-12_dp)call bad('k table')
  p=norm_oc_content(kt(1,1,1),20,0.05_dp,1,'HE');if(abs(p-0.95_dp)>2e-6_dp)call bad('OC P')
  a=norm_oc_alpha(kt(1,1,1),20,0.95_dp,1,'HE');if(abs(a-0.05_dp)>2e-6_dp)call bad('OC alpha')
  if(abs(k_factor(20,0.05_dp,0.95_dp,1,'HE')-2.3960016837521687_dp)>2e-8_dp)call bad('K side1')
  bt=bonf_tol_int(norm_cb,0.01_dp,0.02_dp,0.05_dp)
  if(abs(bt%p-0.97_dp)>1e-12_dp .or. bt%upper<bt%lower)call bad('bonf callback')
  if(fail==0)then;print '(a)','test_support: PASS';else;error stop 1;end if
contains
  function norm_cb(pp,aa) result(out)
    real(dp),intent(in)::pp,aa
    type(tolerance_interval)::out
    out=normtol_int(x,aa,pp,1,'HE')
  end function norm_cb
  subroutine bad(nm);character(len=*),intent(in)::nm;print *,trim(nm);fail=fail+1;end subroutine
end program
