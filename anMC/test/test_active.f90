program test_active
  use anmc
  implicit none
  integer,parameter::n=50,q=10
  real(dp)::e(n,2),mu(n),sig(n,n),pn(n)
  integer,allocatable::idx(:)
  integer::i,m
  type(active_dims_result)::r
  call seed_fortran_rng(2468)
  mu=0.0_dp; sig=0.2_dp
  do i=1,n
    sig(i,i)=1.0_dp
    e(i,1)=real(i-1,dp)/real(n-1,dp)
    e(i,2)=sin(0.2_dp*real(i,dp))
    pn(i)=0.05_dp+0.9_dp*real(i,dp)/real(n,dp)
  end do
  do m=1,5
    idx=select_active_dims(q,e,0.0_dp,mu,sig,pn,m)
    if(size(idx)/=q) error stop 'wrong active-set size'
    if(any(idx<1).or.any(idx>n)) error stop 'active index out of range'
    do i=2,q
      if(idx(i)<=idx(i-1)) error stop 'active indexes not unique/sorted'
    end do
  end do
  r=select_q_dims(e,0.0_dp,mu,sig,pn,method=0,limits=[8,18],reduced_return=.false., &
                  prob_control=genz_bretz(maxpts=50000,batches=12,abseps=1e-4_dp))
  if(.not.r%ok .and. len_trim(r%message)==0) error stop 'select_q_dims failed silently'
  if(size(r%ind_q)<8 .or. size(r%ind_q)>18) error stop 'selected q outside limits'
  if(size(r%eq,1)/=size(r%ind_q)) error stop 'Eq size mismatch'
  if(size(r%k_eq,1)/=size(r%ind_q)) error stop 'KEq size mismatch'
  print *, 'test_active: PASS'
end program test_active
