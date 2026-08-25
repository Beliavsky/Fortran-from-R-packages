program test_rotation
  use gpa_kinds, only: dp
  use gpa_linalg, only: eye, random_orthogonal
  use gpa_criteria, only: criterion_result, criterion_options, evaluate_criterion
  use gpa_rotation
  use gpa_transforms
  use gpa_diagnostics
  implicit none
  real(dp) :: a(9,3), rec(9,3), q(3,3), w(9), rmat(9,9), phi(3,3)
  real(dp) :: fsi(3),mfsi,minfsi,auc(3),aadj(3),mauc,minauc,pct
  integer :: hp(3),hpt,hpm,info,fails
  integer,parameter :: idx(3)=[1,4,7]
  type(rotation_result)::ro,rq,rr,re,rlp
  type(rotation_options)::opt
  type(criterion_result)::c0
  type(criterion_options)::co
  type(fit_stats)::fs
  type(simplicity_stats)::ss
  fails=0
  a(1,:)=[0.82_dp,0.21_dp,0.10_dp]
  a(2,:)=[0.77_dp,0.28_dp,0.12_dp]
  a(3,:)=[0.71_dp,0.18_dp,0.22_dp]
  a(4,:)=[0.15_dp,0.80_dp,0.18_dp]
  a(5,:)=[0.21_dp,0.74_dp,0.20_dp]
  a(6,:)=[0.09_dp,0.70_dp,0.27_dp]
  a(7,:)=[0.12_dp,0.15_dp,0.81_dp]
  a(8,:)=[0.25_dp,0.11_dp,0.73_dp]
  a(9,:)=[0.18_dp,0.24_dp,0.69_dp]
  opt%maxit=1000; opt%eps=1.0e-7_dp; opt%algorithm='bb'; opt%fwindow=10
  call evaluate_criterion('varimax',a,c0,co)
  call gpforth(a,'varimax',ro,co,opt)
  if(.not.ro%converged) then; print *,'varimax did not converge',ro%table(size(ro%table,1),3); fails=fails+1; end if
  if(ro%objective>c0%f+1.0e-9_dp) then; print *,'varimax objective did not improve'; fails=fails+1; end if
  if(maxval(abs(matmul(transpose(ro%th),ro%th)-eye(3)))>2.0e-6_dp) then
    print *,'orthogonality failure'; fails=fails+1
  end if
  rec=matmul(a,ro%th)
  if(maxval(abs(rec-ro%loadings))>1.0e-8_dp) then; print *,'orth reconstruction failure'; fails=fails+1; end if
  call evaluate_criterion('quartimin',a,c0,co)
  call gpfoblq(a,'quartimin',rq,co,opt)
  if(rq%objective>c0%f+1.0e-8_dp) then; print *,'quartimin objective did not improve'; fails=fails+1; end if
  rec=matmul(rq%loadings,transpose(rq%th))
  if(maxval(abs(rec-a))>2.0e-6_dp) then; print *,'oblique reconstruction failure',maxval(abs(rec-a)); fails=fails+1; end if
  if(maxval(abs([(rq%phi(info,info)-1.0_dp,info=1,3)]))>2.0e-6_dp) then
    print *,'oblique phi diagonal failure'; fails=fails+1
  end if
  call random_orthogonal(3,q,info)
  if(info/=0 .or. maxval(abs(matmul(transpose(q),q)-eye(3)))>2.0e-10_dp) then
    print *,'random orthogonal failure'; fails=fails+1
  end if
  call normalizing_weights(a,1,w)
  if(minval(w)<=0.0_dp) then; print *,'normalization failure'; fails=fails+1; end if
  call eiv_rotate(a,idx,rr)
  if(rr%info/=0) then; print *,'eiv failed'; fails=fails+1
  else
    if(maxval(abs(rr%loadings(idx,:)-eye(3)))>1.0e-8_dp) then; print *,'eiv identity failure'; fails=fails+1; end if
  end if
  call echelon_rotate(a,idx,re)
  if(re%info/=0) then; print *,'echelon failed'; fails=fails+1
  else
    if(abs(re%loadings(idx(1),2))+abs(re%loadings(idx(1),3))+abs(re%loadings(idx(2),3))>1.0e-8_dp) then
      print *,'echelon triangular failure'; fails=fails+1
    end if
  end if
  opt%maxit=100; opt%eps=1.0e-5_dp
  call lp_rotate(a,1.0_dp,.true.,rlp,opt)
  if(.not.(rlp%objective<huge(1.0_dp))) then; print *,'lp rotation failure'; fails=fails+1; end if
  phi=eye(3); rmat=matmul(matmul(ro%loadings,phi),transpose(ro%loadings))
  do info=1,9; rmat(info,info)=1.0_dp; end do
  call calc_fitstats(rmat,ro%loadings,phi,500,fs)
  if(fs%srmr>1.0e-10_dp) then; print *,'fit stats SRMR failure',fs%srmr; fails=fails+1; end if
  call calc_fsi(ro%loadings,fsi,mfsi,minfsi)
  call calc_auc(ro%loadings,auc,aadj,mauc,minauc)
  ss=calc_simplicity(ro%loadings)
  call calc_hyperplane(ro%loadings,0.15_dp,hp,hpt,hpm,pct)
  if(minval(fsi)<-1.0e-10_dp .or. minval(auc)<0.0_dp .or. hpt<0 .or. ss%gini<0.0_dp) then
    print *,'simplicity diagnostics failure'; fails=fails+1
  end if
  if(fails/=0) error stop 1
  print *, 'test_rotation: PASS'
end program test_rotation
