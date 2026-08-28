module gpa_api
  use gpa_kinds, only: dp
  use gpa_linalg, only: random_orthogonal
  use gpa_criteria, only: criterion_options
  use gpa_rotation, only: rotation_result, rotation_options, gpforth, gpfoblq
  use gpa_rotation, only: rotate_random_starts, lp_rotate
  use gpa_transforms, only: eiv_rotate, echelon_rotate
  implicit none
  private
  public :: rotation_result, rotation_options, criterion_options
  public :: gpfrsorth, gpfrsoblq, random_start
  public :: oblimin, quartimin, target_t, target_q, pst_t, pst_q
  public :: oblimax, binormamin, entropy, quartimax, varimax, simplimax
  public :: bentler_t, bentler_q, tandem_i, tandem_ii, geomin_t, geomin_q
  public :: bigeomin_t, bigeomin_q, cf_t, cf_q, infomax_t, infomax_q
  public :: mccammon, bifactor_t, bifactor_q, equamax, parsimax, varimin
  public :: lp_t, lp_q, eiv, echelon
contains

  subroutine gpfrsorth(a,method,nstarts,r,crit,opts)
    real(dp),intent(in)::a(:,:)
    character(len=*),intent(in)::method
    integer,intent(in)::nstarts
    type(rotation_result),intent(out)::r
    type(criterion_options),intent(in),optional::crit
    type(rotation_options),intent(in),optional::opts
    call rotate_random_starts(a,method,.true.,nstarts,r,crit,opts)
  end subroutine gpfrsorth

  subroutine gpfrsoblq(a,method,nstarts,r,crit,opts)
    real(dp),intent(in)::a(:,:)
    character(len=*),intent(in)::method
    integer,intent(in)::nstarts
    type(rotation_result),intent(out)::r
    type(criterion_options),intent(in),optional::crit
    type(rotation_options),intent(in),optional::opts
    call rotate_random_starts(a,method,.false.,nstarts,r,crit,opts)
  end subroutine gpfrsoblq

  subroutine random_start(k,t,info)
    integer,intent(in)::k
    real(dp),intent(out)::t(k,k)
    integer,intent(out)::info
    call random_orthogonal(k,t,info)
  end subroutine random_start

  subroutine oblimin(a,r,gam,opts)
    real(dp),intent(in)::a(:,:)
    type(rotation_result),intent(out)::r
    real(dp),intent(in),optional::gam
    type(rotation_options),intent(in),optional::opts
    type(criterion_options)::c
    if(present(gam)) c%gam=gam
    call gpfoblq(a,'oblimin',r,c,opts)
  end subroutine oblimin

  subroutine quartimin(a,r,opts)
    real(dp),intent(in)::a(:,:)
    type(rotation_result),intent(out)::r
    type(rotation_options),intent(in),optional::opts
    type(criterion_options)::c
    call gpfoblq(a,'quartimin',r,c,opts)
  end subroutine quartimin

  subroutine target_t(a,target,r,opts)
    real(dp),intent(in)::a(:,:),target(:,:)
    type(rotation_result),intent(out)::r
    type(rotation_options),intent(in),optional::opts
    type(criterion_options)::c
    c%target=target
    call gpforth(a,'target',r,c,opts)
  end subroutine target_t

  subroutine target_q(a,target,r,opts)
    real(dp),intent(in)::a(:,:),target(:,:)
    type(rotation_result),intent(out)::r
    type(rotation_options),intent(in),optional::opts
    type(criterion_options)::c
    c%target=target
    call gpfoblq(a,'target',r,c,opts)
  end subroutine target_q

  subroutine pst_t(a,w,target,r,opts)
    real(dp),intent(in)::a(:,:),w(:,:),target(:,:)
    type(rotation_result),intent(out)::r
    type(rotation_options),intent(in),optional::opts
    type(criterion_options)::c
    c%weight=w
    c%target=target
    call gpforth(a,'pst',r,c,opts)
  end subroutine pst_t

  subroutine pst_q(a,w,target,r,opts)
    real(dp),intent(in)::a(:,:),w(:,:),target(:,:)
    type(rotation_result),intent(out)::r
    type(rotation_options),intent(in),optional::opts
    type(criterion_options)::c
    c%weight=w
    c%target=target
    call gpfoblq(a,'pst',r,c,opts)
  end subroutine pst_q

  subroutine oblimax(a,r,opts)
    real(dp),intent(in)::a(:,:)
    type(rotation_result),intent(out)::r
    type(rotation_options),intent(in),optional::opts
    type(criterion_options)::c
    call gpfoblq(a,'oblimax',r,c,opts)
  end subroutine oblimax

  subroutine binormamin(a,r,opts)
    real(dp),intent(in)::a(:,:)
    type(rotation_result),intent(out)::r
    type(rotation_options),intent(in),optional::opts
    type(criterion_options)::c
    call gpfoblq(a,'binormamin',r,c,opts)
  end subroutine binormamin

  subroutine entropy(a,r,opts)
    real(dp),intent(in)::a(:,:)
    type(rotation_result),intent(out)::r
    type(rotation_options),intent(in),optional::opts
    type(criterion_options)::c
    call gpforth(a,'entropy',r,c,opts)
  end subroutine entropy

  subroutine quartimax(a,r,opts)
    real(dp),intent(in)::a(:,:)
    type(rotation_result),intent(out)::r
    type(rotation_options),intent(in),optional::opts
    type(criterion_options)::c
    call gpforth(a,'quartimax',r,c,opts)
  end subroutine quartimax

  subroutine varimax(a,r,opts)
    real(dp),intent(in)::a(:,:)
    type(rotation_result),intent(out)::r
    type(rotation_options),intent(in),optional::opts
    type(criterion_options)::c
    call gpforth(a,'varimax',r,c,opts)
  end subroutine varimax

  subroutine simplimax(a,r,kzero,opts)
    real(dp),intent(in)::a(:,:)
    type(rotation_result),intent(out)::r
    integer,intent(in),optional::kzero
    type(rotation_options),intent(in),optional::opts
    type(criterion_options)::c
    if(present(kzero)) c%simplimax_k=kzero
    call gpfoblq(a,'simplimax',r,c,opts)
  end subroutine simplimax

  subroutine bentler_t(a,r,opts)
    real(dp),intent(in)::a(:,:)
    type(rotation_result),intent(out)::r
    type(rotation_options),intent(in),optional::opts
    type(criterion_options)::c
    call gpforth(a,'bentler',r,c,opts)
  end subroutine bentler_t

  subroutine bentler_q(a,r,opts)
    real(dp),intent(in)::a(:,:)
    type(rotation_result),intent(out)::r
    type(rotation_options),intent(in),optional::opts
    type(criterion_options)::c
    call gpfoblq(a,'bentler',r,c,opts)
  end subroutine bentler_q

  subroutine tandem_i(a,r,opts)
    real(dp),intent(in)::a(:,:)
    type(rotation_result),intent(out)::r
    type(rotation_options),intent(in),optional::opts
    type(criterion_options)::c
    call gpforth(a,'tandemI',r,c,opts)
  end subroutine tandem_i

  subroutine tandem_ii(a,r,opts)
    real(dp),intent(in)::a(:,:)
    type(rotation_result),intent(out)::r
    type(rotation_options),intent(in),optional::opts
    type(criterion_options)::c
    call gpforth(a,'tandemII',r,c,opts)
  end subroutine tandem_ii

  subroutine geomin_t(a,r,delta,opts)
    real(dp),intent(in)::a(:,:)
    type(rotation_result),intent(out)::r
    real(dp),intent(in),optional::delta
    type(rotation_options),intent(in),optional::opts
    type(criterion_options)::c
    if(present(delta)) c%delta=delta
    call gpforth(a,'geomin',r,c,opts)
  end subroutine geomin_t

  subroutine geomin_q(a,r,delta,opts)
    real(dp),intent(in)::a(:,:)
    type(rotation_result),intent(out)::r
    real(dp),intent(in),optional::delta
    type(rotation_options),intent(in),optional::opts
    type(criterion_options)::c
    if(present(delta)) c%delta=delta
    call gpfoblq(a,'geomin',r,c,opts)
  end subroutine geomin_q

  subroutine bigeomin_t(a,r,delta,opts)
    real(dp),intent(in)::a(:,:)
    type(rotation_result),intent(out)::r
    real(dp),intent(in),optional::delta
    type(rotation_options),intent(in),optional::opts
    type(criterion_options)::c
    if(present(delta)) c%delta=delta
    call gpforth(a,'bigeomin',r,c,opts)
  end subroutine bigeomin_t

  subroutine bigeomin_q(a,r,delta,opts)
    real(dp),intent(in)::a(:,:)
    type(rotation_result),intent(out)::r
    real(dp),intent(in),optional::delta
    type(rotation_options),intent(in),optional::opts
    type(criterion_options)::c
    if(present(delta)) c%delta=delta
    call gpfoblq(a,'bigeomin',r,c,opts)
  end subroutine bigeomin_q

  subroutine cf_t(a,r,kappa,opts)
    real(dp),intent(in)::a(:,:)
    type(rotation_result),intent(out)::r
    real(dp),intent(in),optional::kappa
    type(rotation_options),intent(in),optional::opts
    type(criterion_options)::c
    if(present(kappa)) c%kappa=kappa
    call gpforth(a,'cf',r,c,opts)
  end subroutine cf_t

  subroutine cf_q(a,r,kappa,opts)
    real(dp),intent(in)::a(:,:)
    type(rotation_result),intent(out)::r
    real(dp),intent(in),optional::kappa
    type(rotation_options),intent(in),optional::opts
    type(criterion_options)::c
    if(present(kappa)) c%kappa=kappa
    call gpfoblq(a,'cf',r,c,opts)
  end subroutine cf_q

  subroutine infomax_t(a,r,opts)
    real(dp),intent(in)::a(:,:)
    type(rotation_result),intent(out)::r
    type(rotation_options),intent(in),optional::opts
    type(criterion_options)::c
    call gpforth(a,'infomax',r,c,opts)
  end subroutine infomax_t

  subroutine infomax_q(a,r,opts)
    real(dp),intent(in)::a(:,:)
    type(rotation_result),intent(out)::r
    type(rotation_options),intent(in),optional::opts
    type(criterion_options)::c
    call gpfoblq(a,'infomax',r,c,opts)
  end subroutine infomax_q

  subroutine mccammon(a,r,opts)
    real(dp),intent(in)::a(:,:)
    type(rotation_result),intent(out)::r
    type(rotation_options),intent(in),optional::opts
    type(criterion_options)::c
    call gpforth(a,'mccammon',r,c,opts)
  end subroutine mccammon

  subroutine bifactor_t(a,r,opts)
    real(dp),intent(in)::a(:,:)
    type(rotation_result),intent(out)::r
    type(rotation_options),intent(in),optional::opts
    type(criterion_options)::c
    call gpforth(a,'bifactor',r,c,opts)
  end subroutine bifactor_t

  subroutine bifactor_q(a,r,opts)
    real(dp),intent(in)::a(:,:)
    type(rotation_result),intent(out)::r
    type(rotation_options),intent(in),optional::opts
    type(criterion_options)::c
    call gpfoblq(a,'bifactor',r,c,opts)
  end subroutine bifactor_q

  subroutine equamax(a,r,opts)
    real(dp),intent(in)::a(:,:)
    type(rotation_result),intent(out)::r
    type(rotation_options),intent(in),optional::opts
    type(criterion_options)::c
    c%kappa=real(size(a,2),dp)/(2.0_dp*real(size(a,1),dp))
    call gpforth(a,'cf',r,c,opts)
  end subroutine equamax

  subroutine parsimax(a,r,opts)
    real(dp),intent(in)::a(:,:)
    type(rotation_result),intent(out)::r
    type(rotation_options),intent(in),optional::opts
    type(criterion_options)::c
    c%kappa=real(size(a,2)-1,dp)/real(size(a,2)+size(a,1)-2,dp)
    call gpforth(a,'cf',r,c,opts)
  end subroutine parsimax

  subroutine varimin(a,r,opts)
    real(dp),intent(in)::a(:,:)
    type(rotation_result),intent(out)::r
    type(rotation_options),intent(in),optional::opts
    type(criterion_options)::c
    call gpforth(a,'varimin',r,c,opts)
  end subroutine varimin

  subroutine lp_t(a,p,r,opts)
    real(dp),intent(in)::a(:,:),p
    type(rotation_result),intent(out)::r
    type(rotation_options),intent(in),optional::opts
    call lp_rotate(a,p,.true.,r,opts)
  end subroutine lp_t

  subroutine lp_q(a,p,r,opts)
    real(dp),intent(in)::a(:,:),p
    type(rotation_result),intent(out)::r
    type(rotation_options),intent(in),optional::opts
    call lp_rotate(a,p,.false.,r,opts)
  end subroutine lp_q

  subroutine eiv(a,identity_rows,r)
    real(dp),intent(in)::a(:,:)
    integer,intent(in)::identity_rows(:)
    type(rotation_result),intent(out)::r
    call eiv_rotate(a,identity_rows,r)
  end subroutine eiv

  subroutine echelon(a,reference_rows,r)
    real(dp),intent(in)::a(:,:)
    integer,intent(in)::reference_rows(:)
    type(rotation_result),intent(out)::r
    call echelon_rotate(a,reference_rows,r)
  end subroutine echelon

end module gpa_api
