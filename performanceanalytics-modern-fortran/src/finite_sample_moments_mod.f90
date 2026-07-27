! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of PerformanceAnalytics.
! Copyright (C) original PerformanceAnalytics authors and translation contributors.
! Distributed under GNU GPL version 2 or, at your option, any later version.
module finite_sample_moments_mod
  use kinds_mod, only: dp
  use statistics_mod, only: mean_value
  use comoments_mod, only: covariance_matrix, coskewness_unique, cokurtosis_unique, &
    m3_vec_to_mat, m4_vec_to_mat, m3_inner_product, m4_inner_product
  implicit none
  private

  type, public :: exact_shrinkage_result
    integer :: order = 0
    integer :: n_targets = 0
    logical :: converged = .false.
    logical :: unbiased_mse = .false.
    real(dp), allocatable :: estimate(:,:)
    real(dp), allocatable :: sample(:,:)
    real(dp), allocatable :: target_vectors(:,:)
    real(dp), allocatable :: lambda(:)
    real(dp), allocatable :: a(:,:)
    real(dp), allocatable :: b(:)
    real(dp) :: objective = huge(1.0_dp)
  end type exact_shrinkage_result

  public :: exact_m2_shrinkage, exact_m3_shrinkage, exact_m4_shrinkage
  public :: exact_vm2_terms, exact_vm3_terms, exact_vm3_kstat_terms, exact_vm4_terms
  public :: solve_shrinkage_qp

contains

  subroutine center_matrix(r, x)
    real(dp), intent(in) :: r(:,:)
    real(dp), allocatable, intent(out) :: x(:,:)
    integer :: j
    allocate(x(size(r,1),size(r,2)))
    do j=1,size(r,2)
      x(:,j)=r(:,j)-mean_value(r(:,j))
    end do
  end subroutine center_matrix

  subroutine pair_moment(x, power_left, power_right, out, sum_scale)
    real(dp), intent(in) :: x(:,:)
    integer, intent(in) :: power_left, power_right
    real(dp), intent(out) :: out(:,:)
    logical, intent(in), optional :: sum_scale
    integer :: i,j
    logical :: raw_sum
    raw_sum=.false.; if(present(sum_scale)) raw_sum=sum_scale
    do i=1,size(x,2)
      do j=1,size(x,2)
        out(i,j)=sum(x(:,i)**power_left*x(:,j)**power_right)
      end do
    end do
    if(.not.raw_sum) out=out/real(max(size(x,1),1),dp)
  end subroutine pair_moment

  pure real(dp) function safe_sqrt_ratio(a,b) result(v)
    real(dp),intent(in)::a,b
    if(a>=0.0_dp .and. b>tiny(1.0_dp))then;v=sqrt(a/b);else;v=0.0_dp;end if
  end function safe_sqrt_ratio

  pure real(dp) function safe_fourth_root_ratio(a,b) result(v)
    real(dp),intent(in)::a,b
    if(a>=0.0_dp .and. b>tiny(1.0_dp))then;v=sqrt(sqrt(a/b));else;v=0.0_dp;end if
  end function safe_fourth_root_ratio

  pure integer function n_unique_m3(p) result(n)
    integer,intent(in)::p
    n=p*(p+1)*(p+2)/6
  end function n_unique_m3

  pure integer function n_unique_m4(p) result(n)
    integer,intent(in)::p
    n=p*(p+1)*(p+2)*(p+3)/24
  end function n_unique_m4


  pure integer function m4_multiplicity(i,j,k,l) result(mult)
    integer,intent(in)::i,j,k,l
    if(i==l)then
      mult=1
    else if((i==k .and. k<l) .or. (i<j .and. j==l))then
      mult=4
    else if(i==j .and. k==l .and. j<k)then
      mult=6
    else if(i==j .or. j==k .or. k==l)then
      mult=12
    else
      mult=24
    end if
  end function m4_multiplicity

  subroutine project_simplex_leq_one(x)
    real(dp),intent(inout)::x(:)
    real(dp)::lo,hi,mid,s
    integer::it
    x=max(x,0.0_dp)
    if(sum(x)<=1.0_dp)return
    lo=minval(x)-1.0_dp;hi=maxval(x)
    do it=1,100
      mid=0.5_dp*(lo+hi);s=sum(max(x-mid,0.0_dp))
      if(s>1.0_dp)then;lo=mid;else;hi=mid;end if
    end do
    x=max(x-hi,0.0_dp)
  end subroutine project_simplex_leq_one

  pure real(dp) function qp_objective(a,b,x) result(value)
    real(dp),intent(in)::a(:,:),b(:),x(:)
    value=0.5_dp*dot_product(x,matmul(a,x))-dot_product(b,x)
  end function qp_objective

  subroutine solve_shrinkage_qp(a,b,lambda,converged,objective)
    real(dp),intent(in)::a(:,:),b(:)
    real(dp),intent(out)::lambda(:)
    logical,intent(out),optional::converged
    real(dp),intent(out),optional::objective
    real(dp),allocatable::grad(:),old(:),trial(:)
    real(dp)::lipschitz,step,obj,trial_obj
    integer::m,it,ls
    logical::ok
    m=size(b);lambda=0.0_dp;ok=.false.
    if(m==0)then
      if(present(converged))converged=.true.
      if(present(objective))objective=0.0_dp
      return
    end if
    if(m==1)then
      if(a(1,1)>tiny(1.0_dp))lambda(1)=max(0.0_dp,min(1.0_dp,b(1)/a(1,1)))
      ok=.true.
      if(present(converged))converged=ok
      if(present(objective))objective=qp_objective(a,b,lambda)
      return
    end if
    allocate(grad(m),old(m),trial(m))
    lipschitz=max(maxval(sum(abs(a),dim=2)),1.0e-12_dp)
    step=1.0_dp/lipschitz
    obj=qp_objective(a,b,lambda)
    do it=1,20000
      old=lambda;grad=matmul(a,lambda)-b
      trial=lambda-step*grad;call project_simplex_leq_one(trial)
      trial_obj=qp_objective(a,b,trial)
      do ls=1,20
        if(trial_obj<=obj+1.0e-14_dp)exit
        step=0.5_dp*step
        trial=lambda-step*grad;call project_simplex_leq_one(trial)
        trial_obj=qp_objective(a,b,trial)
      end do
      lambda=trial;obj=trial_obj
      if(maxval(abs(lambda-old))<=1.0e-11_dp*(1.0_dp+maxval(abs(old))))then
        ok=.true.;exit
      end if
      if(mod(it,50)==0)step=min(1.0_dp/lipschitz,1.1_dp*step)
    end do
    if(present(converged))converged=ok
    if(present(objective))objective=obj
  end subroutine solve_shrinkage_qp

  subroutine exact_vm2_terms(m11,m22,n,terms)
    real(dp),intent(in)::m11(:,:),m22(:,:)
    integer,intent(in)::n
    real(dp),intent(out)::terms(3)
    integer::p,i,j
    real(dp)::temp,rn
    p=size(m11,1);rn=real(n,dp);terms=0.0_dp
    if(n<=0)return
    do i=1,p
      do j=i,p
        if(i==j)then
          temp=m22(i,i)-m11(i,i)*m11(i,i)
          terms(1)=terms(1)+temp/rn
          terms(3)=terms(3)+temp/rn
        else
          terms(1)=terms(1)+2.0_dp*(m22(i,j)-m11(i,j)*m11(i,j))/rn
        end if
      end do
    end do
    terms(2)=terms(3)
    do i=1,p
      do j=i+1,p
        terms(2)=terms(2)+2.0_dp*(m22(i,j)-m11(i,i)*m11(j,j))/rn
      end do
    end do
    terms(2)=terms(2)/real(p,dp)
  end subroutine exact_vm2_terms

  real(dp) function cm2_one_factor(x,fc,fvar,m11,m22) result(value)
    real(dp),intent(in)::x(:,:),fc(:),fvar,m11(:,:),m22(:,:)
    real(dp),allocatable::covxf(:)
    real(dp)::s211,s121,s112,temp_i,temp_j,temp_var,rn,fvar2
    integer::n,p,i,j
    n=size(x,1);p=size(x,2);rn=real(n,dp);value=0.0_dp
    if(n<=0 .or. fvar<=tiny(1.0_dp))return
    allocate(covxf(p));fvar2=fvar*fvar
    do i=1,p;covxf(i)=sum(x(:,i)*fc)/rn;end do
    do i=1,p
      do j=i,p
        if(i==j)then
          value=value+(m22(i,i)-m11(i,i)*m11(i,i))/rn
        else
          s211=sum(x(:,i)**2*x(:,j)*fc)
          s121=sum(x(:,i)*x(:,j)**2*fc)
          s112=sum(x(:,i)*x(:,j)*fc**2)
          temp_i=s211/rn-m11(i,j)*covxf(i)
          temp_j=s121/rn-m11(i,j)*covxf(j)
          temp_var=s112/rn-m11(i,j)*fvar
          value=value+2.0_dp*(covxf(j)*temp_i/fvar+covxf(i)*temp_j/fvar- &
            covxf(i)*covxf(j)*temp_var/fvar2)/rn
        end if
      end do
    end do
  end function cm2_one_factor

  real(dp) function cm2_constant_correlation(x,rcoef,m11,m22) result(value)
    real(dp),intent(in)::x(:,:),rcoef,m11(:,:),m22(:,:)
    real(dp)::s31,s13,temp_i,temp_j,rn
    integer::n,p,i,j
    n=size(x,1);p=size(x,2);rn=real(n,dp);value=0.0_dp
    if(n<=0)return
    do i=1,p
      do j=i,p
        if(i==j)then
          value=value+(m22(i,i)-m11(i,i)*m11(i,i))/rn
        else if(m11(i,i)>tiny(1.0_dp) .and. m11(j,j)>tiny(1.0_dp))then
          s31=sum(x(:,i)**3*x(:,j));s13=sum(x(:,i)*x(:,j)**3)
          temp_i=s31/rn-m11(i,j)*m11(i,i)
          temp_j=s13/rn-m11(i,j)*m11(j,j)
          value=value+rcoef*(sqrt(m11(j,j)/m11(i,i))*temp_i+ &
            sqrt(m11(i,i)/m11(j,j))*temp_j)/rn
        end if
      end do
    end do
  end function cm2_constant_correlation

  subroutine m3_target_independent(margskew,target)
    real(dp),intent(in)::margskew(:)
    real(dp),intent(out)::target(:)
    integer::p,i,j,k,q
    p=size(margskew);q=0
    do i=1,p;do j=i,p;do k=j,p
      q=q+1;target(q)=0.0_dp;if(i==j .and. j==k)target(q)=margskew(i)
    end do;end do;end do
  end subroutine m3_target_independent

  subroutine m3_target_simaan(margroot,target)
    real(dp),intent(in)::margroot(:)
    real(dp),intent(out)::target(:)
    integer::p,i,j,k,q
    p=size(margroot);q=0
    do i=1,p;do j=i,p;do k=j,p
      q=q+1;target(q)=margroot(i)*margroot(j)*margroot(k)
    end do;end do;end do
  end subroutine m3_target_simaan

  subroutine m3_target_one_factor(margskew,beta,fskew,target)
    real(dp),intent(in)::margskew(:),beta(:),fskew
    real(dp),intent(out)::target(:)
    integer::p,i,j,k,q
    p=size(margskew);q=0
    do i=1,p;do j=i,p;do k=j,p
      q=q+1
      if(i==j .and. j==k)then;target(q)=margskew(i)
      else;target(q)=beta(i)*beta(j)*beta(k)*fskew;end if
    end do;end do;end do
  end subroutine m3_target_one_factor

  subroutine m3_cc_coefficients(x,margvar,margkurt,m21,m22,coef)
    real(dp),intent(in)::x(:,:),margvar(:),margkurt(:),m21(:,:),m22(:,:)
    real(dp),intent(out)::coef(3)
    real(dp)::rn,rp,m111,nc
    integer::n,p,i,j,k
    n=size(x,1);p=size(x,2);rn=real(n,dp);rp=real(p,dp);coef=0.0_dp
    if(n<=0 .or. p<2)return
    do i=1,p
      do j=i+1,p
        if(margkurt(i)>0.0_dp .and. margvar(j)>0.0_dp)coef(1)=coef(1)+m21(i,j)/sqrt(margkurt(i)*margvar(j))
        if(margkurt(i)>0.0_dp .and. margkurt(j)>0.0_dp)coef(3)=coef(3)+m22(i,j)/sqrt(margkurt(i)*margkurt(j))
      end do
    end do
    coef(1)=coef(1)*2.0_dp/(rp*(rp-1.0_dp));coef(3)=coef(3)*2.0_dp/(rp*(rp-1.0_dp))
    if(p<3)return
    do i=1,p;do j=i+1,p;do k=j+1,p
      m111=sum(x(:,i)*x(:,j)*x(:,k))/rn
      nc=0.0_dp
      if(margvar(i)>0.0_dp .and. margkurt(j)>0.0_dp .and. margkurt(k)>0.0_dp .and. coef(3)>=0.0_dp) &
        nc=nc+sqrt(margvar(i)*coef(3)*sqrt(margkurt(j)*margkurt(k)))
      if(margvar(j)>0.0_dp .and. margkurt(i)>0.0_dp .and. margkurt(k)>0.0_dp .and. coef(3)>=0.0_dp) &
        nc=nc+sqrt(margvar(j)*coef(3)*sqrt(margkurt(i)*margkurt(k)))
      if(margvar(k)>0.0_dp .and. margkurt(i)>0.0_dp .and. margkurt(j)>0.0_dp .and. coef(3)>=0.0_dp) &
        nc=nc+sqrt(margvar(k)*coef(3)*sqrt(margkurt(i)*margkurt(j)))
      nc=nc/3.0_dp
      if(abs(nc)>tiny(1.0_dp))coef(2)=coef(2)+m111/nc
    end do;end do;end do
    coef(2)=coef(2)*6.0_dp/(rp*(rp-1.0_dp)*(rp-2.0_dp))
  end subroutine m3_cc_coefficients

  subroutine m3_target_cc(margvar,margskew,margkurt,coef,target)
    real(dp),intent(in)::margvar(:),margskew(:),margkurt(:),coef(3)
    real(dp),intent(out)::target(:)
    real(dp)::r2,r4,r5
    integer::p,i,j,k,q
    p=size(margvar);r2=coef(1);r4=coef(2);r5=max(coef(3),0.0_dp);q=0
    do i=1,p;do j=i,p;do k=j,p
      q=q+1
      if(i==j)then
        if(j==k)then;target(q)=margskew(i)
        else;target(q)=r2*sqrt(max(margvar(k)*margkurt(i),0.0_dp));end if
      else if(j==k)then
        target(q)=r2*sqrt(max(margvar(i)*margkurt(j),0.0_dp))
      else
        target(q)=r4*sqrt(r5)*(sqrt(max(margvar(k)*sqrt(max(margkurt(i)*margkurt(j),0.0_dp)),0.0_dp))+ &
          sqrt(max(margvar(j)*sqrt(max(margkurt(i)*margkurt(k),0.0_dp)),0.0_dp))+ &
          sqrt(max(margvar(i)*sqrt(max(margkurt(j)*margkurt(k),0.0_dp)),0.0_dp)))/3.0_dp
      end if
    end do;end do;end do
  end subroutine m3_target_cc


  subroutine exact_vm3_terms(x,m11,m21,m22,m31,m42,m33,terms)
    real(dp),intent(in)::x(:,:),m11(:,:),m21(:,:),m22(:,:),m31(:,:),m42(:,:),m33(:,:)
    real(dp),intent(out)::terms(3)
    real(dp)::rn,s211,s121,s112,s111,s222,temp
    integer::n,p,i,j,k
    n=size(x,1);p=size(x,2);rn=real(n,dp);terms=0.0_dp
    if(n<=0)return
    do i=1,p
      do j=i,p
        do k=j,p
          if(i==j)then
            if(j==k)then
              temp=(m42(i,i)-m21(i,i)**2-6.0_dp*m22(i,i)*m11(i,i)+9.0_dp*m11(i,i)**3)/rn
              terms(1)=terms(1)+temp;terms(3)=terms(3)+temp
            else
              terms(1)=terms(1)+3.0_dp*(m42(i,k)-m21(i,k)**2-4.0_dp*m31(i,k)*m11(i,k)- &
                2.0_dp*m22(i,k)*m11(i,i)+8.0_dp*m11(i,i)*m11(i,k)**2+m11(k,k)*m11(i,i)**2)/rn
            end if
          else if(j==k)then
            terms(1)=terms(1)+3.0_dp*(m42(j,i)-m21(j,i)**2-4.0_dp*m31(j,i)*m11(i,j)- &
              2.0_dp*m22(i,j)*m11(j,j)+8.0_dp*m11(j,j)*m11(i,j)**2+m11(i,i)*m11(j,j)**2)/rn
          else
            s211=sum(x(:,i)**2*x(:,j)*x(:,k));s121=sum(x(:,i)*x(:,j)**2*x(:,k))
            s112=sum(x(:,i)*x(:,j)*x(:,k)**2);s111=sum(x(:,i)*x(:,j)*x(:,k))
            s222=sum(x(:,i)**2*x(:,j)**2*x(:,k)**2)
            terms(1)=terms(1)+6.0_dp*(s222/rn-s111*s111/(rn*rn)-2.0_dp*s211/rn*m11(j,k)- &
              2.0_dp*s121/rn*m11(i,k)-2.0_dp*s112/rn*m11(i,j)+6.0_dp*m11(i,k)*m11(j,k)*m11(i,j)+ &
              m11(i,i)*m11(j,k)**2+m11(j,j)*m11(i,k)**2+m11(k,k)*m11(i,j)**2)/rn
          end if
        end do
      end do
    end do
    terms(2)=terms(3)
    do i=1,p
      do j=i+1,p
        terms(2)=terms(2)+2.0_dp*(m33(i,j)-m21(i,i)*m21(j,j)-3.0_dp*m31(i,j)*m11(j,j)- &
          3.0_dp*m31(j,i)*m11(i,i)+9.0_dp*m11(i,i)*m11(j,j)*m11(i,j))/rn
      end do
    end do
    terms(2)=terms(2)/real(p,dp)
  end subroutine exact_vm3_terms

  subroutine exact_vm3_kstat_terms(x,s11,s21,s22,s31,s42,s33,terms)
    real(dp),intent(in)::x(:,:),s11(:,:),s21(:,:),s22(:,:),s31(:,:),s42(:,:),s33(:,:)
    real(dp),intent(out)::terms(3)
    real(dp)::rn,n2,n3,n4,n5,n122,alpha
    real(dp)::ciii_s2_3,ciii_s4s2,ciii_s3_2,ciii_s6
    real(dp)::ciij_s02s20_2,ciij_s20s11_2,ciij_s20s22,ciij_s31s11,ciij_s40s02,ciij_s21_2,ciij_s12s30
    real(dp)::cijk_s002s110_2,cijk_s011s101s110,cijk_s112s110,cijk_s002s020s200
    real(dp)::cijk_s022s200,cijk_s111_2,cijk_s102s120
    real(dp)::c_s11_3,c_s02s20s11,c_s22s11,c_s13s20,c_s21s12,c_s03s30
    real(dp)::s211i,s211j,s211k,s111,s222,temp
    integer::n,p,i,j,k
    n=size(x,1);p=size(x,2);terms=0.0_dp
    if(n<6)return
    rn=real(n,dp);n2=rn*rn;n3=n2*rn;n4=n3*rn;n5=n4*rn
    n122=(rn-1.0_dp)**2*(rn-2.0_dp)**2
    alpha=rn*n122*(rn-3.0_dp)*(rn-4.0_dp)*(rn-5.0_dp)
    ciii_s2_3=(9.0_dp*n4-72.0_dp*n3+213.0_dp*n2-270.0_dp*rn+120.0_dp)/alpha
    ciii_s4s2=(-6.0_dp*n5+33.0_dp*n4-42.0_dp*n3-75.0_dp*n2+210.0_dp*rn-120.0_dp)/alpha
    ciii_s3_2=(-n5-4.0_dp*n4+41.0_dp*n3-40.0_dp*n2-100.0_dp*rn+80.0_dp)/alpha
    ciii_s6=(n5-5.0_dp*n4+13.0_dp*n3-23.0_dp*n2+22.0_dp*rn-8.0_dp)/ &
      (n122*(rn-3.0_dp)*(rn-4.0_dp)*(rn-5.0_dp))
    ciij_s02s20_2=(n4-8.0_dp*n3+29.0_dp*n2-46.0_dp*rn+24.0_dp)/alpha
    ciij_s20s11_2=(8.0_dp*n4-64.0_dp*n3+184.0_dp*n2-224.0_dp*rn+96.0_dp)/alpha
    ciij_s20s22=(-2.0_dp*n5+10.0_dp*n4-10.0_dp*n3-34.0_dp*n2+84.0_dp*rn-48.0_dp)/alpha
    ciij_s31s11=(-4.0_dp*n5+24.0_dp*n4-36.0_dp*n3-32.0_dp*n2+112.0_dp*rn-64.0_dp)/alpha
    ciij_s40s02=(-n4+4.0_dp*n3-9.0_dp*n2+14.0_dp*rn-8.0_dp)/alpha
    ciij_s21_2=(-n5+25.0_dp*n3-36.0_dp*n2-60.0_dp*rn+48.0_dp)/alpha
    ciij_s12s30=(-4.0_dp*n4+16.0_dp*n3-4.0_dp*n2-40.0_dp*rn+32.0_dp)/alpha
    cijk_s002s110_2=(n4-8.0_dp*n3+25.0_dp*n2-34.0_dp*rn+16.0_dp)/alpha
    cijk_s011s101s110=(6.0_dp*n4-48.0_dp*n3+134.0_dp*n2-156.0_dp*rn+64.0_dp)/alpha
    cijk_s112s110=(-2.0_dp*n5+12.0_dp*n4-18.0_dp*n3-16.0_dp*n2+56.0_dp*rn-32.0_dp)/alpha
    cijk_s002s020s200=(4.0_dp*n2-12.0_dp*rn+8.0_dp)/alpha
    cijk_s022s200=(-n4+4.0_dp*n3-9.0_dp*n2+14.0_dp*rn-8.0_dp)/alpha
    cijk_s111_2=(-n5+2.0_dp*n4+17.0_dp*n3-34.0_dp*n2-40.0_dp*rn+32.0_dp)/alpha
    cijk_s102s120=(-2.0_dp*n4+8.0_dp*n3-2.0_dp*n2-20.0_dp*rn+16.0_dp)/alpha
    do i=1,p
      do j=i,p
        do k=j,p
          if(i==j)then
            if(j==k)then
              temp=ciii_s2_3*s11(i,i)**3+ciii_s4s2*s22(i,i)*s11(i,i)+ciii_s3_2*s21(i,i)**2+ciii_s6*s42(i,i)
              terms(1)=terms(1)+temp;terms(3)=terms(3)+temp
            else
              terms(1)=terms(1)+3.0_dp*(ciij_s02s20_2*s11(k,k)*s11(i,i)**2+ &
                ciij_s20s11_2*s11(i,i)*s11(i,k)**2+ciij_s20s22*s11(i,i)*s22(i,k)+ &
                ciij_s31s11*s31(i,k)*s11(i,k)+ciij_s40s02*s22(i,i)*s11(k,k)+ &
                ciij_s21_2*s21(i,k)**2+ciij_s12s30*s21(k,i)*s21(i,i)+ciii_s6*s42(i,k))
            end if
          else if(j==k)then
            terms(1)=terms(1)+3.0_dp*(ciij_s02s20_2*s11(i,i)*s11(j,j)**2+ &
              ciij_s20s11_2*s11(j,j)*s11(i,j)**2+ciij_s20s22*s11(j,j)*s22(i,j)+ &
              ciij_s31s11*s31(j,i)*s11(i,j)+ciij_s40s02*s22(j,j)*s11(i,i)+ &
              ciij_s21_2*s21(j,i)**2+ciij_s12s30*s21(i,j)*s21(j,j)+ciii_s6*s42(j,i))
          else
            s211i=sum(x(:,i)**2*x(:,j)*x(:,k));s211j=sum(x(:,i)*x(:,j)**2*x(:,k))
            s211k=sum(x(:,i)*x(:,j)*x(:,k)**2);s111=sum(x(:,i)*x(:,j)*x(:,k))
            s222=sum(x(:,i)**2*x(:,j)**2*x(:,k)**2)
            terms(1)=terms(1)+6.0_dp*(cijk_s002s110_2*s11(k,k)*s11(i,j)**2+ &
              (cijk_s011s101s110*s11(j,k)*s11(i,k)+cijk_s112s110*s211k)*s11(i,j)+ &
              cijk_s002s110_2*s11(j,j)*s11(i,k)**2+cijk_s112s110*s211j*s11(i,k)+ &
              cijk_s002s110_2*s11(i,i)*s11(j,k)**2+cijk_s112s110*s211i*s11(j,k)+ &
              (cijk_s002s020s200*s11(j,j)*s11(k,k)+cijk_s022s200*s22(j,k))*s11(i,i)+ &
              cijk_s022s200*s22(i,k)*s11(j,j)+cijk_s022s200*s22(i,j)*s11(k,k)+ &
              cijk_s111_2*s111**2+cijk_s102s120*s21(k,i)*s21(j,i)+ &
              cijk_s102s120*s21(k,j)*s21(i,j)+cijk_s102s120*s21(j,k)*s21(i,k)+ciii_s6*s222)
          end if
        end do
      end do
    end do
    terms(2)=terms(3)
    c_s11_3=(24.0_dp*n2-72.0_dp*rn+48.0_dp)/alpha
    c_s02s20s11=(9.0_dp*n4-72.0_dp*n3+189.0_dp*n2-198.0_dp*rn+72.0_dp)/alpha
    c_s22s11=(-9.0_dp*n4+36.0_dp*n3-81.0_dp*n2+126.0_dp*rn-72.0_dp)/alpha
    c_s13s20=(-3.0_dp*n5+21.0_dp*n4-39.0_dp*n3+3.0_dp*n2+42.0_dp*rn-24.0_dp)/alpha
    c_s21s12=(-9.0_dp*n4+36.0_dp*n3-9.0_dp*n2-90.0_dp*rn+72.0_dp)/alpha
    c_s03s30=(-n5+5.0_dp*n4+5.0_dp*n3-31.0_dp*n2-10.0_dp*rn+8.0_dp)/alpha
    do i=1,p
      do j=i+1,p
        terms(2)=terms(2)+2.0_dp*(c_s11_3*s11(i,j)**3+ &
          (c_s02s20s11*s11(i,i)*s11(j,j)+c_s22s11*s22(i,j))*s11(i,j)+ &
          c_s13s20*s31(j,i)*s11(i,i)+c_s13s20*s31(i,j)*s11(j,j)+ &
          c_s21s12*s21(i,j)*s21(j,i)+c_s03s30*s21(i,i)*s21(j,j)+ciii_s6*s33(i,j))
      end do
    end do
    terms(2)=terms(2)/real(p,dp)
  end subroutine exact_vm3_kstat_terms


  real(dp) function cm3_simaan(x,marg_skew_inv23,m11,m21,m22,m31,m42,m51) result(value)
    real(dp),intent(in)::x(:,:),marg_skew_inv23(:),m11(:,:),m21(:,:),m22(:,:),m31(:,:),m42(:,:),m51(:,:)
    real(dp)::rn,s411,s141,s114,s111,s211,s121,s112,temp_i,temp_j,temp_k
    integer::n,p,i,j,k
    n=size(x,1);p=size(x,2);rn=real(n,dp);value=0.0_dp
    if(n<=0)return
    do i=1,p
      do j=i,p
        do k=j,p
          if(i==j .and. j==k)then
            value=value+(m42(i,i)-m21(i,i)**2-6.0_dp*m22(i,i)*m11(i,i)+9.0_dp*m11(i,i)**3)
          else if(i==j)then
            temp_i=m51(i,k)-m21(i,k)*m21(i,i)-4.0_dp*m31(i,k)*m11(i,i)- &
              2.0_dp*m22(i,i)*m11(i,k)+9.0_dp*m11(i,i)**2*m11(i,k)
            temp_k=m42(k,i)-m21(i,k)*m21(k,k)-3.0_dp*m22(i,k)*m11(k,k)- &
              m22(k,k)*m11(i,i)-2.0_dp*m31(k,i)*m11(i,k)+6.0_dp*m11(k,k)*m11(i,k)**2+ &
              3.0_dp*m11(k,k)**2*m11(i,i)
            value=value+marg_skew_inv23(i)**2*marg_skew_inv23(k)* &
              (2.0_dp*m21(i,i)*m21(k,k)*temp_i+m21(i,i)**2*temp_k)
          else if(j==k)then
            temp_i=m42(i,j)-m21(j,i)*m21(i,i)-3.0_dp*m22(i,j)*m11(i,i)- &
              m22(i,i)*m11(j,j)-2.0_dp*m31(j,i)*m11(i,j)+6.0_dp*m11(i,i)*m11(i,j)**2+ &
              3.0_dp*m11(i,i)**2*m11(j,j)
            temp_j=m51(j,i)-m21(j,i)*m21(j,j)-4.0_dp*m31(j,i)*m11(j,j)- &
              2.0_dp*m22(j,j)*m11(i,j)+9.0_dp*m11(j,j)**2*m11(i,j)
            value=value+marg_skew_inv23(i)*marg_skew_inv23(j)**2* &
              (m21(j,j)**2*temp_i+2.0_dp*m21(i,i)*m21(j,j)*temp_j)
          else
            s411=sum(x(:,i)**4*x(:,j)*x(:,k));s141=sum(x(:,i)*x(:,j)**4*x(:,k))
            s114=sum(x(:,i)*x(:,j)*x(:,k)**4);s111=sum(x(:,i)*x(:,j)*x(:,k))
            s211=sum(x(:,i)**2*x(:,j)*x(:,k));s121=sum(x(:,i)*x(:,j)**2*x(:,k))
            s112=sum(x(:,i)*x(:,j)*x(:,k)**2)
            temp_i=s411/rn-s111/rn*m21(i,i)-3.0_dp*s211/rn*m11(i,i)-m22(i,i)*m11(j,k)- &
              m31(i,j)*m11(i,k)-m31(i,k)*m11(i,j)+6.0_dp*m11(i,i)*m11(i,j)*m11(i,k)+ &
              3.0_dp*m11(i,i)**2*m11(j,k)
            temp_j=s141/rn-s111/rn*m21(j,j)-3.0_dp*s121/rn*m11(j,j)-m22(j,j)*m11(i,k)- &
              m31(j,i)*m11(j,k)-m31(j,k)*m11(i,j)+6.0_dp*m11(j,j)*m11(i,j)*m11(j,k)+ &
              3.0_dp*m11(j,j)**2*m11(i,k)
            temp_k=s114/rn-s111/rn*m21(k,k)-3.0_dp*s112/rn*m11(k,k)-m22(k,k)*m11(i,j)- &
              m31(k,j)*m11(i,k)-m31(k,i)*m11(j,k)+6.0_dp*m11(k,k)*m11(j,k)*m11(i,k)+ &
              3.0_dp*m11(k,k)**2*m11(i,j)
            value=value+2.0_dp*marg_skew_inv23(i)*marg_skew_inv23(j)*marg_skew_inv23(k)* &
              (m21(j,j)*m21(k,k)*temp_i+m21(i,i)*m21(k,k)*temp_j+m21(i,i)*m21(j,j)*temp_k)
          end if
        end do
      end do
    end do
    value=value/rn
  end function cm3_simaan

  real(dp) function cm3_one_factor(x,fc,fvar,fskew,m11,m21,m22,m42) result(value)
    real(dp),intent(in)::x(:,:),fc(:),fvar,fskew,m11(:,:),m21(:,:),m22(:,:),m42(:,:)
    real(dp),allocatable::covxf(:),x1f2(:),x1f3(:),x11f1(:,:)
    real(dp)::rn,fvar3,elem
    real(dp)::s311,s221,s213,s211,s212,s131,s123,s121,s122
    real(dp)::s2111,s1211,s1121,s1112,s1111,s1113,s1110
    real(dp)::temp_i,temp_j,temp_k,temp_fs,temp_fv
    integer::n,p,i,j,k
    n=size(x,1);p=size(x,2);rn=real(n,dp);value=0.0_dp
    if(n<=0 .or. fvar<=tiny(1.0_dp))return
    fvar3=fvar**3
    allocate(covxf(p),x1f2(p),x1f3(p),x11f1(p,p))
    do i=1,p
      covxf(i)=sum(x(:,i)*fc)/rn;x1f2(i)=sum(x(:,i)*fc**2)/rn;x1f3(i)=sum(x(:,i)*fc**3)/rn
      do j=i,p
        elem=sum(x(:,i)*x(:,j)*fc)/rn;x11f1(i,j)=elem;x11f1(j,i)=elem
      end do
    end do
    do i=1,p
      do j=i,p
        do k=j,p
          if(i==j .and. j==k)then
            value=value+m42(i,i)-m21(i,i)**2-6.0_dp*m22(i,i)*m11(i,i)+9.0_dp*m11(i,i)**3
          else if(i==j)then
            s311=sum(x(:,i)**3*x(:,k)*fc);s221=sum(x(:,i)**2*x(:,k)**2*fc)
            s213=sum(x(:,i)**2*x(:,k)*fc**3);s211=sum(x(:,i)**2*x(:,k)*fc)
            s212=sum(x(:,i)**2*x(:,k)*fc**2)
            temp_i=s311/rn-m21(i,k)*covxf(i)-2.0_dp*m11(i,k)*x11f1(i,i)-m11(i,i)*x11f1(i,k)
            temp_k=s221/rn-m21(i,k)*covxf(k)-m11(i,i)*x11f1(k,k)-2.0_dp*m11(i,k)*x11f1(i,k)
            temp_fs=s213/rn-m21(i,k)*fskew-3.0_dp*s211/rn*fvar-2.0_dp*x1f3(i)*m11(i,k)- &
              x1f3(k)*m11(i,i)+3.0_dp*m11(i,i)*covxf(k)*fvar+6.0_dp*m11(i,k)*covxf(i)*fvar
            temp_fv=s212/rn-m21(i,k)*fvar-2.0_dp*x1f2(i)*m11(i,k)-x1f2(k)*m11(i,i)
            value=value+3.0_dp*((2.0_dp*covxf(i)*covxf(k)*temp_i+covxf(i)**2*temp_k)*fskew+ &
              covxf(i)**2*covxf(k)*temp_fs-3.0_dp*covxf(i)**2*covxf(k)*fskew*temp_fv/fvar)/fvar3
          else if(j==k)then
            s221=sum(x(:,i)**2*x(:,j)**2*fc);s131=sum(x(:,i)*x(:,j)**3*fc)
            s123=sum(x(:,i)*x(:,j)**2*fc**3);s121=sum(x(:,i)*x(:,j)**2*fc)
            s122=sum(x(:,i)*x(:,j)**2*fc**2)
            temp_i=s221/rn-m21(j,i)*covxf(i)-m11(j,j)*x11f1(i,i)-2.0_dp*m11(i,j)*x11f1(i,j)
            temp_j=s131/rn-m21(j,i)*covxf(j)-2.0_dp*m11(i,j)*x11f1(j,j)-m11(j,j)*x11f1(i,j)
            temp_fs=s123/rn-m21(j,i)*fskew-3.0_dp*s121/rn*fvar-x1f3(i)*m11(j,j)- &
              2.0_dp*x1f3(j)*m11(i,j)+6.0_dp*m11(i,j)*covxf(j)*fvar+3.0_dp*m11(j,j)*covxf(i)*fvar
            temp_fv=s122/rn-m21(j,i)*fvar-x1f2(i)*m11(j,j)-2.0_dp*x1f2(j)*m11(i,j)
            value=value+3.0_dp*((covxf(j)**2*temp_i+2.0_dp*covxf(i)*covxf(j)*temp_j)*fskew+ &
              covxf(i)*covxf(j)**2*temp_fs-3.0_dp*covxf(i)*covxf(j)**2*fskew*temp_fv/fvar)/fvar3
          else
            s2111=sum(x(:,i)**2*x(:,j)*x(:,k)*fc);s1211=sum(x(:,i)*x(:,j)**2*x(:,k)*fc)
            s1121=sum(x(:,i)*x(:,j)*x(:,k)**2*fc);s1112=sum(x(:,i)*x(:,j)*x(:,k)*fc**2)
            s1113=sum(x(:,i)*x(:,j)*x(:,k)*fc**3);s1111=sum(x(:,i)*x(:,j)*x(:,k)*fc)
            s1110=sum(x(:,i)*x(:,j)*x(:,k))
            temp_i=s2111/rn-s1110/rn*covxf(i)-m11(j,k)*x11f1(i,i)-m11(i,k)*x11f1(i,j)-m11(i,j)*x11f1(i,k)
            temp_j=s1211/rn-s1110/rn*covxf(j)-m11(i,k)*x11f1(j,j)-m11(j,k)*x11f1(i,j)-m11(i,j)*x11f1(j,k)
            temp_k=s1121/rn-s1110/rn*covxf(k)-m11(i,j)*x11f1(k,k)-m11(j,k)*x11f1(i,k)-m11(i,k)*x11f1(j,k)
            temp_fs=s1113/rn-s1110/rn*fskew-3.0_dp*s1111/rn*fvar-x1f3(i)*m11(j,k)- &
              x1f3(j)*m11(i,k)-x1f3(k)*m11(i,j)+3.0_dp*m11(i,j)*covxf(k)*fvar+ &
              3.0_dp*m11(i,k)*covxf(j)*fvar+3.0_dp*m11(j,k)*covxf(i)*fvar
            temp_fv=s1112/rn-s1110/rn*fvar-x1f2(i)*m11(j,k)-x1f2(j)*m11(i,k)-x1f2(k)*m11(i,j)
            value=value+6.0_dp*((covxf(j)*covxf(k)*temp_i+covxf(i)*covxf(k)*temp_j+ &
              covxf(i)*covxf(j)*temp_k)*fskew+covxf(i)*covxf(j)*covxf(k)*temp_fs- &
              3.0_dp*covxf(i)*covxf(j)*covxf(k)*fskew*temp_fv/fvar)/fvar3
          end if
        end do
      end do
    end do
    value=value/rn
  end function cm3_one_factor


  real(dp) function cm3_constant_correlation(x,margvar,margskew,margkurt,marg5,marg6, &
      m11,m21,m31,m32,m41,m61,coef) result(value)
    real(dp),intent(in)::x(:,:),margvar(:),margskew(:),margkurt(:),marg5(:),marg6(:)
    real(dp),intent(in)::m11(:,:),m21(:,:),m31(:,:),m32(:,:),m41(:,:),m61(:,:),coef(3)
    real(dp)::rn,r2,r4,r5,s511,s151,s115,s211,s121,s112,s311,s131,s113,m111
    real(dp)::ti4,tj4,tk4,ti2,tj2,tk2,ci,cj,ck
    integer::n,p,i,j,k
    n=size(x,1);p=size(x,2);rn=real(n,dp);r2=coef(1);r4=coef(2);r5=max(coef(3),0.0_dp);value=0.0_dp
    if(n<=0)return
    do i=1,p
      do j=i,p
        do k=j,p
          if(i==j .and. j==k)then
            value=value+marg6(i)-margskew(i)**2-6.0_dp*margkurt(i)*margvar(i)+9.0_dp*margvar(i)**3
          else if(i==j)then
            ti4=m61(i,k)-m21(i,k)*margkurt(i)-4.0_dp*m31(i,k)*margskew(i)- &
              2.0_dp*m11(i,k)*marg5(i)-margvar(i)*m41(i,k)+12.0_dp*margvar(i)*m11(i,k)*margskew(i)
            tk2=m32(k,i)-margvar(k)*m21(i,k)-margskew(k)*margvar(i)-2.0_dp*m21(k,i)*m11(i,k)
            if(margkurt(i)>0.0_dp .and. margvar(k)>0.0_dp) value=value+3.0_dp*r2* &
              (sqrt(margvar(k)/margkurt(i))*ti4+sqrt(margkurt(i)/margvar(k))*tk2)/2.0_dp
          else if(j==k)then
            ti2=m32(i,j)-margvar(i)*m21(j,i)-margskew(i)*margvar(j)-2.0_dp*m21(i,j)*m11(i,j)
            tj4=m61(j,i)-m21(j,i)*margkurt(j)-4.0_dp*m31(j,i)*margskew(j)- &
              2.0_dp*m11(i,j)*marg5(j)-margvar(j)*m41(j,i)+12.0_dp*margvar(j)*m11(i,j)*margskew(j)
            if(margkurt(j)>0.0_dp .and. margvar(i)>0.0_dp) value=value+3.0_dp*r2* &
              (sqrt(margvar(i)/margkurt(j))*tj4+sqrt(margkurt(j)/margvar(i))*ti2)/2.0_dp
          else
            s511=sum(x(:,i)**5*x(:,j)*x(:,k));s151=sum(x(:,i)*x(:,j)**5*x(:,k));s115=sum(x(:,i)*x(:,j)*x(:,k)**5)
            s211=sum(x(:,i)**2*x(:,j)*x(:,k));s121=sum(x(:,i)*x(:,j)**2*x(:,k));s112=sum(x(:,i)*x(:,j)*x(:,k)**2)
            s311=sum(x(:,i)**3*x(:,j)*x(:,k));s131=sum(x(:,i)*x(:,j)**3*x(:,k));s113=sum(x(:,i)*x(:,j)*x(:,k)**3)
            m111=sum(x(:,i)*x(:,j)*x(:,k))/rn
            ti4=s511/rn-m111*margkurt(i)-4.0_dp*margskew(i)*s211/rn-marg5(i)*m11(j,k)- &
              m41(i,j)*m11(i,k)-m41(i,k)*m11(i,j)+8.0_dp*margskew(i)*m11(i,j)*m11(i,k)+ &
              4.0_dp*margskew(i)*margvar(i)*m11(j,k)
            tj4=s151/rn-m111*margkurt(j)-4.0_dp*margskew(j)*s121/rn-marg5(j)*m11(i,k)- &
              m41(j,k)*m11(i,j)-m41(j,i)*m11(j,k)+8.0_dp*margskew(j)*m11(j,k)*m11(i,j)+ &
              4.0_dp*margskew(j)*margvar(j)*m11(i,k)
            tk4=s115/rn-m111*margkurt(k)-4.0_dp*margskew(k)*s112/rn-marg5(k)*m11(i,j)- &
              m41(k,i)*m11(j,k)-m41(k,j)*m11(i,k)+8.0_dp*margskew(k)*m11(i,k)*m11(j,k)+ &
              4.0_dp*margskew(k)*margvar(k)*m11(i,j)
            ti2=s311/rn-margvar(i)*m111-margskew(i)*m11(j,k)-m21(i,j)*m11(i,k)-m21(i,k)*m11(i,j)
            tj2=s131/rn-margvar(j)*m111-margskew(j)*m11(i,k)-m21(j,i)*m11(j,k)-m21(j,k)*m11(i,j)
            tk2=s113/rn-margvar(k)*m111-margskew(k)*m11(i,j)-m21(k,i)*m11(j,k)-m21(k,j)*m11(i,k)
            ci=0.0_dp;cj=0.0_dp;ck=0.0_dp
            if(margvar(i)>0.0_dp .and. margkurt(j)>0.0_dp .and. margkurt(k)>0.0_dp)then
              ci=sqrt(margvar(i)*sqrt(margkurt(j)/margkurt(k)**3))*tk4/2.0_dp+ &
                sqrt(margvar(i)*sqrt(margkurt(k)/margkurt(j)**3))*tj4/2.0_dp+ &
                sqrt(sqrt(margkurt(j)*margkurt(k))/margvar(i))*ti2
            end if
            if(margvar(j)>0.0_dp .and. margkurt(i)>0.0_dp .and. margkurt(k)>0.0_dp)then
              cj=sqrt(margvar(j)*sqrt(margkurt(i)/margkurt(k)**3))*tk4/2.0_dp+ &
                sqrt(margvar(j)*sqrt(margkurt(k)/margkurt(i)**3))*ti4/2.0_dp+ &
                sqrt(sqrt(margkurt(i)*margkurt(k))/margvar(j))*tj2
            end if
            if(margvar(k)>0.0_dp .and. margkurt(i)>0.0_dp .and. margkurt(j)>0.0_dp)then
              ck=sqrt(margvar(k)*sqrt(margkurt(j)/margkurt(i)**3))*ti4/2.0_dp+ &
                sqrt(margvar(k)*sqrt(margkurt(i)/margkurt(j)**3))*tj4/2.0_dp+ &
                sqrt(sqrt(margkurt(i)*margkurt(j))/margvar(k))*tk2
            end if
            value=value+r4*sqrt(r5)*(ci+cj+ck)
          end if
        end do
      end do
    end do
    value=value/rn
  end function cm3_constant_correlation


  subroutine m4_target_independent(marg_kurt,marg_pair,target)
    real(dp),intent(in)::marg_kurt(:),marg_pair(:)
    real(dp),intent(out)::target(:)
    integer::p,i,j,k,l,q
    p=size(marg_kurt);q=0
    do i=1,p;do j=i,p;do k=j,p;do l=k,p
      q=q+1;target(q)=0.0_dp
      if(i==j .and. j==k .and. k==l)target(q)=marg_kurt(i)
      if(i==j .and. k==l .and. j/=k)target(q)=marg_pair(i)*marg_pair(k)
    end do;end do;end do;end do
  end subroutine m4_target_independent

  subroutine m4_cc_coefficients(x,margvar,margkurt,marg6,m22,m31,coef)
    real(dp),intent(in)::x(:,:),margvar(:),margkurt(:),marg6(:),m22(:,:),m31(:,:)
    real(dp),intent(out)::coef(4)
    real(dp)::rn,rp,m211,m1111,den
    integer::n,p,i,j,k,l
    n=size(x,1);p=size(x,2);rn=real(n,dp);rp=real(p,dp);coef=0.0_dp
    if(n<=0 .or. p<2)return
    do i=1,p;do j=i+1,p
      if(marg6(i)>0.0_dp .and. margvar(j)>0.0_dp)coef(1)=coef(1)+m31(i,j)/sqrt(marg6(i)*margvar(j))
      if(margkurt(i)>0.0_dp .and. margkurt(j)>0.0_dp)coef(2)=coef(2)+m22(i,j)/sqrt(margkurt(i)*margkurt(j))
    end do;end do
    coef(1)=coef(1)*2.0_dp/(rp*(rp-1.0_dp));coef(2)=coef(2)*2.0_dp/(rp*(rp-1.0_dp))
    if(p>=3 .and. coef(2)>0.0_dp)then
      do i=1,p;do j=i+1,p;do k=j+1,p
        m211=sum(x(:,i)**2*x(:,j)*x(:,k))/rn
        den=sqrt(max(margkurt(i)*coef(2)*sqrt(max(margkurt(j)*margkurt(k),0.0_dp)),0.0_dp))
        if(den>tiny(1.0_dp))coef(3)=coef(3)+m211/den
      end do;end do;end do
      coef(3)=coef(3)*6.0_dp/(rp*(rp-1.0_dp)*(rp-2.0_dp))
    end if
    if(p>=4 .and. abs(coef(2))>tiny(1.0_dp))then
      do i=1,p;do j=i+1,p;do k=j+1,p;do l=k+1,p
        m1111=sum(x(:,i)*x(:,j)*x(:,k)*x(:,l))/rn
        den=coef(2)*sqrt(sqrt(max(margkurt(i)*margkurt(j)*margkurt(k)*margkurt(l),0.0_dp)))
        if(abs(den)>tiny(1.0_dp))coef(4)=coef(4)+m1111/den
      end do;end do;end do;end do
      coef(4)=coef(4)*24.0_dp/(rp*(rp-1.0_dp)*(rp-2.0_dp)*(rp-3.0_dp))
    end if
  end subroutine m4_cc_coefficients

  subroutine m4_target_cc(margvar,margkurt,marg6,coef,target)
    real(dp),intent(in)::margvar(:),margkurt(:),marg6(:),coef(4)
    real(dp),intent(out)::target(:)
    real(dp)::r3,r5,r6,r7
    integer::p,i,j,k,l,q
    p=size(margvar);r3=coef(1);r5=coef(2);r6=coef(3);r7=coef(4);q=0
    do i=1,p;do j=i,p;do k=j,p;do l=k,p
      q=q+1
      if(i==j)then
        if(j==k)then
          if(k==l)then;target(q)=margkurt(i)
          else;target(q)=r3*sqrt(max(marg6(i)*margvar(l),0.0_dp));end if
        else if(k==l)then
          target(q)=r5*sqrt(max(margkurt(i)*margkurt(k),0.0_dp))
        else
          target(q)=r6*sqrt(max(margvar(i)*r5*sqrt(max(margkurt(k)*margkurt(l),0.0_dp)),0.0_dp))
        end if
      else if(j==k)then
        if(k==l)then;target(q)=r3*sqrt(max(margvar(i)*marg6(j),0.0_dp))
        else;target(q)=r6*sqrt(max(margvar(j)*r5*sqrt(max(margkurt(i)*margkurt(l),0.0_dp)),0.0_dp));end if
      else if(k==l)then
        target(q)=r6*sqrt(max(margvar(k)*r5*sqrt(max(margkurt(i)*margkurt(j),0.0_dp)),0.0_dp))
      else
        target(q)=r7*r5*sqrt(sqrt(max(margkurt(i)*margkurt(j)*margkurt(k)*margkurt(l),0.0_dp)))
      end if
    end do;end do;end do;end do
  end subroutine m4_target_cc

  subroutine m4_target_one_factor(margkurt,fvar,fkurt,epsvar,beta,target)
    real(dp),intent(in)::margkurt(:),fvar,fkurt,epsvar(:),beta(:)
    real(dp),intent(out)::target(:)
    integer::p,i,j,k,l,q
    p=size(margkurt);q=0
    do i=1,p;do j=i,p;do k=j,p;do l=k,p
      q=q+1
      if(i==j)then
        if(j==k)then
          if(k==l)then
            target(q)=margkurt(i)
          else
            target(q)=beta(i)**3*beta(l)*fkurt+3.0_dp*beta(i)*beta(l)*fvar*epsvar(i)
          end if
        else if(k==l)then
          target(q)=beta(i)**2*beta(k)**2*fkurt+fvar*(beta(i)**2*epsvar(k)+beta(k)**2*epsvar(i))+epsvar(i)*epsvar(k)
        else
          target(q)=beta(i)**2*beta(k)*beta(l)*fkurt+beta(k)*beta(l)*fvar*epsvar(i)
        end if
      else if(j==k)then
        if(k==l)then
          target(q)=beta(i)*beta(j)**3*fkurt+3.0_dp*beta(i)*beta(j)*fvar*epsvar(j)
        else
          target(q)=beta(i)*beta(j)**2*beta(l)*fkurt+beta(i)*beta(l)*fvar*epsvar(j)
        end if
      else if(k==l)then
        target(q)=beta(i)*beta(j)*beta(k)**2*fkurt+beta(i)*beta(j)*fvar*epsvar(k)
      else
        target(q)=beta(i)*beta(j)*beta(k)*beta(l)*fkurt
      end if
    end do;end do;end do;end do
  end subroutine m4_target_one_factor


  subroutine exact_vm4_terms(x,m11,m21,m22,m31,m32,m41,m42,terms)
    real(dp),intent(in)::x(:,:),m11(:,:),m21(:,:),m22(:,:),m31(:,:),m32(:,:),m41(:,:),m42(:,:)
    real(dp),intent(out)::terms(3)
    real(dp)::rn,rp,s8,s62,s44,s422,s311,s221,s212,m211,m111
    real(dp)::s2222,s2111,s1211,s1121,s1112,m1111,m0111,m1011,m1101,m1110
    real(dp)::temp,temp_i,temp_k
    integer::n,p,i,j,k,l,m
    n=size(x,1);p=size(x,2);rn=real(n,dp);rp=real(p,dp);terms=0.0_dp
    if(n<=0)return
    do i=1,p
      do j=i,p
        do k=j,p
          do l=k,p
            if(i==j)then
              if(j==k)then
                if(k==l)then
                  s8=sum(x(:,i)**8)
                  temp=(s8/rn-m31(i,i)**2-8.0_dp*m32(i,i)*m21(i,i)+16.0_dp*m11(i,i)*m21(i,i)**2)/rn
                  terms(1)=terms(1)+temp;terms(2)=terms(2)+temp/rp;terms(3)=terms(3)+temp
                else
                  s62=sum(x(:,i)**6*x(:,l)**2)
                  terms(1)=terms(1)+4.0_dp*(s62/rn-m31(i,l)**2-6.0_dp*m41(i,l)*m21(i,l)- &
                    2.0_dp*m32(i,l)*m21(i,i)+6.0_dp*m11(i,l)*m21(i,l)*m21(i,i)+ &
                    6.0_dp*m11(i,i)*m21(i,l)**2+m21(i,i)**2*m11(l,l)+ &
                    3.0_dp*m21(i,l)**2*m11(i,i))/rn
                end if
              else if(k==l)then
                s44=sum(x(:,i)**4*x(:,k)**4)
                terms(1)=terms(1)+6.0_dp*(s44/rn-m22(i,k)**2-4.0_dp*m32(i,k)*m21(k,i)- &
                  4.0_dp*m32(k,i)*m21(i,k)+8.0_dp*m11(i,k)*m21(k,i)*m21(i,k)+ &
                  4.0_dp*m21(i,k)**2*m11(k,k)+4.0_dp*m21(k,i)**2*m11(i,i))/rn
                temp=0.0_dp
                do m=1,p
                  temp=temp+2.0_dp*m11(m,m)*(sum(x(:,i)**2*x(:,k)**2*x(:,m)**2)/rn- &
                    m11(m,m)*m22(i,k)-2.0_dp*m21(m,i)*m21(k,i)-2.0_dp*m21(m,k)*m21(i,k))/rn
                end do
                terms(2)=terms(2)+6.0_dp*temp/rp
                temp_i=m42(i,k)-m11(i,i)*m22(i,k)-2.0_dp*m21(i,i)*m21(k,i)-2.0_dp*m21(i,k)**2
                temp_k=m42(k,i)-m11(k,k)*m22(k,i)-2.0_dp*m21(k,k)*m21(i,k)-2.0_dp*m21(k,i)**2
                terms(3)=terms(3)+6.0_dp*(m11(k,k)*temp_i+m11(i,i)*temp_k)/rn
              else
                s422=sum(x(:,i)**4*x(:,k)**2*x(:,l)**2);s311=sum(x(:,i)**3*x(:,k)*x(:,l))
                s221=sum(x(:,i)**2*x(:,k)**2*x(:,l));s212=sum(x(:,i)**2*x(:,k)*x(:,l)**2)
                m211=sum(x(:,i)**2*x(:,k)*x(:,l))/rn;m111=sum(x(:,i)*x(:,k)*x(:,l))/rn
                terms(1)=terms(1)+12.0_dp*(s422/rn-m211**2-4.0_dp*s311/rn*m111- &
                  2.0_dp*s221/rn*m21(i,l)-2.0_dp*s212/rn*m21(i,k)+ &
                  4.0_dp*m11(i,l)*m21(i,k)*m111+4.0_dp*m11(i,k)*m111*m21(i,l)+ &
                  2.0_dp*m11(k,l)*m21(i,l)*m21(i,k)+4.0_dp*m11(i,i)*m111**2+ &
                  m21(i,k)**2*m11(l,l)+m21(i,l)**2*m11(k,k))/rn
              end if
            else if(j==k)then
              if(k==l)then
                s62=sum(x(:,j)**6*x(:,i)**2)
                terms(1)=terms(1)+4.0_dp*(s62/rn-m31(j,i)**2-6.0_dp*m41(j,i)*m21(j,i)- &
                  2.0_dp*m32(j,i)*m21(j,j)+6.0_dp*m11(i,j)*m21(j,i)*m21(j,j)+ &
                  6.0_dp*m11(j,j)*m21(j,i)**2+m21(j,j)**2*m11(i,i)+ &
                  3.0_dp*m21(j,i)**2*m11(j,j))/rn
              else
                s422=sum(x(:,j)**4*x(:,i)**2*x(:,l)**2);s311=sum(x(:,j)**3*x(:,i)*x(:,l))
                s221=sum(x(:,j)**2*x(:,i)**2*x(:,l));s212=sum(x(:,j)**2*x(:,i)*x(:,l)**2)
                m211=sum(x(:,j)**2*x(:,i)*x(:,l))/rn;m111=sum(x(:,j)*x(:,i)*x(:,l))/rn
                terms(1)=terms(1)+12.0_dp*(s422/rn-m211**2-4.0_dp*s311/rn*m111- &
                  2.0_dp*s221/rn*m21(j,l)-2.0_dp*s212/rn*m21(j,i)+ &
                  4.0_dp*m11(j,l)*m21(j,i)*m111+4.0_dp*m11(i,j)*m111*m21(j,l)+ &
                  2.0_dp*m11(i,l)*m21(j,l)*m21(j,i)+4.0_dp*m11(j,j)*m111**2+ &
                  m21(j,i)**2*m11(l,l)+m21(j,l)**2*m11(i,i))/rn
              end if
            else if(k==l)then
              s422=sum(x(:,k)**4*x(:,i)**2*x(:,j)**2);s311=sum(x(:,k)**3*x(:,i)*x(:,j))
              s221=sum(x(:,k)**2*x(:,i)**2*x(:,j));s212=sum(x(:,k)**2*x(:,i)*x(:,j)**2)
              m211=sum(x(:,k)**2*x(:,i)*x(:,j))/rn;m111=sum(x(:,k)*x(:,i)*x(:,j))/rn
              terms(1)=terms(1)+12.0_dp*(s422/rn-m211**2-4.0_dp*s311/rn*m111- &
                2.0_dp*s221/rn*m21(k,j)-2.0_dp*s212/rn*m21(k,i)+ &
                4.0_dp*m11(j,k)*m21(k,i)*m111+4.0_dp*m11(i,k)*m111*m21(k,j)+ &
                2.0_dp*m11(i,j)*m21(k,j)*m21(k,i)+4.0_dp*m11(k,k)*m111**2+ &
                m21(k,i)**2*m11(j,j)+m21(k,j)**2*m11(i,i))/rn
            else
              s2222=sum(x(:,i)**2*x(:,j)**2*x(:,k)**2*x(:,l)**2)
              s2111=sum(x(:,i)**2*x(:,j)*x(:,k)*x(:,l));s1211=sum(x(:,i)*x(:,j)**2*x(:,k)*x(:,l))
              s1121=sum(x(:,i)*x(:,j)*x(:,k)**2*x(:,l));s1112=sum(x(:,i)*x(:,j)*x(:,k)*x(:,l)**2)
              m1111=sum(x(:,i)*x(:,j)*x(:,k)*x(:,l))/rn;m0111=sum(x(:,j)*x(:,k)*x(:,l))/rn
              m1011=sum(x(:,i)*x(:,k)*x(:,l))/rn;m1101=sum(x(:,i)*x(:,j)*x(:,l))/rn
              m1110=sum(x(:,i)*x(:,j)*x(:,k))/rn
              terms(1)=terms(1)+24.0_dp*(s2222/rn-m1111**2-2.0_dp*s2111/rn*m0111- &
                2.0_dp*s1211/rn*m1011-2.0_dp*s1121/rn*m1101-2.0_dp*s1112/rn*m1110+ &
                2.0_dp*m11(i,l)*m1110*m0111+2.0_dp*m11(j,l)*m1011*m1110+ &
                2.0_dp*m11(k,l)*m1101*m1110+2.0_dp*m11(i,k)*m0111*m1101+ &
                2.0_dp*m11(j,k)*m1011*m1101+2.0_dp*m11(i,j)*m0111*m1011+ &
                m1110**2*m11(l,l)+m1101**2*m11(k,k)+m1011**2*m11(j,j)+m0111**2*m11(i,i))/rn
            end if
          end do
        end do
      end do
    end do
    do i=1,p
      do j=i+1,p
        s44=sum(x(:,i)**4*x(:,j)**4)
        terms(2)=terms(2)+2.0_dp*(s44/rn-m22(i,i)*m22(j,j)-4.0_dp*m41(i,j)*m21(j,j)- &
          4.0_dp*m41(j,i)*m21(i,i)+16.0_dp*m21(i,i)*m21(j,j)*m11(i,j))/(rn*rp)
      end do
    end do
  end subroutine exact_vm4_terms


  real(dp) function cm4_constant_correlation(x,m11,m21,m22,m31,m32,m33,m41,coef,marg6,marg7) result(value)
    real(dp),intent(in)::x(:,:),m11(:,:),m21(:,:),m22(:,:),m31(:,:),m32(:,:),m33(:,:),m41(:,:)
    real(dp),intent(in)::coef(4),marg6(:),marg7(:)
    real(dp)::rn,r3,r5,r6,r7,s8,s91,s62,s26,s611,s251,s215,s311,s221,s212,m211,m111
    real(dp)::s161,s521,s125,s131,s122,m121,s116,s152,s512,s113,m112
    real(dp)::s5111,s1511,s1151,s1115,s2111,s1211,s1121,s1112,m1111,m0111,m1011,m1101,m1110
    real(dp)::ti,tj,tk,tl,wi,wj,wk,wl,root5
    integer::n,p,i,j,k,l
    n=size(x,1);p=size(x,2);rn=real(n,dp);r3=coef(1);r5=coef(2);r6=coef(3);r7=coef(4)
    root5=sqrt(max(r5,0.0_dp));value=0.0_dp
    if(n<=0)return
    do i=1,p
      do j=i,p
        do k=j,p
          do l=k,p
            if(i==j)then
              if(j==k)then
                if(k==l)then
                  s8=sum(x(:,i)**8)
                  value=value+(s8/rn-m31(i,i)**2-8.0_dp*m41(i,i)*m21(i,i)+ &
                    16.0_dp*m11(i,i)*m21(i,i)**2)/rn
                else
                  s91=sum(x(:,i)**9*x(:,l));s62=sum(x(:,i)**6*x(:,l)**2)
                  ti=s91/rn-m31(i,l)*marg6(i)-3.0_dp*m21(i,l)*marg7(i)-m21(i,i)*s62/rn- &
                    6.0_dp*m41(i,i)*m41(i,l)+18.0_dp*m41(i,i)*m21(i,l)*m11(i,i)+ &
                    6.0_dp*m41(i,i)*m21(i,i)*m11(i,l)
                  tl=m33(i,l)-m31(i,l)*m11(l,l)-3.0_dp*m21(l,i)*m21(i,l)-m21(i,i)*m21(l,l)
                  value=value+2.0_dp*r3*(safe_sqrt_ratio(m11(l,l),marg6(i))*ti+ &
                    safe_sqrt_ratio(marg6(i),m11(l,l))*tl)/rn
                end if
              else if(k==l)then
                s62=sum(x(:,i)**6*x(:,k)**2);s26=sum(x(:,i)**2*x(:,k)**6)
                ti=s62/rn-m22(i,k)*m22(i,i)-4.0_dp*m32(i,k)*m21(i,i)-2.0_dp*m21(k,i)*m32(i,i)- &
                  2.0_dp*m21(i,k)*m41(i,k)+8.0_dp*m21(i,k)*m21(i,i)*m11(i,k)+ &
                  8.0_dp*m21(k,i)*m21(i,i)*m11(i,i)
                tk=s26/rn-m22(k,i)*m22(k,k)-4.0_dp*m32(k,i)*m21(k,k)-2.0_dp*m21(i,k)*m32(k,k)- &
                  2.0_dp*m21(k,i)*m41(k,i)+8.0_dp*m21(k,i)*m21(k,k)*m11(i,k)+ &
                  8.0_dp*m21(i,k)*m21(k,k)*m11(k,k)
                value=value+3.0_dp*r5*(safe_sqrt_ratio(m22(k,k),m22(i,i))*ti+ &
                  safe_sqrt_ratio(m22(i,i),m22(k,k))*tk)/rn
              else
                s611=sum(x(:,i)**6*x(:,k)*x(:,l));s251=sum(x(:,i)**2*x(:,k)**5*x(:,l))
                s215=sum(x(:,i)**2*x(:,k)*x(:,l)**5);s311=sum(x(:,i)**3*x(:,k)*x(:,l))
                s221=sum(x(:,i)**2*x(:,k)**2*x(:,l));s212=sum(x(:,i)**2*x(:,k)*x(:,l)**2)
                m211=sum(x(:,i)**2*x(:,k)*x(:,l))/rn;m111=sum(x(:,i)*x(:,k)*x(:,l))/rn
                ti=s611/rn-m211*m22(i,i)-4.0_dp*s311/rn*m21(i,i)-2.0_dp*m111*m32(i,i)- &
                  m21(i,l)*m41(i,k)-m21(i,k)*m41(i,l)+4.0_dp*m21(i,k)*m21(i,i)*m11(i,l)+ &
                  4.0_dp*m21(i,l)*m21(i,i)*m11(i,k)+8.0_dp*m111*m21(i,i)*m11(i,i)
                tk=s251/rn-m211*m22(k,k)-4.0_dp*s221/rn*m21(k,k)-m21(i,l)*m32(k,k)- &
                  2.0_dp*m111*m41(k,i)-m21(i,k)*m41(k,l)+4.0_dp*m21(i,k)*m21(k,k)*m11(k,l)+ &
                  8.0_dp*m111*m21(k,k)*m11(i,k)+4.0_dp*m21(i,l)*m11(k,k)*m21(k,k)
                tl=s215/rn-m211*m22(l,l)-4.0_dp*s212/rn*m21(l,l)-m21(i,k)*m32(l,l)- &
                  2.0_dp*m111*m41(l,i)-m21(i,l)*m41(l,k)+4.0_dp*m21(i,l)*m21(l,l)*m11(k,l)+ &
                  8.0_dp*m111*m21(l,l)*m11(i,l)+4.0_dp*m21(i,k)*m11(l,l)*m21(l,l)
                wi=safe_sqrt_ratio(sqrt(max(m22(k,k)*m22(l,l),0.0_dp)),m22(i,i))
                wk=0.0_dp;wl=0.0_dp
                if(m22(i,i)>=0.0_dp .and. m22(k,k)>tiny(1.0_dp) .and. m22(l,l)>=0.0_dp) &
                  wk=sqrt(max(m22(i,i)*sqrt(m22(l,l)/m22(k,k)**3),0.0_dp))
                if(m22(i,i)>=0.0_dp .and. m22(l,l)>tiny(1.0_dp) .and. m22(k,k)>=0.0_dp) &
                  wl=sqrt(max(m22(i,i)*sqrt(m22(k,k)/m22(l,l)**3),0.0_dp))
                value=value+3.0_dp*r6*root5*(2.0_dp*wi*ti+wk*tk+wl*tl)/rn
              end if
            else if(j==k)then
              if(k==l)then
                s91=sum(x(:,j)**9*x(:,i));s62=sum(x(:,j)**6*x(:,i)**2)
                tj=s91/rn-m31(j,i)*marg6(j)-3.0_dp*m21(j,i)*marg7(j)-m21(j,j)*s62/rn- &
                  6.0_dp*m41(j,j)*m41(j,i)+18.0_dp*m41(j,j)*m21(j,i)*m11(j,j)+ &
                  6.0_dp*m41(j,j)*m21(j,j)*m11(i,j)
                ti=m33(j,i)-m31(j,i)*m11(i,i)-3.0_dp*m21(i,j)*m21(j,i)-m21(j,j)*m21(i,i)
                value=value+2.0_dp*r3*(safe_sqrt_ratio(m11(i,i),marg6(j))*tj+ &
                  safe_sqrt_ratio(marg6(j),m11(i,i))*ti)/rn
              else
                s161=sum(x(:,i)*x(:,j)**6*x(:,l));s521=sum(x(:,i)**5*x(:,j)**2*x(:,l))
                s125=sum(x(:,i)*x(:,j)**2*x(:,l)**5);s131=sum(x(:,i)*x(:,j)**3*x(:,l))
                s221=sum(x(:,i)**2*x(:,j)**2*x(:,l));s122=sum(x(:,i)*x(:,j)**2*x(:,l)**2)
                m121=sum(x(:,i)*x(:,j)**2*x(:,l))/rn;m111=sum(x(:,i)*x(:,j)*x(:,l))/rn
                tj=s161/rn-m121*m22(j,j)-4.0_dp*s131/rn*m21(j,j)-2.0_dp*m111*m32(j,j)- &
                  m21(j,l)*m41(j,i)-m21(j,i)*m41(j,l)+4.0_dp*m21(j,i)*m21(j,j)*m11(j,l)+ &
                  4.0_dp*m21(j,l)*m21(j,j)*m11(i,j)+8.0_dp*m111*m21(j,j)*m11(j,j)
                ti=s521/rn-m121*m22(i,i)-4.0_dp*s221/rn*m21(i,i)-m21(j,l)*m32(i,i)- &
                  2.0_dp*m111*m41(i,j)-m21(j,i)*m41(i,l)+4.0_dp*m21(j,i)*m21(i,i)*m11(i,l)+ &
                  8.0_dp*m111*m21(i,i)*m11(i,j)+4.0_dp*m21(j,l)*m11(i,i)*m21(i,i)
                tl=s125/rn-m121*m22(l,l)-4.0_dp*s122/rn*m21(l,l)-m21(j,i)*m32(l,l)- &
                  2.0_dp*m111*m41(l,j)-m21(j,l)*m41(l,i)+4.0_dp*m21(j,l)*m21(l,l)*m11(i,l)+ &
                  8.0_dp*m111*m21(l,l)*m11(j,l)+4.0_dp*m21(j,i)*m11(l,l)*m21(l,l)
                wj=safe_sqrt_ratio(sqrt(max(m22(i,i)*m22(l,l),0.0_dp)),m22(j,j))
                wi=0.0_dp;wl=0.0_dp
                if(m22(j,j)>=0.0_dp .and. m22(i,i)>tiny(1.0_dp) .and. m22(l,l)>=0.0_dp) &
                  wi=sqrt(max(m22(j,j)*sqrt(m22(l,l)/m22(i,i)**3),0.0_dp))
                if(m22(j,j)>=0.0_dp .and. m22(l,l)>tiny(1.0_dp) .and. m22(i,i)>=0.0_dp) &
                  wl=sqrt(max(m22(j,j)*sqrt(m22(i,i)/m22(l,l)**3),0.0_dp))
                value=value+3.0_dp*r6*root5*(2.0_dp*wj*tj+wi*ti+wl*tl)/rn
              end if
            else if(k==l)then
              s116=sum(x(:,i)*x(:,j)*x(:,k)**6);s152=sum(x(:,i)*x(:,j)**5*x(:,k)**2)
              s512=sum(x(:,i)**5*x(:,j)*x(:,k)**2);s113=sum(x(:,i)*x(:,j)*x(:,k)**3)
              s122=sum(x(:,i)*x(:,j)**2*x(:,k)**2);s212=sum(x(:,i)**2*x(:,j)*x(:,k)**2)
              m112=sum(x(:,i)*x(:,j)*x(:,k)**2)/rn;m111=sum(x(:,i)*x(:,j)*x(:,k))/rn
              tk=s116/rn-m112*m22(k,k)-4.0_dp*s113/rn*m21(k,k)-2.0_dp*m111*m32(k,k)- &
                m21(k,j)*m41(k,i)-m21(k,i)*m41(k,j)+4.0_dp*m21(k,i)*m21(k,k)*m11(j,k)+ &
                4.0_dp*m21(k,j)*m21(k,k)*m11(i,k)+8.0_dp*m111*m21(k,k)*m11(k,k)
              ti=s512/rn-m112*m22(i,i)-4.0_dp*s212/rn*m21(i,i)-m21(k,j)*m32(i,i)- &
                2.0_dp*m111*m41(i,k)-m21(k,i)*m41(i,j)+4.0_dp*m21(k,i)*m21(i,i)*m11(i,j)+ &
                8.0_dp*m111*m21(i,i)*m11(i,k)+4.0_dp*m21(k,j)*m11(i,i)*m21(i,i)
              tj=s152/rn-m112*m22(j,j)-4.0_dp*s122/rn*m21(j,j)-m21(k,i)*m32(j,j)- &
                2.0_dp*m111*m41(j,k)-m21(k,j)*m41(j,i)+4.0_dp*m21(k,j)*m21(j,j)*m11(i,j)+ &
                8.0_dp*m111*m21(j,j)*m11(j,k)+4.0_dp*m21(k,i)*m11(j,j)*m21(j,j)
              wk=safe_sqrt_ratio(sqrt(max(m22(i,i)*m22(j,j),0.0_dp)),m22(k,k))
              wi=0.0_dp;wj=0.0_dp
              if(m22(k,k)>=0.0_dp .and. m22(i,i)>tiny(1.0_dp) .and. m22(j,j)>=0.0_dp) &
                wi=sqrt(max(m22(k,k)*sqrt(m22(j,j)/m22(i,i)**3),0.0_dp))
              if(m22(k,k)>=0.0_dp .and. m22(j,j)>tiny(1.0_dp) .and. m22(i,i)>=0.0_dp) &
                wj=sqrt(max(m22(k,k)*sqrt(m22(i,i)/m22(j,j)**3),0.0_dp))
              value=value+3.0_dp*r6*root5*(2.0_dp*wk*tk+wi*ti+wj*tj)/rn
            else
              s5111=sum(x(:,i)**5*x(:,j)*x(:,k)*x(:,l));s1511=sum(x(:,i)*x(:,j)**5*x(:,k)*x(:,l))
              s1151=sum(x(:,i)*x(:,j)*x(:,k)**5*x(:,l));s1115=sum(x(:,i)*x(:,j)*x(:,k)*x(:,l)**5)
              s2111=sum(x(:,i)**2*x(:,j)*x(:,k)*x(:,l));s1211=sum(x(:,i)*x(:,j)**2*x(:,k)*x(:,l))
              s1121=sum(x(:,i)*x(:,j)*x(:,k)**2*x(:,l));s1112=sum(x(:,i)*x(:,j)*x(:,k)*x(:,l)**2)
              m1111=sum(x(:,i)*x(:,j)*x(:,k)*x(:,l))/rn;m0111=sum(x(:,j)*x(:,k)*x(:,l))/rn
              m1011=sum(x(:,i)*x(:,k)*x(:,l))/rn;m1101=sum(x(:,i)*x(:,j)*x(:,l))/rn
              m1110=sum(x(:,i)*x(:,j)*x(:,k))/rn
              ! The C source divides m1111 by N a second time in these four terms.  The
              ! influence-function covariance requires m1111 itself; the corrected form is used.
              ti=s5111/rn-m1111*m22(i,i)-4.0_dp*s2111/rn*m21(i,i)-m0111*m41(i,i)- &
                m1011*m41(i,j)-m1101*m41(i,k)-m1110*m41(i,l)+4.0_dp*m1110*m21(i,i)*m11(i,l)+ &
                4.0_dp*m1101*m21(i,i)*m11(i,k)+4.0_dp*m1011*m21(i,i)*m11(i,j)+ &
                4.0_dp*m0111*m11(i,i)*m21(i,i)
              tj=s1511/rn-m1111*m22(j,j)-4.0_dp*s1211/rn*m21(j,j)-m0111*m41(j,i)- &
                m1011*m41(j,j)-m1101*m41(j,k)-m1110*m41(j,l)+4.0_dp*m1110*m21(j,j)*m11(j,l)+ &
                4.0_dp*m1101*m21(j,j)*m11(j,k)+4.0_dp*m1011*m21(j,j)*m11(j,j)+ &
                4.0_dp*m0111*m11(i,j)*m21(j,j)
              tk=s1151/rn-m1111*m22(k,k)-4.0_dp*s1121/rn*m21(k,k)-m0111*m41(k,i)- &
                m1011*m41(k,j)-m1101*m41(k,k)-m1110*m41(k,l)+4.0_dp*m1110*m21(k,k)*m11(k,l)+ &
                4.0_dp*m1101*m21(k,k)*m11(k,k)+4.0_dp*m1011*m21(k,k)*m11(j,k)+ &
                4.0_dp*m0111*m11(i,k)*m21(k,k)
              tl=s1115/rn-m1111*m22(l,l)-4.0_dp*s1112/rn*m21(l,l)-m0111*m41(l,i)- &
                m1011*m41(l,j)-m1101*m41(l,k)-m1110*m41(l,l)+4.0_dp*m1110*m21(l,l)*m11(l,l)+ &
                4.0_dp*m1101*m21(l,l)*m11(k,l)+4.0_dp*m1011*m21(l,l)*m11(j,l)+ &
                4.0_dp*m0111*m11(i,l)*m21(l,l)
              wi=safe_fourth_root_ratio(m22(j,j)*m22(k,k)*m22(l,l),m22(i,i)**3)
              wj=safe_fourth_root_ratio(m22(i,i)*m22(k,k)*m22(l,l),m22(j,j)**3)
              wk=safe_fourth_root_ratio(m22(i,i)*m22(j,j)*m22(l,l),m22(k,k)**3)
              wl=safe_fourth_root_ratio(m22(i,i)*m22(j,j)*m22(k,k),m22(l,l)**3)
              value=value+6.0_dp*r7*r5*(wi*ti+wj*tj+wk*tk+wl*tl)/rn
            end if
          end do
        end do
      end do
    end do
  end function cm4_constant_correlation

  pure real(dp) function m4_one_factor_directional(i,j,k,l,beta,fvar,fkurt,epsvar, &
      dmargkurt,dbeta,dfvar,dfkurt,depsvar) result(dt)
    integer,intent(in)::i,j,k,l
    real(dp),intent(in)::beta(:),fvar,fkurt,epsvar(:)
    real(dp),intent(in)::dmargkurt(:),dbeta(:),dfvar,dfkurt,depsvar(:)
    real(dp)::bi,bj,bk,bl,di,dj,dk,dl
    bi=beta(i);bj=beta(j);bk=beta(k);bl=beta(l)
    di=dbeta(i);dj=dbeta(j);dk=dbeta(k);dl=dbeta(l)
    if(i==j)then
      if(j==k)then
        if(k==l)then
          dt=dmargkurt(i)
        else
          dt=(3.0_dp*bi*bi*di*bl+bi**3*dl)*fkurt+bi**3*bl*dfkurt+ &
            3.0_dp*((di*bl+bi*dl)*fvar*epsvar(i)+bi*bl*dfvar*epsvar(i)+ &
            bi*bl*fvar*depsvar(i))
        end if
      else if(k==l)then
        dt=(2.0_dp*bi*di*bk*bk+2.0_dp*bi*bi*bk*dk)*fkurt+bi*bi*bk*bk*dfkurt+ &
          dfvar*(bi*bi*epsvar(k)+bk*bk*epsvar(i))+fvar*(2.0_dp*bi*di*epsvar(k)+ &
          bi*bi*depsvar(k)+2.0_dp*bk*dk*epsvar(i)+bk*bk*depsvar(i))+ &
          depsvar(i)*epsvar(k)+epsvar(i)*depsvar(k)
      else
        dt=(2.0_dp*bi*di*bk*bl+bi*bi*dk*bl+bi*bi*bk*dl)*fkurt+bi*bi*bk*bl*dfkurt+ &
          (dk*bl+bk*dl)*fvar*epsvar(i)+bk*bl*dfvar*epsvar(i)+bk*bl*fvar*depsvar(i)
      end if
    else if(j==k)then
      if(k==l)then
        dt=(di*bj**3+3.0_dp*bi*bj*bj*dj)*fkurt+bi*bj**3*dfkurt+ &
          3.0_dp*((di*bj+bi*dj)*fvar*epsvar(j)+bi*bj*dfvar*epsvar(j)+ &
          bi*bj*fvar*depsvar(j))
      else
        dt=(di*bj*bj*bl+2.0_dp*bi*bj*dj*bl+bi*bj*bj*dl)*fkurt+bi*bj*bj*bl*dfkurt+ &
          (di*bl+bi*dl)*fvar*epsvar(j)+bi*bl*dfvar*epsvar(j)+bi*bl*fvar*depsvar(j)
      end if
    else if(k==l)then
      dt=(di*bj*bk*bk+bi*dj*bk*bk+2.0_dp*bi*bj*bk*dk)*fkurt+bi*bj*bk*bk*dfkurt+ &
        (di*bj+bi*dj)*fvar*epsvar(k)+bi*bj*dfvar*epsvar(k)+bi*bj*fvar*depsvar(k)
    else
      dt=(di*bj*bk*bl+bi*dj*bk*bl+bi*bj*dk*bl+bi*bj*bk*dl)*fkurt+ &
        bi*bj*bk*bl*dfkurt
    end if
  end function m4_one_factor_directional

  real(dp) function cm4_one_factor(x,fc,fvar,fskew,fkurt,m11) result(value)
    real(dp),intent(in)::x(:,:),fc(:),fvar,fskew,fkurt,m11(:,:)
    real(dp),allocatable::covxf(:),beta(:),epsvar(:),m3diag(:),m4diag(:),m3full(:,:,:)
    real(dp),allocatable::dmarg(:),dbeta(:),deps(:)
    real(dp)::rn,dfvar,dfkurt,if_sample,if_target,mu4
    integer::n,p,i,j,k,l,t,mult
    n=size(x,1);p=size(x,2);rn=real(n,dp);value=0.0_dp
    if(n<=0 .or. fvar<=tiny(1.0_dp))return
    allocate(covxf(p),beta(p),epsvar(p),m3diag(p),m4diag(p),m3full(p,p,p))
    allocate(dmarg(p),dbeta(p),deps(p))
    do i=1,p
      covxf(i)=sum(x(:,i)*fc)/rn
      beta(i)=covxf(i)/fvar
      epsvar(i)=m11(i,i)-covxf(i)*covxf(i)/fvar
      m3diag(i)=sum(x(:,i)**3)/rn
      m4diag(i)=sum(x(:,i)**4)/rn
    end do
    do i=1,p;do j=1,p;do k=1,p
      m3full(i,j,k)=sum(x(:,i)*x(:,j)*x(:,k))/rn
    end do;end do;end do
    do t=1,n
      dfvar=fc(t)*fc(t)-fvar
      dfkurt=fc(t)**4-fkurt-4.0_dp*fc(t)*fskew
      do i=1,p
        dmarg(i)=x(t,i)**4-m4diag(i)-4.0_dp*x(t,i)*m3diag(i)
        dbeta(i)=(x(t,i)*fc(t)-covxf(i))/fvar-covxf(i)*dfvar/(fvar*fvar)
        deps(i)=x(t,i)**2-m11(i,i)-2.0_dp*beta(i)*fvar*dbeta(i)-beta(i)*beta(i)*dfvar
      end do
      do i=1,p;do j=i,p;do k=j,p;do l=k,p
        mu4=sum(x(:,i)*x(:,j)*x(:,k)*x(:,l))/rn
        if_sample=x(t,i)*x(t,j)*x(t,k)*x(t,l)-mu4- &
          x(t,i)*m3full(j,k,l)-x(t,j)*m3full(i,k,l)- &
          x(t,k)*m3full(i,j,l)-x(t,l)*m3full(i,j,k)
        if_target=m4_one_factor_directional(i,j,k,l,beta,fvar,fkurt,epsvar, &
          dmarg,dbeta,dfvar,dfkurt,deps)
        mult=m4_multiplicity(i,j,k,l)
        value=value+real(mult,dp)*if_sample*if_target/(rn*rn)
      end do;end do;end do;end do
    end do
  end function cm4_one_factor

  pure logical function selected_target(targets,code) result(selected)
    integer,intent(in)::targets(:),code
    selected=any(targets==code)
  end function selected_target

  integer function shrink_target_count(targets,max_code,factors,one_factor_code) result(count)
    integer,intent(in)::targets(:),max_code,one_factor_code
    real(dp),intent(in),optional::factors(:,:)
    integer::code
    count=0
    do code=1,max_code
      if(selected_target(targets,code))count=count+1
    end do
    if(selected_target(targets,one_factor_code) .and. present(factors)) &
      count=count+max(size(factors,2)-1,0)
  end function shrink_target_count

  subroutine finalize_shrinkage(sample_vec,targets,p,order,a,b,result)
    real(dp),intent(in)::sample_vec(:),targets(:,:),a(:,:),b(:)
    integer,intent(in)::p,order
    type(exact_shrinkage_result),intent(out)::result
    real(dp),allocatable::estimate_vec(:)
    integer::j
    result%order=order;result%n_targets=size(targets,2)
    allocate(result%lambda(size(b)),result%a(size(a,1),size(a,2)),result%b(size(b)))
    allocate(result%target_vectors(size(targets,1),size(targets,2)))
    result%a=a;result%b=b;result%target_vectors=targets
    call solve_shrinkage_qp(a,b,result%lambda,result%converged,result%objective)
    allocate(estimate_vec(size(sample_vec)));estimate_vec=(1.0_dp-sum(result%lambda))*sample_vec
    do j=1,size(targets,2);estimate_vec=estimate_vec+result%lambda(j)*targets(:,j);end do
    select case(order)
    case(2)
      allocate(result%sample(p,p),result%estimate(p,p))
      result%sample=reshape(sample_vec,[p,p]);result%estimate=reshape(estimate_vec,[p,p])
    case(3)
      allocate(result%sample(p,p*p),result%estimate(p,p*p))
      call m3_vec_to_mat(sample_vec,p,result%sample);call m3_vec_to_mat(estimate_vec,p,result%estimate)
    case(4)
      allocate(result%sample(p,p*p*p),result%estimate(p,p*p*p))
      call m4_vec_to_mat(sample_vec,p,result%sample);call m4_vec_to_mat(estimate_vec,p,result%estimate)
    end select
  end subroutine finalize_shrinkage

  subroutine exact_m2_shrinkage(r,targets,result,factors)
    real(dp),intent(in)::r(:,:)
    integer,intent(in)::targets(:)
    type(exact_shrinkage_result),intent(out)::result
    real(dp),intent(in),optional::factors(:,:)
    real(dp),allocatable::x(:,:),m11(:,:),m22(:,:),sample(:,:),sample_vec(:),tv(:,:),tmat(:,:)
    real(dp),allocatable::margvar(:),fc(:),beta(:)
    real(dp),allocatable::a(:,:),b(:)
    real(dp)::terms(3),fvar,rcoef
    integer::n,p,nt,code,col,jf,i,j
    if(size(r,2)<2)error stop 'exact_m2_shrinkage: at least two variables are required'
    if(size(targets)==0 .or. any(targets<1) .or. any(targets>4))error stop 'exact_m2_shrinkage: invalid target'
    if(selected_target(targets,3) .and. .not.present(factors))error stop 'exact_m2_shrinkage: factors required'
    n=size(r,1);p=size(r,2);rcoef=0.0_dp
    if(present(factors))then
      if(size(factors,1)/=n)error stop 'exact_m2_shrinkage: factor length mismatch'
    end if
    nt=shrink_target_count(targets,4,factors,3)
    call center_matrix(r,x)
    allocate(m11(p,p),m22(p,p),sample(p,p),margvar(p),sample_vec(p*p),tv(p*p,nt),tmat(p,p))
    call pair_moment(x,1,1,m11);call pair_moment(x,2,2,m22)
    sample=m11;margvar=[(m11(i,i),i=1,p)];sample_vec=reshape(sample,[p*p]);col=0
    do code=1,4
      if(.not.selected_target(targets,code))cycle
      col=col+1;tmat=0.0_dp
      select case(code)
      case(1)
        do i=1,p;tmat(i,i)=margvar(i);end do
      case(2)
        do i=1,p;tmat(i,i)=sum(margvar)/real(p,dp);end do
      case(3)
        allocate(fc(n),beta(p));fc=factors(:,1)-mean_value(factors(:,1));fvar=sum(fc*fc)/real(n,dp)
        do i=1,p;beta(i)=sum(x(:,i)*fc)/real(n,dp)/fvar;end do
        do i=1,p;do j=1,p;tmat(i,j)=fvar*beta(i)*beta(j);end do;end do
        do i=1,p;tmat(i,i)=margvar(i);end do
        deallocate(fc,beta)
      case(4)
        rcoef=0.0_dp
        do i=1,p;do j=i+1,p
          if(margvar(i)>0.0_dp .and. margvar(j)>0.0_dp)rcoef=rcoef+m11(i,j)/sqrt(margvar(i)*margvar(j))
        end do;end do
        rcoef=2.0_dp*rcoef/real(p*(p-1),dp)
        do i=1,p;do j=1,p
          if(i==j)then;tmat(i,j)=margvar(i)
          else;tmat(i,j)=rcoef*sqrt(max(margvar(i)*margvar(j),0.0_dp));end if
        end do;end do
      end select
      tv(:,col)=reshape(tmat,[p*p])
    end do
    if(selected_target(targets,3))then
      if(size(factors,2)>1)then
      do jf=2,size(factors,2)
        col=col+1;allocate(fc(n),beta(p));fc=factors(:,jf)-mean_value(factors(:,jf));fvar=sum(fc*fc)/real(n,dp)
        do i=1,p;beta(i)=sum(x(:,i)*fc)/real(n,dp)/fvar;end do
        do i=1,p;do j=1,p;tmat(i,j)=fvar*beta(i)*beta(j);end do;end do
        do i=1,p;tmat(i,i)=margvar(i);end do
        tv(:,col)=reshape(tmat,[p*p]);deallocate(fc,beta)
      end do
      end if
    end if
    allocate(a(nt,nt),b(nt));do i=1,nt;do j=i,nt
      a(i,j)=dot_product(tv(:,i)-sample_vec,tv(:,j)-sample_vec);a(j,i)=a(i,j)
    end do;end do
    call exact_vm2_terms(m11,m22,n,terms);b=terms(1);col=0
    do code=1,4
      if(.not.selected_target(targets,code))cycle
      col=col+1
      select case(code)
      case(1);b(col)=b(col)-terms(3)
      case(2);b(col)=b(col)-terms(2)
      case(3)
        allocate(fc(n));fc=factors(:,1)-mean_value(factors(:,1));fvar=sum(fc*fc)/real(n,dp)
        b(col)=b(col)-cm2_one_factor(x,fc,fvar,m11,m22);deallocate(fc)
      case(4);b(col)=b(col)-cm2_constant_correlation(x,rcoef,m11,m22)
      end select
    end do
    if(selected_target(targets,3))then
      if(size(factors,2)>1)then
      do jf=2,size(factors,2)
        col=col+1;allocate(fc(n));fc=factors(:,jf)-mean_value(factors(:,jf));fvar=sum(fc*fc)/real(n,dp)
        b(col)=b(col)-cm2_one_factor(x,fc,fvar,m11,m22);deallocate(fc)
      end do
      end if
    end if
    call finalize_shrinkage(sample_vec,tv,p,2,a,b,result)
  end subroutine exact_m2_shrinkage

  subroutine exact_m3_shrinkage(r,targets,result,factors,unbiased_mse)
    real(dp),intent(in)::r(:,:)
    integer,intent(in)::targets(:)
    type(exact_shrinkage_result),intent(out)::result
    real(dp),intent(in),optional::factors(:,:)
    logical,intent(in),optional::unbiased_mse
    real(dp),allocatable::x(:,:),m11(:,:),m21(:,:),m22(:,:),m31(:,:),m42(:,:),m33(:,:)
    real(dp),allocatable::m41(:,:),m61(:,:),m32(:,:),m51(:,:),sample_vec(:),tv(:,:),a(:,:),b(:)
    real(dp),allocatable::margvar(:),margskew(:),margkurt(:),marg5(:),marg6(:),root(:),fc(:),beta(:)
    real(dp)::terms(3),coef(3),fvar,fskew,scale
    integer::n,p,nu,nt,code,col,jf,i,j
    logical::unbiased
    if(size(r,2)<2)error stop 'exact_m3_shrinkage: at least two variables are required'
    if(size(targets)==0 .or. any(targets<1) .or. any(targets>6))error stop 'exact_m3_shrinkage: invalid target'
    unbiased=.false.;if(present(unbiased_mse))unbiased=unbiased_mse
    if(unbiased .and. any([(selected_target(targets,code),code=3,5)])) &
      error stop 'exact_m3_shrinkage: unbiased MSE supports targets 1, 2, and 6 only'
    if(selected_target(targets,3) .and. .not.present(factors))error stop 'exact_m3_shrinkage: factors required'
    n=size(r,1);p=size(r,2);if(unbiased .and. n<6)error stop 'exact_m3_shrinkage: six observations required'
    if(present(factors))then;if(size(factors,1)/=n)error stop 'exact_m3_shrinkage: factor length mismatch';end if
    nu=n_unique_m3(p);nt=shrink_target_count(targets,6,factors,3)
    call center_matrix(r,x)
    allocate(m11(p,p),m21(p,p),m22(p,p),m31(p,p),m42(p,p),m33(p,p),sample_vec(nu),tv(nu,nt))
    call pair_moment(x,1,1,m11);call pair_moment(x,2,1,m21);call pair_moment(x,2,2,m22)
    call pair_moment(x,3,1,m31);call pair_moment(x,4,2,m42);call pair_moment(x,3,3,m33)
    call coskewness_unique(r,sample_vec,unbiased)
    allocate(margvar(p),margskew(p),margkurt(p),marg5(p),marg6(p),root(p))
    do i=1,p
      margvar(i)=sum(x(:,i)**2)/real(n,dp);margskew(i)=sum(x(:,i)**3)/real(n,dp)
      margkurt(i)=sum(x(:,i)**4)/real(n,dp);marg5(i)=sum(x(:,i)**5)/real(n,dp)
      marg6(i)=sum(x(:,i)**6)/real(n,dp)
    end do
    if(unbiased)then;scale=real(n*n,dp)/real((n-1)*(n-2),dp);margskew=margskew*scale;end if
    col=0
    do code=1,6
      if(.not.selected_target(targets,code))cycle
      col=col+1
      select case(code)
      case(1);call m3_target_independent(margskew,tv(:,col))
      case(2);call m3_target_independent(spread(sum(margskew)/real(p,dp),1,p),tv(:,col))
      case(3)
        allocate(fc(n),beta(p));fc=factors(:,1)-mean_value(factors(:,1));fvar=sum(fc**2)/real(n,dp)
        fskew=sum(fc**3)/real(n,dp);do i=1,p;beta(i)=sum(x(:,i)*fc)/real(n,dp)/fvar;end do
        call m3_target_one_factor(margskew,beta,fskew,tv(:,col));deallocate(fc,beta)
      case(4);call m3_cc_coefficients(x,margvar,margkurt,m21,m22,coef);call m3_target_cc(margvar,margskew,margkurt,coef,tv(:,col))
      case(5)
        do i=1,p;root(i)=sign(abs(margskew(i))**(1.0_dp/3.0_dp),margskew(i));end do
        call m3_target_simaan(root,tv(:,col))
      case(6);tv(:,col)=0.0_dp
      end select
    end do
    if(selected_target(targets,3))then
      if(size(factors,2)>1)then
      do jf=2,size(factors,2)
        col=col+1;allocate(fc(n),beta(p));fc=factors(:,jf)-mean_value(factors(:,jf));fvar=sum(fc**2)/real(n,dp)
        fskew=sum(fc**3)/real(n,dp);do i=1,p;beta(i)=sum(x(:,i)*fc)/real(n,dp)/fvar;end do
        call m3_target_one_factor(margskew,beta,fskew,tv(:,col));deallocate(fc,beta)
      end do
      end if
    end if
    allocate(a(nt,nt),b(nt));do i=1,nt;do j=i,nt
      a(i,j)=m3_inner_product(tv(:,i)-sample_vec,tv(:,j)-sample_vec,p);a(j,i)=a(i,j)
    end do;end do
    if(unbiased)then
      call exact_vm3_kstat_terms(x,real(n,dp)*m11,real(n,dp)*m21,real(n,dp)*m22, &
        real(n,dp)*m31,real(n,dp)*m42,real(n,dp)*m33,terms)
    else;call exact_vm3_terms(x,m11,m21,m22,m31,m42,m33,terms);end if
    b=terms(1);col=0
    allocate(m41(p,p),m61(p,p),m32(p,p),m51(p,p));call pair_moment(x,4,1,m41)
    call pair_moment(x,6,1,m61);call pair_moment(x,3,2,m32);call pair_moment(x,5,1,m51)
    do code=1,6
      if(.not.selected_target(targets,code))cycle
      col=col+1
      select case(code)
      case(1);b(col)=b(col)-terms(3)
      case(2);b(col)=b(col)-terms(2)
      case(3)
        allocate(fc(n));fc=factors(:,1)-mean_value(factors(:,1));fvar=sum(fc**2)/real(n,dp);fskew=sum(fc**3)/real(n,dp)
        b(col)=b(col)-cm3_one_factor(x,fc,fvar,fskew,m11,m21,m22,m42);deallocate(fc)
      case(4)
        b(col)=b(col)-cm3_constant_correlation(x,margvar,margskew,margkurt,marg5,marg6, &
          m11,m21,m31,m32,m41,m61,coef)
      case(5)
        do i=1,p
          if(abs(margskew(i))>tiny(1.0_dp))then;root(i)=abs(margskew(i))**(-2.0_dp/3.0_dp)
          else;root(i)=0.0_dp;end if
        end do
        b(col)=b(col)-cm3_simaan(x,root,m11,m21,m22,m31,m42,m51)
      case(6)
      end select
    end do
    if(selected_target(targets,3))then
      if(size(factors,2)>1)then
      do jf=2,size(factors,2)
        col=col+1;allocate(fc(n));fc=factors(:,jf)-mean_value(factors(:,jf));fvar=sum(fc**2)/real(n,dp)
        fskew=sum(fc**3)/real(n,dp);b(col)=b(col)-cm3_one_factor(x,fc,fvar,fskew,m11,m21,m22,m42);deallocate(fc)
      end do
      end if
    end if
    call finalize_shrinkage(sample_vec,tv,p,3,a,b,result);result%unbiased_mse=unbiased
  end subroutine exact_m3_shrinkage

  subroutine exact_m4_shrinkage(r,targets,result,factors)
    real(dp),intent(in)::r(:,:)
    integer,intent(in)::targets(:)
    type(exact_shrinkage_result),intent(out)::result
    real(dp),intent(in),optional::factors(:,:)
    real(dp),allocatable::x(:,:),m11(:,:),m21(:,:),m22(:,:),m31(:,:),m32(:,:),m33(:,:),m41(:,:),m42(:,:)
    real(dp),allocatable::sample_vec(:),tv(:,:),a(:,:),b(:),margvar(:),margkurt(:),marg6(:),marg7(:)
    real(dp),allocatable::fc(:),beta(:),epsvar(:),pair(:)
    real(dp)::terms(3),coef(4),fvar,fskew,fkurt
    integer::n,p,nu,nt,code,col,jf,i,j
    if(size(r,2)<2)error stop 'exact_m4_shrinkage: at least two variables are required'
    if(size(targets)==0 .or. any(targets<1) .or. any(targets>4))error stop 'exact_m4_shrinkage: invalid target'
    if(selected_target(targets,3) .and. .not.present(factors))error stop 'exact_m4_shrinkage: factors required'
    n=size(r,1);p=size(r,2);if(present(factors))then
      if(size(factors,1)/=n)error stop 'exact_m4_shrinkage: factor length mismatch'
    end if
    nu=n_unique_m4(p);nt=shrink_target_count(targets,4,factors,3);call center_matrix(r,x)
    allocate(m11(p,p),m21(p,p),m22(p,p),m31(p,p),m32(p,p),m33(p,p),m41(p,p),m42(p,p))
    call pair_moment(x,1,1,m11);call pair_moment(x,2,1,m21);call pair_moment(x,2,2,m22)
    call pair_moment(x,3,1,m31);call pair_moment(x,3,2,m32);call pair_moment(x,3,3,m33)
    call pair_moment(x,4,1,m41);call pair_moment(x,4,2,m42)
    allocate(sample_vec(nu),tv(nu,nt),margvar(p),margkurt(p),marg6(p),marg7(p),pair(p))
    call cokurtosis_unique(r,sample_vec)
    do i=1,p
      margvar(i)=sum(x(:,i)**2)/real(n,dp);margkurt(i)=sum(x(:,i)**4)/real(n,dp)
      marg6(i)=sum(x(:,i)**6)/real(n,dp);marg7(i)=sum(x(:,i)**7)/real(n,dp)
    end do
    call m4_cc_coefficients(x,margvar,margkurt,marg6,m22,m31,coef);col=0
    do code=1,4
      if(.not.selected_target(targets,code))cycle
      col=col+1
      select case(code)
      case(1);call m4_target_independent(margkurt,margvar,tv(:,col))
      case(2)
        pair=sqrt(sum(margvar*margvar)/real(p,dp));call m4_target_independent( &
          spread(sum(margkurt)/real(p,dp),1,p),pair,tv(:,col))
      case(3)
        allocate(fc(n),beta(p),epsvar(p));fc=factors(:,1)-mean_value(factors(:,1));fvar=sum(fc**2)/real(n,dp)
        fkurt=sum(fc**4)/real(n,dp);do i=1,p
          beta(i)=sum(x(:,i)*fc)/real(n,dp)/fvar;epsvar(i)=margvar(i)-beta(i)*beta(i)*fvar
        end do
        call m4_target_one_factor(margkurt,fvar,fkurt,epsvar,beta,tv(:,col));deallocate(fc,beta,epsvar)
      case(4);call m4_target_cc(margvar,margkurt,marg6,coef,tv(:,col))
      end select
    end do
    if(selected_target(targets,3))then
      if(size(factors,2)>1)then
      do jf=2,size(factors,2)
        col=col+1;allocate(fc(n),beta(p),epsvar(p));fc=factors(:,jf)-mean_value(factors(:,jf));fvar=sum(fc**2)/real(n,dp)
        fkurt=sum(fc**4)/real(n,dp);do i=1,p
          beta(i)=sum(x(:,i)*fc)/real(n,dp)/fvar;epsvar(i)=margvar(i)-beta(i)*beta(i)*fvar
        end do
        call m4_target_one_factor(margkurt,fvar,fkurt,epsvar,beta,tv(:,col));deallocate(fc,beta,epsvar)
      end do
      end if
    end if
    allocate(a(nt,nt),b(nt));do i=1,nt;do j=i,nt
      a(i,j)=m4_inner_product(tv(:,i)-sample_vec,tv(:,j)-sample_vec,p);a(j,i)=a(i,j)
    end do;end do
    call exact_vm4_terms(x,m11,m21,m22,m31,m32,m41,m42,terms);b=terms(1);col=0
    do code=1,4
      if(.not.selected_target(targets,code))cycle
      col=col+1
      select case(code)
      case(1);b(col)=b(col)-terms(3)
      case(2);b(col)=b(col)-terms(2)
      case(3)
        allocate(fc(n));fc=factors(:,1)-mean_value(factors(:,1));fvar=sum(fc**2)/real(n,dp)
        fskew=sum(fc**3)/real(n,dp);fkurt=sum(fc**4)/real(n,dp)
        b(col)=b(col)-cm4_one_factor(x,fc,fvar,fskew,fkurt,m11);deallocate(fc)
      case(4);b(col)=b(col)-cm4_constant_correlation(x,m11,m21,m22,m31,m32,m33,m41,coef,marg6,marg7)
      end select
    end do
    if(selected_target(targets,3))then
      if(size(factors,2)>1)then
      do jf=2,size(factors,2)
        col=col+1;allocate(fc(n));fc=factors(:,jf)-mean_value(factors(:,jf));fvar=sum(fc**2)/real(n,dp)
        fskew=sum(fc**3)/real(n,dp);fkurt=sum(fc**4)/real(n,dp)
        b(col)=b(col)-cm4_one_factor(x,fc,fvar,fskew,fkurt,m11);deallocate(fc)
      end do
      end if
    end if
    call finalize_shrinkage(sample_vec,tv,p,4,a,b,result)
  end subroutine exact_m4_shrinkage

end module finite_sample_moments_mod
