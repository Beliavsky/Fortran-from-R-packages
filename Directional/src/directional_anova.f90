module directional_anova
   use directional_kinds, only : dp, pi
   use directional_inference, only : vmf_mle, vmf_mle_result
   use directional_special, only : log_bessel_i
   use directional_tests, only : test_result
   implicit none
   private
   public :: embed_aov, hcf_aov, het_aov, lr_aov
contains
   function embed_aov(x,group,mc_reps) result(res)
      real(dp),intent(in)::x(:,:);integer,intent(in)::group(:);integer,intent(in),optional::mc_reps
      type(test_result)::res;integer::lab(size(group)),b,i,ex
      res%statistic=embed_stat_groups(x,group);res%df=(maxval(group)-1)*(size(x,2)-1)
      b=0;if(present(mc_reps))b=mc_reps
      if(b>0)then;lab=group;ex=0;do i=1,b;call shuffle(lab);if(embed_stat_groups(x,lab)>res%statistic)ex=ex+1;end do;res%p_value=real(ex+1,dp)/real(b+1,dp);end if
   end function
   function hcf_aov(x,group,mc_reps,correction) result(res)
      real(dp),intent(in)::x(:,:);integer,intent(in)::group(:);integer,intent(in),optional::mc_reps;logical,intent(in),optional::correction
      type(test_result)::res;integer::lab(size(group)),b,i,ex;logical::fc
      fc=.true.;if(present(correction))fc=correction;res%statistic=hcf_stat_groups(x,group,fc);res%df=(maxval(group)-1)*(size(x,2)-1)
      b=0;if(present(mc_reps))b=mc_reps
      if(b>0)then;lab=group;ex=0;do i=1,b;call shuffle(lab);if(hcf_stat_groups(x,lab,fc)>res%statistic)ex=ex+1;end do;res%p_value=real(ex+1,dp)/real(b+1,dp);end if
   end function
   function het_aov(x,group,mc_reps) result(res)
      real(dp),intent(in)::x(:,:);integer,intent(in)::group(:);integer,intent(in),optional::mc_reps
      type(test_result)::res;integer::lab(size(group)),b,i,ex
      res%statistic=het_stat_groups(x,group);res%df=(maxval(group)-1)*(size(x,2)-1)
      b=0;if(present(mc_reps))b=mc_reps
      if(b>0)then;lab=group;ex=0;do i=1,b;call shuffle(lab);if(het_stat_groups(x,lab)>res%statistic)ex=ex+1;end do;res%p_value=real(ex+1,dp)/real(b+1,dp);end if
   end function
   function lr_aov(x,group,mc_reps) result(res)
      real(dp),intent(in)::x(:,:);integer,intent(in)::group(:);integer,intent(in),optional::mc_reps
      type(test_result)::res;integer::lab(size(group)),b,i,ex
      res%statistic=lr_stat_groups(x,group);res%df=(maxval(group)-1)*(size(x,2)-1)
      b=0;if(present(mc_reps))b=mc_reps
      if(b>0)then;lab=group;ex=0;do i=1,b;call shuffle(lab);if(lr_stat_groups(x,lab)>res%statistic)ex=ex+1;end do;res%p_value=real(ex+1,dp)/real(b+1,dp);end if
   end function
   real(dp) function embed_stat_groups(x,group) result(v)
      real(dp),intent(in)::x(:,:);integer,intent(in)::group(:);integer::g,n,p,j,i,ng;real(dp)::s(size(x,2)),r,within
      g=maxval(group);n=size(x,1);p=size(x,2);within=0.0_dp
      do j=1,g;s=0;ng=0;do i=1,n;if(group(i)==j)then;s=s+x(i,:);ng=ng+1;end if;end do;if(ng>0)within=within+sum(s*s)/real(ng,dp);end do
      s=sum(x,dim=1);r=sum(s*s)/real(n,dp);v=real(n-g,dp)*(within-r)/max(tiny(1.0_dp),real(n,dp)-within)
   end function
   real(dp) function hcf_stat_groups(x,group,fc) result(v)
      real(dp),intent(in)::x(:,:);integer,intent(in)::group(:);logical,intent(in)::fc
      integer::g,n,p,j,i;real(dp)::s(size(x,2)),ri,r,sumri,k,corr;type(vmf_mle_result)::fit
      g=maxval(group);n=size(x,1);p=size(x,2);sumri=0
      do j=1,g;s=0;do i=1,n;if(group(i)==j)s=s+x(i,:);end do;ri=sqrt(sum(s*s));sumri=sumri+ri;end do
      s=sum(x,dim=1);r=sqrt(sum(s*s));v=real(n-g,dp)*(sumri-r)/max(tiny(1.0_dp),real(g-1,dp)*(real(n,dp)-sumri))
      if(fc)then;fit=vmf_mle(x);k=max(fit%kappa,1.0e-8_dp);corr=1.0_dp;if(p==3)corr=k*(1.0_dp/k-1.0_dp/(5.0_dp*k**3));if(p>3)corr=k*(1.0_dp/k-real(p-3,dp)/(4*k*k)-real(p-3,dp)/(4*k**3));v=v*corr;end if
   end function
   real(dp) function het_stat_groups(x,group) result(v)
      real(dp),intent(in)::x(:,:);integer,intent(in)::group(:);integer::g,n,j,i,ng;real(dp)::s(size(x,2)),tw(size(x,2)),r,k;real(dp),allocatable::z(:,:);type(vmf_mle_result)::fit
      g=maxval(group);n=size(x,1);tw=0;v=0
      do j=1,g;ng=count(group==j);allocate(z(ng,size(x,2)));ng=0;s=0;do i=1,n;if(group(i)==j)then;ng=ng+1;z(ng,:)=x(i,:);s=s+x(i,:);end if;end do;fit=vmf_mle(z);k=fit%kappa;r=sqrt(sum(s*s));v=v+k*r;tw=tw+k*s;deallocate(z);end do
      v=2.0_dp*(v-sqrt(sum(tw*tw)))
   end function
   real(dp) function lr_stat_groups(x,group) result(v)
      real(dp),intent(in)::x(:,:);integer,intent(in)::group(:);integer::g,n,p,j,i;real(dp)::s(size(x,2)),r,sumri,k0,k1,nu,l0,l1
      g=maxval(group);n=size(x,1);p=size(x,2);sumri=0
      do j=1,g;s=0;do i=1,n;if(group(i)==j)s=s+x(i,:);end do;sumri=sumri+sqrt(sum(s*s));end do
      s=sum(x,dim=1);r=sqrt(sum(s*s));k0=kappa_from_rbar(r/real(n,dp),p);k1=kappa_from_rbar(sumri/real(n,dp),p);nu=0.5_dp*p-1.0_dp
      l0=real(n,dp)*(nu*log(max(k0,tiny(1.0_dp)))-0.5_dp*p*log(2*pi)-log_bessel_i(nu,k0))+k0*r
      l1=real(n,dp)*(nu*log(max(k1,tiny(1.0_dp)))-0.5_dp*p*log(2*pi)-log_bessel_i(nu,k1))+k1*sumri
      v=max(0.0_dp,2.0_dp*(l1-l0))
   end function
   real(dp) function kappa_from_rbar(rbar,p) result(k)
      real(dp),intent(in)::rbar;integer,intent(in)::p;real(dp)::a,kn;integer::it
      if(rbar<1.0e-10_dp)then;k=0;return;end if;k=max(1.0e-6_dp,rbar*(real(p,dp)-rbar*rbar)/max(1.0e-12_dp,1-rbar*rbar))
      do it=1,50;a=exp(log_bessel_i(0.5_dp*p,k)-log_bessel_i(0.5_dp*p-1.0_dp,k));kn=k-(a-rbar)/max(1.0e-10_dp,1-a*a-real(p-1,dp)*a/k);kn=max(1.0e-8_dp,kn);if(abs(kn-k)<1e-10_dp*max(1.0_dp,k))exit;k=kn;end do;k=kn
   end function
   subroutine shuffle(a)
      integer,intent(inout)::a(:);integer::i,j,t;real(dp)::u
      do i=size(a),2,-1;call random_number(u);j=min(i,1+int(u*i));t=a(i);a(i)=a(j);a(j)=t;end do
   end subroutine
end module directional_anova
