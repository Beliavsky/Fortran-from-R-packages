module directional_mixtures
   use directional_kinds, only : dp, pi
   use directional_special, only : log_bessel_i
   use directional_distributions, only : dspcauchy, dpkbd
   use directional_random, only : rspcauchy, rpkbd
   implicit none
   private
   public :: mixvmf_mle, mixspcauchy_mle, mixpkbd_mle, dmixspcauchy, dmixpkbd
   public :: rmixspcauchy, rmixpkbd
contains
   subroutine mixvmf_mle(x,g,probs,kappa,mu,loglik,cluster,tol,maxiters)
      real(dp),intent(in)::x(:,:);integer,intent(in)::g
      real(dp),intent(out)::probs(g),kappa(g),mu(g,size(x,2)),loglik
      integer,intent(out),optional::cluster(size(x,1));real(dp),intent(in),optional::tol;integer,intent(in),optional::maxiters
      real(dp)::w(size(x,1),g),ld(size(x,1),g),m(size(x,2)),rs(size(x,1)),r,old,eps,nu,mx
      integer::i,j,it,nit,p,n,cl
      n=size(x,1);p=size(x,2);nu=0.5_dp*p-1;eps=1e-6_dp;if(present(tol))eps=tol;nit=100;if(present(maxiters))nit=maxiters
      ! deterministic farthest-point initialization
      mu=0;mu(1,:)=x(1,:);do j=2,g
         cl=1;mx=-1
         do i=1,n;r=minval(1.0_dp-matmul(mu(1:j-1,:),x(i,:)));if(r>mx)then;mx=r;cl=i;end if;end do;mu(j,:)=x(cl,:)
      end do
      w=0
      do i=1,n;j=maxloc(matmul(mu,x(i,:)),dim=1);w(i,j)=1;end do
      old=-huge(1.0_dp)
      do it=1,nit
         do j=1,g
            probs(j)=max(sum(w(:,j))/n,1e-12_dp);m=matmul(transpose(x),w(:,j));r=sqrt(sum(m*m))/max(sum(w(:,j)),tiny(1.0_dp));mu(j,:)=m/max(sqrt(sum(m*m)),tiny(1.0_dp));kappa(j)=solve_kappa(r,p)
            ld(:,j)=log(probs(j))+nu*log(max(kappa(j),tiny(1.0_dp)))-0.5_dp*p*log(2*pi)-log_bessel_i(nu,max(kappa(j),tiny(1.0_dp)))+kappa(j)*matmul(x,mu(j,:))
         end do
         loglik=0
         do i=1,n;mx=maxval(ld(i,:));rs(i)=sum(exp(ld(i,:)-mx));loglik=loglik+mx+log(rs(i));w(i,:)=exp(ld(i,:)-mx)/rs(i);end do
         if(abs(loglik-old)<eps)exit;old=loglik
      end do
      probs=sum(w,dim=1)/n
      if(present(cluster))then;do i=1,n;cluster(i)=maxloc(w(i,:),dim=1);end do;end if
   end subroutine
   subroutine rmixspcauchy(n,probs,mu,rho,x,id)
      integer,intent(in)::n;real(dp),intent(in)::probs(:),mu(:,:),rho(:);real(dp),intent(out)::x(n,size(mu,2));integer,intent(out),optional::id(n)
      real(dp)::u,cum,tmp(1,size(mu,2));integer::i,j,g
      g=size(probs)
      do i=1,n
         call random_number(u);cum=0.0_dp
         do j=1,g;cum=cum+probs(j);if(u<=cum)exit;end do
         tmp=rspcauchy(1,mu(j,:),rho(j));x(i,:)=tmp(1,:);if(present(id))id(i)=j
      end do
   end subroutine

   subroutine rmixpkbd(n,probs,mu,rho,x,id)
      integer,intent(in)::n;real(dp),intent(in)::probs(:),mu(:,:),rho(:);real(dp),intent(out)::x(n,size(mu,2));integer,intent(out),optional::id(n)
      real(dp)::u,cum,tmp(1,size(mu,2));integer::i,j,g
      g=size(probs)
      do i=1,n
         call random_number(u);cum=0.0_dp
         do j=1,g;cum=cum+probs(j);if(u<=cum)exit;end do
         if(j>g)j=g;tmp=rpkbd(1,mu(j,:),rho(j));x(i,:)=tmp(1,:);if(present(id))id(i)=j
      end do
   end subroutine

   function dmixspcauchy(y,probs,mu,rho,logden) result(v)
      real(dp),intent(in)::y(:,:),probs(:),mu(:,:),rho(:);logical,intent(in),optional::logden
      real(dp)::v(size(y,1)),tmp(size(y,1)),s;integer::i,j;logical::ll
      ll=.false.;if(present(logden))ll=logden;v=0.0_dp
      do j=1,size(probs);tmp=dspcauchy(y,mu(j,:),rho(j));v=v+probs(j)*tmp;end do
      if(ll)then;do i=1,size(v);v(i)=log(max(v(i),tiny(1.0_dp)));end do;end if
   end function

   function dmixpkbd(y,probs,mu,rho,logden) result(v)
      real(dp),intent(in)::y(:,:),probs(:),mu(:,:),rho(:);logical,intent(in),optional::logden
      real(dp)::v(size(y,1)),tmp(size(y,1));integer::i,j;logical::ll
      ll=.false.;if(present(logden))ll=logden;v=0.0_dp
      do j=1,size(probs);tmp=dpkbd(y,mu(j,:),rho(j));v=v+probs(j)*tmp;end do
      if(ll)then;do i=1,size(v);v(i)=log(max(v(i),tiny(1.0_dp)));end do;end if
   end function

   subroutine mixspcauchy_mle(x,g,probs,rho,mu,loglik,cluster,tol,maxiters)
      real(dp),intent(in)::x(:,:);integer,intent(in)::g
      real(dp),intent(out)::probs(g),rho(g),mu(g,size(x,2)),loglik
      integer,intent(out),optional::cluster(size(x,1));real(dp),intent(in),optional::tol;integer,intent(in),optional::maxiters
      real(dp)::w(size(x,1),g),ld(size(x,1),g),mx,rs,old,eps,mutmp(size(x,2))
      integer::i,j,it,nit,n,cl
      n=size(x,1);eps=1e-6_dp;if(present(tol))eps=tol;nit=100;if(present(maxiters))nit=maxiters
      call init_directional_clusters(x,g,w,mu);old=-huge(1.0_dp)
      do it=1,nit
         do j=1,g
            probs(j)=max(sum(w(:,j))/real(n,dp),1e-12_dp)
            mutmp=mu(j,:);call weighted_spcauchy(x,w(:,j),mutmp,rho(j),eps,80);mu(j,:)=mutmp
            ld(:,j)=log(probs(j))+dspcauchy(x,mu(j,:),rho(j),.true.)
         end do
         loglik=0.0_dp
         do i=1,n;mx=maxval(ld(i,:));rs=sum(exp(ld(i,:)-mx));loglik=loglik+mx+log(rs);w(i,:)=exp(ld(i,:)-mx)/rs;end do
         if(abs(loglik-old)<eps)exit;old=loglik
      end do
      probs=sum(w,dim=1)/real(n,dp)
      if(present(cluster))then;do i=1,n;cl=maxloc(w(i,:),dim=1);cluster(i)=cl;end do;end if
   end subroutine

   subroutine mixpkbd_mle(x,g,probs,rho,mu,loglik,cluster,tol,maxiters)
      real(dp),intent(in)::x(:,:);integer,intent(in)::g
      real(dp),intent(out)::probs(g),rho(g),mu(g,size(x,2)),loglik
      integer,intent(out),optional::cluster(size(x,1));real(dp),intent(in),optional::tol;integer,intent(in),optional::maxiters
      real(dp)::w(size(x,1),g),ld(size(x,1),g),mx,rs,old,eps,mutmp(size(x,2))
      integer::i,j,it,nit,n,cl
      n=size(x,1);eps=1e-6_dp;if(present(tol))eps=tol;nit=100;if(present(maxiters))nit=maxiters
      call init_directional_clusters(x,g,w,mu);old=-huge(1.0_dp)
      do it=1,nit
         do j=1,g
            probs(j)=max(sum(w(:,j))/real(n,dp),1e-12_dp)
            mutmp=mu(j,:);call weighted_pkbd(x,w(:,j),mutmp,rho(j),eps,100);mu(j,:)=mutmp
            ld(:,j)=log(probs(j))+dpkbd(x,mu(j,:),rho(j),.true.)
         end do
         loglik=0.0_dp
         do i=1,n;mx=maxval(ld(i,:));rs=sum(exp(ld(i,:)-mx));loglik=loglik+mx+log(rs);w(i,:)=exp(ld(i,:)-mx)/rs;end do
         if(abs(loglik-old)<eps)exit;old=loglik
      end do
      probs=sum(w,dim=1)/real(n,dp)
      if(present(cluster))then;do i=1,n;cl=maxloc(w(i,:),dim=1);cluster(i)=cl;end do;end if
   end subroutine

   subroutine init_directional_clusters(x,g,w,mu)
      real(dp),intent(in)::x(:,:);integer,intent(in)::g;real(dp),intent(out)::w(size(x,1),g),mu(g,size(x,2))
      real(dp)::r,mx;integer::i,j,cl
      mu=0.0_dp;mu(1,:)=x(1,:)
      do j=2,g
         cl=1;mx=-1.0_dp
         do i=1,size(x,1);r=minval(1.0_dp-matmul(mu(1:j-1,:),x(i,:)));if(r>mx)then;mx=r;cl=i;end if;end do
         mu(j,:)=x(cl,:)
      end do
      w=0.0_dp;do i=1,size(x,1);j=maxloc(matmul(mu,x(i,:)),dim=1);w(i,j)=1.0_dp;end do
   end subroutine

   subroutine weighted_spcauchy(x,w,mu,rho,tol,maxit)
      real(dp),intent(in)::x(:,:),w(:),tol;integer,intent(in)::maxit;real(dp),intent(out)::mu(size(x,2)),rho
      real(dp)::m(size(x,2)),down(size(x,1)),a,b,c,d,fc,fd,gr,old,nr,sw
      integer::it
      sw=sum(w);m=matmul(transpose(x),w);nr=sqrt(sum(m*m));if(nr<=tiny(1.0_dp))m=x(1,:);mu=m/max(sqrt(sum(m*m)),tiny(1.0_dp));rho=0.25_dp;old=-huge(1.0_dp)
      do it=1,maxit
         a=1e-8_dp;b=1.0_dp-1e-8_dp;gr=(sqrt(5.0_dp)-1.0_dp)/2.0_dp;c=b-gr*(b-a);d=a+gr*(b-a)
         fc=spc_wll(x,w,mu,c,sw);fd=spc_wll(x,w,mu,d,sw)
         do while(abs(b-a)>tol)
            if(fc>fd)then;b=d;d=c;fd=fc;c=b-gr*(b-a);fc=spc_wll(x,w,mu,c,sw)
            else;a=c;c=d;fc=fd;d=a+gr*(b-a);fd=spc_wll(x,w,mu,d,sw);end if
         end do
         rho=0.5_dp*(a+b);down=1.0_dp+rho*rho-2.0_dp*rho*matmul(x,mu)
         m=matmul(transpose(x),w/down);nr=sqrt(sum(m*m));if(nr>tiny(1.0_dp))mu=m/nr
         if(abs(max(fc,fd)-old)<tol)exit;old=max(fc,fd)
      end do
   end subroutine

   pure real(dp) function spc_wll(x,w,mu,rho,sw) result(v)
      real(dp),intent(in)::x(:,:),w(:),mu(:),rho,sw;integer::d
      d=size(x,2)-1;v=real(d,dp)*sw*log(max(tiny(1.0_dp),1.0_dp-rho*rho))-real(d,dp)*sum(w*log(max(tiny(1.0_dp),1.0_dp+rho*rho-2.0_dp*rho*matmul(x,mu))))
   end function

   subroutine weighted_pkbd(x,w,mu,rho,tol,maxit)
      real(dp),intent(in)::x(:,:),w(:),tol;integer,intent(in)::maxit;real(dp),intent(out)::mu(size(x,2)),rho
      real(dp)::mes(size(x,2)),grad(size(x,2)),cand(size(x,2)),h,val,base,step,gamma,com,nr
      integer::j,it
      mes=matmul(transpose(x),w)/max(sum(w),tiny(1.0_dp));base=wpk_mes_ll(x,w,mes);step=0.2_dp
      do it=1,maxit
         do j=1,size(mes)
            h=1e-5_dp*max(1.0_dp,abs(mes(j)));cand=mes;cand(j)=cand(j)+h;val=wpk_mes_ll(x,w,cand);cand(j)=mes(j)-h
            grad(j)=(val-wpk_mes_ll(x,w,cand))/(2.0_dp*h)
         end do
         nr=sqrt(sum(grad*grad));if(nr<=tiny(1.0_dp))exit;cand=mes+step*grad/nr;val=wpk_mes_ll(x,w,cand)
         if(val>base)then;if(abs(val-base)<tol)then;mes=cand;base=val;exit;end if;mes=cand;base=val;step=min(1.0_dp,1.1_dp*step)
         else;step=0.5_dp*step;if(step<1e-10_dp)exit;end if
      end do
      gamma=sqrt(sum(mes*mes));com=sqrt(gamma*gamma+1.0_dp);rho=(com-1.0_dp)/max(gamma,tiny(1.0_dp));mu=mes/max(gamma,tiny(1.0_dp))
   end subroutine

   pure real(dp) function wpk_mes_ll(x,w,mes) result(v)
      real(dp),intent(in)::x(:,:),w(:),mes(:);real(dp)::g2,com,sw;integer::d
      d=size(x,2)-1;sw=sum(w);g2=max(sum(mes*mes),1e-16_dp);com=sqrt(g2+1.0_dp)
      v=-0.5_dp*real(d+1,dp)*sum(w*log(max(tiny(1.0_dp),com-matmul(x,mes))))-0.5_dp*sw*real(d-1,dp)*(log(com-1.0_dp)-log(g2))
   end function

   real(dp) function solve_kappa(r,p) result(k)
      real(dp),intent(in)::r;integer,intent(in)::p;real(dp)::a,k2,nu;integer::it
      if(r<1e-10_dp)then;k=0;return;end if;k=max(1e-8_dp,r*(p-r*r)/max(1e-10_dp,1-r*r));nu=0.5_dp*p-1
      do it=1,100;a=exp(log_bessel_i(nu+1,k)-log_bessel_i(nu,k));k2=max(1e-10_dp,k-(a-r)/max(1e-12_dp,1-a*a-(p-1)*a/k));if(abs(k2-k)<1e-9_dp*max(1.0_dp,k))exit;k=k2;end do;k=k2
   end function
end module directional_mixtures
