module urca_restrictions
   use urca_kinds, only : dp
   use urca_types, only : johansen_result, restriction_result
   use urca_linalg, only : invert_matrix, invert_spd, chol_lower, symmetric_eigen, orthogonal_complement
   use urca_distributions, only : chi_square_sf
   implicit none
   private
   public :: beta_restriction_test, alpha_restriction_test, alpha_beta_restriction_test
   public :: partly_known_beta_test, iterated_partly_known_beta_test, linear_trend_lr_test
contains
   subroutine moments(z,m00,m11,mkk,m01,m0k,mk0,m10,m1k,mk1,s00,s0k,sk0,skk,m11i,info)
      type(johansen_result),intent(in)::z
      real(dp),allocatable,intent(out)::m00(:,:),m11(:,:),mkk(:,:),m01(:,:),m0k(:,:),mk0(:,:),m10(:,:),m1k(:,:), &
         & mk1(:,:)
      real(dp),allocatable,intent(out)::s00(:,:),s0k(:,:),sk0(:,:),skk(:,:),m11i(:,:)
      integer,intent(out)::info
      integer::n
      n=size(z%z0,1)
      info=0
      m00=matmul(transpose(z%z0),z%z0)/real(n,dp)
      m11=matmul(transpose(z%z1),z%z1)/real(n,dp)
      mkk=matmul(transpose(z%zk),z%zk)/real(n,dp)
      m01=matmul(transpose(z%z0),z%z1)/real(n,dp)
      m10=transpose(m01)
      m0k=matmul(transpose(z%z0),z%zk)/real(n,dp)
      mk0=transpose(m0k)
      m1k=matmul(transpose(z%z1),z%zk)/real(n,dp)
      mk1=transpose(m1k)
      call invert_spd(m11,m11i,info)
      if(info/=0)call invert_matrix(m11,m11i,info)
      if(info/=0)return
      s00=m00-matmul(m01,matmul(m11i,m10))
      s0k=m0k-matmul(m01,matmul(m11i,m1k))
      sk0=transpose(s0k)
      skk=mkk-matmul(mk1,matmul(m11i,m1k))
   end subroutine moments

   subroutine whiten_eigen(m,b,eval,evec,info)
! Solve the eigenproblem M^{-1/2} B M^{-1/2}' for SPD M.
      real(dp),intent(in)::m(:,:),b(:,:)
      real(dp),allocatable,intent(out)::eval(:),evec(:,:)
      integer,intent(out)::info
      real(dp),allocatable::l(:,:),li(:,:),a(:,:)
      call chol_lower(m,l,info)
      if(info/=0)return
      call invert_matrix(l,li,info)
      if(info/=0)return
      a=matmul(li,matmul(b,transpose(li)))
      call symmetric_eigen(a,eval,evec,info,.true.)
   end subroutine whiten_eigen

   subroutine finish_beta_like(s00,s0k,sk0,skk,m01,mk1,m11i,v,vorg,out,info)
      real(dp),intent(in)::s00(:,:),s0k(:,:),sk0(:,:),skk(:,:),m01(:,:),mk1(:,:),m11i(:,:),v(:,:),vorg(:,:)
      type(restriction_result),intent(inout)::out
      integer,intent(out)::info
      real(dp),allocatable::vv(:,:),vvi(:,:)
      vv=matmul(transpose(v),matmul(skk,v))
      call invert_matrix(vv,vvi,info)
      if(info/=0)return
      out%w=matmul(s0k,matmul(v,vvi))
      out%pi=matmul(out%w,transpose(v))
      out%delta=s00-matmul(s0k,matmul(v,matmul(vvi,matmul(transpose(v),sk0))))
      out%gamma=matmul(m01,m11i)-matmul(out%pi,matmul(mk1,m11i))
      out%v=v
      out%vorg=vorg
   end subroutine finish_beta_like

   function beta_restriction_test(z,h,r) result(out)
      type(johansen_result),intent(in)::z
      real(dp),intent(in)::h(:,:)
      integer,intent(in)::r
      type(restriction_result)::out
      real(dp),allocatable::m00(:,:),m11(:,:),mkk(:,:),m01(:,:),m0k(:,:),mk0(:,:),m10(:,:),m1k(:,:),mk1(:,:),s00(:, &
         & :),s0k(:,:),sk0(:,:),skk(:,:),m11i(:,:)
      real(dp),allocatable::hh(:,:),l(:,:),li(:,:),s00i(:,:),a(:,:),ev(:),e(:,:),vorg(:,:),v(:,:),b(:,:)
      integer::info,j,q,n
      call moments(z,m00,m11,mkk,m01,m0k,mk0,m10,m1k,mk1,s00,s0k,sk0,skk,m11i,info)
      if(info/=0)then
      out%info=info
      return
      end if
      q=size(z%zk,2)
      n=size(z%z0,1)
      if(r<1.or.r>=z%p.or.size(h,1)/=q)then
      out%info=-1
      return
      end if
      hh=matmul(transpose(h),matmul(skk,h))
      call chol_lower(hh,l,info)
      if(info/=0)then
      out%info=100+info
      return
      end if
      call invert_matrix(l,li,info)
      if(info/=0)then
      out%info=200+info
      return
      end if
      call invert_spd(s00,s00i,info)
      if(info/=0)call invert_matrix(s00,s00i,info)
      if(info/=0)then
      out%info=300+info
      return
      end if
      a=matmul(li,matmul(transpose(h),matmul(sk0,matmul(s00i,matmul(s0k,matmul(h,transpose(li)))))))
      call symmetric_eigen(a,ev,e,info,.true.)
      if(info/=0)then
      out%info=400+info
      return
      end if
      vorg=matmul(h,matmul(transpose(li),e))
      v=vorg
      do j=1,size(v,2)
      if(abs(v(1,j))>1e-14_dp)v(:,j)=v(:,j)/v(1,j)
      end do
      call finish_beta_like(s00,s0k,sk0,skk,m01,mk1,m11i,v,vorg,out,info)
      if(info/=0)then
      out%info=500+info
      return
      end if
      out%lambda=ev
      out%statistic=real(n,dp)*sum(log(max(1e-15_dp,1.0_dp-ev(1:r))/max(1e-15_dp,1.0_dp-z%lambda(1:r))))
      out%df=r*(q-size(h,2))
      out%p_value=chi_square_sf(max(0.0_dp,out%statistic),real(out%df,dp))
      out%info=0
   end function beta_restriction_test

   function alpha_restriction_test(z,a,r) result(out)
      type(johansen_result),intent(in)::z
      real(dp),intent(in)::a(:,:)
      integer,intent(in)::r
      type(restriction_result)::out
      real(dp),allocatable::m00(:,:),m11(:,:),mkk(:,:),m01(:,:),m0k(:,:),mk0(:,:),m10(:,:),m1k(:,:),mk1(:,:),s00(:, &
         & :),s0k(:,:),sk0(:,:),skk(:,:),m11i(:,:)
      real(dp),allocatable::b(:,:),sab(:,:),skb(:,:),sbb(:,:),sbbi(:,:),ra(:,:),rkc(:,:),saa(:,:),sak(:,:),ska(:,:), &
         & skkc(:,:),sai(:,:),ev(:),e(:,:),l(:,:),li(:,:),vorg(:,:),v(:,:),ata(:,:),atai(:,:),phi(:,:),alpha(:,:)
      integer::info,j,n,p,q
      call moments(z,m00,m11,mkk,m01,m0k,mk0,m10,m1k,mk1,s00,s0k,sk0,skk,m11i,info)
      if(info/=0)then
      out%info=info
      return
      end if
      n=size(z%z0,1)
      p=z%p
      q=size(z%zk,2)
      if(r<1.or.r>=p.or.size(a,1)/=p.or.size(a,2)>=p)then
      out%info=-1
      return
      end if
      call orthogonal_complement(a,b,info)
      if(info/=0)then
      out%info=100+info
      return
      end if
      sab=matmul(transpose(a),matmul(s00,b))
      skb=matmul(transpose(s0k),b)
      sbb=matmul(transpose(b),matmul(s00,b))
      call invert_spd(sbb,sbbi,info)
      if(info/=0)call invert_matrix(sbb,sbbi,info)
      if(info/=0)then
      out%info=200+info
      return
      end if
      ra=matmul(z%r0,a)-matmul(z%r0,matmul(b,matmul(sbbi,transpose(sab))))
      rkc=z%rk-matmul(z%r0,matmul(b,matmul(sbbi,transpose(skb))))
      saa=matmul(transpose(ra),ra)/real(n,dp)
      sak=matmul(transpose(ra),rkc)/real(n,dp)
      ska=transpose(sak)
      skkc=matmul(transpose(rkc),rkc)/real(n,dp)
      call chol_lower(skkc,l,info)
      if(info/=0)then
      out%info=300+info
      return
      end if
      call invert_matrix(l,li,info)
      if(info/=0)then
      out%info=301+info
      return
      end if
      call invert_spd(saa,sai,info)
      if(info/=0)call invert_matrix(saa,sai,info)
      if(info/=0)then
      out%info=400+info
      return
      end if
      call symmetric_eigen(matmul(li,matmul(ska,matmul(sai,matmul(sak,transpose(li))))),ev,e,info,.true.)
      if(info/=0)then
      out%info=500+info
      return
      end if
      vorg=matmul(transpose(li),e(:,1:r))
      v=vorg
      do j=1,r
      if(abs(v(1,j))>1e-14_dp)v(:,j)=v(:,j)/v(1,j)
      end do
      ata=matmul(transpose(a),a)
      call invert_spd(ata,atai,info)
      if(info/=0)call invert_matrix(ata,atai,info)
      if(info/=0)then
      out%info=600+info
      return
      end if
      phi=matmul(atai,matmul(sak,vorg))
      alpha=matmul(a,phi)
      do j=1,r
      alpha(:,j)=alpha(:,j)*vorg(1,j)
      end do
      out%w=alpha
      out%pi=matmul(alpha,transpose(v))
      out%gamma=matmul(m01,m11i)-matmul(out%pi,matmul(mk1,m11i))
      out%delta_bb=sbb
      out%delta_ab=sab-matmul(transpose(a),matmul(alpha,matmul(transpose(v),skb)))
      out%delta_aa_b=saa-matmul(transpose(a),matmul(alpha,matmul(transpose(alpha),a)))
      out%v=v
      out%vorg=vorg
      out%lambda=ev
      out%statistic=real(n,dp)*sum(log(max(1e-15_dp,1.0_dp-ev(1:r))/max(1e-15_dp,1.0_dp-z%lambda(1:r))))
      out%df=r*(p-size(a,2))
      out%p_value=chi_square_sf(max(0.0_dp,out%statistic),real(out%df,dp))
      out%info=0
   end function alpha_restriction_test

   function alpha_beta_restriction_test(z,h,a,r) result(out)
      type(johansen_result),intent(in)::z
      real(dp),intent(in)::h(:,:),a(:,:)
      integer,intent(in)::r
      type(restriction_result)::out
      real(dp),allocatable::m00(:,:),m11(:,:),mkk(:,:),m01(:,:),m0k(:,:),mk0(:,:),m10(:,:),m1k(:,:),mk1(:,:),s00(:, &
         & :),s0k(:,:),sk0(:,:),skk(:,:),m11i(:,:)
      real(dp),allocatable::b(:,:),sab(:,:),skb(:,:),sbb(:,:),sbbi(:,:),ra(:,:),rkc(:,:),saa(:,:),sak(:,:),ska(:,:), &
         & skkc(:,:),sai(:,:),hh(:,:),l(:,:),li(:,:),ev(:),e(:,:),vorg(:,:),v(:,:),ata(:,:),atai(:,:),phi(:,:), &
         & alpha(:,:)
      integer::info,j,n,p,q
      call moments(z,m00,m11,mkk,m01,m0k,mk0,m10,m1k,mk1,s00,s0k,sk0,skk,m11i,info)
      if(info/=0)then
      out%info=info
      return
      end if
      n=size(z%z0,1)
      p=z%p
      q=size(z%zk,2)
      if(r<1.or.r>=p.or.size(a,1)/=p.or.size(a,2)>=p.or.size(h,1)/=q)then
      out%info=-1
      return
      end if
      call orthogonal_complement(a,b,info)
      if(info/=0)then
      out%info=100+info
      return
      end if
      sab=matmul(transpose(a),matmul(s00,b))
      skb=matmul(transpose(s0k),b)
      sbb=matmul(transpose(b),matmul(s00,b))
      call invert_spd(sbb,sbbi,info)
      if(info/=0)call invert_matrix(sbb,sbbi,info)
      if(info/=0)then
      out%info=200+info
      return
      end if
      ra=matmul(z%r0,a)-matmul(z%r0,matmul(b,matmul(sbbi,transpose(sab))))
      rkc=z%rk-matmul(z%r0,matmul(b,matmul(sbbi,transpose(skb))))
      saa=matmul(transpose(ra),ra)/real(n,dp)
      sak=matmul(transpose(ra),rkc)/real(n,dp)
      ska=transpose(sak)
      skkc=matmul(transpose(rkc),rkc)/real(n,dp)
      hh=matmul(transpose(h),matmul(skkc,h))
      call chol_lower(hh,l,info)
      if(info/=0)then
      out%info=300+info
      return
      end if
      call invert_matrix(l,li,info)
      if(info/=0)then
      out%info=301+info
      return
      end if
      call invert_spd(saa,sai,info)
      if(info/=0)call invert_matrix(saa,sai,info)
      if(info/=0)then
      out%info=400+info
      return
      end if
      call symmetric_eigen(matmul(li,matmul(transpose(h),matmul(ska,matmul(sai,matmul(sak,matmul(h, &
         & transpose(li))))))),ev,e,info,.true.)
      if(info/=0)then
      out%info=500+info
      return
      end if
      vorg=matmul(h,matmul(transpose(li),e(:,1:r)))
      v=vorg
      do j=1,r
      if(abs(v(1,j))>1e-14_dp)v(:,j)=v(:,j)/v(1,j)
      end do
      ata=matmul(transpose(a),a)
      call invert_spd(ata,atai,info)
      if(info/=0)call invert_matrix(ata,atai,info)
      if(info/=0)then
      out%info=600+info
      return
      end if
      phi=matmul(atai,matmul(sak,vorg))
      alpha=matmul(a,phi)
      do j=1,r
      alpha(:,j)=alpha(:,j)*vorg(1,j)
      end do
      out%w=alpha
      out%pi=matmul(alpha,transpose(v))
      out%gamma=matmul(m01,m11i)-matmul(out%pi,matmul(mk1,m11i))
      out%delta_bb=sbb
      out%delta_ab=sab-matmul(transpose(a),matmul(alpha,matmul(transpose(v),skb)))
      out%delta_aa_b=saa-matmul(transpose(a),matmul(alpha,matmul(transpose(alpha),a)))
      out%v=v
      out%vorg=vorg
      out%lambda=ev
      out%statistic=real(n,dp)*sum(log(max(1e-15_dp,1.0_dp-ev(1:r))/max(1e-15_dp,1.0_dp-z%lambda(1:r))))
      out%df=r*(p-size(a,2))+r*(p-size(h,2))
      out%p_value=chi_square_sf(max(0.0_dp,out%statistic),real(out%df,dp))
      out%info=0
   end function alpha_beta_restriction_test

   subroutine residualize_on_beta(s00,s0k,sk0,skk,beta,s00b,s0kb,sk0b,skkb,info)
      real(dp),intent(in)::s00(:,:),s0k(:,:),sk0(:,:),skk(:,:),beta(:,:)
      real(dp),allocatable,intent(out)::s00b(:,:),s0kb(:,:),sk0b(:,:),skkb(:,:)
      integer,intent(out)::info
      real(dp),allocatable::g(:,:),gi(:,:)
      g=matmul(transpose(beta),matmul(skk,beta))
      call invert_matrix(g,gi,info)
      if(info/=0)return
      s00b=s00-matmul(s0k,matmul(beta,matmul(gi,matmul(transpose(beta),sk0))))
      s0kb=s0k-matmul(s0k,matmul(beta,matmul(gi,matmul(transpose(beta),skk))))
      sk0b=sk0-matmul(skk,matmul(beta,matmul(gi,matmul(transpose(beta),sk0))))
      skkb=skk-matmul(skk,matmul(beta,matmul(gi,matmul(transpose(beta),skk))))
   end subroutine residualize_on_beta

   subroutine psd_whitener(a,nkeep,c,info)
      real(dp),intent(in)::a(:,:)
      integer,intent(in)::nkeep
      real(dp),allocatable,intent(out)::c(:,:)
      integer,intent(out)::info
      real(dp),allocatable::ev(:),e(:,:)
      integer::j
      call symmetric_eigen(a,ev,e,info,.true.)
      if(info/=0)return
      if(nkeep<1.or.nkeep>size(ev).or.any(ev(1:nkeep)<=1e-12_dp))then
      info=-1
      return
      end if
      allocate(c(size(a,1),nkeep))
      c=e(:,1:nkeep)
      do j=1,nkeep
      c(:,j)=c(:,j)/sqrt(ev(j))
      end do
   end subroutine psd_whitener

   function partly_known_beta_test(z,h,r) result(out)
      type(johansen_result),intent(in)::z
      real(dp),intent(in)::h(:,:)
      integer,intent(in)::r
      type(restriction_result)::out
      real(dp),allocatable::m00(:,:),m11(:,:),mkk(:,:),m01(:,:),m0k(:,:),mk0(:,:),m10(:,:),m1k(:,:),mk1(:,:),s00(:, &
         & :),s0k(:,:),sk0(:,:),skk(:,:),m11i(:,:)
      real(dp),allocatable::s00h(:,:),s0kh(:,:),sk0h(:,:),skkh(:,:),c(:,:),s00hi(:,:),ev(:),e(:,:),psi(:,:),vorg(:, &
         & :),v(:,:),hh(:,:),l(:,:),li(:,:),s00i(:,:),rho(:),re(:,:)
      integer::info,j,n,q,r1,r2
      call moments(z,m00,m11,mkk,m01,m0k,mk0,m10,m1k,mk1,s00,s0k,sk0,skk,m11i,info)
      if(info/=0)then
      out%info=info
      return
      end if
      n=size(z%z0,1)
      q=size(z%zk,2)
      r1=size(h,2)
      r2=r-r1
      if(r<1.or.r>=z%p.or.size(h,1)/=q.or.r1<1.or.r2<1)then
      out%info=-1
      return
      end if
      call residualize_on_beta(s00,s0k,sk0,skk,h,s00h,s0kh,sk0h,skkh,info)
      if(info/=0)then
      out%info=100+info
      return
      end if
      call psd_whitener(skkh,q-r1,c,info)
      if(info/=0)then
      out%info=200+info
      return
      end if
      call invert_spd(s00h,s00hi,info)
      if(info/=0)call invert_matrix(s00h,s00hi,info)
      if(info/=0)then
      out%info=300+info
      return
      end if
      call symmetric_eigen(matmul(transpose(c),matmul(sk0h,matmul(s00hi,matmul(s0kh,c)))),ev,e,info,.true.)
      if(info/=0)then
      out%info=400+info
      return
      end if
      psi=matmul(c,e(:,1:r2))
      allocate(vorg(q,r))
      vorg(:,1:r1)=h
      vorg(:,r1+1:r)=psi
      v=vorg
      do j=r1+1,r
      if(abs(v(1,j))>1e-14_dp)v(:,j)=v(:,j)/v(1,j)
      end do
      call finish_beta_like(s00,s0k,sk0,skk,m01,mk1,m11i,v,vorg,out,info)
      if(info/=0)then
      out%info=500+info
      return
      end if
      hh=matmul(transpose(h),matmul(skk,h))
      call chol_lower(hh,l,info)
      if(info/=0)then
      out%info=600+info
      return
      end if
      call invert_matrix(l,li,info)
      if(info/=0)then
      out%info=601+info
      return
      end if
      call invert_spd(s00,s00i,info)
      if(info/=0)call invert_matrix(s00,s00i,info)
      if(info/=0)then
      out%info=700+info
      return
      end if
      call symmetric_eigen(matmul(li,matmul(transpose(h),matmul(sk0,matmul(s00i,matmul(s0k,matmul(h, &
         & transpose(li))))))),rho,re,info,.true.)
      if(info/=0)then
      out%info=800+info
      return
      end if
      out%lambda=ev
      out%statistic=real(n,dp)*(sum(log(max(1e-15_dp,1.0_dp-rho(1:r1))))+sum(log(max(1e-15_dp, &
         & 1.0_dp-ev(1:r2))))-sum(log(max(1e-15_dp,1.0_dp-z%lambda(1:r)))))
      out%df=(q-r)*r1
      out%p_value=chi_square_sf(max(0.0_dp,out%statistic),real(out%df,dp))
      out%info=0
   end function partly_known_beta_test

   function iterated_partly_known_beta_test(z,h,r,r1,conv,max_iter) result(out)
      type(johansen_result),intent(in)::z
      real(dp),intent(in)::h(:,:)
      integer,intent(in)::r,r1
      real(dp),intent(in),optional::conv
      integer,intent(in),optional::max_iter
      type(restriction_result)::out
      real(dp),allocatable::m00(:,:),m11(:,:),mkk(:,:),m01(:,:),m0k(:,:),mk0(:,:),m10(:,:),m1k(:,:),mk1(:,:),s00(:, &
         & :),s0k(:,:),sk0(:,:),skk(:,:),m11i(:,:)
      real(dp),allocatable::hh(:,:),l(:,:),li(:,:),s00i(:,:),ev(:),e(:,:),beta1(:,:),beta2(:,:),s00b(:,:),s0kb(:,:), &
         & sk0b(:,:),skkb(:,:),c(:,:),eval2(:),evec2(:,:),last(:),vorg(:,:),v(:,:),rho(:),re(:,:)
      real(dp)::tol,diff
      integer::mi,iter,n,q,s,r2,info,j
      call moments(z,m00,m11,mkk,m01,m0k,mk0,m10,m1k,mk1,s00,s0k,sk0,skk,m11i,info)
      if(info/=0)then
      out%info=info
      return
      end if
      n=size(z%z0,1)
      q=size(z%zk,2)
      s=size(h,2)
      r2=r-r1
      tol=1e-4_dp
      if(present(conv))tol=conv
      mi=50
      if(present(max_iter))mi=max_iter
      if(r<1.or.r>=z%p.or.r1<1.or.r2<1.or.size(h,1)/=q.or.r1>s)then
      out%info=-1
      return
      end if
      hh=matmul(transpose(h),matmul(skk,h))
      call chol_lower(hh,l,info)
      if(info/=0)then
      out%info=100+info
      return
      end if
      call invert_matrix(l,li,info)
      if(info/=0)then
      out%info=101+info
      return
      end if
      call invert_spd(s00,s00i,info)
      if(info/=0)call invert_matrix(s00,s00i,info)
      if(info/=0)then
      out%info=200+info
      return
      end if
      call symmetric_eigen(matmul(li,matmul(transpose(h),matmul(sk0,matmul(s00i,matmul(s0k,matmul(h, &
         & transpose(li))))))),ev,e,info,.true.)
      if(info/=0)then
      out%info=300+info
      return
      end if
      beta1=matmul(h,e(:,1:r1))
      allocate(last(q-r1))
      last=1.0_dp
      diff=huge(1.0_dp)
      do iter=1,mi
         call residualize_on_beta(s00,s0k,sk0,skk,beta1,s00b,s0kb,sk0b,skkb,info)
         if(info/=0)exit
         call psd_whitener(skkb,q-r1,c,info)
         if(info/=0)exit
         call invert_spd(s00b,s00i,info)
         if(info/=0)call invert_matrix(s00b,s00i,info)
         if(info/=0)exit
         call symmetric_eigen(matmul(transpose(c),matmul(sk0b,matmul(s00i,matmul(s0kb,c)))),eval2,evec2,info,.true.)
         if(info/=0)exit
         diff=sum((eval2-last)**2)
         last=eval2
         beta2=matmul(c,evec2(:,1:r2))
         call residualize_on_beta(s00,s0k,sk0,skk,beta2,s00b,s0kb,sk0b,skkb,info)
         if(info/=0)exit
         hh=matmul(transpose(h),matmul(skkb,h))
         call psd_whitener(hh,s,c,info)
         if(info/=0)exit
         call invert_spd(s00b,s00i,info)
         if(info/=0)call invert_matrix(s00b,s00i,info)
         if(info/=0)exit
         call symmetric_eigen(matmul(transpose(c),matmul(transpose(h),matmul(sk0b,matmul(s00i,matmul(s0kb,matmul(h, &
            & c)))))),ev,e,info,.true.)
         if(info/=0)exit
         beta1=matmul(h,e(:,1:r1))
         if(diff<=tol)exit
      end do
      if(info/=0)then
      out%info=400+info
      return
      end if
      allocate(vorg(q,r))
      vorg(:,1:r1)=beta1
      vorg(:,r1+1:r)=beta2
      v=vorg
      do j=1,r
      if(abs(v(1,j))>1e-14_dp)v(:,j)=v(:,j)/v(1,j)
      end do
      call finish_beta_like(s00,s0k,sk0,skk,m01,mk1,m11i,v,vorg,out,info)
      if(info/=0)then
      out%info=500+info
      return
      end if
      hh=matmul(transpose(beta1),matmul(skk,beta1))
      call chol_lower(hh,l,info)
      if(info/=0)then
      out%info=600+info
      return
      end if
      call invert_matrix(l,li,info)
      if(info/=0)then
      out%info=601+info
      return
      end if
      call invert_spd(s00,s00i,info)
      if(info/=0)call invert_matrix(s00,s00i,info)
      if(info/=0)then
      out%info=700+info
      return
      end if
      call symmetric_eigen(matmul(li,matmul(transpose(beta1),matmul(sk0,matmul(s00i,matmul(s0k,matmul(beta1, &
         & transpose(li))))))),rho,re,info,.true.)
      if(info/=0)then
      out%info=800+info
      return
      end if
      out%lambda=eval2
      out%statistic=real(n,dp)*(sum(log(max(1e-15_dp,1.0_dp-rho(1:r1))))+sum(log(max(1e-15_dp, &
         & 1.0_dp-eval2(1:r2))))-sum(log(max(1e-15_dp,1.0_dp-z%lambda(1:r)))))
      out%df=(q-s-r2)*r1
      out%p_value=chi_square_sf(max(0.0_dp,out%statistic),real(out%df,dp))
      out%info=merge(1,0,iter>=mi.and.diff>tol)
   end function iterated_partly_known_beta_test

   function linear_trend_lr_test(x,k,r,season,dumvar) result(out)
      use urca_cointegration, only : johansen_test, JO_TRACE, JO_CONST, JO_NONE, JO_LONGRUN
      real(dp),intent(in)::x(:,:)
      integer,intent(in)::k,r
      integer,intent(in),optional::season
      real(dp),intent(in),optional::dumvar(:,:)
      type(restriction_result)::out
      type(johansen_result)::zc,zn
      integer::p,n,idx
      zc=johansen_test(x,JO_TRACE,JO_CONST,k,JO_LONGRUN,season,dumvar)
      zn=johansen_test(x,JO_TRACE,JO_NONE,k,JO_LONGRUN,season,dumvar)
      if(zc%info/=0.or.zn%info/=0)then
      out%info=merge(zc%info,zn%info,zc%info/=0)
      return
      end if
      p=size(x,2)
      if(r<1.or.r>=p)then
      out%info=-1
      return
      end if
      idx=r+1
      n=size(zc%z0,1)
      out%statistic=-real(n,dp)*sum(log(max(1e-15_dp,1.0_dp-zc%lambda(idx:p))/max(1e-15_dp,1.0_dp-zn%lambda(idx:p))))
      out%df=p-r
      out%p_value=chi_square_sf(max(0.0_dp,out%statistic),real(out%df,dp))
      out%info=0
   end function linear_trend_lr_test
end module urca_restrictions
