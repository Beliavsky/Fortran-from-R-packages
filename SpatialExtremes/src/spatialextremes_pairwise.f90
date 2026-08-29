module spatialextremes_pairwise
   use spatialextremes_base, only: dp,pair_count,pair_indices,is_finite,neg_huge
   use r_compat, only: dnorm,normal_cdf,dt,pt
   implicit none
   private
   public :: lplik_smith,lplik_schlather,lplik_schlather_ind,lplik_extremalt
   public :: lplik_smith_contributions,lplik_schlather_contributions
   public :: lplik_schlather_ind_contributions,lplik_extremalt_contributions
contains
   function lplik_schlather_contributions(data,rho,logjac,weights) result(cn)
      real(dp),intent(in)::data(:,:),rho(:),logjac(:,:)
      real(dp),intent(in),optional::weights(:)
      real(dp)::cn(size(data,1),pair_count(size(data,2)))
      integer::k,i,j,t
      real(dp)::x,y,c,lF,d1,d2,dm,w
      cn=0.0_dp
      do k=1,pair_count(size(data,2))
         call pair_indices(k,size(data,2),i,j)
         w=1.0_dp
         if(present(weights))w=weights(k)
         if(w==0.0_dp)cycle
         do t=1,size(data,1)
            x=data(t,i)
            y=data(t,j)
            if(.not.(is_finite(x).and.is_finite(y)))cycle
            if(rho(k)>.99999996_dp)then
               if(x>=y)then
               cn(t,k)=w*(-2*log(y)-1/y+logjac(t,i)+logjac(t,j))
               else
               cn(t,k)=w*(-2*log(x)-1/x+logjac(t,i)+logjac(t,j))
               end if
            else
               c=sqrt(x*x+y*y-2*x*y*rho(k))
               lF=-(x+y+c)/(2*x*y)
               d1=-(rho(k)*x-c-y)/(2*c*x*x)
               d2=-(rho(k)*y-c-x)/(2*c*y*y)
               dm=(1-rho(k)*rho(k))/(2*c**3)+d1*d2
               if(dm<=0)then
               cn(t,k)=neg_huge
               else
               cn(t,k)=w*(log(dm)+lF+logjac(t,i)+logjac(t,j))
               end if
            end if
         end do
      end do
   end function lplik_schlather_contributions

   real(dp) function lplik_schlather(data,rho,logjac,weights) result(ll)
      real(dp),intent(in)::data(:,:),rho(:),logjac(:,:)
      real(dp),intent(in),optional::weights(:)
      real(dp)::cn(size(data,1),pair_count(size(data,2)))
      cn=lplik_schlather_contributions(data,rho,logjac,weights)
      if(any(cn<=0.5_dp*neg_huge))then
      ll=neg_huge
      else
      ll=sum(cn)
      end if
   end function lplik_schlather

   function lplik_smith_contributions(data,a,logjac,weights) result(cn)
      real(dp),intent(in)::data(:,:),a(:),logjac(:,:)
      real(dp),intent(in),optional::weights(:)
      real(dp)::cn(size(data,1),pair_count(size(data,2)))
      integer::k,i,j,t
      real(dp)::x,y,ia,c1,c2,ix,iy,mix,phi1,phi2,P1,P2,lF,d1,d2,dm,w
      cn=0.0_dp
      do k=1,pair_count(size(data,2))
         call pair_indices(k,size(data,2),i,j)
         w=1.0_dp
         if(present(weights))w=weights(k)
         if(w==0.0_dp)cycle
         ia=1.0_dp/a(k)
         do t=1,size(data,1)
            x=data(t,i)
            y=data(t,j)
            if(.not.(is_finite(x).and.is_finite(y)))cycle
            ix=1/x
            iy=1/y
            mix=ix*iy*ia
            c1=log(y*ix)*ia+0.5_dp*a(k)
            c2=a(k)-c1
            if(c1>38 .and. c2< -38)then
            cn(t,k)=w*(2*log(ix)-ix+logjac(t,i)+logjac(t,j))
            cycle
            end if
            if(c1< -38 .and. c2>38)then
            cn(t,k)=w*(2*log(iy)-iy+logjac(t,i)+logjac(t,j))
            cycle
            end if
            if(c1>38 .and. c2>38)then
            cn(t,k)=w*(2*log(ix*iy)-ix-iy+logjac(t,i)+logjac(t,j))
            cycle
            end if
            phi1=dnorm(c1,0.0_dp,1.0_dp,.false.)
            phi2=dnorm(c2,0.0_dp,1.0_dp,.false.)
            P1=normal_cdf(c1)
            P2=normal_cdf(c2)
            lF=-P1*ix-P2*iy
            d1=(phi1*ia+P1)*ix*ix-phi2*mix
            d2=(phi2*ia+P2)*iy*iy-phi1*mix
            dm=d1*d2+(y*c2*phi1+x*c1*phi2)*mix*mix
            if(dm<=0)then
            cn(t,k)=neg_huge
            else
            cn(t,k)=w*(log(dm)+lF+logjac(t,i)+logjac(t,j))
            end if
         end do
      end do
   end function lplik_smith_contributions

   real(dp) function lplik_smith(data,a,logjac,weights) result(ll)
      real(dp),intent(in)::data(:,:),a(:),logjac(:,:)
      real(dp),intent(in),optional::weights(:)
      real(dp)::cn(size(data,1),pair_count(size(data,2)))
      cn=lplik_smith_contributions(data,a,logjac,weights)
      if(any(cn<=0.5_dp*neg_huge))then
      ll=neg_huge
      else
      ll=sum(cn)
      end if
   end function lplik_smith

   function lplik_schlather_ind_contributions(data,alpha,rho,logjac,weights) result(cn)
      real(dp),intent(in)::data(:,:),alpha,rho(:),logjac(:,:)
      real(dp),intent(in),optional::weights(:)
      real(dp)::cn(size(data,1),pair_count(size(data,2)))
      integer::k,i,j,t
      real(dp)::x,y,c,lF,d1,d2,dm,w
      if(alpha<=epsilon(1.0_dp))then
      cn=lplik_schlather_contributions(data,rho,logjac,weights)
      return
      end if
      cn=0.0_dp
      do k=1,pair_count(size(data,2))
         call pair_indices(k,size(data,2),i,j)
         w=1.0_dp
         if(present(weights))w=weights(k)
         if(w==0.0_dp)cycle
         do t=1,size(data,1)
            x=data(t,i)
            y=data(t,j)
            if(.not.(is_finite(x).and.is_finite(y)))cycle
            if(alpha>=1.0_dp-epsilon(1.0_dp))then
               cn(t,k)=w*(-1/x-1/y-2*log(x*y)+logjac(t,i)+logjac(t,j))
               cycle
            end if
            c=sqrt(x*x+y*y-2*x*y*rho(k))
            lF=((-alpha-1)*(x+y)+(alpha-1)*c)/(2*x*y)
            d1=(alpha-1)*(rho(k)*x-c-y)/(2*c*x*x)+alpha/(x*x)
            d2=(alpha-1)*(rho(k)*y-c-x)/(2*c*y*y)+alpha/(y*y)
            dm=(1-alpha)*(1-rho(k)*rho(k))/(2*c**3)+d1*d2
            if(dm<=0)then
            cn(t,k)=neg_huge
            else
            cn(t,k)=w*(log(dm)+lF+logjac(t,i)+logjac(t,j))
            end if
         end do
      end do
   end function lplik_schlather_ind_contributions

   real(dp) function lplik_schlather_ind(data,alpha,rho,logjac,weights) result(ll)
      real(dp),intent(in)::data(:,:),alpha,rho(:),logjac(:,:)
      real(dp),intent(in),optional::weights(:)
      real(dp)::cn(size(data,1),pair_count(size(data,2)))
      cn=lplik_schlather_ind_contributions(data,alpha,rho,logjac,weights)
      if(any(cn<=0.5_dp*neg_huge))then
      ll=neg_huge
      else
      ll=sum(cn)
      end if
   end function lplik_schlather_ind

   function lplik_extremalt_contributions(data,rho,nu,logjac,weights) result(cn)
      real(dp),intent(in)::data(:,:),rho(:),nu,logjac(:,:)
      real(dp),intent(in),optional::weights(:)
      real(dp)::cn(size(data,1),pair_count(size(data,2)))
      integer::k,i,j,t
      real(dp)::x,y,ix,iy,r21,r12,c1,c2,tc1,tc2,T1,T2,a,idf,df1,lF,d1,d2,dc1,dc2,mix,dm,w
      idf=1/nu
      df1=nu+1
      cn=0.0_dp
      do k=1,pair_count(size(data,2))
         call pair_indices(k,size(data,2),i,j)
         w=1.0_dp
         if(present(weights))w=weights(k)
         if(w==0.0_dp)cycle
         a=sqrt(df1/(1-rho(k)*rho(k)))
         do t=1,size(data,1)
            x=data(t,i)
            y=data(t,j)
            if(.not.(is_finite(x).and.is_finite(y)))cycle
            ix=1/x
            iy=1/y
            r21=(y*ix)**idf
            r12=1/r21
            c1=(r21-rho(k))*a
            c2=(r12-rho(k))*a
            tc1=dt(c1,df1,.false.)
            tc2=dt(c2,df1,.false.)
            T1=pt(c1,df1)
            T2=pt(c2,df1)
            if(T1<=0)then
            cn(t,k)=w*(2*log(iy)-iy+logjac(t,j))
            cycle
            end if
            if(T2<=0)then
            cn(t,k)=w*(2*log(ix)-ix+logjac(t,i))
            cycle
            end if
            lF=-T1*ix-T2*iy
            d1=ix*(ix*T1+a*idf*(ix*r21*tc1-iy*r12*tc2))
            d2=iy*(iy*T2+a*idf*(iy*r12*tc2-ix*r21*tc1))
            dc1=-(1+1/df1)*c1/(1+c1*c1/df1)*tc1
            dc2=-(1+1/df1)*c2/(1+c2*c2/df1)*tc2
            mix=ix*iy*idf*idf*a*(df1*(r21*ix*tc1+r12*iy*tc2)+(r21*r21*ix*dc1+r12*r12*iy*dc2)*a)
            dm=mix+d1*d2
            if(dm<=0)then
            cn(t,k)=neg_huge
            else
            cn(t,k)=w*(log(dm)+lF+logjac(t,i)+logjac(t,j))
            end if
         end do
      end do
   end function lplik_extremalt_contributions

   real(dp) function lplik_extremalt(data,rho,nu,logjac,weights) result(ll)
      real(dp),intent(in)::data(:,:),rho(:),nu,logjac(:,:)
      real(dp),intent(in),optional::weights(:)
      real(dp)::cn(size(data,1),pair_count(size(data,2)))
      cn=lplik_extremalt_contributions(data,rho,nu,logjac,weights)
      if(any(cn<=0.5_dp*neg_huge))then
      ll=neg_huge
      else
      ll=sum(cn)
      end if
   end function lplik_extremalt
end module spatialextremes_pairwise
