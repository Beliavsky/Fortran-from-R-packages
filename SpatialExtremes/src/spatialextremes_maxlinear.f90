module spatialextremes_maxlinear
   use spatialextremes_base, only: dp,pi,exp_rand
   use r_compat, only: runif1
   implicit none
   private
   public :: max_linear,simulate_max_linear,maxlinear_design,maxlinear_design_1d,maxlinear_design_2d
   public :: conditional_max_linear_latent,simulate_conditional_max_linear

   interface maxlinear_design
      module procedure maxlinear_design_1d
      module procedure maxlinear_design_2d
   end interface
contains
   pure function max_linear(a,z) result(x)
      real(dp),intent(in)::a(:,:),z(:)
      real(dp)::x(size(a,1))
      integer::i
      if(size(a,2)/=size(z)) then
         x=-huge(1.0_dp)
         return
      end if
      do i=1,size(a,1)
         x(i)=maxval(a(i,:)*z)
      end do
   end function max_linear

   function simulate_max_linear(n,a) result(x)
      integer,intent(in)::n
      real(dp),intent(in)::a(:,:)
      real(dp)::x(n,size(a,1)),z(size(a,2))
      integer::i,j
      do i=1,n
         do j=1,size(z)
            z(j)=1.0_dp/exp_rand()
         end do
         x(i,:)=max_linear(a,z)
      end do
   end function simulate_max_linear

   function conditional_max_linear_latent(data,a,n,info,tol) result(z)
      ! Conditional simulation of the independent unit-Frechet latent
      ! variables of a max-linear process, following upstream rcondMaxLin
      ! and the hitting-scenario construction of Stoev and Wang.
      real(dp),intent(in)::data(:),a(:,:)
      integer,intent(in)::n
      integer,intent(out),optional::info
      real(dp),intent(in),optional::tol
      real(dp)::z(n,size(a,2))
      real(dp)::bound(size(a,2)),eps,u,total
      real(dp),allocatable::weights(:)
      logical::hit(size(a,1),size(a,2)),in_class(size(a,1)),assigned(size(a,1))
      logical,allocatable::js(:),jbar(:)
      integer::i,j,r,seed,chosen,njs,istat
      logical::changed,connects

      istat=0
      z=0.0_dp
      eps=1.0e-14_dp
      if(present(tol))eps=tol
      if(size(a,1)/=size(data) .or. n<0 .or. any(data<=0.0_dp))then
         istat=1
         if(present(info))info=istat
         return
      end if

      do j=1,size(a,2)
         bound(j)=huge(1.0_dp)
         do i=1,size(a,1)
            if(a(i,j)>0.0_dp)bound(j)=min(bound(j),data(i)/a(i,j))
         end do
      end do
      hit=.false.
      do j=1,size(a,2)
         if(bound(j)>=0.5_dp*huge(1.0_dp))cycle
         do i=1,size(a,1)
            if(a(i,j)>0.0_dp)hit(i,j)=abs(a(i,j)*bound(j)-data(i))<=eps*max(1.0_dp,abs(data(i)))
         end do
      end do

      ! Draw every non-extremal latent variable from its Frechet law
      ! truncated above at its hitting bound.  This also handles columns
      ! that are zero at every conditioning site but active elsewhere.
      do r=1,n
         do j=1,size(a,2)
            u=max(tiny(1.0_dp),runif1())
            if(bound(j)>=0.5_dp*huge(1.0_dp))then
               z(r,j)=1.0_dp/(-log(u))
            else
               z(r,j)=1.0_dp/(1.0_dp/bound(j)-log(u))
            end if
         end do
      end do

      assigned=.false.
      do seed=1,size(a,1)
         if(assigned(seed))cycle
         in_class=.false.
         in_class(seed)=.true.
         changed=.true.
         do while(changed)
            changed=.false.
            do i=1,size(a,1)
               if(in_class(i).or.assigned(i))cycle
               connects=.false.
               do j=1,size(a,2)
                  if(hit(i,j).and.any(hit(:,j).and.in_class))then
                     connects=.true.
                     exit
                  end if
               end do
               if(connects)then
               in_class(i)=.true.
               changed=.true.
               end if
            end do
         end do
         assigned=assigned.or.in_class
         allocate(js(size(a,2)),jbar(size(a,2)))
         do j=1,size(a,2)
            jbar(j)=any(hit(:,j).and.in_class)
            js(j)=jbar(j).and.all(pack(hit(:,j),in_class))
         end do
         njs=count(js)
         if(njs<=0)then
            istat=2
            deallocate(js,jbar)
            exit
         end if
         allocate(weights(njs))
         r=0
         total=0.0_dp
         do j=1,size(a,2)
            if(js(j))then
               r=r+1
               weights(r)=1.0_dp/bound(j)
               total=total+weights(r)
            end if
         end do
         if(total<=0.0_dp)then
            istat=3
            deallocate(weights,js,jbar)
            exit
         end if
         weights=weights/total
         do r=1,n
            u=runif1()
            total=0.0_dp
            chosen=0
            njs=0
            do j=1,size(a,2)
               if(js(j))then
                  njs=njs+1
                  total=total+weights(njs)
                  if(chosen==0.and.u<=total)chosen=j
               end if
            end do
            if(chosen==0)then
               do j=size(a,2),1,-1
               if(js(j))then
               chosen=j
               exit
               end if
               end do
            end if
            z(r,chosen)=bound(chosen)
         end do
         deallocate(weights,js,jbar)
      end do
      if(present(info))info=istat
   end function conditional_max_linear_latent

   function simulate_conditional_max_linear(data,a_cond,a_sim,n,info,tol) result(x)
      real(dp),intent(in)::data(:),a_cond(:,:),a_sim(:,:)
      integer,intent(in)::n
      integer,intent(out),optional::info
      real(dp),intent(in),optional::tol
      real(dp)::x(n,size(a_sim,1))
      real(dp),allocatable::z(:,:)
      integer::istat,r
      if(size(a_cond,2)/=size(a_sim,2))then
         x=-huge(1.0_dp)
         if(present(info))info=1
         return
      end if
      if(present(tol))then
      z=conditional_max_linear_latent(data,a_cond,n,istat,tol)
      else
      z=conditional_max_linear_latent(data,a_cond,n,istat)
      end if
      if(istat/=0)then
      x=-huge(1.0_dp)
      if(present(info))info=istat
      return
      end if
      do r=1,n
      x(r,:)=max_linear(a_sim,z(r,:))
      end do
      if(present(info))info=0
   end function simulate_conditional_max_linear

   function maxlinear_design_1d(coord,grid,variance,pixel_area) result(a)
      ! Exact numerical kernel of upstream maxLinDsgnMat for the
      ! one-dimensional discretized Smith model.
      real(dp),intent(in)::coord(:),grid(:),variance,pixel_area
      real(dp)::a(size(coord),size(grid)),ivar,cst,d
      integer::i,j
      if(variance<=0.0_dp .or. pixel_area<=0.0_dp) then
         a=0.0_dp
         return
      end if
      ivar=1.0_dp/variance
      cst=pixel_area*sqrt(ivar/(2.0_dp*pi))
      do j=1,size(grid)
         do i=1,size(coord)
            d=coord(i)-grid(j)
            a(i,j)=exp(-0.5_dp*d*d*ivar)*cst
         end do
      end do
      where(a<=1.0e-8_dp)a=0.0_dp
   end function maxlinear_design_1d

   function maxlinear_design_2d(coord,grid,cov,pixel_area) result(a)
      ! Exact numerical kernel of upstream maxLinDsgnMat for the
      ! two-dimensional discretized Smith model.  coord/grid are
      ! (n,2)/(p,2); cov is the 2x2 storm covariance matrix.
      real(dp),intent(in)::coord(:,:),grid(:,:),cov(:,:),pixel_area
      real(dp)::a(size(coord,1),size(grid,1)),det,idet,cst,d1,d2,q
      integer::i,j
      if(size(coord,2)/=2 .or. size(grid,2)/=2 .or. size(cov,1)/=2 .or. size(cov,2)/=2) then
         a=0.0_dp
         return
      end if
      det=cov(1,1)*cov(2,2)-cov(1,2)*cov(2,1)
      if(det<=0.0_dp .or. pixel_area<=0.0_dp) then
         a=0.0_dp
         return
      end if
      idet=1.0_dp/det
      cst=pixel_area*sqrt(idet)/(2.0_dp*pi)
      do j=1,size(grid,1)
         do i=1,size(coord,1)
            d1=coord(i,1)-grid(j,1)
            d2=coord(i,2)-grid(j,2)
            q=(cov(2,2)*d1*d1-(cov(1,2)+cov(2,1))*d1*d2+cov(1,1)*d2*d2)*idet
            a(i,j)=exp(-0.5_dp*q)*cst
         end do
      end do
      where(a<=1.0e-8_dp)a=0.0_dp
   end function maxlinear_design_2d
end module spatialextremes_maxlinear
