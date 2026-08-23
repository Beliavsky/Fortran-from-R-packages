module directional_permutation
   use directional_kinds, only : dp
   use directional_inference, only : vmf_mle, vmf_mle_result
   use directional_tests, only : test_result
   use directional_geometry, only : rotation_matrix
   implicit none
   private
   public :: embed_perm, hcf_perm, het_perm
   public :: embed_boot, hcf_boot, het_boot
contains
   function embed_perm(x1,x2,b) result(res)
      real(dp),intent(in)::x1(:,:),x2(:,:);integer,intent(in),optional::b
      type(test_result)::res
      real(dp),allocatable::x(:,:);integer,allocatable::lab(:);integer::n1,n2,n,bb,i,ex
      real(dp)::obs,st
      n1=size(x1,1);n2=size(x2,1);n=n1+n2;bb=999;if(present(b))bb=b
      allocate(x(n,size(x1,2)),lab(n));x(1:n1,:)=x1;x(n1+1:n,:)=x2;lab=[(1,i=1,n1),(2,i=1,n2)]
      obs=embed_stat(x,lab,n1,n2);ex=0
      do i=1,bb;call shuffle_int(lab);st=embed_stat(x,lab,n1,n2);if(st>obs)ex=ex+1;end do
      res%statistic=obs;res%p_value=real(ex+1,dp)/real(bb+1,dp)
   end function

   function hcf_perm(x1,x2,b) result(res)
      real(dp),intent(in)::x1(:,:),x2(:,:);integer,intent(in),optional::b
      type(test_result)::res
      real(dp),allocatable::x(:,:);integer,allocatable::lab(:);integer::n1,n2,n,bb,i,ex
      real(dp)::obs,st
      n1=size(x1,1);n2=size(x2,1);n=n1+n2;bb=999;if(present(b))bb=b
      allocate(x(n,size(x1,2)),lab(n));x(1:n1,:)=x1;x(n1+1:n,:)=x2;lab=[(1,i=1,n1),(2,i=1,n2)]
      obs=hcf_stat(x,lab,n1,n2);ex=0
      do i=1,bb;call shuffle_int(lab);st=hcf_stat(x,lab,n1,n2);if(st>obs)ex=ex+1;end do
      res%statistic=obs;res%p_value=real(ex+1,dp)/real(bb+1,dp)
   end function

   function het_perm(x1,x2,b) result(res)
      real(dp),intent(in)::x1(:,:),x2(:,:);integer,intent(in),optional::b
      type(test_result)::res
      real(dp),allocatable::x(:,:);integer,allocatable::lab(:);integer::n1,n2,n,bb,i,ex
      real(dp)::obs,st
      n1=size(x1,1);n2=size(x2,1);n=n1+n2;bb=999;if(present(b))bb=b
      allocate(x(n,size(x1,2)),lab(n));x(1:n1,:)=x1;x(n1+1:n,:)=x2;lab=[(1,i=1,n1),(2,i=1,n2)]
      obs=het_stat(x,lab,n1,n2);ex=0
      do i=1,bb;call shuffle_int(lab);st=het_stat(x,lab,n1,n2);if(st>obs)ex=ex+1;end do
      res%statistic=2.0_dp*obs;res%p_value=real(ex+1,dp)/real(bb+1,dp)
   end function


   function embed_boot(x1,x2,b) result(res)
      real(dp),intent(in)::x1(:,:),x2(:,:);integer,intent(in),optional::b
      type(test_result)::res
      real(dp),allocatable::y1(:,:),y2(:,:),yb(:,:)
      integer,allocatable::lab(:);integer::n1,n2,n,bb,i,j,ex,idx
      real(dp)::m1(size(x1,2)),m2(size(x1,2)),m(size(x1,2)),r1,r2,r,st
      real(dp)::rot1(size(x1,2),size(x1,2)),rot2(size(x1,2),size(x1,2)),u
      n1=size(x1,1);n2=size(x2,1);n=n1+n2;bb=999;if(present(b))bb=b
      allocate(y1(n1,size(x1,2)),y2(n2,size(x1,2)),yb(n,size(x1,2)),lab(n))
      lab=[(1,i=1,n1),(2,i=1,n2)]
      m1=sum(x1,dim=1);r1=sqrt(sum(m1*m1));m1=m1/max(r1,tiny(1.0_dp))
      m2=sum(x2,dim=1);r2=sqrt(sum(m2*m2));m2=m2/max(r2,tiny(1.0_dp))
      m=sum(x1,dim=1)+sum(x2,dim=1);r=sqrt(sum(m*m));m=m/max(r,tiny(1.0_dp))
      rot1=rotation_matrix(m1,m);rot2=rotation_matrix(m2,m)
      y1=matmul(x1,rot1);y2=matmul(x2,rot2)
      yb(1:n1,:)=x1;yb(n1+1:n,:)=x2;res%statistic=embed_stat(yb,lab,n1,n2);ex=0
      do i=1,bb
         do j=1,n1;call random_number(u);idx=min(n1,1+int(u*n1));yb(j,:)=y1(idx,:);end do
         do j=1,n2;call random_number(u);idx=min(n2,1+int(u*n2));yb(n1+j,:)=y2(idx,:);end do
         st=embed_stat(yb,lab,n1,n2);if(st>res%statistic)ex=ex+1
      end do
      res%p_value=real(ex+1,dp)/real(bb+1,dp)
   end function

   function hcf_boot(x1,x2,b) result(res)
      real(dp),intent(in)::x1(:,:),x2(:,:);integer,intent(in),optional::b
      type(test_result)::res
      real(dp),allocatable::y1(:,:),y2(:,:),yb(:,:)
      integer,allocatable::lab(:);integer::n1,n2,n,bb,i,j,ex,idx
      real(dp)::m1(size(x1,2)),m2(size(x1,2)),m(size(x1,2)),r1,r2,r,st,u
      real(dp)::rot1(size(x1,2),size(x1,2)),rot2(size(x1,2),size(x1,2))
      n1=size(x1,1);n2=size(x2,1);n=n1+n2;bb=999;if(present(b))bb=b
      allocate(y1(n1,size(x1,2)),y2(n2,size(x1,2)),yb(n,size(x1,2)),lab(n));lab=[(1,i=1,n1),(2,i=1,n2)]
      m1=sum(x1,dim=1);r1=sqrt(sum(m1*m1));m1=m1/max(r1,tiny(1.0_dp));m2=sum(x2,dim=1);r2=sqrt(sum(m2*m2));m2=m2/max(r2,tiny(1.0_dp))
      m=sum(x1,dim=1)+sum(x2,dim=1);r=sqrt(sum(m*m));m=m/max(r,tiny(1.0_dp));rot1=rotation_matrix(m1,m);rot2=rotation_matrix(m2,m)
      y1=matmul(x1,rot1);y2=matmul(x2,rot2);yb(1:n1,:)=x1;yb(n1+1:n,:)=x2;res%statistic=hcf_stat(yb,lab,n1,n2);ex=0
      do i=1,bb
         do j=1,n1;call random_number(u);idx=min(n1,1+int(u*n1));yb(j,:)=y1(idx,:);end do
         do j=1,n2;call random_number(u);idx=min(n2,1+int(u*n2));yb(n1+j,:)=y2(idx,:);end do
         st=hcf_stat(yb,lab,n1,n2);if(st>res%statistic)ex=ex+1
      end do
      res%p_value=real(ex+1,dp)/real(bb+1,dp)
   end function

   function het_boot(x1,x2,b) result(res)
      real(dp),intent(in)::x1(:,:),x2(:,:);integer,intent(in),optional::b
      type(test_result)::res
      real(dp),allocatable::y1(:,:),y2(:,:),yb(:,:)
      integer,allocatable::lab(:);integer::n1,n2,n,bb,i,j,ex,idx
      real(dp)::m1(size(x1,2)),m2(size(x1,2)),m(size(x1,2)),r1,r2,r,st,u
      real(dp)::rot1(size(x1,2),size(x1,2)),rot2(size(x1,2),size(x1,2))
      n1=size(x1,1);n2=size(x2,1);n=n1+n2;bb=999;if(present(b))bb=b
      allocate(y1(n1,size(x1,2)),y2(n2,size(x1,2)),yb(n,size(x1,2)),lab(n));lab=[(1,i=1,n1),(2,i=1,n2)]
      m1=sum(x1,dim=1);r1=sqrt(sum(m1*m1));m1=m1/max(r1,tiny(1.0_dp));m2=sum(x2,dim=1);r2=sqrt(sum(m2*m2));m2=m2/max(r2,tiny(1.0_dp))
      m=sum(x1,dim=1)+sum(x2,dim=1);r=sqrt(sum(m*m));m=m/max(r,tiny(1.0_dp));rot1=rotation_matrix(m1,m);rot2=rotation_matrix(m2,m)
      y1=matmul(x1,rot1);y2=matmul(x2,rot2);yb(1:n1,:)=x1;yb(n1+1:n,:)=x2;res%statistic=2.0_dp*het_stat(yb,lab,n1,n2);ex=0
      do i=1,bb
         do j=1,n1;call random_number(u);idx=min(n1,1+int(u*n1));yb(j,:)=y1(idx,:);end do
         do j=1,n2;call random_number(u);idx=min(n2,1+int(u*n2));yb(n1+j,:)=y2(idx,:);end do
         st=2.0_dp*het_stat(yb,lab,n1,n2);if(st>res%statistic)ex=ex+1
      end do
      res%p_value=real(ex+1,dp)/real(bb+1,dp)
   end function
   pure real(dp) function embed_stat(x,lab,n1,n2) result(v)
      real(dp),intent(in)::x(:,:);integer,intent(in)::lab(:),n1,n2
      real(dp)::s1(size(x,2)),s2(size(x,2)),sa(size(x,2)),r1,r2,rb;integer::i,n
      s1=0.0_dp;s2=0.0_dp;do i=1,size(x,1);if(lab(i)==1)then;s1=s1+x(i,:);else;s2=s2+x(i,:);end if;end do
      n=n1+n2;r1=sqrt(sum((s1/real(n1,dp))**2));r2=sqrt(sum((s2/real(n2,dp))**2));sa=(s1+s2)/real(n,dp);rb=sqrt(sum(sa*sa))
      v=real(n-2,dp)*(real(n1,dp)*r1*r1+real(n2,dp)*r2*r2-real(n,dp)*rb*rb)/max(tiny(1.0_dp),real(n,dp)-real(n1,dp)*r1*r1-real(n2,dp)*r2*r2)
   end function

   pure real(dp) function hcf_stat(x,lab,n1,n2) result(v)
      real(dp),intent(in)::x(:,:);integer,intent(in)::lab(:),n1,n2
      real(dp)::s1(size(x,2)),s2(size(x,2)),r1,r2,r;integer::i,n
      s1=0.0_dp;s2=0.0_dp;do i=1,size(x,1);if(lab(i)==1)then;s1=s1+x(i,:);else;s2=s2+x(i,:);end if;end do
      n=n1+n2;r1=sqrt(sum(s1*s1));r2=sqrt(sum(s2*s2));r=sqrt(sum((s1+s2)**2))
      v=real(n-2,dp)*(r1+r2-r)/max(tiny(1.0_dp),real(n,dp)-r1-r2)
   end function

   real(dp) function het_stat(x,lab,n1,n2) result(v)
      real(dp),intent(in)::x(:,:);integer,intent(in)::lab(:),n1,n2
      real(dp)::g1(n1,size(x,2)),g2(n2,size(x,2)),m1(size(x,2)),m2(size(x,2)),tw(size(x,2));integer::i,j1,j2
      type(vmf_mle_result)::a,b
      j1=0;j2=0;do i=1,size(x,1);if(lab(i)==1)then;j1=j1+1;g1(j1,:)=x(i,:);else;j2=j2+1;g2(j2,:)=x(i,:);end if;end do
      a=vmf_mle(g1);b=vmf_mle(g2);m1=sum(g1,dim=1)/real(n1,dp);m2=sum(g2,dim=1)/real(n2,dp)
      tw=a%kappa*real(n1,dp)*m1+b%kappa*real(n2,dp)*m2
      v=a%kappa*real(n1,dp)*sqrt(sum(m1*m1))+b%kappa*real(n2,dp)*sqrt(sum(m2*m2))-sqrt(sum(tw*tw))
   end function

   subroutine shuffle_int(a)
      integer,intent(inout)::a(:);integer::i,j,t;real(dp)::u
      do i=size(a),2,-1;call random_number(u);j=1+int(u*real(i,dp));if(j>i)j=i;t=a(i);a(i)=a(j);a(j)=t;end do
   end subroutine
end module directional_permutation
