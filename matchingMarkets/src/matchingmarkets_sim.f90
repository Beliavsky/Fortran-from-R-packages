module matchingmarkets_sim
   use matchingmarkets_kinds, only : dp, i8
   use matchingmarkets_rng, only : rng_t
   use matchingmarkets_mechanisms, only : iaa
   use matchingmarkets_types, only : assignment_result_t
   implicit none
   private
   public :: one_sided_sim_t, two_sided_sim_t, stabsim, stabsim2
   type :: one_sided_sim_t
      integer, allocatable :: market(:), group(:)
      real(dp), allocatable :: x1(:), x2(:), outcome(:)
   end type
   type :: two_sided_sim_t
      integer, allocatable :: market(:), student(:), college(:), matched(:)
      real(dp), allocatable :: sx(:), cx(:), surplus(:), outcome(:)
   end type
contains
   function stabsim(nmarkets,group_size,groups_per_market,seed) result(out)
      integer,intent(in)::nmarkets,group_size,groups_per_market
      integer(i8),intent(in),optional::seed
      type(one_sided_sim_t)::out
      type(rng_t)::rng
      integer::n,i,m,g,k
      real(dp)::latent
      if(present(seed))call rng%seed(seed)
      n=nmarkets*group_size*groups_per_market
      allocate(out%market(n),out%group(n),out%x1(n),out%x2(n),out%outcome(n));k=0
      do m=1,nmarkets
         do g=1,groups_per_market
            latent=rng%normal()
            do i=1,group_size
               k=k+1;out%market(k)=m;out%group(k)=g
               out%x1(k)=rng%normal();out%x2(k)=rng%normal()
               out%outcome(k)=0.5_dp+0.7_dp*out%x1(k)-0.3_dp*out%x2(k)+0.6_dp*latent+0.4_dp*rng%normal()
            end do
         end do
      end do
   end function stabsim

   function stabsim2(nmarkets,nstudents,ncolleges,slots,seed) result(out)
      integer,intent(in)::nmarkets,nstudents,ncolleges,slots(:)
      integer(i8),intent(in),optional::seed
      type(two_sided_sim_t)::out
      type(rng_t)::rng
      integer::m,s,c,k,n
      real(dp),allocatable::su(:),cu(:)
      integer,allocatable::sp(:,:),cp(:,:)
      type(assignment_result_t)::mat
      if(size(slots)/=ncolleges)error stop 'stabsim2: slots mismatch'
      if(present(seed))call rng%seed(seed)
      n=nmarkets*nstudents*ncolleges
      allocate(out%market(n),out%student(n),out%college(n),out%matched(n),out%sx(n),out%cx(n),out%surplus(n),out%outcome(n))
      k=0
      do m=1,nmarkets
         allocate(su(nstudents),cu(ncolleges),sp(ncolleges,nstudents),cp(nstudents,ncolleges))
         do s=1,nstudents;su(s)=rng%normal();end do
         do c=1,ncolleges;cu(c)=rng%normal();end do
         do s=1,nstudents
            call rank_colleges(su(s),cu,rng,sp(:,s))
         end do
         do c=1,ncolleges
            call rank_students(cu(c),su,rng,cp(:,c))
         end do
         mat=iaa(sp,cp,slots,'deferred')
         do s=1,nstudents
            do c=1,ncolleges
               k=k+1;out%market(k)=m;out%student(k)=s;out%college(k)=c
               out%sx(k)=su(s);out%cx(k)=cu(c)
               out%surplus(k)=su(s)+cu(c)+0.3_dp*rng%normal()
               out%matched(k)=merge(1,0,mat%assignment(s)==c)
               out%outcome(k)=1.0_dp+0.5_dp*su(s)+0.4_dp*cu(c)+0.3_dp*rng%normal()
            end do
         end do
         deallocate(su,cu,sp,cp)
      end do
   contains
      subroutine rank_colleges(xs,xc,r,ord)
         real(dp),intent(in)::xs,xc(:);type(rng_t),intent(inout)::r;integer,intent(out)::ord(:)
         real(dp)::u(size(xc));integer::i,j,t
         do i=1,size(xc);u(i)=0.7_dp*xs+xc(i)+0.2_dp*r%normal();ord(i)=i;end do
         do i=1,size(xc)-1;do j=i+1,size(xc);if(u(ord(j))>u(ord(i)))then;t=ord(i);ord(i)=ord(j);ord(j)=t;end if;end do;end do
      end subroutine
      subroutine rank_students(xc,xs,r,ord)
         real(dp),intent(in)::xc,xs(:);type(rng_t),intent(inout)::r;integer,intent(out)::ord(:)
         real(dp)::u(size(xs));integer::i,j,t
         do i=1,size(xs);u(i)=0.7_dp*xc+xs(i)+0.2_dp*r%normal();ord(i)=i;end do
         do i=1,size(xs)-1;do j=i+1,size(xs);if(u(ord(j))>u(ord(i)))then;t=ord(i);ord(i)=ord(j);ord(j)=t;end if;end do;end do
      end subroutine
   end function stabsim2
end module matchingmarkets_sim
