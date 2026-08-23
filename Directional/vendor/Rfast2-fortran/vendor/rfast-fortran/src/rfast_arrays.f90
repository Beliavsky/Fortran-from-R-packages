module rfast_arrays
   use rfast_special, only : dp, nan_r
   implicit none
   private
   public :: mean_r, variance_r, median_r, mad_r, skewness_r, kurtosis_r, harmonic_mean, gini_r
   public :: sort_real, sort_integer, order_real, rank_average, nth_value, range_real
   public :: colmeans, rowmeans, colvars, rowvars, colsums, rowsums, colprods, rowprods
   public :: colmins, colmaxs, rowmins, rowmaxs, colmedians, rowmedians
   public :: colskewness, colkurtosis, rowcvs, colcvs
   public :: cumulative_sum, cumulative_prod, cumulative_min, cumulative_max
   public :: standardise_cols, standardise_vector, count_value_real, count_value_int
   public :: unique_real, unique_int, tabulate_int, binary_search_real
   public :: rep_row, rep_col, lower_tri_values, upper_tri_values
   public :: is_symmetric_matrix, is_integer_value, is_element_real
   public :: pmin_vec, pmax_vec, pmin_pmax_vec
   public :: colmads, rowmads, colranges, rowranges, colhameans, rowhameans
   public :: colsort, rowsort, colorder, roworder, colranks, rowranks, colnth, rownth
   public :: coltrue, rowtrue, colfalse, rowfalse, colall, rowall, colany, rowany
   public :: coltruefalse, rowtruefalse

contains

   pure real(dp) function mean_r(x) result(m)
      real(dp), intent(in) :: x(:)
      if (size(x)==0) then
         m=nan_r()
      else
         m=sum(x)/real(size(x),dp)
      end if
   end function mean_r

   pure real(dp) function variance_r(x, population) result(v)
      real(dp), intent(in) :: x(:)
      logical, intent(in), optional :: population
      real(dp) :: m, den
      logical :: pop
      pop=.false.; if(present(population))pop=population
      if(size(x)<2 .and. .not.pop) then
         v=nan_r(); return
      end if
      if(size(x)==0) then
         v=nan_r(); return
      end if
      m=mean_r(x)
      if(pop)then; den=real(size(x),dp); else; den=real(size(x)-1,dp); end if
      v=sum((x-m)**2)/den
   end function variance_r

   subroutine sort_real(x, decreasing)
      real(dp), intent(inout) :: x(:)
      logical, intent(in), optional :: decreasing
      logical :: dec
      dec=.false.;if(present(decreasing))dec=decreasing
      if(size(x)>1)call quicksort_real(x,1,size(x))
      if(dec)x=x(size(x):1:-1)
   end subroutine sort_real

   recursive subroutine quicksort_real(x,left,right)
      real(dp),intent(inout)::x(:);integer,intent(in)::left,right
      integer::i,j;real(dp)::pivot,tmp
      if(left>=right)return
      i=left;j=right;pivot=x(left+(right-left)/2)
      do
         if(i>j)exit
         do
            if(i>right)exit
            if(x(i)>=pivot)exit
            i=i+1
         end do
         do
            if(j<left)exit
            if(x(j)<=pivot)exit
            j=j-1
         end do
         if(i>j)exit
         tmp=x(i);x(i)=x(j);x(j)=tmp;i=i+1;j=j-1
      end do
      if(left<j)call quicksort_real(x,left,j)
      if(i<right)call quicksort_real(x,i,right)
   end subroutine quicksort_real

   subroutine sort_integer(x, decreasing)
      integer, intent(inout) :: x(:)
      logical, intent(in), optional :: decreasing
      logical :: dec
      dec=.false.;if(present(decreasing))dec=decreasing
      if(size(x)>1)call quicksort_integer(x,1,size(x))
      if(dec)x=x(size(x):1:-1)
   end subroutine sort_integer

   recursive subroutine quicksort_integer(x,left,right)
      integer,intent(inout)::x(:);integer,intent(in)::left,right
      integer::i,j,pivot,tmp
      if(left>=right)return
      i=left;j=right;pivot=x(left+(right-left)/2)
      do
         if(i>j)exit
         do
            if(i>right)exit
            if(x(i)>=pivot)exit
            i=i+1
         end do
         do
            if(j<left)exit
            if(x(j)<=pivot)exit
            j=j-1
         end do
         if(i>j)exit
         tmp=x(i);x(i)=x(j);x(j)=tmp;i=i+1;j=j-1
      end do
      if(left<j)call quicksort_integer(x,left,j)
      if(i<right)call quicksort_integer(x,i,right)
   end subroutine quicksort_integer

   function order_real(x, decreasing) result(idx)
      real(dp), intent(in) :: x(:)
      logical, intent(in), optional :: decreasing
      integer, allocatable :: idx(:)
      logical :: dec
      integer :: i
      dec=.false.;if(present(decreasing))dec=decreasing
      allocate(idx(size(x)));idx=[(i,i=1,size(x))]
      if(size(x)>1)call quicksort_index(x,idx,1,size(x))
      if(dec)idx=idx(size(idx):1:-1)
   end function order_real

   recursive subroutine quicksort_index(x,idx,left,right)
      real(dp),intent(in)::x(:);integer,intent(inout)::idx(:);integer,intent(in)::left,right
      integer::i,j,pivot,tmp
      if(left>=right)return
      i=left;j=right;pivot=idx(left+(right-left)/2)
      do
         if(i>j)exit
         do
            if(i>right)exit
            if(.not.index_less(x,idx(i),pivot))exit
            i=i+1
         end do
         do
            if(j<left)exit
            if(.not.index_less(x,pivot,idx(j)))exit
            j=j-1
         end do
         if(i>j)exit
         tmp=idx(i);idx(i)=idx(j);idx(j)=tmp;i=i+1;j=j-1
      end do
      if(left<j)call quicksort_index(x,idx,left,j)
      if(i<right)call quicksort_index(x,idx,i,right)
   end subroutine quicksort_index

   pure logical function index_less(x,i,j) result(less)
      real(dp),intent(in)::x(:);integer,intent(in)::i,j
      if(x(i)<x(j))then
         less=.true.
      else if(x(i)>x(j))then
         less=.false.
      else
         less=i<j
      end if
   end function index_less

   function rank_average(x) result(r)
      real(dp), intent(in) :: x(:)
      real(dp), allocatable :: r(:)
      integer, allocatable :: ord(:)
      integer :: i,j,k,n
      real(dp) :: av
      n=size(x); allocate(r(n)); ord=order_real(x)
      i=1
      do while(i<=n)
         j=i
         do while(j<n)
            if(x(ord(j+1))<x(ord(i)).or.x(ord(j+1))>x(ord(i)))exit
            j=j+1
         end do
         av=0.5_dp*real(i+j,dp)
         do k=i,j; r(ord(k))=av; end do
         i=j+1
      end do
   end function rank_average

   function median_r(x) result(m)
      real(dp), intent(in) :: x(:)
      real(dp) :: m
      real(dp), allocatable :: y(:)
      integer :: n
      n=size(x)
      if(n==0)then; m=nan_r(); return; end if
      y=x; call sort_real(y)
      if(mod(n,2)==1)then
         m=y((n+1)/2)
      else
         m=0.5_dp*(y(n/2)+y(n/2+1))
      end if
   end function median_r

   function nth_value(x, k) result(v)
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: k
      real(dp) :: v
      real(dp), allocatable :: y(:)
      if(k<1 .or. k>size(x))then; v=nan_r(); return; end if
      y=x; call sort_real(y); v=y(k)
   end function nth_value

   function mad_r(x, constant) result(v)
      real(dp), intent(in) :: x(:)
      real(dp), intent(in), optional :: constant
      real(dp) :: v,c,m
      c=1.482602218505602_dp; if(present(constant))c=constant
      m=median_r(x); v=c*median_r(abs(x-m))
   end function mad_r

   pure real(dp) function skewness_r(x, corrected) result(s)
      real(dp), intent(in) :: x(:)
      logical, intent(in), optional :: corrected
      real(dp) :: m,m2,m3,n
      logical :: corr
      corr=.false.; if(present(corrected))corr=corrected
      n=real(size(x),dp); m=mean_r(x)
      m2=sum((x-m)**2)/n; m3=sum((x-m)**3)/n
      if(m2<=0.0_dp)then; s=0.0_dp; return; end if
      s=m3/m2**1.5_dp
      if(corr .and. size(x)>2)s=sqrt(n*(n-1.0_dp))/(n-2.0_dp)*s
   end function skewness_r

   pure real(dp) function kurtosis_r(x, excess) result(k)
      real(dp), intent(in) :: x(:)
      logical, intent(in), optional :: excess
      logical :: ex
      real(dp) :: m,m2,m4,n
      ex=.true.; if(present(excess))ex=excess
      n=real(size(x),dp); m=mean_r(x)
      m2=sum((x-m)**2)/n; m4=sum((x-m)**4)/n
      if(m2<=0.0_dp)then; k=0.0_dp; return; end if
      k=m4/(m2*m2)
      if(ex)k=k-3.0_dp
   end function kurtosis_r

   pure real(dp) function harmonic_mean(x) result(h)
      real(dp), intent(in) :: x(:)
      if(any(abs(x)<=0.0_dp))then
         h=0.0_dp
      else
         h=real(size(x),dp)/sum(1.0_dp/x)
      end if
   end function harmonic_mean

   function gini_r(x) result(g)
      real(dp), intent(in) :: x(:)
      real(dp) :: g, sx
      real(dp), allocatable :: y(:)
      integer :: i,n
      n=size(x); if(n==0)then; g=nan_r(); return; end if
      y=x; call sort_real(y); sx=sum(y)
      if(abs(sx)<=0.0_dp)then; g=0.0_dp; return; end if
      g=0.0_dp
      do i=1,n; g=g+real(2*i-n-1,dp)*y(i); end do
      g=g/(real(n,dp)*sx)
   end function gini_r

   pure function range_real(x) result(r)
      real(dp), intent(in) :: x(:)
      real(dp) :: r(2)
      r=[minval(x),maxval(x)]
   end function range_real

   pure function colsums(x) result(v)
      real(dp), intent(in) :: x(:,:)
      real(dp) :: v(size(x,2))
      v=sum(x,dim=1)
   end function colsums

   pure function rowsums(x) result(v)
      real(dp), intent(in) :: x(:,:)
      real(dp) :: v(size(x,1))
      v=sum(x,dim=2)
   end function rowsums

   pure function colmeans(x) result(v)
      real(dp), intent(in) :: x(:,:)
      real(dp) :: v(size(x,2))
      v=sum(x,dim=1)/real(size(x,1),dp)
   end function colmeans

   pure function rowmeans(x) result(v)
      real(dp), intent(in) :: x(:,:)
      real(dp) :: v(size(x,1))
      v=sum(x,dim=2)/real(size(x,2),dp)
   end function rowmeans

   pure function colvars(x, population) result(v)
      real(dp), intent(in) :: x(:,:)
      logical, intent(in), optional :: population
      real(dp) :: v(size(x,2)),m(size(x,2)),den
      logical :: pop
      pop=.false.; if(present(population))pop=population
      m=colmeans(x)
      if(pop)then; den=real(size(x,1),dp); else; den=real(size(x,1)-1,dp); end if
      v=sum((x-spread(m,1,size(x,1)))**2,dim=1)/den
   end function colvars

   pure function rowvars(x, population) result(v)
      real(dp), intent(in) :: x(:,:)
      logical, intent(in), optional :: population
      real(dp) :: v(size(x,1)),m(size(x,1)),den
      logical :: pop
      pop=.false.; if(present(population))pop=population
      m=rowmeans(x)
      if(pop)then; den=real(size(x,2),dp); else; den=real(size(x,2)-1,dp); end if
      v=sum((x-spread(m,2,size(x,2)))**2,dim=2)/den
   end function rowvars

   pure function colprods(x) result(v)
      real(dp), intent(in) :: x(:,:)
      real(dp) :: v(size(x,2))
      integer :: j
      do j=1,size(x,2); v(j)=product(x(:,j)); end do
   end function colprods

   pure function rowprods(x) result(v)
      real(dp), intent(in) :: x(:,:)
      real(dp) :: v(size(x,1))
      integer :: i
      do i=1,size(x,1); v(i)=product(x(i,:)); end do
   end function rowprods

   pure function colmins(x) result(v)
      real(dp),intent(in)::x(:,:); real(dp)::v(size(x,2)); v=minval(x,dim=1)
   end function colmins
   pure function colmaxs(x) result(v)
      real(dp),intent(in)::x(:,:); real(dp)::v(size(x,2)); v=maxval(x,dim=1)
   end function colmaxs
   pure function rowmins(x) result(v)
      real(dp),intent(in)::x(:,:); real(dp)::v(size(x,1)); v=minval(x,dim=2)
   end function rowmins
   pure function rowmaxs(x) result(v)
      real(dp),intent(in)::x(:,:); real(dp)::v(size(x,1)); v=maxval(x,dim=2)
   end function rowmaxs

   function colmedians(x) result(v)
      real(dp),intent(in)::x(:,:); real(dp)::v(size(x,2)); integer::j
      do j=1,size(x,2); v(j)=median_r(x(:,j)); end do
   end function colmedians
   function rowmedians(x) result(v)
      real(dp),intent(in)::x(:,:); real(dp)::v(size(x,1)); integer::i
      do i=1,size(x,1); v(i)=median_r(x(i,:)); end do
   end function rowmedians

   pure function colskewness(x) result(v)
      real(dp),intent(in)::x(:,:); real(dp)::v(size(x,2)); integer::j
      do j=1,size(x,2); v(j)=skewness_r(x(:,j)); end do
   end function colskewness
   pure function colkurtosis(x) result(v)
      real(dp),intent(in)::x(:,:); real(dp)::v(size(x,2)); integer::j
      do j=1,size(x,2); v(j)=kurtosis_r(x(:,j)); end do
   end function colkurtosis

   pure function colcvs(x) result(v)
      real(dp),intent(in)::x(:,:); real(dp)::v(size(x,2)); real(dp)::m(size(x,2))
      m=colmeans(x); v=sqrt(colvars(x))/m
   end function colcvs
   pure function rowcvs(x) result(v)
      real(dp),intent(in)::x(:,:); real(dp)::v(size(x,1)); real(dp)::m(size(x,1))
      m=rowmeans(x); v=sqrt(rowvars(x))/m
   end function rowcvs

   pure function cumulative_sum(x) result(v)
      real(dp),intent(in)::x(:); real(dp)::v(size(x)); integer::i
      if(size(x)==0)return; v(1)=x(1); do i=2,size(x);v(i)=v(i-1)+x(i);end do
   end function cumulative_sum
   pure function cumulative_prod(x) result(v)
      real(dp),intent(in)::x(:); real(dp)::v(size(x)); integer::i
      if(size(x)==0)return; v(1)=x(1); do i=2,size(x);v(i)=v(i-1)*x(i);end do
   end function cumulative_prod
   pure function cumulative_min(x) result(v)
      real(dp),intent(in)::x(:); real(dp)::v(size(x)); integer::i
      if(size(x)==0)return; v(1)=x(1); do i=2,size(x);v(i)=min(v(i-1),x(i));end do
   end function cumulative_min
   pure function cumulative_max(x) result(v)
      real(dp),intent(in)::x(:); real(dp)::v(size(x)); integer::i
      if(size(x)==0)return; v(1)=x(1); do i=2,size(x);v(i)=max(v(i-1),x(i));end do
   end function cumulative_max

   pure function standardise_vector(x) result(z)
      real(dp),intent(in)::x(:); real(dp)::z(size(x)); real(dp)::m,s
      m=mean_r(x); s=sqrt(variance_r(x)); if(s>0)then; z=(x-m)/s; else; z=0; end if
   end function standardise_vector
   pure function standardise_cols(x) result(z)
      real(dp),intent(in)::x(:,:); real(dp)::z(size(x,1),size(x,2)); integer::j
      do j=1,size(x,2);z(:,j)=standardise_vector(x(:,j));end do
   end function standardise_cols

   pure integer function count_value_real(x,value) result(n)
      real(dp),intent(in)::x(:),value; n=count(x<=value.and.x>=value)
   end function count_value_real
   pure integer function count_value_int(x,value) result(n)
      integer,intent(in)::x(:),value; n=count(x<=value.and.x>=value)
   end function count_value_int

   function unique_real(x) result(u)
      real(dp),intent(in)::x(:); real(dp),allocatable::u(:),y(:); integer::i,n
      if(size(x)==0)then;allocate(u(0));return;end if
      y=x;call sort_real(y);n=1
      do i=2,size(y);if(y(i)<y(i-1).or.y(i)>y(i-1))n=n+1;end do
      allocate(u(n));u(1)=y(1);n=1
      do i=2,size(y);if(y(i)<y(i-1).or.y(i)>y(i-1))then;n=n+1;u(n)=y(i);end if;end do
   end function unique_real

   function unique_int(x) result(u)
      integer,intent(in)::x(:); integer,allocatable::u(:),y(:); integer::i,n
      if(size(x)==0)then;allocate(u(0));return;end if
      y=x;call sort_integer(y);n=1
      do i=2,size(y);if(y(i)<y(i-1).or.y(i)>y(i-1))n=n+1;end do
      allocate(u(n));u(1)=y(1);n=1
      do i=2,size(y);if(y(i)<y(i-1).or.y(i)>y(i-1))then;n=n+1;u(n)=y(i);end if;end do
   end function unique_int

   function tabulate_int(x, nbins) result(tab)
      integer,intent(in)::x(:),nbins; integer,allocatable::tab(:); integer::i
      allocate(tab(nbins));tab=0
      do i=1,size(x);if(x(i)>=1.and.x(i)<=nbins)tab(x(i))=tab(x(i))+1;end do
   end function tabulate_int

   integer function binary_search_real(x,value) result(pos)
      real(dp),intent(in)::x(:),value; integer::lo,hi,mid
      lo=1;hi=size(x);pos=0
      do while(lo<=hi)
         mid=(lo+hi)/2
         if(x(mid)<=value.and.x(mid)>=value)then;pos=mid;return
         else if(x(mid)<value)then;lo=mid+1
         else;hi=mid-1;end if
      end do
   end function binary_search_real

   function rep_row(x,n) result(a)
      real(dp),intent(in)::x(:);integer,intent(in)::n;real(dp)::a(n,size(x));integer::i
      do i=1,n;a(i,:)=x;end do
   end function rep_row
   function rep_col(x,n) result(a)
      real(dp),intent(in)::x(:);integer,intent(in)::n;real(dp)::a(size(x),n);integer::j
      do j=1,n;a(:,j)=x;end do
   end function rep_col

   function lower_tri_values(a,diag) result(v)
      real(dp),intent(in)::a(:,:);logical,intent(in),optional::diag
      real(dp),allocatable::v(:);logical::d;integer::i,j,k,n
      d=.false.;if(present(diag))d=diag;n=0
      do j=1,size(a,2);do i=1,size(a,1);if(i>j.or.(d.and.i==j))n=n+1;end do;end do
      allocate(v(n));k=0
      do j=1,size(a,2);do i=1,size(a,1);if(i>j.or.(d.and.i==j))then;k=k+1;v(k)=a(i,j);end if;end do;end do
   end function lower_tri_values
   function upper_tri_values(a,diag) result(v)
      real(dp),intent(in)::a(:,:);logical,intent(in),optional::diag
      real(dp),allocatable::v(:);logical::d;integer::i,j,k,n
      d=.false.;if(present(diag))d=diag;n=0
      do j=1,size(a,2);do i=1,size(a,1);if(i<j.or.(d.and.i==j))n=n+1;end do;end do
      allocate(v(n));k=0
      do j=1,size(a,2);do i=1,size(a,1);if(i<j.or.(d.and.i==j))then;k=k+1;v(k)=a(i,j);end if;end do;end do
   end function upper_tri_values

   pure logical function is_symmetric_matrix(a,tol) result(ok)
      real(dp),intent(in)::a(:,:);real(dp),intent(in),optional::tol;real(dp)::t
      t=1e-12_dp;if(present(tol))t=tol
      ok=size(a,1)==size(a,2);if(ok)ok=maxval(abs(a-transpose(a)))<=t
   end function is_symmetric_matrix
   pure elemental logical function is_integer_value(x,tol) result(ok)
      real(dp),intent(in)::x;real(dp),intent(in),optional::tol;real(dp)::t
      t=1e-12_dp;if(present(tol))t=tol;ok=abs(x-anint(x))<=t
   end function is_integer_value
   pure logical function is_element_real(x,value) result(ok)
      real(dp),intent(in)::x(:),value;ok=any(x<=value.and.x>=value)
   end function is_element_real

   pure function pmin_vec(x,y) result(z)
      real(dp),intent(in)::x(:),y(:);real(dp)::z(size(x));z=min(x,y)
   end function pmin_vec
   pure function pmax_vec(x,y) result(z)
      real(dp),intent(in)::x(:),y(:);real(dp)::z(size(x));z=max(x,y)
   end function pmax_vec
   pure subroutine pmin_pmax_vec(x,y,mn,mx)
      real(dp),intent(in)::x(:),y(:);real(dp),intent(out)::mn(size(x)),mx(size(x));mn=min(x,y);mx=max(x,y)
   end subroutine pmin_pmax_vec


   function colmads(x) result(v)
      real(dp),intent(in)::x(:,:);real(dp)::v(size(x,2));integer::j
      do j=1,size(x,2);v(j)=mad_r(x(:,j));end do
   end function colmads

   function rowmads(x) result(v)
      real(dp),intent(in)::x(:,:);real(dp)::v(size(x,1));integer::i
      do i=1,size(x,1);v(i)=mad_r(x(i,:));end do
   end function rowmads

   pure function colranges(x) result(v)
      real(dp),intent(in)::x(:,:);real(dp)::v(size(x,2))
      v=maxval(x,dim=1)-minval(x,dim=1)
   end function colranges

   pure function rowranges(x) result(v)
      real(dp),intent(in)::x(:,:);real(dp)::v(size(x,1))
      v=maxval(x,dim=2)-minval(x,dim=2)
   end function rowranges

   pure function colhameans(x) result(v)
      real(dp),intent(in)::x(:,:);real(dp)::v(size(x,2));integer::j
      do j=1,size(x,2);v(j)=harmonic_mean(x(:,j));end do
   end function colhameans

   pure function rowhameans(x) result(v)
      real(dp),intent(in)::x(:,:);real(dp)::v(size(x,1));integer::i
      do i=1,size(x,1);v(i)=harmonic_mean(x(i,:));end do
   end function rowhameans

   function colsort(x,decreasing) result(z)
      real(dp),intent(in)::x(:,:);logical,intent(in),optional::decreasing
      real(dp)::z(size(x,1),size(x,2));integer::j
      z=x;do j=1,size(x,2);call sort_real(z(:,j),decreasing);end do
   end function colsort

   function rowsort(x,decreasing) result(z)
      real(dp),intent(in)::x(:,:);logical,intent(in),optional::decreasing
      real(dp)::z(size(x,1),size(x,2));integer::i
      z=x;do i=1,size(x,1);call sort_real(z(i,:),decreasing);end do
   end function rowsort

   function colorder(x,decreasing) result(ind)
      real(dp),intent(in)::x(:,:);logical,intent(in),optional::decreasing
      integer::ind(size(x,1),size(x,2)),j
      do j=1,size(x,2);ind(:,j)=order_real(x(:,j),decreasing);end do
   end function colorder

   function roworder(x,decreasing) result(ind)
      real(dp),intent(in)::x(:,:);logical,intent(in),optional::decreasing
      integer::ind(size(x,1),size(x,2)),i
      do i=1,size(x,1);ind(i,:)=order_real(x(i,:),decreasing);end do
   end function roworder

   function colranks(x) result(r)
      real(dp),intent(in)::x(:,:);real(dp)::r(size(x,1),size(x,2));integer::j
      do j=1,size(x,2);r(:,j)=rank_average(x(:,j));end do
   end function colranks

   function rowranks(x) result(r)
      real(dp),intent(in)::x(:,:);real(dp)::r(size(x,1),size(x,2));integer::i
      do i=1,size(x,1);r(i,:)=rank_average(x(i,:));end do
   end function rowranks

   function colnth(x,k) result(v)
      real(dp),intent(in)::x(:,:);integer,intent(in)::k;real(dp)::v(size(x,2));integer::j
      do j=1,size(x,2);v(j)=nth_value(x(:,j),k);end do
   end function colnth

   function rownth(x,k) result(v)
      real(dp),intent(in)::x(:,:);integer,intent(in)::k;real(dp)::v(size(x,1));integer::i
      do i=1,size(x,1);v(i)=nth_value(x(i,:),k);end do
   end function rownth

   pure function coltrue(x) result(v)
      logical,intent(in)::x(:,:);integer::v(size(x,2));v=count(x,dim=1)
   end function coltrue
   pure function rowtrue(x) result(v)
      logical,intent(in)::x(:,:);integer::v(size(x,1));v=count(x,dim=2)
   end function rowtrue
   pure function colfalse(x) result(v)
      logical,intent(in)::x(:,:);integer::v(size(x,2));v=size(x,1)-count(x,dim=1)
   end function colfalse
   pure function rowfalse(x) result(v)
      logical,intent(in)::x(:,:);integer::v(size(x,1));v=size(x,2)-count(x,dim=2)
   end function rowfalse
   pure function colall(x) result(v)
      logical,intent(in)::x(:,:);logical::v(size(x,2));v=all(x,dim=1)
   end function colall
   pure function rowall(x) result(v)
      logical,intent(in)::x(:,:);logical::v(size(x,1));v=all(x,dim=2)
   end function rowall
   pure function colany(x) result(v)
      logical,intent(in)::x(:,:);logical::v(size(x,2));v=any(x,dim=1)
   end function colany
   pure function rowany(x) result(v)
      logical,intent(in)::x(:,:);logical::v(size(x,1));v=any(x,dim=2)
   end function rowany
   pure function coltruefalse(x) result(v)
      logical,intent(in)::x(:,:);integer::v(2,size(x,2))
      v(2,:)=count(x,dim=1);v(1,:)=size(x,1)-v(2,:)
   end function coltruefalse
   pure function rowtruefalse(x) result(v)
      logical,intent(in)::x(:,:);integer::v(2,size(x,1))
      v(2,:)=count(x,dim=2);v(1,:)=size(x,2)-v(2,:)
   end function rowtruefalse

end module rfast_arrays
