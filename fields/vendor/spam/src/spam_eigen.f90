module spam_eigen
use spam_kinds, only: dp
use spam_types, only: csr_matrix, eigen_result
implicit none
private
public :: spam_eigen_symmetric, spam_eigen_general
interface
 subroutine ds_eigen_f(maxnev,ncv,maxitr,n,iwhich,na,a,ja,ia,v,d,iparam)
  integer maxnev,ncv,maxitr,n,iwhich,na,ja(*),ia(na+1),iparam(8)
  double precision a(*),v(n,ncv),d(maxnev)
 end subroutine
 subroutine dn_eigen_f(maxnev,ncv,maxitr,n,iwhich,na,a,ja,ia,v,dr,di,iparam)
  integer maxnev,ncv,maxitr,n,iwhich,na,ja(*),ia(na+1),iparam(8)
  double precision a(*),v(n,ncv),dr(maxnev+1),di(maxnev+1)
 end subroutine
end interface
contains
function spam_eigen_symmetric(a,nev,which,ncv,maxitr) result(res)
type(csr_matrix),intent(in)::a
integer,intent(in)::nev
character(len=*),intent(in),optional::which
integer,intent(in),optional::ncv,maxitr
type(eigen_result)::res
integer::nc,mi,iw,n,iparam(8),nv
real(dp),allocatable::v(:,:),d(:)
if(a%nrow/=a%ncol) error stop 'spam_eigen_symmetric: square matrix required'
n=a%nrow;if(nev<1 .or. nev>=n)error stop 'spam_eigen_symmetric: require 1 <= nev < n'
nc=min(n,max(nev+2,2*nev+1));if(present(ncv))nc=min(n,max(nev+1,ncv))
mi=300;if(present(maxitr))mi=maxitr;iw=which_sym(which)
allocate(v(n,nc),d(nev));v=0;d=0;iparam=0
call ds_eigen_f(nev,nc,mi,n,iw,n,a%entries,a%colindices,a%rowpointers,v,d,iparam)
nv=max(0,min(nev,iparam(5)))
allocate(res%values(nv),res%vectors(n,nv));if(nv>0)then;res%values=d(:nv);res%vectors=v(:,:nv);end if
res%nconv=nv;res%niter=iparam(3);res%info=merge(0,1,nv==nev)
end function

function spam_eigen_general(a,nev,which,ncv,maxitr) result(res)
type(csr_matrix),intent(in)::a
integer,intent(in)::nev
character(len=*),intent(in),optional::which
integer,intent(in),optional::ncv,maxitr
type(eigen_result)::res
integer::nc,mi,iw,n,iparam(8),nv
real(dp),allocatable::v(:,:),dr(:),di(:)
if(a%nrow/=a%ncol) error stop 'spam_eigen_general: square matrix required'
n=a%nrow;if(nev<1 .or. nev>=n)error stop 'spam_eigen_general: require 1 <= nev < n'
nc=min(n,max(nev+2,2*nev+1));if(present(ncv))nc=min(n,max(nev+1,ncv))
mi=300;if(present(maxitr))mi=maxitr;iw=which_gen(which)
allocate(v(n,nc),dr(nev+1),di(nev+1));v=0;dr=0;di=0;iparam=0
call dn_eigen_f(nev,nc,mi,n,iw,n,a%entries,a%colindices,a%rowpointers,v,dr,di,iparam)
nv=max(0,min(nev,iparam(5)))
allocate(res%values(nv),res%imag_values(nv),res%vectors(n,nv))
if(nv>0)then;res%values=dr(:nv);res%imag_values=di(:nv);res%vectors=v(:,:nv);end if
res%nconv=nv;res%niter=iparam(3);res%info=merge(0,1,nv==nev)
end function

integer function which_sym(w) result(i)
character(len=*),intent(in),optional::w;character(len=2)::s
s='LM';if(present(w))s=upper2(w)
select case(s);case('LM');i=1;case('SM');i=2;case('LA');i=7;case('SA');i=8;case('BE');i=9;case default;i=1;end select
end function
integer function which_gen(w) result(i)
character(len=*),intent(in),optional::w;character(len=2)::s
s='LM';if(present(w))s=upper2(w)
select case(s);case('LM');i=1;case('SM');i=2;case('LR');i=3;case('SR');i=4;case('LI');i=5;case('SI');i=6;case default;i=1;end select
end function
pure function upper2(s) result(t)
character(len=*),intent(in)::s;character(len=2)::t;integer::j,c
t='  ';t(1:min(2,len_trim(s)))=s(1:min(2,len_trim(s)))
do j=1,2;c=iachar(t(j:j));if(c>=97.and.c<=122)t(j:j)=achar(c-32);end do
end function
end module spam_eigen
