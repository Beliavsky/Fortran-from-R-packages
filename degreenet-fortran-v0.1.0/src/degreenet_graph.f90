! SPDX-License-Identifier: GPL-3.0-or-later
! Based on 'statnet' project software (statnet.org).
module degreenet_graph
  use degreenet_kinds, only : dp
  use degreenet_simulation, only : sample_model
  use degreenet_models, only : MODEL_YULE
  implicit none
  private
  type, public :: edge_list
    integer :: nvertex=0, nedge=0
    integer, allocatable :: tail(:), head(:)
  end type edge_list
  public :: reed_molloy, yule_graph

contains
  subroutine reed_molloy(degree,g,ok)
    integer,intent(in)::degree(:)
    type(edge_list),intent(out)::g
    logical,intent(out),optional::ok
    integer::n,m,i,j,k,v,ne
    integer,allocatable::d(:),idx(:),tail(:),head(:)
    logical::good
    n=size(degree);good=.true.
    if(any(degree<0).or.any(degree>=n).or.mod(sum(degree),2)/=0)good=.false.
    allocate(d(n),idx(n));d=degree;idx=[(i,i=1,n)]
    m=sum(degree)/2;allocate(tail(max(1,m)),head(max(1,m)));ne=0
    do while(good.and.maxval(d)>0)
      call sort_desc(d,idx)
      k=d(1);v=idx(1)
      if(k>n-1.or.k<0)then;good=.false.;exit;end if
      d(1)=0
      do j=2,k+1
        if(j>n.or.d(j)<=0)then;good=.false.;exit;end if
        ne=ne+1;if(ne>m)then;good=.false.;exit;end if
        tail(ne)=v;head(ne)=idx(j);d(j)=d(j)-1
      end do
    end do
    g%nvertex=n
    if(good)then
      g%nedge=ne;allocate(g%tail(ne),g%head(ne));g%tail=tail(1:ne);g%head=head(1:ne)
    else
      g%nedge=0;allocate(g%tail(0),g%head(0))
    end if
    if(present(ok))ok=good
  contains
    subroutine sort_desc(a,ii)
      integer,intent(inout)::a(:),ii(:)
      integer::q,r,ta,ti
      do q=2,size(a)
        ta=a(q);ti=ii(q);r=q-1
        do while(r>=1)
          if(a(r)>=ta)exit
          a(r+1)=a(r);ii(r+1)=ii(r);r=r-1
        end do
        a(r+1)=ta;ii(r+1)=ti
      end do
    end subroutine sort_desc
  end subroutine reed_molloy

  subroutine yule_graph(n,rho,g,ok,maxdeg,maxit)
    integer,intent(in)::n
    real(dp),intent(in)::rho
    type(edge_list),intent(out)::g
    logical,intent(out),optional::ok
    integer,intent(in),optional::maxdeg,maxit
    integer::md,mit,it
    integer,allocatable::deg(:)
    logical::good
    md=n-1;if(present(maxdeg))md=min(maxdeg,n-1);mit=100;if(present(maxit))mit=maxit
    allocate(deg(n));good=.false.
    do it=1,mit
      call sample_model(MODEL_YULE,[rho],n,deg,1,md)
      if(mod(sum(deg),2)/=0)cycle
      call reed_molloy(deg,g,good)
      if(good)exit
    end do
    if(present(ok))ok=good
  end subroutine yule_graph
end module degreenet_graph
