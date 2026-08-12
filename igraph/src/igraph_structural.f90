module igraph_structural
  use igraph_kinds, only : dp
  use igraph_graph, only : graph_t, are_adjacent, degree
  implicit none
  private
  public :: edge_density, reciprocity, triangle_count, local_triangle_count
  public :: transitivity_global, transitivity_local, assortativity_degree

contains

  real(dp) function edge_density(g, loops) result(x)
    type(graph_t),intent(in)::g
    logical,intent(in),optional::loops
    logical::lp
    real(dp)::den
    lp=.false.; if(present(loops)) lp=loops
    if(g%n<=1 .and. .not.lp) then; x=0.0_dp; return; end if
    if(g%directed) then
      den=real(g%n,dp)*real(g%n-1,dp); if(lp) den=real(g%n,dp)**2
    else
      den=real(g%n,dp)*real(g%n-1,dp)/2.0_dp; if(lp) den=den+real(g%n,dp)
    end if
    if(den>0.0_dp) then; x=real(g%m,dp)/den; else; x=0.0_dp; end if
  end function edge_density

  real(dp) function reciprocity(g) result(r)
    type(graph_t),intent(in)::g
    integer::k,u,v,mut,total
    if(.not.g%directed) then; r=1.0_dp; return; end if
    mut=0; total=0
    do k=1,g%m
      u=g%edge(1,k); v=g%edge(2,k); if(u==v) cycle
      total=total+1; if(are_adjacent(g,v,u)) mut=mut+1
    end do
    if(total==0) then; r=0.0_dp; else; r=real(mut,dp)/real(total,dp); end if
  end function reciprocity

  integer function triangle_count(g) result(t)
    type(graph_t),intent(in)::g
    integer::i,j,k
    t=0
    do i=1,g%n-2
      do j=i+1,g%n-1
        if(.not. und_adj(g,i,j)) cycle
        do k=j+1,g%n
          if(und_adj(g,i,k) .and. und_adj(g,j,k)) t=t+1
        end do
      end do
    end do
  contains
    logical function und_adj(h,u,v)
      type(graph_t),intent(in)::h; integer,intent(in)::u,v
      und_adj=are_adjacent(h,u,v)
      if(h%directed .and. .not.und_adj) und_adj=are_adjacent(h,v,u)
    end function und_adj
  end function triangle_count

  function local_triangle_count(g) result(t)
    type(graph_t),intent(in)::g
    integer,allocatable::t(:)
    integer::i,j,k
    allocate(t(g%n),source=0)
    do i=1,g%n-2
      do j=i+1,g%n-1
        if(.not.und_adj(g,i,j)) cycle
        do k=j+1,g%n
          if(und_adj(g,i,k) .and. und_adj(g,j,k)) then
            t(i)=t(i)+1; t(j)=t(j)+1; t(k)=t(k)+1
          end if
        end do
      end do
    end do
  contains
    logical function und_adj(h,u,v)
      type(graph_t),intent(in)::h; integer,intent(in)::u,v
      und_adj=are_adjacent(h,u,v)
      if(h%directed .and. .not.und_adj) und_adj=are_adjacent(h,v,u)
    end function und_adj
  end function local_triangle_count

  real(dp) function transitivity_global(g) result(c)
    type(graph_t),intent(in)::g
    integer,allocatable::d(:)
    integer::triples,t
    d=degree(g,'all',.false.); triples=sum(d*(d-1)/2); t=triangle_count(g)
    if(triples==0) then; c=0.0_dp; else; c=3.0_dp*real(t,dp)/real(triples,dp); end if
  end function transitivity_global

  function transitivity_local(g) result(c)
    type(graph_t),intent(in)::g
    real(dp),allocatable::c(:)
    integer,allocatable::d(:),t(:)
    integer::i,den
    d=degree(g,'all',.false.); t=local_triangle_count(g); allocate(c(g%n))
    do i=1,g%n
      den=d(i)*(d(i)-1)/2
      if(den>0) then; c(i)=real(t(i),dp)/real(den,dp); else; c(i)=0.0_dp; end if
    end do
  end function transitivity_local

  real(dp) function assortativity_degree(g) result(r)
    type(graph_t),intent(in)::g
    integer,allocatable::d(:)
    real(dp)::sx,sy,sxx,syy,sxy,m,dx,dy
    integer::k,u,v
    d=degree(g,'all',.false.); sx=0;sy=0;sxx=0;syy=0;sxy=0;m=0
    do k=1,g%m
      u=g%edge(1,k); v=g%edge(2,k); if(u==v) cycle
      dx=real(d(u),dp); dy=real(d(v),dp)
      sx=sx+dx; sy=sy+dy; sxx=sxx+dx*dx; syy=syy+dy*dy; sxy=sxy+dx*dy; m=m+1
      if(.not.g%directed) then
        sx=sx+dy; sy=sy+dx; sxx=sxx+dy*dy; syy=syy+dx*dx; sxy=sxy+dx*dy; m=m+1
      end if
    end do
    if(m<=1) then; r=0; return; end if
    sx=sx/m; sy=sy/m
    sxx=sxx/m-sx*sx; syy=syy/m-sy*sy; sxy=sxy/m-sx*sy
    if(sxx<=0 .or. syy<=0) then; r=0; else; r=sxy/sqrt(sxx*syy); end if
  end function assortativity_degree

end module igraph_structural
