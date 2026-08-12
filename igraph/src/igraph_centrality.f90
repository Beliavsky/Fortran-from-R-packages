module igraph_centrality
  use igraph_kinds, only : dp, igraph_inf
  use igraph_graph, only : graph_t, degree
  use igraph_traversal, only : distances_unweighted
  implicit none
  private
  public :: closeness_centrality, betweenness_centrality
  public :: pagerank, eigenvector_centrality

contains

  function closeness_centrality(g, normalized) result(c)
    type(graph_t),intent(in)::g
    logical,intent(in),optional::normalized
    real(dp),allocatable::c(:)
    integer,allocatable::d(:,:)
    logical::norm
    integer::i,j,reach
    real(dp)::s
    norm=.true.; if(present(normalized)) norm=normalized
    d=distances_unweighted(g,'out'); allocate(c(g%n),source=0.0_dp)
    do i=1,g%n
      s=0.0_dp; reach=0
      do j=1,g%n
        if(i==j) cycle
        if(d(i,j)>=0) then; s=s+real(d(i,j),dp); reach=reach+1; end if
      end do
      if(s>0.0_dp) then
        c(i)=real(reach,dp)/s
        if(norm .and. g%n>1) c(i)=c(i)*real(reach,dp)/real(g%n-1,dp)
      end if
    end do
  end function closeness_centrality

  function betweenness_centrality(g, normalized) result(cb)
    type(graph_t),intent(in)::g
    logical,intent(in),optional::normalized
    real(dp),allocatable::cb(:)
    integer,allocatable::q(:),stack(:),dist(:),pred_count(:),pred(:,:),pred_len(:)
    real(dp),allocatable::sigma(:),delta(:)
    integer::s,v,p,head,tail,sp,k,w,maxpred
    logical::norm
    norm=.false.; if(present(normalized)) norm=normalized
    maxpred=max(1,g%m*2+1)
    allocate(cb(g%n),q(g%n),stack(g%n),dist(g%n),sigma(g%n),delta(g%n))
    allocate(pred_count(g%n),pred_len(g%n),pred(maxpred,g%n)); cb=0.0_dp
    do s=1,g%n
      dist=-1; sigma=0.0_dp; delta=0.0_dp; pred_count=0; pred_len=0
      head=1; tail=1; q(1)=s; dist(s)=0; sigma(s)=1.0_dp; sp=0
      do while(head<=tail)
        v=q(head); head=head+1; sp=sp+1; stack(sp)=v
        do p=g%out_ptr(v),g%out_ptr(v+1)-1
          w=g%out_nei(p)
          if(dist(w)<0) then
            tail=tail+1; q(tail)=w; dist(w)=dist(v)+1
          end if
          if(dist(w)==dist(v)+1) then
            sigma(w)=sigma(w)+sigma(v)
            pred_len(w)=pred_len(w)+1
            if(pred_len(w)>maxpred) error stop 'betweenness: predecessor storage overflow'
            pred(pred_len(w),w)=v
          end if
        end do
      end do
      do k=sp,1,-1
        w=stack(k)
        do p=1,pred_len(w)
          v=pred(p,w)
          if(sigma(w)>0.0_dp) delta(v)=delta(v)+(sigma(v)/sigma(w))*(1.0_dp+delta(w))
        end do
        if(w/=s) cb(w)=cb(w)+delta(w)
      end do
    end do
    if(.not.g%directed) cb=cb/2.0_dp
    if(norm .and. g%n>2) then
      if(g%directed) then
        cb=cb/(real(g%n-1,dp)*real(g%n-2,dp))
      else
        cb=2.0_dp*cb/(real(g%n-1,dp)*real(g%n-2,dp))
      end if
    end if
  end function betweenness_centrality

  subroutine pagerank(g, score, damping, tol, maxiter, iterations)
    type(graph_t),intent(in)::g
    real(dp),allocatable,intent(out)::score(:)
    real(dp),intent(in),optional::damping,tol
    integer,intent(in),optional::maxiter
    integer,intent(out),optional::iterations
    real(dp)::damp,eps,dang,diff
    integer::mi,it,u,v,p,outdeg
    real(dp),allocatable::old(:),new(:)
    damp=0.85_dp; if(present(damping)) damp=damping
    eps=1.0e-12_dp; if(present(tol)) eps=tol
    mi=1000; if(present(maxiter)) mi=maxiter
    allocate(old(g%n),new(g%n),score(g%n))
    if(g%n==0) then; if(present(iterations)) iterations=0; return; end if
    old=1.0_dp/real(g%n,dp)
    do it=1,mi
      dang=0.0_dp
      do u=1,g%n
        outdeg=g%out_ptr(u+1)-g%out_ptr(u)
        if(outdeg==0) dang=dang+old(u)
      end do
      new=(1.0_dp-damp)/real(g%n,dp)+damp*dang/real(g%n,dp)
      do u=1,g%n
        outdeg=g%out_ptr(u+1)-g%out_ptr(u); if(outdeg==0) cycle
        do p=g%out_ptr(u),g%out_ptr(u+1)-1
          v=g%out_nei(p); new(v)=new(v)+damp*old(u)/real(outdeg,dp)
        end do
      end do
      diff=maxval(abs(new-old)); old=new
      if(diff<eps) exit
    end do
    score=old/sum(old); if(present(iterations)) iterations=it
  end subroutine pagerank

  subroutine eigenvector_centrality(g, score, value, tol, maxiter)
    type(graph_t),intent(in)::g
    real(dp),allocatable,intent(out)::score(:)
    real(dp),intent(out),optional::value
    real(dp),intent(in),optional::tol
    integer,intent(in),optional::maxiter
    real(dp)::eps,nrm,lambda,diff
    integer::mi,it,u,v,p
    real(dp),allocatable::x(:),y(:)
    eps=1.0e-12_dp; if(present(tol)) eps=tol
    mi=1000; if(present(maxiter)) mi=maxiter
    allocate(x(g%n),y(g%n),score(g%n))
    if(g%n==0) then; if(present(value)) value=0.0_dp; return; end if
    x=1.0_dp/sqrt(real(g%n,dp)); lambda=0.0_dp
    do it=1,mi
      y=0.0_dp
      if(g%directed) then
        do u=1,g%n
          do p=g%out_ptr(u),g%out_ptr(u+1)-1
            v=g%out_nei(p); y(v)=y(v)+g%out_w(p)*x(u)
          end do
        end do
      else
        do u=1,g%n
          do p=g%out_ptr(u),g%out_ptr(u+1)-1
            v=g%out_nei(p); y(u)=y(u)+g%out_w(p)*x(v)
          end do
        end do
      end if
      nrm=sqrt(sum(y*y)); if(nrm<=tiny(1.0_dp)) exit
      y=y/nrm; diff=maxval(abs(y-x)); x=y
      if(diff<eps) exit
    end do
    y=0.0_dp
    do u=1,g%n
      do p=g%out_ptr(u),g%out_ptr(u+1)-1
        v=g%out_nei(p)
        if(g%directed) then; y(v)=y(v)+g%out_w(p)*x(u)
        else; y(u)=y(u)+g%out_w(p)*x(v); end if
      end do
    end do
    lambda=dot_product(x,y); score=x
    if(maxval(abs(score))>0.0_dp) score=score/maxval(abs(score))
    if(present(value)) value=lambda
  end subroutine eigenvector_centrality

end module igraph_centrality
