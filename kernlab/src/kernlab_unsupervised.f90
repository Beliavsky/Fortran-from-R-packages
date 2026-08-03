! SPDX-License-Identifier: GPL-2.0-only
module kernlab_unsupervised
  use kernlab_kinds
  use kernlab_types
  use kernlab_kernels, only: kernel_matrix
  use kernlab_linalg
  use kernlab_core, only: inchol
  implicit none
  private
  public :: kpca, kpca_predict, kcca, kha, kha_predict, kfa, kfa_predict
  public :: kkmeans, kkmeans_from_kernel, specc, specc_from_kernel
  public :: ranking, ranking_from_kernel, csi

contains

  subroutine kpca(x,kernel,result,features,threshold)
    real(dp),intent(in)::x(:,:)
    type(kernel_spec),intent(in)::kernel
    type(kpca_result),intent(out)::result
    integer,intent(in),optional::features
    real(dp),intent(in),optional::threshold
    real(dp),allocatable::k(:,:),kc(:,:),eval(:),evec(:,:),rowm(:)
    real(dp)::gm,th
    integer::st,n,nf,i
    call kernel_matrix(kernel,x,k,st);result%status=st;if(st/=KL_SUCCESS)return
    call center_kernel(k,kc,rowm,gm);call symmetric_eigen(kc,eval,evec,st)
    result%status=st;if(st/=KL_SUCCESS)return
    n=size(x,1);th=1.0e-4_dp;if(present(threshold))th=threshold
    nf=count(eval>th*max(1.0_dp,eval(1)));if(present(features))then;if(features>0)nf=min(nf,features);end if
    nf=max(1,nf)
    allocate(result%train(size(x,1),size(x,2)),result%coefficients(n,nf),result%rotated(n,nf), &
             result%eigenvalues(nf),result%kernel_row_mean(n))
    result%train=x;result%kernel=kernel;result%eigenvalues=eval(1:nf);result%kernel_row_mean=rowm;result%kernel_grand_mean=gm
    do i=1,nf
      result%coefficients(:,i)=evec(:,i)/sqrt(max(eval(i),tiny(1.0_dp)))
    end do
    result%rotated=matmul(kc,result%coefficients);result%status=KL_SUCCESS
  end subroutine kpca

  subroutine kpca_predict(model,x,scores,status)
    type(kpca_result),intent(in)::model
    real(dp),intent(in)::x(:,:)
    real(dp),allocatable,intent(out)::scores(:,:)
    integer,intent(out)::status
    real(dp),allocatable::kn(:,:),kc(:,:)
    real(dp),allocatable::newmean(:)
    integer::i,j,n,m
    call kernel_matrix(model%kernel,x,kn,status,model%train);if(status/=KL_SUCCESS)then;allocate(scores(0,0));return;end if
    n=size(x,1);m=size(model%train,1);allocate(kc(n,m),newmean(n));newmean=sum(kn,dim=2)/real(m,dp)
    do i=1,n;do j=1,m
      kc(i,j)=kn(i,j)-newmean(i)-model%kernel_row_mean(j)+model%kernel_grand_mean
    end do;end do
    allocate(scores(n,size(model%coefficients,2)));scores=matmul(kc,model%coefficients)
  end subroutine kpca_predict

  subroutine kcca(x,y,kernel,result,gamma,ncomps)
    real(dp),intent(in)::x(:,:),y(:,:)
    type(kernel_spec),intent(in)::kernel
    type(kcca_result),intent(out)::result
    real(dp),intent(in),optional::gamma
    integer,intent(in),optional::ncomps
    real(dp),allocatable::kx(:,:),ky(:,:),lh(:,:),rh(:,:),bs(:,:),tmp(:,:),val(:),vec(:,:)
    real(dp)::g
    integer::n,nc,st
    result%status=KL_INVALID_ARGUMENT;if(size(x,1)/=size(y,1))return
    n=size(x,1);g=0.1_dp;if(present(gamma))g=gamma;nc=min(10,n);if(present(ncomps))nc=min(n,ncomps)
    call kernel_matrix(kernel,x,kx,st);if(st/=KL_SUCCESS)return
    call kernel_matrix(kernel,y,ky,st);if(st/=KL_SUCCESS)return
    allocate(lh(2*n,2*n),rh(2*n,2*n));lh=0.0_dp;rh=0.0_dp
    lh(1:n,n+1:2*n)=matmul(kx,ky);lh(n+1:2*n,1:n)=transpose(lh(1:n,n+1:2*n))
    rh(1:n,1:n)=matmul(kx+diagmat(n,g),kx)+diagmat(n,1.0e-6_dp)
    rh(n+1:2*n,n+1:2*n)=matmul(ky+diagmat(n,g),ky)+diagmat(n,1.0e-6_dp)
    call matrix_inverse_sqrt(rh,bs,st);if(st/=KL_SUCCESS)return
    tmp=matmul(bs,matmul(lh,bs));call symmetric_eigen(tmp,val,vec,st);if(st/=KL_SUCCESS)return
    vec=matmul(bs,vec)
    allocate(result%xtrain(size(x,1),size(x,2)), result%ytrain(size(y,1),size(y,2)), &
      result%xcoef(n,nc), result%ycoef(n,nc), result%correlations(nc))
    result%xtrain=x;result%ytrain=y;result%kernel=kernel
    result%xcoef=vec(1:n,1:nc);result%ycoef=vec(n+1:2*n,1:nc);result%correlations=val(1:nc);result%status=KL_SUCCESS
  contains
    pure function diagmat(n,v) result(a)
      integer,intent(in)::n;real(dp),intent(in)::v;real(dp)::a(n,n);integer::q
      a=0.0_dp;do q=1,n;a(q,q)=v;end do
    end function diagmat
  end subroutine kcca

  subroutine kha(x,kernel,result,features,eta,threshold,maxiter)
    real(dp),intent(in)::x(:,:)
    type(kernel_spec),intent(in)::kernel
    type(kpca_result),intent(out)::result
    integer,intent(in),optional::features,maxiter
    real(dp),intent(in),optional::eta,threshold
    ! Deterministic batch equivalent: the converged KHA subspace equals kernel PCA.
    if(present(eta))then
      if(eta<0.0_dp) continue
    end if
    if(present(maxiter))then
      if(maxiter<0) continue
    end if
    call kpca(x,kernel,result,features,threshold)
  end subroutine kha

  subroutine kha_predict(model,x,scores,status)
    type(kpca_result),intent(in)::model;real(dp),intent(in)::x(:,:)
    real(dp),allocatable,intent(out)::scores(:,:);integer,intent(out)::status
    call kpca_predict(model,x,scores,status)
  end subroutine kha_predict

  subroutine kfa(x,kernel,result,features,subset,normalize)
    real(dp),intent(in)::x(:,:)
    type(kernel_spec),intent(in)::kernel
    type(kernel_model),intent(out)::result
    integer,intent(in),optional::features,subset
    logical,intent(in),optional::normalize
    type(inchol_result)::ic
    real(dp),allocatable::kbb(:,:),eye(:,:),inv(:,:)
    integer::nf,ss,st,i
    if(present(normalize))then
      if(.not.normalize) continue
    end if
    ss=size(x,1)
    if(present(subset))ss=min(ss,subset)
    nf=ss
    if(present(features))then
      if(features>0)nf=min(ss,features)
    end if
    call inchol(x,kernel,ic,tol=1.0e-10_dp,maxiter=nf);result%status=ic%status;if(ic%rank<1)return
    allocate(result%train(ic%rank,size(x,2)));do i=1,ic%rank;result%train(i,:)=x(ic%pivots(i),:);end do
    call kernel_matrix(kernel,result%train,kbb,st);if(st/=KL_SUCCESS)return
    allocate(eye(ic%rank,ic%rank));eye=0.0_dp;do i=1,ic%rank;eye(i,i)=1.0_dp;end do
    call solve_linear(kbb+1.0e-10_dp*eye,eye,inv,st);if(st/=KL_SUCCESS)return
    allocate(result%coefficients(ic%rank,ic%rank),result%bias(ic%rank));result%coefficients=inv;result%bias=0.0_dp
    result%kernel=kernel;result%model_type=MODEL_REGRESSION;result%status=KL_SUCCESS
  end subroutine kfa

  subroutine kfa_predict(model,x,scores,status)
    type(kernel_model),intent(in)::model;real(dp),intent(in)::x(:,:)
    real(dp),allocatable,intent(out)::scores(:,:);integer,intent(out)::status
    real(dp),allocatable::k(:,:)
    call kernel_matrix(model%kernel,x,k,status,model%train);if(status/=KL_SUCCESS)then;allocate(scores(0,0));return;end if
    allocate(scores(size(x,1),size(model%coefficients,2)));scores=matmul(k,model%coefficients)
    scores=scores-spread(sum(scores,dim=1)/real(size(scores,1),dp),1,size(scores,1))
  end subroutine kfa_predict

  subroutine kkmeans(x,centers,kernel,result,maxiter,initial_labels)
    real(dp),intent(in)::x(:,:)
    integer,intent(in)::centers
    type(kernel_spec),intent(in)::kernel
    type(cluster_result),intent(out)::result
    integer,intent(in),optional::maxiter
    integer,intent(in),optional::initial_labels(:)
    real(dp),allocatable::k(:,:)
    integer::st
    call kernel_matrix(kernel,x,k,st);if(st/=KL_SUCCESS)then;result%status=st;return;end if
    call kkmeans_from_kernel(k,centers,result,maxiter,initial_labels)
    if(result%status==KL_SUCCESS)call input_centers(x,result)
  end subroutine kkmeans

  subroutine kkmeans_from_kernel(k,centers,result,maxiter,initial_labels)
    real(dp),intent(in)::k(:,:)
    integer,intent(in)::centers
    type(cluster_result),intent(out)::result
    integer,intent(in),optional::maxiter
    integer,intent(in),optional::initial_labels(:)
    integer::n,c,i,j,it,limit,best,count
    integer,allocatable::lab(:),old(:),cnt(:),seeds(:)
    real(dp),allocatable::dist(:,:),within(:),cluster_sum(:),mindist(:)
    result%status=KL_INVALID_ARGUMENT;n=size(k,1);if(size(k,2)/=n.or.centers<1.or.centers>n)return
    limit=200;if(present(maxiter))limit=maxiter
    allocate(lab(n),old(n),cnt(centers),dist(n,centers),within(centers),cluster_sum(centers))
    if(present(initial_labels))then
      if(size(initial_labels)/=n)return
      lab=initial_labels
    else
      allocate(seeds(centers),mindist(n))
      seeds(1)=1
      mindist=huge(1.0_dp)
      do c=2,centers
        do i=1,n
          mindist(i)=min(mindist(i), k(i,i)+k(seeds(c-1),seeds(c-1)) &
            -2.0_dp*k(i,seeds(c-1)))
        end do
        seeds(c)=maxloc(mindist,dim=1)
      end do
      do i=1,n
        best=1
        mindist(i)=k(i,i)+k(seeds(1),seeds(1))-2.0_dp*k(i,seeds(1))
        do c=2,centers
          if(k(i,i)+k(seeds(c),seeds(c))-2.0_dp*k(i,seeds(c)) &
             <mindist(i))then
            best=c
            mindist(i)=k(i,i)+k(seeds(c),seeds(c))-2.0_dp*k(i,seeds(c))
          end if
        end do
        lab(i)=best
      end do
    end if
    do it=1,limit
      old=lab;cnt=0;cluster_sum=0.0_dp
      do c=1,centers
        cnt(c)=count(lab==c)
        if(cnt(c)>0)then
          do i=1,n;do j=1,n;if(lab(i)==c.and.lab(j)==c)cluster_sum(c)=cluster_sum(c)+k(i,j);end do;end do
          cluster_sum(c)=cluster_sum(c)/real(cnt(c)*cnt(c),dp)
        end if
      end do
      do i=1,n
        do c=1,centers
          if(cnt(c)>0)then
            dist(i,c)=k(i,i)-2.0_dp*sum(k(i,:),mask=lab==c)/real(cnt(c),dp)+cluster_sum(c)
          else
            dist(i,c)=huge(1.0_dp)
          end if
        end do
        best=minloc(dist(i,:),dim=1);lab(i)=best
      end do
      if(all(lab==old))exit
    end do
    within=0.0_dp
    do i=1,n;within(lab(i))=within(lab(i))+dist(i,lab(i));end do
    allocate(result%labels(n),result%withinss(centers));result%labels=lab;result%withinss=within
    result%iterations=min(it,limit);result%status=KL_SUCCESS
  end subroutine kkmeans_from_kernel

  subroutine input_centers(x,result)
    real(dp),intent(in)::x(:,:);type(cluster_result),intent(inout)::result
    integer::c,i,k,n,p,countc
    n=size(x,1);p=size(x,2);k=maxval(result%labels);allocate(result%centers(k,p));result%centers=0.0_dp
    do c=1,k
      countc=count(result%labels==c)
      do i=1,n;if(result%labels(i)==c)result%centers(c,:)=result%centers(c,:)+x(i,:);end do
      if(countc>0)result%centers(c,:)=result%centers(c,:)/real(countc,dp)
    end do
  end subroutine input_centers

  subroutine specc(x,centers,kernel,result,iterations)
    real(dp),intent(in)::x(:,:);integer,intent(in)::centers;type(kernel_spec),intent(in)::kernel
    type(cluster_result),intent(out)::result;integer,intent(in),optional::iterations
    real(dp),allocatable::k(:,:);integer::st
    call kernel_matrix(kernel,x,k,st);if(st/=KL_SUCCESS)then;result%status=st;return;end if
    call specc_from_kernel(k,centers,result,iterations)
    if(result%status==KL_SUCCESS)call input_centers(x,result)
  end subroutine specc

  subroutine specc_from_kernel(k,centers,result,iterations)
    real(dp),intent(in)::k(:,:);integer,intent(in)::centers;type(cluster_result),intent(out)::result
    integer,intent(in),optional::iterations
    real(dp),allocatable::a(:,:),deg(:),val(:),vec(:,:),emb(:,:),cen(:,:),within(:)
    integer,allocatable::lab(:)
    integer::n,i,j,st,it
    n=size(k,1);result%status=KL_INVALID_ARGUMENT;if(size(k,2)/=n.or.centers<1.or.centers>n)return
    allocate(a(n,n),deg(n));a=k;do i=1,n;a(i,i)=0.0_dp;end do;deg=sum(a,dim=2)
    do i=1,n;do j=1,n
      if(deg(i)>0.0_dp.and.deg(j)>0.0_dp)a(i,j)=a(i,j)/sqrt(deg(i)*deg(j))
    end do;end do
    call symmetric_eigen(a,val,vec,st);if(st/=KL_SUCCESS)then;result%status=st;return;end if
    allocate(emb(n,centers));emb=vec(:,1:centers)
    do i=1,n;if(vec_norm(emb(i,:))>0.0_dp)emb(i,:)=emb(i,:)/vec_norm(emb(i,:));end do
    call kmeans_dense(emb,centers,lab,cen,within,it,st,iterations)
    allocate(result%labels(n),result%embedding(n,centers),result%withinss(centers))
    result%labels=lab;result%embedding=emb;result%withinss=within
    result%iterations=it;result%status=st
  end subroutine specc_from_kernel

  subroutine ranking(x,y,kernel,result,alpha,iterations,edgegraph)
    real(dp),intent(in)::x(:,:),y(:);type(kernel_spec),intent(in)::kernel;type(ranking_result),intent(out)::result
    real(dp),intent(in),optional::alpha;integer,intent(in),optional::iterations;logical,intent(in),optional::edgegraph
    real(dp),allocatable::k(:,:);integer::st
    if(present(edgegraph))then
      if(edgegraph .neqv. edgegraph) continue
    end if
    call kernel_matrix(kernel,x,k,st)
    if(st/=KL_SUCCESS)then;result%status=st;return;end if
    call ranking_from_kernel(k,y,result,alpha,iterations)
  end subroutine ranking

  subroutine ranking_from_kernel(kin,y,result,alpha,iterations)
    real(dp),intent(in)::kin(:,:),y(:);type(ranking_result),intent(out)::result
    real(dp),intent(in),optional::alpha;integer,intent(in),optional::iterations
    real(dp),allocatable::k(:,:),d(:),base(:),score(:),newscore(:),ranks(:)
    logical,allocatable::labelled(:)
    real(dp)::a
    integer::n,i,j,it,limit
    result%status=KL_INVALID_ARGUMENT;n=size(kin,1);if(size(kin,2)/=n.or.size(y)/=n)return
    a=0.99_dp;if(present(alpha))a=alpha;limit=600;if(present(iterations))limit=iterations
    allocate(k(n,n),d(n),base(n),score(n),newscore(n),labelled(n));k=kin;do i=1,n;k(i,i)=0.0_dp;end do
    d=sum(k,dim=2);do i=1,n;do j=1,n;if(d(i)>0.0_dp.and.d(j)>0.0_dp)k(i,j)=k(i,j)/sqrt(d(i)*d(j));end do;end do
    labelled=abs(y)>tiny(1.0_dp);if(.not.any(labelled))return;base=matmul(k,merge(1.0_dp,0.0_dp,labelled));score=base
    allocate(result%convergence(n,limit));result%convergence(:,1)=score
    do it=2,limit
      newscore=base+a*matmul(k,merge(0.0_dp,score,.not.labelled));score=newscore;result%convergence(:,it)=score
    end do
    call vector_ranks(-score,ranks)
    allocate(result%score(n),result%rank(n));result%score=score;result%rank=ranks
    result%status=KL_SUCCESS
  end subroutine ranking_from_kernel

  subroutine csi(x,y,kernel,rank,result,centering,kappa,delta,tol)
    real(dp),intent(in)::x(:,:),y(:,:);type(kernel_spec),intent(in)::kernel;integer,intent(in)::rank
    type(csi_result),intent(out)::result;logical,intent(in),optional::centering
    real(dp),intent(in),optional::kappa,tol;integer,intent(in),optional::delta
    type(inchol_result)::ic
    real(dp),allocatable::g(:,:),q(:,:),r(:,:),v(:)
    real(dp)::threshold,nrm
    integer::n,m,i,j
    logical :: do_center
    threshold=1.0e-5_dp;if(present(tol))threshold=tol
    if(present(kappa))then
      if(kappa<0.0_dp) continue
    end if
    if(present(delta))then
      if(delta<0) continue
    end if
    call inchol(x,kernel,ic,threshold,rank)
    result%status=ic%status
    if(ic%rank<1)return
    n=size(x,1);m=ic%rank;allocate(g(n,m));g=ic%factor
    do_center=.true.
    if(present(centering))do_center=centering
    if(do_center)g=g-spread(sum(g,dim=1)/real(n,dp),1,n)
    allocate(q(n,m),r(m,m));q=0.0_dp;r=0.0_dp
    do j=1,m
      v=g(:,j)
      do i=1,j-1;r(i,j)=dot_product(q(:,i),v);v=v-r(i,j)*q(:,i);end do
      nrm=vec_norm(v);if(nrm>tiny(1.0_dp))then;r(j,j)=nrm;q(:,j)=v/nrm;end if
    end do
    allocate(result%g(n,m),result%q(n,m),result%r(m,m),result%pivots(m),result%predicted_gain(m),result%true_gain(m))
    result%g=ic%factor;result%q=q;result%r=r;result%pivots=ic%pivots;result%predicted_gain=ic%max_residuals
    do i=1,m;result%true_gain(i)=sum(matmul(transpose(q(:,1:i)),y)**2)/max(sum(y*y),tiny(1.0_dp));end do
    result%rank=m;result%status=KL_SUCCESS
  end subroutine csi

end module kernlab_unsupervised
