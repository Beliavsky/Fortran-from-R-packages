program test_v04_features
    use ecos
    use ecos_sparse, only : sparse_triplet_builder, sparse_ldl_factor, triplet_to_csc
    implicit none
    integer, parameter :: ncache=20, nstar=25
    type(ecos_dims) :: d
    type(ecos_problem) :: p
    type(ecos_workspace) :: ws
    type(ecos_result) :: r1,r2,r3,r
    type(ecos_csc_matrix) :: gc,astar
    type(sparse_triplet_builder) :: tb
    type(sparse_ldl_factor) :: famd,fid
    real(dp) :: gd(ncache,ncache),c(ncache),c2(ncache),h(ncache)
    real(dp) :: illg(2,2),illc(2),illh(2)
    real(dp) :: g1(1,1),h1(1),cu(1),g2(2,1),h2(2),ci(1)
    real(dp) :: gbb(2,2),cbb(2),hbb(2)
    real(dp) :: gexp(3,1),aexp(1,1),cexp(1),hexp(3),bexp(1)
    type(ecos_csc_matrix) :: aec
    integer :: i,ie,info,bv(1)

    ! Persistent symbolic and warm-start reuse.  A value-only matrix update
    ! keeps symbolic analysis but invalidates the scaled warm point.
    gd=0.0_dp
    do i=1,ncache
        gd(i,i)=-1.0_dp
    end do
    c=1.0_dp; h=-1.0_dp; d%l=ncache; d%e=0
    if (allocated(d%q)) deallocate(d%q)
    call make_csc_matrix(gd,gc)
    call setup_problem_csc(p,c,gc,h,d,ierr=ie)
    call check(ie==0,'workspace setup')
    call ecos_setup(ws,p)
    call ecos_solve(ws,r1)
    call check(r1%exitflag==ECOS_OPTIMAL,'workspace first solve')
    call check(r1%symbolic_analyses==1 .and. r1%cached_symbolic_reuses==0,'first symbolic analysis')
    c2=1.0_dp; c2(1)=2.0_dp
    call ecos_update(ws,c=c2,ierr=ie)
    call ecos_solve(ws,r2)
    call check(r2%exitflag==ECOS_OPTIMAL,'workspace objective update')
    call check(r2%symbolic_analyses==0 .and. r2%cached_symbolic_reuses>=1,'symbolic cache reuse')
    call check(r2%cached_warm_starts>=1,'warm-start reuse')
    gd(1,1)=-2.0_dp
    call make_csc_matrix(gd,gc)
    call ecos_update(ws,g_csc=gc,ierr=ie)
    call ecos_solve(ws,r3)
    call check(r3%exitflag==ECOS_OPTIMAL,'workspace value-only matrix update')
    call check(r3%symbolic_analyses==0 .and. r3%cached_symbolic_reuses>=1,'matrix-value symbolic reuse')
    call check(r3%cached_warm_starts==0,'matrix-value warm invalidation')

    ! Cone-preserving equilibration on an extremely ill-scaled diagonal LP.
    illg=0.0_dp; illg(1,1)=-1.0e-8_dp; illg(2,2)=-1.0e8_dp
    illh=[-1.0e-8_dp,-1.0e8_dp]; illc=[1.0_dp,1.0_dp]; d%l=2
    call make_csc_matrix(illg,gc)
    call ecos_csolve(illc,gc,illh,d,r,control=ecos_control(maxit=100,equilibrate=.true.))
    call check(r%exitflag==ECOS_OPTIMAL,'equilibrated LP status')
    call check(maxval(abs(r%x-1.0_dp))<2.0e-7_dp,'equilibrated LP solution')
    call check(r%max_col_scale/r%min_col_scale>1.0e5_dp,'nontrivial column equilibration')
    call check(r%max_row_scale/r%min_row_scale>1.0e5_dp,'nontrivial row equilibration')

    ! AMD-style minimum degree avoids the catastrophic identity fill of a
    ! star graph whose center is numbered first.
    call tb%init(nstar,nstar,2*nstar)
    do i=1,nstar
        call tb%add(i,i,2.0_dp)
    end do
    do i=2,nstar
        call tb%add(1,i,-0.1_dp)
    end do
    call triplet_to_csc(tb,astar,.true.)
    call famd%analyze(astar,.false.,info,.true.)
    call check(info==0,'AMD analysis')
    call fid%analyze(astar,.false.,info,.false.)
    call check(info==0,'identity analysis')
    call check(famd%symbolic_nnz<fid%symbolic_nnz/4,'AMD fill reduction')

    ! Sparse primal-unbounded certificate.
    d%l=1; g1(1,1)=-1.0_dp; h1=0.0_dp; cu=-1.0_dp
    call make_csc_matrix(g1,gc)
    call ecos_csolve(cu,gc,h1,d,r,control=ecos_control(maxit=30))
    call check(r%exitflag==ECOS_DINF .and. r%dual_certificate_valid,'sparse unbounded certificate')
    call check(allocated(r%dual_certificate),'unbounded certificate vector')
    call check(dot_product(cu,r%dual_certificate)<-0.5_dp,'unbounded certificate normalization')

    ! Sparse primal-infeasible certificate.
    d%l=2; g2(:,1)=[-1.0_dp,1.0_dp]; h2=[-1.0_dp,0.0_dp]; ci=0.0_dp
    call make_csc_matrix(g2,gc)
    call ecos_csolve(ci,gc,h2,d,r,control=ecos_control(maxit=30))
    call check(r%exitflag==ECOS_PINF .and. r%primal_certificate_valid,'sparse infeasible certificate')
    call check(allocated(r%primal_certificate),'infeasible certificate vector')

    ! Infeasible exponential-cone case exercises the explicit dual-exp map
    ! (-v,-u,e*w) in K_exp used by the sparse certificate builder.
    gexp=0.0_dp; aexp=1.0_dp; cexp=0.0_dp
    hexp=[1.0_dp,0.0_dp,0.0_dp]; bexp=0.0_dp
    d%l=0; d%e=1
    call make_csc_matrix(gexp,gc); call make_csc_matrix(aexp,aec)
    call ecos_csolve(cexp,gc,hexp,d,r,aec,bexp, &
        control=ecos_control(maxit=20,certificate_maxit=120))
    call check(r%exitflag==ECOS_PINF .and. r%primal_certificate_valid, &
        'dual exponential certificate map')

    ! Sparse ECOS_BB path with fixed bound-row structure and symbolic reuse.
    d%l=2; d%e=0
    cbb=[-1.1_dp,-1.0_dp]
    gbb=reshape([2.0_dp,3.0_dp,1.0_dp,4.0_dp],[2,2])
    hbb=[4.0_dp,12.0_dp]; bv=[1]
    call make_csc_matrix(gbb,gc)
    call ecos_csolve(cbb,gc,hbb,d,r,bool_vars=bv)
    call check(r%exitflag==ECOS_OPTIMAL,'sparse BB status')
    call check(maxval(abs(r%x-[1.0_dp,2.0_dp]))<2.0e-6_dp,'sparse BB solution')
    call check(r%sparse_backend_used,'sparse BB retained sparse backend')
    call check(r%bb_symbolic_reuses>=1,'sparse BB symbolic reuse')

    ! ECOS inaccurate-optimal status offset when the best iterate meets the
    ! deliberately loose inaccurate tolerances but the exact tolerances do not.
    d%l=1; g1(1,1)=-1.0_dp; h1=-1.0_dp; cu=1.0_dp
    call make_csc_matrix(g1,gc)
    call ecos_csolve(cu,gc,h1,d,r,control=ecos_control(maxit=0, &
        feastol_inacc=1.0e3_dp,abstol_inacc=1.0e3_dp,reltol_inacc=1.0e3_dp))
    call check(r%exitflag==ECOS_OPTIMAL+ECOS_INACC_OFFSET,'inaccurate optimal status')

    print '(a)', 'PASS test_v04_features'
contains
    subroutine check(ok,msg)
        logical,intent(in)::ok
        character(*),intent(in)::msg
        if(.not.ok) then
            print '(a,1x,a)','FAIL',msg
            error stop 1
        end if
    end subroutine check
end program test_v04_features
