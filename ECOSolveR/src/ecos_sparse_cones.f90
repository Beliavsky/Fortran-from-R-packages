! SPDX-License-Identifier: GPL-3.0-or-later
module ecos_sparse_cones
    use ecos_types, only : dp, ecos_problem, ecos_csr_matrix, ecos_csc_matrix
    use ecos_sparse, only : sparse_triplet_builder, triplet_to_csr, triplet_to_csc, csr_matvec
    implicit none
    private
    public :: sparse_cone_linearize, sparse_cone_values
    public :: sparse_cone_slack, sparse_limit_exp_step

contains

    subroutine sparse_cone_slack(prob,x,s)
        type(ecos_problem), intent(in) :: prob
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: s(:)
        real(dp), allocatable :: gx(:)
        if(size(s)==0) return
        allocate(gx(size(s)))
        call csr_matvec(prob%g_csr,x,gx)
        s=prob%h-gx
    end subroutine sparse_cone_slack

    subroutine add_scaled_row(tb,outrow,g,row,scale)
        type(sparse_triplet_builder), intent(inout) :: tb
        integer, intent(in) :: outrow,row
        type(ecos_csr_matrix), intent(in) :: g
        real(dp), intent(in) :: scale
        integer :: k
        do k=g%rowptr(row),g%rowptr(row+1)-1
            call tb%add(outrow,g%colind(k),scale*g%values(k))
        end do
    end subroutine add_scaled_row

    subroutine add_row_outer(tb,g,row,coef)
        type(sparse_triplet_builder), intent(inout) :: tb
        type(ecos_csr_matrix), intent(in) :: g
        integer, intent(in) :: row
        real(dp), intent(in) :: coef
        integer :: ka,kb,ia,ib
        real(dp) :: va,vb
        do ka=g%rowptr(row),g%rowptr(row+1)-1
            ia=g%colind(ka); va=g%values(ka)
            do kb=ka,g%rowptr(row+1)-1
                ib=g%colind(kb); vb=g%values(kb)
                call tb%add(min(ia,ib),max(ia,ib),coef*va*vb)
            end do
        end do
    end subroutine add_row_outer

    subroutine add_dense_outer_touched(tb,idx,val,nt,coef)
        type(sparse_triplet_builder), intent(inout) :: tb
        integer, intent(in) :: idx(:),nt
        real(dp), intent(in) :: val(:),coef
        integer :: a,b,ia,ib
        do a=1,nt
            ia=idx(a)
            do b=a,nt
                ib=idx(b)
                call tb%add(min(ia,ib),max(ia,ib),coef*val(ia)*val(ib))
            end do
        end do
    end subroutine add_dense_outer_touched

    subroutine add_cross_rows(tb,g,rowa,rowb,coef,worka,workb,touched,mark,stamp)
        type(sparse_triplet_builder), intent(inout) :: tb
        type(ecos_csr_matrix), intent(in) :: g
        integer, intent(in) :: rowa,rowb
        real(dp), intent(in) :: coef
        real(dp), intent(inout) :: worka(:),workb(:)
        integer, intent(inout) :: touched(:),mark(:),stamp
        integer :: k,nt,i,j,a,b
        stamp=stamp+1; nt=0
        do k=g%rowptr(rowa),g%rowptr(rowa+1)-1
            i=g%colind(k)
            if(mark(i)/=stamp) then
                nt=nt+1; touched(nt)=i; mark(i)=stamp; worka(i)=0.0_dp; workb(i)=0.0_dp
            end if
            worka(i)=worka(i)+g%values(k)
        end do
        do k=g%rowptr(rowb),g%rowptr(rowb+1)-1
            i=g%colind(k)
            if(mark(i)/=stamp) then
                nt=nt+1; touched(nt)=i; mark(i)=stamp; worka(i)=0.0_dp; workb(i)=0.0_dp
            end if
            workb(i)=workb(i)+g%values(k)
        end do
        do a=1,nt
            i=touched(a)
            do b=a,nt
                j=touched(b)
                call tb%add(min(i,j),max(i,j),coef*(worka(i)*workb(j)+workb(i)*worka(j)))
            end do
        end do
    end subroutine add_cross_rows

    subroutine sparse_cone_values(prob,x,gval)
        type(ecos_problem), intent(in) :: prob
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: gval(:)
        real(dp), allocatable :: s(:)
        real(dp) :: t,nr,bb,cc,bs,cs
        integer :: row,idx,ir,iq,qd
        allocate(s(prob%ncone()))
        call sparse_cone_slack(prob,x,s)
        row=0; idx=0
        do ir=1,prob%dims%l
            row=row+1; idx=idx+1; gval(idx)=-s(row)
        end do
        if(allocated(prob%dims%q)) then
            do iq=1,size(prob%dims%q)
                qd=prob%dims%q(iq); idx=idx+1
                if(qd==1) then
                    gval(idx)=-s(row+1)
                else
                    t=s(row+1); nr=sqrt(dot_product(s(row+2:row+qd),s(row+2:row+qd)))
                    gval(idx)=nr-t
                end if
                row=row+qd
            end do
        end if
        do ir=1,prob%dims%e
            bb=s(row+2); cc=s(row+3); bs=max(bb,1.0e-12_dp); cs=max(cc,1.0e-12_dp)
            idx=idx+1; gval(idx)=-bb
            idx=idx+1; gval(idx)=-cc
            idx=idx+1; gval(idx)=s(row+1)-cc*log(bs/cs)
            row=row+3
        end do
    end subroutine sparse_cone_values

    subroutine sparse_cone_linearize(prob,x,lambda,gval,jac,hlag)
        type(ecos_problem), intent(in) :: prob
        real(dp), intent(in) :: x(:),lambda(:)
        real(dp), intent(out) :: gval(:)
        type(ecos_csr_matrix), intent(out) :: jac
        type(ecos_csc_matrix), intent(out) :: hlag
        type(sparse_triplet_builder) :: jb,hb
        real(dp), allocatable :: s(:),work(:),workb(:)
        integer, allocatable :: touched(:),mark(:)
        real(dp) :: t,nr,bb,cc,bs,cs,lg,coef
        integer :: row,idx,ir,iq,qd,k,kk,col,nt,stamp,n
        n=prob%nvar()
        call jb%init(size(gval),n,max(64,size(prob%g_csr%values)))
        call hb%init(n,n,max(64,size(prob%g_csr%values)))
        allocate(s(prob%ncone()),work(n),workb(n),touched(n),mark(n))
        work=0.0_dp; workb=0.0_dp; mark=0; stamp=0
        call sparse_cone_slack(prob,x,s)
        row=0; idx=0
        do ir=1,prob%dims%l
            row=row+1; idx=idx+1; gval(idx)=-s(row)
            call add_scaled_row(jb,idx,prob%g_csr,row,1.0_dp)
        end do
        if(allocated(prob%dims%q)) then
            do iq=1,size(prob%dims%q)
                qd=prob%dims%q(iq); idx=idx+1
                if(qd==1) then
                    gval(idx)=-s(row+1)
                    call add_scaled_row(jb,idx,prob%g_csr,row+1,1.0_dp)
                else
                    t=s(row+1)
                    nr=max(sqrt(dot_product(s(row+2:row+qd),s(row+2:row+qd))),1.0e-12_dp)
                    gval(idx)=nr-t
                    call add_scaled_row(jb,idx,prob%g_csr,row+1,1.0_dp)
                    do k=1,qd-1
                        call add_scaled_row(jb,idx,prob%g_csr,row+1+k,-s(row+1+k)/nr)
                    end do
                    if(lambda(idx)>0.0_dp) then
                        do k=1,qd-1
                            call add_row_outer(hb,prob%g_csr,row+1+k,lambda(idx)/nr)
                        end do
                        stamp=stamp+1; nt=0
                        do k=1,qd-1
                            do kk=prob%g_csr%rowptr(row+1+k),prob%g_csr%rowptr(row+2+k)-1
                                col=prob%g_csr%colind(kk)
                                if(mark(col)/=stamp) then
                                    nt=nt+1; touched(nt)=col; mark(col)=stamp; work(col)=0.0_dp
                                end if
                                work(col)=work(col)+s(row+1+k)*prob%g_csr%values(kk)
                            end do
                        end do
                        call add_dense_outer_touched(hb,touched,work,nt,-lambda(idx)/(nr**3))
                    end if
                end if
                row=row+qd
            end do
        end if
        do ir=1,prob%dims%e
            bb=s(row+2); cc=s(row+3); bs=max(bb,1.0e-12_dp); cs=max(cc,1.0e-12_dp)
            idx=idx+1; gval(idx)=-bb
            call add_scaled_row(jb,idx,prob%g_csr,row+2,1.0_dp)
            idx=idx+1; gval(idx)=-cc
            call add_scaled_row(jb,idx,prob%g_csr,row+3,1.0_dp)
            idx=idx+1; lg=log(bs/cs); gval(idx)=s(row+1)-cc*lg
            call add_scaled_row(jb,idx,prob%g_csr,row+1,-1.0_dp)
            call add_scaled_row(jb,idx,prob%g_csr,row+2,cc/bs)
            call add_scaled_row(jb,idx,prob%g_csr,row+3,-1.0_dp+lg)
            if(lambda(idx)>0.0_dp) then
                coef=lambda(idx)*cc/(bs*bs)
                call add_row_outer(hb,prob%g_csr,row+2,coef)
                coef=lambda(idx)/cs
                call add_row_outer(hb,prob%g_csr,row+3,coef)
                call add_cross_rows(hb,prob%g_csr,row+2,row+3,-lambda(idx)/bs, &
                                    work,workb,touched,mark,stamp)
            end if
            row=row+3
        end do
        call triplet_to_csr(jb,jac)
        call triplet_to_csc(hb,hlag,.true.)
    end subroutine sparse_cone_linearize

    subroutine sparse_limit_exp_step(prob,x,dx,alpha)
        type(ecos_problem), intent(in) :: prob
        real(dp), intent(in) :: x(:),dx(:)
        real(dp), intent(inout) :: alpha
        real(dp), allocatable :: s(:),gdx(:)
        real(dp) :: aa,ds
        integer :: row,k
        if(prob%dims%e<=0) return
        allocate(s(prob%ncone()),gdx(prob%ncone()))
        call sparse_cone_slack(prob,x,s)
        call csr_matvec(prob%g_csr,dx,gdx)
        row=prob%dims%l
        if(allocated(prob%dims%q)) row=row+sum(prob%dims%q)
        aa=alpha
        do k=1,prob%dims%e
            ds=-gdx(row+2)
            if(ds<0.0_dp) aa=min(aa,0.99_dp*max(s(row+2),1.0e-14_dp)/(-ds))
            ds=-gdx(row+3)
            if(ds<0.0_dp) aa=min(aa,0.99_dp*max(s(row+3),1.0e-14_dp)/(-ds))
            row=row+3
        end do
        alpha=max(1.0e-8_dp,min(alpha,aa))
    end subroutine sparse_limit_exp_step

end module ecos_sparse_cones
