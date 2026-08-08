! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of the ManifoldOptim/ROPTLIB computational interface.
! See LICENSE and UPSTREAM_PROVENANCE.md for upstream copyright and provenance.
module manifoldoptim_manifolds
  use manifoldoptim_kinds, only : dp
  use manifoldoptim_types
  use manifoldoptim_linalg
  implicit none
  private
  public :: project_tangent, retract_point, transport_vector, inverse_transport_vector
  public :: cotangent_vector
  public :: manifold_metric, manifold_beta, differentiated_retraction, metric_dual
  public :: euclidean_to_riemannian_gradient, euclidean_hess_to_riemannian
  public :: random_manifold_point
  public :: orthonorm, point_is_valid

contains

  subroutine project_tangent(domain, x, v, p)
    type(manifold_domain), intent(in) :: domain
    real(dp), intent(in) :: x(:), v(:)
    real(dp), intent(out) :: p(:)
    integer :: ic, j, pos, len

    if (size(x) /= domain%length() .or. size(v) /= size(x) .or. size(p) /= size(x)) &
      error stop 'project_tangent: bad size'
    pos = 1
    do ic = 1, size(domain%component)
      len = domain%component(ic)%block_length()
      do j = 1, domain%component(ic)%numofmani
        call project_block(domain%component(ic), x(pos:pos+len-1), &
          v(pos:pos+len-1), p(pos:pos+len-1))
        pos = pos + len
      end do
    end do
  end subroutine project_tangent

  subroutine project_block(c, x, v, p)
    type(manifold_component), intent(in) :: c
    real(dp), intent(in) :: x(:), v(:)
    real(dp), intent(out) :: p(:)
    integer :: n, m, r, l1, l2
    real(dp), allocatable :: xm(:,:), vm(:,:), pm(:,:), tmp(:,:), sm(:,:)
    real(dp), allocatable :: um(:,:), dm(:,:), wm(:,:), du(:,:), dd(:,:), dv(:,:)
    real(dp), allocatable :: ou(:,:), ov(:,:), duperp(:,:), dvperp(:,:), dc(:,:)

    select case (c%kind)
    case (MANI_EUCLIDEAN)
      p = v
    case (MANI_SPHERE)
      p = v - dot_product(x,v) * x
    case (MANI_STIEFEL, MANI_ORTHGROUP)
      n = c%n
      if (c%kind == MANI_ORTHGROUP) then
        r = n
      else
        r = c%p
      end if
      allocate(xm(n,r), vm(n,r), pm(n,r), tmp(r,r), sm(r,r))
      xm = reshape(x,[n,r])
      vm = reshape(v,[n,r])
      tmp = matmul(transpose(xm),vm)
      sm = 0.5_dp*(tmp+transpose(tmp))
      pm = vm - matmul(xm,sm)
      p = reshape(pm,[n*r])
    case (MANI_GRASSMANN)
      n = c%n
      r = c%p
      allocate(xm(n,r), vm(n,r), pm(n,r), tmp(r,r))
      xm = reshape(x,[n,r])
      vm = reshape(v,[n,r])
      tmp = matmul(transpose(xm),vm)
      pm = vm - matmul(xm,tmp)
      p = reshape(pm,[n*r])
    case (MANI_SPD)
      n = c%n
      allocate(vm(n,n),pm(n,n))
      vm = reshape(v,[n,n])
      pm = 0.5_dp*(vm+transpose(vm))
      p = reshape(pm,[n*n])
    case (MANI_LOWRANK)
      n = c%n
      m = c%m
      r = c%p
      l1 = n*r
      l2 = l1+r*r
      allocate(um(n,r), dm(r,r), wm(m,r), du(n,r), dd(r,r), dv(m,r))
      allocate(ou(r,r), ov(r,r), duperp(n,r), dvperp(m,r), dc(r,r))
      um = reshape(x(1:l1),[n,r])
      dm = reshape(x(l1+1:l2),[r,r])
      wm = reshape(x(l2+1:),[m,r])
      du = reshape(v(1:l1),[n,r])
      dd = reshape(v(l1+1:l2),[r,r])
      dv = reshape(v(l2+1:),[m,r])
      call skew_part(matmul(transpose(um),du),ou)
      call skew_part(matmul(transpose(wm),dv),ov)
      duperp = du - matmul(um,matmul(transpose(um),du))
      dvperp = dv - matmul(wm,matmul(transpose(wm),dv))
      dc = dd + matmul(ou,dm) + matmul(dm,transpose(ov))
      p(1:l1) = reshape(duperp,[l1])
      p(l1+1:l2) = reshape(dc,[r*r])
      p(l2+1:) = reshape(dvperp,[m*r])
    case default
      p = v
    end select
  end subroutine project_block

  subroutine euclidean_to_riemannian_gradient(domain, x, eg, rg)
    type(manifold_domain), intent(in) :: domain
    real(dp), intent(in) :: x(:), eg(:)
    real(dp), intent(out) :: rg(:)
    integer :: ic, j, pos, len, n
    real(dp), allocatable :: xm(:,:), gm(:,:), sm(:,:), tmp(:,:)

    pos = 1
    do ic = 1, size(domain%component)
      len = domain%component(ic)%block_length()
      do j = 1, domain%component(ic)%numofmani
        if (domain%component(ic)%kind == MANI_SPD) then
          n = domain%component(ic)%n
          allocate(xm(n,n),gm(n,n),sm(n,n),tmp(n,n))
          xm = reshape(x(pos:pos+len-1),[n,n])
          gm = reshape(eg(pos:pos+len-1),[n,n])
          sm = 0.5_dp*(gm+transpose(gm))
          tmp = matmul(xm,matmul(sm,xm))
          rg(pos:pos+len-1) = reshape(0.5_dp*(tmp+transpose(tmp)),[len])
          deallocate(xm,gm,sm,tmp)
        else
          call project_block(domain%component(ic), x(pos:pos+len-1), &
            eg(pos:pos+len-1), rg(pos:pos+len-1))
        end if
        pos = pos + len
      end do
    end do
  end subroutine euclidean_to_riemannian_gradient

  subroutine euclidean_hess_to_riemannian(domain, x, eta, egrad, eh, hv)
    type(manifold_domain), intent(in) :: domain
    real(dp), intent(in) :: x(:), eta(:), egrad(:), eh(:)
    real(dp), intent(out) :: hv(:)
    integer :: ic, j, pos, len, n, m, r, l1, l2
    real(dp), allocatable :: xm(:,:), em(:,:), gm(:,:), hm(:,:), corr(:,:)
    real(dp), allocatable :: work(:)

    pos = 1
    do ic = 1, size(domain%component)
      len = domain%component(ic)%block_length()
      do j = 1, domain%component(ic)%numofmani
        select case(domain%component(ic)%kind)
        case(MANI_SPHERE)
          hv(pos:pos+len-1) = eh(pos:pos+len-1) - &
            dot_product(x(pos:pos+len-1),egrad(pos:pos+len-1))*eta(pos:pos+len-1)
          hv(pos:pos+len-1) = hv(pos:pos+len-1) - &
            dot_product(x(pos:pos+len-1),hv(pos:pos+len-1))*x(pos:pos+len-1)
        case(MANI_STIEFEL,MANI_ORTHGROUP)
          n = domain%component(ic)%n
          if (domain%component(ic)%kind == MANI_ORTHGROUP) then
            r = n
          else
            r = domain%component(ic)%p
          end if
          allocate(xm(n,r),em(n,r),gm(n,r),hm(n,r),corr(r,r),work(len))
          xm = reshape(x(pos:pos+len-1),[n,r])
          em = reshape(eta(pos:pos+len-1),[n,r])
          gm = reshape(egrad(pos:pos+len-1),[n,r])
          hm = reshape(eh(pos:pos+len-1),[n,r])
          corr = matmul(transpose(xm),gm)
          corr = 0.5_dp*(corr+transpose(corr))
          work = reshape(hm-matmul(em,corr),[len])
          call project_block(domain%component(ic),x(pos:pos+len-1),work,hv(pos:pos+len-1))
          deallocate(xm,em,gm,hm,corr,work)
        case(MANI_GRASSMANN)
          n = domain%component(ic)%n
          r = domain%component(ic)%p
          allocate(xm(n,r),em(n,r),gm(n,r),hm(n,r),corr(r,r),work(len))
          xm = reshape(x(pos:pos+len-1),[n,r])
          em = reshape(eta(pos:pos+len-1),[n,r])
          gm = reshape(egrad(pos:pos+len-1),[n,r])
          hm = reshape(eh(pos:pos+len-1),[n,r])
          corr = matmul(transpose(xm),gm)
          work = reshape(hm-matmul(em,corr),[len])
          call project_block(domain%component(ic),x(pos:pos+len-1),work,hv(pos:pos+len-1))
          deallocate(xm,em,gm,hm,corr,work)
        case(MANI_LOWRANK)
          n = domain%component(ic)%n
          m = domain%component(ic)%m
          r = domain%component(ic)%p
          l1 = n*r
          l2 = l1+r*r
          allocate(work(len))
          work = eh(pos:pos+len-1)
          block
            real(dp) :: xu(n,r), eu(n,r), gu(n,r), hu(n,r), cu(r,r)
            real(dp) :: xv(m,r), ev(m,r), gv(m,r), hvv(m,r), cv(r,r)
            xu = reshape(x(pos:pos+l1-1),[n,r])
            eu = reshape(eta(pos:pos+l1-1),[n,r])
            gu = reshape(egrad(pos:pos+l1-1),[n,r])
            hu = reshape(eh(pos:pos+l1-1),[n,r])
            cu = matmul(transpose(xu),gu)
            cu = 0.5_dp*(cu+transpose(cu))
            work(1:l1) = reshape(hu-matmul(eu,cu),[l1])
            xv = reshape(x(pos+l2:pos+len-1),[m,r])
            ev = reshape(eta(pos+l2:pos+len-1),[m,r])
            gv = reshape(egrad(pos+l2:pos+len-1),[m,r])
            hvv = reshape(eh(pos+l2:pos+len-1),[m,r])
            cv = matmul(transpose(xv),gv)
            cv = 0.5_dp*(cv+transpose(cv))
            work(l2+1:len) = reshape(hvv-matmul(ev,cv),[m*r])
          end block
          call project_block(domain%component(ic),x(pos:pos+len-1),work,hv(pos:pos+len-1))
          deallocate(work)
        case default
          call project_block(domain%component(ic),x(pos:pos+len-1), &
            eh(pos:pos+len-1),hv(pos:pos+len-1))
        end select
        pos = pos+len
      end do
    end do
  end subroutine euclidean_hess_to_riemannian

  subroutine retract_point(domain, x, eta, y, ok)
    type(manifold_domain), intent(in) :: domain
    real(dp), intent(in) :: x(:), eta(:)
    real(dp), intent(out) :: y(:)
    logical, intent(out), optional :: ok
    integer :: ic, j, pos, len
    logical :: lok, allok

    pos = 1
    allok = .true.
    do ic = 1, size(domain%component)
      len = domain%component(ic)%block_length()
      do j = 1, domain%component(ic)%numofmani
        call retract_block(domain%component(ic), x(pos:pos+len-1), &
          eta(pos:pos+len-1), y(pos:pos+len-1), lok)
        allok = allok .and. lok
        pos = pos + len
      end do
    end do
    if (present(ok)) ok = allok
  end subroutine retract_point

  subroutine retract_block(c, x, eta, y, ok)
    type(manifold_component), intent(in) :: c
    real(dp), intent(in) :: x(:), eta(:)
    real(dp), intent(out) :: y(:)
    logical, intent(out) :: ok
    integer :: n, m, r, l1, l2, info, i
    real(dp) :: nr
    real(dp), allocatable :: a(:,:), q(:,:), xm(:,:), em(:,:), isx(:,:), ym(:,:)
    real(dp), allocatable :: d(:), ev(:,:)

    select case(c%kind)
    case(MANI_EUCLIDEAN)
      y = x + eta
      ok = .true.
    case(MANI_SPHERE)
      if (c%param_set == 2) then
        nr = vecnorm(eta)
        if (nr <= 100.0_dp*epsilon(1.0_dp)) then
          y = x
        else
          y = cos(nr)*x + (sin(nr)/nr)*eta
          y = y/max(vecnorm(y),tiny(1.0_dp))
        end if
        ok = .true.
      else
        y = x + eta
        nr = vecnorm(y)
        ok = nr > 100.0_dp*epsilon(1.0_dp)
        if (ok) y = y/nr
      end if
    case(MANI_STIEFEL)
      n = c%n
      r = c%p
      if (c%param_set == 2) then
        call stiefel_constructed_retraction(n,r,x,eta,y,ok)
      else
        allocate(a(n,r),q(n,r))
        a = reshape(x+eta,[n,r])
        call mgs_orthonormalize(a,q,ok)
        y = reshape(q,[n*r])
      end if
    case(MANI_GRASSMANN,MANI_ORTHGROUP)
      n = c%n
      if (c%kind == MANI_ORTHGROUP) then
        r = n
      else
        r = c%p
      end if
      allocate(a(n,r),q(n,r))
      a = reshape(x+eta,[n,r])
      call mgs_orthonormalize(a,q,ok)
      y = reshape(q,[n*r])
    case(MANI_SPD)
      n = c%n
      allocate(xm(n,n),em(n,n),isx(n,n),ym(n,n),d(n),ev(n,n))
      xm = reshape(x,[n,n])
      em = reshape(eta,[n,n])
      em = 0.5_dp*(em+transpose(em))
      call matrix_inverse(xm,isx,info)
      if (info /= 0) then
        ok = .false.
        y = x
        return
      end if
      ym = xm + em + 0.5_dp*matmul(em,matmul(isx,em))
      ym = 0.5_dp*(ym+transpose(ym))
      call symmetric_eigen(ym,d,ev,info)
      ok = info == 0 .and. minval(d) > 0.0_dp
      if (.not. ok .and. info == 0) then
        do i = 1, n
          d(i) = max(d(i),1.0e-12_dp)
        end do
        ym = matmul(ev,matmul(diag_matrix(d),transpose(ev)))
        ok = .true.
      end if
      y = reshape(ym,[n*n])
    case(MANI_LOWRANK)
      n = c%n
      m = c%m
      r = c%p
      l1 = n*r
      l2 = l1+r*r
      allocate(a(max(n,m),r),q(max(n,m),r))
      a = 0.0_dp
      q = 0.0_dp
      a(1:n,:) = reshape(x(1:l1)+eta(1:l1),[n,r])
      call mgs_orthonormalize(a(1:n,:),q(1:n,:),ok)
      if (.not. ok) return
      y(1:l1) = reshape(q(1:n,:),[l1])
      y(l1+1:l2) = x(l1+1:l2)+eta(l1+1:l2)
      a = 0.0_dp
      q = 0.0_dp
      a(1:m,:) = reshape(x(l2+1:)+eta(l2+1:),[m,r])
      call mgs_orthonormalize(a(1:m,:),q(1:m,:),ok)
      if (.not. ok) return
      y(l2+1:) = reshape(q(1:m,:),[m*r])
    case default
      y = x + eta
      ok = .true.
    end select
  end subroutine retract_block

  subroutine stiefel_constructed_retraction(n, p, x, eta, y, ok)
    integer, intent(in) :: n, p
    real(dp), intent(in) :: x(:), eta(:)
    real(dp), intent(out) :: y(:)
    logical, intent(out) :: ok
    integer :: nc, info
    real(dp), allocatable :: xm(:,:), em(:,:), xp(:,:), frame(:,:)
    real(dp), allocatable :: a(:,:), k(:,:), omega(:,:), eomega(:,:), ym(:,:)

    nc = n-p
    allocate(xm(n,p),em(n,p),xp(n,nc),frame(n,n))
    allocate(a(p,p),k(nc,p),omega(n,n),eomega(n,n),ym(n,p))
    xm = reshape(x,[n,p])
    em = reshape(eta,[n,p])
    call orthogonal_complement(xm,xp,ok)
    if (.not. ok) then
      y = x
      return
    end if
    a = matmul(transpose(xm),em)
    a = 0.5_dp*(a-transpose(a))
    if (nc > 0) k = matmul(transpose(xp),em)
    omega = 0.0_dp
    omega(1:p,1:p) = a
    if (nc > 0) then
      omega(1:p,p+1:n) = -transpose(k)
      omega(p+1:n,1:p) = k
    end if
    call matrix_exponential(omega,eomega,info)
    if (info /= 0) then
      ok = .false.
      y = x
      return
    end if
    frame(:,1:p) = xm
    if (nc > 0) frame(:,p+1:n) = xp
    ym = matmul(frame,eomega(:,1:p))
    y = reshape(ym,[n*p])
    ok = .true.
  end subroutine stiefel_constructed_retraction

  pure function diag_matrix(d) result(a)
    real(dp), intent(in) :: d(:)
    real(dp) :: a(size(d),size(d))
    integer :: i
    a = 0.0_dp
    do i = 1, size(d)
      a(i,i) = 1.0_dp*d(i)
    end do
  end function diag_matrix

  subroutine transport_vector(domain, x, y, v, tv)
    type(manifold_domain), intent(in) :: domain
    real(dp), intent(in) :: x(:), y(:), v(:)
    real(dp), intent(out) :: tv(:)
    integer :: ic, j, pos, len

    if (size(x) /= size(y) .or. size(v) /= size(x) .or. size(tv) /= size(x)) &
      error stop 'transport_vector: incompatible sizes'
    pos = 1
    do ic = 1, size(domain%component)
      len = domain%component(ic)%block_length()
      do j = 1, domain%component(ic)%numofmani
        call transport_block(domain%component(ic),x(pos:pos+len-1), &
          y(pos:pos+len-1),v(pos:pos+len-1),tv(pos:pos+len-1))
        pos = pos + len
      end do
    end do
  end subroutine transport_vector

  subroutine inverse_transport_vector(domain, x, y, v, tv)
    type(manifold_domain), intent(in) :: domain
    real(dp), intent(in) :: x(:), y(:), v(:)
    real(dp), intent(out) :: tv(:)
    call transport_vector(domain,y,x,v,tv)
  end subroutine inverse_transport_vector


  subroutine cotangent_vector(domain, x, eta, y, z, cz)
    type(manifold_domain), intent(in) :: domain
    real(dp), intent(in) :: x(:), eta(:), y(:), z(:)
    real(dp), intent(out) :: cz(:)
    integer :: ic, j, pos, len

    pos = 1
    do ic = 1, size(domain%component)
      len = domain%component(ic)%block_length()
      do j = 1, domain%component(ic)%numofmani
        call cotangent_block(domain%component(ic),x(pos:pos+len-1), &
          eta(pos:pos+len-1),y(pos:pos+len-1),z(pos:pos+len-1), &
          cz(pos:pos+len-1))
        pos = pos+len
      end do
    end do
  end subroutine cotangent_vector

  subroutine cotangent_block(c, x, eta, y, z, cz)
    type(manifold_component), intent(in) :: c
    real(dp), intent(in) :: x(:), eta(:), y(:), z(:)
    real(dp), intent(out) :: cz(:)
    real(dp) :: nr, ztx, zte, sn, cs

    select case(c%kind)
    case(MANI_EUCLIDEAN,MANI_SPD)
      cz = z
    case(MANI_SPHERE)
      if (c%param_set == 2) then
        nr = vecnorm(eta)
        if (nr <= 100.0_dp*epsilon(1.0_dp)) then
          cz = z-dot_product(x,z)*x
        else
          ztx = dot_product(x,z)
          zte = dot_product(z,eta)
          sn = sin(nr)
          cs = cos(nr)
          cz = (sn/nr)*z
          cz = cz + (zte*cs/nr-ztx*sn-zte*sn/(nr*nr))/nr*eta
          cz = cz-dot_product(x,cz)*x
        end if
      else
        call qf_cotangent_matrix(c%n,1,.true.,x,eta,y,z,cz)
      end if
    case(MANI_STIEFEL)
      if (c%param_set == 2) then
        call frame_transport(c%n,c%p,y,x,z,cz)
      else
        call qf_cotangent_matrix(c%n,c%p,.true.,x,eta,y,z,cz)
      end if
    case(MANI_GRASSMANN)
      call qf_cotangent_matrix(c%n,c%p,.false.,x,eta,y,z,cz)
    case(MANI_ORTHGROUP)
      call qf_cotangent_matrix(c%n,c%n,.true.,x,eta,y,z,cz)
    case(MANI_LOWRANK)
      ! LowRank has a specialized upstream cotangent formula. The intrinsic
      ! coordinate transport is its stable fallback for the public flat API.
      call transport_block(c,y,x,z,cz)
    case default
      cz = z
    end select
  end subroutine cotangent_block

  subroutine qf_cotangent_matrix(n, p, stiefel_case, x, eta, y, z, cz)
    integer, intent(in) :: n, p
    logical, intent(in) :: stiefel_case
    real(dp), intent(in) :: x(:), eta(:), y(:), z(:)
    real(dp), intent(out) :: cz(:)
    integer :: i, j, info
    real(dp) :: xm(n,p), ym(n,p), zm(n,p), am(n,p), rr(p,p), ir(p,p)
    real(dp) :: tt(p,p), wm(n,p), pm(n,p)

    xm = reshape(x,[n,p])
    ym = reshape(y,[n,p])
    zm = reshape(z,[n,p])
    am = reshape(x+eta,[n,p])
    rr = matmul(transpose(ym),am)
    call matrix_inverse(rr,ir,info)
    if (info /= 0) then
      if (stiefel_case) then
        call stiefel_project_matrix(xm,zm,pm)
      else
        pm = zm-matmul(xm,matmul(transpose(xm),zm))
      end if
      cz = reshape(pm,[n*p])
      return
    end if
    tt = matmul(transpose(ym),zm)
    if (stiefel_case) then
      do j = 1, p
        do i = 1, j
          tt(i,j) = -tt(i,j)
        end do
      end do
    end if
    wm = zm+matmul(ym,tt)
    wm = matmul(wm,transpose(ir))
    if (stiefel_case) then
      call stiefel_project_matrix(xm,wm,pm)
    else
      pm = wm-matmul(xm,matmul(transpose(xm),wm))
    end if
    cz = reshape(pm,[n*p])
  end subroutine qf_cotangent_matrix

  subroutine transport_block(c, x, y, v, tv)
    type(manifold_component), intent(in) :: c
    real(dp), intent(in) :: x(:), y(:), v(:)
    real(dp), intent(out) :: tv(:)
    integer :: n, m, r, nc, l1, l2, info
    real(dp) :: den, coeff
    logical :: okx, oky
    real(dp), allocatable :: xm(:,:), ym(:,:), vm(:,:), xp(:,:), yp(:,:), k(:,:), tm(:,:)
    real(dp), allocatable :: px(:), um(:,:), dm(:,:), wm(:,:), dun(:,:), dd(:,:), dvn(:,:)
    real(dp), allocatable :: uy(:,:), dy(:,:), wy(:,:), up(:,:), uyp(:,:), wp(:,:), wyp(:,:)
    real(dp), allocatable :: ia(:,:), ib(:,:), dinv(:,:), dtinv(:,:)

    select case(c%kind)
    case(MANI_EUCLIDEAN)
      tv = v
    case(MANI_SPHERE)
      if (c%param_set >= 2 .and. c%param_set <= 4) then
        den = dot_product(x+y,x+y)
        if (den <= tiny(1.0_dp)) then
          tv = v-dot_product(y,v)*y
        else
          coeff = 2.0_dp*dot_product(v,y)/den
          tv = v-coeff*(x+y)
          tv = tv-dot_product(y,tv)*y
        end if
      else
        call frame_transport(c%n,1,x,y,v,tv)
      end if
    case(MANI_STIEFEL)
      call frame_transport(c%n,c%p,x,y,v,tv)
    case(MANI_GRASSMANN)
      n = c%n
      r = c%p
      nc = n-r
      allocate(xm(n,r),ym(n,r),vm(n,r),xp(n,nc),yp(n,nc),k(nc,r),tm(n,r))
      xm = reshape(x,[n,r])
      ym = reshape(y,[n,r])
      vm = reshape(v,[n,r])
      call orthogonal_complement(xm,xp,okx)
      call orthogonal_complement(ym,yp,oky)
      if (okx .and. oky) then
        if (nc > 0) then
          k = matmul(transpose(xp),vm)
          tm = matmul(yp,k)
        else
          tm = 0.0_dp
        end if
        tv = reshape(tm,[n*r])
      else
        call project_block(c,y,v,tv)
      end if
    case(MANI_SPD)
      call project_block(c,y,v,tv)
    case(MANI_ORTHGROUP)
      call frame_transport(c%n,c%n,x,y,v,tv)
    case(MANI_LOWRANK)
      n = c%n
      m = c%m
      r = c%p
      l1 = n*r
      l2 = l1+r*r
      allocate(px(size(v)))
      call project_block(c,x,v,px)
      allocate(um(n,r),dm(r,r),wm(m,r),dun(n,r),dd(r,r),dvn(m,r))
      allocate(uy(n,r),dy(r,r),wy(m,r),up(n,n-r),uyp(n,n-r))
      allocate(wp(m,m-r),wyp(m,m-r),ia(n-r,r),ib(m-r,r),dinv(r,r),dtinv(r,r))
      um = reshape(x(1:l1),[n,r])
      dm = reshape(x(l1+1:l2),[r,r])
      wm = reshape(x(l2+1:),[m,r])
      dun = reshape(px(1:l1),[n,r])
      dd = reshape(px(l1+1:l2),[r,r])
      dvn = reshape(px(l2+1:),[m,r])
      uy = reshape(y(1:l1),[n,r])
      dy = reshape(y(l1+1:l2),[r,r])
      wy = reshape(y(l2+1:),[m,r])
      call orthogonal_complement(um,up,okx)
      call orthogonal_complement(uy,uyp,oky)
      if (.not. (okx .and. oky)) then
        call project_block(c,y,v,tv)
        return
      end if
      call orthogonal_complement(wm,wp,okx)
      call orthogonal_complement(wy,wyp,oky)
      if (.not. (okx .and. oky)) then
        call project_block(c,y,v,tv)
        return
      end if
      call matrix_inverse(dy,dinv,info)
      if (info /= 0) then
        call project_block(c,y,v,tv)
        return
      end if
      call matrix_inverse(transpose(dy),dtinv,info)
      if (info /= 0) then
        call project_block(c,y,v,tv)
        return
      end if
      if (n-r > 0) ia = matmul(transpose(up),dun)
      if (m-r > 0) ib = matmul(transpose(wp),dvn)
      if (n-r > 0) dun = matmul(uyp,matmul(matmul(ia,dm),dinv))
      if (n-r == 0) dun = 0.0_dp
      if (m-r > 0) dvn = matmul(wyp,matmul(matmul(ib,transpose(dm)),dtinv))
      if (m-r == 0) dvn = 0.0_dp
      tv(1:l1) = reshape(dun,[l1])
      tv(l1+1:l2) = reshape(dd,[r*r])
      tv(l2+1:) = reshape(dvn,[m*r])
    case default
      call project_block(c,y,v,tv)
    end select
  end subroutine transport_block

  subroutine frame_transport(n, p, x, y, v, tv)
    integer, intent(in) :: n, p
    real(dp), intent(in) :: x(:), y(:), v(:)
    real(dp), intent(out) :: tv(:)
    integer :: nc
    logical :: okx, oky
    real(dp), allocatable :: xm(:,:),ym(:,:),vm(:,:),xp(:,:),yp(:,:),a(:,:),k(:,:),tm(:,:)

    nc = n-p
    allocate(xm(n,p),ym(n,p),vm(n,p),xp(n,nc),yp(n,nc),a(p,p),k(nc,p),tm(n,p))
    xm = reshape(x,[n,p])
    ym = reshape(y,[n,p])
    vm = reshape(v,[n,p])
    call orthogonal_complement(xm,xp,okx)
    call orthogonal_complement(ym,yp,oky)
    if (.not. (okx .and. oky)) then
      call stiefel_project_matrix(ym,vm,tm)
      tv = reshape(tm,[n*p])
      return
    end if
    a = matmul(transpose(xm),vm)
    a = 0.5_dp*(a-transpose(a))
    if (nc > 0) k = matmul(transpose(xp),vm)
    tm = matmul(ym,a)
    if (nc > 0) tm = tm+matmul(yp,k)
    tv = reshape(tm,[n*p])
  end subroutine frame_transport

  subroutine stiefel_project_matrix(x, v, p)
    real(dp), intent(in) :: x(:,:), v(:,:)
    real(dp), intent(out) :: p(:,:)
    real(dp) :: tmp(size(x,2),size(x,2))
    tmp = matmul(transpose(x),v)
    p = v-matmul(x,0.5_dp*(tmp+transpose(tmp)))
  end subroutine stiefel_project_matrix

  subroutine differentiated_retraction(domain, x, eta, xi, y, dxi)
    type(manifold_domain), intent(in) :: domain
    real(dp), intent(in) :: x(:), eta(:), xi(:), y(:)
    real(dp), intent(out) :: dxi(:)
    integer :: ic, j, pos, len

    pos = 1
    do ic = 1, size(domain%component)
      len = domain%component(ic)%block_length()
      do j = 1, domain%component(ic)%numofmani
        call diff_retract_block(domain%component(ic),x(pos:pos+len-1), &
          eta(pos:pos+len-1),xi(pos:pos+len-1),y(pos:pos+len-1), &
          dxi(pos:pos+len-1))
        pos = pos+len
      end do
    end do
  end subroutine differentiated_retraction

  subroutine diff_retract_block(c, x, eta, xi, y, dxi)
    type(manifold_component), intent(in) :: c
    real(dp), intent(in) :: x(:), eta(:), xi(:), y(:)
    real(dp), intent(out) :: dxi(:)
    real(dp) :: nr, er, etaxi

    select case(c%kind)
    case(MANI_EUCLIDEAN)
      dxi = xi
    case(MANI_SPHERE)
      if (c%param_set == 2) then
        nr = vecnorm(eta)
        if (nr <= 100.0_dp*epsilon(1.0_dp)) then
          dxi = xi
        else
          etaxi = dot_product(eta,xi)
          dxi = -sin(nr)*etaxi/nr*x + sin(nr)/nr*xi
          dxi = dxi + (cos(nr)-sin(nr)/nr)*etaxi/(nr*nr)*eta
        end if
      else
        er = max(vecnorm(x+eta),tiny(1.0_dp))
        dxi = (xi-y*dot_product(y,xi))/er
      end if
    case default
      ! For the remaining domains the transported vector is the stable public
      ! representation of the differentiated retraction used by the solvers.
      call transport_block(c,x,y,xi,dxi)
    end select
  end subroutine diff_retract_block

  real(dp) function manifold_beta(domain, x, eta) result(beta)
    type(manifold_domain), intent(in) :: domain
    real(dp), intent(in) :: x(:), eta(:)
    real(dp), allocatable :: y(:), d(:)
    real(dp) :: a, b
    logical :: ok
    integer :: ic, j, pos, len
    logical :: need_beta

    need_beta = .false.
    do ic = 1, size(domain%component)
      do j = 1, domain%component(ic)%numofmani
        if (domain%component(ic)%kind == MANI_SPHERE .and. &
            domain%component(ic)%param_set == 4) need_beta = .true.
      end do
    end do
    if (.not. need_beta) then
      beta = 1.0_dp
      return
    end if
    allocate(y(size(x)),d(size(x)))
    call retract_point(domain,x,eta,y,ok)
    if (.not. ok) then
      beta = 1.0_dp
      return
    end if
    d = 0.0_dp
    pos = 1
    do ic = 1, size(domain%component)
      len = domain%component(ic)%block_length()
      do j = 1, domain%component(ic)%numofmani
        if (domain%component(ic)%kind == MANI_SPHERE .and. &
            domain%component(ic)%param_set == 4) then
          call diff_retract_block(domain%component(ic),x(pos:pos+len-1), &
            eta(pos:pos+len-1),eta(pos:pos+len-1),y(pos:pos+len-1), &
            d(pos:pos+len-1))
        else
          call transport_block(domain%component(ic),x(pos:pos+len-1), &
            y(pos:pos+len-1),eta(pos:pos+len-1),d(pos:pos+len-1))
        end if
        pos = pos+len
      end do
    end do
    a = manifold_metric(domain,x,eta,eta)
    b = manifold_metric(domain,y,d,d)
    if (a > tiny(1.0_dp) .and. b > tiny(1.0_dp)) then
      beta = sqrt(a/b)
    else
      beta = 1.0_dp
    end if
  end function manifold_beta

  subroutine metric_dual(domain, x, u, du)
    type(manifold_domain), intent(in) :: domain
    real(dp), intent(in) :: x(:), u(:)
    real(dp), intent(out) :: du(:)
    integer :: ic, j, pos, len, n, m, r, l1, l2, info
    real(dp), allocatable :: xm(:,:), im(:,:), um(:,:), dm(:,:), vm(:,:)

    pos = 1
    do ic = 1, size(domain%component)
      len = domain%component(ic)%block_length()
      do j = 1, domain%component(ic)%numofmani
        select case(domain%component(ic)%kind)
        case(MANI_SPD)
          n = domain%component(ic)%n
          allocate(xm(n,n),im(n,n),um(n,n))
          xm = reshape(x(pos:pos+len-1),[n,n])
          um = reshape(u(pos:pos+len-1),[n,n])
          call matrix_inverse(xm,im,info)
          if (info == 0) then
            um = matmul(im,matmul(um,im))
            du(pos:pos+len-1) = reshape(0.5_dp*(um+transpose(um)),[len])
          else
            du(pos:pos+len-1) = u(pos:pos+len-1)
          end if
          deallocate(xm,im,um)
        case(MANI_LOWRANK)
          n = domain%component(ic)%n
          m = domain%component(ic)%m
          r = domain%component(ic)%p
          l1 = n*r
          l2 = l1+r*r
          allocate(um(n,r),dm(r,r),vm(m,r))
          um = reshape(u(pos:pos+l1-1),[n,r])
          dm = reshape(x(pos+l1:pos+l2-1),[r,r])
          vm = reshape(u(pos+l2:pos+len-1),[m,r])
          um = matmul(um,matmul(dm,transpose(dm)))
          vm = matmul(vm,matmul(transpose(dm),dm))
          du(pos:pos+l1-1) = reshape(um,[l1])
          du(pos+l1:pos+l2-1) = u(pos+l1:pos+l2-1)
          du(pos+l2:pos+len-1) = reshape(vm,[m*r])
        case default
          du(pos:pos+len-1) = u(pos:pos+len-1)
        end select
        pos = pos+len
      end do
    end do
  end subroutine metric_dual

  real(dp) function manifold_metric(domain, x, u, v) result(val)
    type(manifold_domain), intent(in) :: domain
    real(dp), intent(in) :: x(:), u(:), v(:)
    integer :: ic, j, pos, len, n, m, r, l1, l2, info
    real(dp), allocatable :: xm(:,:), im(:,:), um(:,:), vm(:,:), a(:,:), b(:,:)
    real(dp), allocatable :: dm(:,:), uu(:,:), uv(:,:), vu(:,:), vv(:,:)

    val = 0.0_dp
    pos = 1
    do ic = 1, size(domain%component)
      len = domain%component(ic)%block_length()
      do j = 1, domain%component(ic)%numofmani
        select case(domain%component(ic)%kind)
        case(MANI_SPD)
          n = domain%component(ic)%n
          allocate(xm(n,n),im(n,n),um(n,n),vm(n,n),a(n,n),b(n,n))
          xm = reshape(x(pos:pos+len-1),[n,n])
          um = reshape(u(pos:pos+len-1),[n,n])
          vm = reshape(v(pos:pos+len-1),[n,n])
          call matrix_inverse(xm,im,info)
          if (info == 0) then
            a = matmul(im,um)
            b = matmul(im,vm)
            val = val+trace_matrix(matmul(a,b))
          else
            val = val+dot_product(u(pos:pos+len-1),v(pos:pos+len-1))
          end if
          deallocate(xm,im,um,vm,a,b)
        case(MANI_LOWRANK)
          n = domain%component(ic)%n
          m = domain%component(ic)%m
          r = domain%component(ic)%p
          l1 = n*r
          l2 = l1+r*r
          allocate(dm(r,r),uu(n,r),uv(n,r),vu(m,r),vv(m,r))
          dm = reshape(x(pos+l1:pos+l2-1),[r,r])
          uu = reshape(u(pos:pos+l1-1),[n,r])
          uv = reshape(v(pos:pos+l1-1),[n,r])
          vu = reshape(u(pos+l2:pos+len-1),[m,r])
          vv = reshape(v(pos+l2:pos+len-1),[m,r])
          val = val+sum((matmul(uu,dm))*(matmul(uv,dm)))
          val = val+dot_product(u(pos+l1:pos+l2-1),v(pos+l1:pos+l2-1))
          val = val+sum((matmul(vu,transpose(dm)))*(matmul(vv,transpose(dm))))
        case default
          val = val+dot_product(u(pos:pos+len-1),v(pos:pos+len-1))
        end select
        pos = pos+len
      end do
    end do
  end function manifold_metric

  subroutine random_manifold_point(domain, x)
    type(manifold_domain), intent(in) :: domain
    real(dp), intent(out) :: x(:)
    integer :: ic, j, pos, len, n, r, m, l1, l2, k
    logical :: ok

    if (size(x) /= domain%length()) error stop 'random_manifold_point: bad size'
    call random_number(x)
    x = 2.0_dp*x-1.0_dp
    pos = 1
    do ic = 1, size(domain%component)
      len = domain%component(ic)%block_length()
      do j = 1, domain%component(ic)%numofmani
        select case(domain%component(ic)%kind)
        case(MANI_SPHERE)
          x(pos:pos+len-1) = x(pos:pos+len-1)/max(vecnorm(x(pos:pos+len-1)),tiny(1.0_dp))
        case(MANI_STIEFEL,MANI_GRASSMANN,MANI_ORTHGROUP)
          n = domain%component(ic)%n
          if (domain%component(ic)%kind == MANI_ORTHGROUP) then
            r = n
          else
            r = domain%component(ic)%p
          end if
          call orthonorm(x(pos:pos+len-1),n,r,x(pos:pos+len-1),ok)
        case(MANI_SPD)
          n = domain%component(ic)%n
          block
            real(dp) :: aa(n,n), ss(n,n)
            aa = reshape(x(pos:pos+len-1),[n,n])
            ss = matmul(aa,transpose(aa))
            do k = 1, n
              ss(k,k) = ss(k,k)+1.0_dp
            end do
            x(pos:pos+len-1) = reshape(ss,[n*n])
          end block
        case(MANI_LOWRANK)
          n = domain%component(ic)%n
          m = domain%component(ic)%m
          r = domain%component(ic)%p
          l1 = n*r
          l2 = l1+r*r
          call orthonorm(x(pos:pos+l1-1),n,r,x(pos:pos+l1-1),ok)
          x(pos+l1:pos+l2-1) = 0.0_dp
          do k = 1, r
            x(pos+l1-1+(k-1)*(r+1)+1) = 1.0_dp
          end do
          call orthonorm(x(pos+l2:pos+len-1),m,r,x(pos+l2:pos+len-1),ok)
        end select
        pos = pos+len
      end do
    end do
  end subroutine random_manifold_point

  subroutine orthonorm(x, n, p, qflat, ok)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: n, p
    real(dp), intent(out) :: qflat(:)
    logical, intent(out), optional :: ok
    real(dp) :: a(n,p), q(n,p)
    logical :: lok

    a = reshape(x,[n,p])
    call mgs_orthonormalize(a,q,lok)
    qflat = reshape(q,[n*p])
    if (present(ok)) ok = lok
  end subroutine orthonorm

  logical function point_is_valid(domain, x, tol) result(ok)
    type(manifold_domain), intent(in) :: domain
    real(dp), intent(in) :: x(:)
    real(dp), intent(in), optional :: tol
    real(dp) :: t
    integer :: ic, j, pos, len, n, r, m, l1, l2
    real(dp), allocatable :: a(:,:), gram(:,:)

    t = 1.0e-8_dp
    if (present(tol)) t = tol
    ok = size(x) == domain%length()
    if (.not. ok) return
    pos = 1
    do ic = 1, size(domain%component)
      len = domain%component(ic)%block_length()
      do j = 1, domain%component(ic)%numofmani
        select case(domain%component(ic)%kind)
        case(MANI_SPHERE)
          if (abs(vecnorm(x(pos:pos+len-1))-1.0_dp) > t) ok = .false.
        case(MANI_STIEFEL,MANI_GRASSMANN,MANI_ORTHGROUP)
          n = domain%component(ic)%n
          if (domain%component(ic)%kind == MANI_ORTHGROUP) then
            r = n
          else
            r = domain%component(ic)%p
          end if
          allocate(a(n,r),gram(r,r))
          a = reshape(x(pos:pos+len-1),[n,r])
          gram = matmul(transpose(a),a)-eye_matrix(r)
          if (maxval(abs(gram)) > t) ok = .false.
          deallocate(a,gram)
        case(MANI_SPD)
          n = domain%component(ic)%n
          allocate(a(n,n))
          a = reshape(x(pos:pos+len-1),[n,n])
          if (.not. is_spd(a,t)) ok = .false.
          deallocate(a)
        case(MANI_LOWRANK)
          n = domain%component(ic)%n
          m = domain%component(ic)%m
          r = domain%component(ic)%p
          l1 = n*r
          l2 = l1+r*r
          allocate(a(n,r),gram(r,r))
          a = reshape(x(pos:pos+l1-1),[n,r])
          gram = matmul(transpose(a),a)-eye_matrix(r)
          if (maxval(abs(gram)) > t) ok = .false.
          deallocate(a,gram)
          allocate(a(m,r),gram(r,r))
          a = reshape(x(pos+l2:pos+len-1),[m,r])
          gram = matmul(transpose(a),a)-eye_matrix(r)
          if (maxval(abs(gram)) > t) ok = .false.
          deallocate(a,gram)
        end select
        if (.not. ok) return
        pos = pos+len
      end do
    end do
  end function point_is_valid

end module manifoldoptim_manifolds
