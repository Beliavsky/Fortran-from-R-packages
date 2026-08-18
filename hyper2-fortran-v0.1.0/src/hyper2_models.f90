! SPDX-License-Identifier: GPL-2.0-or-later
module hyper2_models
    use hyper2_kinds, only : dp, name_len
    use hyper2_types, only : hyper2_model, hyper3_model, hyper2_create, hyper3_create, player_index
    implicit none
    private

    public :: dirichlet, generalized_dirichlet, gd_wong, dirichlet3
    public :: hyper2_matrix, hyper3_matrix
    public :: zipf, trial, pick, pass_model, beats, rankvec_likelihood, ordervec2supp
    public :: ordervec2supp3, ordervec2supp3a, pairwise, home_away, home_away3, home_draw_away3
    public :: zermelo, hyper3_to_hyper2, as_hyper3, setweight, pwa, pwa3, pwa23
    public :: keep_players, discard_players, substitute_players, balance

contains


    function hyper2_matrix(mask, powers, names) result(h)
        logical, intent(in) :: mask(:,:)
        real(dp), intent(in) :: powers(:)
        character(len=*), intent(in) :: names(:)
        type(hyper2_model) :: h
        character(len=name_len), allocatable :: bracket(:)
        integer :: i, j, n

        if (size(mask,2) /= size(names) .or. size(mask,1) /= size(powers)) then
            error stop "hyper2_matrix: size mismatch"
        end if
        h = hyper2_create(names)
        do i = 1, size(mask,1)
            n = count(mask(i,:))
            if (n == 0 .or. abs(powers(i)) <= 0.0_dp) cycle
            allocate(bracket(n))
            n = 0
            do j = 1, size(names)
                if (mask(i,j)) then
                    n = n + 1
                    bracket(n) = names(j)
                end if
            end do
            call h%add_term(bracket, powers(i))
            deallocate(bracket)
        end do
    end function hyper2_matrix

    function hyper3_matrix(weights, powers, names, stripzeros) result(h)
        real(dp), intent(in) :: weights(:,:), powers(:)
        character(len=*), intent(in) :: names(:)
        logical, intent(in), optional :: stripzeros
        type(hyper3_model) :: h
        character(len=name_len), allocatable :: bracket(:)
        real(dp), allocatable :: w(:)
        logical :: strip
        integer :: i, j, n

        if (size(weights,2) /= size(names) .or. size(weights,1) /= size(powers)) then
            error stop "hyper3_matrix: size mismatch"
        end if
        if (any(weights < 0.0_dp)) error stop "hyper3_matrix: negative weights"
        strip = .true.
        if (present(stripzeros)) strip = stripzeros
        h = hyper3_create(names)
        do i = 1, size(weights,1)
            n = count(weights(i,:) > 0.0_dp)
            if (strip .and. n == 0) cycle
            if (n == 0 .or. abs(powers(i)) <= 0.0_dp) cycle
            allocate(bracket(n),w(n))
            n = 0
            do j = 1, size(names)
                if (weights(i,j) > 0.0_dp) then
                    n = n + 1
                    bracket(n) = names(j)
                    w(n) = weights(i,j)
                end if
            end do
            call h%add_term(bracket,w,powers(i))
            deallocate(bracket,w)
        end do
    end function hyper3_matrix

    function zipf(n) result(p)
        integer, intent(in) :: n
        real(dp), allocatable :: p(:)
        integer :: i
        allocate(p(max(0,n)))
        do i = 1, n
            p(i) = 1.0_dp / real(i, dp)
        end do
        if (n > 0) p = p / sum(p)
    end function zipf

    function dirichlet(names, powers, alpha) result(h)
        character(len=*), intent(in) :: names(:)
        real(dp), intent(in), optional :: powers(:), alpha(:)
        type(hyper2_model) :: h
        real(dp), allocatable :: pw(:)
        integer :: i
        if (present(powers) .eqv. present(alpha)) error stop "dirichlet: supply exactly one of powers or alpha"
        if (present(powers)) then
            if (size(powers) /= size(names)) error stop "dirichlet size mismatch"
            pw = powers
        else
            if (size(alpha) /= size(names)) error stop "dirichlet size mismatch"
            pw = alpha - 1.0_dp
        end if
        h = hyper2_create(names)
        do i = 1, size(names)
            call h%add_term(names(i:i), pw(i))
        end do
        if (size(names) > 0) call h%add_term(names, -sum(pw))
    end function dirichlet

    function generalized_dirichlet(names, alpha, beta, beta0) result(h)
        character(len=*), intent(in) :: names(:)
        real(dp), intent(in) :: alpha(:), beta(:)
        real(dp), intent(in), optional :: beta0
        type(hyper2_model) :: h
        real(dp) :: b0
        integer :: k, i
        k = size(alpha)
        if (size(names) /= k .or. size(beta) /= k - 1) error stop "GD size mismatch"
        b0 = 0.0_dp
        if (present(beta0)) b0 = beta0
        h = dirichlet(names, powers=alpha - 1.0_dp)
        if (k >= 3) then
            do i = 2, k - 1
                call h%set_term(names(i:k), beta(i-1) - (alpha(i) + beta(i)))
            end do
        end if
        if (k >= 2) call h%set_term(names(k:k), beta(k-1) - 1.0_dp)
        if (k >= 2) call h%set_term(names, b0 - (alpha(1) + beta(1)))
    end function generalized_dirichlet

    function gd_wong(names, alpha, beta) result(h)
        character(len=*), intent(in) :: names(:)
        real(dp), intent(in) :: alpha(:), beta(:)
        type(hyper2_model) :: h
        real(dp), allocatable :: gamma(:)
        integer :: k, i
        k = size(alpha)
        if (size(names) /= k .or. size(beta) /= k) error stop "GD_wong size mismatch"
        allocate(gamma(k))
        if (k > 1) gamma(1:k-1) = beta(1:k-1) - alpha(2:k) - beta(2:k)
        gamma(k) = beta(k) - 1.0_dp
        h = dirichlet(names, powers=alpha - 1.0_dp)
        call h%set_term(names(k:k), gamma(k))
        do i = 1, k - 1
            call h%set_term(names(i+1:k), gamma(i))
        end do
    end function gd_wong

    function dirichlet3(names, powers, lambda) result(h)
        character(len=*), intent(in) :: names(:)
        real(dp), intent(in) :: powers(:)
        real(dp), intent(in), optional :: lambda(:)
        type(hyper3_model) :: h
        real(dp), allocatable :: lam(:)
        integer :: i
        if (size(names) /= size(powers)) error stop "dirichlet3 size mismatch"
        allocate(lam(size(names)))
        lam = 1.0_dp
        if (present(lambda)) then
            if (size(lambda) > size(names)) error stop "dirichlet3 lambda too long"
            lam(1:size(lambda)) = lambda
        end if
        h = hyper3_create(names)
        do i = 1, size(names)
            call h%add_term(names(i:i), lam(i:i), powers(i))
        end do
        if (size(names) > 0) call h%add_term(names, lam, -sum(powers))
    end function dirichlet3

    function trial(winners, players, value) result(h)
        character(len=*), intent(in) :: winners(:), players(:)
        real(dp), intent(in), optional :: value
        type(hyper2_model) :: h
        real(dp) :: v
        v = 1.0_dp
        if (present(value)) v = value
        h = hyper2_create(players)
        call h%add_term(winners, v)
        call h%add_term(players, -v)
    end function trial

    function pick(winners, allplayers) result(h)
        character(len=*), intent(in) :: winners(:), allplayers(:)
        type(hyper2_model) :: h
        h = trial(winners, allplayers)
    end function pick

    function pass_model(losers, allplayers) result(h)
        character(len=*), intent(in) :: losers(:), allplayers(:)
        type(hyper2_model) :: h
        character(len=name_len), allocatable :: winners(:)
        integer :: i, n
        allocate(winners(size(allplayers)))
        n = 0
        do i = 1, size(allplayers)
            if (.not. contains_name(losers, allplayers(i))) then
                n = n + 1
                winners(n) = allplayers(i)
            end if
        end do
        h = trial(winners(1:n), allplayers)
    end function pass_model

    function beats(winners, losers) result(h)
        character(len=*), intent(in) :: winners(:), losers(:)
        type(hyper2_model) :: h
        character(len=name_len), allocatable :: allp(:)
        integer :: n
        n = size(winners) + size(losers)
        allocate(allp(n))
        if (size(winners)>0) allp(1:size(winners)) = winners
        if (size(losers)>0) allp(size(winners)+1:n) = losers
        h = trial(winners, allp)
    end function beats

    logical function contains_name(names, name)
        character(len=*), intent(in) :: names(:), name
        integer :: i
        contains_name = .false.
        do i=1,size(names)
            if(trim(names(i))==trim(name)) then
                contains_name=.true.
                return
            end if
        end do
    end function contains_name

    function unique_names(v, extra) result(names)
        character(len=*), intent(in) :: v(:)
        character(len=*), intent(in), optional :: extra(:)
        character(len=name_len), allocatable :: names(:), tmp(:)
        integer :: i, n
        allocate(names(0))
        do i=1,size(v)
            if(.not.contains_name(names,v(i))) then
                n=size(names)
                allocate(tmp(n+1))
                if(n>0) tmp(1:n)=names
                tmp(n+1)=v(i)
                call move_alloc(tmp,names)
            end if
        end do
        if(present(extra)) then
            do i=1,size(extra)
                if(.not.contains_name(names,extra(i))) then
                    n=size(names)
                    allocate(tmp(n+1))
                    if(n>0) tmp(1:n)=names
                    tmp(n+1)=extra(i)
                    call move_alloc(tmp,names)
                end if
            end do
        end if
    end function unique_names

    function rankvec_likelihood(v, nonfinishers) result(h)
        character(len=*), intent(in) :: v(:)
        character(len=*), intent(in), optional :: nonfinishers(:)
        type(hyper2_model) :: h
        character(len=name_len), allocatable :: names(:), denom(:)
        integer :: i, n, m
        if (present(nonfinishers)) then
            names = unique_names(v, nonfinishers)
        else
            names = unique_names(v)
        end if
        if (size(names) < size(v)) error stop "rankvec_likelihood: repeated finisher"
        h = hyper2_create(names)
        do i = 1, size(v)
            call h%add_term(v(i:i), 1.0_dp)
            n = size(v) - i + 1
            m = 0
            if (present(nonfinishers)) m = size(nonfinishers)
            allocate(denom(n+m))
            denom(1:n)=v(i:size(v))
            if(m>0) denom(n+1:n+m)=nonfinishers
            call h%add_term(denom, -1.0_dp)
            deallocate(denom)
        end do
    end function rankvec_likelihood

    function ordervec2supp(v, nonfinishers) result(h)
        character(len=*), intent(in) :: v(:)
        character(len=*), intent(in), optional :: nonfinishers(:)
        type(hyper2_model) :: h
        if (present(nonfinishers)) then
            h = rankvec_likelihood(v, nonfinishers)
        else
            h = rankvec_likelihood(v)
        end if
    end function ordervec2supp

    function ordervec2supp3(v, nonfinishers) result(h)
        character(len=*), intent(in) :: v(:)
        character(len=*), intent(in), optional :: nonfinishers(:)
        type(hyper3_model) :: h
        character(len=name_len), allocatable :: names(:), denom_names(:), un(:)
        real(dp), allocatable :: denom_weights(:)
        integer :: i, j, n, m, cnt
        if (present(nonfinishers)) then
            names = unique_names(v, nonfinishers)
        else
            names = unique_names(v)
        end if
        h = hyper3_create(names)
        un = unique_names(v)
        do i=1,size(un)
            cnt=0
            do j=1,size(v)
                if(trim(v(j))==trim(un(i))) cnt=cnt+1
            end do
            call h%add_term(un(i:i), [1.0_dp], real(cnt,dp))
        end do
        do i=1,size(v)
            n=size(v)-i+1
            m=0
            if(present(nonfinishers)) m=size(nonfinishers)
            block
                character(len=name_len), allocatable :: raw(:)
                allocate(raw(n+m))
                raw(1:n)=v(i:size(v))
                if(m>0) raw(n+1:n+m)=nonfinishers
                denom_names=unique_names(raw)
                allocate(denom_weights(size(denom_names)))
                do j=1,size(denom_names)
                    cnt=0
                    do m=1,size(raw)
                        if(trim(raw(m))==trim(denom_names(j))) cnt=cnt+1
                    end do
                    denom_weights(j)=real(cnt,dp)
                end do
            end block
            call h%add_term(denom_names,denom_weights,-1.0_dp)
            deallocate(denom_names,denom_weights)
        end do
    end function ordervec2supp3

    function ordervec2supp3a(v, helped, lambda, nonfinishers) result(h)
        character(len=*), intent(in) :: v(:), helped(:)
        real(dp), intent(in) :: lambda
        character(len=*), intent(in), optional :: nonfinishers(:)
        type(hyper3_model) :: h
        character(len=name_len), allocatable :: names(:), dn(:), raw(:), un(:)
        real(dp), allocatable :: w(:)
        integer :: i,j,k,cnt,n,m
        if(present(nonfinishers)) then
            names=unique_names(v,nonfinishers)
        else
            names=unique_names(v)
        end if
        h=hyper3_create(names)
        un=unique_names(v)
        do i=1,size(un)
            cnt = 0
            do j=1,size(v)
                if(trim(v(j))==trim(un(i))) cnt=cnt+1
            end do
            if(contains_name(helped,un(i))) then
                call h%add_term(un(i:i),[lambda],real(cnt,dp))
            else
                call h%add_term(un(i:i),[1.0_dp],real(cnt,dp))
            end if
        end do
        do i=1,size(v)
            n=size(v)-i+1
            m=0
            if(present(nonfinishers)) m=size(nonfinishers)
            allocate(raw(n+m))
            raw(1:n)=v(i:size(v))
            if(m>0)raw(n+1:n+m)=nonfinishers
            dn=unique_names(raw)
            allocate(w(size(dn)))
            do j=1,size(dn)
                cnt=0
                do k=1,size(raw)
                    if(trim(raw(k))==trim(dn(j))) cnt=cnt+1
                end do
                w(j)=real(cnt,dp)
                if(contains_name(helped,dn(j))) w(j)=lambda*w(j)
            end do
            call h%add_term(dn,w,-1.0_dp)
            deallocate(raw,dn,w)
        end do
    end function ordervec2supp3a

    function pairwise(wins, names) result(h)
        real(dp), intent(in) :: wins(:,:)
        character(len=*), intent(in) :: names(:)
        type(hyper2_model) :: h
        integer :: i,j,n
        real(dp) :: t, rt
        n=size(names)
        if(size(wins,1)/=n .or. size(wins,2)/=n) error stop "pairwise size mismatch"
        h=hyper2_create(names)
        do i=1,n
            rt = sum(wins(i,:))
            if(abs(rt) > 0.0_dp) call h%add_term(names(i:i),rt)
        end do
        do i=1,n-1
            do j=i+1,n
                t=wins(i,j)+wins(j,i)
                if(abs(t) > 0.0_dp) call h%add_term([names(i),names(j)],-t)
            end do
        end do
    end function pairwise

    function home_away(home_wins, away_wins, teams, monster) result(h)
        real(dp), intent(in) :: home_wins(:,:), away_wins(:,:)
        character(len=*), intent(in) :: teams(:), monster
        type(hyper2_model) :: h
        character(len=name_len), allocatable :: names(:)
        integer :: i,j,n
        real(dp)::hw,aw
        n=size(teams)
        allocate(names(n+1))
        names(1:n)=teams
        names(n+1)=monster
        h=hyper2_create(names)
        do i=1,n
        do j=1,n
            if(i/=j) then
                hw=home_wins(i,j)
                aw=away_wins(i,j)
                call h%add_term([teams(i),monster],hw)
                call h%add_term(teams(j:j),aw)
                call h%add_term([teams(i),teams(j),monster],-(hw+aw))
            end if
        end do
        end do
    end function home_away

    function home_away3(home_wins,away_wins,teams,lambda) result(h)
        real(dp),intent(in)::home_wins(:,:),away_wins(:,:),lambda
        character(len=*),intent(in)::teams(:)
        type(hyper3_model)::h
        integer::i,j,n
        real(dp)::hw,aw
        n=size(teams)
        h=hyper3_create(teams)
        do i=1,n
        do j=1,n
            if(i/=j)then
                hw=home_wins(i,j)
                aw=away_wins(i,j)
                call h%add_term(teams(i:i),[lambda],hw)
                call h%add_term(teams(j:j),[1.0_dp],aw)
                call h%add_term([teams(i),teams(j)],[lambda,1.0_dp],-(hw+aw))
            end if
        end do
        end do
    end function home_away3

    function home_draw_away3(home_wins,draws,away_wins,teams,lambda,d) result(h)
        real(dp),intent(in)::home_wins(:,:),draws(:,:),away_wins(:,:),lambda,d
        character(len=*),intent(in)::teams(:)
        type(hyper3_model)::h
        integer::i,j,n
        real(dp)::hw,aw,dr,t
        n=size(teams)
        h=hyper3_create(teams)
        do i=1,n
        do j=1,n
            if(i/=j)then
                hw=home_wins(i,j)
                aw=away_wins(i,j)
                dr=draws(i,j)
                t=hw+aw+dr
                call h%add_term(teams(i:i),[lambda],hw)
                call h%add_term(teams(j:j),[1.0_dp],aw)
                call h%add_term([teams(i),teams(j)],[d,d],dr)
                call h%add_term([teams(i),teams(j)],[d+lambda,d+1.0_dp],-t)
            end if
        end do
        end do
    end function home_draw_away3

    function zermelo(m_in,maxit,start,tol) result(p)
        real(dp),intent(in)::m_in(:,:)
        integer,intent(in),optional::maxit
        real(dp),intent(in),optional::start(:),tol
        real(dp),allocatable::p(:),pn(:),m(:,:),r(:)
        real(dp)::tt,den
        integer::n,it,mi,i,j
        n=size(m_in,1)
        if(size(m_in,2)/=n) error stop "zermelo square matrix required"
        mi=100
        if(present(maxit))mi=maxit
        tt=0.0_dp
        if(present(tol))tt=tol
        m=m_in
        do i=1,n
        m(i,i)=0.0_dp
        end do
        allocate(r(n))
        r=sum(m,dim=2)
        m=m+transpose(m)
        allocate(p(n),pn(n))
        if(present(start))then
        p=start
        else
        p=1.0_dp/real(n,dp)
        end if
        do it=1,mi
            do i=1,n
                den=0.0_dp
                do j=1,n
                    if(m(j,i)>0.0_dp .and. p(j)+p(i)>0.0_dp) den=den+m(j,i)/(p(j)+p(i))
                end do
                if(den>0.0_dp)then
                pn(i)=r(i)/den
                else
                pn(i)=0.0_dp
                end if
            end do
            if(sum(pn)>0.0_dp)pn=pn/sum(pn)
            if(all(abs(p-pn)<=tt))exit
            p=pn
        end do
        p=pn
    end function zermelo

    function as_hyper3(h2) result(h3)
        type(hyper2_model),intent(in)::h2
        type(hyper3_model)::h3
        character(len=name_len),allocatable::bn(:)
        real(dp),allocatable::w(:)
        integer::i,j
        h3=hyper3_create(h2%pnames)
        do i=1,size(h2%terms)
            allocate(bn(size(h2%terms(i)%ids)),w(size(h2%terms(i)%ids)))
            w=1.0_dp
            do j=1,size(bn)
            bn(j)=h2%pnames(h2%terms(i)%ids(j))
            end do
            call h3%add_term(bn,w,h2%terms(i)%power)
            deallocate(bn,w)
        end do
    end function as_hyper3

    function hyper3_to_hyper2(h3) result(h2)
        type(hyper3_model),intent(in)::h3
        type(hyper2_model)::h2
        character(len=name_len),allocatable::bn(:)
        integer::i,j
        h2=hyper2_create(h3%pnames)
        do i=1,size(h3%terms)
            allocate(bn(size(h3%terms(i)%ids)))
            do j=1,size(bn)
            bn(j)=h3%pnames(h3%terms(i)%ids(j))
            end do
            call h2%add_term(bn,h3%terms(i)%power)
            deallocate(bn)
        end do
    end function hyper3_to_hyper2

    subroutine setweight(h,players,values)
        type(hyper3_model),intent(inout)::h
        character(len=*),intent(in)::players(:)
        real(dp),intent(in)::values(:)
        integer::i,j,k,id
        if(size(players)/=size(values))error stop "setweight size mismatch"
        do k=1,size(players)
            id=player_index(h%pnames,players(k))
            if(id==0)cycle
            do i=1,size(h%terms)
                do j=1,size(h%terms(i)%ids)
                    if(h%terms(i)%ids(j)==id)h%terms(i)%weights(j)=values(k)
                end do
            end do
        end do
    end subroutine setweight

    function pwa3(h_in,players,lambda) result(h)
        type(hyper3_model),intent(in)::h_in
        character(len=*),intent(in)::players(:)
        real(dp),intent(in)::lambda(:)
        type(hyper3_model)::h
        h=h_in
        call setweight(h,players,lambda)
    end function pwa3

    function pwa23(h2,player,chameleon) result(h3)
        type(hyper2_model),intent(in)::h2
        character(len=*),intent(in)::player,chameleon
        type(hyper3_model)::h3
        character(len=name_len),allocatable::bn(:),names(:)
        real(dp),allocatable::w(:)
        integer::i,j,n,id
        names=h2%pnames
        if(.not.contains_name(names,chameleon))then
            n=size(names)
            block
                character(len=name_len),allocatable::tmp(:)
                allocate(tmp(n+1))
                if(n>0)tmp(1:n)=names
                tmp(n+1)=chameleon
                call move_alloc(tmp,names)
            end block
        end if
        h3=hyper3_create(names)
        id=player_index(h2%pnames,player)
        do i=1,size(h2%terms)
            n=size(h2%terms(i)%ids)
            if(any(h2%terms(i)%ids==id))then
                allocate(bn(n+1),w(n+1))
                w=1.0_dp
                do j=1,n
                bn(j)=h2%pnames(h2%terms(i)%ids(j))
                end do
                bn(n+1)=chameleon
            else
                allocate(bn(n),w(n))
                w=1.0_dp
                do j=1,n
                bn(j)=h2%pnames(h2%terms(i)%ids(j))
                end do
            end if
            call h3%add_term(bn,w,h2%terms(i)%power)
            deallocate(bn,w)
        end do
    end function pwa23

    function pwa(h,player,chameleon) result(out)
        type(hyper2_model),intent(in)::h
        character(len=*),intent(in)::player,chameleon
        type(hyper2_model)::out
        character(len=name_len),allocatable::bn(:)
        integer::i,j,n,id
        out=hyper2_create()
        id=player_index(h%pnames,player)
        do i=1,size(h%terms)
            n=size(h%terms(i)%ids)
            if(any(h%terms(i)%ids==id))then
                allocate(bn(n+1))
                do j=1,n
                bn(j)=h%pnames(h%terms(i)%ids(j))
                end do
                bn(n+1)=chameleon
            else
                allocate(bn(n))
                do j=1,n
                bn(j)=h%pnames(h%terms(i)%ids(j))
                end do
            end if
            call out%add_term(bn,h%terms(i)%power)
            deallocate(bn)
        end do
    end function pwa

    function keep_players(h,wanted) result(out)
        type(hyper2_model),intent(in)::h
        character(len=*),intent(in)::wanted(:)
        type(hyper2_model)::out
        character(len=name_len),allocatable::bn(:)
        integer::i,j,n
        out=hyper2_create(h%pnames)
        do i=1,size(h%terms)
            allocate(bn(size(h%terms(i)%ids)))
            n=0
            do j=1,size(h%terms(i)%ids)
                if(contains_name(wanted,h%pnames(h%terms(i)%ids(j))))then
                    n=n+1
                    bn(n)=h%pnames(h%terms(i)%ids(j))
                end if
            end do
            if(n>0)call out%add_term(bn(1:n),h%terms(i)%power)
            deallocate(bn)
        end do
    end function keep_players

    function discard_players(h,unwanted) result(out)
        type(hyper2_model),intent(in)::h
        character(len=*),intent(in)::unwanted(:)
        type(hyper2_model)::out
        character(len=name_len),allocatable::wanted(:)
        integer::i,n
        allocate(wanted(size(h%pnames)))
        n=0
        do i=1,size(h%pnames)
            if(.not.contains_name(unwanted,h%pnames(i)))then
            n=n+1
            wanted(n)=h%pnames(i)
            end if
        end do
        out=keep_players(h,wanted(1:n))
    end function discard_players

    function substitute_players(h,from,to) result(out)
        type(hyper2_model),intent(in)::h
        character(len=*),intent(in)::from(:),to(:)
        type(hyper2_model)::out
        character(len=name_len),allocatable::bn(:)
        integer::i,j,k
        if(size(from)/=size(to))error stop "substitute_players size mismatch"
        out=hyper2_create()
        do i=1,size(h%terms)
            allocate(bn(size(h%terms(i)%ids)))
            do j=1,size(bn)
                bn(j)=h%pnames(h%terms(i)%ids(j))
                do k=1,size(from)
                    if(trim(bn(j))==trim(from(k)))bn(j)=to(k)
                end do
            end do
            call out%add_term(bn,h%terms(i)%power)
            deallocate(bn)
        end do
    end function substitute_players

    function balance(h) result(out)
        type(hyper2_model),intent(in)::h
        type(hyper2_model)::out
        real(dp)::s
        integer::i
        out=h
        do i=1,size(out%pnames)
        call out%set_term(out%pnames(i:i),0.0_dp)
        end do
        s=0.0_dp
        do i=1,size(out%terms)
        s=s+out%terms(i)%power
        end do
        call out%set_term(out%pnames,-s)
    end function balance

end module hyper2_models
