program test_random_fmatch
   use optmatch_kinds, only : dp
   use optmatch_types, only : distance_spec, match_result
   use optmatch_matching, only : fmatch_core
   implicit none
   integer :: trial, nt, nc, i, j
   integer, allocatable :: seed(:)
   real(dp) :: r, mn, mx, brute_cost, omf
   logical :: brute_feasible, has_omit
   type(distance_spec) :: d
   type(match_result) :: m

   call random_seed(size=i)
   allocate(seed(i))
   seed = 918273 + 37 * [(j, j=1,i)]
   call random_seed(put=seed)

   do trial = 1, 400
      nt = 1 + mod(trial, 3)
      nc = 1 + mod(trial / 3, 3)
      allocate(d%value(nt,nc), d%allowed(nt,nc))
      do j = 1, nc
         do i = 1, nt
            call random_number(r)
            d%allowed(i,j) = r > 0.20_dp
            call random_number(r)
            d%value(i,j) = 0.1_dp + 9.9_dp*r
         end do
      end do
      select case(mod(trial, 4))
      case(0)
         mn = 1.0_dp / real(nt,dp); mx = real(nc,dp)
      case(1)
         mn = 1.0_dp; mx = real(max(1,nc),dp)
      case(2)
         mn = 1.0_dp / real(min(2,nt),dp); mx = 1.0_dp
      case default
         mn = 1.0_dp; mx = 1.0_dp
      end select
      has_omit = mod(trial, 3) == 0
      if (has_omit) then
         omf = real(mod(trial, nc + 1), dp) / real(nc, dp)
      else
         omf = 0.0_dp
      end if
      m = fmatch_core(d, mn, mx, has_omit, omf)
      call brute_fmatch(d, mn, mx, has_omit, omf, brute_feasible, brute_cost)
      if (m%feasible .neqv. brute_feasible) then
         write(*,'(a,i0)') 'feasibility mismatch trial ', trial
         error stop 1
      end if
      if (m%feasible) then
         if (abs(m%objective-brute_cost) > 1.0e-8_dp) then
            write(*,'(a,i0,2(1x,es18.10))') 'objective mismatch trial ', trial, m%objective, brute_cost
            error stop 1
         end if
      end if
      deallocate(d%value,d%allowed)
   end do
   print '(a)', 'Randomized fmatch brute-force validation passed (400 cases).'

contains

subroutine brute_fmatch(d, min_cpt, max_cpt, has_omit, omit_fraction, feasible, best)
   type(distance_spec), intent(in) :: d
   real(dp), intent(in) :: min_cpt, max_cpt, omit_fraction
   logical, intent(in) :: has_omit
   logical, intent(out) :: feasible
   real(dp), intent(out) :: best
   integer :: nt, nc, mxc, mnc, mxr, nmc, ne, mask, k, i, j
   integer, allocatable :: ei(:), ej(:), rd(:), cd(:)
   real(dp) :: cost
   logical :: ok
   nt=size(d%value,1); nc=size(d%value,2)
   mxc=ceiling(max_cpt); mnc=max(1,floor(min_cpt)); mxr=ceiling(1.0_dp/min_cpt)
   if(mnc>1) mxr=1
   if (has_omit) then
      nmc = nint(real(nc,dp) * (1.0_dp-omit_fraction))
   else
      nmc = nc
   end if
   ne=count(d%allowed)
   allocate(ei(ne),ej(ne),rd(nt),cd(nc))
   k=0
   do j=1,nc
      do i=1,nt
         if(.not.d%allowed(i,j)) cycle
         k=k+1; ei(k)=i; ej(k)=j
      end do
   end do
   feasible=.false.; best=huge(1.0_dp)
   if(ne>20) error stop 'brute test too large'
   do mask=0,2**ne-1
      rd=0; cd=0; cost=0.0_dp
      do k=1,ne
         if(.not.btest(mask,k-1)) cycle
         rd(ei(k))=rd(ei(k))+1
         cd(ej(k))=cd(ej(k))+1
         cost=cost+d%value(ei(k),ej(k))
      end do
      ok=all(rd>=mnc .and. rd<=mxc)
      if(.not.ok) cycle
      if(any(cd>mxr)) cycle
      ! There must exist sink flows s_j in {0,1}, with s_j<=cd_j,
      ! cd_j-s_j <= mxr-1, and sum(s_j)=nmc.
      if(count(cd==mxr)>nmc) cycle
      if(count(cd>=1)<nmc) cycle
      feasible=.true.
      best=min(best,cost)
   end do
end subroutine brute_fmatch

end program test_random_fmatch
