! coda-fortran: computational translation of the R package coda.
! Original coda license: GPL (>= 2). This translation is GPL-2.0-or-later.
module coda_types
   use coda_kinds, only : dp
   implicit none
   private

   type, public :: mcmc_chain
      real(dp), allocatable :: x(:,:)
      integer :: start = 1
      integer :: finish = 0
      integer :: thin = 1
      character(len=:), allocatable :: var_names(:)
   contains
      procedure :: niter => chain_niter
      procedure :: nvar => chain_nvar
   end type mcmc_chain

   type, public :: mcmc_list
      type(mcmc_chain), allocatable :: chain(:)
   contains
      procedure :: nchain => list_nchain
      procedure :: niter => list_niter
      procedure :: nvar => list_nvar
   end type mcmc_list

   public :: make_mcmc, make_mcmc_list, window_mcmc, pool_chains

contains

   function make_mcmc(data, start, finish, thin, var_names) result(chain)
      real(dp), intent(in) :: data(:,:)
      integer, intent(in), optional :: start, finish, thin
      character(len=*), intent(in), optional :: var_names(:)
      type(mcmc_chain) :: chain
      integer :: s, e, t, j, maxlen, nobs, ndata

      ndata = size(data,1)
      t = 1
      if (present(thin)) t = thin
      if (t <= 0) error stop "make_mcmc: thin must be positive"

      if (present(start)) then
         s = start
      else if (present(finish)) then
         s = finish - (ndata - 1) * t
      else
         s = 1
      end if

      if (present(finish)) then
         e = finish
         nobs = floor(real(e-s,dp)/real(t,dp) + 1.0_dp)
         if (nobs < 1 .or. ndata < nobs) error stop "make_mcmc: start, finish and thin incompatible with data"
      else
         nobs = ndata
      end if
      e = s + (nobs - 1) * t

      allocate(chain%x(nobs,size(data,2)))
      chain%x = data(1:nobs,:)
      chain%start = s
      chain%thin = t
      chain%finish = e

      if (present(var_names)) then
         if (size(var_names) /= size(data, 2)) error stop "make_mcmc: wrong number of variable names"
         maxlen = 1
         do j = 1, size(var_names)
            maxlen = max(maxlen, len_trim(var_names(j)))
         end do
         allocate(character(len=maxlen) :: chain%var_names(size(var_names)))
         chain%var_names = var_names
      end if
   end function make_mcmc

   function make_mcmc_list(chains) result(lst)
      type(mcmc_chain), intent(in) :: chains(:)
      type(mcmc_list) :: lst
      integer :: i

      if (size(chains) < 1) error stop "make_mcmc_list: at least one chain is required"
      do i = 2, size(chains)
         if (chains(i)%niter() /= chains(1)%niter()) error stop "make_mcmc_list: niter mismatch"
         if (chains(i)%nvar() /= chains(1)%nvar()) error stop "make_mcmc_list: nvar mismatch"
         if (chains(i)%start /= chains(1)%start) error stop "make_mcmc_list: start mismatch"
         if (chains(i)%finish /= chains(1)%finish) error stop "make_mcmc_list: end mismatch"
         if (chains(i)%thin /= chains(1)%thin) error stop "make_mcmc_list: thin mismatch"
      end do
      lst%chain = chains
   end function make_mcmc_list

   integer function chain_niter(self) result(n)
      class(mcmc_chain), intent(in) :: self
      if (allocated(self%x)) then
         n = size(self%x, 1)
      else
         n = 0
      end if
   end function chain_niter

   integer function chain_nvar(self) result(n)
      class(mcmc_chain), intent(in) :: self
      if (allocated(self%x)) then
         n = size(self%x, 2)
      else
         n = 0
      end if
   end function chain_nvar

   integer function list_nchain(self) result(n)
      class(mcmc_list), intent(in) :: self
      if (allocated(self%chain)) then
         n = size(self%chain)
      else
         n = 0
      end if
   end function list_nchain

   integer function list_niter(self) result(n)
      class(mcmc_list), intent(in) :: self
      if (allocated(self%chain) .and. size(self%chain) > 0) then
         n = self%chain(1)%niter()
      else
         n = 0
      end if
   end function list_niter

   integer function list_nvar(self) result(n)
      class(mcmc_list), intent(in) :: self
      if (allocated(self%chain) .and. size(self%chain) > 0) then
         n = self%chain(1)%nvar()
      else
         n = 0
      end if
   end function list_nvar

   function window_mcmc(chain, start, finish, thin) result(out)
      type(mcmc_chain), intent(in) :: chain
      integer, intent(in), optional :: start, finish, thin
      type(mcmc_chain) :: out
      integer :: s, e, t, first, last, stride, n, k, idx

      s = chain%start
      if (present(start)) s = max(start, chain%start)
      e = chain%finish
      if (present(finish)) e = min(finish, chain%finish)
      t = chain%thin
      if (present(thin)) then
         if (mod(thin, chain%thin) == 0) t = thin
      end if
      if (s > e) error stop "window_mcmc: start after finish"

      first = ceiling(real(s - chain%start, dp) / real(chain%thin, dp)) + 1
      first = max(1, first)
      last = floor(real(e - chain%start, dp) / real(chain%thin, dp)) + 1
      last = min(chain%niter(), last)
      stride = max(1, t / chain%thin)
      n = 1 + (last - first) / stride
      allocate(out%x(n, chain%nvar()))
      do k = 1, n
         idx = first + (k - 1) * stride
         out%x(k,:) = chain%x(idx,:)
      end do
      out%start = chain%start + (first - 1) * chain%thin
      out%thin = t
      out%finish = out%start + (n - 1) * t
      if (allocated(chain%var_names)) out%var_names = chain%var_names
   end function window_mcmc

   function pool_chains(lst) result(x)
      type(mcmc_list), intent(in) :: lst
      real(dp), allocatable :: x(:,:)
      integer :: i, n, p

      n = lst%niter()
      p = lst%nvar()
      allocate(x(n * lst%nchain(), p))
      do i = 1, lst%nchain()
         x((i - 1) * n + 1:i * n,:) = lst%chain(i)%x
      end do
   end function pool_chains

end module coda_types
