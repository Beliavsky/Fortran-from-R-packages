module lmomco_named
   use lmomco_kinds, only : dp
   use lmomco_types, only : lmomco_params
   use lmomco_distributions, only : lmomco_pdf, lmomco_cdf, lmomco_quantile
   use lmomco_moments, only : theoretical_lmoments, fit_lmoments
   implicit none
   private
   public :: pdfexp, cdfexp, quaexp
   public :: pdfnor, cdfnor, quanor
   public :: pdfgev, cdfgev, quagev
   public :: pdfgum, cdfgum, quagum
   public :: pdfgpa, cdfgpa, quagpa
   public :: pdflap, cdflap, qualap
   public :: pdfcau, cdfcau, quacau
   public :: pdfwei, cdfwei, quawei
   public :: lmomexp
   public :: lmomnor
   public :: lmomgev
   public :: lmomgum
   public :: lmomgpa
   public :: lmomwei
   public :: parexp
   public :: parnor
   public :: pargev
   public :: pargum
contains

   real(dp) function pdfexp(x, par) result(value)
      real(dp), intent(in) :: x
      type(lmomco_params), intent(in) :: par
      value = lmomco_pdf(x, par)
   end function pdfexp

   real(dp) function cdfexp(x, par) result(value)
      real(dp), intent(in) :: x
      type(lmomco_params), intent(in) :: par
      value = lmomco_cdf(x, par)
   end function cdfexp

   real(dp) function quaexp(f, par) result(value)
      real(dp), intent(in) :: f
      type(lmomco_params), intent(in) :: par
      value = lmomco_quantile(f, par)
   end function quaexp

   real(dp) function pdfnor(x, par) result(value)
      real(dp), intent(in) :: x
      type(lmomco_params), intent(in) :: par
      value = lmomco_pdf(x, par)
   end function pdfnor

   real(dp) function cdfnor(x, par) result(value)
      real(dp), intent(in) :: x
      type(lmomco_params), intent(in) :: par
      value = lmomco_cdf(x, par)
   end function cdfnor

   real(dp) function quanor(f, par) result(value)
      real(dp), intent(in) :: f
      type(lmomco_params), intent(in) :: par
      value = lmomco_quantile(f, par)
   end function quanor

   real(dp) function pdfgev(x, par) result(value)
      real(dp), intent(in) :: x
      type(lmomco_params), intent(in) :: par
      value = lmomco_pdf(x, par)
   end function pdfgev

   real(dp) function cdfgev(x, par) result(value)
      real(dp), intent(in) :: x
      type(lmomco_params), intent(in) :: par
      value = lmomco_cdf(x, par)
   end function cdfgev

   real(dp) function quagev(f, par) result(value)
      real(dp), intent(in) :: f
      type(lmomco_params), intent(in) :: par
      value = lmomco_quantile(f, par)
   end function quagev

   real(dp) function pdfgum(x, par) result(value)
      real(dp), intent(in) :: x
      type(lmomco_params), intent(in) :: par
      value = lmomco_pdf(x, par)
   end function pdfgum

   real(dp) function cdfgum(x, par) result(value)
      real(dp), intent(in) :: x
      type(lmomco_params), intent(in) :: par
      value = lmomco_cdf(x, par)
   end function cdfgum

   real(dp) function quagum(f, par) result(value)
      real(dp), intent(in) :: f
      type(lmomco_params), intent(in) :: par
      value = lmomco_quantile(f, par)
   end function quagum

   real(dp) function pdfgpa(x, par) result(value)
      real(dp), intent(in) :: x
      type(lmomco_params), intent(in) :: par
      value = lmomco_pdf(x, par)
   end function pdfgpa

   real(dp) function cdfgpa(x, par) result(value)
      real(dp), intent(in) :: x
      type(lmomco_params), intent(in) :: par
      value = lmomco_cdf(x, par)
   end function cdfgpa

   real(dp) function quagpa(f, par) result(value)
      real(dp), intent(in) :: f
      type(lmomco_params), intent(in) :: par
      value = lmomco_quantile(f, par)
   end function quagpa

   real(dp) function pdflap(x, par) result(value)
      real(dp), intent(in) :: x
      type(lmomco_params), intent(in) :: par
      value = lmomco_pdf(x, par)
   end function pdflap

   real(dp) function cdflap(x, par) result(value)
      real(dp), intent(in) :: x
      type(lmomco_params), intent(in) :: par
      value = lmomco_cdf(x, par)
   end function cdflap

   real(dp) function qualap(f, par) result(value)
      real(dp), intent(in) :: f
      type(lmomco_params), intent(in) :: par
      value = lmomco_quantile(f, par)
   end function qualap

   real(dp) function pdfcau(x, par) result(value)
      real(dp), intent(in) :: x
      type(lmomco_params), intent(in) :: par
      value = lmomco_pdf(x, par)
   end function pdfcau

   real(dp) function cdfcau(x, par) result(value)
      real(dp), intent(in) :: x
      type(lmomco_params), intent(in) :: par
      value = lmomco_cdf(x, par)
   end function cdfcau

   real(dp) function quacau(f, par) result(value)
      real(dp), intent(in) :: f
      type(lmomco_params), intent(in) :: par
      value = lmomco_quantile(f, par)
   end function quacau

   real(dp) function pdfwei(x, par) result(value)
      real(dp), intent(in) :: x
      type(lmomco_params), intent(in) :: par
      value = lmomco_pdf(x, par)
   end function pdfwei

   real(dp) function cdfwei(x, par) result(value)
      real(dp), intent(in) :: x
      type(lmomco_params), intent(in) :: par
      value = lmomco_cdf(x, par)
   end function cdfwei

   real(dp) function quawei(f, par) result(value)
      real(dp), intent(in) :: f
      type(lmomco_params), intent(in) :: par
      value = lmomco_quantile(f, par)
   end function quawei

   subroutine lmomexp(par, nmom, lmom)
      type(lmomco_params), intent(in) :: par
      integer, intent(in) :: nmom
      real(dp), intent(out) :: lmom(nmom)
      call theoretical_lmoments(par, nmom, lmom)
   end subroutine lmomexp

   subroutine lmomnor(par, nmom, lmom)
      type(lmomco_params), intent(in) :: par
      integer, intent(in) :: nmom
      real(dp), intent(out) :: lmom(nmom)
      call theoretical_lmoments(par, nmom, lmom)
   end subroutine lmomnor

   subroutine lmomgev(par, nmom, lmom)
      type(lmomco_params), intent(in) :: par
      integer, intent(in) :: nmom
      real(dp), intent(out) :: lmom(nmom)
      call theoretical_lmoments(par, nmom, lmom)
   end subroutine lmomgev

   subroutine lmomgum(par, nmom, lmom)
      type(lmomco_params), intent(in) :: par
      integer, intent(in) :: nmom
      real(dp), intent(out) :: lmom(nmom)
      call theoretical_lmoments(par, nmom, lmom)
   end subroutine lmomgum

   subroutine lmomgpa(par, nmom, lmom)
      type(lmomco_params), intent(in) :: par
      integer, intent(in) :: nmom
      real(dp), intent(out) :: lmom(nmom)
      call theoretical_lmoments(par, nmom, lmom)
   end subroutine lmomgpa

   subroutine lmomwei(par, nmom, lmom)
      type(lmomco_params), intent(in) :: par
      integer, intent(in) :: nmom
      real(dp), intent(out) :: lmom(nmom)
      call theoretical_lmoments(par, nmom, lmom)
   end subroutine lmomwei

   function parexp(lmom) result(par)
      real(dp), intent(in) :: lmom(:)
      type(lmomco_params) :: par
      par = fit_lmoments('exp', lmom)
   end function parexp

   function parnor(lmom) result(par)
      real(dp), intent(in) :: lmom(:)
      type(lmomco_params) :: par
      par = fit_lmoments('nor', lmom)
   end function parnor

   function pargev(lmom) result(par)
      real(dp), intent(in) :: lmom(:)
      type(lmomco_params) :: par
      par = fit_lmoments('gev', lmom)
   end function pargev

   function pargum(lmom) result(par)
      real(dp), intent(in) :: lmom(:)
      type(lmomco_params) :: par
      par = fit_lmoments('gum', lmom)
   end function pargum

end module lmomco_named
