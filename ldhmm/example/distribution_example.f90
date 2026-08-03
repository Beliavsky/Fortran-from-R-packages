! SPDX-License-Identifier: Artistic-2.0
program distribution_example
   use ldhmm
   implicit none
   type(ecld_type) :: distribution
   real(dp) :: x

   distribution = ecld_create(lambda=3.0_dp,sigma=0.02_dp,mu=0.001_dp)
   x = 0.03_dp
   write(*,'(a,es14.6)') 'pdf:      ', ecld_pdf(distribution,x)
   write(*,'(a,es14.6)') 'cdf:      ', ecld_cdf(distribution,x)
   write(*,'(a,es14.6)') 'mean:     ', ecld_mean(distribution)
   write(*,'(a,es14.6)') 'sd:       ', ecld_sd(distribution)
   write(*,'(a,es14.6)') 'kurtosis: ', ecld_kurtosis(distribution)
end program distribution_example
