program fit_csv
   use rugarch
   implicit none
   character(len=1024) :: filename, line
   real(dp), allocatable :: x(:), tmp(:)
   type(garch_fit_result) :: fit
   integer :: unit, ios, n, capacity
   real(dp) :: value

   if (command_argument_count()<1) then
      write(*,'(a)') 'usage: fpm run --example fit_csv -- FILE.csv'
      stop 1
   end if
   call get_command_argument(1,filename)
   capacity=1024
   allocate(x(capacity))
   n=0
   open(newunit=unit,file=trim(filename),status='old',action='read',iostat=ios)
   if (ios/=0) error stop 'cannot open input file'
   do
      read(unit,'(a)',iostat=ios) line
      if (ios/=0) exit
      read(line,*,iostat=ios) value
      if (ios/=0) cycle
      if (n==capacity) then
         allocate(tmp(2*capacity))
         tmp(1:n)=x(1:n)
         call move_alloc(tmp,x)
         capacity=size(x)
      end if
      n=n+1
      x(n)=value
   end do
   close(unit)
   if (n<30) error stop 'need at least 30 numeric observations'

   fit=fit_garch11(x(1:n),dist_std,fit_shape=.true.,max_iterations=1500)
   print '(a,i0)', 'observations: ', n
   print '(a,f12.6)', 'mean:          ', fit%spec%mean
   print '(a,es14.6)', 'omega:         ', fit%spec%omega
   print '(a,f12.6)', 'alpha:         ', fit%spec%alpha(1)
   print '(a,f12.6)', 'beta:          ', fit%spec%beta(1)
   print '(a,f12.6)', 'shape:         ', fit%spec%shape
   print '(a,f14.3)', 'log likelihood:', fit%log_likelihood
   print '(a,f14.3)', 'AIC:           ', fit%aic
end program fit_csv
