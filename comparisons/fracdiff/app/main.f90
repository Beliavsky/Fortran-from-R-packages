program compare_fracdiff
 use fracdiff, only: dp,diffseries,fd_gph,fd_sperio,fractional_d_estimate
 implicit none
 integer,parameter::n=4096; real(dp)::x(n),dx(n),v,t0,t1; integer::i,j,u,st
 type(fractional_d_estimate)::est; character(512)::out
 call get_command_argument(1,out); if(len_trim(out)==0)out='fortran_results.csv'
 open(newunit=u,file=trim(out),status='replace'); write(u,'(a)')'case,value,seconds,abs_tol,rel_tol'
 do i=1,n; x(i)=sin(.017_dp*i)+.35_dp*cos(.0031_dp*i)+.0002_dp*i; end do
 call cpu_time(t0); do j=1,100; call diffseries(x,.25_dp,dx,st); v=wcheck(dx); end do; call cpu_time(t1); call emit('diffseries_d025',v,t1-t0,1e-5_dp)
 call cpu_time(t0); do j=1,100; call diffseries(x,.70_dp,dx,st); v=wcheck(dx); end do; call cpu_time(t1); call emit('diffseries_d070',v,t1-t0,1e-5_dp)
 call cpu_time(t0); do j=1,30; est=fd_gph(x); v=est%d; end do; call cpu_time(t1); call emit('gph_d',v,t1-t0,1e-7_dp)
 call cpu_time(t0); do j=1,30; est=fd_sperio(x); v=est%d; end do; call cpu_time(t1); call emit('sperio_d',v,t1-t0,1e-7_dp)
 close(u)
contains
 real(dp) function wcheck(q); real(dp),intent(in)::q(:); integer::k; wcheck=sum(q*[(real(k,dp),k=1,size(q))]); end function
 subroutine emit(name,value,secs,atol); character(*),intent(in)::name; real(dp),intent(in)::value,secs,atol
  write(u,'(a,",",es24.16,",",es16.8,",",es12.4,",",es12.4)')trim(name),value,secs,atol,1e-8_dp
 end subroutine
end program
