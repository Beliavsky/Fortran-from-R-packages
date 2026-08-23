program compare_cluster
 use cluster, only: dp,partition_result,hierarchy_result,silhouette_result,daisy,pam,silhouette,agnes,diana
 implicit none
 real(dp)::x(8,2),v,t0,t1; real(dp),allocatable::d(:,:); integer::labels(8),i,u,st
 type(partition_result)::part; type(hierarchy_result)::hier; type(silhouette_result)::sil
 character(512)::out; character(:),allocatable::msg
 x=reshape([0._dp,.2_dp,-.1_dp,.1_dp,5._dp,5.2_dp,4.9_dp,5.1_dp,0._dp,-.1_dp,.2_dp,.1_dp,5._dp,4.8_dp,5.2_dp,5.1_dp],[8,2])
 labels=[1,1,1,1,2,2,2,2]; call get_command_argument(1,out); if(len_trim(out)==0)out='fortran_results.csv'
 open(newunit=u,file=trim(out),status='replace'); write(u,'(a)')'case,value,seconds,abs_tol,rel_tol'
 call cpu_time(t0); do i=1,50000; call daisy(x,d,status=st,message=msg); v=sum(d); end do; call cpu_time(t1); call emit('daisy_sum',v,t1-t0,1e-9_dp)
 call cpu_time(t0); do i=1,50000; call pam(x,2,part); v=part%objective; end do; call cpu_time(t1); call emit('pam_objective',v,t1-t0,1e-9_dp)
 call daisy(x,d,status=st,message=msg); call cpu_time(t0); do i=1,50000; call silhouette(labels,d,sil); v=sil%average_width; end do; call cpu_time(t1); call emit('silhouette_average',v,t1-t0,1e-9_dp)
 call cpu_time(t0); do i=1,50000; call agnes(x,hier,'average'); v=hier%coefficient; end do; call cpu_time(t1); call emit('agnes_coefficient',v,t1-t0,1e-8_dp)
 call cpu_time(t0); do i=1,50000; call diana(x,hier); v=hier%coefficient; end do; call cpu_time(t1); call emit('diana_coefficient',v,t1-t0,1e-8_dp)
 close(u)
contains
 subroutine emit(name,value,secs,atol); character(*),intent(in)::name; real(dp),intent(in)::value,secs,atol
  write(u,'(a,",",es24.16,",",es16.8,",",es12.4,",",es12.4)')trim(name),value,secs,atol,1e-8_dp
 end subroutine
end program
