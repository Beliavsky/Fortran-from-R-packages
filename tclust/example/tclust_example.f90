program tclust_example
    use tclust_mod, only : dp, simulation_result, tclust_result, simulate_tclust, tclust
    implicit none
    type(simulation_result) :: sample
    type(tclust_result) :: fit
    integer :: j

    call simulate_tclust(120, sample, p=2, k=3, type_id=1, balanced=1, seed=2468)
    call tclust(sample%x, 3, fit, alpha=0.10_dp, restr_fact=12.0_dp, &
                nstart=40, niter1=3, niter2=20, nkeep=5, seed=1357)

    write (*, '(a,f14.6)') 'Trimmed classification log-likelihood: ', fit%obj
    write (*, '(a,i0)') 'Trimmed observations: ', count(fit%cluster == 0)
    do j = 1, fit%k
        write (*, '(a,i0,a,f8.4,a,*(1x,f10.4))') 'Cluster ', j, ', weight=', fit%weights(j), &
                                                ', center=', fit%centers(:, j)
    end do
end program tclust_example
