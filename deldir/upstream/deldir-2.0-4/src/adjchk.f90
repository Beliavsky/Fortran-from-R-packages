subroutine adjchk(i,j,adj,nadj,madj,ntot)
! Check if vertices i and j are adjacent.
! Called by insrt, delet, trifnd, swap, delseg, dirseg.

dimension :: nadj(-3:ntot,0:madj)
logical :: adj
dimension :: ndi(1)

! Set dummy integer for call to intpr(...).
ndi(1) = 0

! Check if j is in the adjacency list of i.
adj = .false.
ni  = nadj(i,0)
if(ni>0) then
    do k = 1,ni
        if(j==nadj(i,k)) then
            adj = .true.
            exit
        endif
    enddo
endif

! Check if i is in the adjacency list of j.
nj = nadj(j,0)
if(nj>0) then
    do k = 1,nj
        if(i==nadj(j,k)) then
            if(adj) then
                return ! Have j in i's list and i in j's.
            else
                call intpr("Contradictory adjacency lists.",-1,ndi,0)
                call rexit("Bailing out of adjchk.")
            endif
        endif
    enddo
endif

! If we get to here i is not in j's list.
if(adj) then ! If adj is true, then j IS in i's list.
    call intpr("Contradictory adjacency lists.",-1,ndi,0)
    call rexit("Bailing out of adjchk.")
endif

end subroutine adjchk
