module deldir_kernel
use deldir_kinds, only: dp
use iso_fortran_env, only: error_unit
implicit none
private
public :: master, binsrt, intri, mnnd, stoke
contains

subroutine intpr(message, nchar, values, n)
character(len=*), intent(in) :: message
integer, intent(in) :: nchar, n
integer, intent(in) :: values(:)
integer :: i
if (nchar < -huge(1)) write(error_unit,*) nchar
write(error_unit,'(a)') trim(message)
if (n > 0) write(error_unit,'(*(i0,1x))') (values(i), i=1,min(n,size(values)))
end subroutine intpr

subroutine dblepr(message, nchar, values, n)
character(len=*), intent(in) :: message
integer, intent(in) :: nchar, n
real(dp), intent(in) :: values(:)
integer :: i
if (nchar < -huge(1)) write(error_unit,*) nchar
write(error_unit,'(a)') trim(message)
if (n > 0) write(error_unit,'(*(es24.16,1x))') (values(i), i=1,min(n,size(values)))
end subroutine dblepr

subroutine rexit(message)
character(len=*), intent(in) :: message
error stop message
end subroutine rexit

! ---- upstream acchk.f90 ----
subroutine acchk(i,j,k,anticl,x,y,ntot,eps)
! Check whether vertices i, j, k, are in anti-clockwise order.
! Called by locn, qtest, qtest1.
implicit none
integer :: i, i1, ijk, j, j1, k, k1, ntot
real(dp) :: cprd, eps
real(dp) :: x(-3:ntot), y(-3:ntot), xt(3), yt(3)
logical :: anticl
if(i<=0) then
    i1 = 1
else
    i1 = 0
endif
if(j<=0) then
    j1 = 1
else
    j1 = 0
endif
if(k<=0) then
    k1 = 1
else
    k1 = 0
endif
ijk = i1*4+j1*2+k1

! Get the coordinates of vertices i, j, and k. (Pseudo-coordinates for
! any ideal points.)
xt(1) = x(i)
yt(1) = y(i)
xt(2) = x(j)
yt(2) = y(j)
xt(3) = x(k)
yt(3) = y(k)

! Get the ``normalized'' cross product.
call cross(xt,yt,ijk,cprd)

! If cprd is positive then (ij-cross-ik) is directed ***upwards*** 
! and so i, j, k, are in anti-clockwise order; else not.
if(cprd > eps) then
    anticl = .true.
else
    anticl = .false.
endif
end subroutine acchk

! ---- upstream addpt.f90 ----
subroutine addpt(j,nadj,madj,x,y,ntot,eps,ntri,incAdj)
! Add point j to the triangulation.
! Called by master, dirseg.
implicit none
integer :: incadj, j, madj, n, ngap, now, ntot, ntri, nxt
real(dp) :: eps
integer :: nadj(-3:ntot,0:madj)
real(dp) :: x(-3:ntot), y(-3:ntot)
logical :: didswp
! enclosing triangle.
call initad(j,nadj,madj,x,y,ntot,eps,ntri,incAdj)

! Look at each `gap', i.e. pair of adjacent segments
! emanating from the new point; they form two sides of a
! quadrilateral; see whether the extant diagonal of this
! quadrilateral should be swapped with its alternative
! (according to the LOP: local optimality principle).
now = nadj(j,1)
nxt = nadj(j,2)
ngap = 0

do
    call swap(j,now,nxt,didswp,nadj,madj,x,y,ntot,eps,incAdj)
    n = nadj(j,0)
    if(.not.didswp) then  ! If no swap of diagonals
        now  = nxt    ! move to the next gap.
        ngap = ngap+1
    endif
    call succ(nxt,j,now,nadj,madj,ntot)
    if(ngap==n) exit
enddo
end subroutine addpt

! ---- upstream adjchk.f90 ----
subroutine adjchk(i,j,adj,nadj,madj,ntot)
implicit none
integer :: i, j, k, madj, ni, nj, ntot
! Check if vertices i and j are adjacent.
! Called by insrt, delet, trifnd, swap, delseg, dirseg.
integer :: nadj(-3:ntot,0:madj)
logical :: adj
integer :: ndi(1)

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

! ---- upstream binsrt.f90 ----
subroutine binsrt(x,y,rw,nn,ind,rind,tx,ty,ilst)
! Sort the data points into bins.
! Called by master.
implicit none
integer :: i, ink, ix, jy, k, kc, kdiv, kx, ky, nn
real(dp) :: dh, dw, h, w, xkdiv, xmax, xmin, xt, ymax, ymin, yt
real(dp) :: x(nn), y(nn), tx(nn), ty(nn)
integer :: rind(nn)
integer :: ind(nn), ilst(nn)
real(dp) :: rw(4)
integer :: ndi(1)

! Set dummy integer for call to intpr(...).
ndi(1) = 0

kdiv   = int(1+dble(nn)**0.25) ! Round high.
xkdiv  = dble(kdiv)

! Dig out the corners of the rectangular window.
xmin = rw(1)
xmax = rw(2)
ymin = rw(3)
ymax = rw(4)

w = xmax-xmin
h = ymax-ymin

! Number of bins is to be approx. sqrt(nn); thus number of subdivisions
! on each side of rectangle is approx. nn**(1/4).
dw  = w/xkdiv
dh  = h/xkdiv

! The width of each bin is dw; the height is dh.  We shall move across
! the rectangle from left to right, then up, then back from right to
! left, then up, ....  Note that kx counts the divisions from the left,
! ky counts the divisions from the bottom; kx is incremented by ink, which
! is +/- 1 and switches sign when ky is incremented; ky is always
! incremented by 1.
kx   = 1
ky   = 1
ink  = 1
k    = 0
do i = 1,nn     ! Keeps a list of those points already added
    ilst(i) = 0 ! to the new list.
enddo
do
    do i = 1,nn
        if(ilst(i) .ne. 1) then
! If the i-th point is in the current bin, add it to the list.
            xt = x(i)
            yt = y(i)
            ix = int(1+(xt-xmin)/dw)
            if(ix>kdiv) ix = kdiv
            jy = int(1+(yt-ymin)/dh)
            if(jy>kdiv) jy = kdiv
            if(ix==kx .and. jy==ky) then
                k       = k+1
                ind(i)  = k   ! Index i is the pos'n. of (x,y) in the
                rind(k) = i   ! old list; k is its pos'n. in the new one.
                tx(k)   = xt
                ty(k)   = yt
                ilst(i) = 1   ! Cross the i-th point off the old list.
            endif
        endif
    enddo
! Move to the next bin.
    kc = kx+ink
    if((1<=kc) .and. (kc<=kdiv)) then
        kx = kc
    else
        ky  = ky+1
        ink = -ink
    endif
    if(ky > kdiv) exit
enddo
! Check that all points from old list have been added to the new,
! with no spurious additions.
if(k .ne. nn) then
    call intpr("Mismatch between number of points",-1,ndi,0)
    call intpr("and number of sorted points.",-1,ndi,0)
    call rexit("Bailing out of binsrt.")
endif

! Copy the new sorted vectors back on top of the old ones.
do i = 1,nn
    x(i)   = tx(i)
    y(i)   = ty(i)
enddo

end  subroutine binsrt

! ---- upstream circen.f90 ----
subroutine circen(i,j,k,x0,y0,x,y,ntot,eps,collin)
! Find the circumcentre (x0,y0) of the triangle with
! vertices (x(i),y(i)), (x(j),y(j)), (x(k),y(k)).
! Called by qtest1, dirseg, dirout.
implicit none
integer :: i, ijk, j, k, ntot
real(dp) :: a, alpha, b, c, c1, c2, cprd, crss, d, eps, x0, y0
real(dp) :: x(-3:ntot), y(-3:ntot), xt(3), yt(3)
integer :: indv(3) ! To facillitate a lucid error message.
real(dp) :: xtmp(1)
integer :: ndi(1)
logical :: collin
ndi(1) = 0

! Get the coordinates.
xt(1) = x(i)
yt(1) = y(i)
xt(2) = x(j)
yt(2) = y(j)
xt(3) = x(k)
yt(3) = y(k)

! Check for collinearity
ijk = 0
call cross(xt,yt,ijk,cprd)
if(abs(cprd) < eps) then
    collin = .true.
else
    collin = .false.
endif

! Form the vector u from i to j, and the vector v from i to k,
! and normalize them.
a  = x(j) - x(i)
b  = y(j) - y(i)
c  = x(k) - x(i)
d  = y(k) - y(i)
c1 = sqrt(a*a+b*b)
c2 = sqrt(c*c+d*d)
a  = a/c1
b  = b/c1
c  = c/c2
d  = d/c2

! If the points are collinear, make sure that they're in the right
! order --- i between j and k.
if(collin) then
    alpha = a*c+b*d
! If they're not in the right order, bring things to
! a shuddering halt.
    if(alpha>0) then
        indv(1) = i
        indv(2) = j
        indv(3) = k
        call intpr("Point numbers:",-1,indv,3)
        xtmp(1) = alpha
        call dblepr("Test value:",-1,xtmp,1)
        call intpr("Points are collinear but in the wrong order.",-1,ndi,0)
        call rexit("Bailing out of circen.")
    endif
! Collinear, but in the right order; think of this as meaning
! that the circumcircle in question has infinite radius.
return
endif

! Not collinear; go ahead, make my circumcentre.  (First, form
! the cross product of the ***unit*** vectors, instead of the
! ``normalized'' cross product produced by ``cross''.)
crss = a*d - b*c
x0   = x(i) + 0.5*(c1*d - c2*b)/crss
y0   = y(i) + 0.5*(c2*a - c1*c)/crss

end subroutine circen

! ---- upstream cross.f90 ----
subroutine cross(x,y,ijk,cprd)
implicit none
integer :: i, ijk, ip
real(dp) :: a, b, c, cn, cprd, d, four, one, s, smin, two, zero
real(dp) :: x(3), y(3)
! Calculates a ``normalized'' cross product of the vectors joining
! [x(1),y(1)] to [x(2),y(2)] and to [x(3),y(3)] respectively.
! The normalization consists in dividing by the square of the
! shortest of the three sides of the triangle.  This normalization is
! for the purposes of testing for collinearity; if the result is less
! than epsilon, then the smallest of the sines of the angles is less than
! epsilon.

! Set constants
zero = 0.e0_dp
one  = 1.e0_dp
two  = 2.e0_dp
four = 4.e0_dp

! Adjust the coordinates depending upon which points are ideal,
! and calculate the squared length of the shortest side.

! case 0: No ideal points; no adjustment necessary.
if(ijk==0) then
    smin = -one
    do i = 1,3
        ip = i+1
        if(ip==4) ip = 1
        a = x(ip) - x(i)
        b = y(ip) - y(i)
        s = a*a+b*b
        if(smin < zero .or. s < smin) smin = s
    enddo
endif

! case 1: Only k ideal.
if(ijk==1) then
    x(2) = x(2) - x(1)
    y(2) = y(2) - y(1)
    x(1) = zero
    y(1) = zero
    cn   = sqrt(x(2)**2+y(2)**2)
    x(2) = x(2)/cn
    y(2) = y(2)/cn
    smin  = one
endif

! case 2: Only j ideal.
if(ijk==2) then
    x(3) = x(3) - x(1)
    y(3) = y(3) - y(1)
    x(1) = zero
    y(1) = zero
    cn   = sqrt(x(3)**2+y(3)**2)
    x(3) = x(3)/cn
    y(3) = y(3)/cn
    smin = one
endif

! case 3: Both j and k ideal (i not).
if(ijk==3) then
        x(1) = zero
        y(1) = zero
        smin = two
endif

! case 4: Only i ideal.
if(ijk==4) then
    x(3) = x(3) - x(2)
    y(3) = y(3) - y(2)
    x(2) = zero
    y(2) = zero
    cn   = sqrt(x(3)**2+y(3)**2)
    x(3) = x(3)/cn
    y(3) = y(3)/cn
    smin = one
endif

! case 5: Both i and k ideal (j not).
if(ijk==5) then
    x(2) = zero
    y(2) = zero
    smin = two
endif

! case 6: Both i and j ideal (k not).
if(ijk==6) then
    x(3) = zero
    y(3) = zero
    smin = two
endif

! case 7: All three points ideal; no adjustment necessary.
if(ijk==7) then
    smin = four
endif

a = x(2)-x(1)
b = y(2)-y(1)
c = x(3)-x(1)
d = y(3)-y(1)

cprd = (a*d - b*c)/smin
end subroutine cross

! ---- upstream delet.f90 ----
subroutine delet(i,j,nadj,madj,ntot)
! Delete i and j from each other's adjacency lists.
! Called by initad, swap.
implicit none
integer :: i, j, madj, ntot
integer :: nadj(-3:ntot,0:madj)
logical :: adj
call adjchk(i,j,adj,nadj,madj,ntot)

! Then do the actual deletion if they are.
if(adj) then
    call delet1(i,j,nadj,madj,ntot)
    call delet1(j,i,nadj,madj,ntot)
endif

return
end

! ---- upstream delet1.f90 ----
subroutine delet1(i,j,nadj,madj,ntot)
! Delete j from the adjacency list of i.
! Called by delet.
implicit none
integer :: i, j, k, kk, madj, n, ntot
integer :: nadj(-3:ntot,0:madj)

n    = nadj(i,0)
do k = 1,n
    if(nadj(i,k)==j) then ! Find j in the list;
                          ! then move everything back one notch.
        do kk = k,n-1
             nadj(i,kk) = nadj(i,kk+1)
        enddo
        nadj(i,n) = -99   ! Changed from the confusing 0 value 25/7/2011.
        nadj(i,0) = n-1
        return
    endif
enddo

end subroutine delet1

! ---- upstream delout.f90 ----
subroutine delout(delsum,nadj,madj,x,y,ntot,nn)

! Put a summary of the Delaunay triangles with a vertex at point i,
! for i = 1, ..., nn, into the array delsum.  Do this in the original
! order of the points, not the order into which they have been
! bin-sorted.
! Called by master.
implicit none
integer :: i, j, j1, k, kp, madj, nn, np, npt, ntot
real(dp) :: area, tmp, xi, xj, xk, yi, yj, yk
integer :: nadj(-3:ntot,0:madj)
real(dp) :: x(-3:ntot), y(-3:ntot)
real(dp) :: delsum(nn,4)

do i = 1,nn
    area = 0.e0_dp ! Initialize area of polygon consisting of triangles
                ! with a vertex at point i.
! Get the coordinates of the point and the number of
! (real) triangles emanating from it.
    np = nadj(i,0)
    xi = x(i)
    yi = y(i)
    npt = np
    do k = 1,np
        kp = k+1
        if(kp>np) kp = 1
        if(nadj(i,k)<=0 .or. nadj(i,kp)<=0) npt = npt-1
    enddo

! For each point in the adjacency list of point i, find its
! successor, and the area of the triangle determined by these
! three points.
    do j1 = 1,np
        j = nadj(i,j1)
        if(j<=0) cycle
        xj = x(j)
        yj = y(j)
        call succ(k,i,j,nadj,madj,ntot)
        if(k<=0) cycle
        xk = x(k)
        yk = y(k)
        call triar(xi,yi,xj,yj,xk,yk,tmp)
! Downweight the area by 1/3, since each
! triangle eventually appears 3 times over.
        area = area+tmp/3.e0_dp
    enddo
    delsum(i,1) = xi
    delsum(i,2) = yi
    delsum(i,3) = npt
    delsum(i,4) = area
enddo

end subroutine delout

! ---- upstream delseg.f90 ----
subroutine delseg(delsgs,ndel,nadj,madj,nn,x,y,ntot,incSeg)

! Output the endpoints of the line segments joining the
! vertices of the Delaunay triangles.
! Called by master.
implicit none
integer :: i, incseg, j, kseg, madj, ndel, nn, ntot
logical :: value
integer :: nadj(-3:ntot,0:madj)
real(dp) :: x(-3:ntot), y(-3:ntot)
real(dp) :: delsgs(6,ndel)

! Initialise incSeg
incSeg = 0

! For each distinct pair of points i and j, if they are adjacent
! then put their endpoints into the output array.
nn = ntot-4
kseg = 0
do i = 2,nn
    do j = 1,i-1
        call adjchk(i,j,value,nadj,madj,ntot)
        if(value) then
            kseg = kseg+1
            if(kseg > ndel) then
                incSeg = 1
                return
            endif
            delsgs(1,kseg) = x(i)
            delsgs(2,kseg) = y(i)
            delsgs(3,kseg) = x(j)
            delsgs(4,kseg) = y(j)
            delsgs(5,kseg) = i
            delsgs(6,kseg) = j
        endif
    enddo
enddo
ndel = kseg

end subroutine delseg

! ---- upstream dirout.f90 ----
subroutine dirout(dirsum,nadj,madj,x,y,ntot,nn,rw,eps)

! Output the description of the Dirichlet tile centred at point
! i for i = 1, ..., nn.  Do this in the original order of the
! points, not in the order into which they have been bin-sorted.
! Called by master.
implicit none
integer :: i, j, j1, k, l, madj, nbpt, nedge, nn, np, npt, ntot
real(dp) :: a, ai, area, b, bi, c, ci, d, di, eps, slope, sn, tmp, xi, xj, xm, xmax, xmin, yi, yj, ym, ymax, &
    & ymin
integer :: nadj(-3:ntot,0:madj)
real(dp) :: x(-3:ntot), y(-3:ntot)
real(dp) :: dirsum(nn,3), rw(4)
integer :: ndi(1)
logical :: collin, intfnd, bptab, bptcd, rwu
ndi(1) = 0

! Note that at this point some Delaunay neighbours may be
! `spurious'; they are the corners of a `large' rectangle in which
! the rectangular window of interest has been suspended.  This
! large window was brought in simply to facilitate output concerning
! the Dirichlet tesselation.  They were added to the triangulation
! in the routine `dirseg' which ***must*** therefore be called before
! this routine (`dirout') is called.  (Likewise `dirseg' must be called
! ***after*** `delseg' and `delout' have been called.)

! Dig out the corners of the rectangular window.
xmin = rw(1)
xmax = rw(2)
ymin = rw(3)
ymax = rw(4)

do i = 1,nn
    area = 0.e0_dp ! Initialize the area of the ith tile to zero.
    nbpt = 0    ! Initialize the number of boundary points of
                ! the ith tile to zero.
    npt  = 0    ! Initialize the number of tile boundaries to zero.
    np   = nadj(i,0)

! Output the point number, its coordinates, and the number of
! its Delaunay neighbours == the number of boundary segments in
! its Dirichlet tile.
! For each Delaunay neighbour, find the circumcentres of the
! triangles on each side of the segment joining point i to that
! neighbour.
    !call dblepr("rw:",-1,rw,4)
    do j1 = 1,np
        j = nadj(i,j1)
        call pred(k,i,j,nadj,madj,ntot)
        call succ(l,i,j,nadj,madj,ntot)
        call circen(i,k,j,a,b,x,y,ntot,eps,collin)
        if(collin) then
            call intpr("Vertices of triangle are collinear.",-1,ndi,0)
            call rexit("Bailing out of dirout.")
        endif
        call circen(i,j,l,c,d,x,y,ntot,eps,collin)
        if(collin) then
            call intpr("Vertices of triangle are collinear.",-1,ndi,0)
            call rexit("Bailing out of dirout.")
        endif

! Increment the area of the current Dirichlet
! tile (intersected with the rectangular window) by applying
! Stokes' Theorem to the segment of tile boundary joining
! (a,b) to (c,d).  (Note that the direction is anti-clockwise.)
        call stoke(a,b,c,d,rw,tmp,sn,eps)
        area = area+sn*tmp

! If a circumcentre is outside the rectangular window, replace
! it with the intersection of the rectangle boundary with the
! line joining the two circumcentres.  Then output
! the number of the current Delaunay neighbour and
! the two circumcentres (or the points with which
! they have been replaced).
! Note: rwu = "right way up".
        xi = x(i)
        xj = x(j)
        yi = y(i)
        yj = y(j)
        if(abs(yi-yj) > eps) then
            slope = (xi - xj)/(yj - yi)
            rwu = .true.
        else
            slope = 0.e0_dp
            rwu = .false.
        endif
        call dldins(a,b,slope,rwu,ai,bi,rw,intfnd,bptab,nedge)
        if(intfnd) then
            call dldins(c,d,slope,rwu,ci,di,rw,intfnd,bptcd,nedge)
            if(.not.intfnd) then
                call intpr("Line from midpoint to circumcenter",-1,ndi,0)
                call intpr("does not intersect rectangle boundary!",-1,ndi,0)
                call intpr("But it HAS to!!!",-1,ndi,0)
                call rexit("Bailing out of dirout.")
            endif
            if(bptab .and. bptcd) then
                xm = 0.5*(ai+ci)
                ym = 0.5*(bi+di)
                if(xmin<xm .and. xm<xmax .and. ymin<ym .and. ym<ymax) then
                    nbpt = nbpt+2
                    npt  = npt+1
                endif
            else
                npt = npt + 1
                if(bptab .or. bptcd) nbpt = nbpt+1
            endif
        endif
    enddo
    dirsum(i,1) = npt
    dirsum(i,2) = nbpt
    dirsum(i,3) = area
enddo

end subroutine dirout

! ---- upstream dirseg.f90 ----
subroutine dirseg(dirsgs,ndir,nadj,madj,nn,x,y,ntot,rw,eps,ntri,incAdj,incSeg)

! Output the endpoints of the segments of boundaries of Dirichlet
! tiles.  (Do it economically; each such segment once and only once.)
! Called by master.
implicit none
integer :: i, incadj, incseg, j, k, kseg, l, madj, ndir, nedgeab, nedgecd, nn, nstt, ntot, ntri
real(dp) :: a, ai, b, bi, c, ci, d, di, eps, slope, xi, xj, xm, xmax, xmin, yi, yj, ym, ymax, ymin
logical :: collin, adjace, intfnd, bptab, bptcd, goferit, rwu
integer :: nadj(-3:ntot,0:madj)
real(dp) :: x(-3:ntot), y(-3:ntot)
real(dp) :: dirsgs(10,ndir), rw(4)
integer :: ndi(1)

! Set dummy integer for call to intpr(...).
ndi(1) = 0

! Initialise incSeg
incSeg = 0

! Add in some dummy corner points, outside the actual window.
! Far enough out so that no resulting tile boundaries intersect the
! window.

! Note that these dummy corners are needed by the routine `dirout'
! but will screw things up for `delseg' and `delout'.  Therefore
! this routine (`dirseg') must be called ***before*** dirout, and
! ***after*** delseg and delout.

! Dig out the corners of the rectangular window.
xmin = rw(1)
xmax = rw(2)
ymin = rw(3)
ymax = rw(4)

a = xmax-xmin
b = ymax-ymin
c = sqrt(a*a+b*b)

nn  = ntot-4
nstt = nn+1
i = nstt
x(i) = xmin-c
y(i) = ymin-c
i = i+1
x(i) = xmax+c
y(i) = ymin-c
i = i+1
x(i) = xmax+c
y(i) = ymax+c
i = i+1
x(i) = xmin-c
y(i) = ymax+c

do j = nstt,ntot
    call addpt(j,nadj,madj,x,y,ntot,eps,ntri,incAdj)
    if(incAdj==1) return
    ntri = ntri + 3
enddo

! Put the segments into the array dirsgs.

! For each distinct pair of (genuine) data points, find out if they are
! adjacent.  If so, find the circumcentres of the triangles lying on each
! side of the segment joining them.
kseg = 0
do i = 2,nn
    do j = 1,i-1
        call adjchk(i,j,adjace,nadj,madj,ntot)
        if(adjace) then
            call pred(k,i,j,nadj,madj,ntot)
            call circen(i,k,j,a,b,x,y,ntot,eps,collin)
            if(collin) then
                call intpr("Vertices of triangle are collinear.",-1,ndi,0)
                call rexit("Bailing out of dirseg.")
            endif
            call succ(l,i,j,nadj,madj,ntot)
            call circen(i,j,l,c,d,x,y,ntot,eps,collin)
            if(collin) then
                call intpr("Vertices of triangle are collinear.",-1,ndi,0)
                call rexit("Bailing out of dirseg.")
            endif
! If a circumcentre is outside the rectangular window
! of interest, draw a line joining it to the other
! circumcentre.  Find the intersection of this line with
! the boundary of the window; for (a,b) and call the point
! of intersection (ai,bi).  For (c,d), call it (ci,di).
! Note: rwu = "right way up".
            xi = x(i)
            xj = x(j)
            yi = y(i)
            yj = y(j)
            if(abs(yi-yj) > eps) then
                slope = (xi - xj)/(yj - yi)
                rwu   = .true.
            else
                slope = 0.e0_dp
                rwu   = .false.
            endif
            call dldins(a,b,slope,rwu,ai,bi,rw,intfnd,bptab,nedgeab)
            if(.not.intfnd) then
                call intpr("Line from midpoint to circumcenter",-1,ndi,0)
                call intpr("does not intersect rectangle boundary!",-1,ndi,0)
                call intpr("But it HAS to!!!",-1,ndi,0)
                call rexit("Bailing out of dirseg.")
            endif
            call dldins(c,d,slope,rwu,ci,di,rw,intfnd,bptcd,nedgecd)
            if(.not.intfnd) then
                call intpr("Line from midpoint to circumcenter",-1,ndi,0)
                call intpr("does not intersect rectangle boundary!",-1,ndi,0)
                call intpr("But it HAS to!!!",-1,ndi,0)
                call rexit("Bailing out of dirseg.")
            endif
            goferit = .false.
            if(bptab .and. bptcd) then
                xm = 0.5*(ai+ci)
                ym = 0.5*(bi+di)
                if(xmin<xm.and.xm<xmax.and.ymin<ym.and.ym<ymax) then
                    goferit = .true.
                endif
            endif
            if((.not.bptab) .or. (.not.bptcd)) goferit = .true.
            if(goferit) then
                kseg = kseg + 1
                if(kseg > ndir) then
                    incSeg = 1
                    return
                endif
                dirsgs(1,kseg) = ai
                dirsgs(2,kseg) = bi
                dirsgs(3,kseg) = ci
                dirsgs(4,kseg) = di
                dirsgs(5,kseg) = i
                dirsgs(6,kseg) = j
                if(bptab) then
                     dirsgs(7,kseg) = 1.e0_dp
                else
                     dirsgs(7,kseg) = 0.e0_dp
                endif
                if(bptcd) then
                    dirsgs(8,kseg) = 1.e0_dp
                else
                    dirsgs(8,kseg) = 0.e0_dp
                endif
                if(bptab) then
                    dirsgs(9,kseg) = -nedgeab
                else
                    dirsgs(9,kseg) = k
                endif
                if(bptcd) then
                    dirsgs(10,kseg) = -nedgecd
                else
                    dirsgs(10,kseg) = l
                endif
            endif
        endif
    enddo
enddo
ndir = kseg

end subroutine dirseg

! ---- upstream dldins.f90 ----
subroutine dldins(a,b,slope,rwu,ai,bi,rw,intfnd,bpt,nedge)

! Get a point ***inside*** the rectangular window on the ray from
! one circumcentre to the next one.  I.e. if the `next one' is
! inside, then that's it; else find the intersection of this ray with
! the boundary of the rectangle.
! Called by dirseg, dirout.
implicit none
integer :: nedge
real(dp) :: a, ai, b, bi, slope, xmax, xmin, ymax, ymin
real(dp) :: rw(4)
logical :: intfnd, bpt, rwu
! and that slope is the slope of the ray joining (a,b) to the
! corresponding circumcentre on the opposite side of an edge of that
! triangle.  When `dldins' is called by `dirout' it is possible
! for the ray not to intersect the window at all.  (The Delaunay
! edge between the two circumcentres might be connected to a `fake
! outer corner', added to facilitate constructing a tiling that
! completely covers the actual window.)  The variable `intfnd' acts
! as an indicator as to whether such an intersection has been found.

! The variable `bpt' acts as an indicator as to whether the returned
! point (ai,bi) is a true circumcentre, inside the window (bpt == .false.),
! or is the intersection of a ray with the boundary of the window
! (bpt = .true.).


intfnd = .true.
bpt = .true.

! Dig out the corners of the rectangular window.
xmin = rw(1)
xmax = rw(2)
ymin = rw(3)
ymax = rw(4)

! Check if (a,b) is inside the rectangle.
if(xmin<=a .and. a<=xmax .and. ymin<=b .and. b<=ymax) then
    ai = a
    bi = b
    bpt = .false.
    nedge = 0
    return
endif

! Look for appropriate intersections with the four lines forming
! the sides of the rectangular window.

! If not "the right way up" then the line joining the two
! circumcentres is vertical.

if(.not.rwu) then
    if(b < ymin) then
        ai = a
        bi = ymin
        nedge = 1
        if(xmin<=ai .and. ai<=xmax) return
    endif
    if(b > ymax) then
        ai = a
        bi = ymax
        nedge = 3
        if(xmin<=ai .and. ai<=xmax) return
    endif
    intfnd = .false.
    return
endif

! Line 1: x = xmin.
if(a<xmin) then
    ai = xmin
    bi = b + slope*(ai-a)
    nedge = 2
    if(ymin<=bi .and. bi<=ymax) return
endif

! Line 2: y = ymin.
if(b<ymin) then
    bi = ymin
    ai = a + (bi-b)/slope
    nedge = 1
    if(xmin<=ai .and. ai<=xmax) return
endif

! Line 3: x = xmax.
if(a>xmax) then
    ai = xmax
    bi = b + slope*(ai-a)
    nedge = 4
    if(ymin<=bi .and. bi<=ymax) return
endif

! Line 4: y = ymax.
if(b>ymax) then
    bi = ymax
    ai = a + (bi-b)/slope
    nedge = 3
    if(xmin<=ai .and. ai<=xmax) return
endif

intfnd = .false.
end subroutine dldins

! ---- upstream initad.f90 ----
subroutine initad(j,nadj,madj,x,y,ntot,eps,ntri,incAdj)

! Initial adding-in of a new point j.
! Called by addpt.
implicit none
integer :: i, incadj, ip, j, k, kk, madj, nedge, ntot, ntri
real(dp) :: eps
integer :: nadj(-3:ntot,0:madj)
real(dp) :: x(-3:ntot), y(-3:ntot)
integer :: tau(3)
call trifnd(j,tau,nedge,nadj,madj,x,y,ntot,eps,ntri)

! If the new point is on the edge of a triangle, detach the two
! vertices of that edge from each other.  Also join j to the vertex
! of the triangle on the reverse side of that edge from the `found'
! triangle (defined by tau) -- given that there ***is*** such a triangle.
if(nedge .ne. 0) then
    ip = nedge
    i  = ip-1
    if(i==0) i = 3 ! Arithmetic modulo 3.
    call pred(k,tau(i),tau(ip),nadj,madj,ntot)
    call succ(kk,tau(ip),tau(i),nadj,madj,ntot)
    call delet(tau(i),tau(ip),nadj,madj,ntot)
    if(k==kk) call insrt(j,k,nadj,madj,x,y,ntot,eps,incAdj)
    if(incAdj==1) return
endif

! Join the new point to each of the three vertices.
do i = 1,3
    call insrt(j,tau(i),nadj,madj,x,y,ntot,eps,incAdj)
enddo
end subroutine initad

! ---- upstream insrt.f90 ----
subroutine insrt(i,j,nadj,madj,x,y,ntot,eps,incAdj)
! Insert i and j into each other's adjacency list.
! Called by master, initad, swap.
implicit none
integer :: i, incadj, j, ki, kj, madj, ntot
real(dp) :: eps
integer :: nadj(-3:ntot,0:madj)
real(dp) :: x(-3:ntot), y(-3:ntot)
logical :: adj
call adjchk(i,j,adj,nadj,madj,ntot)
if(adj) return

! If not, find where in each list they should respectively be.
call locn(i,j,kj,nadj,madj,x,y,ntot,eps)
call locn(j,i,ki,nadj,madj,x,y,ntot,eps)

! Put them in each other's lists in the appropriate position.
call insrt1(i,j,kj,nadj,madj,ntot,incAdj)
if(incAdj==1) return
call insrt1(j,i,ki,nadj,madj,ntot,incAdj)
if(incAdj==1) return ! This seems unnecessary; check on this.

end subroutine insrt

! ---- upstream insrt1.f90 ----
subroutine insrt1(i,j,kj,nadj,madj,ntot,incAdj)

! Insert j into the adjacency list of i.
! Called by insrt.
implicit none
integer :: i, incadj, j, kj, kk, madj, n, ntot
integer :: nadj(-3:ntot,0:madj)

! Initialise incAdj.
incAdj = 0

! Variable  kj is the index which j ***will***
! have when it is inserted into the adjacency list of i in
! the appropriate position.

! If the adjacency list of i had no points just stick j into the list.
n = nadj(i,0)
if(n==0) then
    nadj(i,0) = 1
    nadj(i,1) = j
    return
endif

! If the adjacency list had some points, move everything ahead of the
! kj-th place one place forward, and put j in position kj.
kk = n+1

if(kk>madj) then ! Watch out for over-writing!!!
    incAdj = 1
    return
endif

do
    nadj(i,kk) = nadj(i,kk-1)
    kk = kk-1
    if(kk <= kj) exit
enddo

nadj(i,kj) = j
nadj(i,0)  = n+1

end subroutine insrt1

! ---- upstream intri.f90 ----
subroutine intri(x,y,u,v,n,okay)
!
! Test whether any of the points (u(i),v(i)) are inside the triangle
! whose vertices are specified by the vectors x and y.
! Called by .Fortran() from triang.list.R.
!
implicit none
integer :: i, j, jp, n
real(dp) :: a, b, c, cp, d, s, zero
real(dp) :: x(3), y(3), u(n), v(n)
integer :: okay
logical :: inside

zero = 0.e0_dp

! Check on order (clockwise or anticlockwise).
s = 1.e0_dp
a = x(2) - x(1)
b = y(2) - y(1)
c = x(3) - x(1)
d = y(3) - y(1)
cp = a*d - b*c
if(cp < 0) s = -s
do i = 1,n
    inside = .true.
    do j = 1,3
        jp = j+1
        if(jp==4) jp = 1 ! Take addition modulo 3.
        a  = x(jp) - x(j)
        b  = y(jp) - y(j)
        c  = u(i)  - x(j)
        d  = v(i)  - y(j)
        cp = s*(a*d - b*c)
        if(cp <= zero) then
            inside = .false.
            exit
        endif
    enddo
    if(inside) then
        okay = 0
        return
    endif
enddo
okay = 1
end subroutine intri

! ---- upstream locn.f90 ----
subroutine locn(i,j,kj,nadj,madj,x,y,ntot,eps)

! Find the appropriate location for j in the adjacency list
! of i.  This is the index which j ***will*** have when
! it is inserted into the adjacency list of i in the
! appropriate place.  Called by insrt.
implicit none
integer :: i, j, k, kj, km, ks, madj, n, ntot
real(dp) :: eps
integer :: nadj(-3:ntot,0:madj)
real(dp) :: x(-3:ntot), y(-3:ntot)
logical :: before

n = nadj(i,0)

! If there is nothing already adjacent to i, then j will have place 1.
if(n==0) then
    kj = 1
    return
endif

! Run through i's list, checking if j should come before each element
! of that list.  (I.e. if i, j, and k are in anti-clockwise order.)
! If j comes before the kj-th item, but not before the (kj-1)-st, then
! j should have place kj.
do ks = 1,n
    kj = ks
    k = nadj(i,kj)
    call acchk(i,j,k,before,x,y,ntot,eps)
    if(before) then
        km = kj-1
        if(km==0) km = n
        k = nadj(i,km)
        call acchk(i,j,k,before,x,y,ntot,eps)
! Got here.
        if(.not.before) then
! If j is before 1 and after n, then it should
! have place n+1.
            if(kj==1) kj = n+1
            return
        endif
    endif
enddo

! We've gone right through the list and haven't been before
! the kj-th item ***and*** after the (kj-1)-st on any occasion.
! Therefore j is before everything (==> place 1) or after
! everything (==> place n+1).
if(before) then
    kj = 1
else
    kj = n+1
endif
end subroutine locn

! ---- upstream master.f90 ----
subroutine master(x,y,rw,nn,ntot,nadj,madj,eps,delsgs,ndel,delsum,&
                  dirsgs,ndir,dirsum,incAdj,incSeg)

! Master subroutine:
! One subroutine to rule them all,
! One subroutine to find them.
! One subroutine to bring them all in,
! And in the darkness bind them.

! Note: "incAdj" <--> increase size of adjacency list.
!       "incSeg" <--> increase size of storage for segments.
implicit none
integer :: i, incadj, incseg, j, k, madj, ndel, ndir, nn, ntot, ntri
real(dp) :: eps, one
real(dp) :: x(-3:ntot), y(-3:ntot)
integer :: nadj(-3:ntot,0:madj)
real(dp) :: rw(4)
real(dp) :: delsgs(6,ndel), dirsgs(10,ndir)
real(dp) :: delsum(nn,4), dirsum(nn,3)

! Define one.
one = 1.e0_dp

! Initialize the adjacency list; counts to 0, other entries to -99.
do i = -3,ntot
    nadj(i,0) = 0
    do j = 1,madj
        nadj(i,j) = -99
    enddo
enddo

! Put the four ideal points into x and y and the adjacency list.
! The ideal points are given pseudo-coordinates
! (-1,-1), (1,-1), (1,1), and (-1,1).  They are numbered as
!    0       -1      -2         -3
! i.e. the numbers decrease anticlockwise from the
! `bottom left corner'.
x(-3) = -one
y(-3) =  one
x(-2) =  one
y(-2) =  one
x(-1) =  one
y(-1) = -one
x(0)  = -one
y(0)  = -one

do i = 1,4
    j = i-4
    k = j+1
    if(k>0) k = -3
    call insrt(j,k,nadj,madj,x,y,ntot,eps,incAdj)
    if(incAdj==1) return
enddo

! Put in the first of the point set into the adjacency list.
do i = 1,4
    j = i-4
    call insrt(1,j,nadj,madj,x,y,ntot,eps,incAdj)
    if(incAdj==1) return
enddo
ntri = 4

! Now add the rest of the point set
do j = 2,nn
    call addpt(j,nadj,madj,x,y,ntot,eps,ntri,incAdj)
    if(incAdj==1) return
    ntri = ntri + 3
enddo

! Obtain the description of the triangulation.
call delseg(delsgs,ndel,nadj,madj,nn,x,y,ntot,incSeg)
if(incSeg==1) return

call delout(delsum,nadj,madj,x,y,ntot,nn)

call dirseg(dirsgs,ndir,nadj,madj,nn,x,y,ntot,rw,eps,ntri,incAdj,incSeg)
if(incAdj==1 .or. incSeg==1) return
call dirout(dirsum,nadj,madj,x,y,ntot,nn,rw,eps)
end subroutine master

! ---- upstream mnnd.f90 ----
subroutine mnnd(x,y,n,dminbig,dminav)
!
! Mean nearest neighbour distance.  Called by .Fortran()
! from mnnd.R.
!
implicit none
integer :: i, j, n
real(dp) :: d, dmin, dminav, dminbig
real(dp) :: x(n), y(n)

dminav = 0.e0_dp
do i = 1,n
    dmin = dminbig
    do j = 1,n
        if(i .ne. j) then
            d = (x(i)-x(j))**2 + (y(i)-y(j))**2
            if(d < dmin) dmin = d
        endif
    enddo
    dminav = dminav + sqrt(dmin)
enddo

dminav = dminav/n

return
end

! ---- upstream pred.f90 ----
subroutine pred(kpr,i,j,nadj,madj,ntot)

! Find the predecessor of j in the adjacency list of i.
! Called by initad, trifnd, swap, dirseg, dirout.
implicit none
integer :: i, j, k, km, kpr, madj, n, ntot
integer :: nadj(-3:ntot,0:madj)
integer :: ndi(1)

! Set dummy integer for call to intpr(...).
ndi(1) = 0

n = nadj(i,0)

! If the adjacency list of i is empty, then clearly j has no predecessor
! in this adjacency list. Something's wrong; stop.
if(n==0) then
    call intpr("Adjacency list of i is empty, and so cannot contain j.",-1,ndi,0)
    call rexit("Bailing out of pred.")
endif

! The adjacency list of i is non-empty; search through it until j is found;
! subtract 1 from the location of j, and find the contents of this new location
do k = 1,n
    if(j==nadj(i,k)) then
        km = k-1
        if(km<1) km = n         ! Take km modulo n. (The adjacency list
        kpr = nadj(i,km)        ! is circular.)
        return
    endif
enddo

! The adjacency list for i doesn't contain j.  Something's wrong; stop.
    call intpr("Adjacency list of i does not contain j.",-1,ndi,0)
    call rexit("Bailing out of pred.")
end subroutine pred

! ---- upstream qtest.f90 ----
subroutine qtest(h,i,j,k,shdswp,x,y,ntot,eps)

! Test whether the LOP is satisified; i.e. whether vertex j
! is outside the circumcircle of vertices h, i, and k of the
! quadrilateral.  Vertex h is the vertex being added; i and k are
! the vertices of the quadrilateral which are currently joined;
! j is the vertex which is ``opposite'' the vertex being added.
! If the LOP is not satisfied, then shdswp ("should-swap") is true,
! i.e. h and j should be joined, rather than i and k.  I.e. if j
! is outside the circumcircle of h, i, and k then all is well as-is;
! *don't* swap ik for hj.  If j is inside the circumcircle of h,
! i, and k then change is needed so swap ik for hj.
! Called by swap.
implicit none
integer :: i, ii, ijk, j, jj, k, kk, ntot
real(dp) :: eps, ss, test, xh, xi, xk, yh, yi, yk
real(dp) :: x(-3:ntot), y(-3:ntot)
integer :: ndi(1)
integer :: h
logical :: shdswp
ndi(1) = 0

! Look for ideal points.
if(i<=0) then
    ii = 1
else
    ii = 0
endif

if(j<=0) then
    jj = 1
else
    jj = 0
endif

if(k<=0) then
    kk = 1
else
    kk = 0
endif
ijk = ii*4+jj*2+kk

! All three corners other than h (the point currently being
! added) are ideal --- so h, i, and k are co-linear; so 
! i and k shouldn't be joined, and h should be joined to j.
! So swap.  (But this can't happen, anyway!!!)
! case 7:
if(ijk==7) then
    shdswp = .true.
    return
endif

! If i and j are ideal, find out which of h and k is closer to the
! intersection point of the two diagonals, and choose the diagonal
! emanating from that vertex.  (I.e. if h is closer, swap.)
! Unless swapping yields a non-convex quadrilateral!!!
! case 6:
if(ijk==6) then
    xh = x(h)
    yh = y(h)
    xk = x(k)
    yk = y(k)
    ss = 1 - 2*mod(-j,2)
    test = (xh*yk+xk*yh-xh*yh-xk*yk)*ss
    if(test>0.e0_dp) then
        shdswp = .true.
    else
        shdswp = .false.
    endif
! Check for convexity:
    if(shdswp) call acchk(j,k,h,shdswp,x,y,ntot,eps)
    return
endif

! Vertices i and k are ideal --- can't happen, but if it did, we'd
! increase the minimum angle ``from 0 to more than 2*0'' by swapping ...
!
! 24/7/2011 --- I now think that the forgoing comment is misleading,
! although it doesn't matter since it can't happen anyway.  The
! ``2*0'' is wrong.  The ``new minimum angle would be min(alpha,beta)
! where alpha and beta are the angles made by the line joining h
! to j with (any) line with slope = -1.  This will be greater than
! 0 unless the line from h to j has slope = - 1.  In this case h,
! i, j, and k are all co-linear, so i and k should not be joined
! (and h and j should be) so swapping is called for.  If h, i,
! j and j are not co-linear then the quadrilateral is definitely
! convex whence swapping is OK.  So let's say swap.
! case 5:
if(ijk==5) then
    shdswp = .true.
    return
endif

! If i is ideal we'd increase the minimum angle ``from 0 to more than
! 2*0'' by swapping, so just check for convexity:
! case 4:
if(ijk==4) then
    call acchk(j,k,h,shdswp,x,y,ntot,eps)
    return
endif

! If j and k are ideal, this is like unto case 6.
! case 3:
if(ijk==3) then
    xi = x(i)
    yi = y(i)
    xh = x(h)
    yh = y(h)
    ss = 1 - 2*mod(-j,2)
    test = (xh*yi+xi*yh-xh*yh-xi*yi)*ss
    if(test>0.e0_dp) then
        shdswp = .true.
    else
        shdswp = .false.
    endif
! Check for convexity:
    if(shdswp) call acchk(h,i,j,shdswp,x,y,ntot,eps)
    return
endif

! If j is ideal we'd decrease the minimum angle ``from more than 2*0
! to 0'', by swapping; so don't swap.
! case 2:
if(ijk==2) then
    shdswp = .false.
    return
endif

! If k is ideal, this is like unto case 4.
! case 1:
if(ijk==1) then
    call acchk(h,i,j,shdswp,x,y,ntot,eps) ! This checks
                                          ! for convexity.
                                          ! (Was i,j,h,...)
    return
endif

! If none of the `other' three corners are ideal, do the Lee-Schacter
! test for the LOP.
! case 0:
if(ijk==0) then
    call qtest1(h,i,j,k,x,y,ntot,eps,shdswp)
    return
endif

! default:  ! This CAN'T happen!
call intpr("Indicator ijk is out of range.",-1,ndi,0)
call intpr("This CAN'T happen!",-1,ndi,0)
call rexit("Bailing out of qtest.")
end subroutine qtest

! ---- upstream qtest1.f90 ----
subroutine qtest1(h,i,j,k,x,y,ntot,eps,shdswp)

! The Lee-Schacter test for the LOP (all points are real,
! i.e. non-ideal).  If the LOP is ***not*** satisfied (i.e. if
! vertex j is inside the circumcircle of vertices h, i, and k) then the
! diagonals should be swapped, i.e. shdswp ("should-swap") is true.
! Called by qtest.
implicit none
integer :: i, j, k, nid, ntot
real(dp) :: a, alpha, b, c, c1, c2, ch, cprd, d, eps, r2, x0, xh, xj, y0, yh, yj
integer :: indv(3)
real(dp) :: x(-3:ntot), y(-3:ntot), xt(3), yt(3)
integer :: itmp(1)
real(dp) :: xtmp(1)
integer :: ndi(1)
integer :: h
logical :: shdswp, collin
ndi(1) = 0

! The vertices of the quadrilateral are labelled
! h, i, j, k in the anticlockwise direction, h
! being the point of central interest.

! Make sure the quadrilateral is convex, so that
! it makes sense to swap the diagonal.
! call acchk(i,j,k,shdswp,x,y,ntot,eps)
! if(!shdswp) return
!
! 23 July 2011:
! The foregoing test is a load of dingoes' kidneys.  (1) It is
! unnecessary, and (2) it is wrong!  (1) If the LOP is not satisfied
! (the only circumstance under which there should be a swap) then the
! quadrilateral ***must*** be convex, and so swapping can sensibly
! take place.  (2) The vertices i, j, k in will ***always*** be in
! anticlockwise order, since the vertices h, i, j, k of the quadrilateral
! are in such order and i is connected to k, whence j can't be inside
! the triangle ihk.  So the test does nothing.  But then it didn't need
! to do anything.

! Check for collinearity of points h, i and k.
xt(1) = x(h)
yt(1) = y(h)
xt(2) = x(i)
yt(2) = y(i)
xt(3) = x(k)
yt(3) = y(k)
nid = 0  ! nid = number of ideal points.
call cross(xt,yt,nid,cprd)
collin = (abs(cprd) < eps) ! Does this work???

! If the points are collinear, make sure that they're in the right
! order --- h between i and k.
if(collin) then
! Form the vector u from h to i, and the vector v from h to k,
! and normalize them.
    a  = xt(2) - xt(1)
    b  = yt(2) - yt(1)
    c  = xt(3) - xt(1)
    d  = yt(3) - yt(1)
    c1 = sqrt(a*a+b*b)
    c2 = sqrt(c*c+d*d)
    a  = a/c1
    b  = b/c1
    c  = c/c2
    d  = d/c2
    alpha = a*c+b*d
! If they're not in the right order, bring things to
! a shuddering halt.
    if(alpha>0) then
        itmp(1) = 1
        indv(1) = i
        indv(2) = j
        indv(3) = k
        itmp(1) = h
        call intpr("Point being added, h:",-1,itmp,1)
        call intpr("now, other vertex, nxt:",-1,indv,3)
        xtmp(1) = alpha
        call dblepr("Test value:",-1,xtmp,1)
        call intpr("Points are collinear but h is not between i and k.",-1,ndi,0)
        call rexit("Bailing out of qtest1.")
    endif
! Collinear, and in the right order; think of this as meaning
! that the circumcircle in question has infinite radius.
    shdswp = .true.
endif

! Get the coordinates of vertices h and j.
xh = x(h)
yh = y(h)
xj = x(j)
yj = y(j)

! Find the centre of the circumcircle of vertices h, i, k.
call circen(h,i,k,x0,y0,x,y,ntot,eps,shdswp)
if(shdswp) return ! The points h, i, and k are colinear, so
                  ! the circumcircle has `infinite radius', so
                  ! (xj,yj) is definitely inside!

! Check whether (xj,yj) is inside the circle of centre
! (x0,y0) and radius r = dist[(x0,y0),(xh,yh)]

a  = x0-xh
b  = y0-yh
r2 = a*a+b*b
a  = x0-xj
b  = y0-yj
ch = a*a + b*b
if(ch<r2) then
    shdswp = .true.
else
    shdswp = .false.
endif
end subroutine qtest1

! ---- upstream stoke.f90 ----
subroutine stoke(x1,y1,x2,y2,rw,area,s1,eps)

! Apply Stokes' theorem to find the area of a polygon;
! we are looking at the boundary segment from (x1,y1)
! to (x2,y2), travelling anti-clockwise.  We find the
! area between this segment and the horizontal base-line
! y = ymin, and attach a sign s1.  (Positive if the
! segment is right-to-left, negative if left to right.)
! The area of the polygon is found by summing the result
! over all boundary segments.

! Just in case you thought this wasn't complicated enough,
! what we really want is the area of the intersection of
! the polygon with the rectangular window that we're using.

! Called by dirout.
implicit none
real(dp) :: area, eps, s1, slope, tmp, w, w1, w2, x, x1, x2, xib, xit, xl, xmax, xmin, xr, y, y1, y2, ybot, yl, &
    & ymax, ymin, yr, ytop, zero
real(dp) :: rw(4)
integer :: ndi(1)
logical :: value
ndi(1) = 0

zero   = 0.e0_dp

! If the segment is vertical, the area is zero.
call testeq(x1,x2,eps,value)
if(value) then
    area = 0.e0_dp
    s1   = 0.e0_dp
    return
endif

! Find which is the right-hand end, and which is the left.
if(x1<x2) then
    xl = x1
    yl = y1
    xr = x2
    yr = y2
    s1 = -1.e0_dp
else
    xl = x2
    yl = y2
    xr = x1
    yr = y1
    s1 = 1.e0_dp
endif

! Dig out the corners of the rectangular window.
xmin = rw(1)
xmax = rw(2)
ymin = rw(3)
ymax = rw(4)

! Now start intersecting with the rectangular window.
! Truncate the segment in the horizontal direction at
! the edges of the rectangle.
slope = (yl-yr)/(xl-xr)
x  = max(xl,xmin)
y  = yl+slope*(x-xl)
xl = x
yl = y

x  = min(xr,xmax)
y  = yr+slope*(x-xr)
xr = x
yr = y

if(xr<=xmin .or. xl>=xmax) then
    area = 0.e0_dp
    return
endif

! We're now looking at a trapezoidal region which may or may
! not protrude above or below the horizontal strip bounded by
! y = ymax and y = ymin.
ybot = min(yl,yr)
ytop = max(yl,yr)

! Case 1; ymax <= ybot:
! The `roof' of the trapezoid is entirely above the
! horizontal strip.
if(ymax<=ybot) then
    area = (xr-xl)*(ymax-ymin)
    return
endif

! Case 2; ymin <= ybot <= ymax <= ytop:
! The `roof' of the trapezoid intersects the top of the
! horizontal strip (y = ymax) but not the bottom (y = ymin).
if(ymin<=ybot .and. ymax<=ytop) then
    call testeq(slope,zero,eps,value)
    if(value) then
        w1 = 0.e0_dp
        w2 = xr-xl
    else
        xit = xl+(ymax-yl)/slope
        w1 = xit-xl
        w2 = xr-xit
        if(slope<0.e0_dp) then
            tmp = w1
            w1  = w2
            w2  = tmp
        endif
    endif
    area = 0.5*w1*((ybot-ymin)+(ymax-ymin))+w2*(ymax-ymin)
    return
endif

! Case 3; ybot <= ymin <= ymax <= ytop:
! The `roof' intersects both the top (y = ymax) and
! the bottom (y = ymin) of the horizontal strip.
if(ybot<=ymin .and. ymax<=ytop) then
    xit = xl+(ymax-yl)/slope
    xib = xl+(ymin-yl)/slope
    if(slope>0.e0_dp) then
            w1 = xit-xib
            w2 = xr-xit
    else
        w1 = xib-xit
        w2 = xit-xl
    endif
    area = 0.5e0_dp*w1*(ymax-ymin)+w2*(ymax-ymin)
    return
endif

! Case 4; ymin <= ybot <= ytop <= ymax:
! The `roof' is ***between*** the bottom (y = ymin) and
! the top (y = ymax) of the horizontal strip.
if(ymin<=ybot .and. ytop<=ymax) then
    area = 0.5e0_dp*(xr-xl)*((ytop-ymin)+(ybot-ymin))
    return
endif

! Case 5; ybot <= ymin <= ytop <= ymax:
! The `roof' intersects the bottom (y = ymin) but not
! the top (y = ymax) of the horizontal strip.
if(ybot<=ymin .and. ymin<=ytop) then
    call testeq(slope,zero,eps,value)
    if(value) then
        area = 0.
        return
    endif
    xib = xl+(ymin-yl)/slope
    if(slope>0.e0_dp) then
        w = xr-xib
    else
        w = xib-xl
    endif
    area = 0.5*w*(ytop-ymin)
    return
endif

! Case 6; ytop <= ymin:
! The `roof' is entirely below the bottom (y = ymin), so
! there is no area contribution at all.
if(ytop<=ymin) then
        area = 0.
        return
endif

! Default; all stuffed up:
call intpr("Fell through all six cases.",-1,ndi,0)
call intpr("Something is totally stuffed up!",-1,ndi,0)
call intpr("Chaos and havoc in stoke.",-1,ndi,0)
call rexit("Bailing out of stoke.")
end subroutine stoke

! ---- upstream succ.f90 ----
subroutine  succ(ksc,i,j,nadj,madj,ntot)

! Find the successor of j in the adjacency list of i.
! Called by addpt, initad, trifnd, swap, delout, dirseg, dirout.
implicit none
integer :: i, j, k, kp, ksc, madj, n, ntot
integer :: nadj(-3:ntot,0:madj)
integer :: ndi(1)

! Set dummy integer for call to intpr(...).
ndi(1) = 0

n = nadj(i,0)

! If the adjacency list of i is empty, then clearly j has no successor
! in this adjacency list.  Something's wrong; stop.
if(n==0) then
    call intpr("Adjacency list of i is empty, and so cannot contain j.",-1,ndi,0)
    call rexit("Bailing out of succ.")
endif

! The adjacency list of i is non-empty; search through it until j is found;
! add 1 to the location of j, and find the contents of this new location.
do k = 1,n
    if(j==nadj(i,k)) then
        kp = k+1
        if(kp>n) kp = 1         ! Take kp modulo n. (The adjacency list
        ksc = nadj(i,kp)        ! is circular.)
        return
    endif
enddo

! The adjacency list doesn't contain j.  Something's wrong.
ndi(1) = i
call intpr("i =",-1,ndi,1)
ndi(1) = j
call intpr("j =",-1,ndi,1)
call intpr("Adjacency list of i does not contain j.",-1,ndi,0)
call rexit("Bailing out of succ.")
end subroutine succ

! ---- upstream swap.f90 ----
subroutine swap(j,k1,k2,shdswp,nadj,madj,x,y,ntot,eps,incAdj)

! The segment k1->k2 is a diagonal of a quadrilateral
! with a vertex at j (the point being added to the
! triangulation).  If the LOP is not satisfied, swap
! it for the other diagonal.
! Called by addpt.
implicit none
integer :: incadj, j, k, k1, k2, kk, madj, ntot
real(dp) :: eps
integer :: nadj(-3:ntot,0:madj)
real(dp) :: x(-3:ntot), y(-3:ntot)
logical :: shdswp
! This could happen if vertices j, k1, and k2 were colinear, but shouldn't.
call adjchk(k1,k2,shdswp,nadj,madj,ntot)
if(.not.shdswp) return

! Get the other vertex of the quadrilateral.
call pred(k,k1,k2,nadj,madj,ntot)  ! If these aren't the same, then
call succ(kk,k2,k1,nadj,madj,ntot) ! there is no other vertex.
if(kk .ne. k) then
    shdswp = .false.
    return
endif

! Check whether the LOP is satisified; i.e. whether
! vertex k is outside the circumcircle of vertices j, k1, and k2
call qtest(j,k1,k,k2,shdswp,x,y,ntot,eps)

! Do the actual swapping.
if(shdswp) then
    call delet(k1,k2,nadj,madj,ntot)
    call insrt(j,k,nadj,madj,x,y,ntot,eps,incAdj)
    if(incAdj==1) return
endif
end subroutine swap

! ---- upstream testeq.f90 ----
subroutine testeq(a,b,eps,value)

! Test for the equality of a and b in a fairly
! robust way.
! Called by trifnd, circen, stoke.
implicit none
real(dp) :: a, b, c, eps, one, ten
logical :: value
one = 1.e0_dp
ten = 1.e10_dp

! If b is essentially 0, check whether a is essentially zero also.
! The following is very sloppy!  Must fix it!
if(abs(b)<=eps) then
    if(abs(a)<=eps) then
        value = .true.
    else
        value = .false.
    endif
    return
endif

! Test if a is a `lot different' from b.  (If it is
! they're obviously not equal.)  This avoids under/overflow
! problems in dividing a by b.
if(abs(a)>ten*abs(b) .or. abs(a)<one*abs(b)) then
    value = .false.
    return
endif

! They're non-zero and fairly close; compare their ratio with 1.
c = a/b
if(abs(c-1.e0_dp)<=eps) then
    value = .true.
else
    value = .false.
endif
end subroutine testeq

! ---- upstream triar.f90 ----
subroutine triar(x0,y0,x1,y1,x2,y2,area)

! Calculate the area of a triangle with given vertices.  Called
! by delout (so that the vertices are presented in the anticlockwise
! direction).
implicit none
real(dp) :: area, half, x0, x1, x2, y0, y1, y2
half = 0.5e0_dp

area = half*((x1-x0)*(y2-y0)-(x2-x0)*(y1-y0))
end subroutine triar

! ---- upstream trifnd.f90 ----
subroutine trifnd(j,tau,nedge,nadj,madj,x,y,ntot,eps,ntri)

! Find the triangle of the extant triangulation in which
! lies the point currently being added.
! Called by initad.
implicit none
integer :: i, i1, ijk, ip, ivtmp, j, j1, k1, ktri, madj, nedge, ntau, ntot, ntri
real(dp) :: cprd, eps
integer :: nadj(-3:ntot,0:madj)
real(dp) :: x(-3:ntot), y(-3:ntot), xt(3), yt(3)
integer :: ndi(1)
integer :: tau(3)
logical :: adjace, anticl
! calling trifnd.
if(j==1) then
    call intpr("No triangles to find.",-1,ndi,0)
    call rexit("Bailing out of trifnd.")
endif

! Get the previous triangle:
j1     = j-1
tau(1) = j1
tau(3) = nadj(j1,1)
call pred(tau(2),j1,tau(3),nadj,madj,ntot)
call adjchk(tau(2),tau(3),adjace,nadj,madj,ntot)
if(.not.adjace) then
    tau(3) = tau(2)
    call pred(tau(2),j1,tau(3),nadj,madj,ntot)
endif

! Move to the adjacent triangle in the direction of the new
! point, until the new point lies in this triangle.
ktri = 0

do
! Check that the vertices of the triangle listed in tau are
! in anticlockwise order.  (If they aren't then reverse the order;
! if they are *still*  not in anticlockwise order, theh alles
! upgefucken ist; throw an error.)
call acchk(tau(1),tau(2),tau(3),anticl,x,y,ntot,eps)
if(.not.anticl) then
    call acchk(tau(3),tau(2),tau(1),anticl,x,y,ntot,eps)
    if(.not.anticl) then
        ndi(1) = j
        call intpr("Point number =",-1,ndi,1)
        call intpr("Previous triangle:",-1,tau,3)
        call intpr("Both vertex orderings are clockwise.",-1,ndi,0)
        call intpr("See help for deldir.",-1,ndi,0)
        call rexit("Bailing out of trifnd.")
    else
        ivtmp  = tau(3)
        tau(3) = tau(1)
        tau(1) = ivtmp
    endif
endif

ntau  = 0 ! This number will identify the triangle to be moved to.
nedge = 0 ! If the point lies on an edge, this number will identify that edge.
do i = 1,3
    ip = i+1
    if(ip==4) ip = 1 ! Take addition modulo 3.

! Get the coordinates of the vertices of the current side,
! and of the point j which is being added:
    xt(1) = x(tau(i))
    yt(1) = y(tau(i))
    xt(2) = x(tau(ip))
    yt(2) = y(tau(ip))
    xt(3) = x(j)
    yt(3) = y(j)

! Create indicator telling which of tau(i), tau(ip), and j
! are ideal points.  (The point being added, j, is ***never*** ideal.)
    if(tau(i)<=0) then
        i1 = 1
    else
        i1 = 0
    endif
    if(tau(ip)<=0) then
        j1 = 1
    else
        j1 = 0
    endif
    k1 = 0
    ijk = i1*4+j1*2+k1

! Calculate the ``normalized'' cross product; if this is positive
! then the point being added is to the left (as we move along the
! edge in an anti-clockwise direction).  If the test value is positive
! for all three edges, then the point is inside the triangle.  Note
! that if the test value is very close to zero, we might get negative
! values for it on both sides of an edge, and hence go into an
! infinite loop.
    call cross(xt,yt,ijk,cprd)
    if(cprd >=  eps) then
        continue
    else
        if(cprd > -eps) then
            nedge = ip
        else
            ntau = ip
            exit
        endif
    endif
enddo

! We've played ring-around-the-triangle; now figure out the
! next move:

! case 0: All tests >= 0.; the point is inside; return.
if(ntau==0) return

! The point is not inside; work out the vertices of the triangle to which
! to move.  Notation: Number the vertices of the current triangle from 1 to 3,
! anti-clockwise. Then "triangle i+1" is adjacent to the side from vertex i to
! vertex i+1, where i+1 is taken modulo 3 (i.e. "3+1 = 1").

! case 1: Move to "triangle 1"
if(ntau==1) then
    tau(2)  = tau(3)
    call succ(tau(3),tau(1),tau(2),nadj,madj,ntot)
endif

! case 2: Move to "triangle 2"
if(ntau==2) then
    tau(3)  = tau(2)
    call pred(tau(2),tau(1),tau(3),nadj,madj,ntot)
endif

! case 3: Move to "triangle 3"
if(ntau==3) then
    tau(1)  = tau(3)
    call succ(tau(3),tau(1),tau(2),nadj,madj,ntot)
endif

! We've moved to a new triangle; check if the point being added lies
! inside this one.
ktri = ktri + 1
if(ktri > ntri) then
    ndi(1) = j
    call intpr("Point being added:",-1,ndi,1)
    call intpr("Cannot find an enclosing triangle.",-1,ndi,0)
    call intpr("See help for deldir.",-1,ndi,0)
    call rexit("Bailing out of trifnd.")
endif
enddo

end subroutine trifnd

end module deldir_kernel
