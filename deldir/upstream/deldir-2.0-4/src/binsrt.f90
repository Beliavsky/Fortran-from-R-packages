subroutine binsrt(x,y,rw,nn,ind,rind,tx,ty,ilst)
! Sort the data points into bins.
! Called by master.

implicit double precision(a-h,o-z)
dimension :: x(nn), y(nn), tx(nn), ty(nn)
integer :: rind(nn)
dimension :: ind(nn), ilst(nn)
dimension :: rw(4)
dimension :: ndi(1)

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
