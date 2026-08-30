subroutine testeq(a,b,eps,value)

! Test for the equality of a and b in a fairly
! robust way.
! Called by trifnd, circen, stoke.

implicit double precision(a-h,o-z)
logical :: value

! Define constants.
one = 1.d0
ten = 1.d10

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
if(abs(c-1.d0)<=eps) then
    value = .true.
else
    value = .false.
endif
end subroutine testeq
