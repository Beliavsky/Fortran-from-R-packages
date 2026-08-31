subroutine updatefactor( m,nnzd, &
& d,jd,id, invp,perm, &
& lindx,xlindx, nsuper,lnz,xlnz, &
& snode, xsuper, &
& cachesize,ierr)

implicit none
integer m,nnzd
integer nsuper,tmpsiz, &
& ierr, &
& jd(nnzd),cachesize, &
& id(m+1),lindx(*),xlindx(*), &
& invp(m),perm(m),xlnz(m+1), &
& snode(m),xsuper(m+1)
double precision d(nnzd),lnz(*)

! temp and working stuff, loops, etc
integer iwork(7*m+3)
integer split(m)
!
! Clean L
call cleanlnz(nsuper,xsuper,xlnz,lnz)
!
! Input numerical values into data structures of L
call inpnv(id,jd,d,perm,invp,nsuper,xsuper,xlindx,lindx, &
& xlnz,lnz,iwork)
!
! Initialization for block factorization
call bfinit(m,nsuper,xsuper,snode,xlindx,lindx,cachesize,tmpsiz, &
& split)
!
! Numerical factorization
call blkfc2(nsuper,xsuper,snode,split,xlindx,lindx,xlnz, &
& lnz,iwork(1),iwork(nsuper+1),iwork(2*nsuper+1), &
& iwork(2*nsuper+m+1),tmpsiz,ierr)
if (ierr .eq. -1) then
ierr = 1
go to 100
elseif (ierr .eq. -2) then
ierr = 3
go to 100
endif
100 continue
return
end




subroutine cholstepwise(  m,nnzd,  d,jd,id, &
& doperm,invp,perm, &
& nsub,nsubmax, &
& lindx,xlindx,nsuper,nnzlmax,lnz,xlnz, &
& snode,xsuper, &
& cachsz,ierr)
!     Modified chol routine
!
!
! Sparse least squares solver via Ng-Peyton's sparse Cholesky
!    factorization for sparse symmetric positive definite
! INPUT:
!     m -- the number of column in the matrix A
!     d -- an nnzd-vector of non-zero values of A
!     jd -- an nnzd-vector of indices in d
!     id -- an (m+1)-vector of pointers to the begining of each
!           row in d and jd
!     nsub  --  [output] length of entries in lindx
!     nsubmax -- upper bound of the dimension of lindx
!     lindx -- an nsub-vector of integer which contains, in
!           column major oder, the row subscripts of the nonzero
!           entries in L in a compressed storage format
!     xlindx -- an nsuper-vector of integer of pointers for lindx
!     nsuper -- the length of xlindx   ???
!     nnzlmax -- the upper bound of the non-zero entries in
!                L stored in lnz, including the diagonal entries
!     lnz -- First contains the non-zero entries of d; later
!            contains the entries of the Cholesky factor
!     xlnz -- column pointer for L stored in lnz
!     invp -- an n-vector of integer of inverse permutation
!             vector
!     perm -- an n-vector of integer of permutation vector
!     colcnt -- array of length m, containing the number of
!               non-zeros in each column of the factor, including
!               the diagonal entries
!     snode -- array of length m for recording supernode
!              membership
!     xsuper -- array of length m+1 containing the supernode
!               partitioning
!     split -- an m-vector with splitting of supernodes so that
!              they fit into cache
!     tmpmax -- upper bound of the dimension of tmpvec
!     tmpvec -- a tmpmax-vector of temporary vector
!     cachsz -- size of the cache (in kilobytes) on the target
!               machine
!     ierr -- error flag
!       1 -- insufficient work space in call to extract
!       2 -- (not defined)
!       3 -- insufficient storage in iwork when calling sfinit;
!       4 -- nnzl > nnzlmax when calling sfinit
!       5 -- nsub > nsubmax when calling sfinit
!       6 -- insufficient work space in iwork when calling symfct
!       7 -- inconsistancy in input when calling symfct
!       8 -- tmpsiz > tmpmax when calling symfct; increase tmpmax
!       9 -- nonpositive diagonal encountered when calling
!            blkfct
!       10 -- insufficient work storage in tmpvec when calling
!            blkfct
!       11 -- insufficient work storage in iwork when calling
!            blkfct
! WORK ARRAYS:
!     adjncy -- the indices of non diag elements
!     iwsiz -- set at 7*m+3
!     iwork -- an iwsiz-vector of integer as work space
!
! FIXME: difference between nsub and nnzl ???
implicit none
integer m,nnzd,doperm
integer nsub,nsuper,nnzl, &
& nnzlmax,nsubmax,cachsz,ierr, &
& jd(nnzd), &
& id(m+1),lindx(nsubmax),xlindx(m+1), &
& invp(m),perm(m),xlnz(m+1), &
& snode(m),xsuper(m+1),split(m)
double precision d(nnzd),lnz(nnzlmax)

! local variables:
!  fix introduced in 29-3
!     &        adj(m+1),adjncy(nnzd-m)
integer adj(m+1), adjncy(nnzd-m+1), colcnt(m), &
& iwsiz,tmpsiz
! temp and working stuff, loops, etc
integer i,j,k,  nnzadj, jtmp, intmp1, intmp2
integer iwork(7*m+3)

! iwsiz is used temporalily
iwsiz=0

! Create the adjacency matrix: eliminate the diagonal elements from
!    (d,id,jd) and make two copies: (*,xlindx,lindx),(*,adj,adjncy)
!    Also to lindx and xlindx, because the matrix structure is destroyed
!    by the minimum degree ordering routine.

nsub = 0

! the adj matrix has m elements less than d
nnzadj = nnzd - m

k=1
do i=1,m
! copy id, but ajust for the missing diagonal.
xlindx(i) = id(i)-i+1
adj(i) = xlindx(i)

! now cycle over all rows
do j=id(i),id(i+1)-1
jtmp=jd(j)
if (jtmp.ne.i) then
lindx(k) = jtmp
adjncy(k) = jtmp
k=k+1
else
if ( d(j) .le. 0) then
ierr = 1
return
endif
iwsiz = iwsiz + 1
endif
enddo
enddo
jtmp=m+1
xlindx(jtmp) = id(jtmp)-m
adj(jtmp) = xlindx(jtmp)

! check if we actually had m elements on the diagonal...
if ( iwsiz .lt. m) then
ierr = 1
return
endif

! initialize iwsiz to the later used value...
iwsiz=7*m+3

!
!
! reorder the matrix using minimum degree ordering routine.
! we call the genmmd function directly.

if (doperm.eq.1) then

!       delta  - tolerance value for multiple elimination.
!                set to 0 below
!       maxint - maximum machine representable (short) integer
!                (any smaller estimate will do) for marking
!                nodes.
!                set to 32767 below
intmp1 = 0
intmp2 = 32767
call genmmd  (  m, xlindx,lindx, invp,perm,intmp1, &
& iwork(1),  iwork(m+1), iwork(2*m+1), iwork(3*m+1) , &
& intmp2, nsub   )
endif
if (doperm.eq.2) then
call genrcm ( m, nnzadj, xlindx,lindx, perm )
do i=1,m
invp(perm(i))=i
enddo
endif
if (doperm.eq.0) then
do i=1,m
invp(perm(i))=i
enddo
endif

!
! Call sfinit: Symbolic factorization initialization
!   to compute supernode partition and storage requirements
!   for symbolic factorization. New ordering is a postordering
!   of the nodal elimination tree
!
call sfinit(m,nnzadj,adj(1),adjncy(1),perm, &
& invp,colcnt,nnzl,nsub,nsuper,snode,xsuper,iwsiz, &
& iwork,ierr)
! we do not have to test ierr, as we have hardwired iwsiz to 7*m+3
if (nnzl .gt. nnzlmax) then
ierr = 4
go to 100
endif
if (nsub .gt. nsubmax) then
ierr = 5
go to 100
endif
!
! Call symfct: Perform supernodal symbolic factorization
!
iwsiz = nsuper + 2 * m + 1

call symfc2(m,nnzadj,adj(1),adjncy(1),perm,invp, &
& colcnt,nsuper,xsuper,snode,nsub,xlindx,lindx, &
& xlnz, &
& iwork(1), iwork(nsuper+1), iwork(nsuper+m+2) ,ierr)
! ierr = -2 "inconsistency in the input"
if (ierr .eq. -2) then
ierr = 6
go to 100
endif
!
! Input numerical values into data structures of L
call inpnv(id,jd,d,perm,invp,nsuper,xsuper,xlindx,lindx, &
& xlnz,lnz,iwork)
!
! Initialization for block factorization
call bfinit(m,nsuper,xsuper,snode,xlindx,lindx,cachsz,tmpsiz, &
& split)
!
! Numerical factorization
call blkfc2(nsuper,xsuper,snode,split,xlindx,lindx,xlnz, &
& lnz,iwork(1),iwork(nsuper+1),iwork(2*nsuper+1), &
& iwork(2*nsuper+m+1),tmpsiz,ierr)
if (ierr .eq. -1) then
ierr = 1
go to 100
elseif (ierr .eq. -2) then
ierr = 3
go to 100
endif
100 continue
return
end





!***********************************************************************
!***********************************************************************
!
!   Authors:        Reinhard Furrer, based on inpnv
!
!
!***********************************************************************
!***********************************************************************
!
!     ------------------------------------------------------
!     Clean the array lnz
!     ------------------------------------------------------
!
SUBROUTINE  CLEANLNZ (NSUPER, XSUPER, XLNZ, LNZ)
!
IMPLICIT NONE

INTEGER             NSUPER
INTEGER             XSUPER(*), XLNZ(*)
DOUBLE PRECISION    LNZ(*)
!
INTEGER             II, J, JSUPER
!
DO  500  JSUPER = 1, NSUPER
DO  400  J = XSUPER(JSUPER), XSUPER(JSUPER+1)-1
DO  200  II = XLNZ(J), XLNZ(J+1)-1
LNZ(II) = 0.0
200 CONTINUE
400 CONTINUE
!
500 CONTINUE
RETURN
END




!***********************************************************************
!***********************************************************************
!
!   Version:        0.4
!   Last modified:  December 27, 1994
!   Authors:        Esmond G. Ng and Barry W. Peyton
!
!   Mathematical Sciences Section, Oak Ridge National Laboratory
!
!***********************************************************************
!***********************************************************************
!************     ASSMB .... INDEXED ASSEMBLY OPERATION     ************
!***********************************************************************
!***********************************************************************
!
!   PURPOSE:
!       THIS ROUTINE PERFORMS AN INDEXED ASSEMBLY (I.E., SCATTER-ADD)
!       OPERATION, ASSUMING DATA STRUCTURES USED IN SOME OF OUR SPARSE
!       CHOLESKY CODES.
!
!   INPUT PARAMETERS:
!       M               -   NUMBER OF ROWS IN Y.
!       Q               -   NUMBER OF COLUMNS IN Y.
!       Y               -   BLOCK UPDATE TO BE INCORPORATED INTO FACTOR
!                           STORAGE.
!       RELIND          -   RELATIVE INDICES FOR MAPPING THE UPDATES
!                           ONTO THE TARGET COLUMNS.
!       XLNZ            -   POINTERS TO THE START OF EACH COLUMN IN THE
!                           TARGET MATRIX.
!
!   OUTPUT PARAMETERS:
!       LNZ             -   CONTAINS COLUMNS MODIFIED BY THE UPDATE
!                           MATRIX.
!
!***********************************************************************
!
SUBROUTINE  ASSMB  (  M     , Q     , Y     , RELIND, XLNZ  , &
& LNZ   , LDA                             )
!
!***********************************************************************
implicit none
!
!       -----------
!       PARAMETERS.
!       -----------
!
INTEGER             LDA   , M     , Q
INTEGER             XLNZ(*)
INTEGER             RELIND(*)
DOUBLE PRECISION    LNZ(*)        , Y(*)
!
!       ----------------
!       LOCAL VARIABLES.
!       ----------------
!
INTEGER             ICOL  , IL1   , IR    , IY1   , LBOT1 , &
& YCOL  , YOFF1
!
!***********************************************************************
!
!
YOFF1 = 0
IY1 = 0
DO  200  ICOL = 1, Q
YCOL = LDA - RELIND(ICOL)
LBOT1 = XLNZ(YCOL+1) - 1
!DIR$ IVDEP
DO  100  IR = ICOL, M
IL1 = LBOT1 - RELIND(IR)
IY1 = YOFF1 + IR
LNZ(IL1) = LNZ(IL1) + Y(IY1)
Y(IY1) = 0.0D0
100 CONTINUE
YOFF1 = IY1 - ICOL
200 CONTINUE
!
RETURN
END
!***********************************************************************
!***********************************************************************
!
!   Version:        0.4
!   Last modified:  December 27, 1994
!   Authors:        Joseph W.H. Liu
!
!   Mathematical Sciences Section, Oak Ridge National Laboratory
!
!***********************************************************************
!***********************************************************************
!******     BETREE ..... BINARY TREE REPRESENTATION OF ETREE     *******
!***********************************************************************
!***********************************************************************
!
!   WRITTEN BY JOSEPH LIU (JUL 17, 1985)
!
!   PURPOSE:
!       TO DETERMINE THE BINARY TREE REPRESENTATION OF THE ELIMINATION
!       TREE GIVEN BY THE PARENT VECTOR.  THE RETURNED REPRESENTATION
!       WILL BE GIVEN BY THE FIRST-SON AND BROTHER VECTORS.  THE ROOT
!       OF THE BINARY TREE IS ALWAYS NEQNS.
!
!   INPUT PARAMETERS:
!       NEQNS           -   NUMBER OF EQUATIONS.
!       PARENT          -   THE PARENT VECTOR OF THE ELIMINATION TREE.
!                           IT IS ASSUMED THAT PARENT(I) > I EXCEPT OF
!                           THE ROOTS.
!
!   OUTPUT PARAMETERS:
!       FSON            -   THE FIRST SON VECTOR.
!       BROTHR          -   THE BROTHER VECTOR.
!
!***********************************************************************
!
SUBROUTINE  BETREE (  NEQNS , PARENT, FSON  , BROTHR          )
!
!***********************************************************************
!
implicit none

INTEGER           BROTHR(*)     , FSON(*)       , &
& PARENT(*)
!
INTEGER           NEQNS
!
!***********************************************************************
!
INTEGER           LROOT , NODE  , NDPAR
!
!***********************************************************************
!
IF  ( NEQNS .LE. 0 )  RETURN
!
DO  100  NODE = 1, NEQNS
FSON(NODE) = 0
BROTHR(NODE) = 0
100 CONTINUE
LROOT = NEQNS
!       ------------------------------------------------------------
!       FOR EACH NODE := NEQNS-1 STEP -1 DOWNTO 1, DO THE FOLLOWING.
!       ------------------------------------------------------------
IF  ( NEQNS .LE. 1 )  RETURN
DO  300  NODE = NEQNS-1, 1, -1
NDPAR = PARENT(NODE)
IF  ( NDPAR .LE. 0  .OR.  NDPAR .EQ. NODE )  THEN
!               -------------------------------------------------
!               NODE HAS NO PARENT.  GIVEN STRUCTURE IS A FOREST.
!               SET NODE TO BE ONE OF THE ROOTS OF THE TREES.
!               -------------------------------------------------
BROTHR(LROOT) = NODE
LROOT = NODE
ELSE
!               -------------------------------------------
!               OTHERWISE, BECOMES FIRST SON OF ITS PARENT.
!               -------------------------------------------
BROTHR(NODE) = FSON(NDPAR)
FSON(NDPAR) = NODE
ENDIF
300 CONTINUE
BROTHR(LROOT) = 0
!
RETURN
END
!***********************************************************************
!***********************************************************************
!
!   Version:        0.4
!   Last modified:  December 27, 1994
!   Authors:        Esmond G. Ng and Barry W. Peyton
!
!   Mathematical Sciences Section, Oak Ridge National Laboratory
!
!***********************************************************************
!***********************************************************************
!******     BFINIT ..... INITIALIZATION FOR BLOCK FACTORIZATION   ******
!***********************************************************************
!***********************************************************************
!
!   PURPOSE:
!       THIS SUBROUTINE COMPUTES ITEMS NEEDED BY THE LEFT-LOOKING
!       BLOCK-TO-BLOCK CHOLESKY FACTORITZATION ROUTINE BLKFCT.
!
!   INPUT PARAMETERS:
!       NEQNS           -   NUMBER OF EQUATIONS.
!       NSUPER          -   NUMBER OF SUPERNODES.
!       XSUPER          -   INTEGER ARRAY OF SIZE (NSUPER+1) CONTAINING
!                           THE SUPERNODE PARTITIONING.
!       SNODE           -   SUPERNODE MEMBERSHIP.
!       (XLINDX,LINDX)  -   ARRAYS DESCRIBING THE SUPERNODAL STRUCTURE.
!       CACHSZ          -   CACHE SIZE (IN KBYTES).
!
!   OUTPUT PARAMETERS:
!       TMPSIZ          -   SIZE OF WORKING STORAGE REQUIRED BY BLKFCT.
!       SPLIT           -   SPLITTING OF SUPERNODES SO THAT THEY FIT
!                           INTO CACHE.
!
!***********************************************************************
!
SUBROUTINE  BFINIT  ( NEQNS , NSUPER, XSUPER, SNODE , XLINDX, &
& LINDX , CACHSZ, TMPSIZ, SPLIT           )
!
!***********************************************************************
implicit none
!
INTEGER     CACHSZ, NEQNS , NSUPER, TMPSIZ
INTEGER     XLINDX(*)       , XSUPER(*)
INTEGER     LINDX (*)       , SNODE (*)   , &
& SPLIT(*)
!
!***********************************************************************
!
!       ---------------------------------------------------
!       DETERMINE FLOATING POINT WORKING SPACE REQUIREMENT.
!       ---------------------------------------------------
CALL  FNTSIZ (  NSUPER, XSUPER, SNODE , XLINDX, LINDX , &
& TMPSIZ                                  )
!
!       -------------------------------
!       PARTITION SUPERNODES FOR CACHE.
!       -------------------------------
CALL  FNSPLT (  NEQNS , NSUPER, XSUPER, XLINDX, CACHSZ, &
& SPLIT                                   )
!
RETURN
END
!***********************************************************************
!***********************************************************************
!
!   Version:        0.3
!   Last modified:  March 6, 1995
!   Authors:        Esmond G. Ng and Barry W. Peyton
!                   RF eliminated dependence on SMXPY and MMPY
!
!   Mathematical Sciences Section, Oak Ridge National Laboratoy
!
!***********************************************************************
!***********************************************************************
!*********     BLKFC2 .....  BLOCK GENERAL SPARSE CHOLESKY     *********
!***********************************************************************
!***********************************************************************
!
!   PURPOSE:
!       THIS SUBROUTINE FACTORS A SPARSE POSITIVE DEFINITE MATRIX.
!       THE COMPUTATION IS ORGANIZED AROUND KERNELS THAT PERFORM
!       SUPERNODE-TO-SUPERNODE UPDATES, I.E., BLOCK-TO-BLOCK UPDATES.
!
!   INPUT PARAMETERS:
!       NSUPER          -   NUMBER OF SUPERNODES.
!       XSUPER          -   SUPERNODE PARTITION.
!       SNODE           -   MAPS EACH COLUMN TO THE SUPERNODE CONTAINING
!                           IT.
!       SPLIT           -   SPLITTING OF SUPERNODES SO THAT THEY FIT
!                           INTO CACHE.
!       (XLINDX,LINDX)  -   ROW INDICES FOR EACH SUPERNODE (INCLUDING
!                           THE DIAGONAL ELEMENTS).
!       (XLNZ,LNZ)      -   ON INPUT, CONTAINS MATRIX TO BE FACTORED.
!       TMPSIZ          -   SIZE OF TEMPORARY WORKING STORAGE.
!
!   OUTPUT PARAMETERS:
!       LNZ             -   ON OUTPUT, CONTAINS CHOLESKY FACTOR.
!       IFLAG           -   ERROR FLAG.
!                               0: SUCCESSFUL FACTORIZATION.
!                              -1: NONPOSITIVE DIAGONAL ENCOUNTERED,
!                                  MATRIX IS NOT POSITIVE DEFINITE.
!                              -2: INSUFFICIENT WORKING STORAGE
!                                  [TEMP(*)].
!
!   WORKING PARAMETERS:
!       LINK            -   LINKS TOGETHER THE SUPERNODES IN A SUPERNODE
!                           ROW.
!       LENGTH          -   LENGTH OF THE ACTIVE PORTION OF EACH
!                           SUPERNODE.
!       INDMAP          -   VECTOR OF SIZE NEQNS INTO WHICH THE GLOBAL
!                           INDICES ARE SCATTERED.
!       RELIND          -   MAPS LOCATIONS IN THE UPDATING COLUMNS TO
!                           THE CORRESPONDING LOCATIONS IN THE UPDATED
!                           COLUMNS.  (RELIND IS GATHERED FROM INDMAP).
!       TEMP            -   REAL VECTOR FOR ACCUMULATING UPDATES.  MUST
!                           ACCOMODATE ALL COLUMNS OF A SUPERNODE.
!
!***********************************************************************
!
SUBROUTINE  BLKFC2 (  NSUPER, XSUPER, SNODE , SPLIT , XLINDX, &
& LINDX , XLNZ  , LNZ   , LINK  , LENGTH, &
& INDMAP, RELIND, TMPSIZ, IFLAG )
!
!*********************************************************************
implicit none
!
!       -----------
!       PARAMETERS.
!       -----------
!
INTEGER             XLINDX(*)     , XLNZ(*)
INTEGER             INDMAP(*)     , LENGTH(*)     , &
& LINDX(*)      , LINK(*)       , &
& RELIND(*)     , SNODE(*)      , &
& SPLIT(*)      , XSUPER(*)
INTEGER             IFLAG , NSUPER, TMPSIZ
DOUBLE PRECISION    LNZ(*)
!
!       ----------------
!       LOCAL VARIABLES.
!       ----------------
!
INTEGER             FJCOL , FKCOL , I     , ILEN  , ILPNT , &
& INDDIF, JLEN  , JLPNT , JSUP  , JXPNT , &
& KFIRST, KLAST , KLEN  , KLPNT , KSUP  , &
& KXPNT , LJCOL , NCOLUP, NJCOLS, NKCOLS, &
& NXKSUP, NXTCOL, NXTSUP, STORE
DOUBLE PRECISION    TEMP(TMPSIZ)
!     RF: put TEMP(*) into a local variable

DOUBLE PRECISION MXDIAG
INTEGER NTINY
!*********************************************************************
!
IFLAG = 0
NTINY = 0
NXTCOL = 0
!
!       -----------------------------------------------------------
!       INITIALIZE EMPTY ROW LISTS IN LINK(*) AND ZERO OUT TEMP(*).
!       -----------------------------------------------------------
DO  100  JSUP = 1, NSUPER
LINK(JSUP) = 0
100 CONTINUE
DO  200  I = 1, TMPSIZ
TEMP(I) = 0.0D+00
200 CONTINUE

!       COMPUTE MAXIMUM DIAGONAL ELEMENT IN INPUT MATRIX
MXDIAG = 0.D0
DO 201 I = 1, XSUPER(NSUPER+1)-1
FJCOL = XLNZ(I)
MXDIAG = MAX(MXDIAG, LNZ(FJCOL))
201 CONTINUE
!
!       ---------------------------
!       FOR EACH SUPERNODE JSUP ...
!       ---------------------------
DO  600  JSUP = 1, NSUPER
!
!           ------------------------------------------------
!           FJCOL  ...  FIRST COLUMN OF SUPERNODE JSUP.
!           LJCOL  ...  LAST COLUMN OF SUPERNODE JSUP.
!           NJCOLS ...  NUMBER OF COLUMNS IN SUPERNODE JSUP.
!           JLEN   ...  LENGTH OF COLUMN FJCOL.
!           JXPNT  ...  POINTER TO INDEX OF FIRST
!                       NONZERO IN COLUMN FJCOL.
!           ------------------------------------------------
FJCOL  = XSUPER(JSUP)
NJCOLS = XSUPER(JSUP+1) - FJCOL
LJCOL  = FJCOL + NJCOLS - 1
JLEN   = XLNZ(FJCOL+1) - XLNZ(FJCOL)
JXPNT  = XLINDX(JSUP)

!            print *, 'Super Node: ', JSUP, ' first: ', FJCOL,
!     .           ' last: ', LJCOL


!
!
!           -----------------------------------------------------
!           SET UP INDMAP(*) TO MAP THE ENTRIES IN UPDATE COLUMNS
!           TO THEIR CORRESPONDING POSITIONS IN UPDATED COLUMNS,
!           RELATIVE THE THE BOTTOM OF EACH UPDATED COLUMN.
!           -----------------------------------------------------
CALL  LDINDX ( JLEN, LINDX(JXPNT), INDMAP )
!
!           -----------------------------------------
!           FOR EVERY SUPERNODE KSUP IN ROW(JSUP) ...
!           -----------------------------------------
KSUP = LINK(JSUP)
300 IF  ( KSUP .GT. 0 )  THEN
NXKSUP = LINK(KSUP)
!
!               -------------------------------------------------------
!               GET INFO ABOUT THE CMOD(JSUP,KSUP) UPDATE.
!
!               FKCOL  ...  FIRST COLUMN OF SUPERNODE KSUP.
!               NKCOLS ...  NUMBER OF COLUMNS IN SUPERNODE KSUP.
!               KLEN   ...  LENGTH OF ACTIVE PORTION OF COLUMN FKCOL.
!               KXPNT  ...  POINTER TO INDEX OF FIRST NONZERO IN ACTIVE
!                           PORTION OF COLUMN FJCOL.
!               -------------------------------------------------------
FKCOL = XSUPER(KSUP)
NKCOLS = XSUPER(KSUP+1) - FKCOL
KLEN = LENGTH(KSUP)
KXPNT = XLINDX(KSUP+1) - KLEN
!
!               -------------------------------------------
!               PERFORM CMOD(JSUP,KSUP), WITH SPECIAL CASES
!               HANDLED DIFFERENTLY.
!               -------------------------------------------
!
IF  ( KLEN .NE. JLEN )  THEN
!
!                   -------------------------------------------
!                   SPARSE CMOD(JSUP,KSUP).
!
!                   NCOLUP ... NUMBER OF COLUMNS TO BE UPDATED.
!                   -------------------------------------------
!
DO  400  I = 0, KLEN-1
NXTCOL = LINDX(KXPNT+I)
IF  ( NXTCOL .GT. LJCOL )  GO TO 500
400 CONTINUE
I = KLEN
500 CONTINUE
NCOLUP = I
!
IF  ( NKCOLS .EQ. 1 )  THEN
!
!                       ----------------------------------------------
!                       UPDATING TARGET SUPERNODE BY TRIVIAL
!                       SUPERNODE (WITH ONE COLUMN).
!
!                       KLPNT  ...  POINTER TO FIRST NONZERO IN ACTIVE
!                                   PORTION OF COLUMN FKCOL.
!                       ----------------------------------------------
KLPNT = XLNZ(FKCOL+1) - KLEN
CALL  MMPYI ( KLEN, NCOLUP, LINDX(KXPNT), &
& LNZ(KLPNT), XLNZ, LNZ, INDMAP )
!
ELSE
!
!                       --------------------------------------------
!                       KFIRST ...  FIRST INDEX OF ACTIVE PORTION OF
!                                   SUPERNODE KSUP (FIRST COLUMN TO
!                                   BE UPDATED).
!                       KLAST  ...  LAST INDEX OF ACTIVE PORTION OF
!                                   SUPERNODE KSUP.
!                       --------------------------------------------
!
KFIRST = LINDX(KXPNT)
KLAST  = LINDX(KXPNT+KLEN-1)
INDDIF = INDMAP(KFIRST) - INDMAP(KLAST)
!
IF  ( INDDIF .LT. KLEN )  THEN
!
!                           ---------------------------------------
!                           DENSE CMOD(JSUP,KSUP).
!
!                           ILPNT  ...  POINTER TO FIRST NONZERO IN
!                                       COLUMN KFIRST.
!                           ILEN   ...  LENGTH OF COLUMN KFIRST.
!                           ---------------------------------------
ILPNT = XLNZ(KFIRST)
ILEN = XLNZ(KFIRST+1) - ILPNT
CALL  MMPY ( KLEN, NKCOLS, NCOLUP, &
& SPLIT(FKCOL), XLNZ(FKCOL), &
& LNZ, LNZ(ILPNT), ILEN  )
!
ELSE
!
!                           -------------------------------
!                           GENERAL SPARSE CMOD(JSUP,KSUP).
!                           COMPUTE CMOD(JSUP,KSUP) UPDATE
!                           IN WORK STORAGE.
!                           -------------------------------
STORE = KLEN * NCOLUP - NCOLUP * &
& (NCOLUP-1) / 2
IF  ( STORE .GT. TMPSIZ )  THEN
IFLAG = -2
RETURN
ENDIF
CALL  MMPY ( KLEN, NKCOLS, NCOLUP, &
& SPLIT(FKCOL), XLNZ(FKCOL), &
& LNZ, TEMP, KLEN  )
!                           ----------------------------------------
!                           GATHER INDICES OF KSUP RELATIVE TO JSUP.
!                           ----------------------------------------
CALL  IGATHR ( KLEN, LINDX(KXPNT), &
& INDMAP, RELIND )
!                           --------------------------------------
!                           INCORPORATE THE CMOD(JSUP,KSUP) BLOCK
!                           UPDATE INTO THE TO APPROPRIATE COLUMNS
!                           OF L.
!                           --------------------------------------
CALL  ASSMB ( KLEN, NCOLUP, TEMP, RELIND, &
& XLNZ(FJCOL), LNZ, JLEN )
!
ENDIF
!
ENDIF
!
ELSE
!
!                   ----------------------------------------------
!                   DENSE CMOD(JSUP,KSUP).
!                   JSUP AND KSUP HAVE IDENTICAL STRUCTURE.
!
!                   JLPNT  ...  POINTER TO FIRST NONZERO IN COLUMN
!                               FJCOL.
!                   ----------------------------------------------
JLPNT = XLNZ(FJCOL)
CALL  MMPY ( KLEN, NKCOLS, NJCOLS, SPLIT(FKCOL), &
& XLNZ(FKCOL), LNZ, LNZ(JLPNT), JLEN)
NCOLUP = NJCOLS
IF  ( KLEN .GT. NJCOLS )  THEN
NXTCOL = LINDX(JXPNT+NJCOLS)
ENDIF
!
ENDIF
!
!               ------------------------------------------------
!               LINK KSUP INTO LINKED LIST OF THE NEXT SUPERNODE
!               IT WILL UPDATE AND DECREMENT KSUP'S ACTIVE
!               LENGTH.
!               ------------------------------------------------
IF  ( KLEN .GT. NCOLUP )  THEN
NXTSUP = SNODE(NXTCOL)
LINK(KSUP) = LINK(NXTSUP)
LINK(NXTSUP) = KSUP
LENGTH(KSUP) = KLEN - NCOLUP
ELSE
LENGTH(KSUP) = 0
ENDIF
!
!               -------------------------------
!               NEXT UPDATING SUPERNODE (KSUP).
!               -------------------------------
KSUP = NXKSUP
GO TO 300
!
ENDIF
!
!           ----------------------------------------------
!           APPLY PARTIAL CHOLESKY TO THE COLUMNS OF JSUP.
!           ----------------------------------------------
!xPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPC
CALL CHLSUP ( JLEN, NJCOLS, SPLIT(FJCOL), XLNZ(FJCOL), LNZ, &
& MXDIAG, NTINY)
!     &                    MXDIAG, NTINY, IFLAG )
IF  ( IFLAG .NE. 0 )  THEN
IFLAG = -1
RETURN
ENDIF
!
!           -----------------------------------------------
!           INSERT JSUP INTO LINKED LIST OF FIRST SUPERNODE
!           IT WILL UPDATE.
!           -----------------------------------------------
IF  ( JLEN .GT. NJCOLS )  THEN
NXTCOL = LINDX(JXPNT+NJCOLS)
NXTSUP = SNODE(NXTCOL)
LINK(JSUP) = LINK(NXTSUP)
LINK(NXTSUP) = JSUP
LENGTH(JSUP) = JLEN - NJCOLS
ELSE
LENGTH(JSUP) = 0
ENDIF
!
600 CONTINUE
!
!xPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPC
!        IF(NTINY .NE. 0) WRITE(6,699) NTINY
! 699    FORMAT(1X,' FOUND ',I6,' TINY DIAGONALS; REPLACED WITH INF')
!
! SET IFLAG TO -1 TO INDICATE PRESENCE OF TINY DIAGONALS
!
IF(NTINY .NE. 0) IFLAG = -1
!xPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPC
RETURN
END
!***********************************************************************
!***********************************************************************
!
!   Written:        October 6, 1996 by SJW. Based on routine BLKSLV of
!                   Esmond G. Ng and Barry W. Peyton.
!
!   Modified:       Sept 30, 1999 to improve efficiency in the case
!                   in which the right-hand side and solution are both
!                   expected to be sparse. Happens a lot in "dense"
!                   column handling.
!
!***********************************************************************
!***********************************************************************
!*********     BLKSLB ... BACK TRIANGULAR SUBSTITUTION        **********
!***********************************************************************
!***********************************************************************
!
!   PURPOSE:
!       GIVEN THE CHOLESKY FACTORIZATION OF A SPARSE SYMMETRIC
!       POSITIVE DEFINITE MATRIX, THIS SUBROUTINE PERFORMS THE
!       BACKWARD TRIANGULAR SUBSTITUTION.  IT USES OUTPUT FROM BLKFCT.
!
!   INPUT PARAMETERS:
!       NSUPER          -   NUMBER OF SUPERNODES.
!       XSUPER          -   SUPERNODE PARTITION.
!       (XLINDX,LINDX)  -   ROW INDICES FOR EACH SUPERNODE.
!       (XLNZ,LNZ)      -   CHOLESKY FACTOR.
!
!   UPDATED PARAMETERS:
!       RHS             -   ON INPUT, CONTAINS THE RIGHT HAND SIDE.  ON
!                           OUTPUT, CONTAINS THE SOLUTION.
!
!***********************************************************************
!
SUBROUTINE  BLKSLB (  NSUPER, XSUPER, XLINDX, LINDX , XLNZ  , &
& LNZ   , RHS                             )
!
!***********************************************************************
implicit none
!
INTEGER             NSUPER
INTEGER             LINDX(*)      , XSUPER(*)
INTEGER             XLINDX(*)     , XLNZ(*)
DOUBLE PRECISION    LNZ(*)        , RHS(*)
!
!***********************************************************************
!
INTEGER             FJCOL , I     , IPNT  , IX    , IXSTOP, &
& IXSTRT, JCOL  , JPNT  , JSUP  , LJCOL
DOUBLE PRECISION    T
!
!***********************************************************************
!
IF  ( NSUPER .LE. 0 )  RETURN
!       -------------------------
!       BACKWARD SUBSTITUTION ...
!       -------------------------
LJCOL = XSUPER(NSUPER+1) - 1
DO  600  JSUP = NSUPER, 1, -1
FJCOL  = XSUPER(JSUP)
IXSTOP = XLNZ(LJCOL+1) - 1
JPNT   = XLINDX(JSUP) + (LJCOL - FJCOL)
DO  500  JCOL = LJCOL, FJCOL, -1
IXSTRT = XLNZ(JCOL)
IPNT   = JPNT + 1
T      = RHS(JCOL)
!DIR$           IVDEP
DO  400  IX = IXSTRT+1, IXSTOP
I = LINDX(IPNT)
IF (abs(RHS(I)) .GT. 0.D0) T = T - LNZ(IX)*RHS(I)
IPNT = IPNT + 1
400 CONTINUE

IF(abs(T) .GT. 0.D0) THEN
RHS(JCOL) = T/LNZ(IXSTRT)
ELSE
RHS(JCOL) = 0.D0
ENDIF

IXSTOP    = IXSTRT - 1
JPNT      = JPNT - 1
500 CONTINUE

LJCOL = FJCOL - 1
600 CONTINUE
!
RETURN
END
!***********************************************************************
!***********************************************************************
!
!   Written:        October 6, 1996 by SJW. Based on routine BLKSLV of
!                   Esmond G. Ng and Barry W. Peyton.
!
!   Modified:       Sept 30, 1999 to improve efficiency in the case
!                   in which the right-hand side and solution are both
!                   expected to be sparse. Happens a lot in "dense"
!                   column handling.
!
!***********************************************************************
!***********************************************************************
!*********     BLKSLF ... FORWARD TRIANGULAR SUBSTITUTION     **********
!***********************************************************************
!***********************************************************************
!
!   PURPOSE:
!       GIVEN THE CHOLESKY FACTORIZATION OF A SPARSE SYMMETRIC
!       POSITIVE DEFINITE MATRIX, THIS SUBROUTINE PERFORMS THE
!       FORWARD TRIANGULAR SUBSTITUTIOn.  IT USES OUTPUT FROM BLKFCT.
!
!   INPUT PARAMETERS:
!       NSUPER          -   NUMBER OF SUPERNODES.
!       XSUPER          -   SUPERNODE PARTITION.
!       (XLINDX,LINDX)  -   ROW INDICES FOR EACH SUPERNODE.
!       (XLNZ,LNZ)      -   CHOLESKY FACTOR.
!
!   UPDATED PARAMETERS:
!       RHS             -   ON INPUT, CONTAINS THE RIGHT HAND SIDE.  ON
!                           OUTPUT, CONTAINS THE SOLUTION.
!
!***********************************************************************
!
SUBROUTINE  BLKSLF (  NSUPER, XSUPER, XLINDX, LINDX , XLNZ  , &
& LNZ   , RHS                             )
!
!***********************************************************************
!
implicit none

INTEGER             NSUPER
INTEGER             LINDX(*)      , XSUPER(*)
INTEGER             XLINDX(*)     , XLNZ(*)
DOUBLE PRECISION    LNZ(*)        , RHS(*)
!
!***********************************************************************
!
INTEGER             FJCOL , I     , IPNT  , IX    , IXSTOP, &
& IXSTRT, JCOL  , JPNT  , JSUP  , LJCOL
DOUBLE PRECISION    T
!
!***********************************************************************
!
IF  ( NSUPER .LE. 0 )  RETURN
!
!       ------------------------
!       FORWARD SUBSTITUTION ...
!       ------------------------
FJCOL = XSUPER(1)
DO  300  JSUP = 1, NSUPER
LJCOL  = XSUPER(JSUP+1) - 1
IXSTRT = XLNZ(FJCOL)
JPNT   = XLINDX(JSUP)
DO  200  JCOL = FJCOL, LJCOL
IXSTOP    = XLNZ(JCOL+1) - 1

IF (abs(RHS(JCOL)) .GT. 0.D0) THEN
T         = RHS(JCOL)/LNZ(IXSTRT)
RHS(JCOL) = T
IPNT      = JPNT + 1
!DIR$           IVDEP
DO  100  IX = IXSTRT+1, IXSTOP
I      = LINDX(IPNT)
RHS(I) = RHS(I) - T*LNZ(IX)
IPNT   = IPNT + 1
100 CONTINUE
ENDIF

IXSTRT = IXSTOP + 1
JPNT   = JPNT + 1
200 CONTINUE
FJCOL = LJCOL + 1
300 CONTINUE
!
RETURN
END
!***********************************************************************
!***********************************************************************
!
!   Version:        0.4
!   Last modified:  December 27, 1994
!   Authors:        Esmond G. Ng and Barry W. Peyton
!
!   Mathematical Sciences Section, Oak Ridge National Laboratory
!
!   Modified:       Sept 30, 1999 to improve efficiency in the case
!                   in which the right-hand side and solution are both
!                   expected to be sparse. Happens a lot in "dense"
!                   column handling.
!
!***********************************************************************
!***********************************************************************
!*********     BLKSLV ... BLOCK TRIANGULAR SOLUTIONS          **********
!***********************************************************************
!***********************************************************************
!
!   PURPOSE:
!       GIVEN THE CHOLESKY FACTORIZATION OF A SPARSE SYMMETRIC
!       POSITIVE DEFINITE MATRIX, THIS SUBROUTINE PERFORMS THE
!       TRIANGULAR SOLUTION.  IT USES OUTPUT FROM BLKFCT.
!
!   INPUT PARAMETERS:
!       NSUPER          -   NUMBER OF SUPERNODES.
!       XSUPER          -   SUPERNODE PARTITION.
!       (XLINDX,LINDX)  -   ROW INDICES FOR EACH SUPERNODE.
!       (XLNZ,LNZ)      -   CHOLESKY FACTOR.
!
!   UPDATED PARAMETERS:
!       RHS             -   ON INPUT, CONTAINS THE RIGHT HAND SIDE.  ON
!                           OUTPUT, CONTAINS THE SOLUTION.
!
!***********************************************************************
!
SUBROUTINE  BLKSLV (  NSUPER, XSUPER, XLINDX, LINDX , XLNZ  , &
& LNZ   , RHS                             )
!
!***********************************************************************
!
implicit none

INTEGER             NSUPER
INTEGER             LINDX(*)      , XSUPER(*)
INTEGER             XLINDX(*)     , XLNZ(*)
DOUBLE PRECISION    LNZ(*)        , RHS(*)
!
!***********************************************************************
!
INTEGER             FJCOL , I     , IPNT  , IX    , IXSTOP, &
& IXSTRT, JCOL  , JPNT  , JSUP  , LJCOL
DOUBLE PRECISION    T
!
!***********************************************************************
!
IF  ( NSUPER .LE. 0 )  RETURN
!
!       ------------------------
!       FORWARD SUBSTITUTION ...
!       ------------------------
FJCOL = XSUPER(1)
DO  300  JSUP = 1, NSUPER
LJCOL  = XSUPER(JSUP+1) - 1
IXSTRT = XLNZ(FJCOL)
JPNT   = XLINDX(JSUP)
DO  200  JCOL = FJCOL, LJCOL
IXSTOP    = XLNZ(JCOL+1) - 1

IF (abs(RHS(JCOL)) .GT. 0.D0) THEN
T         = RHS(JCOL)/LNZ(IXSTRT)
RHS(JCOL) = T
IPNT      = JPNT + 1
!DIR$           IVDEP
DO  100  IX = IXSTRT+1, IXSTOP
I      = LINDX(IPNT)
RHS(I) = RHS(I) - T*LNZ(IX)
IPNT   = IPNT + 1
100 CONTINUE
ENDIF

IXSTRT = IXSTOP + 1
JPNT   = JPNT + 1
200 CONTINUE
FJCOL = LJCOL + 1
300 CONTINUE
!
!       -------------------------
!       BACKWARD SUBSTITUTION ...
!       -------------------------
LJCOL = XSUPER(NSUPER+1) - 1
DO  600  JSUP = NSUPER, 1, -1
FJCOL  = XSUPER(JSUP)
IXSTOP = XLNZ(LJCOL+1) - 1
JPNT   = XLINDX(JSUP) + (LJCOL - FJCOL)
DO  500  JCOL = LJCOL, FJCOL, -1
IXSTRT = XLNZ(JCOL)
IPNT   = JPNT + 1
T      = RHS(JCOL)
!DIR$           IVDEP
DO  400  IX = IXSTRT+1, IXSTOP
I    = LINDX(IPNT)
IF(abs(RHS(I)) .GT. 0.D0) T = T - LNZ(IX)*RHS(I)
IPNT = IPNT + 1
400 CONTINUE

IF (abs(T) .GT. 0.D0) THEN
RHS(JCOL) = T/LNZ(IXSTRT)
ELSE
RHS(JCOL) = 0.D0
ENDIF

IXSTOP    = IXSTRT - 1
JPNT      = JPNT - 1
500 CONTINUE
LJCOL = FJCOL - 1
600 CONTINUE
!
RETURN
END





!***********************************************************************
!***********************************************************************
!
!   Version:        0.4
!   Last modified:  January 12, 1995
!   Authors:        Esmond G. Ng and Barry W. Peyton
!
!   Mathematical Sciences Section, Oak Ridge National Laboratory
!
!***********************************************************************
!***********************************************************************
!******     BTREE2 ..... BINARY TREE REPRESENTATION OF ETREE     *******
!***********************************************************************
!***********************************************************************
!
!   PURPOSE:
!       TO DETERMINE A BINARY TREE REPRESENTATION OF THE ELIMINATION
!       TREE, FOR WHICH EVERY "LAST CHILD" HAS THE MAXIMUM POSSIBLE
!       COLUMN NONZERO COUNT IN THE FACTOR.  THE RETURNED REPRESENTATION
!       WILL BE GIVEN BY THE FIRST-SON AND BROTHER VECTORS.  THE ROOT OF
!       THE BINARY TREE IS ALWAYS NEQNS.
!
!   INPUT PARAMETERS:
!       NEQNS           -   NUMBER OF EQUATIONS.
!       PARENT          -   THE PARENT VECTOR OF THE ELIMINATION TREE.
!                           IT IS ASSUMED THAT PARENT(I) > I EXCEPT OF
!                           THE ROOTS.
!       COLCNT          -   COLUMN NONZERO COUNTS OF THE FACTOR.
!
!   OUTPUT PARAMETERS:
!       FSON            -   THE FIRST SON VECTOR.
!       BROTHR          -   THE BROTHER VECTOR.
!
!   WORKING PARAMETERS:
!       LSON            -   LAST SON VECTOR.
!
!***********************************************************************
!
SUBROUTINE  BTREE2 (  NEQNS , PARENT, COLCNT, FSON  , BROTHR, &
& LSON    )
!
implicit none
!***********************************************************************
!
INTEGER             BROTHR(*)     , COLCNT(*)     , &
& FSON(*)       , LSON(*)       , &
& PARENT(*)
!
INTEGER             NEQNS
!
!***********************************************************************
!
INTEGER           LROOT , NODE  , NDLSON, NDPAR
!
!***********************************************************************
!
IF  ( NEQNS .LE. 0 )  RETURN
!
DO  100  NODE = 1, NEQNS
FSON(NODE) = 0
BROTHR(NODE) = 0
LSON(NODE) = 0
100 CONTINUE
LROOT = NEQNS
!       ------------------------------------------------------------
!       FOR EACH NODE := NEQNS-1 STEP -1 DOWNTO 1, DO THE FOLLOWING.
!       ------------------------------------------------------------
IF  ( NEQNS .LE. 1 )  RETURN
DO  300  NODE = NEQNS-1, 1, -1
NDPAR = PARENT(NODE)
IF  ( NDPAR .LE. 0  .OR.  NDPAR .EQ. NODE )  THEN
!               -------------------------------------------------
!               NODE HAS NO PARENT.  GIVEN STRUCTURE IS A FOREST.
!               SET NODE TO BE ONE OF THE ROOTS OF THE TREES.
!               -------------------------------------------------
BROTHR(LROOT) = NODE
LROOT = NODE
ELSE
!               -------------------------------------------
!               OTHERWISE, BECOMES FIRST SON OF ITS PARENT.
!               -------------------------------------------
NDLSON = LSON(NDPAR)
IF  ( NDLSON .NE. 0 )  THEN
IF  ( COLCNT(NODE) .GE. COLCNT(NDLSON) )  THEN
BROTHR(NODE) = FSON(NDPAR)
FSON(NDPAR) = NODE
ELSE
BROTHR(NDLSON) = NODE
LSON(NDPAR) = NODE
ENDIF
ELSE
FSON(NDPAR) = NODE
LSON(NDPAR) = NODE
ENDIF
ENDIF
300 CONTINUE
BROTHR(LROOT) = 0
!
RETURN
END
!***********************************************************************
!***********************************************************************
!
!   Version:        0.3
!   Last modified:  December 27, 1994
!   Authors:        Esmond G. Ng and Barry W. Peyton
!
!   Modified by RF:   Eliminated the  MMPYN, SMXPY as arguments
!
!   Mathematical Sciences Section, Oak Ridge National Laboratoy
!
!***********************************************************************
!***********************************************************************
!******     CHLSUP .... DENSE CHOLESKY WITHIN SUPERNODE   **************
!***********************************************************************
!***********************************************************************
!
!     PURPOSE - THIS ROUTINE PERFORMS CHOLESKY
!               FACTORIZATION ON THE COLUMNS OF A SUPERNODE
!               THAT HAVE RECEIVED ALL UPDATES FROM COLUMNS
!               EXTERNAL TO THE SUPERNODE.
!
!     INPUT PARAMETERS -
!        M      - NUMBER OF ROWS (LENGTH OF THE FIRST COLUMN).
!        N      - NUMBER OF COLUMNS IN THE SUPERNODE.
!        XPNT   - XPNT(J+1) POINTS ONE LOCATION BEYOND THE END
!                 OF THE J-TH COLUMN OF THE SUPERNODE.
!        X(*)   - CONTAINS THE COLUMNS OF OF THE SUPERNODE TO
!                 BE FACTORED.
!
!     EXTERNAL ROUTINES -
!        MMPY8  -  MATRIX-MATRIX MULTIPLY WITH 8 LOOP UNROLLING.
!
!     OUTPUT PARAMETERS -
!        X(*)   - ON OUTPUT, CONTAINS THE FACTORED COLUMNS OF
!                 THE SUPERNODE.
!        IFLAG  - UNCHANGED IF THERE IS NO ERROR.
!                 =1 IF NONPOSITIVE DIAGONAL ENTRY IS ENCOUNTERED.
!              RF: removed!
!
!***********************************************************************
!
SUBROUTINE  CHLSUP  ( M, N, SPLIT, XPNT, X, MXDIAG, NTINY)
!     &                      IFLAG )
!
!***********************************************************************
!
implicit none
!     -----------
!     PARAMETERS.
!     -----------
!
EXTERNAL            MMPY8
!
INTEGER             M, N
!
INTEGER             XPNT(*), SPLIT(*)
!
DOUBLE PRECISION    X(*), MXDIAG
INTEGER             NTINY
!
!     ----------------
!     LOCAL VARIABLES.
!     ----------------
!
INTEGER             FSTCOL, JBLK  , JPNT  , MM    , NN    , &
& NXTCOL, Q
!
!***********************************************************************
!
JBLK = 0
FSTCOL = 1
MM = M
JPNT = XPNT(FSTCOL)
!
!       ----------------------------------------
!       FOR EACH BLOCK JBLK IN THE SUPERNODE ...
!       ----------------------------------------
100 CONTINUE
IF  ( FSTCOL .LE. N )  THEN
JBLK = JBLK + 1
NN = SPLIT(JBLK)
!           ------------------------------------------
!           ... PERFORM PARTIAL CHOLESKY FACTORIZATION
!               ON THE BLOCK.
!           ------------------------------------------
CALL PCHOL ( MM, NN, XPNT(FSTCOL), X, MXDIAG, NTINY)
!           ----------------------------------------------
!           ... APPLY THE COLUMNS IN JBLK TO ANY COLUMNS
!               OF THE SUPERNODE REMAINING TO BE COMPUTED.
!           ----------------------------------------------
NXTCOL = FSTCOL + NN
Q = N - NXTCOL + 1
MM = MM - NN
JPNT = XPNT(NXTCOL)
IF  ( Q .GT. 0 )  THEN
CALL  MMPY8( MM, NN, Q, XPNT(FSTCOL), X, X(JPNT), MM )
ENDIF
FSTCOL = NXTCOL
GO TO 100
ENDIF
!
RETURN
END
!***********************************************************************
!***********************************************************************
!
!   Version:        0.4
!   Last modified:  December 27, 1994
!   Authors:        Esmond G. Ng and Barry W. Peyton
!
!   Mathematical Sciences Section, Oak Ridge National Laboratory
!
!***********************************************************************
!***********************************************************************
!**********     CHORDR ..... CHILD REORDERING                ***********
!***********************************************************************
!***********************************************************************
!
!   PURPOSE:
!       REARRANGE THE CHILDREN OF EACH VERTEX SO THAT THE LAST ONE
!       MAXIMIZES (AMONG THE CHILDREN) THE NUMBER OF NONZEROS IN THE
!       CORRESPONDING COLUMN OF L.  ALSO DETERMINE AN NEW POSTORDERING
!       BASED ON THE STRUCTURE OF THE MODIFIED ELIMINATION TREE.
!
!   INPUT PARAMETERS:
!       NEQNS           -   NUMBER OF EQUATIONS.
!
!   UPDATED PARAMETERS:
!       (PERM,INVP)     -   ON INPUT, THE GIVEN PERM AND INVERSE PERM
!                           VECTORS.  ON OUTPUT, THE NEW PERM AND
!                           INVERSE PERM VECTORS OF THE NEW
!                           POSTORDERING.
!       COLCNT          -   COLUMN COUNTS IN L UNDER INITIAL ORDERING;
!                           MODIFIED TO REFLECT THE NEW ORDERING.
!
!   OUTPUT PARAMETERS:
!       PARENT          -   THE PARENT VECTOR OF THE ELIMINATION TREE
!                           ASSOCIATED WITH THE NEW ORDERING.
!
!   WORKING PARAMETERS:
!       FSON            -   THE FIRST SON VECTOR.
!       BROTHR          -   THE BROTHER VECTOR.
!       INVPOS          -   THE INVERSE PERM VECTOR FOR THE
!                           POSTORDERING.
!
!   PROGRAM SUBROUTINES:
!       BTREE2, EPOST2, INVINV.
!
!***********************************************************************
!
SUBROUTINE  CHORDR (  NEQNS , PERM  , INVP  , &
& COLCNT, PARENT, FSON  , BROTHR, INVPOS  )
!
!***********************************************************************
!
INTEGER             BROTHR(*)     , &
& COLCNT(*)     , FSON(*)       , &
& INVP(*)       , INVPOS(*)     , &
& PARENT(*)     , PERM(*)
!
INTEGER             NEQNS
!
!***********************************************************************
!
!       ----------------------------------------------------------
!       COMPUTE A BINARY REPRESENTATION OF THE ELIMINATION TREE,
!       SO THAT EACH "LAST CHILD" MAXIMIZES AMONG ITS SIBLINGS THE
!       NUMBER OF NONZEROS IN THE CORRESPONDING COLUMNS OF L.
!       ----------------------------------------------------------
CALL  BTREE2  ( NEQNS , PARENT, COLCNT, FSON  , BROTHR, &
& INVPOS                                  )
!
!       ----------------------------------------------------
!       POSTORDER THE ELIMINATION TREE (USING THE NEW BINARY
!       REPRESENTATION.
!       ----------------------------------------------------
CALL  EPOST2  ( NEQNS , FSON  , BROTHR, INVPOS, PARENT, &
& COLCNT, PERM                            )
!
!       --------------------------------------------------------
!       COMPOSE THE ORIGINAL ORDERING WITH THE NEW POSTORDERING.
!       --------------------------------------------------------
CALL  INVINV  ( NEQNS , INVP  , INVPOS, PERM    )
!
RETURN
END
!***********************************************************************
!***********************************************************************
!
!   Version:        0.4
!   Last modified:  December 27, 1994
!   Authors:        Esmond G. Ng and Barry W. Peyton
!
!   Mathematical Sciences Section, Oak Ridge National Laboratory
!
!***********************************************************************
!***********************************************************************
!******     DSCAL1 .... SCALE A VECTOR                     **************
!***********************************************************************
!***********************************************************************
!
!     PURPOSE - THIS ROUTINE COMPUTES A <-- AX, WHERE A IS A
!               SCALAR AND X IS A VECTOR.
!
!     INPUT PARAMETERS -
!        N - LENGTH OF THE VECTOR X.
!        A - SCALAR MULIPLIER.
!        X - VECTOR TO BE SCALED.
!
!     OUTPUT PARAMETERS -
!        X - REPLACED BY THE SCALED VECTOR, AX.
!
!***********************************************************************
!
SUBROUTINE  DSCAL1 ( N, A, X )
!
!***********************************************************************
!
!     -----------
!     PARAMETERS.
!     -----------
INTEGER             N
DOUBLE PRECISION    A, X(N)
!
!     ----------------
!     LOCAL VARIABLES.
!     ----------------
INTEGER             I
!
!***********************************************************************
!
DO  100  I = 1, N
X(I) = A * X(I)
100 CONTINUE
RETURN
END
!***********************************************************************
!***********************************************************************
!
!   Version:        0.4
!   Last modified:  December 27, 1994
!   Authors:        Esmond G. Ng and Barry W. Peyton
!
!   Mathematical Sciences Section, Oak Ridge National Laboratory
!
!***********************************************************************
!***********************************************************************
!***************     EPOST2 ..... ETREE POSTORDERING #2  ***************
!***********************************************************************
!***********************************************************************
!
!   PURPOSE:
!       BASED ON THE BINARY REPRESENTATION (FIRST-SON,BROTHER) OF THE
!       ELIMINATION TREE, A POSTORDERING IS DETERMINED. THE
!       CORRESPONDING PARENT AND COLCNT VECTORS ARE ALSO MODIFIED TO
!       REFLECT THE REORDERING.
!
!   INPUT PARAMETERS:
!       ROOT            -   ROOT OF THE ELIMINATION TREE (USUALLY IT
!                           IS NEQNS).
!       FSON            -   THE FIRST SON VECTOR.
!       BROTHR          -   THE BROTHR VECTOR.
!
!   UPDATED PARAMETERS:
!       PARENT          -   THE PARENT VECTOR.
!       COLCNT          -   COLUMN NONZERO COUNTS OF THE FACTOR.
!
!   OUTPUT PARAMETERS:
!       INVPOS          -   INVERSE PERMUTATION FOR THE POSTORDERING.
!
!   WORKING PARAMETERS:
!       STACK           -   THE STACK FOR POSTORDER TRAVERSAL OF THE
!                           TREE.
!
!***********************************************************************
!
SUBROUTINE  EPOST2 (  ROOT  , FSON  , BROTHR, INVPOS, PARENT, &
& COLCNT, STACK                           )
!
!***********************************************************************
!
INTEGER           BROTHR(*)     , COLCNT(*)     , &
& FSON(*)       , INVPOS(*)     , &
& PARENT(*)     , STACK(*)
!
INTEGER           ROOT
!
!***********************************************************************
!
INTEGER           ITOP  , NDPAR , NODE  , NUM   , NUNODE
!
!***********************************************************************
!
NUM = 0
ITOP = 0
NODE = ROOT
!       -------------------------------------------------------------
!       TRAVERSE ALONG THE FIRST SONS POINTER AND PUSH THE TREE NODES
!       ALONG THE TRAVERSAL INTO THE STACK.
!       -------------------------------------------------------------
100 CONTINUE
ITOP = ITOP + 1
STACK(ITOP) = NODE
NODE = FSON(NODE)
IF  ( NODE .GT. 0 )  GO TO 100
!           ----------------------------------------------------------
!           IF POSSIBLE, POP A TREE NODE FROM THE STACK AND NUMBER IT.
!           ----------------------------------------------------------
200 CONTINUE
IF  ( ITOP .LE. 0 )  GO TO 300
NODE = STACK(ITOP)
ITOP = ITOP - 1
NUM = NUM + 1
INVPOS(NODE) = NUM
!               ----------------------------------------------------
!               THEN, TRAVERSE TO ITS YOUNGER BROTHER IF IT HAS ONE.
!               ----------------------------------------------------
NODE = BROTHR(NODE)
IF  ( NODE .LE. 0 )  GO TO 200
GO TO 100
!
300 CONTINUE
!       ------------------------------------------------------------
!       DETERMINE THE NEW PARENT VECTOR OF THE POSTORDERING.  BROTHR
!       IS USED TEMPORARILY FOR THE NEW PARENT VECTOR.
!       ------------------------------------------------------------
DO  400  NODE = 1, NUM
NUNODE = INVPOS(NODE)
NDPAR = PARENT(NODE)
IF  ( NDPAR .GT. 0 )  NDPAR = INVPOS(NDPAR)
BROTHR(NUNODE) = NDPAR
400 CONTINUE
!
DO  500  NUNODE = 1, NUM
PARENT(NUNODE) = BROTHR(NUNODE)
500 CONTINUE
!
!       ----------------------------------------------
!       PERMUTE COLCNT(*) TO REFLECT THE NEW ORDERING.
!       ----------------------------------------------
DO  600  NODE = 1, NUM
NUNODE = INVPOS(NODE)
STACK(NUNODE) = COLCNT(NODE)
600 CONTINUE
!
DO  700  NODE = 1, NUM
COLCNT(NODE) = STACK(NODE)
700 CONTINUE
!
RETURN
END
!***********************************************************************
!***********************************************************************
!
!   Version:        0.4
!   Last modified:  December 27, 1994
!   Authors:        Joseph W.H. Liu
!
!   Mathematical Sciences Section, Oak Ridge National Laboratory
!
!***********************************************************************
!***********************************************************************
!**********     ETORDR ..... ELIMINATION TREE REORDERING     ***********
!***********************************************************************
!***********************************************************************
!
!   WRITTEN BY JOSEPH LIU (JUL 17, 1985)
!
!   PURPOSE:
!       TO DETERMINE AN EQUIVALENT REORDERING BASED ON THE STRUCTURE OF
!       THE ELIMINATION TREE.  A POSTORDERING OF THE GIVEN ELIMINATION
!       TREE IS RETURNED.
!
!   INPUT PARAMETERS:
!       NEQNS           -   NUMBER OF EQUATIONS.
!       (XADJ,ADJNCY)   -   THE ADJACENCY STRUCTURE.
!
!   UPDATED PARAMETERS:
!       (PERM,INVP)     -   ON INPUT, THE GIVEN PERM AND INVERSE PERM
!                           VECTORS.  ON OUTPUT, THE NEW PERM AND
!                           INVERSE PERM VECTORS OF THE EQUIVALENT
!                           ORDERING.
!
!   OUTPUT PARAMETERS:
!       PARENT          -   THE PARENT VECTOR OF THE ELIMINATION TREE
!                           ASSOCIATED WITH THE NEW ORDERING.
!
!   WORKING PARAMETERS:
!       FSON            -   THE FIRST SON VECTOR.
!       BROTHR          -   THE BROTHER VECTOR.
!       INVPOS          -   THE INVERSE PERM VECTOR FOR THE
!                           POSTORDERING.
!
!   PROGRAM SUBROUTINES:
!       BETREE, ETPOST, ETREE , INVINV.
!
!***********************************************************************
!
SUBROUTINE  ETORDR (  NEQNS , XADJ  , ADJNCY, PERM  , INVP  , &
& PARENT, FSON  , BROTHR, INVPOS          )
!
!***********************************************************************
!
INTEGER           ADJNCY(*)     , BROTHR(*)     , &
& FSON(*)       , INVP(*)       , &
& INVPOS(*)     , PARENT(*)     , &
& PERM(*)
!
INTEGER           XADJ(*)
INTEGER           NEQNS
!
!***********************************************************************
!
!       -----------------------------
!       COMPUTE THE ELIMINATION TREE.
!       -----------------------------
CALL  ETREE ( NEQNS, XADJ, ADJNCY, PERM, INVP, PARENT, INVPOS )
!
!       --------------------------------------------------------
!       COMPUTE A BINARY REPRESENTATION OF THE ELIMINATION TREE.
!       --------------------------------------------------------
CALL  BETREE ( NEQNS, PARENT, FSON, BROTHR )
!
!       -------------------------------
!       POSTORDER THE ELIMINATION TREE.
!       -------------------------------
CALL  ETPOST ( NEQNS, FSON, BROTHR, INVPOS, PARENT, PERM )
!
!       --------------------------------------------------------
!       COMPOSE THE ORIGINAL ORDERING WITH THE NEW POSTORDERING.
!       --------------------------------------------------------
CALL  INVINV ( NEQNS, INVP, INVPOS, PERM )
!
RETURN
END
!***********************************************************************
!***********************************************************************
!
!   Version:        0.4
!   Last modified:  December 27, 1994
!   Authors:        Joseph W.H. Liu
!
!   Mathematical Sciences Section, Oak Ridge National Laboratory
!
!***********************************************************************
!***********************************************************************
!***************     ETPOST ..... ETREE POSTORDERING     ***************
!***********************************************************************
!***********************************************************************
!
!   WRITTEN BY JOSEPH LIU (SEPT 17, 1986)
!
!   PURPOSE:
!       BASED ON THE BINARY REPRESENTATION (FIRST-SON,BROTHER) OF
!       THE ELIMINATION TREE, A POSTORDERING IS DETERMINED. THE
!       CORRESPONDING PARENT VECTOR IS ALSO MODIFIED TO REFLECT
!       THE REORDERING.
!
!   INPUT PARAMETERS:
!       ROOT            -   ROOT OF THE ELIMINATION TREE (USUALLY IT
!                           IS NEQNS).
!       FSON            -   THE FIRST SON VECTOR.
!       BROTHR          -   THE BROTHR VECTOR.
!
!   UPDATED PARAMETERS:
!       PARENT          -   THE PARENT VECTOR.
!
!   OUTPUT PARAMETERS:
!       INVPOS          -   INVERSE PERMUTATION FOR THE POSTORDERING.
!
!   WORKING PARAMETERS:
!       STACK           -   THE STACK FOR POSTORDER TRAVERSAL OF THE
!                           TREE.
!
!***********************************************************************
!
SUBROUTINE  ETPOST (  ROOT  , FSON  , BROTHR, INVPOS, PARENT, &
& STACK                                   )
!
!***********************************************************************
!
INTEGER           BROTHR(*)     , FSON(*)       , &
& INVPOS(*)     , PARENT(*)     , &
& STACK(*)
!
INTEGER           ROOT
!
!***********************************************************************
!
INTEGER           ITOP  , NDPAR , NODE  , NUM   , NUNODE
!
!***********************************************************************
!
NUM = 0
ITOP = 0
NODE = ROOT
!       -------------------------------------------------------------
!       TRAVERSE ALONG THE FIRST SONS POINTER AND PUSH THE TREE NODES
!       ALONG THE TRAVERSAL INTO THE STACK.
!       -------------------------------------------------------------
100 CONTINUE
ITOP = ITOP + 1
STACK(ITOP) = NODE
NODE = FSON(NODE)
IF  ( NODE .GT. 0 )  GO TO 100
!           ----------------------------------------------------------
!           IF POSSIBLE, POP A TREE NODE FROM THE STACK AND NUMBER IT.
!           ----------------------------------------------------------
200 CONTINUE
IF  ( ITOP .LE. 0 )  GO TO 300
NODE = STACK(ITOP)
ITOP = ITOP - 1
NUM = NUM + 1
INVPOS(NODE) = NUM
!               ----------------------------------------------------
!               THEN, TRAVERSE TO ITS YOUNGER BROTHER IF IT HAS ONE.
!               ----------------------------------------------------
NODE = BROTHR(NODE)
IF  ( NODE .LE. 0 )  GO TO 200
GO TO 100
!
300 CONTINUE
!       ------------------------------------------------------------
!       DETERMINE THE NEW PARENT VECTOR OF THE POSTORDERING.  BROTHR
!       IS USED TEMPORARILY FOR THE NEW PARENT VECTOR.
!       ------------------------------------------------------------
DO  400  NODE = 1, NUM
NUNODE = INVPOS(NODE)
NDPAR = PARENT(NODE)
IF  ( NDPAR .GT. 0 )  NDPAR = INVPOS(NDPAR)
BROTHR(NUNODE) = NDPAR
400 CONTINUE
!
DO  500  NUNODE = 1, NUM
PARENT(NUNODE) = BROTHR(NUNODE)
500 CONTINUE
!
RETURN
END
!***********************************************************************
!***********************************************************************
!
!   Version:        0.4
!   Last modified:  December 27, 1994
!   Authors:        Joseph W.H. Liu
!
!   Mathematical Sciences Section, Oak Ridge National Laboratory
!
!***********************************************************************
!***********************************************************************
!****************     ETREE ..... ELIMINATION TREE     *****************
!***********************************************************************
!***********************************************************************
!
!   WRITTEN BY JOSEPH LIU (JUL 17, 1985)
!
!   PURPOSE:
!       TO DETERMINE THE ELIMINATION TREE FROM A GIVEN ORDERING AND
!       THE ADJACENCY STRUCTURE.  THE PARENT VECTOR IS RETURNED.
!
!   INPUT PARAMETERS:
!       NEQNS           -   NUMBER OF EQUATIONS.
!       (XADJ,ADJNCY)   -   THE ADJACENCY STRUCTURE.
!       (PERM,INVP)     -   PERMUTATION AND INVERSE PERMUTATION VECTORS
!
!   OUTPUT PARAMETERS:
!       PARENT          -   THE PARENT VECTOR OF THE ELIMINATION TREE.
!
!   WORKING PARAMETERS:
!       ANCSTR          -   THE ANCESTOR VECTOR.
!
!***********************************************************************
!
SUBROUTINE  ETREE (   NEQNS , XADJ  , ADJNCY, PERM  , INVP  , &
& PARENT, ANCSTR                          )
!
!***********************************************************************
!
INTEGER           ADJNCY(*)     , ANCSTR(*)     , &
& INVP(*)       , PARENT(*)     , &
& PERM(*)
!
INTEGER           NEQNS
INTEGER           XADJ(*)
!
!***********************************************************************
!
INTEGER           I     , J     , JSTOP , JSTRT , NBR   , &
& NEXT  , NODE
!
!***********************************************************************
!
IF  ( NEQNS .LE. 0 )  RETURN
!
DO  400  I = 1, NEQNS
PARENT(I) = 0
ANCSTR(I) = 0
NODE = PERM(I)
!
JSTRT = XADJ(NODE)
JSTOP = XADJ(NODE+1) - 1
IF  ( JSTRT .LE. JSTOP )  THEN
DO  300  J = JSTRT, JSTOP
NBR = ADJNCY(J)
NBR = INVP(NBR)
IF  ( NBR .LT. I )  THEN
!                       -------------------------------------------
!                       FOR EACH NBR, FIND THE ROOT OF ITS CURRENT
!                       ELIMINATION TREE.  PERFORM PATH COMPRESSION
!                       AS THE SUBTREE IS TRAVERSED.
!                       -------------------------------------------
100 CONTINUE
IF  ( ANCSTR(NBR) .EQ. I )  GO TO 300
IF  ( ANCSTR(NBR) .GT. 0 )  THEN
NEXT = ANCSTR(NBR)
ANCSTR(NBR) = I
NBR = NEXT
GO TO 100
ENDIF
!                       --------------------------------------------
!                       NOW, NBR IS THE ROOT OF THE SUBTREE.  MAKE I
!                       THE PARENT NODE OF THIS ROOT.
!                       --------------------------------------------
PARENT(NBR) = I
ANCSTR(NBR) = I
ENDIF
300 CONTINUE
ENDIF
400 CONTINUE
!
RETURN
END
!***********************************************************************
!***********************************************************************
!
!   Version:        0.4
!   Last modified:  January 12, 1995
!   Authors:        Esmond G. Ng and Barry W. Peyton
!
!   Mathematical Sciences Section, Oak Ridge National Laboratory
!
!***********************************************************************
!***********************************************************************
!**************     FCNTHN  ..... FIND NONZERO COUNTS    ***************
!***********************************************************************
!***********************************************************************
!
!   PURPOSE:
!       THIS SUBROUTINE DETERMINES THE ROW COUNTS AND COLUMN COUNTS IN
!       THE CHOLESKY FACTOR.  IT USES A DISJOINT SET UNION ALGORITHM.
!
!       TECHNIQUES:
!       1) SUPERNODE DETECTION.
!       2) PATH HALVING.
!       3) NO UNION BY RANK.
!
!   NOTES:
!       1) ASSUMES A POSTORDERING OF THE ELIMINATION TREE.
!
!   INPUT PARAMETERS:
!       (I) NEQNS       -   NUMBER OF EQUATIONS.
!       (I) ADJLEN      -   LENGTH OF ADJACENCY STRUCTURE.
!       (I) XADJ(*)     -   ARRAY OF LENGTH NEQNS+1, CONTAINING POINTERS
!                           TO THE ADJACENCY STRUCTURE.
!       (I) ADJNCY(*)   -   ARRAY OF LENGTH XADJ(NEQNS+1)-1, CONTAINING
!                           THE ADJACENCY STRUCTURE.
!       (I) PERM(*)     -   ARRAY OF LENGTH NEQNS, CONTAINING THE
!                           POSTORDERING.
!       (I) INVP(*)     -   ARRAY OF LENGTH NEQNS, CONTAINING THE
!                           INVERSE OF THE POSTORDERING.
!       (I) ETPAR(*)    -   ARRAY OF LENGTH NEQNS, CONTAINING THE
!                           ELIMINATION TREE OF THE POSTORDERED MATRIX.
!
!   OUTPUT PARAMETERS:
!       (I) ROWCNT(*)   -   ARRAY OF LENGTH NEQNS, CONTAINING THE NUMBER
!                           OF NONZEROS IN EACH ROW OF THE FACTOR,
!                           INCLUDING THE DIAGONAL ENTRY.
!       (I) COLCNT(*)   -   ARRAY OF LENGTH NEQNS, CONTAINING THE NUMBER
!                           OF NONZEROS IN EACH COLUMN OF THE FACTOR,
!                           INCLUDING THE DIAGONAL ENTRY.
!       (I) NLNZ        -   NUMBER OF NONZEROS IN THE FACTOR, INCLUDING
!                           THE DIAGONAL ENTRIES.
!
!   WORK PARAMETERS:
!       (I) SET(*)      -   ARRAY OF LENGTH NEQNS USED TO MAINTAIN THE
!                           DISJOINT SETS (I.E., SUBTREES).
!       (I) PRVLF(*)    -   ARRAY OF LENGTH NEQNS USED TO RECORD THE
!                           PREVIOUS LEAF OF EACH ROW SUBTREE.
!       (I) LEVEL(*)    -   ARRAY OF LENGTH NEQNS+1 CONTAINING THE LEVEL
!                           (DISTANCE FROM THE ROOT).
!       (I) WEIGHT(*)   -   ARRAY OF LENGTH NEQNS+1 CONTAINING WEIGHTS
!                           USED TO COMPUTE COLUMN COUNTS.
!       (I) FDESC(*)    -   ARRAY OF LENGTH NEQNS+1 CONTAINING THE
!                           FIRST (I.E., LOWEST-NUMBERED) DESCENDANT.
!       (I) NCHILD(*)   -   ARRAY OF LENGTH NEQNS+1 CONTAINING THE
!                           NUMBER OF CHILDREN.
!       (I) PRVNBR(*)   -   ARRAY OF LENGTH NEQNS USED TO RECORD THE
!                           PREVIOUS ``LOWER NEIGHBOR'' OF EACH NODE.
!
!   FIRST CREATED ON    APRIL 12, 1990.
!   LAST UPDATED ON     JANUARY 12, 1995.
!
!***********************************************************************
!
SUBROUTINE FCNTHN  (  NEQNS , ADJLEN, XADJ  , ADJNCY, PERM  , &
& INVP  , ETPAR , ROWCNT, COLCNT, NLNZ  , &
& SET   , PRVLF , LEVEL , WEIGHT, FDESC , &
& NCHILD, PRVNBR                          )
!
!       -----------
!       PARAMETERS.
!       -----------
INTEGER             ADJLEN, NEQNS , NLNZ
INTEGER             ADJNCY(ADJLEN)  , COLCNT(NEQNS) , &
& ETPAR(NEQNS)    , FDESC(0:NEQNS), &
& INVP(NEQNS)     , LEVEL(0:NEQNS), &
& NCHILD(0:NEQNS) , PERM(NEQNS)   , &
& PRVLF(NEQNS)    , PRVNBR(NEQNS) , &
& ROWCNT(NEQNS)   , SET(NEQNS)    , &
& WEIGHT(0:NEQNS)
INTEGER             XADJ(*)
!
!       ----------------
!       LOCAL VARIABLES.
!       ----------------
INTEGER             HINBR , IFDESC, J     , JSTOP , JSTRT , &
& K     , LAST1 , LAST2 , LCA   , LFLAG , &
& LOWNBR, OLDNBR, PARENT, PLEAF , TEMP  , &
& XSUP
!
!***********************************************************************
!
!       --------------------------------------------------
!       COMPUTE LEVEL(*), FDESC(*), NCHILD(*).
!       INITIALIZE XSUP, ROWCNT(*), COLCNT(*),
!                  SET(*), PRVLF(*), WEIGHT(*), PRVNBR(*).
!       --------------------------------------------------
XSUP = 1
LEVEL(0) = 0
DO  100  K = NEQNS, 1, -1
ROWCNT(K) = 1
COLCNT(K) = 0
SET(K) = K
PRVLF(K) = 0
LEVEL(K) = LEVEL(ETPAR(K)) + 1
WEIGHT(K) = 1
FDESC(K) = K
NCHILD(K) = 0
PRVNBR(K) = 0
100 CONTINUE
NCHILD(0) = 0
FDESC(0) = 0
DO  200  K = 1, NEQNS
PARENT = ETPAR(K)
WEIGHT(PARENT) = 0
NCHILD(PARENT) = NCHILD(PARENT) + 1
IFDESC = FDESC(K)
IF  ( IFDESC .LT. FDESC(PARENT) )  THEN
FDESC(PARENT) = IFDESC
ENDIF
200 CONTINUE
!       ------------------------------------
!       FOR EACH ``LOW NEIGHBOR'' LOWNBR ...
!       ------------------------------------
DO  600  LOWNBR = 1, NEQNS
LFLAG = 0
IFDESC = FDESC(LOWNBR)
OLDNBR = PERM(LOWNBR)
JSTRT = XADJ(OLDNBR)
JSTOP = XADJ(OLDNBR+1) - 1
!           -----------------------------------------------
!           FOR EACH ``HIGH NEIGHBOR'', HINBR OF LOWNBR ...
!           -----------------------------------------------
DO  500  J = JSTRT, JSTOP
HINBR = INVP(ADJNCY(J))
IF  ( HINBR .GT. LOWNBR )  THEN
IF  ( IFDESC .GT. PRVNBR(HINBR) )  THEN
!                       -------------------------
!                       INCREMENT WEIGHT(LOWNBR).
!                       -------------------------
WEIGHT(LOWNBR) = WEIGHT(LOWNBR) + 1
PLEAF = PRVLF(HINBR)
!                       -----------------------------------------
!                       IF HINBR HAS NO PREVIOUS ``LOW NEIGHBOR''
!                       THEN ...
!                       -----------------------------------------
IF  ( PLEAF .EQ. 0 )  THEN
!                           -----------------------------------------
!                           ... ACCUMULATE LOWNBR-->HINBR PATH LENGTH
!                               IN ROWCNT(HINBR).
!                           -----------------------------------------
ROWCNT(HINBR) = ROWCNT(HINBR) + &
& LEVEL(LOWNBR) - LEVEL(HINBR)
ELSE
!                           -----------------------------------------
!                           ... OTHERWISE, LCA <-- FIND(PLEAF), WHICH
!                               IS THE LEAST COMMON ANCESTOR OF PLEAF
!                               AND LOWNBR.
!                               (PATH HALVING.)
!                           -----------------------------------------
LAST1 = PLEAF
LAST2 = SET(LAST1)
LCA = SET(LAST2)
300 CONTINUE
IF  ( LCA .NE. LAST2 )  THEN
SET(LAST1) = LCA
LAST1 = LCA
LAST2 = SET(LAST1)
LCA = SET(LAST2)
GO TO 300
ENDIF
!                           -------------------------------------
!                           ACCUMULATE PLEAF-->LCA PATH LENGTH IN
!                           ROWCNT(HINBR).
!                           DECREMENT WEIGHT(LCA).
!                           -------------------------------------
ROWCNT(HINBR) = ROWCNT(HINBR) &
& + LEVEL(LOWNBR) - LEVEL(LCA)
WEIGHT(LCA) = WEIGHT(LCA) - 1
ENDIF
!                       ----------------------------------------------
!                       LOWNBR NOW BECOMES ``PREVIOUS LEAF'' OF HINBR.
!                       ----------------------------------------------
PRVLF(HINBR) = LOWNBR
LFLAG = 1
ENDIF
!                   --------------------------------------------------
!                   LOWNBR NOW BECOMES ``PREVIOUS NEIGHBOR'' OF HINBR.
!                   --------------------------------------------------
PRVNBR(HINBR) = LOWNBR
ENDIF
500 CONTINUE
!           ----------------------------------------------------
!           DECREMENT WEIGHT ( PARENT(LOWNBR) ).
!           SET ( P(LOWNBR) ) <-- SET ( P(LOWNBR) ) + SET(XSUP).
!           ----------------------------------------------------
PARENT = ETPAR(LOWNBR)
WEIGHT(PARENT) = WEIGHT(PARENT) - 1
IF  ( LFLAG .EQ. 1     .OR. &
& NCHILD(LOWNBR) .GE. 2 )  THEN
XSUP = LOWNBR
ENDIF
SET(XSUP) = PARENT
600 CONTINUE
!       ---------------------------------------------------------
!       USE WEIGHTS TO COMPUTE COLUMN (AND TOTAL) NONZERO COUNTS.
!       ---------------------------------------------------------
NLNZ = 0
DO  700  K = 1, NEQNS
TEMP = COLCNT(K) + WEIGHT(K)
COLCNT(K) = TEMP
NLNZ = NLNZ + TEMP
PARENT = ETPAR(K)
IF  ( PARENT .NE. 0 )  THEN
COLCNT(PARENT) = COLCNT(PARENT) + TEMP
ENDIF
700 CONTINUE
!
RETURN
END
!***********************************************************************
!***********************************************************************
!
!   Version:        0.4
!   Last modified:  May 26, 1995
!   Authors:        Esmond G. Ng and Barry W. Peyton
!
!   Mathematical Sciences Section, Oak Ridge National Laboratory
!
!***********************************************************************
!***********************************************************************
!****     FNSPLT ..... COMPUTE FINE PARTITIONING OF SUPERNODES     *****
!***********************************************************************
!***********************************************************************
!
!   PURPOSE:
!       THIS SUBROUTINE DETERMINES A FINE PARTITIONING OF SUPERNODES
!       WHEN THERE IS A CACHE AVAILABLE ON THE MACHINE.  THE FINE
!       PARTITIONING IS CHOSEN SO THAT DATA RE-USE IS MAXIMIZED.
!
!   INPUT PARAMETERS:
!       NEQNS           -   NUMBER OF EQUATIONS.
!       NSUPER          -   NUMBER OF SUPERNODES.
!       XSUPER          -   INTEGER ARRAY OF SIZE (NSUPER+1) CONTAINING
!                           THE SUPERNODE PARTITIONING.
!       XLINDX          -   INTEGER ARRAY OF SIZE (NSUPER+1) CONTAINING
!                           POINTERS IN THE SUPERNODE INDICES.
!       CACHSZ          -   CACHE SIZE IN KILO BYTES.
!                           IF THERE IS NO CACHE, SET CACHSZ = 0.
!
!   OUTPUT PARAMETERS:
!       SPLIT           -   INTEGER ARRAY OF SIZE NEQNS CONTAINING THE
!                           FINE PARTITIONING.
!
!***********************************************************************
!
SUBROUTINE  FNSPLT ( NEQNS , NSUPER, XSUPER, XLINDX, &
& CACHSZ, SPLIT )
!
!***********************************************************************
!
!       -----------
!       PARAMETERS.
!       -----------
INTEGER         CACHSZ, NEQNS , NSUPER
INTEGER         XSUPER(*), SPLIT(*)
INTEGER         XLINDX(*)
!
!       ----------------
!       LOCAL VARIABLES.
!       ----------------
INTEGER         CACHE , CURCOL, FSTCOL, HEIGHT, KCOL  , &
& KSUP  , LSTCOL, NCOLS , NXTBLK, USED  , &
& WIDTH
!
! *******************************************************************
!
!       --------------------------------------------
!       COMPUTE THE NUMBER OF 8-BYTE WORDS IN CACHE.
!       --------------------------------------------
IF  ( CACHSZ .LE. 0 )  THEN
CACHE = 2000000000
ELSE
CACHE = CACHSZ * 116!INT( ( FLOAT(CACHSZ) * 1024. / 8. ) * 0.9 )
ENDIF
!
!       ---------------
!       INITIALIZATION.
!       ---------------
DO  100  KCOL = 1, NEQNS
SPLIT(KCOL) = 0
100 CONTINUE
!
!       ---------------------------
!       FOR EACH SUPERNODE KSUP ...
!       ---------------------------
DO  1000  KSUP = 1, NSUPER
!           -----------------------
!           ... GET SUPERNODE INFO.
!           -----------------------
HEIGHT = XLINDX(KSUP+1) - XLINDX(KSUP)
FSTCOL = XSUPER(KSUP)
LSTCOL = XSUPER(KSUP+1) - 1
WIDTH = LSTCOL - FSTCOL + 1
NXTBLK = FSTCOL
!           --------------------------------------
!           ... UNTIL ALL COLUMNS OF THE SUPERNODE
!               HAVE BEEN PROCESSED ...
!           --------------------------------------
CURCOL = FSTCOL - 1
200 CONTINUE
!               -------------------------------------------
!               ... PLACE THE FIRST COLUMN(S) IN THE CACHE.
!               -------------------------------------------
CURCOL = CURCOL + 1
IF  ( CURCOL .LT. LSTCOL )  THEN
CURCOL = CURCOL + 1
NCOLS = 2
USED = 4 * HEIGHT - 1
HEIGHT = HEIGHT - 2
ELSE
NCOLS = 1
USED = 3 * HEIGHT
HEIGHT = HEIGHT - 1
ENDIF
!
!               --------------------------------------
!               ... WHILE THE CACHE IS NOT FILLED AND
!                   THERE ARE COLUMNS OF THE SUPERNODE
!                   REMAINING TO BE PROCESSED ...
!               --------------------------------------
300 CONTINUE
IF  ( USED+HEIGHT .LT. CACHE  .AND. &
& CURCOL      .LT. LSTCOL       )  THEN
!                   --------------------------------
!                   ... ADD ANOTHER COLUMN TO CACHE.
!                   --------------------------------
CURCOL = CURCOL + 1
NCOLS = NCOLS + 1
USED = USED + HEIGHT
HEIGHT = HEIGHT - 1
GO TO 300
ENDIF
!               -------------------------------------
!               ... RECORD THE NUMBER OF COLUMNS THAT
!                   FILLED THE CACHE.
!               -------------------------------------
SPLIT(NXTBLK) = NCOLS
NXTBLK = NXTBLK + 1
!               --------------------------
!               ... GO PROCESS NEXT BLOCK.
!               --------------------------
IF  ( CURCOL .LT. LSTCOL )  GO TO 200
1000 CONTINUE
!
RETURN
END
!***********************************************************************
!***********************************************************************
!
!   Version:        0.4
!   Last modified:  December 27, 1994
!   Authors:        Esmond G. Ng and Barry W. Peyton
!
!   Mathematical Sciences Section, Oak Ridge National Laboratory
!
!***********************************************************************
!***********************************************************************
!******     FNTSIZ ..... COMPUTE WORK STORAGE SIZE FOR BLKFCT     ******
!***********************************************************************
!***********************************************************************
!
!   PURPOSE:
!       THIS SUBROUTINE DETERMINES THE SIZE OF THE WORKING STORAGE
!       REQUIRED BY BLKFCT.
!
!   INPUT PARAMETERS:
!       NSUPER          -   NUMBER OF SUPERNODES.
!       XSUPER          -   INTEGER ARRAY OF SIZE (NSUPER+1) CONTAINING
!                           THE SUPERNODE PARTITIONING.
!       SNODE           -   SUPERNODE MEMBERSHIP.
!       (XLINDX,LINDX)  -   ARRAYS DESCRIBING THE SUPERNODAL STRUCTURE.
!
!   OUTPUT PARAMETERS:
!       TMPSIZ          -   SIZE OF WORKING STORAGE REQUIRED BY BLKFCT.
!
!***********************************************************************
!
SUBROUTINE  FNTSIZ ( NSUPER, XSUPER, SNODE , XLINDX, &
& LINDX , TMPSIZ  )
!
!***********************************************************************
!
INTEGER     NSUPER, TMPSIZ
INTEGER     XLINDX(*)       , XSUPER(*)
INTEGER     LINDX (*)       , SNODE (*)
!
INTEGER     BOUND , CLEN  , CURSUP, I     , IBEGIN, IEND  , &
& KSUP  , LENGTH, NCOLS , NXTSUP, &
& TSIZE , WIDTH
!
!***********************************************************************
!
!       RETURNS SIZE OF TEMP ARRAY USED BY BLKFCT FACTORIZATION ROUTINE.
!       NOTE THAT THE VALUE RETURNED IS AN ESTIMATE, THOUGH IT IS USUALLY
!       TIGHT.
!
!       ----------------------------------------
!       COMPUTE SIZE OF TEMPORARY STORAGE VECTOR
!       NEEDED BY BLKFCT.
!       ----------------------------------------
TMPSIZ = 0
DO  500  KSUP = NSUPER, 1, -1
NCOLS = XSUPER(KSUP+1) - XSUPER(KSUP)
IBEGIN = XLINDX(KSUP) + NCOLS
IEND = XLINDX(KSUP+1) - 1
LENGTH = IEND - IBEGIN + 1
BOUND = LENGTH * (LENGTH + 1) / 2
IF  ( BOUND .GT. TMPSIZ )  THEN
CURSUP = SNODE(LINDX(IBEGIN))
CLEN = XLINDX(CURSUP+1) - XLINDX(CURSUP)
WIDTH = 0
DO  400  I = IBEGIN, IEND
NXTSUP = SNODE(LINDX(I))
IF  ( NXTSUP .EQ. CURSUP )  THEN
WIDTH = WIDTH + 1
IF  ( I .EQ. IEND )  THEN
IF  ( CLEN .GT. LENGTH )  THEN
TSIZE = LENGTH * WIDTH - &
& (WIDTH - 1) * WIDTH / 2
TMPSIZ = MAX ( TSIZE , TMPSIZ )
ENDIF
ENDIF
ELSE
IF  ( CLEN .GT. LENGTH )  THEN
TSIZE = LENGTH * WIDTH - &
& (WIDTH - 1) * WIDTH / 2
TMPSIZ = MAX ( TSIZE , TMPSIZ )
ENDIF
LENGTH = LENGTH - WIDTH
BOUND = LENGTH * (LENGTH + 1) / 2
IF  ( BOUND .LE. TMPSIZ )  GO TO 500
WIDTH = 1
CURSUP = NXTSUP
CLEN = XLINDX(CURSUP+1) - XLINDX(CURSUP)
ENDIF
400 CONTINUE
ENDIF
500 CONTINUE
!
RETURN
END
!***********************************************************************
!***********************************************************************
!
!   Version:        0.4
!   Last modified:  December 27, 1994
!   Authors:        Esmond G. Ng and Barry W. Peyton
!
!   Mathematical Sciences Section, Oak Ridge National Laboratory
!
!***********************************************************************
!***********************************************************************
!****************    FSUP1 ..... FIND SUPERNODES #1    *****************
!***********************************************************************
!***********************************************************************
!
!   PURPOSE:
!       THIS SUBROUTINE IS THE FIRST OF TWO ROUTINES FOR FINDING A
!       MAXIMAL SUPERNODE PARTITION.  IT RETURNS ONLY THE NUMBER OF
!       SUPERNODES NSUPER AND THE SUPERNODE MEMBERSHIP VECTOR SNODE(*),
!       WHICH IS OF LENGTH NEQNS.  THE VECTORS OF LENGTH NSUPER ARE
!       COMPUTED SUBSEQUENTLY BY THE COMPANION ROUTINE FSUP2.
!
!   METHOD AND ASSUMPTIONS:
!       THIS ROUTINE USES THE ELIMINATION TREE AND THE FACTOR COLUMN
!       COUNTS TO COMPUTE THE SUPERNODE PARTITION; IT ALSO ASSUMES A
!       POSTORDERING OF THE ELIMINATION TREE.
!
!   INPUT PARAMETERS:
!       (I) NEQNS       -   NUMBER OF EQUATIONS.
!       (I) ETPAR(*)    -   ARRAY OF LENGTH NEQNS, CONTAINING THE
!                           ELIMINATION TREE OF THE POSTORDERED MATRIX.
!       (I) COLCNT(*)   -   ARRAY OF LENGTH NEQNS, CONTAINING THE
!                           FACTOR COLUMN COUNTS: I.E., THE NUMBER OF
!                           NONZERO ENTRIES IN EACH COLUMN OF L
!                           (INCLUDING THE DIAGONAL ENTRY).
!
!   OUTPUT PARAMETERS:
!       (I) NOFSUB      -   NUMBER OF SUBSCRIPTS.
!       (I) NSUPER      -   NUMBER OF SUPERNODES (<= NEQNS).
!       (I) SNODE(*)    -   ARRAY OF LENGTH NEQNS FOR RECORDING
!                           SUPERNODE MEMBERSHIP.
!
!   FIRST CREATED ON    JANUARY 18, 1992.
!   LAST UPDATED ON     NOVEMBER 11, 1994.
!
!***********************************************************************
!
SUBROUTINE  FSUP1  (  NEQNS , ETPAR , COLCNT, NOFSUB, NSUPER, &
& SNODE                                   )
!
!***********************************************************************
!
!       -----------
!       PARAMETERS.
!       -----------
INTEGER             NEQNS , NOFSUB, NSUPER
INTEGER             COLCNT(*)     , ETPAR(*)      , &
& SNODE(*)
!
!       ----------------
!       LOCAL VARIABLES.
!       ----------------
INTEGER             KCOL
!
!***********************************************************************
!
!       --------------------------------------------
!       COMPUTE THE FUNDAMENTAL SUPERNODE PARTITION.
!       --------------------------------------------
NSUPER = 1
SNODE(1) = 1
NOFSUB = COLCNT(1)
DO  300  KCOL = 2, NEQNS
IF  ( ETPAR(KCOL-1) .EQ. KCOL )  THEN
IF  ( COLCNT(KCOL-1) .EQ. COLCNT(KCOL)+1 )  THEN
SNODE(KCOL) = NSUPER
GO TO 300
ENDIF
ENDIF
NSUPER = NSUPER + 1
SNODE(KCOL) = NSUPER
NOFSUB = NOFSUB + COLCNT(KCOL)
300 CONTINUE
!
RETURN
END
!***********************************************************************
!***********************************************************************
!
!   Version:        0.4
!   Last modified:  December 27, 1994
!   Authors:        Esmond G. Ng and Barry W. Peyton
!
!   Mathematical Sciences Section, Oak Ridge National Laboratory
!
!***********************************************************************
!***********************************************************************
!****************    FSUP2  ..... FIND SUPERNODES #2   *****************
!***********************************************************************
!***********************************************************************
!
!   PURPOSE:
!       THIS SUBROUTINE IS THE SECOND OF TWO ROUTINES FOR FINDING A
!       MAXIMAL SUPERNODE PARTITION.  IT'S SOLE PURPOSE IS TO
!       CONSTRUCT THE NEEDED VECTOR OF LENGTH NSUPER: XSUPER(*).  THE
!       FIRST ROUTINE FSUP1 COMPUTES THE NUMBER OF SUPERNODES AND THE
!       SUPERNODE MEMBERSHIP VECTOR SNODE(*), WHICH IS OF LENGTH NEQNS.
!
!
!   ASSUMPTIONS:
!       THIS ROUTINE ASSUMES A POSTORDERING OF THE ELIMINATION TREE.  IT
!       ALSO ASSUMES THAT THE OUTPUT FROM FSUP1 IS AVAILABLE.
!
!   INPUT PARAMETERS:
!       (I) NEQNS       -   NUMBER OF EQUATIONS.
!       (I) NSUPER      -   NUMBER OF SUPERNODES (<= NEQNS).
!       (I) SNODE(*)    -   ARRAY OF LENGTH NEQNS FOR RECORDING
!                           SUPERNODE MEMBERSHIP.
!
!   OUTPUT PARAMETERS:
!       (I) XSUPER(*)   -   ARRAY OF LENGTH NSUPER+1, CONTAINING THE
!                           SUPERNODE PARTITIONING.
!
!   FIRST CREATED ON    JANUARY 18, 1992.
!   LAST UPDATED ON     NOVEMEBER 22, 1994.
!
!***********************************************************************
!
SUBROUTINE  FSUP2  (  NEQNS , NSUPER, SNODE , XSUPER  )
!
!***********************************************************************
!
!       -----------
!       PARAMETERS.
!       -----------
INTEGER             NEQNS , NSUPER
INTEGER             SNODE(*)      , &
& XSUPER(*)
!
!       ----------------
!       LOCAL VARIABLES.
!       ----------------
INTEGER             KCOL  , KSUP  , LSTSUP
!
!***********************************************************************
!
!       -------------------------------------------------
!       COMPUTE THE SUPERNODE PARTITION VECTOR XSUPER(*).
!       -------------------------------------------------
LSTSUP = NSUPER + 1
DO  100  KCOL = NEQNS, 1, -1
KSUP = SNODE(KCOL)
IF  ( KSUP .NE. LSTSUP )  THEN
XSUPER(LSTSUP) = KCOL + 1
ENDIF
LSTSUP = KSUP
100 CONTINUE
XSUPER(1) = 1
!
RETURN
END
!***********************************************************************
!***********************************************************************
!
!   Version:        0.4
!   Last modified:  December 27, 1994
!   Authors:        Joseph W.H. Liu
!
!   Mathematical Sciences Section, Oak Ridge National Laboratory
!
!***********************************************************************
!***********************************************************************
!--- SPARSPAK-A (ANSI FORTRAN) RELEASE III --- NAME = GENMMD
!  (C)  UNIVERSITY OF WATERLOO   JANUARY 1984
!***********************************************************************
!***********************************************************************
!****     GENMMD ..... MULTIPLE MINIMUM EXTERNAL DEGREE     ************
!***********************************************************************
!***********************************************************************
!
!     PURPOSE - THIS ROUTINE IMPLEMENTS THE MINIMUM DEGREE
!        ALGORITHM.  IT MAKES USE OF THE IMPLICIT REPRESENTATION
!        OF ELIMINATION GRAPHS BY QUOTIENT GRAPHS, AND THE
!        NOTION OF INDISTINGUISHABLE NODES.  IT ALSO IMPLEMENTS
!        THE MODIFICATIONS BY MULTIPLE ELIMINATION AND MINIMUM
!        EXTERNAL DEGREE.
!        ---------------------------------------------
!        CAUTION - THE ADJACENCY VECTOR ADJNCY WILL BE
!        DESTROYED.
!        ---------------------------------------------
!
!     INPUT PARAMETERS -
!        NEQNS  - NUMBER OF EQUATIONS.
!        (XADJ,ADJNCY) - THE ADJACENCY STRUCTURE.
!        DELTA  - TOLERANCE VALUE FOR MULTIPLE ELIMINATION.
!        MAXINT - MAXIMUM MACHINE REPRESENTABLE (SHORT) INTEGER
!                 (ANY SMALLER ESTIMATE WILL DO) FOR MARKING
!                 NODES.
!
!     OUTPUT PARAMETERS -
!        PERM   - THE MINIMUM DEGREE ORDERING.
!        INVP   - THE INVERSE OF PERM.
!        NOFSUB - AN UPPER BOUND ON THE NUMBER OF NONZERO
!                 SUBSCRIPTS FOR THE COMPRESSED STORAGE SCHEME.
!
!     WORKING PARAMETERS -
!        DHEAD  - VECTOR FOR HEAD OF DEGREE LISTS.
!        INVP   - USED TEMPORARILY FOR DEGREE FORWARD LINK.
!        PERM   - USED TEMPORARILY FOR DEGREE BACKWARD LINK.
!        QSIZE  - VECTOR FOR SIZE OF SUPERNODES.
!        LLIST  - VECTOR FOR TEMPORARY LINKED LISTS.
!        MARKER - A TEMPORARY MARKER VECTOR.
!
!     PROGRAM SUBROUTINES -
!        MMDELM, MMDINT, MMDNUM, MMDUPD.
!
!***********************************************************************
!
SUBROUTINE  GENMMD ( NEQNS, XADJ, ADJNCY, INVP, PERM, &
& DELTA, DHEAD, QSIZE, LLIST, MARKER, &
& MAXINT, NOFSUB )
!
!***********************************************************************
!
INTEGER    ADJNCY(*), DHEAD(*) , INVP(*)  , LLIST(*) , &
& MARKER(*), PERM(*)  , QSIZE(*)
INTEGER    XADJ(*)
INTEGER    DELTA , EHEAD , I     , MAXINT, MDEG  , &
& MDLMT , MDNODE, NEQNS , NEXTMD, NOFSUB, &
& NUM, TAG
!
!***********************************************************************
!
IF  ( NEQNS .LE. 0 )  RETURN
!
!        ------------------------------------------------
!        INITIALIZATION FOR THE MINIMUM DEGREE ALGORITHM.
!        ------------------------------------------------
NOFSUB = 0
CALL  MMDINT ( NEQNS, XADJ, DHEAD, INVP, PERM, &
& QSIZE, LLIST, MARKER )
!
!        ----------------------------------------------
!        NUM COUNTS THE NUMBER OF ORDERED NODES PLUS 1.
!        ----------------------------------------------
NUM = 1
!
!        -----------------------------
!        ELIMINATE ALL ISOLATED NODES.
!        -----------------------------
NEXTMD = DHEAD(1)
100 CONTINUE
IF  ( NEXTMD .LE. 0 )  GO TO 200
MDNODE = NEXTMD
NEXTMD = INVP(MDNODE)
MARKER(MDNODE) = MAXINT
INVP(MDNODE) = - NUM
NUM = NUM + 1
GO TO 100
!
200 CONTINUE
!        ----------------------------------------
!        SEARCH FOR NODE OF THE MINIMUM DEGREE.
!        MDEG IS THE CURRENT MINIMUM DEGREE;
!        TAG IS USED TO FACILITATE MARKING NODES.
!        ----------------------------------------
IF  ( NUM .GT. NEQNS )  GO TO 1000
TAG = 1
DHEAD(1) = 0
MDEG = 2
300 CONTINUE
IF  ( DHEAD(MDEG) .GT. 0 )  GO TO 400
MDEG = MDEG + 1
GO TO 300
400 CONTINUE
!            -------------------------------------------------
!            USE VALUE OF DELTA TO SET UP MDLMT, WHICH GOVERNS
!            WHEN A DEGREE UPDATE IS TO BE PERFORMED.
!            -------------------------------------------------
MDLMT = MDEG + DELTA
EHEAD = 0
!
500 CONTINUE
MDNODE = DHEAD(MDEG)
IF  ( MDNODE .GT. 0 )  GO TO 600
MDEG = MDEG + 1
IF  ( MDEG .GT. MDLMT )  GO TO 900
GO TO 500
600 CONTINUE
!                ----------------------------------------
!                REMOVE MDNODE FROM THE DEGREE STRUCTURE.
!                ----------------------------------------
NEXTMD = INVP(MDNODE)
DHEAD(MDEG) = NEXTMD
IF  ( NEXTMD .GT. 0 )  PERM(NEXTMD) = - MDEG
INVP(MDNODE) = - NUM
NOFSUB = NOFSUB + MDEG + QSIZE(MDNODE) - 2
IF  ( NUM+QSIZE(MDNODE) .GT. NEQNS )  GO TO 1000
!                ----------------------------------------------
!                ELIMINATE MDNODE AND PERFORM QUOTIENT GRAPH
!                TRANSFORMATION.  RESET TAG VALUE IF NECESSARY.
!                ----------------------------------------------
TAG = TAG + 1
IF  ( TAG .LT. MAXINT )  GO TO 800
TAG = 1
DO  700  I = 1, NEQNS
IF  ( MARKER(I) .LT. MAXINT )  MARKER(I) = 0
700 CONTINUE
800 CONTINUE
CALL  MMDELM ( MDNODE, XADJ, ADJNCY, DHEAD, INVP, &
& PERM, QSIZE, LLIST, MARKER, MAXINT, &
& TAG )
NUM = NUM + QSIZE(MDNODE)
LLIST(MDNODE) = EHEAD
EHEAD = MDNODE
IF  ( DELTA .GE. 0 )  GO TO 500
900 CONTINUE
!            -------------------------------------------
!            UPDATE DEGREES OF THE NODES INVOLVED IN THE
!            MINIMUM DEGREE NODES ELIMINATION.
!            -------------------------------------------
IF  ( NUM .GT. NEQNS )  GO TO 1000
CALL  MMDUPD ( EHEAD, NEQNS, XADJ, ADJNCY, DELTA, MDEG, &
& DHEAD, INVP, PERM, QSIZE, LLIST, MARKER, &
& MAXINT, TAG )
GO TO 300
!
1000 CONTINUE
CALL  MMDNUM ( NEQNS, PERM, INVP, QSIZE )
RETURN
!
END
!***********************************************************************
!***********************************************************************
!
!   Version:        0.4
!   Last modified:  December 27, 1994
!   Authors:        Esmond G. Ng and Barry W. Peyton
!
!   Mathematical Sciences Section, Oak Ridge National Laboratory
!
!***********************************************************************
!***********************************************************************
!******         IGATHR .... INTEGER GATHER OPERATION      **************
!***********************************************************************
!***********************************************************************
!
!     PURPOSE - THIS ROUTINE PERFORMS A STANDARD INTEGER GATHER
!               OPERATION.
!
!     INPUT PARAMETERS -
!        KLEN   - LENGTH OF THE LIST OF GLOBAL INDICES.
!        LINDX  - LIST OF GLOBAL INDICES.
!        INDMAP - INDEXED BY GLOBAL INDICES, IT CONTAINS THE
!                 REQUIRED RELATIVE INDICES.
!
!     OUTPUT PARAMETERS -
!        RELIND - LIST RELATIVE INDICES.
!
!***********************************************************************
!
SUBROUTINE  IGATHR ( KLEN  , LINDX, INDMAP, RELIND )
!
!***********************************************************************
!
!     -----------
!     PARAMETERS.
!     -----------
INTEGER             KLEN
INTEGER             INDMAP(*), LINDX (*), RELIND(*)
!
!     ----------------
!     LOCAL VARIABLES.
!     ----------------
INTEGER             I
!
!***********************************************************************
!
!DIR$ IVDEP
DO  100  I = 1, KLEN
RELIND(I) = INDMAP(LINDX(I))
100 CONTINUE
RETURN
END
!***********************************************************************
!***********************************************************************
!
!   Version:        0.4
!   Last modified:  December 27, 1994
!   Authors:        Esmond G. Ng and Barry W. Peyton
!
!   Mathematical Sciences Section, Oak Ridge National Laboratory
!
!***********************************************************************
!***********************************************************************
!
!     ------------------------------------------------------
!     INPUT NUMERICAL VALUES INTO SPARSE DATA STRUCTURES ...
!     ------------------------------------------------------
!
SUBROUTINE  INPNV  (  XADJF, ADJF, ANZF, PERM, INVP, &
& NSUPER, XSUPER, XLINDX, LINDX, &
& XLNZ, LNZ, OFFSET )
!
INTEGER             XADJF(*), ADJF(*)
DOUBLE PRECISION    ANZF(*)
INTEGER             PERM(*), INVP(*)
INTEGER             NSUPER
INTEGER             XSUPER(*), XLINDX(*), LINDX(*)
INTEGER             XLNZ(*)
DOUBLE PRECISION    LNZ(*)
INTEGER             OFFSET(*)
!
INTEGER             I, II, J, JLEN, JSUPER, LAST, OLDJ
!
DO  500  JSUPER = 1, NSUPER
!
!           ----------------------------------------
!           FOR EACH SUPERNODE, DO THE FOLLOWING ...
!           ----------------------------------------
!
!           -----------------------------------------------
!           FIRST GET OFFSET TO FACILITATE NUMERICAL INPUT.
!           -----------------------------------------------
JLEN = XLINDX(JSUPER+1) - XLINDX(JSUPER)
DO  100  II = XLINDX(JSUPER), XLINDX(JSUPER+1)-1
I = LINDX(II)
JLEN = JLEN - 1
OFFSET(I) = JLEN
100 CONTINUE
!
DO  400  J = XSUPER(JSUPER), XSUPER(JSUPER+1)-1
!               -----------------------------------------
!               FOR EACH COLUMN IN THE CURRENT SUPERNODE,
!               FIRST INITIALIZE THE DATA STRUCTURE.
!               -----------------------------------------
DO  200  II = XLNZ(J), XLNZ(J+1)-1
LNZ(II) = 0.0
200 CONTINUE
!     re-commented the previous lines because of vector_dc!
!   C        The previous lines are not required as R initializes the arrays
!   C        Reinhard Furrer, Nov 19, 2007
!
!               -----------------------------------
!               NEXT INPUT THE INDIVIDUAL NONZEROS.
!               -----------------------------------
OLDJ = PERM(J)
LAST = XLNZ(J+1) - 1
DO  300  II = XADJF(OLDJ), XADJF(OLDJ+1)-1
I = INVP(ADJF(II))
IF  ( I .GE. J )  THEN
LNZ(LAST-OFFSET(I)) = ANZF(II)
ENDIF
300 CONTINUE
400 CONTINUE
!
500 CONTINUE
RETURN
END
!***********************************************************************
!***********************************************************************
!
!   Version:        0.4
!   Last modified:  December 27, 1994
!   Authors:        Joseph W.H. Liu
!
!   Mathematical Sciences Section, Oak Ridge National Laboratory
!
!***********************************************************************
!***********************************************************************
!***********     INVINV ..... CONCATENATION OF TWO INVP     ************
!***********************************************************************
!***********************************************************************
!
!   WRITTEN BY JOSEPH LIU (JUL 17, 1985)
!
!   PURPOSE:
!       TO PERFORM THE MAPPING OF
!           ORIGINAL-INVP --> INTERMEDIATE-INVP --> NEW INVP
!       AND THE RESULTING ORDERING REPLACES INVP.  THE NEW PERMUTATION
!       VECTOR PERM IS ALSO COMPUTED.
!
!   INPUT PARAMETERS:
!       NEQNS           -   NUMBER OF EQUATIONS.
!       INVP2           -   THE SECOND INVERSE PERMUTATION VECTOR.
!
!   UPDATED PARAMETERS:
!       INVP            -   THE FIRST INVERSE PERMUTATION VECTOR.  ON
!                           OUTPUT, IT CONTAINS THE NEW INVERSE
!                           PERMUTATION.
!
!   OUTPUT PARAMETER:
!       PERM            -   NEW PERMUTATION VECTOR (CAN BE THE SAME AS
!                           INVP2).
!
!***********************************************************************
!
SUBROUTINE  INVINV (  NEQNS , INVP  , INVP2 , PERM            )
!
!***********************************************************************
!
INTEGER           INVP(*)       , INVP2(*)      , &
& PERM(*)
!
INTEGER           NEQNS
!
!***********************************************************************
!
INTEGER           I     , INTERM, NODE
!
!***********************************************************************
!
DO  100  I = 1, NEQNS
INTERM = INVP(I)
INVP(I) = INVP2(INTERM)
100 CONTINUE
!
DO  200  I = 1, NEQNS
NODE = INVP(I)
PERM(NODE) = I
200 CONTINUE
!
RETURN
END
!***********************************************************************
!***********************************************************************
!
!   Version:        0.4
!   Last modified:  December 27, 1994
!   Authors:        Esmond G. Ng and Barry W. Peyton
!
!   Mathematical Sciences Section, Oak Ridge National Laboratory
!
!***********************************************************************
!***********************************************************************
!******         LDINDX .... LOAD INDEX VECTOR             **************
!***********************************************************************
!***********************************************************************
!
!     PURPOSE - THIS ROUTINE COMPUTES THE SECOND INDEX VECTOR
!               USED TO IMPLEMENT THE DOUBLY-INDIRECT SAXPY-LIKE
!               LOOPS THAT ALLOW US TO ACCUMULATE UPDATE
!               COLUMNS DIRECTLY INTO FACTOR STORAGE.
!
!     INPUT PARAMETERS -
!        JLEN   - LENGTH OF THE FIRST COLUMN OF THE SUPERNODE,
!                 INCLUDING THE DIAGONAL ENTRY.
!        LINDX  - THE OFF-DIAGONAL ROW INDICES OF THE SUPERNODE,
!                 I.E., THE ROW INDICES OF THE NONZERO ENTRIES
!                 LYING BELOW THE DIAGONAL ENTRY OF THE FIRST
!                 COLUMN OF THE SUPERNODE.
!
!     OUTPUT PARAMETERS -
!        INDMAP - THIS INDEX VECTOR MAPS EVERY GLOBAL ROW INDEX
!                 OF NONZERO ENTRIES IN THE FIRST COLUMN OF THE
!                 SUPERNODE TO ITS POSITION IN THE INDEX LIST
!                 RELATIVE TO THE LAST INDEX IN THE LIST.  MORE
!                 PRECISELY, IT GIVES THE DISTANCE OF EACH INDEX
!                 FROM THE LAST INDEX IN THE LIST.
!
!***********************************************************************
!
SUBROUTINE  LDINDX ( JLEN, LINDX, INDMAP )
!
!***********************************************************************
!
!     -----------
!     PARAMETERS.
!     -----------
INTEGER             JLEN
INTEGER             LINDX(*), INDMAP(*)
!
!     ----------------
!     LOCAL VARIABLES.
!     ----------------
INTEGER             CURLEN, J, JSUB
!
!***********************************************************************
!

CURLEN = JLEN

DO  200  J = 1, JLEN
JSUB = LINDX(J)
CURLEN = CURLEN - 1
INDMAP(JSUB) = CURLEN
200 CONTINUE
RETURN
END
!***********************************************************************
!***********************************************************************
!
!   Version:        0.4
!   Last modified:  December 27, 1994
!   Authors:        Joseph W.H. Liu
!
!   Mathematical Sciences Section, Oak Ridge National Laboratory
!
!***********************************************************************
!***********************************************************************
!--- SPARSPAK-A (ANSI FORTRAN) RELEASE III --- NAME = MMDELM
!  (C)  UNIVERSITY OF WATERLOO   JANUARY 1984
!***********************************************************************
!***********************************************************************
!**     MMDELM ..... MULTIPLE MINIMUM DEGREE ELIMINATION     ***********
!***********************************************************************
!***********************************************************************
!
!     PURPOSE - THIS ROUTINE ELIMINATES THE NODE MDNODE OF
!        MINIMUM DEGREE FROM THE ADJACENCY STRUCTURE, WHICH
!        IS STORED IN THE QUOTIENT GRAPH FORMAT.  IT ALSO
!        TRANSFORMS THE QUOTIENT GRAPH REPRESENTATION OF THE
!        ELIMINATION GRAPH.
!
!     INPUT PARAMETERS -
!        MDNODE - NODE OF MINIMUM DEGREE.
!        MAXINT - ESTIMATE OF MAXIMUM REPRESENTABLE (SHORT)
!                 INTEGER.
!        TAG    - TAG VALUE.
!
!     UPDATED PARAMETERS -
!        (XADJ,ADJNCY) - UPDATED ADJACENCY STRUCTURE.
!        (DHEAD,DFORW,DBAKW) - DEGREE DOUBLY LINKED STRUCTURE.
!        QSIZE  - SIZE OF SUPERNODE.
!        MARKER - MARKER VECTOR.
!        LLIST  - TEMPORARY LINKED LIST OF ELIMINATED NABORS.
!
!***********************************************************************
!
SUBROUTINE  MMDELM ( MDNODE, XADJ, ADJNCY, DHEAD, DFORW, &
& DBAKW, QSIZE, LLIST, MARKER, MAXINT, &
& TAG )
!
!***********************************************************************
!
INTEGER    ADJNCY(*), DBAKW(*) , DFORW(*) , DHEAD(*) , &
& LLIST(*) , MARKER(*), QSIZE(*)
INTEGER    XADJ(*)
INTEGER    ELMNT , I     , ISTOP , ISTRT , J     , &
& JSTOP , JSTRT , LINK  , MAXINT, MDNODE, &
& NABOR , NODE  , NPV   , NQNBRS, NXNODE, &
& PVNODE, RLMT  , RLOC  , RNODE , TAG   , &
& XQNBR
!
!***********************************************************************
!
!        -----------------------------------------------
!        FIND REACHABLE SET AND PLACE IN DATA STRUCTURE.
!        -----------------------------------------------
MARKER(MDNODE) = TAG
ISTRT = XADJ(MDNODE)
ISTOP = XADJ(MDNODE+1) - 1
!        -------------------------------------------------------
!        ELMNT POINTS TO THE BEGINNING OF THE LIST OF ELIMINATED
!        NABORS OF MDNODE, AND RLOC GIVES THE STORAGE LOCATION
!        FOR THE NEXT REACHABLE NODE.
!        -------------------------------------------------------
ELMNT = 0
RLOC = ISTRT
RLMT = ISTOP
DO  200  I = ISTRT, ISTOP
NABOR = ADJNCY(I)
IF  ( NABOR .EQ. 0 )  GO TO 300
IF  ( MARKER(NABOR) .GE. TAG )  GO TO 200
MARKER(NABOR) = TAG
IF  ( DFORW(NABOR) .LT. 0 )  GO TO 100
ADJNCY(RLOC) = NABOR
RLOC = RLOC + 1
GO TO 200
100 CONTINUE
LLIST(NABOR) = ELMNT
ELMNT = NABOR
200 CONTINUE
300 CONTINUE
!            -----------------------------------------------------
!            MERGE WITH REACHABLE NODES FROM GENERALIZED ELEMENTS.
!            -----------------------------------------------------
IF  ( ELMNT .LE. 0 )  GO TO 1000
ADJNCY(RLMT) = - ELMNT
LINK = ELMNT
400 CONTINUE
JSTRT = XADJ(LINK)
JSTOP = XADJ(LINK+1) - 1
DO  800  J = JSTRT, JSTOP
NODE = ADJNCY(J)
LINK = - NODE
if ( NODE < 0 )  then
GO TO 400
else if ( NODE == 0 ) then
GO TO 900
else
GO TO 500
end if

500 CONTINUE
IF  ( MARKER(NODE) .GE. TAG  .OR. &
& DFORW(NODE) .LT. 0 )  GO TO 800
MARKER(NODE) = TAG
!                            ---------------------------------
!                            USE STORAGE FROM ELIMINATED NODES
!                            IF NECESSARY.
!                            ---------------------------------
600 CONTINUE
IF  ( RLOC .LT. RLMT )  GO TO 700
LINK = - ADJNCY(RLMT)
RLOC = XADJ(LINK)
RLMT = XADJ(LINK+1) - 1
GO TO 600
700 CONTINUE
ADJNCY(RLOC) = NODE
RLOC = RLOC + 1
800 CONTINUE
900 CONTINUE
ELMNT = LLIST(ELMNT)
GO TO 300
1000 CONTINUE
IF  ( RLOC .LE. RLMT )  ADJNCY(RLOC) = 0
!        --------------------------------------------------------
!        FOR EACH NODE IN THE REACHABLE SET, DO THE FOLLOWING ...
!        --------------------------------------------------------
LINK = MDNODE
1100 CONTINUE
ISTRT = XADJ(LINK)
ISTOP = XADJ(LINK+1) - 1
DO  1700  I = ISTRT, ISTOP
RNODE = ADJNCY(I)
LINK = - RNODE

if ( RNODE < 0 )  then
GO TO 1100
else if ( RNODE == 0 ) then
GO TO 1800
else
GO TO 1200
end if

1200 CONTINUE
!                --------------------------------------------
!                IF RNODE IS IN THE DEGREE LIST STRUCTURE ...
!                --------------------------------------------
PVNODE = DBAKW(RNODE)
IF  ( PVNODE .EQ. 0  .OR. &
& PVNODE .EQ. (-MAXINT) )  GO TO 1300
!                    -------------------------------------
!                    THEN REMOVE RNODE FROM THE STRUCTURE.
!                    -------------------------------------
NXNODE = DFORW(RNODE)
IF  ( NXNODE .GT. 0 )  DBAKW(NXNODE) = PVNODE
IF  ( PVNODE .GT. 0 )  DFORW(PVNODE) = NXNODE
NPV = - PVNODE
IF  ( PVNODE .LT. 0 )  DHEAD(NPV) = NXNODE
1300 CONTINUE
!                ----------------------------------------
!                PURGE INACTIVE QUOTIENT NABORS OF RNODE.
!                ----------------------------------------
JSTRT = XADJ(RNODE)
JSTOP = XADJ(RNODE+1) - 1
XQNBR = JSTRT
DO  1400  J = JSTRT, JSTOP
NABOR = ADJNCY(J)
IF  ( NABOR .EQ. 0 )  GO TO 1500
IF  ( MARKER(NABOR) .GE. TAG )  GO TO 1400
ADJNCY(XQNBR) = NABOR
XQNBR = XQNBR + 1
1400 CONTINUE
1500 CONTINUE
!                ----------------------------------------
!                IF NO ACTIVE NABOR AFTER THE PURGING ...
!                ----------------------------------------
NQNBRS = XQNBR - JSTRT
IF  ( NQNBRS .GT. 0 )  GO TO 1600
!                    -----------------------------
!                    THEN MERGE RNODE WITH MDNODE.
!                    -----------------------------
QSIZE(MDNODE) = QSIZE(MDNODE) + QSIZE(RNODE)
QSIZE(RNODE) = 0
MARKER(RNODE) = MAXINT
DFORW(RNODE) = - MDNODE
DBAKW(RNODE) = - MAXINT
GO TO 1700
1600 CONTINUE
!                --------------------------------------
!                ELSE FLAG RNODE FOR DEGREE UPDATE, AND
!                ADD MDNODE AS A NABOR OF RNODE.
!                --------------------------------------
DFORW(RNODE) = NQNBRS + 1
DBAKW(RNODE) = 0
ADJNCY(XQNBR) = MDNODE
XQNBR = XQNBR + 1
IF  ( XQNBR .LE. JSTOP )  ADJNCY(XQNBR) = 0
!
1700 CONTINUE
1800 CONTINUE
RETURN
!
END
!***********************************************************************
!***********************************************************************
!
!   Version:        0.4
!   Last modified:  December 27, 1994
!   Authors:        Joseph W.H. Liu
!
!   Mathematical Sciences Section, Oak Ridge National Laboratory
!
!***********************************************************************
!***********************************************************************
!--- SPARSPAK-A (ANSI FORTRAN) RELEASE III --- NAME = MMDINT
!  (C)  UNIVERSITY OF WATERLOO   JANUARY 1984
!***********************************************************************
!***********************************************************************
!***     MMDINT ..... MULT MINIMUM DEGREE INITIALIZATION     ***********
!***********************************************************************
!***********************************************************************
!
!     PURPOSE - THIS ROUTINE PERFORMS INITIALIZATION FOR THE
!        MULTIPLE ELIMINATION VERSION OF THE MINIMUM DEGREE
!        ALGORITHM.
!
!     INPUT PARAMETERS -
!        NEQNS  - NUMBER OF EQUATIONS.
!        (XADJ,ADJNCY) - ADJACENCY STRUCTURE.
!
!     OUTPUT PARAMETERS -
!        (DHEAD,DFORW,DBAKW) - DEGREE DOUBLY LINKED STRUCTURE.
!        QSIZE  - SIZE OF SUPERNODE (INITIALIZED TO ONE).
!        LLIST  - LINKED LIST.
!        MARKER - MARKER VECTOR.
!
!***********************************************************************
!
SUBROUTINE  MMDINT ( NEQNS, XADJ, DHEAD, DFORW, &
& DBAKW, QSIZE, LLIST, MARKER )
!
!***********************************************************************
!
INTEGER    DBAKW(*) , DFORW(*) , DHEAD(*) , &
& LLIST(*) , MARKER(*), QSIZE(*)
INTEGER    XADJ(*)
INTEGER    FNODE , NDEG  , NEQNS , NODE
!
!***********************************************************************
!
DO  100  NODE = 1, NEQNS
DHEAD(NODE) = 0
QSIZE(NODE) = 1
MARKER(NODE) = 0
LLIST(NODE) = 0
100 CONTINUE
!        ------------------------------------------
!        INITIALIZE THE DEGREE DOUBLY LINKED LISTS.
!        ------------------------------------------
DO  200  NODE = 1, NEQNS
NDEG = XADJ(NODE+1) - XADJ(NODE) + 1
FNODE = DHEAD(NDEG)
DFORW(NODE) = FNODE
DHEAD(NDEG) = NODE
IF  ( FNODE .GT. 0 )  DBAKW(FNODE) = NODE
DBAKW(NODE) = - NDEG
200 CONTINUE
RETURN
!
END
!***********************************************************************
!***********************************************************************
!
!   Version:        0.4
!   Last modified:  December 27, 1994
!   Authors:        Joseph W.H. Liu
!
!   Mathematical Sciences Section, Oak Ridge National Laboratory
!
!***********************************************************************
!***********************************************************************
!--- SPARSPAK-A (ANSI FORTRAN) RELEASE III --- NAME = MMDNUM
!  (C)  UNIVERSITY OF WATERLOO   JANUARY 1984
!***********************************************************************
!***********************************************************************
!*****     MMDNUM ..... MULTI MINIMUM DEGREE NUMBERING     *************
!***********************************************************************
!***********************************************************************
!
!     PURPOSE - THIS ROUTINE PERFORMS THE FINAL STEP IN
!        PRODUCING THE PERMUTATION AND INVERSE PERMUTATION
!        VECTORS IN THE MULTIPLE ELIMINATION VERSION OF THE
!        MINIMUM DEGREE ORDERING ALGORITHM.
!
!     INPUT PARAMETERS -
!        NEQNS  - NUMBER OF EQUATIONS.
!        QSIZE  - SIZE OF SUPERNODES AT ELIMINATION.
!
!     UPDATED PARAMETERS -
!        INVP   - INVERSE PERMUTATION VECTOR.  ON INPUT,
!                 IF QSIZE(NODE)=0, THEN NODE HAS BEEN MERGED
!                 INTO THE NODE -INVP(NODE); OTHERWISE,
!                 -INVP(NODE) IS ITS INVERSE LABELLING.
!
!     OUTPUT PARAMETERS -
!        PERM   - THE PERMUTATION VECTOR.
!
!***********************************************************************
!
SUBROUTINE  MMDNUM ( NEQNS, PERM, INVP, QSIZE )
!
!***********************************************************************
!
INTEGER    INVP(*)  , PERM(*)  , QSIZE(*)
INTEGER    FATHER, NEQNS , NEXTF , NODE  , NQSIZE, &
& NUM   , ROOT
!
!***********************************************************************
!
DO  100  NODE = 1, NEQNS
NQSIZE = QSIZE(NODE)
IF  ( NQSIZE .LE. 0 )  PERM(NODE) = INVP(NODE)
IF  ( NQSIZE .GT. 0 )  PERM(NODE) = - INVP(NODE)
100 CONTINUE
!        ------------------------------------------------------
!        FOR EACH NODE WHICH HAS BEEN MERGED, DO THE FOLLOWING.
!        ------------------------------------------------------
DO  500  NODE = 1, NEQNS
IF  ( PERM(NODE) .GT. 0 )  GO TO 500
!                -----------------------------------------
!                TRACE THE MERGED TREE UNTIL ONE WHICH HAS
!                NOT BEEN MERGED, CALL IT ROOT.
!                -----------------------------------------
FATHER = NODE
200 CONTINUE
IF  ( PERM(FATHER) .GT. 0 )  GO TO 300
FATHER = - PERM(FATHER)
GO TO 200
300 CONTINUE
!                -----------------------
!                NUMBER NODE AFTER ROOT.
!                -----------------------
ROOT = FATHER
NUM = PERM(ROOT) + 1
INVP(NODE) = - NUM
PERM(ROOT) = NUM
!                ------------------------
!                SHORTEN THE MERGED TREE.
!                ------------------------
FATHER = NODE
400 CONTINUE
NEXTF = - PERM(FATHER)
IF  ( NEXTF .LE. 0 )  GO TO 500
PERM(FATHER) = - ROOT
FATHER = NEXTF
GO TO 400
500 CONTINUE
!        ----------------------
!        READY TO COMPUTE PERM.
!        ----------------------
DO  600  NODE = 1, NEQNS
NUM = - INVP(NODE)
INVP(NODE) = NUM
PERM(NUM) = NODE
600 CONTINUE
RETURN
!
END
!***********************************************************************
!***********************************************************************
!
!   Version:        0.4
!   Last modified:  December 27, 1994
!   Authors:        Joseph W.H. Liu
!
!   Mathematical Sciences Section, Oak Ridge National Laboratory
!
!***********************************************************************
!***********************************************************************
!--- SPARSPAK-A (ANSI FORTRAN) RELEASE III --- NAME = MMDUPD
!  (C)  UNIVERSITY OF WATERLOO   JANUARY 1984
!***********************************************************************
!***********************************************************************
!*****     MMDUPD ..... MULTIPLE MINIMUM DEGREE UPDATE     *************
!***********************************************************************
!***********************************************************************
!
!     PURPOSE - THIS ROUTINE UPDATES THE DEGREES OF NODES
!        AFTER A MULTIPLE ELIMINATION STEP.
!
!     INPUT PARAMETERS -
!        EHEAD  - THE BEGINNING OF THE LIST OF ELIMINATED
!                 NODES (I.E., NEWLY FORMED ELEMENTS).
!        NEQNS  - NUMBER OF EQUATIONS.
!        (XADJ,ADJNCY) - ADJACENCY STRUCTURE.
!        DELTA  - TOLERANCE VALUE FOR MULTIPLE ELIMINATION.
!        MAXINT - MAXIMUM MACHINE REPRESENTABLE (SHORT)
!                 INTEGER.
!
!     UPDATED PARAMETERS -
!        MDEG   - NEW MINIMUM DEGREE AFTER DEGREE UPDATE.
!        (DHEAD,DFORW,DBAKW) - DEGREE DOUBLY LINKED STRUCTURE.
!        QSIZE  - SIZE OF SUPERNODE.
!        LLIST  - WORKING LINKED LIST.
!        MARKER - MARKER VECTOR FOR DEGREE UPDATE.
!        TAG    - TAG VALUE.
!
!***********************************************************************
!
SUBROUTINE  MMDUPD ( EHEAD, NEQNS, XADJ, ADJNCY, DELTA, &
& MDEG, DHEAD, DFORW, DBAKW, QSIZE, &
& LLIST, MARKER, MAXINT, TAG )
!
!***********************************************************************
!
INTEGER    ADJNCY(*), DBAKW(*) , DFORW(*) , DHEAD(*) , &
& LLIST(*) , MARKER(*), QSIZE(*)
INTEGER    XADJ(*)
INTEGER    DEG   , DEG0  , DELTA , EHEAD , ELMNT , &
& ENODE , FNODE , I     , IQ2   , ISTOP , &
& ISTRT , J     , JSTOP , JSTRT , LINK  , &
& MAXINT, MDEG  , MDEG0 , MTAG  , NABOR , &
& NEQNS , NODE  , Q2HEAD, QXHEAD, TAG
!
!***********************************************************************
!
MDEG0 = MDEG + DELTA
ELMNT = EHEAD
100 CONTINUE
!            -------------------------------------------------------
!            FOR EACH OF THE NEWLY FORMED ELEMENT, DO THE FOLLOWING.
!            (RESET TAG VALUE IF NECESSARY.)
!            -------------------------------------------------------
IF  ( ELMNT .LE. 0 )  RETURN
MTAG = TAG + MDEG0
IF  ( MTAG .LT. MAXINT )  GO TO 300
TAG = 1
DO  200  I = 1, NEQNS
IF  ( MARKER(I) .LT. MAXINT )  MARKER(I) = 0
200 CONTINUE
MTAG = TAG + MDEG0
300 CONTINUE
!            ---------------------------------------------
!            CREATE TWO LINKED LISTS FROM NODES ASSOCIATED
!            WITH ELMNT: ONE WITH TWO NABORS (Q2HEAD) IN
!            ADJACENCY STRUCTURE, AND THE OTHER WITH MORE
!            THAN TWO NABORS (QXHEAD).  ALSO COMPUTE DEG0,
!            NUMBER OF NODES IN THIS ELEMENT.
!            ---------------------------------------------
Q2HEAD = 0
QXHEAD = 0
DEG0 = 0
LINK = ELMNT
400 CONTINUE
ISTRT = XADJ(LINK)
ISTOP = XADJ(LINK+1) - 1
DO  700  I = ISTRT, ISTOP
ENODE = ADJNCY(I)
LINK = - ENODE

if ( ENODE < 0 )  then
GO TO 400
else if ( ENODE == 0 ) then
GO TO 800
else
GO TO 500
end if
!
500 CONTINUE
IF  ( QSIZE(ENODE) .EQ. 0 )  GO TO 700
DEG0 = DEG0 + QSIZE(ENODE)
MARKER(ENODE) = MTAG
!                        ----------------------------------
!                        IF ENODE REQUIRES A DEGREE UPDATE,
!                        THEN DO THE FOLLOWING.
!                        ----------------------------------
IF  ( DBAKW(ENODE) .NE. 0 )  GO TO 700
!                            ---------------------------------------
!                            PLACE EITHER IN QXHEAD OR Q2HEAD LISTS.
!                            ---------------------------------------
IF  ( DFORW(ENODE) .EQ. 2 )  GO TO 600
LLIST(ENODE) = QXHEAD
QXHEAD = ENODE
GO TO 700
600 CONTINUE
LLIST(ENODE) = Q2HEAD
Q2HEAD = ENODE
700 CONTINUE
800 CONTINUE
!            --------------------------------------------
!            FOR EACH ENODE IN Q2 LIST, DO THE FOLLOWING.
!            --------------------------------------------
ENODE = Q2HEAD
IQ2 = 1
900 CONTINUE
IF  ( ENODE .LE. 0 )  GO TO 1500
IF  ( DBAKW(ENODE) .NE. 0 )  GO TO 2200
TAG = TAG + 1
DEG = DEG0
!                    ------------------------------------------
!                    IDENTIFY THE OTHER ADJACENT ELEMENT NABOR.
!                    ------------------------------------------
ISTRT = XADJ(ENODE)
NABOR = ADJNCY(ISTRT)
IF  ( NABOR .EQ. ELMNT )  NABOR = ADJNCY(ISTRT+1)
!                    ------------------------------------------------
!                    IF NABOR IS UNELIMINATED, INCREASE DEGREE COUNT.
!                    ------------------------------------------------
LINK = NABOR
IF  ( DFORW(NABOR) .LT. 0 )  GO TO 1000
DEG = DEG + QSIZE(NABOR)
GO TO 2100
1000 CONTINUE
!                        --------------------------------------------
!                        OTHERWISE, FOR EACH NODE IN THE 2ND ELEMENT,
!                        DO THE FOLLOWING.
!                        --------------------------------------------
ISTRT = XADJ(LINK)
ISTOP = XADJ(LINK+1) - 1
DO  1400  I = ISTRT, ISTOP
NODE = ADJNCY(I)
LINK = - NODE
IF  ( NODE .EQ. ENODE )  GO TO 1400
!.ARITH.if                       IF  ( NODE )  1000, 2100, 1100
if ( NODE < 0 )  then
GO TO 1000
else if ( NODE == 0 ) then
GO TO 2100
else
GO TO 1100
end if
!
1100 CONTINUE
IF  ( QSIZE(NODE) .EQ. 0 )  GO TO 1400
IF  ( MARKER(NODE) .GE. TAG )  GO TO 1200
!                                -------------------------------------
!                                CASE WHEN NODE IS NOT YET CONSIDERED.
!                                -------------------------------------
MARKER(NODE) = TAG
DEG = DEG + QSIZE(NODE)
GO TO 1400
1200 CONTINUE
!                            ----------------------------------------
!                            CASE WHEN NODE IS INDISTINGUISHABLE FROM
!                            ENODE.  MERGE THEM INTO A NEW SUPERNODE.
!                            ----------------------------------------
IF  ( DBAKW(NODE) .NE. 0 )  GO TO 1400
IF  ( DFORW(NODE) .NE. 2 )  GO TO 1300
QSIZE(ENODE) = QSIZE(ENODE) + &
& QSIZE(NODE)
QSIZE(NODE) = 0
MARKER(NODE) = MAXINT
DFORW(NODE) = - ENODE
DBAKW(NODE) = - MAXINT
GO TO 1400
1300 CONTINUE
!                            --------------------------------------
!                            CASE WHEN NODE IS OUTMATCHED BY ENODE.
!                            --------------------------------------
IF  ( DBAKW(NODE) .EQ.0 ) &
& DBAKW(NODE) = - MAXINT
1400 CONTINUE
GO TO 2100
1500 CONTINUE
!                ------------------------------------------------
!                FOR EACH ENODE IN THE QX LIST, DO THE FOLLOWING.
!                ------------------------------------------------
ENODE = QXHEAD
IQ2 = 0
1600 CONTINUE
IF  ( ENODE .LE. 0 )  GO TO 2300
IF  ( DBAKW(ENODE) .NE. 0 )  GO TO 2200
TAG = TAG + 1
DEG = DEG0
!                        ---------------------------------
!                        FOR EACH UNMARKED NABOR OF ENODE,
!                        DO THE FOLLOWING.
!                        ---------------------------------
ISTRT = XADJ(ENODE)
ISTOP = XADJ(ENODE+1) - 1
DO  2000  I = ISTRT, ISTOP
NABOR = ADJNCY(I)
IF  ( NABOR .EQ. 0 )  GO TO 2100
IF  ( MARKER(NABOR) .GE. TAG )  GO TO 2000
MARKER(NABOR) = TAG
LINK = NABOR
!                                ------------------------------
!                                IF UNELIMINATED, INCLUDE IT IN
!                                DEG COUNT.
!                                ------------------------------
IF  ( DFORW(NABOR) .LT. 0 )  GO TO 1700
DEG = DEG + QSIZE(NABOR)
GO TO 2000
1700 CONTINUE
!                                    -------------------------------
!                                    IF ELIMINATED, INCLUDE UNMARKED
!                                    NODES IN THIS ELEMENT INTO THE
!                                    DEGREE COUNT.
!                                    -------------------------------
JSTRT = XADJ(LINK)
JSTOP = XADJ(LINK+1) - 1
DO  1900  J = JSTRT, JSTOP
NODE = ADJNCY(J)
LINK = - NODE
!
if ( NODE < 0 )  then
GO TO 1700
else if ( NODE == 0 ) then
GO TO 2000
else
GO TO 1800
end if
!
1800 CONTINUE
IF  ( MARKER(NODE) .GE. TAG ) &
& GO TO 1900
MARKER(NODE) = TAG
DEG = DEG + QSIZE(NODE)
1900 CONTINUE
2000 CONTINUE
2100 CONTINUE
!                    -------------------------------------------
!                    UPDATE EXTERNAL DEGREE OF ENODE IN DEGREE
!                    STRUCTURE, AND MDEG (MIN DEG) IF NECESSARY.
!                    -------------------------------------------
DEG = DEG - QSIZE(ENODE) + 1
FNODE = DHEAD(DEG)
DFORW(ENODE) = FNODE
DBAKW(ENODE) = - DEG
IF  ( FNODE .GT. 0 )  DBAKW(FNODE) = ENODE
DHEAD(DEG) = ENODE
IF  ( DEG .LT. MDEG )  MDEG = DEG
2200 CONTINUE
!                    ----------------------------------
!                    GET NEXT ENODE IN CURRENT ELEMENT.
!                    ----------------------------------
ENODE = LLIST(ENODE)
IF  ( IQ2 .EQ. 1 )  GO TO 900
GO TO 1600
2300 CONTINUE
!            -----------------------------
!            GET NEXT ELEMENT IN THE LIST.
!            -----------------------------
TAG = MTAG
ELMNT = LLIST(ELMNT)
GO TO 100
!
END
!***********************************************************************
!***********************************************************************
!
!   Version:        0.4
!   Last modified:  December 27, 1994
!   Authors:        Esmond G. Ng and Barry W. Peyton
!                   RF: modified mmpy8 dependence
!
!   Mathematical Sciences Section, Oak Ridge National Laboratory
!
!***********************************************************************
!***********************************************************************
!**************     MMPY  .... MATRIX-MATRIX MULTIPLY     **************
!***********************************************************************
!***********************************************************************
!
!   PURPOSE -
!       THIS ROUTINE PERFORMS A MATRIX-MATRIX MULTIPLY, Y = Y + XA,
!       ASSUMING DATA STRUCTURES USED IN SOME OF OUR SPARSE CHOLESKY
!       CODES.
!
!   INPUT PARAMETERS -
!       M               -   NUMBER OF ROWS IN X AND IN Y.
!       N               -   NUMBER OF COLUMNS IN X AND NUMBER OF ROWS
!                           IN A.
!       Q               -   NUMBER OF COLUMNS IN A AND Y.
!       SPLIT(*)        -   BLOCK PARTITIONING OF X.
!       XPNT(*)         -   XPNT(J+1) POINTS ONE LOCATION BEYOND THE
!                           END OF THE J-TH COLUMN OF X.  XPNT IS ALSO
!                           USED TO ACCESS THE ROWS OF A.
!       X(*)            -   CONTAINS THE COLUMNS OF X AND THE ROWS OF A.
!       LDY             -   LENGTH OF FIRST COLUMN OF Y.
!
!   EXTERNAL ROUTINES:
!       MMPYN           -   MATRIX-MATRIX MULTIPLY,
!                           WITH LEVEL 8 LOOP UNROLLING.
!
!   UPDATED PARAMETERS -
!       Y(*)            -   ON OUTPUT, Y = Y + AX.
!
!***********************************************************************
!
SUBROUTINE  MMPY   (  M     , N     , Q     , SPLIT , XPNT  , &
& X     , Y     , LDY    )
!
!***********************************************************************
!
!       -----------
!       PARAMETERS.
!       -----------
!
EXTERNAL            MMPY8
INTEGER             LDY   , M     , N     , Q
INTEGER             SPLIT(*)      , XPNT(*)
DOUBLE PRECISION    X(*)          , Y(*)
!
!       ----------------
!       LOCAL VARIABLES.
!       ----------------
!
INTEGER             BLK   , FSTCOL, NN
!
!***********************************************************************
!
BLK = 1
FSTCOL = 1
100 CONTINUE
IF  ( FSTCOL .LE. N )  THEN
NN = SPLIT(BLK)
CALL  MMPY8 ( M, NN, Q, XPNT(FSTCOL), X, Y, LDY )
FSTCOL = FSTCOL + NN
BLK = BLK + 1
GO TO 100
ENDIF
RETURN
!
END
!***********************************************************************
!***********************************************************************
!
!   Version:        0.4
!   Last modified:  May 26, 1995
!   Authors:        Esmond G. Ng, Barry W. Peyton, and Guodong Zhang
!
!   Mathematical Sciences Section, Oak Ridge National Laboratory
!
!***********************************************************************
!***********************************************************************
!*************     MMPY8  .... MATRIX-MATRIX MULTIPLY     **************
!***********************************************************************
!***********************************************************************
!
!   PURPOSE -
!       THIS ROUTINE PERFORMS A MATRIX-MATRIX MULTIPLY, Y = Y + XA,
!       ASSUMING DATA STRUCTURES USED IN SOME OF OUR SPARSE CHOLESKY
!       CODES.
!
!       LOOP UNROLLING: LEVEL 8 UPDATING TWO COLUMNS AT A TIME
!
!   INPUT PARAMETERS -
!       M               -   NUMBER OF ROWS IN X AND IN Y.
!       N               -   NUMBER OF COLUMNS IN X AND NUMBER OF ROWS
!                           IN A.
!       Q               -   NUMBER OF COLUMNS IN A AND Y.
!       XPNT(*)         -   XPNT(J+1) POINTS ONE LOCATION BEYOND THE
!                           END OF THE J-TH COLUMN OF X.  XPNT IS ALSO
!                           USED TO ACCESS THE ROWS OF A.
!       X(*)            -   CONTAINS THE COLUMNS OF X AND THE ROWS OF A.
!       LDY             -   LENGTH OF FIRST COLUMN OF Y.
!
!   UPDATED PARAMETERS -
!       Y(*)            -   ON OUTPUT, Y = Y + AX.
!
!***********************************************************************
!
SUBROUTINE  MMPY8  (  M     , N     , Q     , XPNT  , X     , &
& Y     , LDY                             )
!
!***********************************************************************
!
!       -----------
!       PARAMETERS.
!       -----------
!
INTEGER               LDY   , M     , N     , Q
INTEGER               XPNT(*)
DOUBLE PRECISION      X(*)          , Y(*)
!
!       ----------------
!       LOCAL VARIABLES.
!       ----------------
!
INTEGER               I     , J     , K     , QQ
INTEGER               I1    , I2    , I3    , I4    , I5    , &
& I6    , I7    , I8
INTEGER               IYBEG , IYBEG1, IYBEG2, LENY  , MM
DOUBLE PRECISION      A1    , A2    , A3    , A4    , A5    , &
& A6    , A7    , A8    , A9    , A10   , &
& A11   , A12   , A13   , A14   , A15   , &
& A16
DOUBLE PRECISION      B1    , B2    , B3    , B4    , B5    , &
& B6    , B7    , B8    , Y1    , Y2
!
!***********************************************************************
!
!       ----------------------------------------------------
!       COMPUTE EACH DIAGONAL ENTRY OF THE ODD COLUMNS OF Y.
!       ----------------------------------------------------
!
MM = M
QQ = MIN(M,Q)
IYBEG = 1
LENY = LDY - 1
do J = 1, QQ-1 , 2
!DIR$   IVDEP
do I = 1, N
I1 = XPNT(I+1) - MM
A1 = X(I1)
Y(IYBEG) = Y(IYBEG) - A1*A1
end do
IYBEG = IYBEG + 2*LENY + 1
LENY = LENY - 2
MM = MM - 2
end do
!
!       -------------------------------------------------------
!       UPDATE TWO COLUMNS OF Y AT A TIME,  EXCEPT THE DIAGONAL
!       ELEMENT.
!       NOTE: THE DIAGONAL ELEMENT OF THE ODD COLUMN HAS
!             BEEN COMPUTED, SO WE COMPUTE THE SAME NUMBER OF
!             ELEMENTS FOR THE TWO COLUMNS.
!       -------------------------------------------------------
!
MM = M
IYBEG = 1
LENY = LDY - 1
!
do  J = 1, QQ-1, 2
!
IYBEG1 = IYBEG
IYBEG2 = IYBEG + LENY
!
do  K = 1, N-7, 8
!
!               -----------------------------------
!               EIGHT COLUMNS UPDATING TWO COLUMNS.
!               -----------------------------------
!
I1 = XPNT(K+1) - MM
I2 = XPNT(K+2) - MM
I3 = XPNT(K+3) - MM
I4 = XPNT(K+4) - MM
I5 = XPNT(K+5) - MM
I6 = XPNT(K+6) - MM
I7 = XPNT(K+7) - MM
I8 = XPNT(K+8) - MM
A1 = X(I1)
A2 = X(I2)
A3 = X(I3)
A4 = X(I4)
A5 = X(I5)
A6 = X(I6)
A7 = X(I7)
A8 = X(I8)
A9  = X(I1+1)
A10 = X(I2+1)
A11 = X(I3+1)
A12 = X(I4+1)
A13 = X(I5+1)
A14 = X(I6+1)
A15 = X(I7+1)
A16 = X(I8+1)
!
Y(IYBEG1+1) =  Y(IYBEG1+1) - &
& A1*A9 - A2*A10 - A3*A11 - A4*A12 - A5*A13 - &
& A6*A14 - A7*A15 - A8*A16
!
Y(IYBEG2+1) =  Y(IYBEG2+1) - &
& A9*A9 - A10*A10 - A11*A11 - A12*A12 - A13*A13 - &
& A14*A14 - A15*A15 - A16*A16
!
do  I = 2, MM-1
Y1 = Y(IYBEG1+I)
B1 = X(I1+I)
Y1 =  Y1 - B1 * A1
Y2 = Y(IYBEG2+I)
B2 = X(I2+I)
Y2 =  Y2 - B1 * A9
Y1 =  Y1 - B2 * A2
B3 = X(I3+I)
Y2 =  Y2 - B2 * A10
Y1 =  Y1 - B3 * A3
B4 = X(I4+I)
Y2 =  Y2 - B3 * A11
Y1 =  Y1 - B4 * A4
B5 = X(I5+I)
Y2 =  Y2 - B4 * A12
Y1 =  Y1 - B5 * A5
B6 = X(I6+I)
Y2 =  Y2 - B5 * A13
Y1 =  Y1 - B6 * A6
B7 = X(I7+I)
Y2 =  Y2 - B6 * A14
Y1 = Y1 - B7 * A7
B8 = X(I8+I)
Y2 = Y2 - B7 * A15
Y1 = Y1 - B8 * A8
Y(IYBEG1+I) = Y1
Y2 =  Y2 - B8 * A16
Y(IYBEG2+I) = Y2
end do
!
end do
!
!           -----------------------------
!           BOUNDARY CODE FOR THE K LOOP.
!           -----------------------------
!
SELECT CASE( N-K+2 )
CASE(1)
GO TO 2000
CASE(2)
GO TO 1700
CASE(3)
GO TO 1500
CASE(4)
GO TO 1300
CASE(5)
GO TO 1100
CASE(6)
GO TO 900
CASE(7)
GO TO 700
CASE(8)
GO TO 500
CASE DEFAULT
GO TO 2000
END SELECT
!
500 CONTINUE
!
!               -----------------------------------
!               SEVEN COLUMNS UPDATING TWO COLUMNS.
!               -----------------------------------
!
I1 = XPNT(K+1) - MM
I2 = XPNT(K+2) - MM
I3 = XPNT(K+3) - MM
I4 = XPNT(K+4) - MM
I5 = XPNT(K+5) - MM
I6 = XPNT(K+6) - MM
I7 = XPNT(K+7) - MM
A1 = X(I1)
A2 = X(I2)
A3 = X(I3)
A4 = X(I4)
A5 = X(I5)
A6 = X(I6)
A7 = X(I7)
A9  = X(I1+1)
A10 = X(I2+1)
A11 = X(I3+1)
A12 = X(I4+1)
A13 = X(I5+1)
A14 = X(I6+1)
A15 = X(I7+1)
!
Y(IYBEG1+1) =  Y(IYBEG1+1) - &
& A1*A9 - A2*A10 - A3*A11 - A4*A12 - A5*A13 - &
& A6*A14 - A7*A15
!
Y(IYBEG2+1) =  Y(IYBEG2+1) - &
& A9*A9 - A10*A10 - A11*A11 - A12*A12 - A13*A13 - &
& A14*A14 - A15*A15
!
do  I = 2, MM-1
Y1 = Y(IYBEG1+I)
B1 = X(I1+I)
Y1 =  Y1 - B1 * A1
Y2 = Y(IYBEG2+I)
B2 = X(I2+I)
Y2 =  Y2 - B1 * A9
Y1 =  Y1 - B2 * A2
B3 = X(I3+I)
Y2 =  Y2 - B2 * A10
Y1 =  Y1 - B3 * A3
B4 = X(I4+I)
Y2 =  Y2 - B3 * A11
Y1 =  Y1 - B4 * A4
B5 = X(I5+I)
Y2 =  Y2 - B4 * A12
Y1 =  Y1 - B5 * A5
B6 = X(I6+I)
Y2 =  Y2 - B5 * A13
Y1 =  Y1 - B6 * A6
B7 = X(I7+I)
Y2 =  Y2 - B6 * A14
Y1 = Y1 - B7 * A7
Y(IYBEG1+I) = Y1
Y2 = Y2 - B7 * A15
Y(IYBEG2+I) = Y2
end do
!
GO TO 2000
!
700 CONTINUE
!
!               ---------------------------------
!               SIX COLUMNS UPDATING TWO COLUMNS.
!               ---------------------------------
!
I1 = XPNT(K+1) - MM
I2 = XPNT(K+2) - MM
I3 = XPNT(K+3) - MM
I4 = XPNT(K+4) - MM
I5 = XPNT(K+5) - MM
I6 = XPNT(K+6) - MM
A1 = X(I1)
A2 = X(I2)
A3 = X(I3)
A4 = X(I4)
A5 = X(I5)
A6 = X(I6)
A9  = X(I1+1)
A10 = X(I2+1)
A11 = X(I3+1)
A12 = X(I4+1)
A13 = X(I5+1)
A14 = X(I6+1)
!
Y(IYBEG1+1) =  Y(IYBEG1+1) - &
& A1*A9 - A2*A10 - A3*A11 - A4*A12 - A5*A13 - &
& A6*A14
!
Y(IYBEG2+1) =  Y(IYBEG2+1) - &
& A9*A9 - A10*A10 - A11*A11 - A12*A12 - A13*A13 - &
& A14*A14
!
do  I = 2, MM-1
Y1 = Y(IYBEG1+I)
B1 = X(I1+I)
Y1 =  Y1 - B1 * A1
Y2 = Y(IYBEG2+I)
B2 = X(I2+I)
Y2 =  Y2 - B1 * A9
Y1 =  Y1 - B2 * A2
B3 = X(I3+I)
Y2 =  Y2 - B2 * A10
Y1 =  Y1 - B3 * A3
B4 = X(I4+I)
Y2 =  Y2 - B3 * A11
Y1 =  Y1 - B4 * A4
B5 = X(I5+I)
Y2 =  Y2 - B4 * A12
Y1 =  Y1 - B5 * A5
B6 = X(I6+I)
Y2 =  Y2 - B5 * A13
Y1 =  Y1 - B6 * A6
Y(IYBEG1+I) = Y1
Y2 =  Y2 - B6 * A14
Y(IYBEG2+I) = Y2
end do
!
GO TO 2000
!
900 CONTINUE
!
!               ----------------------------------
!               FIVE COLUMNS UPDATING TWO COLUMNS.
!               ----------------------------------
!
I1 = XPNT(K+1) - MM
I2 = XPNT(K+2) - MM
I3 = XPNT(K+3) - MM
I4 = XPNT(K+4) - MM
I5 = XPNT(K+5) - MM
A1 = X(I1)
A2 = X(I2)
A3 = X(I3)
A4 = X(I4)
A5 = X(I5)
A9  = X(I1+1)
A10 = X(I2+1)
A11 = X(I3+1)
A12 = X(I4+1)
A13 = X(I5+1)
!
Y(IYBEG1+1) =  Y(IYBEG1+1) - &
& A1*A9 - A2*A10 - A3*A11 - A4*A12 - A5*A13
!
Y(IYBEG2+1) =  Y(IYBEG2+1) - &
& A9*A9 - A10*A10 - A11*A11 - A12*A12 - A13*A13
!
DO  1000  I = 2, MM-1
Y1 = Y(IYBEG1+I)
B1 = X(I1+I)
Y1 =  Y1 - B1 * A1
Y2 = Y(IYBEG2+I)
B2 = X(I2+I)
Y2 =  Y2 - B1 * A9
Y1 =  Y1 - B2 * A2
B3 = X(I3+I)
Y2 =  Y2 - B2 * A10
Y1 =  Y1 - B3 * A3
B4 = X(I4+I)
Y2 =  Y2 - B3 * A11
Y1 =  Y1 - B4 * A4
B5 = X(I5+I)
Y2 =  Y2 - B4 * A12
Y1 =  Y1 - B5 * A5
Y(IYBEG1+I) = Y1
Y2 =  Y2 - B5 * A13
Y(IYBEG2+I) = Y2
1000 CONTINUE
!
GO TO 2000
!
1100 CONTINUE
!
!               ----------------------------------
!               FOUR COLUMNS UPDATING TWO COLUMNS.
!               ----------------------------------
!
I1 = XPNT(K+1) - MM
I2 = XPNT(K+2) - MM
I3 = XPNT(K+3) - MM
I4 = XPNT(K+4) - MM
A1 = X(I1)
A2 = X(I2)
A3 = X(I3)
A4 = X(I4)
A9  = X(I1+1)
A10 = X(I2+1)
A11 = X(I3+1)
A12 = X(I4+1)
!
Y(IYBEG1+1) =  Y(IYBEG1+1) - &
& A1*A9 - A2*A10 - A3*A11 - A4*A12
!
Y(IYBEG2+1) =  Y(IYBEG2+1) - &
& A9*A9 - A10*A10 - A11*A11 - A12*A12
!
DO  1200  I = 2, MM-1
Y1 = Y(IYBEG1+I)
B1 = X(I1+I)
Y1 =  Y1 - B1 * A1
Y2 = Y(IYBEG2+I)
B2 = X(I2+I)
Y2 =  Y2 - B1 * A9
Y1 =  Y1 - B2 * A2
B3 = X(I3+I)
Y2 =  Y2 - B2 * A10
Y1 =  Y1 - B3 * A3
B4 = X(I4+I)
Y2 =  Y2 - B3 * A11
Y1 =  Y1 - B4 * A4
Y(IYBEG1+I) = Y1
Y2 =  Y2 - B4 * A12
Y(IYBEG2+I) = Y2
1200 CONTINUE
!
GO TO 2000
!
1300 CONTINUE
!
!               -----------------------------------
!               THREE COLUMNS UPDATING TWO COLUMNS.
!               -----------------------------------
!
I1 = XPNT(K+1) - MM
I2 = XPNT(K+2) - MM
I3 = XPNT(K+3) - MM
A1 = X(I1)
A2 = X(I2)
A3 = X(I3)
A9  = X(I1+1)
A10 = X(I2+1)
A11 = X(I3+1)
!
Y(IYBEG1+1) =  Y(IYBEG1+1) - &
& A1*A9 - A2*A10 - A3*A11
!
Y(IYBEG2+1) =  Y(IYBEG2+1) - &
& A9*A9 - A10*A10 - A11*A11
!
DO  1400  I = 2, MM-1
Y1 = Y(IYBEG1+I)
B1 = X(I1+I)
Y1 =  Y1 - B1 * A1
Y2 = Y(IYBEG2+I)
B2 = X(I2+I)
Y2 =  Y2 - B1 * A9
Y1 =  Y1 - B2 * A2
B3 = X(I3+I)
Y2 =  Y2 - B2 * A10
Y1 =  Y1 - B3 * A3
Y(IYBEG1+I) = Y1
Y2 =  Y2 - B3 * A11
Y(IYBEG2+I) = Y2
1400 CONTINUE
!
GO TO 2000
!
1500 CONTINUE
!
!     ---------------------------------
!     TWO COLUMNS UPDATING TWO COLUMNS.
!     ---------------------------------
!
I1 = XPNT(K+1) - MM
I2 = XPNT(K+2) - MM
A1 = X(I1)
A2 = X(I2)
A9  = X(I1+1)
A10 = X(I2+1)
!
Y(IYBEG1+1) =  Y(IYBEG1+1) - &
& A1*A9 - A2*A10
!
Y(IYBEG2+1) =  Y(IYBEG2+1) - &
& A9*A9 - A10*A10
!
do  I = 2, MM-1
Y1 = Y(IYBEG1+I)
B1 = X(I1+I)
Y1 =  Y1 - B1 * A1
Y2 = Y(IYBEG2+I)
B2 = X(I2+I)
Y2 =  Y2 - B1 * A9
Y1 =  Y1 - B2 * A2
Y(IYBEG1+I) = Y1
Y2 =  Y2 - B2 * A10
Y(IYBEG2+I) = Y2
end do
!
GO TO 2000
!
1700 CONTINUE
!
!               --------------------------------
!               ONE COLUMN UPDATING TWO COLUMNS.
!               --------------------------------
!
I1 = XPNT(K+1) - MM
A1 = X(I1)
A9  = X(I1+1)
!
Y(IYBEG1+1) =  Y(IYBEG1+1) - &
& A1*A9
!
Y(IYBEG2+1) =  Y(IYBEG2+1) - &
& A9*A9
!
do  I = 2, MM-1
Y1 = Y(IYBEG1+I)
B1 = X(I1+I)
Y1 =  Y1 - B1 * A1
Y2 = Y(IYBEG2+I)
Y(IYBEG1+I) = Y1
Y2 =  Y2 - B1 * A9
Y(IYBEG2+I) = Y2
end do
!
GO TO 2000
!
!           -----------------------------------------------
!           PREPARE FOR NEXT PAIR OF COLUMNS TO BE UPDATED.
!           -----------------------------------------------
!
2000 CONTINUE
MM = MM - 2
IYBEG = IYBEG2 + LENY + 1
LENY = LENY - 2
!
end do
!
!       -----------------------------------------------------
!       BOUNDARY CODE FOR J LOOP:  EXECUTED WHENVER Q IS ODD.
!       -----------------------------------------------------
!
IF  ( J .EQ. QQ )  THEN
CALL  SMXPY8  ( MM, N, Y(IYBEG), XPNT, X )
ENDIF
!
RETURN
END
!***********************************************************************
!***********************************************************************
!
!   Version:        0.4
!   Last modified:  December 27, 1994
!   Authors:        Esmond G. Ng and Barry W. Peyton
!
!   Mathematical Sciences Section, Oak Ridge National Laboratory
!
!***********************************************************************
!***********************************************************************
!*************     MMPYI  .... MATRIX-MATRIX MULTIPLY     **************
!***********************************************************************
!***********************************************************************
!
!   PURPOSE -
!       THIS ROUTINE PERFORMS A MATRIX-MATRIX MULTIPLY, Y = Y + XA,
!       ASSUMING DATA STRUCTURES USED IN SOME OF OUR SPARSE CHOLESKY
!       CODES.
!
!       MATRIX X HAS ONLY 1 COLUMN.
!
!   INPUT PARAMETERS -
!       M               -   NUMBER OF ROWS IN X AND IN Y.
!       Q               -   NUMBER OF COLUMNS IN A AND Y.
!       XPNT(*)         -   XPNT(J+1) POINTS ONE LOCATION BEYOND THE
!                           END OF THE J-TH COLUMN OF X.  XPNT IS ALSO
!                           USED TO ACCESS THE ROWS OF A.
!       X(*)            -   CONTAINS THE COLUMNS OF X AND THE ROWS OF A.
!       IY(*)           -   IY(COL) POINTS TO THE BEGINNING OF COLUMN
!       RELIND(*)       -   RELATIVE INDICES.
!
!   UPDATED PARAMETERS -
!       Y(*)            -   ON OUTPUT, Y = Y + AX.
!
!***********************************************************************
!
SUBROUTINE  MMPYI  (  M     , Q     , XPNT  , X     , IY    , &
& Y     , RELIND                          )
!
!***********************************************************************
!
!       -----------
!       PARAMETERS.
!       -----------
!
INTEGER             M     , Q
INTEGER             IY(*)         , RELIND(*)     , &
& XPNT(*)
DOUBLE PRECISION    X(*)      , Y(*)
!
!       ----------------
!       LOCAL VARIABLES.
!       ----------------
!
INTEGER             COL   , I     , ISUB  , K     , YLAST
DOUBLE PRECISION    A
!
!***********************************************************************
!
DO  200  K = 1, Q
COL = XPNT(K)
YLAST = IY(COL+1) - 1
A = - X(K)
!DIR$   IVDEP
DO  100  I = K, M
ISUB = XPNT(I)
ISUB = YLAST - RELIND(ISUB)
Y(ISUB) = Y(ISUB) + A*X(I)
100 CONTINUE
200 CONTINUE
RETURN
!
END
!***********************************************************************
!***********************************************************************
!
!   Version:        0.3
!   Last modified:  December 27, 1994
!   Authors:        Esmond G. Ng and Barry W. Peyton
!
!   Mathematical Sciences Section, Oak Ridge National Laboratoy
!
!***********************************************************************
!***********************************************************************
!******     PCHOL .... DENSE PARTIAL CHOLESKY             **************
!***********************************************************************
!***********************************************************************
!
!     PURPOSE - THIS ROUTINE PERFORMS CHOLESKY
!               FACTORIZATION ON THE COLUMNS OF A SUPERNODE
!               THAT HAVE RECEIVED ALL UPDATES FROM COLUMNS
!               EXTERNAL TO THE SUPERNODE.
!
!     INPUT PARAMETERS -
!        M      - NUMBER OF ROWS (LENGTH OF THE FIRST COLUMN).
!        N      - NUMBER OF COLUMNS IN THE SUPERNODE.
!        XPNT   - XPNT(J+1) POINTS ONE LOCATION BEYOND THE END
!                 OF THE J-TH COLUMN OF THE SUPERNODE.
!        X(*)   - CONTAINS THE COLUMNS OF OF THE SUPERNODE TO
!                 BE FACTORED.
!
!     EXTERNAL ROUTINE:
!        SMXPY8 -  MATRIX-VECTOR MULTIPLY WITH 8 LOOP UNROLLING.
!
!     OUTPUT PARAMETERS -
!        X(*)   - ON OUTPUT, CONTAINS THE FACTORED COLUMNS OF
!                 THE SUPERNODE.
!
!***********************************************************************
!
SUBROUTINE  PCHOL  ( M, N, XPNT, X, MXDIAG, NTINY )
!
!***********************************************************************
!
!     -----------
!     PARAMETERS.
!     -----------
!
EXTERNAL            SMXPY8
!
INTEGER             M, N
!
INTEGER             XPNT(*)
!
!xPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPC
DOUBLE PRECISION    X(*), MXDIAG
INTEGER             NTINY
!
!     ----------------
!     LOCAL VARIABLES.
!     ----------------
!
INTEGER             JPNT  , JCOL  , MM
!
DOUBLE PRECISION    DIAG
!
!***********************************************************************
!
!       ------------------------------------------
!       FOR EVERY COLUMN JCOL IN THE SUPERNODE ...
!       ------------------------------------------
MM     = M
JPNT = XPNT(1)
DO  100  JCOL = 1, N
!
!           ----------------------------------
!           UPDATE JCOL WITH PREVIOUS COLUMNS.
!           ----------------------------------
IF  ( JCOL .GT. 1 )  THEN
CALL SMXPY8 ( MM, JCOL-1, X(JPNT), XPNT, X )
ENDIF
!
!           ---------------------------
!           COMPUTE THE DIAGONAL ENTRY.
!           ---------------------------
DIAG = X(JPNT)
!xPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPC
IF (DIAG .LE. 1.0D-30*MXDIAG) THEN
DIAG = 1.0D+128
NTINY = NTINY+1
ENDIF
!xPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPC

DIAG = SQRT ( DIAG )
X(JPNT) = DIAG
DIAG = 1.0D+00 / DIAG
!
!           ----------------------------------------------------
!           SCALE COLUMN JCOL WITH RECIPROCAL OF DIAGONAL ENTRY.
!           ----------------------------------------------------
MM = MM - 1
JPNT = JPNT + 1
CALL DSCAL1 ( MM, DIAG, X(JPNT) )
JPNT = JPNT + MM
!
100 CONTINUE
!
RETURN
END
!***********************************************************************
!***********************************************************************
!
!   Version:        0.4
!   Last modified:  January 12, 1995
!   Authors:        Esmond G. Ng and Barry W. Peyton
!
!   Mathematical Sciences Section, Oak Ridge National Laboratory
!
!***********************************************************************
!***********************************************************************
!**************    SFINIT  ..... SET UP FOR SYMB. FACT.     ************
!***********************************************************************
!***********************************************************************
!
!   PURPOSE:
!       THIS SUBROUTINE COMPUTES THE STORAGE REQUIREMENTS AND SETS UP
!       PRELIMINARY DATA STRUCTURES FOR THE SYMBOLIC FACTORIZATION.
!
!   NOTE:
!       THIS VERSION PRODUCES THE MAXIMAL SUPERNODE PARTITION (I.E.,
!       THE ONE WITH THE FEWEST POSSIBLE SUPERNODES).
!
!   INPUT PARAMETERS:
!       NEQNS       -   NUMBER OF EQUATIONS.
!       NNZA        -   LENGTH OF ADJACENCY STRUCTURE.
!       XADJ(*)     -   ARRAY OF LENGTH NEQNS+1, CONTAINING POINTERS
!                       TO THE ADJACENCY STRUCTURE.
!       ADJNCY(*)   -   ARRAY OF LENGTH XADJ(NEQNS+1)-1, CONTAINING
!                       THE ADJACENCY STRUCTURE.
!       PERM(*)     -   ARRAY OF LENGTH NEQNS, CONTAINING THE
!                       POSTORDERING.
!       INVP(*)     -   ARRAY OF LENGTH NEQNS, CONTAINING THE
!                       INVERSE OF THE POSTORDERING.
!       IWSIZ       -   SIZE OF INTEGER WORKING STORAGE.
!
!   OUTPUT PARAMETERS:
!       COLCNT(*)   -   ARRAY OF LENGTH NEQNS, CONTAINING THE NUMBER
!                       OF NONZEROS IN EACH COLUMN OF THE FACTOR,
!                       INCLUDING THE DIAGONAL ENTRY.
!       NNZL        -   NUMBER OF NONZEROS IN THE FACTOR, INCLUDING
!                       THE DIAGONAL ENTRIES.
!       NSUB        -   NUMBER OF SUBSCRIPTS.
!       NSUPER      -   NUMBER OF SUPERNODES (<= NEQNS).
!       SNODE(*)    -   ARRAY OF LENGTH NEQNS FOR RECORDING
!                       SUPERNODE MEMBERSHIP.
!       XSUPER(*)   -   ARRAY OF LENGTH NEQNS+1, CONTAINING THE
!                       SUPERNODE PARTITIONING.
!       IFLAG(*)    -   ERROR FLAG.
!                          0: SUCCESSFUL SF INITIALIZATION.
!                         -1: INSUFFICENT WORKING STORAGE
!                             [IWORK(*)].
!
!   WORK PARAMETERS:
!       IWORK(*)    -   INTEGER WORK ARRAY OF LENGTH 7*NEQNS+3.
!
!   FIRST CREATED ON    NOVEMEBER 14, 1994.
!   LAST UPDATED ON     January 12, 1995.
!
!***********************************************************************
!
SUBROUTINE  SFINIT (  NEQNS , NNZA  , XADJ  , ADJNCY, PERM  , &
& INVP  , COLCNT, NNZL  , NSUB  , NSUPER, &
& SNODE , XSUPER, IWSIZ , IWORK , IFLAG   )
!
!       -----------
!       PARAMETERS.
!       -----------
INTEGER             IFLAG , IWSIZ , NNZA  , NEQNS , NNZL  , &
& NSUB  , NSUPER
INTEGER             ADJNCY(NNZA)    , COLCNT(NEQNS)   , &
& INVP(NEQNS)     , IWORK(7*NEQNS+3), &
& PERM(NEQNS)     , SNODE(NEQNS)    , &
& XADJ(NEQNS+1)   , XSUPER(NEQNS+1)
!
!***********************************************************************
!
!       --------------------------------------------------------
!       RETURN IF THERE IS INSUFFICIENT INTEGER WORKING STORAGE.
!       --------------------------------------------------------
IFLAG = 0
IF  ( IWSIZ .LT. 7*NEQNS+3 )  THEN
IFLAG = -1
RETURN
ENDIF
!
!       ------------------------------------------
!       COMPUTE ELIMINATION TREE AND POSTORDERING.
!       ------------------------------------------
CALL  ETORDR (  NEQNS , XADJ  , ADJNCY, PERM  , INVP  , &
& IWORK(1)              , &
& IWORK(NEQNS+1)        , &
& IWORK(2*NEQNS+1)      , &
& IWORK(3*NEQNS+1)        )
!
!       ---------------------------------------------
!       COMPUTE ROW AND COLUMN FACTOR NONZERO COUNTS.
!       ---------------------------------------------
CALL  FCNTHN (  NEQNS , NNZA  , XADJ  , ADJNCY, PERM  , &
& INVP  , IWORK(1)      , SNODE , COLCNT, &
& NNZL  , &
& IWORK(NEQNS+1)        , &
& IWORK(2*NEQNS+1)      , &
& XSUPER                , &
& IWORK(3*NEQNS+1)      , &
& IWORK(4*NEQNS+2)      , &
& IWORK(5*NEQNS+3)      , &
& IWORK(6*NEQNS+4)        )
!
!       ---------------------------------------------------------
!       REARRANGE CHILDREN SO THAT THE LAST CHILD HAS THE MAXIMUM
!       NUMBER OF NONZEROS IN ITS COLUMN OF L.
!       ---------------------------------------------------------
CALL  CHORDR (  NEQNS , PERM  , INVP  , &
& COLCNT, &
& IWORK(1)              , &
& IWORK(NEQNS+1)        , &
& IWORK(2*NEQNS+1)      , &
& IWORK(3*NEQNS+1)        )
!
!       ----------------
!       FIND SUPERNODES.
!       ----------------
CALL  FSUP1  (  NEQNS , IWORK(1)      , COLCNT, NSUB  , &
& NSUPER, SNODE                           )
CALL  FSUP2  (  NEQNS , NSUPER,  SNODE,     XSUPER      )
!
RETURN
END
!***********************************************************************
!***********************************************************************
!
!   Version:        0.4
!   Last modified:  December 27, 1994
!   Authors:        Esmond G. Ng and Barry W. Peyton
!
!   Mathematical Sciences Section, Oak Ridge National Laboratory
!
!***********************************************************************
!***********************************************************************
!******     SMXPY8 .... MATRIX-VECTOR MULTIPLY            **************
!***********************************************************************
!***********************************************************************
!
!     PURPOSE - THIS ROUTINE PERFORMS A MATRIX-VECTOR MULTIPLY,
!               Y = Y + AX, ASSUMING DATA STRUCTURES USED IN
!               RECENTLY DEVELOPED SPARSE CHOLESKY CODES.  THE
!               '8' SIGNIFIES LEVEL 8 LOOP UNROLLING.
!
!     INPUT PARAMETERS -
!        M      - NUMBER OF ROWS.
!        N      - NUMBER OF COLUMNS.
!        Y      - M-VECTOR TO WHICH AX WILL BE ADDED.
!        APNT   - INDEX VECTOR FOR A.  APNT(I) POINTS TO THE
!                 FIRST NONZERO IN COLUMN I OF A.
!        Y      - ON OUTPUT, CONTAINS Y = Y + AX.
!
!***********************************************************************
!
SUBROUTINE  SMXPY8 ( M, N, Y, APNT, A )
!
!***********************************************************************
!
!     -----------
!     PARAMETERS.
!     -----------
!
INTEGER             M, N, LEVEL
!
INTEGER             APNT(*)
!
DOUBLE PRECISION    Y(*), A(*)
!
PARAMETER           ( LEVEL = 8 )
!
!     ----------------
!     LOCAL VARIABLES.
!     ----------------
!
INTEGER             I, I1, I2, I3, I4, I5, I6, I7, I8, &
& J, REMAIN
!
DOUBLE PRECISION    A1, A2, A3, A4, A5, A6, A7, A8
!
!***********************************************************************
!
REMAIN = MOD ( N, LEVEL )
!
SELECT CASE ( REMAIN+1 )
CASE(1)
GO TO 2000
CASE(2)
GO TO 100
CASE(3)
GO TO 200
CASE(4)
GO TO 300
CASE(5)
GO TO 400
CASE(6)
GO TO 500
CASE(7)
GO TO 600
CASE(8)
GO TO 700
CASE DEFAULT
GO TO 2000
END SELECT
!
100 CONTINUE
I1 = APNT(1+1) - M
A1 = - A(I1)
DO  150  I = 1, M
Y(I) = Y(I) + A1*A(I1)
I1 = I1 + 1
150 CONTINUE
GO TO 2000
!
200 CONTINUE
I1 = APNT(1+1) - M
I2 = APNT(1+2) - M
A1 = - A(I1)
A2 = - A(I2)
DO  250  I = 1, M
Y(I) = ( (Y(I)) &
& + A1*A(I1)) + A2*A(I2)
I1 = I1 + 1
I2 = I2 + 1
250 CONTINUE
GO TO 2000
!
300 CONTINUE
I1 = APNT(1+1) - M
I2 = APNT(1+2) - M
I3 = APNT(1+3) - M
A1 = - A(I1)
A2 = - A(I2)
A3 = - A(I3)
DO  350  I = 1, M
Y(I) = (( (Y(I)) &
& + A1*A(I1)) + A2*A(I2)) &
& + A3*A(I3)
I1 = I1 + 1
I2 = I2 + 1
I3 = I3 + 1
350 CONTINUE
GO TO 2000
!
400 CONTINUE
I1 = APNT(1+1) - M
I2 = APNT(1+2) - M
I3 = APNT(1+3) - M
I4 = APNT(1+4) - M
A1 = - A(I1)
A2 = - A(I2)
A3 = - A(I3)
A4 = - A(I4)
DO  450  I = 1, M
Y(I) = ((( (Y(I)) &
& + A1*A(I1)) + A2*A(I2)) &
& + A3*A(I3)) + A4*A(I4)
I1 = I1 + 1
I2 = I2 + 1
I3 = I3 + 1
I4 = I4 + 1
450 CONTINUE
GO TO 2000
!
500 CONTINUE
I1 = APNT(1+1) - M
I2 = APNT(1+2) - M
I3 = APNT(1+3) - M
I4 = APNT(1+4) - M
I5 = APNT(1+5) - M
A1 = - A(I1)
A2 = - A(I2)
A3 = - A(I3)
A4 = - A(I4)
A5 = - A(I5)
DO  550  I = 1, M
Y(I) = (((( (Y(I)) &
& + A1*A(I1)) + A2*A(I2)) &
& + A3*A(I3)) + A4*A(I4)) &
& + A5*A(I5)
I1 = I1 + 1
I2 = I2 + 1
I3 = I3 + 1
I4 = I4 + 1
I5 = I5 + 1
550 CONTINUE
GO TO 2000
!
600 CONTINUE
I1 = APNT(1+1) - M
I2 = APNT(1+2) - M
I3 = APNT(1+3) - M
I4 = APNT(1+4) - M
I5 = APNT(1+5) - M
I6 = APNT(1+6) - M
A1 = - A(I1)
A2 = - A(I2)
A3 = - A(I3)
A4 = - A(I4)
A5 = - A(I5)
A6 = - A(I6)
DO  650  I = 1, M
Y(I) = ((((( (Y(I)) &
& + A1*A(I1)) + A2*A(I2)) &
& + A3*A(I3)) + A4*A(I4)) &
& + A5*A(I5)) + A6*A(I6)
I1 = I1 + 1
I2 = I2 + 1
I3 = I3 + 1
I4 = I4 + 1
I5 = I5 + 1
I6 = I6 + 1
650 CONTINUE
GO TO 2000
!
700 CONTINUE
I1 = APNT(1+1) - M
I2 = APNT(1+2) - M
I3 = APNT(1+3) - M
I4 = APNT(1+4) - M
I5 = APNT(1+5) - M
I6 = APNT(1+6) - M
I7 = APNT(1+7) - M
A1 = - A(I1)
A2 = - A(I2)
A3 = - A(I3)
A4 = - A(I4)
A5 = - A(I5)
A6 = - A(I6)
A7 = - A(I7)
DO  750  I = 1, M
Y(I) = (((((( (Y(I)) &
& + A1*A(I1)) + A2*A(I2)) &
& + A3*A(I3)) + A4*A(I4)) &
& + A5*A(I5)) + A6*A(I6)) &
& + A7*A(I7)
I1 = I1 + 1
I2 = I2 + 1
I3 = I3 + 1
I4 = I4 + 1
I5 = I5 + 1
I6 = I6 + 1
I7 = I7 + 1
750 CONTINUE
GO TO 2000
!
2000 CONTINUE
DO  4000  J = REMAIN+1, N, LEVEL
I1 = APNT(J+1) - M
I2 = APNT(J+2) - M
I3 = APNT(J+3) - M
I4 = APNT(J+4) - M
I5 = APNT(J+5) - M
I6 = APNT(J+6) - M
I7 = APNT(J+7) - M
I8 = APNT(J+8) - M
A1 = - A(I1)
A2 = - A(I2)
A3 = - A(I3)
A4 = - A(I4)
A5 = - A(I5)
A6 = - A(I6)
A7 = - A(I7)
A8 = - A(I8)
DO  3000  I = 1, M
Y(I) = ((((((( (Y(I)) &
& + A1*A(I1)) + A2*A(I2)) &
& + A3*A(I3)) + A4*A(I4)) &
& + A5*A(I5)) + A6*A(I6)) &
& + A7*A(I7)) + A8*A(I8)
I1 = I1 + 1
I2 = I2 + 1
I3 = I3 + 1
I4 = I4 + 1
I5 = I5 + 1
I6 = I6 + 1
I7 = I7 + 1
I8 = I8 + 1
3000 CONTINUE
4000 CONTINUE
!
RETURN
END
!***********************************************************************
!***********************************************************************
!
!   Version:        0.4
!   Last modified:  February 13, 1995
!   Authors:        Esmond G. Ng and Barry W. Peyton
!
!   Mathematical Sciences Section, Oak Ridge National Laboratory
!
!***********************************************************************
!***********************************************************************
!*************     SYMFC2 ..... SYMBOLIC FACTORIZATION    **************
!***********************************************************************
!***********************************************************************
!
!   PURPOSE:
!       THIS ROUTINE PERFORMS SUPERNODAL SYMBOLIC FACTORIZATION ON A
!       REORDERED LINEAR SYSTEM.  IT ASSUMES ACCESS TO THE COLUMNS
!       COUNTS, SUPERNODE PARTITION, AND SUPERNODAL ELIMINATION TREE
!       ASSOCIATED WITH THE FACTOR MATRIX L.
!
!   INPUT PARAMETERS:
!       (I) NEQNS       -   NUMBER OF EQUATIONS
!       (I) ADJLEN      -   LENGTH OF THE ADJACENCY LIST.
!       (I) XADJ(*)     -   ARRAY OF LENGTH NEQNS+1 CONTAINING POINTERS
!                           TO THE ADJACENCY STRUCTURE.
!       (I) ADJNCY(*)   -   ARRAY OF LENGTH XADJ(NEQNS+1)-1 CONTAINING
!                           THE ADJACENCY STRUCTURE.
!       (I) PERM(*)     -   ARRAY OF LENGTH NEQNS CONTAINING THE
!                           POSTORDERING.
!       (I) INVP(*)     -   ARRAY OF LENGTH NEQNS CONTAINING THE
!                           INVERSE OF THE POSTORDERING.
!       (I) COLCNT(*)   -   ARRAY OF LENGTH NEQNS, CONTAINING THE NUMBER
!                           OF NONZEROS IN EACH COLUMN OF THE FACTOR,
!                           INCLUDING THE DIAGONAL ENTRY.
!       (I) NSUPER      -   NUMBER OF SUPERNODES.
!       (I) XSUPER(*)   -   ARRAY OF LENGTH NSUPER+1, CONTAINING THE
!                           FIRST COLUMN OF EACH SUPERNODE.
!       (I) SNODE(*)    -   ARRAY OF LENGTH NEQNS FOR RECORDING
!                           SUPERNODE MEMBERSHIP.
!       (I) NOFSUB      -   NUMBER OF SUBSCRIPTS TO BE STORED IN
!                           LINDX(*).
!
!   OUTPUT PARAMETERS:
!       (I) XLINDX      -   ARRAY OF LENGTH NEQNS+1, CONTAINING POINTERS
!                           INTO THE SUBSCRIPT VECTOR.
!       (I) LINDX       -   ARRAY OF LENGTH MAXSUB, CONTAINING THE
!                           COMPRESSED SUBSCRIPTS.
!       (I) XLNZ        -   COLUMN POINTERS FOR L.
!       (I) FLAG        -   ERROR FLAG:
!                               0 - NO ERROR.
!                               1 - INCONSISTANCY IN THE INPUT.
!
!   WORKING PARAMETERS:
!       (I) MRGLNK      -   ARRAY OF LENGTH NSUPER, CONTAINING THE
!                           CHILDREN OF EACH SUPERNODE AS A LINKED LIST.
!       (I) RCHLNK      -   ARRAY OF LENGTH NEQNS+1, CONTAINING THE
!                           CURRENT LINKED LIST OF MERGED INDICES (THE
!                           "REACH" SET).
!       (I) MARKER      -   ARRAY OF LENGTH NEQNS USED TO MARK INDICES
!                           AS THEY ARE INTRODUCED INTO EACH SUPERNODE'S
!                           INDEX SET.
!
!***********************************************************************
!
SUBROUTINE  SYMFC2 (  NEQNS , ADJLEN, XADJ  , ADJNCY, PERM  , &
& INVP  , COLCNT, NSUPER, XSUPER, SNODE , &
& NOFSUB, XLINDX, LINDX , XLNZ  , MRGLNK, &
& RCHLNK, MARKER, FLAG    )
!
!***********************************************************************
!
!       -----------
!       PARAMETERS.
!       -----------
INTEGER             ADJLEN, FLAG  , NEQNS , NOFSUB, NSUPER
INTEGER             ADJNCY(ADJLEN), COLCNT(NEQNS) , &
& INVP(NEQNS)   , MARKER(NEQNS) , &
& MRGLNK(NSUPER), LINDX(NOFSUB) , &
& PERM(NEQNS)   , RCHLNK(0:NEQNS), &
& SNODE(NEQNS)  , XSUPER(NSUPER+1)
INTEGER             XADJ(NEQNS+1) , XLINDX(NSUPER+1), &
& XLNZ(NEQNS+1)
!
!       ----------------
!       LOCAL VARIABLES.
!       ----------------
INTEGER             FSTCOL, HEAD  , I     , JNZBEG, JNZEND, &
& JPTR  , JSUP  , JWIDTH, KNZ   , KNZBEG, &
& KNZEND, KPTR  , KSUP  , LENGTH, LSTCOL, &
& NEWI  , NEXTI , NODE  , NZBEG , NZEND , &
& PCOL  , PSUP  , POINT , TAIL  , WIDTH
!
!***********************************************************************
!
FLAG = 0
IF  ( NEQNS .LE. 0 )  RETURN
!
!       ---------------------------------------------------
!       INITIALIZATIONS ...
!           NZEND  : POINTS TO THE LAST USED SLOT IN LINDX.
!           TAIL   : END OF LIST INDICATOR
!                    (IN RCHLNK(*), NOT MRGLNK(*)).
!           MRGLNK : CREATE EMPTY LISTS.
!           MARKER : "UNMARK" THE INDICES.
!       ---------------------------------------------------
NZEND = 0
HEAD = 0
TAIL = NEQNS + 1
POINT = 1
DO  50  I = 1, NEQNS
MARKER(I) = 0
XLNZ(I) = POINT
POINT = POINT + COLCNT(I)
50 CONTINUE
XLNZ(NEQNS+1) = POINT
POINT = 1
DO  100  KSUP = 1, NSUPER
MRGLNK(KSUP) = 0
FSTCOL = XSUPER(KSUP)
XLINDX(KSUP) = POINT
POINT = POINT + COLCNT(FSTCOL)
100 CONTINUE
XLINDX(NSUPER+1) = POINT
!
!       ---------------------------
!       FOR EACH SUPERNODE KSUP ...
!       ---------------------------
DO  1000  KSUP = 1, NSUPER
!
!           ---------------------------------------------------------
!           INITIALIZATIONS ...
!               FSTCOL : FIRST COLUMN OF SUPERNODE KSUP.
!               LSTCOL : LAST COLUMN OF SUPERNODE KSUP.
!               KNZ    : WILL COUNT THE NONZEROS OF L IN COLUMN KCOL.
!               RCHLNK : INITIALIZE EMPTY INDEX LIST FOR KCOL.
!           ---------------------------------------------------------
FSTCOL = XSUPER(KSUP)
LSTCOL = XSUPER(KSUP+1) - 1
WIDTH  = LSTCOL - FSTCOL + 1
LENGTH = COLCNT(FSTCOL)
KNZ = 0
RCHLNK(HEAD) = TAIL
JSUP = MRGLNK(KSUP)
!
!           -------------------------------------------------
!           IF KSUP HAS CHILDREN IN THE SUPERNODAL E-TREE ...
!           -------------------------------------------------
IF  ( JSUP .GT. 0 )  THEN
!               ---------------------------------------------
!               COPY THE INDICES OF THE FIRST CHILD JSUP INTO
!               THE LINKED LIST, AND MARK EACH WITH THE VALUE
!               KSUP.
!               ---------------------------------------------
JWIDTH = XSUPER(JSUP+1) - XSUPER(JSUP)
JNZBEG = XLINDX(JSUP) + JWIDTH
JNZEND = XLINDX(JSUP+1) - 1
DO  200  JPTR = JNZEND, JNZBEG, -1
NEWI = LINDX(JPTR)
KNZ = KNZ+1
MARKER(NEWI) = KSUP
RCHLNK(NEWI) = RCHLNK(HEAD)
RCHLNK(HEAD) = NEWI
200 CONTINUE
!               ------------------------------------------
!               FOR EACH SUBSEQUENT CHILD JSUP OF KSUP ...
!               ------------------------------------------
JSUP = MRGLNK(JSUP)
300 CONTINUE
IF  ( JSUP .NE. 0  .AND.  KNZ .LT. LENGTH )  THEN
!                   ----------------------------------------
!                   MERGE THE INDICES OF JSUP INTO THE LIST,
!                   AND MARK NEW INDICES WITH VALUE KSUP.
!                   ----------------------------------------
JWIDTH = XSUPER(JSUP+1) - XSUPER(JSUP)
JNZBEG = XLINDX(JSUP) + JWIDTH
JNZEND = XLINDX(JSUP+1) - 1
NEXTI = HEAD
DO  500  JPTR = JNZBEG, JNZEND
NEWI = LINDX(JPTR)
400 CONTINUE
I = NEXTI
NEXTI = RCHLNK(I)
IF  ( NEWI .GT. NEXTI )  GO TO 400
IF  ( NEWI .LT. NEXTI )  THEN
KNZ = KNZ+1
RCHLNK(I) = NEWI
RCHLNK(NEWI) = NEXTI
MARKER(NEWI) = KSUP
NEXTI = NEWI
ENDIF
500 CONTINUE
JSUP = MRGLNK(JSUP)
GO TO 300
ENDIF
ENDIF
!           ---------------------------------------------------
!           STRUCTURE OF A(*,FSTCOL) HAS NOT BEEN EXAMINED YET.
!           "SORT" ITS STRUCTURE INTO THE LINKED LIST,
!           INSERTING ONLY THOSE INDICES NOT ALREADY IN THE
!           LIST.
!           ---------------------------------------------------
IF  ( KNZ .LT. LENGTH )  THEN
NODE = PERM(FSTCOL)
KNZBEG = XADJ(NODE)
KNZEND = XADJ(NODE+1) - 1
DO  700  KPTR = KNZBEG, KNZEND
NEWI = ADJNCY(KPTR)
NEWI = INVP(NEWI)
IF  ( NEWI .GT. FSTCOL  .AND. &
& MARKER(NEWI) .NE. KSUP )  THEN
!                       --------------------------------
!                       POSITION AND INSERT NEWI IN LIST
!                       AND MARK IT WITH KCOL.
!                       --------------------------------
NEXTI = HEAD
600 CONTINUE
I = NEXTI
NEXTI = RCHLNK(I)
IF  ( NEWI .GT. NEXTI )  GO TO 600
KNZ = KNZ + 1
RCHLNK(I) = NEWI
RCHLNK(NEWI) = NEXTI
MARKER(NEWI) = KSUP
ENDIF
700 CONTINUE
ENDIF
!           ------------------------------------------------------------
!           IF KSUP HAS NO CHILDREN, INSERT FSTCOL INTO THE LINKED LIST.
!           ------------------------------------------------------------
IF  ( RCHLNK(HEAD) .NE. FSTCOL )  THEN
RCHLNK(FSTCOL) = RCHLNK(HEAD)
RCHLNK(HEAD) = FSTCOL
KNZ = KNZ + 1
ENDIF
!
!           --------------------------------------------
!           COPY INDICES FROM LINKED LIST INTO LINDX(*).
!           --------------------------------------------
NZBEG = NZEND + 1
NZEND = NZEND + KNZ
IF  ( NZEND+1 .NE. XLINDX(KSUP+1) )  GO TO 8000
I = HEAD
DO  800  KPTR = NZBEG, NZEND
I = RCHLNK(I)
LINDX(KPTR) = I
800 CONTINUE
!
!           ---------------------------------------------------
!           IF KSUP HAS A PARENT, INSERT KSUP INTO ITS PARENT'S
!           "MERGE" LIST.
!           ---------------------------------------------------
IF  ( LENGTH .GT. WIDTH )  THEN
PCOL = LINDX ( XLINDX(KSUP) + WIDTH )
PSUP = SNODE(PCOL)
MRGLNK(KSUP) = MRGLNK(PSUP)
MRGLNK(PSUP) = KSUP
ENDIF
!
1000 CONTINUE
!
RETURN
!
!       -----------------------------------------------
!       INCONSISTENCY IN DATA STRUCTURE WAS DISCOVERED.
!       -----------------------------------------------
8000 CONTINUE
FLAG = -2
RETURN
!
END













subroutine genrcm ( node_num, adj_num, adj_row, adj, perm )

!*****************************************************************************80
!
!! GENRCM finds the reverse Cuthill-Mckee ordering for a general graph.
!
!  Discussion:
!
!    For each connected component in the graph, the routine obtains
!    an ordering by calling RCM.
!
!  Modified:
!
!    04 January 2003
!
!  Author:
!
!    Alan George, Joseph Liu
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, integer ADJ_NUM, the number of adjacency entries.
!
!    Input, integer ADJ_ROW(NODE_NUM+1).  Information about row I is stored
!    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
!    Input, integer ADJ(ADJ_NUM), the adjacency structure.
!    For each row, it contains the column indices of the nonzero entries.
!
!    Output, integer PERM(NODE_NUM), the RCM ordering.
!
!  Local Parameters:
!
!    Local, integer LEVEL_ROW(NODE_NUM+1), the index vector for a level
!    structure.  The level structure is stored in the currently unused
!    spaces in the permutation vector PERM.
!
!    Local, integer MASK(NODE_NUM), marks variables that have been numbered.
!
implicit none

integer adj_num,node_num

integer adj(adj_num)
integer adj_row(node_num+1)
integer i
integer iccsze
integer mask(node_num)
integer level_num
integer level_row(node_num+1)
integer num
integer perm(node_num)
integer root


do i=1,node_num
mask(i) = 1
enddo
num = 1

do i = 1, node_num
!
!  For each masked connected component...
!
if ( mask(i).ne. 0 ) then

root = i

!
!  Find a pseudo-peripheral node ROOT.  The level structure found by
!  ROOT_FIND is stored starting at PERM(NUM).
!
call root_find ( root, adj_num, adj_row, adj, mask, &
& level_num,    level_row, perm(num), node_num )
!
!  RCM orders the component using ROOT as the starting node.
!
call rcm ( root, adj_num, adj_row, adj, mask, perm(num), &
& iccsze,    node_num )

num = num + iccsze
!
!  We can stop once every node is in one of the connected components.
!
if ( node_num .lt. num ) then
return
endif

endif

enddo

return
end

subroutine rcm ( root, adj_num, adj_row, adj, mask, perm, iccsze, &
& node_num )

!*****************************************************************************80
!
!! RCM renumbers a connected component by the reverse Cuthill McKee algorithm.
!
!  Discussion:
!
!    The connected component is specified by a node ROOT and a mask.
!    The numbering starts at the root node.
!
!    An outline of the algorithm is as follows:
!
!    X(1) = ROOT.
!
!    for ( I = 1 to N-1)
!      Find all unlabeled neighbors of X(I),
!      assign them the next available labels, in order of increasing degree.
!
!    When done, reverse the ordering.
!
!  Modified:
!
!    02 January 2007
!
!  Author:
!
!    Alan George, Joseph Liu
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Input, integer ROOT, the node that defines the connected component.
!    It is used as the starting point for the RCM ordering.
!
!    Input, integer ADJ_NUM, the number of adjacency entries.
!
!    Input, integer ADJ_ROW(NODE_NUM+1).  Information about row I is stored
!    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
!    Input, integer ADJ(ADJ_NUM), the adjacency structure.
!    For each row, it contains the column indices of the nonzero entries.
!
!    Input/output, integer MASK(NODE_NUM), a mask for the nodes.  Only
!    those nodes with nonzero input mask values are considered by the
!    routine.  The nodes numbered by RCM will have their mask values
!    set to zero.
!
!    Output, integer PERM(NODE_NUM), the RCM ordering.
!
!    Output, integer ICCSZE, the size of the connected component
!    that has been numbered.
!
!    Input, integer NODE_NUM, the number of nodes.
!
!  Local Parameters:
!
!    Workspace, integer DEG(NODE_NUM), a temporary vector used to hold
!    the degree of the nodes in the section graph specified by mask and root.
!
implicit none

integer adj_num
integer node_num

integer adj(adj_num)
integer adj_row(node_num+1)
integer deg(node_num)
integer fnbr
integer i
integer iccsze
integer j
integer jstop
integer jstrt
integer k
integer l
integer lbegin
integer lnbr
integer lperm
integer lvlend
integer mask(node_num)
integer nbr
integer node
integer perm(node_num)
integer root
!
!  Find the degrees of the nodes in the component specified by MASK and ROOT.
!
call degree ( root, adj_num, adj_row, adj, mask, deg, iccsze, &
& perm, node_num )

mask(root) = 0

if ( iccsze .le. 1 ) then
return
end if

lvlend = 0
lnbr = 1
!
!  LBEGIN and LVLEND point to the beginning and
!  the end of the current level respectively.
!
do while ( lvlend .lt. lnbr )

lbegin = lvlend + 1
lvlend = lnbr

do i = lbegin, lvlend
!
!  For each node in the current level...
!
node = perm(i)
jstrt = adj_row(node)
jstop = adj_row(node+1) - 1
!
!  Find the unnumbered neighbors of NODE.
!
!  FNBR and LNBR point to the first and last neighbors
!  of the current node in PERM.
!
fnbr = lnbr + 1

do j = jstrt, jstop

nbr = adj(j)

if ( mask(nbr) .ne. 0 ) then
lnbr = lnbr + 1
mask(nbr) = 0
perm(lnbr) = nbr
end if

end do
!
!  If no neighbors, skip to next node in this level.
!
!c            if ( lnbr .le. fnbr ) then
!c               cycle
!c            end if
if ( lnbr .gt. fnbr ) then
!
!  Sort the neighbors of NODE in increasing order by degree.
!  Linear insertion is used.
!
k = fnbr

do while ( k .lt. lnbr )

l = k
k = k + 1
nbr = perm(k)

do while ( fnbr .lt. l )

lperm = perm(l)

if ( deg(lperm) .le. deg(nbr) ) then
exit
end if

perm(l+1) = lperm
l = l - 1

end do

perm(l+1) = nbr

end do
end if

end do

end do
!
!  We now have the Cuthill-McKee ordering.  Reverse it.
!
k=iccsze/2
l=iccsze
do i=1,k
lperm=perm(l)
perm(l)=perm(i)
perm(i)=lperm
l=l-1
enddo

return
end
subroutine root_find ( root, adj_num, adj_row, adj, mask, &
& level_num,   level_row, level, node_num )

!*****************************************************************************80
!
!! ROOT_FIND finds a pseudo-peripheral node.
!
!  Discussion:
!
!    The diameter of a graph is the maximum distance (number of edges)
!    between any two nodes of the graph.
!
!    The eccentricity of a node is the maximum distance between that
!    node and any other node of the graph.
!
!    A peripheral node is a node whose eccentricity equals the
!    diameter of the graph.
!
!    A pseudo-peripheral node is an approximation to a peripheral node;
!    it may be a peripheral node, but all we know is that we tried our
!    best.
!
!    The routine is given a graph, and seeks pseudo-peripheral nodes,
!    using a modified version of the scheme of Gibbs, Poole and
!    Stockmeyer.  It determines such a node for the section subgraph
!    specified by MASK and ROOT.
!
!    The routine also determines the level structure associated with
!    the given pseudo-peripheral node; that is, how far each node
!    is from the pseudo-peripheral node.  The level structure is
!    returned as a list of nodes LS, and pointers to the beginning
!    of the list of nodes that are at a distance of 0, 1, 2, ...,
!    NODE_NUM-1 from the pseudo-peripheral node.
!
!  Modified:
!
!    28 October 2003
!
!  Author:
!
!    Alan George, Joseph Liu
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!    Norman Gibbs, William Poole, Paul Stockmeyer,
!    An Algorithm for Reducing the Bandwidth and Profile of a Sparse Matrix,
!    SIAM Journal on Numerical Analysis,
!    Volume 13, pages 236-250, 1976.
!
!    Norman Gibbs,
!    Algorithm 509: A Hybrid Profile Reduction Algorithm,
!    ACM Transactions on Mathematical Software,
!    Volume 2, pages 378-387, 1976.
!
!  Parameters:
!
!    Input/output, integer ROOT.  On input, ROOT is a node in the
!    the component of the graph for which a pseudo-peripheral node is
!    sought.  On output, ROOT is the pseudo-peripheral node obtained.
!
!    Input, integer ADJ_NUM, the number of adjacency entries.
!
!    Input, integer ADJ_ROW(NODE_NUM+1).  Information about row I is stored
!    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
!    Input, integer ADJ(ADJ_NUM), the adjacency structure.
!    For each row, it contains the column indices of the nonzero entries.
!
!    Input, integer MASK(NODE_NUM), specifies a section subgraph.  Nodes
!    for which MASK is zero are ignored by FNROOT.
!
!    Output, integer LEVEL_NUM, is the number of levels in the level structure
!    rooted at the node ROOT.
!
!    Output, integer LEVEL_ROW(NODE_NUM+1), LEVEL(NODE_NUM), the
!    level structure array pair containing the level structure found.
!
!    Input, integer NODE_NUM, the number of nodes.
!
implicit none

integer adj_num
integer node_num

integer adj(adj_num)
integer adj_row(node_num+1)
integer iccsze
integer j
integer jstrt
integer k
integer kstop
integer kstrt
integer level(node_num)
integer level_num
integer level_num2
integer level_row(node_num+1)
integer mask(node_num)
integer mindeg
integer nabor
integer ndeg
integer node
integer root
!
!  Determine the level structure rooted at ROOT.
!
call level_set ( root, adj_num, adj_row, adj, mask, level_num, &
& level_row, level, node_num )
!
!  Count the number of nodes in this level structure.
!
iccsze = level_row(level_num+1) - 1
!
!  Extreme case:
!    A complete graph has a level set of only a single level.
!    Every node is equally good (or bad).
!
if ( level_num .eq. 1 ) then
return
end if
!
!  Extreme case:
!    A "line graph" 0--0--0--0--0 has every node in its only level.
!    By chance, we've stumbled on the ideal root.
!
if ( level_num .eq. iccsze ) then
return
end if
!
!  Pick any node from the last level that has minimum degree
!  as the starting point to generate a new level set.
!
do

mindeg = iccsze

jstrt = level_row(level_num)
root = level(jstrt)

if ( jstrt .lt. iccsze ) then

do j = jstrt, iccsze

node = level(j)
ndeg = 0
kstrt = adj_row(node)
kstop = adj_row(node+1) - 1

do k = kstrt, kstop
nabor = adj(k)
if ( 0 .lt. mask(nabor) ) then
ndeg = ndeg + 1
end if
end do

if ( ndeg .lt. mindeg ) then
root = node
mindeg = ndeg
end if

end do

end if
!
!  Generate the rooted level structure associated with this node.
!
call level_set ( root, adj_num, adj_row, adj, mask, &
& level_num2,   level_row, level, node_num )
!
!  If the number of levels did not increase, accept the new ROOT.
!
if ( level_num2 .le. level_num ) then
exit
end if

level_num = level_num2
!
!  In the unlikely case that ROOT is one endpoint of a line graph,
!  we can exit now.
!
if ( iccsze .le. level_num ) then
exit
end if

end do

return
end

subroutine level_set ( root, adj_num, adj_row, adj, mask, &
& level_num,   level_row, level, node_num )

!*****************************************************************************80
!
!! LEVEL_SET generates the connected level structure rooted at a given node.
!
!  Discussion:
!
!    Only nodes for which MASK is nonzero will be considered.
!
!    The root node chosen by the user is assigned level 1, and masked.
!    All (unmasked) nodes reachable from a node in level 1 are
!    assigned level 2 and masked.  The process continues until there
!    are no unmasked nodes adjacent to any node in the current level.
!    The number of levels may vary between 2 and NODE_NUM.
!
!  Modified:
!
!    28 October 2003
!
!  Author:
!
!    Alan George, Joseph Liu
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Input, integer ROOT, the node at which the level structure
!    is to be rooted.
!
!    Input, integer ADJ_NUM, the number of adjacency entries.
!
!    Input, integer ADJ_ROW(NODE_NUM+1).  Information about row I is stored
!    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
!    Input, integer ADJ(ADJ_NUM), the adjacency structure.
!    For each row, it contains the column indices of the nonzero entries.
!
!    Input/output, integer MASK(NODE_NUM).  On input, only nodes with nonzero
!    MASK are to be processed.  On output, those nodes which were included
!    in the level set have MASK set to 1.
!
!    Output, integer LEVEL_NUM, the number of levels in the level
!    structure.  ROOT is in level 1.  The neighbors of ROOT
!    are in level 2, and so on.
!
!    Output, integer LEVEL_ROW(NODE_NUM+1), LEVEL(NODE_NUM), the rooted
!    level structure.
!
!    Input, integer NODE_NUM, the number of nodes.
!
implicit none

integer adj_num
integer node_num

integer adj(adj_num)
integer adj_row(node_num+1)
integer i
integer iccsze
integer j
integer jstop
integer jstrt
integer lbegin
integer level_num
integer level_row(node_num+1)
integer level(node_num)
integer lvlend
integer lvsize
integer mask(node_num)
integer nbr
integer node
integer root

mask(root) = 0
level(1) = root
level_num = 0
lvlend = 0
iccsze = 1
!
!  LBEGIN is the pointer to the beginning of the current level, and
!  LVLEND points to the end of this level.
!
do

lbegin = lvlend + 1
lvlend = iccsze
level_num = level_num + 1
level_row(level_num) = lbegin
!
!  Generate the next level by finding all the masked neighbors of nodes
!  in the current level.
!
do i = lbegin, lvlend

node = level(i)
jstrt = adj_row(node)
jstop = adj_row(node+1) - 1

do j = jstrt, jstop

nbr = adj(j)

if ( mask(nbr) .ne. 0 ) then
iccsze = iccsze + 1
level(iccsze) = nbr
mask(nbr) = 0
end if

end do

end do
!
!  Compute the current level width (the number of nodes encountered.)
!  If it is positive, generate the next level.
!
lvsize = iccsze - lvlend

if ( lvsize .le. 0 ) then
exit
end if

end do

level_row(level_num+1) = lvlend + 1
!
!  Reset MASK to 1 for the nodes in the level structure.
!
do i =1 ,iccsze
mask(level(i)) = 1
enddo
return
end


subroutine degree ( root, adj_num, adj_row, adj, mask, deg, &
& iccsze, ls,  node_num )

!*****************************************************************************80
!
!! DEGREE computes the degrees of the nodes in the connected component.
!
!  Discussion:
!
!    The connected component is specified by MASK and ROOT.
!    Nodes for which MASK is zero are ignored.
!
!  Modified:
!
!    05 January 2003
!
!  Author:
!
!    Alan George, Joseph Liu
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Input, integer ROOT, the node that defines the connected component.
!
!    Input, integer ADJ_NUM, the number of adjacency entries.
!
!    Input, integer ADJ_ROW(NODE_NUM+1).  Information about row I is stored
!    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
!    Input, integer ADJ(ADJ_NUM), the adjacency structure.
!    For each row, it contains the column indices of the nonzero entries.
!
!    Input, integer MASK(NODE_NUM), is nonzero for those nodes which are
!    to be considered.
!
!    Output, integer DEG(NODE_NUM), contains, for each  node in the connected
!    component, its degree.
!
!    Output, integer ICCSIZE, the number of nodes in the connected component.
!
!    Output, integer LS(NODE_NUM), stores in entries 1 through ICCSIZE the nodes
!    in the connected component, starting with ROOT, and proceeding
!    by levels.
!
!    Input, integer NODE_NUM, the number of nodes.
!
implicit none

integer adj_num
integer node_num

integer adj(adj_num)
integer adj_row(node_num+1)
integer deg(node_num)
integer i
integer iccsze
integer ideg
integer j
integer jstop
integer jstrt
integer lbegin
integer ls(node_num)
integer lvlend
integer lvsize
integer mask(node_num)
integer nbr
integer node
integer root
!
!  The sign of ADJ_ROW(I) is used to indicate if node I has been considered.
ls(1) = root
adj_row(root) = -adj_row(root)
lvlend = 0
iccsze = 1
!
!  LBEGIN is the pointer to the beginning of the current level, and
!  LVLEND points to the end of this level.
do

lbegin = lvlend + 1
lvlend = iccsze
!
!  Find the degrees of nodes in the current level,
!  and at the same time, generate the next level.
do i = lbegin, lvlend

node = ls(i)
jstrt = -adj_row(node)
jstop = abs ( adj_row(node+1) ) - 1
ideg = 0

do j = jstrt, jstop

nbr = adj(j)

if ( mask(nbr) .ne. 0 ) then

ideg = ideg + 1

if ( 0 .le. adj_row(nbr) ) then
adj_row(nbr) = -adj_row(nbr)
iccsze = iccsze + 1
ls(iccsze) = nbr
end if

end if

end do

deg(node) = ideg

end do
!
!  Compute the current level width.
lvsize = iccsze - lvlend
!
!  If the current level width is nonzero, generate another level.
if ( lvsize .eq. 0 ) then
exit
end if

end do
!
!  Reset ADJ_ROW to its correct sign and return.
do i = 1, iccsze
node = ls(i)
adj_row(node) = -adj_row(node)
end do

return
end
