\name{deldir}
\alias{deldir}
\title{
   Delaunay triangulation and Dirichlet tessellation
}
\description{
This function computes the Delaunay triangulation (and hence the
Dirichlet or Voronoi tessellation) of a planar point set according
to the second (iterative) algorithm of Lee and Schacter --- see
\bold{References}.  The triangulation is made to be with respect to
the whole plane by \dQuote{suspending} it from so-called ideal points
(-Inf,-Inf), (Inf,-Inf) (Inf,Inf), and (-Inf,Inf).  The triangulation
is also enclosed in a finite rectangular window.
}
\usage{
deldir(x, y, z=NULL, rw=NULL, eps=1e-09, sort=TRUE, plot=FALSE,
       round=TRUE,digits=6, id=NULL, \dots)
}
\arguments{
\item{x,y}{
These arguments specify the coordinates of the point set being
triangulated/tessellated.  Argument \code{x} may be a numeric
vector or it may be a data structure consisting of a matrix, a
data frame, a generic list, or an object of class \code{"ppp"}.
(See package \code{spatstat}.)  Argument \code{y}, if specified, is
always a numeric vector.  The \code{"x"} and \code{"y"} coordinates
are extracted from arguments \code{x} and \code{y}  according to
(what can be, in some instances) a rather complicated protocol.
See \bold{Notes on extracting the coordinates} for details of
this protocol.
}
\item{z}{
Optional argument specifying \dQuote{auxiliary} values or
\dQuote{tags} associated with the respective points.  (See
\bold{Notes on \dQuote{tags}}.) This argument may be a vector or
factor whose entries constitute these tags, or it may be a text
string naming such a vector or factor.  If \code{z}, or the object
named by \code{z} is a vector (rather than a factor) it may be
of any mode (numeric, character, logical, etc.).  See \bold{Notes
on extracting \code{z}} for how \code{z} is handled when argument
\code{x} is a data structure (rather than a numeric vector).
}
\item{rw}{
The coordinates of the corners of the rectangular window enclosing
the triangulation, in the order (xmin, xmax, ymin, ymax).  Any data
points outside this window are discarded.  If this argument is
omitted, it defaults to values given by the range of the data,
plus and minus 10 percent.
}
\item{eps}{
A value of epsilon used in testing whether a quantity is zero, mainly
in the context of whether points are collinear.  If anomalous errors
arise, it is possible that these may averted by adjusting the value
of eps upward or downward.
}
\item{sort}{
Logical argument; if \code{TRUE} (the default) the data are sorted
into a sequence of \dQuote{bins} prior to triangulation; this makes
the algorithm slightly more efficient.  Normally one would set
\code{sort} equal to \code{FALSE} only if one wished to observe
some of the fine detail of the way in which adding a point to a
data set affected the triangulation, and therefore wished to make
sure that the point in question was added last.  Essentially this
argument would get used only in a debugging process.
}
\item{plot}{
Logical argument; if \code{TRUE} a plot is produced.  The nature
of the plot may be controlled by using the \code{\dots} argument
to pass appropriate arguments to \code{plot.deldir()}.  Without
\dQuote{further instruction} a plot of the points being triangulated
and of both the triangulation and the tessellation is produced;
}
\item{round}{
Logical scalar.  Should the data stored in the returned value
be rounded to \code{digits} decimal digits?  This is essentially
for cosmetic purposes.  This argument is a \dQuote{new addition}
to \code{deldir()}, as of version 0.1-26.  Previously rounding
was done willy-nilly.  The change was undertaken when Kodi
Arfer pointed out that the rounding might have unwanted effects
upon \dQuote{downstream} operations.
}
\item{digits}{
The number of decimal places to which all numeric values in the
returned list should be rounded.  Defaults to 6. Ignored if
\code{round} (see above) is set to \code{FALSE}.
}
\item{id}{
Optional vector of the same length as \code{x} and
\code{y}, providing identifiers for the points that are being
triangulated/tessellated.  This vector will be coerced to mode
\code{"character"} internally.  If \code{id} is not specified,
it defaults to the sequence of indices of the points.  Note that
\code{id} must specify identifiers for \emph{all} points, even
though some of these points may be dropped from consideration if
they do not lie inside the \dQuote{rectangular window} \code{rw}.
Note also that \code{id} must contain no duplicate values.

The \code{id} argument should not be confused with the \code{z}
argument.  The two are entirely different and serve entirely different
purposes.
}

\item{...}{
Auxiliary arguments \code{add}, \code{wlines},
\code{labelPts}, \code{lex}, \code{col}, \code{lty}, \code{pch},
\code{xlim}, and \code{ylim} (and possibly other plotting parameters)
may be passed to \code{\link{plot.deldir}()} through \code{\dots}
if \code{plot=TRUE}.
}
}
\value{
A list (of class \code{deldir}), invisible if \code{plot=TRUE}, with components:

\item{delsgs}{
A data frame with 6 columns.  The first 4 entries of each row are the
coordinates of the points joined by an edge of a Delaunay triangle,
in the order \code{(x1,y1,x2,y2)}.  The last two entries are the
identifiers (defaulting to the indices, if \code{id} is not supplied)
of the two points which are joined.
}
\item{dirsgs}{
A data frame with 10 columns.  The first 4 entries of each
row are the coordinates of the endpoints of one the edges of a
Dirichlet tile, in the order \code{(x1,y1,x2,y2)}.  The fifth and
sixth entries, in the columns named \code{ind1} and \code{ind2},
are the identifiers (defaulting to the indices, if \code{id} is
not supplied)  of the two points, in the set being triangulated,
which are separated by that edge. The seventh and eighth entries,
in the columns named \code{bp1} and \code{bp2} are logical values.
The entry in column \code{bp1} indicates whether the first endpoint
of the corresponding edge of a Dirichlet tile is a boundary point
(a point on the boundary of the rectangular window).  Likewise for
the entry in column \code{bp2} and the second endpoint of the edge.

The ninth and tenth entries, in columns named \code{thirdv1}
and \code{thirdv2} are the identifiers of the respective third
vertices of the Delaunay triangles whose circumcentres constitute
the corresponding endpoints of the edge under consideration.
(The other two vertices of the triangle in question are identified by
the entries of columns \code{ind1} and \code{ind2}.)

The entries of columns \code{thirdv1} and \code{thirdv2} may (also)
take the values $-1, -2, -3$, and $-4$. (These will be coerced to
mode character if \code{id} was supplied.)  This will be the case
if the circumcentre in question lies outside of the rectangular
window \code{rw}.  In these circumstances the corresponding
endpoint of the tile edge is the intersection of the line joining
the two circumcentres with the boundary of \code{rw}, and the
numeric value of the entry of column \dQuote{thirdv1} (respectively
\dQuote{thirdv2}) indicates which side.  The numbering follows the
convention for numbering the sides of a plot region in \code{R}:
1 for the bottom side, 2 for the left hand side, 3 for the top side
and 4 for the right hand side.

Note that the entry in column \code{thirdv1} will be negative if and
only if the corresponding entry in column \code{bp1} is \code{TRUE}.
Similarly for columns \code{thirdv2} and \code{bp2}.
}

\item{summary}{
a data frame with 9, 10 or 11 columns and \code{n.data} rows
(see below).  The rows correspond to the points in the set being
triangulated.  Note that the row names are the indices of the points
in the original sequence of points being triangulated/tessellated.
Usually these will be the sequence 1, 2, ..., n.  However if there
were \emph{duplicated} points then the row name corresponding to a
point is the \emph{first} of the indices of the set of duplicated
points in which the given point appears.  The columns are:

\itemize{
\item \code{x} (the \eqn{x}-coordinate of the point)
\item \code{y} (the \eqn{y}-coordinate of the point)
\item \code{id} (the identifers of the points; present only if
\code{id} was specified)
\item \code{z} (the auxiliary values or \dQuote{tags}; present
only if these were specified)
\item \code{n.tri} (the number of Delaunay triangles emanating from
the point)
\item \code{del.area} (1/3 of the total area of all the Delaunay
triangles emanating from the point)
\item \code{del.wts} (the corresponding entry of the \code{del.area}
column divided by the sum of this column)
\item \code{n.tside} (the number of sides --- within the rectangular
window --- of the Dirichlet tile surrounding the point)
\item \code{nbpt} (the number of points in which the Dirichlet tile
intersects the boundary of the rectangular window)
\item \code{dir.area} (the area of the Dirichlet tile surrounding
the point)
\item \code{dir.wts} (the corresponding entry of the \code{dir.area}
column divided by the sum of this column).
}

Note that the factor of 1/3 associated with the del.area column
arises because each triangle occurs three times --- once for
each corner.
}
\item{n.data}{
the number of points in the set which was triangulated, with any
duplicate points eliminated.  It is the same as the number of rows
of \code{summary}.
}
\item{del.area}{
the area of the convex hull of the set of points being triangulated,
as formed by summing the \code{del.area} column of \code{summary}.
}
\item{dir.area}{
the area of the rectangular window enclosing the points being triangulated,
as formed by summing the \code{dir.area} column of \code{summary}.
}
\item{rw}{
the specification of the corners of the rectangular window enclosing
the data, in the order (xmin, xmax, ymin, ymax).
}
\item{ind.orig}{A vector of the indices of the points (x,y) in the
set of coordinates initially supplied to \code{deldir()} before
duplicate points (if any) were removed.  These indices are used
by \code{\link{triang.list}()}.
}}
\section{Side Effects}{
If \code{plot=TRUE} a plot of the triangulation and/or tessellation
is produced or added to an existing plot.
}

\details{
This package had its origins (way back in the mists of time!)  as an
Splus library section named \dQuote{delaunay}.  That library section
in turn was a re-write of a stand-alone Fortran program written
in 1987/88 while the author was with the Division of Mathematics
and Statistics, CSIRO, Sydney, Australia.  This program was an
implementation of the second (iterative) Lee-Schacter algorithm.
The stand-alone Fortran program was re-written as an Splus function
(which called upon dynamically loaded Fortran code) while the author
was visiting the University of Western Australia, May, 1995.

Further revisions were made December 1996. The author gratefully
acknowledges the contributions, assistance, and guidance of Mark
Berman, of D.M.S., CSIRO, in collaboration with whom this project
was originally undertaken.  The author also acknowledges much useful
advice from Adrian Baddeley, formerly of D.M.S., CSIRO (now
Professor of Statistics at Curtin University).
Daryl Tingley of the Department of Mathematics
and Statistics, University of New Brunswick, provided some helpful
insight.  Special thanks are extended to Alan Johnson, of the
Alaska Fisheries Science Centre, who supplied two data sets which
were extremely valuable in tracking down some errors in the code.

Don MacQueen, of Lawrence Livermore National Lab, wrote an Splus
driver function for the old stand-alone version of this software.
That driver, which was available on Statlib, was deprecated in favour
of the Statlib package \dQuote{delaunay}.  Don also collaborated
in the preparation of the latter package.  It is not clear to me
whether the \dQuote{delaunay} package, or indeed Statlib (or indeed
Splus!) still exist.

See the \code{ChangeLog} for information about further revisions
and bug-fixes.
}

\section{Notes on extracting the coordinates}{
The protocol for extracting the \eqn{x} and \eqn{y} coordinates
from the arguments \code{x} and \code{y} is a bit complicated
and confusing.  It is designed to handle a number of different
desiderata and to accommodate various feature-requests that users
have made over the years.  Basically the protocol is:
\itemize{
   \item If \code{x} is a numeric vector and \code{y} is a numeric
   vector then \code{x} is used as the \eqn{x}-coordinates and
   \code{y} is used as the \eqn{y}-coordinates.

   \item If \code{x} is a matrix, a data frame, or a generic list),
   and \code{y} is a numeric vector, then the \eqn{x}-coordinates
   are sought amongst the components of \code{x} and \code{y}
   is used as the \eqn{y}-coordinates.

   \item If \code{x} is a matrix, a data frame, or a generic list
   and \code{y} is not specified or cannot be found, then both the
   \eqn{x}-coordinates and \eqn{y}-coordinates are sought amongst
   the components of \code{x}.

   \item If \code{x} an object of class \code{"ppp"} then both the
   \eqn{x}-coordinates and \eqn{y}-coordinates are taken from the
   components of \code{x}.   If \code{y} is specified, it is ignored
   (with a warning).

   \item If \code{x} is a numeric vector and \code{y} is not specified
   or cannot be found, then an error is thrown.
}
A few more details:
\itemize{
  \item If \code{x} is of class \code{"ppp"} then it will definitely
  have components named \code{"x"} and \code{"y"}.

  \item If \code{x} is a generic list, it \emph{must} have a
  component named \code{"x"} (otherwise an error is thrown),
  and the \eqn{x}-coordinates are set equal to this component.
  If \code{y} is not specified or cannot be found, then a \code{"y"}
  component of \code{x} is sought.  If such a component exists
  then the \eqn{y}-coordinates are set equal to this component.
  Otherwise an error is thrown).

  \item If \code{x} is a matrix or a data frame, the protocol gets
  a bit more intricate.
  \itemize{
     \item If \code{x} has a column named \code{"x"} then this column
     is taken to be the \eqn{x}-coordinates.
     \item Otherwise the \eqn{x}-coordinates are taken to be the
     \emph{first} column of \code{x} that is not named \code{"y"} or
     \code{znm} (where \code{znm} is the name of the object providing
     the \dQuote{tags}, if \dQuote{tags} have been specified).
     \item If there is no such first column (e.g. if there are only
     two columns and these have names \code{"y"} and \code{znm})
     then an error is thrown.
     \item If \code{y} is not specified or cannot be found, and
     if \code{x} has a column named \code{"y"} then this column is
     taken to be the \eqn{y}-coordinates.
     \item Otherwise, in this situation, the \eqn{y}-coordinates
     are taken to be the \emph{first} column of \code{x} that is
     not named \code{"x"} or \code{znm} and is not equal to the
     column previously selected to be the \code{x}-coordinates.
     \item If there is no such first column (e.g. if there are only
     two columns and these have names \code{"x"} and \code{znm}),
     then an error is thrown.
  }
}

Got all that? \code{:-) } If these instructions seem rather
bewildering (and indeed they are!) just keep things simple and make
calls like \code{deldir(x,y)} where \code{x} and \code{y} are numeric
vectors that have been previously assigned in the global environment.
}

\section{Notes on extracting \code{z}}{
If argument \code{x} is a data structure (rather than a numeric
vector) and is \emph{not} an object of class \code{"ppp"} then
\code{z}, if specified and not found, is searched for in \code{x}.
If \code{x} \emph{is} of class \code{"ppp"} then what happens depends
on whether \code{z} was specified or left to take its default value
of \code{NULL}.  In the former case, \code{z} takes the specified
value.  In the latter case the value of \code{"z"} is taken from
the marks of \code{x} provided that \code{x} is indeed a marked
point pattern and that the marks are \emph{atomic} (essentially
provided that the marks are not a data frame). Otherwise \code{z}
is left \code{NULL}, i.e. there are no \dQuote{tags} associated
with the points.
}

\section{Notes on \dQuote{tags}}{
The \dQuote{tags} are simply values that are associated in some way
with the data points and hence with the tiles of the tessellation
produced.  They \bold{DO NOT} affect the tessellation.  In previous
versions of this package (0.2-10 and earlier) the entries of \code{z}
were referred to as \dQuote{weights}.  This terminology has been
changed since it is misleading.  The tessellation produced when
a \code{z} argument is supplied is the same as is it would be
if there were no \code{z} argument (i.e. no \dQuote{weights}).
The \code{deldir} package \bold{DOES NOT do weighted tessellation}.
}

\section{Notes on Memory Allocation}{
It is difficult-to-impossible to determine \emph{a priori} how much
memory needs to be allocated (in the Fortran code) for storing the
edges of the Delaunay triangles and Dirichlet tiles, and for storing
the \dQuote{adjacency list} used by the Lee-Schacter algorithm.
In the code, an attempt is made to allocate sufficient storage.
If, during the course of running the algorithm, the amount of
storage turns out to be inadequate, the algorithm is halted, the
storage is incremented, and the algorithm is restarted (with an
informative message).  This message may be suppressed by wrapping
the call to \code{deldir()} in \code{\link{suppressMessages}()}.
}

\section{Notes on error messages}{

In previous versions of this package, error traps were set in
the underlying Fortran code for 17 different errors. These were
identified by an error number which was passed back up the call stack
and finally printed out by \code{deldir()} which then invisibly
returned a \code{NULL} value.  A glossary of the meanings of the
values of was provided in a file to be found in a file located in the
\code{inst} directory (\dQuote{folder} if you are a Windoze weenie).

This was a pretty shaganappi system.  Consequently, as of version
1.2-1, conversion to \dQuote{proper} error trapping was implemented.
Such error trapping is effected via the \code{rexit()} subroutine
which is now available in \code{R}.  (See \dQuote{Writing R
Extensions}, section 6.2.1.)

Note that when an error is detected, \code{deldir()} now exits with
a genuine error, rather than returning \code{NULL}.  The glossary
of the meanings of \dQuote{error numbers} is now irrelevant and
has been removed from the \code{inst} directory.

An error trap that merits particular mention was introduced in
version \code{0.1-16} of \code{deldir}.  This error trap relates to
\dQuote{triangle problems}.  It was drawn to my attention by Adam
Dadvar (on 18 December, 2018) that in some data sets collinearity
problems may cause the \dQuote{triangle finding} procedure, used
by the algorithm to successively add new points to a tessellation,
to go into an infinite loop.  A symptom of the collinearity is
that the vertices of a putative triangle appear \emph{not} to be
in anticlockwise order irrespective of whether they are presented
in the order \code{i, j, k} or \code{k, j, i}.  The result of this
anomaly is that the procedure keeps alternating between moving to
\dQuote{triangle} \code{i, j, k} and moving to \dQuote{triangle}
\code{k, j, i}, forever.

The error trap in question is set in \code{trifnd}, the triangle
finding subroutine.  It detects such occurrences of \dQuote{clockwise
in either orientation} vertices.  The trap causes the \code{deldir()}
function to throw an error rather than disappearing into a black
hole.

When an error of the \dQuote{triangle problems} nature occurs, a
\emph{possible} remedy is to increase the value of the \code{eps}
argument of \code{deldir()}. (See the \bold{Examples}.)  There may
conceivably be other problems that lead to infinite loops and so I
put in another error trap to detect whether the procedure has
inspected more triangles than actually exist, and if so to throw
an error.

Note that the strategy of increasing the value of \code{eps}
is \emph{probably} the appropriate response in most (if not all)
of the cases where errors of this nature arise.  Similarly this
strategy is \emph{probably} the appropriate response to the errors
\itemize{
    \item Vertices of \dQuote{triangle} are collinear and vertex 2 is not
          between 1 and 3.  Error in circen.
    \item Vertices of triangle are collinear.  Error in dirseg.
    \item Vertices of triangle are collinear.  Error in dirout.
}

However it is impossible to be sure.  The intricacy and numerical
delicacy of triangulations is too great for anyone to be able to
foresee all the possibilities that could arise.

If there is any doubt as to the appropriateness of the
\dQuote{increase \code{eps}} strategy, users are advised to do
their best to explore the data set, graphically or by other means,
and thereby determine what is actually going on and why problems
are occurring.
}

\section{Warnings}{
\enumerate{
\item The process for determining if points are duplicated
changed between versions 0.1-9 and 0.1-10.  Previously there
was an argument \code{frac} for this function, which defaulted
to 0.0001.  Points were deemed to be duplicates if the difference
in \code{x}-coordinates was less than \code{frac} times the width
of \code{rw} and \code{y}-coordinates was less than \code{frac}
times the height of \code{rw}.  This process has been changed to
one which uses \code{\link{duplicated}()} on the data frame whose
columns are \code{x} and \code{y}.

As a result it may happen that points which were previously
eliminated as duplicates will no longer be eliminated.  (And
possibly vice-versa.)

\item The components \code{delsgs} and \code{summary} of the value
returned by \code{deldir()} are now \emph{data frames} rather than
matrices.  The component \code{summary} was changed to allow the
\dQuote{auxiliary} values \code{z} to be of arbitrary mode (i.e.
not necessarily numeric).  The component \code{delsgs} was then
changed for consistency.  Note that the other \dQuote{matrix-like}
component \code{dirsgs} has been a data frame since time immemorial.
}
}
\section{Acknowledgement}{
I would like to express my most warm and sincere thanks to Duncan
Murdoch (Emeritus Professor of Statistics, Western University) for
helping me, with incredible patience and forbearance, to straighten
out my thinking in respect of adjustments that I recently (October
2021) made to the argument processing protocol in the \code{deldir()}
function.  Duncan provided numerous simple examples to demonstrate
when and how things were going wrong, and patiently explained to
me how I was getting one aspect of the protocol backwards.
}
\author{\rolf}

\references{
Lee, D. T. and Schacter, B. J. (1980). Two algorithms for
constructing a Delaunay triangulation, \emph{International Journal
of Computer and Information Sciences} \bold{9} (3), pp. 219 -- 242.

Ahuja, N. and Schacter, B. J. (1983).  Pattern Models.  New York: Wiley.
}
\seealso{
  \code{\link{plot.deldir}()}, \code{\link{tile.list}()}, \code{\link{triang.list}()}
}
\examples{
x <- c(2.3,3.0,7.0,1.0,3.0,8.0)
y <- c(2.3,3.0,2.0,5.0,8.0,9.0)

# Let deldir() choose the rectangular window.
dxy1 <- deldir(x,y)

# User chooses the rectangular window.
dxy2 <- deldir(x,y,rw=c(0,10,0,10))

# Put "dummy" points at the corners of the rectangular
# window, i.e. at (0,0), (10,0), (10,10), and (0,10)
xx <- c(x,0,10,10,0)
yy <- c(y,0,0,10,10)
dxy3 <- deldir(xx,yy,rw=c(0,10,0,10))

# Plot the triangulation created (but not the tessellation).
dxy2 <- deldir(x,y,rw=c(0,10,0,10),plot=TRUE,wl="tr")

# Example of collinearity error.
\dontrun{
    dniP <- deldir(niProperties) # Throws an error
}
    dniP <- deldir(niProperties,eps=1e-8) # No error.

# Example of using data stored in a data frame.
dsw  <- deldir(seaweed)
# Example where "data" is of class "ppp".
dtoy <- deldir(toyPattern)
# The "tags", in dtoy$summary$z, are the marks of the toy ppp
# object which consists of the letters "a","b","c" and "d".

# Artificial example of tags.
set.seed(42)
trees1  <- sample(c("spruce","birch","poplar","shoe"),20,TRUE)
trees2  <- sample(c("fir","maple","larch","palm"),20,TRUE)
egDat   <- data.frame(x=runif(20),y=runif(20),species=trees1)
tagNm   <- "species"
species <- trees2
dd1     <- deldir(egDat) # No tags.
dd2     <- deldir(egDat,z=species)  # Uses trees1 as the tags.
dd3     <- deldir(egDat,z="species") # Same as dd2.
dd4     <- deldir(egDat,z=tagNm) # Same as dd2 and dd3.
spec    <- species
dd5     <- deldir(egDat,z=spec) # Uses trees2 as the tags.

# Duncan Murdoch's examples.  The deldir() function was not
# handling such examples correctly until Duncan kindly set
# me on the right path.
set.seed(123)
dd6 <- deldir(rnorm(32),rnorm(32),rnorm(32))
#
x <- cbind(x = 1:10, junk = 11:20)
y <- 21:30
z <- 31:40
d7 <- deldir(x=x, y=y, z=z)
#
# print(d7$summary) reveals that x is 1:10, y is 21:30
# and z is 31:40; x[,"junk"] is ignored as it should be.
x <- cbind(x = 1:10, "rnorm(10)" = 11:20)
y <- 21:30
z <- 41:50
d8 <- deldir(x=x, y=y, z=rnorm(10))  
#
# print(d8$summary) reveals that x is 1:10, y is 21:30 and z is a
# vector of standard normal values.  Even though x has a column with
# the name of the z argument i.e. "rnorm(10)" (!!!)  the specified
# value of z takes precedence over this column (and, as per the usual
# R syntax) over the object named "z" in the global environment.

# Artificial example of the use of the "id" argument.
set.seed(42)
x  <- runif(50)
y  <- runif(50)
ll <- expand.grid(a=letters[1:10],b=letters[1:10])
aa <- sample(paste0(ll[["a"]],ll[["b"]]),50)
dxy.wid <- deldir(x,y,id=aa)
}

\keyword{spatial}
