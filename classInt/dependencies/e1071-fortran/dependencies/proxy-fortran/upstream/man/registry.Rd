\name{pr_DB}
\alias{pr_DB}
\alias{registry}
\alias{summary.pr_DB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Registry of proximities}
\description{
  Registry containing similarities and distances.
}
\usage{
pr_DB
\special{pr_DB$get_field(name)}
\special{pr_DB$get_fields()}
\special{pr_DB$get_field_names()}
\special{pr_DB$set_field(name, default = NA, type = NA, is_mandatory = FALSE,
                is_modifiable = TRUE, validity_FUN = NULL)}

\special{pr_DB$entry_exists(name)}
\special{pr_DB$get_entry(name)}
\special{pr_DB$get_entries(name = NULL, pattern = NULL)}
\special{pr_DB$get_entry_names(name)}
\special{pr_DB$set_entry(...)}
\special{pr_DB$modify_entry(...)}
\special{pr_DB$delete_entry(name)}

\method{summary}{pr_DB}(object, verbosity = c("short", "long"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{name}{character string representing the name of an
    entry (case-insensitive).}
  \item{pattern}{regular expression to be matched to all fields of class
    \code{"character"} in all entries.}
  \item{default}{optional default value for the field.}
  \item{type}{optional character string specifying the class to be
    required for this field. If \code{type} is a character vector with more
    than two elements, the entries will be used as fixed set of
    alternatives. If \code{type} is not a character string or vector, the
    class will be inferred from the argument given.}
  \item{is_mandatory}{logical specifying whether new entries are required
    to have a value for this field.}
  \item{is_modifiable}{logical specifying whether entries can be changed
    with respect to that field.}
  \item{validity_FUN}{optional function or character string with the name of a
    function that checks the validity of a field entry. Such a function
    gets the value to be investigated as argument, and should stop with an
    error message if the value is not correct.}
  \item{object}{a registry object.}
  \item{verbosity}{controlling the verbosity of the output of the
    summary method for the registry. \code{"short"} gives just a list, \code{"long"}
    also gives the formulas.}
  \item{\dots}{for \code{pr_DB$set_entry} and \code{pr_DB$modify_entry}:
    named list of fields to be modified in or added to the registry (see details).
    This must include the index field (\code{"names"}).}
}
\details{
  \code{pr_DB} represents the registry of all proximity measures
  available. For each
  measure, it comprises meta-information that can be queried and
  extended. Also, new measures can be added. This is done using
  the following accessor functions of the \code{pr_DB} object:
  
  \code{get_field_names()} returns a character
  vector with all field names. \code{get_field()} returns the information
  for a specific field as a list with components named as described
  above. \code{get_fields()} returns a list with all field
  entries. \code{set_field()} is used to create new fields in the
  repository (the default value will be set in all
  entries).

  \code{get_entry_names()} returns a character vector with (the first
  alias of) all entries. \code{entry_exists()} is a predicate checking
  if an entry with the specified alias exists in the
  registry. \code{get_entry()} returns the specified entry if it exists (and, by
  default, gives an error if it does not). \code{get_entries()} is used to
  query more than one entry: either those matching \code{name} exactly, or
  those where the regular expression in \code{pattern} matches \emph{any}
  character field in an entry. By default, all values are
  returned. \code{delete_entry} removes an existing entry from the
  registry (note that only user-provided entries can be deleted).
  \code{set_entry} and \code{modify_entry} require a named list
  of arguments used as field entries.
  At least the \code{names} index field is required. \code{set_entry}
  will check for all other mandatory fields. If specified in the field
  meta data, each field entry and the entry as a whole is checked for
  validity. Note that only user-specified fields and/or entries can be
  modified, the data shipped with the package are read-only.

  The registry fields currently available are as follows:
  
\describe{
  \item{FUN}{Function to register (see below).}
  \item{names}{Character vector with an alias(es) for the measure.}
  \item{PREFUN}{Optional function (or function name) for preprocessing
    code (see below).}
  \item{POSTFUN}{Optional function (or function name) for postprocessing
    code (see below).}
  \item{distance}{logical indicating whether this measure is a distance (\code{TRUE})
    or similarity (\code{FALSE}).}
  \item{convert}{Optional Function or function name for converting
    between similarities and distances when needed.}
  \item{type}{Optional, the scale the measure applies to
    (\code{"metric"}, \code{"ordinal"}, \code{"nominal"},
  \code{"binary"}, or \code{"other"}). If
    \code{NULL}, it is assumed to apply to some other unknown scale.}
  \item{loop}{logical indicating whether \code{FUN} is just a measure,
    and therefore, if \code{\link{dist}} shall do the loop over all pairs of
    observations/variables, or if \code{FUN} does the loop on its own.}
  \item{\code{C_FUN}}{logical indicating whether \code{FUN} is a C function.}
  \item{abcd}{logical; if \code{TRUE} and binary data (or data to be
    interpreted as such) are supplied, the number of concordant and
    discordant pairs is precomputed for every two binary data vectors and
    supplied to the measure function.}
  \item{formula}{Optional character string with the symbolic representation of
    the formula.}
  \item{reference}{Optional reference (character).}
  \item{description}{Optional description (character). Ideally,
    describes the context in which the measure can be applied.}  
}

A function specified as \code{FUN} parameter has mandatory arguments
\code{x} and \code{y} (if \code{abcd} is \code{FALSE}), and \code{a},
\code{b}, \code{c}, \code{d}, \code{n} otherwise. Additionally, it gets
all optional parameters specified by the user in the \code{\dots}
argument of the \code{\link{dist}} and \code{\link{simil}} functions, possibly
changed and/or complemented by the corresponding (optional)
\code{PREFUN} function. It must return the
(diss-)similarity value computed from the arguments.
\code{x} and \code{y} are two vectors from the
data matrix (matrices) supplied. If \code{abcd} is \code{FALSE}, it is
assumed that binary measures will be used, and the number of all
\code{n} concordant and discordant pairs (x_k, y_k)
precomputed and supplied instead of \code{x} and
\code{y}. \code{a}, \code{b}, \code{c}, and \code{d} are the counts of
all (TRUE, TRUE), (TRUE, FALSE), (FALSE, TRUE), and (FALSE, FALSE)
pairs, respectively.

A function specified as \code{PREFUN} parameter has mandatory arguments
\code{x}, \code{y}, \code{p}, and \code{reg_entry}, with \code{y} and
\code{p} possibly being \code{NULL} depending on the task at
hand. \code{x} and \code{y} are the data objects, \code{p} is a
(possibly empty) list with all specified proximity parameters, and
\code{reg_entry} is the registry entry (a named list containing all
information specified in \code{reg_add}).
The preprocessing function is allowed to change all these
information, and if so, is required to return *all* arguments
as a named list in the same order.

A function specified as \code{POSTFUN} parameter has two mandatory
arguments: \code{result} and \code{p}. \code{result} will contain the
computed raw data, i.e. a vector of length \eqn{n * (n - 1) / 2} for
auto-distances (see \code{\link[stats]{dist}} for details on
\code{dist} objects), or a matrix for cross-distances. \code{p} contains
the specified proximity parameters. Post-processing functions need to
return the \code{result} object (even if unmodified).

A function specified as \code{convert} parameter should preserve the
type of its argument.
}
\author{David Meyer \email{David.Meyer@R-project.org}}

\seealso{\code{\link{dist}}}
\examples{
## create a new distance measure
mydist <- function(x,y) x * y

## create a new entry in the registry with two aliases
pr_DB$set_entry(FUN = mydist, names = c("test", "mydist"))

## look it up (index is case insensitive):
pr_DB$get_entry("TEST")

## modify the content of the description field in the new entry
pr_DB$modify_entry(names = "test", description = "foo function")

## create a new field
pr_DB$set_field("New")

## look up the test entry again (two ways)
pr_DB$get_entry("test")
pr_DB[["test"]]

## show total number of entries
length(pr_DB)

## show all entries (short list)
pr_DB$get_entries(pattern = "foo")

## show more details
summary(pr_DB, "long")

## get all entries in a list (and extract first two ones)
pr_DB$get_entries()[1:2]

## get all entries as a data frame (select first 3 fields)
as.data.frame(pr_DB)[,1:3]

## delete test entry
pr_DB$delete_entry("test")

## check if it is really gone
pr_DB$entry_exists("test")

}
\keyword{cluster}% __ONLY ONE__ keyword per line
