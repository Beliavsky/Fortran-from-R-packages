# Rcsdp 0.1.57

* Fix remaining CRAN check issues (CC variable in Solaris - thanks Brian Ripley, ignore function return value)

# Rcsdp 0.1.56

* Added a `NEWS.md` file to track changes to the package.
* Exposes function 'csdp.control' through NAMESPACE
* Adds function 'csdp_minimal' (contributed by Don van den Bergh - https://github.com/vandenman, https://github.com/hcorrada/rcsdp/pull/1)
* Fixes 'significant warning' from use of `R CMD config` to get
CPP variable

# Rcsdp 0.1.52

* Updated CSDP to version 6.1.1
* Fixed bug in constraint matrix conversion from default behavior of mapply
* Updated configure script to fix build errors

# Rcsdp 0.1.50

* Changed version naming
* Added configure.ac to source package
* Fixed issues in Makevars.in requested by Brian Ripley

# Rcsdp 0.1-5

* Fixed bug reported by Dustin Lennon, where all-zero constraint vectors on linear blocks were not properly initialized.

# Rcsdp 0.1-4

* Fixed bug reported by Jacques-Olivier Moussafir <msfr@mac.com> where an error is thrown on matrices having only one non-zero entry
