
<!-- README.md is generated from README.Rmd. Please edit that file -->

# <img src="man/figures/logo.png" align="right" width="150px"/> eratosthenes: Tools for Archaeological Synchronism

<!-- badges: start -->
<!-- badges: end -->

The `R` package `eratosthenes` aims to provide a coherent foundation for
archaeological chronology-building by incorporating, computationally,
all relevant sources of information on uncertain archaeological or
historical dates. Archaeological dates are often subject to relational
conditions (via seriation or stratigraphic relationships) and absolute
constraints (such as radiocarbon dates, datable artifacts, or other
known historical events, as *termini post* or *ante quem*), which prompt
the use of a joint conditional probability density to convey those
relationships. The date of any one event can then be marginalized from
that full, joint conditional distribution, which is achieved using a
Gibbs sampler to draw estimates uniformly between potential earliest and
latest bounds. Ancillary functions include checking for discrepancies in
sequences of events and constraining optimal seriations to known
sequences.

Dates are estimated for the following types of events:

- **deposition**: the marginal density of the date of the final
  deposition of a context or find-type.
- **externals**: the marginal density of date of any *terminus post
  quem* or *terminus ante quem*, as affected by depositional variates in
  the joint conditional distribution.
- **production**: the marginal density of the production date of a given
  type or class of artifact, given the stipulation that the type’s
  earliest date of production lies before its earliest date of
  deposition and after the depositional date of the context immediately
  prior.

See vignettes for more information on the package functionality. The
package requires `Rcpp` for faster Gibbs sampling.

The package is named after Eratosthenes of Cyrene, author of the
*Chronographiai*.

## Installation

To obtain the current development version of `eratosthenes` from GitHub,
install the package in the `R` command line with `devtools`:

``` r
library(devtools)
install_github("scollinselliott/eratosthenes", dependencies = TRUE, build_vignettes = TRUE) 
```
