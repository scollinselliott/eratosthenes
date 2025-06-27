
<!-- README.md is generated from README.Rmd. Please edit that file -->

# <img src="man/figures/logo.png" align="right" width="150px"/> eratosthenes: Archaeological Synchronism

<!-- badges: start -->
<!-- badges: end -->

The `R` package `eratosthenes` aims to provide a coherent foundation for
archaeological chronology-building by incorporating, computationally,
all relevant sources of information on uncertain archaeological or
historical dates. Archaeological dates are often subject to relational
conditions (via seriation or stratigraphic relationships) and absolute
constraints (such as radiocarbon dates, datable artifacts, or other
known historical events, as *termini post* or *ante quos*), which prompt
the use of a joint conditional probability density to convey those
relationships. The date of any one event can then be marginalized from
that full, joint conditional distribution, which is achieved using a
two-stage Gibbs sampler to draw estimates uniformly between potential
earliest and latest bounds. Ancillary functions include checking for
discrepancies in sequences of events and constraining optimal seriations
to known sequences.

While software exists for calibrating and conditioning radiocarbon dates
upon relative constraints, such as [BCal](https://bcal.shef.ac.uk/)
(Buck, Christen, and James 1999) and
[OxCal](https://c14.arch.ox.ac.uk/oxcal.html) (Bronk Ramsey 2009), the
aim of `eratosthenes` is to extend the application of probability theory
more generally to dating all archaeological phenomena, especially the
production dates of artifact types. The method of sampling employed in
`eratosthenes` involves a two-step process, which entails first
performing iterative routines of Gibbs sampling to determine the initial
value for the main sampler, which uses consistent batch means (CBM) and
Monte Carlo standard errors (MCSE) to determine convergence (Flegal,
Haran, and Jones 2008). Furthermore, `eratosthenes` provides tools for
analyzing the impact of events on each other with the conditional
structure stipulated by the investigator, by implementing a
jackknife-style estimator of squared displacement (how much the date of
one event shifts when another is omitted, squared). `Rcpp` is required
for faster Gibbs sampling.

The package is named after Eratosthenes of Cyrene, author of the
*Chronographiai*.

## Installation

To obtain the current development version of `eratosthenes` from GitHub,
install the package in the `R` command line with `devtools`:

``` r
library(devtools)
install_github("scollinselliott/eratosthenes", dependencies = TRUE, build_vignettes = TRUE) 
```

## Usage

The following comments are intended as a general introduction. See
vignettes for more information on package functionality.

The basic objects of interest in `eratosthenes` are:

- **sequences** of relative events, typically stratigraphic deposits,
  but also isolated contexts such as may be part of a frequency or
  contextual seriation
- **finds**, elements which belong to those events, typically artifacts
- **absolute constraints**, as either *termini post* or *ante quos*,
  expressed as samples from a probability density

Information related to these three items must be formatted in objects of
a `list` class, as follows.

### Sequences

Relative sequences should run in order from left (earliest) to right
(latest). All sequences should consist of vectors, contained in a
`list`. In the following example, the object `contexts` contains three
sequences of events.

``` r
x <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J")
y <- c("B", "D", "G", "H", "K")
z <- c("F", "K", "L", "M")
contexts <- list(x, y, z)
```

See also the section [Evaluating Sequences](#evaluating-sequences)
below.

### Finds

Finds should be formatted as a `list` of `list`s, each of which contains
the entries of the following:

- `id` : a unique identifcation number or code
- `assoc` : an element in the sequences `list` to which that find or
  element pertains
- `type` : optional – one or more types, attributes, features, or
  aspects that pertain to that find (`NULL` if none)

In the following example, the `artifacts` object contains six artifacts
which pertain to elements of the sequences contained in `contexts`
above”

``` r
f1 <- list(id = "find01", assoc = "D", type = c("type1", "form1"))
f2 <- list(id = "find02", assoc = "E", type = c("type1", "form2"))
f3 <- list(id = "find03", assoc = "G", type = c("type1", "form1"))
f4 <- list(id = "find04", assoc = "H", type = c("type2", "form1"))
f5 <- list(id = "find05", assoc = "I", type = "type2")
f6 <- list(id = "find06", assoc = "H", type = NULL)
artifacts <- list(f1, f2, f3, f4, f5, f6)
```

### Absolute Constraints

Constraints should be given as two separate `list`s , one for *termini
post quos* and the other for *termini ante quos*. These take the same
form as the finds object, as a `list` of `list`s, with the same
headings, but include one additional heading of `samples` which contains
the absolute dates pertinent to that *t.p.q.* or *t.a.q*.

``` r
coin1 <- list(id = "coin1", assoc = "B", type = NULL, samples = runif(100, -320, -300))
coin2 <- list(id = "coin2", assoc = "G", type = NULL, samples = runif(100, 37, 41))
destr <- list(id = "destr", assoc = "J", type = NULL, samples = 79)

tpq_info <- list(coin1, coin2)
taq_info <- list(destr)
```

Absolute dates can take any form:

- Single dates, e.g., `79` for 79 CE.
- Samples between two potential dates for a date range, e.g., `-91:-88`,
  `seq(-91, -88, length = 10^5)`, or `runif(10^5, -91, -88)` for 91-88
  BCE.
- Samples from a bespoke density, e.g., from a calibrated radiocarbon
  date. `eratosthenes` does not provide functionality for calibrating
  dates, which can be accomplished using preexisting software or
  directly from a calibration curve. The `R` package `Bchron` (Haslett
  and Parnell 2008) provides functions for calibrating dates. As a brief
  example, given an uncalibrated date and its standard deviation, a
  crude sample of calibrated dates can be drawn from the IntCal20 curve
  data, available from IntCal
  [here](https://www.intcal.org/curves/intcal20.14c) (Reimer et al.
  2020), using the following script:

``` r
intcal20 <- read.csv("../path/to/intcal20.14c")

# 14c date mean and st.dev.
mu <- 2040  
sigma <- 30

# samples of 14c dates
uncalib <- round(rnorm(10^5, mu, sigma))

calib <- c()

for (i in 1:length(uncalib)) {
  x <- intcal20$CAL.BP[ intcal20$X14C.age == uncalib[i] ] 
  #g <- intcal20$Sigma[ intcal20$X14C.age == uncalib[i] ]

  if (length(x) > 0) {
    for (j in 1:length(x)) {
      calib <- c(calib, x[j])  
    }
  }
}

# samples of cal BC date
calBC <- 1950 - calib
hist(calBC, breaks = 100)
```

## Estimating Dates

The core approach of `eratosthenes` is a Gibbs sampler, a common Markov
Chain Monte Carlo (MCMC) technique used for dating archaeological
events, above all radiocarbon dates (Geman and Geman 1984; Buck,
Cavanagh, and Litton 1996; Bronk Ramsey 2009). Gibbs sampling however
can take a number of different forms, and so it is worthwhile to
describe explicitly how it is conducted in `eratosthenes`. The precise
method is as follows:

- To initialize, the earliest possible *t.p.q.* and latest possible
  *t.a.q* dates are selected.
- Relative events are put into a single sequence for the purpose of
  sampling.
- To select the initial date for each relative event, a sample is drawn
  uniformly at random between its upper and lower constraints (absolute
  and relative).
  - For each initial date, a subroutine of Gibbs sampling is performed
    in order to avoid catastrophic collapse of dates due to floating
    point errors (e.g., if one has a high number of events compressed
    into a brief span of time).
- After all dates are initialized, the main Gibbs sampler is performed
  for a specified maximum number of samples, which will stop
  automatically if convergence in distribution has been achieved.
  - Given that dates have already been initialized via Gibbs
    subroutines, the need to discard initial samples due to burn-in is
    obviated (or reduced).
  - Convergence is determined using consistent batch means (CBM), which
    divides the samples into batches. For all events (i.e., variates),
    the Monte Carlo standard errors (MCSE) of their batch means are
    computed. If the mean MCSE falls below a specified criterion (by
    default 0.5, to determine the date of an event +/- 1 year), the main
    Gibbs sampler will stop. See Jones et al. (2006) and Flegal, Haran,
    and Jones (2008) for details.
  - Given that this is the mean MCSE of all events, certain events will
    have higher or lower MCSE, and so each event’s MCSE should be
    reported.

There are two functions in `eratosthenes` for estimating dates:

- `gibbs_ad()` is the primary function, which will yield samples for
  dates of deposition, production, and any absolute constraints
  themselves (that is, the density of that extrinsic date as impacted by
  all other events in the joint distribution).
- `gibbs_ad_use()` is the secondary function, for estimating the date of
  use of a find given its production and depositional date.

See the section [Evaluating Displacement](#evaluating-displacement)
below for tools on assessing the effective influence of events upon each
other within the joint conditional density.

### Dates of Production and Deposition

The function `gibbs_ad()` takes as inputs the following objects:

- `sequences` : A `list` of relative sequences of contexts or events.
- `finds` : A `list` of any elements which belong to a context or event,
  which may be assigned a given type.
- `max_samples` : The maximum number of samples to run, which will stop
  the main sampling routine even if convergence has not been achieved
  (default is `10^5`).
- `size` : How many samples to take between each check for convergence
  (default is `10^3`).
- `mcse_crit` : The criterion of the mean MCSE at which to stop the
  sampler (default is `0.5`)
- `tpq` and `taq`: Separate `lists` that indicate any elements that
  provide extrinsic (i.e., absolute) chronological information, as
  *termini post* and *ante quos*.
- `alpha_` and `omega_`: lowest and highest bounds within which to
  sample.
- `trim`: whether to remove contexts from the output that are before or
  after user-provided *t.p.q.* and *t.a.q.* (i.e., those which depend on
  `alpha_` and `omega_`).
- `rule`: the rule for determining the earliest date of production of an
  artifact type. Initial threshold boundaries are first established
  between the earliest depositional context containing an artifact of
  that type and the next earliest context which lacks it. Then, the
  following rules will sample a date accordingly:
  - `naive`: samples are drawn between the initial threshold sample and
    the depositional date of that artifact
  - `earliest`: samples are drawn within the initial threshold
    boundaries

For example, to sample from the sequences, finds, and constraints given
above, the following inputs are entered into the `gibbs_ad()` function:

``` r
result <- gibbs_ad(contexts, finds = artifacts, tpq = tpq_info, taq = taq_info)
```

The output is a `list` object of class `marginals` containing the
following objects:

- `deposition` : a `list` of the marginal densities of the date of the
  final deposition of contexts.
- `externals` : a `list` of the marginal densities of date of any
  *terminus post quem* or *terminus ante quem*, as affected by
  depositional variates in the joint conditional distribution.
- `production` : a `list` of the marginal densities of the production
  date of a given type or class of artifact, given the rule stipulated
  in the input.
- `mcse` : a vector of the MCSE of all events.

Information on the `marginals` object can be accessed with `print()` and
`summary()`. Density plots and density histograms of one more events can
be produced using `plot()` and `hist()` respectively (see documentation
for details).

### Dates of Use

Determining the use date proceeds along the same method of Gibbs
sampling discussed above, using consistent batch means to determine
convergence. To compute a density of the use date of an artifact or
artifact type, an object of class `marginals` is necessary (i.e., the
estimation of the production and depositional dates must first be
performed). Then, the `use_dates()` function will return a density
conditional on the production and depositional dates.

Only one type at a time can be estimated with `use_dates()`. A type can
be defined on the basis of:

- One or more `id`s in the finds list.
- One or more `type`s in the finds list.

That is, one can pool together multiple finds as a type on the basis of
their `id`, even if they were not so explicitly given a `type` in the
finds object. Similarly, one can pool together more than one type of
artifact, e.g., if one is dealing with multiple subtypes and one wants
to evaluate them as a single type (e.g., pooling the labels of “Late
Greco-Italic amphora”, “MGS V amphora”, “MGS VI amphora” into a single
type).

The `gibbs_ad_use()` function, takes the following inputs, similar to
`gibbs_ad()`, but with a field for either `id` or `type` (only one or
the other field must be used):

- `gibbs` : A `list` object of class `marginals`, as computed via
  `gibbs_ad()`.
- `finds` : Either the `list` object of finds originally used as input
  to produce `gibbs`, or a `data.frame` of two columns, the first column
  listing the context and the second the incidence of the id or type in
  that context.
- `id` : A vector of the `id` of one or more specific finds whose use
  date is to be estimated. The values of `id` must match those in the
  `list` of `finds`. If `type` is used, `id` is ignored.
- `type` : A vector of one or more types to estimate a use density for.
  Must contain a value if `id` is left as `NULL`.
- `max_samples`, `size`, `mcse_crit` : The same information used for
  determining the maximum length of the Gibbs sampler and when
  convergence has been achieved.

Using the `result` object above, the densities of the use dates of the
following types is computed using the `use_dates()` function as follows:

``` r
# use dates by specifying ids
gibbs_ad_use(result, artifacts, id = c("find04", "find05"))

# use dates by speciifying types
gibbs_ad_use(result, artifacts, type = "type1")
```

The result is an object of the class `use_marginals`, which contains
information on the density of the date of use as well as the MCSE of the
type specified, in the same fashion as the result of the `gibbs_ad()`
function.

## Graphics

Base R graphics are provided by `eratosthenes` to examine traceplots of
the results of `gibbs_ad()`, as well as produce density histograms of
the results of `gibbs_ad()` and `gibbs_ad_use()`. For `gibbs_ad()`,
histograms may contain up to 12 distinct events. For `gibbs_ad_use()`,
the production, use, and deposition of the stipulated artifact type are
shown.

## Evaluating Sequences

Managing and evaluating the validity of relative sequences consists of
checking multiple partial sequences against one another. Not all
relative sequences are of the same informational validity, and not all
sequences will contain the same elements. An investigator may seek to
constrain one sequence against another, i.e., keeping elements of
sequence as close as possible to one another while reordering only some
of the elements.

Some functions related to relative sequences:

- `seq_check()` sees whether partial sequences agree in their relative
  ordering of elements.
- `seq_adj()` provides the means to coerce an “input” sequence to a
  discrepant “target” sequence which contains fewer elements. E.g., if
  one has obtained an optimal seriation of contexts (of both single,
  unrelated deposits and stratigraphic deposits) as determined by the
  presence/absence of find-types, which conflicts with a sequence
  obtained from a stratigraphic sequence whose physical relationships
  are certain, this function will reorder the optimal seriation, fitting
  any single deposits missing from the stratigraphic sequence
  accordingly.

The package `eratosthenes` does not have functionality to produce
seriations or ordinations, as packages `seriation`, `vegan`, and
`lakhesis` can perform this task already.

## Evaluating Displacement

As real-world joint condition densities will comprise hundreds of events
or more, it is easy for an investigator to loose track of which
relative/absolute events are determinative or influential upon others,
in terms of the estimation of their date. `eratosthenes` assess such
influence within the conditional structure via the estimation of
“displacement.” That is, given the omission of an event $j$ (either a
depositional event or an absolute constraint) from the set of all
events, how much does the estimation of the date of another event
change?

The squared displacement $\delta^2(i,j)$ of a target event $i$ caused by
the omission of $j$ is computed as follows. Let $x_i$ be the estimated
marginalized Monte Carlo mean date using all events within the full
joint conditional, and then let $\tilde{x}_i^{(-j)}$ be the “jackknife”
estimated date, when event $j$ has been omitted from all sequences and
absolute constraints. Squared displacement of $j$ upon $i$ is then:
\begin{equation\*} ^2(i,j) = (\_i^{(-j)} - x_i)^2 \end{equation} If
squared displacement is high, then the omission of $j$ has greatly
shifted the date of $i$. If squared displacement is low, then the
omission of $j$ has not altered the date of $i$ much. Squared
displacement is measured in continuous time, whichever scale the
investigator is using (typically years).

Conversely, one can estimate the effective influence of an event $j$
upon all others by taking the mean squared displacement (MSD). This
involves taking the mean of the squared displacements of all others
events when $j$ is omitted. Where $\Theta$ represents the set of all
relative and absolute events, the MSD is defined as \begin{equation\*} =
\_{i , i j} ^2 (i,j) \end{equation}

The squared displacement and MSD are computed in `eratosthenes` after
running the `gibbs_ad()` function, as follows. Note that squared
displacement may be computed for any event $i$ that represents a
relative or absolute constraint, as well as an artifact production date,
while $j$ can only be a relative event or absolute constraint (it would
make no sense to omit an artifact production date, since these are
conditional upon relaive/absolute dates to begin with). Similarly, MSD
can only be computed for relative/absolute events.

Objects in the example below are provided from the section
[Usage](#usage) above. As these routines are fairly intensive,
computational time can be reduced by lowering the values of
`max_samples` and/or raising `mcse_crit`.

``` r
# run gibbs_ad() first
result <- gibbs_ad(contexts, finds = artifacts, tpq = tpq_info, taq = taq_info)

# squared displacement is estimated for a target event ("j" above) and all other events

# squared displacement for depositional context "E"
sq_disp(result, target = "E", sequences = contexts, 
        max_samples = 20000, mcse_crit = 2, tpq = tpq_info, taq = taq_info)

# squared displacement for production of artifact type "type1"
sq_disp(result, target = "type1", sequences = contexts, finds = artifacts,
        max_samples = 20000, mcse_crit = 2, tpq = tpq_info, taq = taq_info)
        
# mean squared displacement (MSD) is estimated for all relative and absolute dates
result_msd <- msd(result, contexts, finds = artifacts,
                  mcse_crit = 1, tpq = tpq_info, taq = taq_info)
```

## Bibliography

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-bronk_ramsey_bayesian_2009" class="csl-entry">

Bronk Ramsey, C. 2009. “Bayesian Analysis of Radiocarbon Dates.”
*Radiocarbon* 51: 337–60.

</div>

<div id="ref-buck_bayesian_1996" class="csl-entry">

Buck, C. E., W. G. Cavanagh, and C. D. Litton. 1996. *Bayesian Approach
to Interpreting Archaeological Data*. Chichester: John Wiley; Sons.

</div>

<div id="ref-buck_bcal_1999" class="csl-entry">

Buck, C. E., J. A. Christen, and G. N. James. 1999. “BCal: An On-Line
Bayesian Radiocarbon Calibration Tool.” *Internet Archaeology* 7.
<https://intarch.ac.uk/journal/issue7/buck/>.

</div>

<div id="ref-flegal_markov_2008" class="csl-entry">

Flegal, J. M., M. Haran, and G. L. Jones. 2008. “Markov Chain Monte
Carlo: Can We Trust the Third Significant Figure?” *Statistical Science*
23: 250–60. <https://doi.org/10.1214/08-STS257>.

</div>

<div id="ref-geman_stochastic_1984" class="csl-entry">

Geman, S., and D. Geman. 1984. “Stochastic Relaxation, Gibbs
Distributions, and the Bayesian Restoration of Images.” *IEEE
Transactions on Pattern Analysis and Machine Intelligence* 6: 721–41.

</div>

<div id="ref-haslett_simple_2008" class="csl-entry">

Haslett, J., and A. C. Parnell. 2008. “A Simple Monotone Process with
Application to Radiocarbon-Dated Depth Chronologies.” *Journal of the
Royal Statistical Society: Series C (Applied Statistics)* 57: 399–418.
<https://doi.org/10.1111/j.1467-9876.2008.00623.x>.

</div>

<div id="ref-jones_fixed-width_2006" class="csl-entry">

Jones, G. L., M. Haran, B. S. Caffo, and R. Neath. 2006. “Fixed-Width
Output Analysis for Markov Chain Monte Carlo.” *Journal of the American
Statistical Association* 101: 1537–47.
<https://doi.org/10.1198/016214506000000492>.

</div>

<div id="ref-reimer_intcal20_2020" class="csl-entry">

Reimer, P. J., W. E. N. Austin, E. Bard, A. Bayliss, P. G. Blackwell, C.
Bronk Ramsey, M. Butzin, et al. 2020. “The IntCal20 Northern Hemisphere
Radiocarbon Age Calibration Curve (0–55 Cal
<span class="nocase">kBP</span>).” *Radiocarbon* 62: 725–57.
<https://doi.org/10.1017/RDC.2020.41>.

</div>

</div>
