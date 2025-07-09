
<!-- README.md is generated from README.Rmd. Please edit that file -->

# <img src="man/figures/logo.png" align="right" width="150px"/> eratosthenes: Archaeological Synchronism

<!-- badges: start -->
<!-- badges: end -->

The `R` package `eratosthenes` aims to provide a general, flexible
toolkit for archaeological chronology-building by incorporating,
computationally, all relevant sources of information on uncertain
archaeological or historical dates. Archaeological dates are subject to
relational conditions (via seriation or stratigraphic relationships) and
absolute constraints (such as radiocarbon dates, datable artifacts, or
other known historical events, as *termini post* or *ante quos*), which
prompt the use of a joint conditional probability density to convey
those relationships. The date of any one event can then be marginalized
from that full, joint conditional distribution.

While software exists for calibrating and conditioning radiocarbon dates
upon relative constraints, such as `OxCal` (Bronk Ramsey 2009) and
`BCal` (Buck, Christen, and James 1999), as well as R packages `oxcAAR`
(Hinz et al. 2021), `Bchron` (Haslett and Parnell 2008), and `rcarbon`
(Crema, Bevan, and Shennan 2017), along with software for general
chronological modeling like `Chronomodel` (Lanos and Philippe 2017) and
`ChronoLog` (Levy et al. 2021), formal methods for dating artifacts and
artifact types are lacking. One of the major goals of `eratosthenes` is
to advance the synchronism of chronologies and the crafting of
large-scale chronological models that rely heavily upon artifact
typologies. The package therefore facilitates the marginalization of
dates of a type’s production, use, and deposition. The method of
sampling employed in `eratosthenes` involves a two-step process of Gibbs
sampling, using consistent batch means (CBM) and Monte Carlo standard
errors (MCSE) to determine convergence (Jones et al. 2006; Flegal,
Haran, and Jones 2008). Finaly, `eratosthenes` provides tools for
analyzing the impact of events on each other with the conditional
structure stipulated by the investigator, by implementing a
jackknife-style estimator of squared displacement (how much the date of
one event shifts when another is omitted). Ancillary functions include
checking for discrepancies in sequences of events and constraining
optimal seriations to known sequences.

The package is motivated by a philosophy of generalism and minimalism,
eschewing the following:

- intervals or durative events. If desired, such instances can be
  asserted as two separate point events in sequences, e.g.,
  `"X - Start"` and `"X - End"`.
- periods and phases. Periods and phases are not actual material or
  behavioral events, but ideal (and often contested) constructs used to
  make sense of the past. If desired, an investigator can always enter
  period-related events, e.g., `"Archaic Period - Start"`, into their
  list of sequences.
- discretization of time into intervals. Samples are drawn along the
  continuum.
- overly cumbersome chronological relationships. As `eratosthenes`
  samples points along the continuum, there is only before and after. If
  desired, overlaping events can be expressed in sequence construction:
  e.g., for sequences `c("A", "B", "C")` and `c("A", "D", "E", "C")`,
  events `"B"` and `"D", "E"` will overlap with each other.

The focus of the package is on the structure of the joint conditional,
rather than specific probability models. Hence, `eratosthenes` relies on
the continuous uniform for estimating relative events. Any model can
however be used for absolute constraints, from single points to
customized densities.

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

The basic items of interest in `eratosthenes` are:

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

Finds should be formatted as a `list` of `lists`, each of which contains
the entries of the following:

- `id` : a unique identification number or code
- `assoc` : an element in the sequences `list` to which that find or
  element pertains
- `type` : optional – one or more types, attributes, features, or
  aspects that pertain to that find (`NULL` if none)
- `residual` : optional – if `TRUE`, it means that the find is
  considered residual to the context (the sequential event) in which it
  was found, and will not be considered when estimating aspects of the
  date of a type (production, use, and deposition).

In the following example, the `artifacts` object contains six artifacts
which pertain to elements of the sequences contained in `contexts`
above”

``` r
f1 <- list(id = "find01", assoc = "D", type = c("type1", "form1"))
f2 <- list(id = "find02", assoc = "E", type = c("type1", "form2"))
f3 <- list(id = "find03", assoc = "G", type = c("type1", "form1"), residual = TRUE)
f4 <- list(id = "find04", assoc = "H", type = c("type2", "form1"))
f5 <- list(id = "find05", assoc = "I", type = "type2")
f6 <- list(id = "find06", assoc = "H", type = NULL)
artifacts <- list(f1, f2, f3, f4, f5, f6)
```

### Absolute Constraints

Constraints should be given as two separate `lists`, one for *termini
post quos* and the other for *termini ante quos*. These take the same
form as the finds object, as a `list` of `lists`, with the same
headings, but include one additional heading of `samples` which contains
the absolute dates pertinent to that *t.p.q.* or *t.a.q*.

``` r
coin1 <- list(id = "coin1", assoc = "B", type = NULL, samples = runif(100, -320, -300))
coin2 <- list(id = "coin2", assoc = "G", type = NULL, samples = runif(100, 37, 41))
destr <- list(id = "destr", assoc = "J", type = NULL, samples = 79)

tpq_info <- list(coin1, coin2)
taq_info <- list(destr)
```

It can be noted that absolute constraints can belong to a type. Any
artifact which carries absolute dating information (i.e., extrinsic to
the joint conditional density) should be assigned as an absolute
constraint. It assumed that if a *t.p.q* has a type, it refers to the
artifact’s date of production, and is treated as such (see the section
[Dates of the Production, Use, and Deposition of a
Type](#dates-of-the-production-use-and-deposition%20of-a-type) below).

Absolute dates can take any form:

- Single dates, e.g., `79` for 79 CE.
- Samples between two potential dates for a date range, e.g., `-91:-88`,
  `seq(-91, -88, length = 10^5)`, or `runif(10^5, -91, -88)` for 91-88
  BCE.
- Samples from a bespoke density, e.g., from a calibrated radiocarbon
  date. `eratosthenes` does not provide functionality for calibrating
  dates, which can be accomplished using preexisting software or
  directly from a calibration curve. As a brief example, given an
  uncalibrated date and its standard deviation, a crude sample of
  calibrated dates can be drawn from the IntCal20 curve data, available
  from IntCal [here](https://www.intcal.org/curves/intcal20.14c) (Reimer
  et al. 2020), using the following script:

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
- Relative events are indexed along a single sequence for the purpose of
  sampling (this does not change their conditional relationships).
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

- `gibbs_ad()` estimates the marginal density of the date of events in
  sequences and absolute constraints (*t.p./a.q*).
- `gibbs_ad_type()` estimates densities of the date of the production,
  use, and deposition of a specified artifact type, given seuqences and
  cosntraints.

See the section [Evaluating Displacement](#evaluating-displacement)
below for tools on assessing the effective influence of events upon each
other within the joint conditional density.

### Dates of Events in Sequences and Absolute Constraints

The function `gibbs_ad()` takes as inputs the following objects:

- `sequences` : A `list` of relative sequences of contexts or events.
- `max_samples` : The maximum number of samples to run, which will stop
  the main sampling routine even if convergence has not been achieved
  (default is `10^5`).
- `size` : How many samples to take between each check for convergence
  (default is `10^3`).
- `mcse_crit` : The criterion of the mean MCSE at which to stop the
  sampler (default is `0.5`)
- `tpq` and `taq`: Separate `lists` that indicate any elements that
  provide extrinsic (i.e., absolute) chronological information, as
  *termini post* and *ante quos*. Format must follow that illustrated in
  the Section above on [Absolute Constraints](#absolute-constraints).
- `alpha_` and `omega_`: lowest and highest bounds within which to
  sample.
- `trim`: whether to remove contexts from the output that are before or
  after user-provided *t.p.q.* and *t.a.q.* (i.e., those which depend on
  `alpha_` and `omega_`).

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
- `mcse` : a vector of the MCSE of all events.

Information on the `marginals` object can be accessed with `print()` and
`summary()`. Density plots and density histograms of one more events can
be produced using `plot()` and `histogram()` respectively (see packag
documentation for details).

### Dates of the Production, Use, and Deposition of a Type

Determining the date of the production, use, and deposition of an
artifact type uses the same method of Gibbs sampling discussed above,
i.e., consistent batch means to determine convergence. Given that types
are ideal constructs used to categorize artifacts, the notion of a
“type” here has flexibility. While only one “type” at a time can be
estimated with `gibbs_ad_type()`, here, a “type” can be defined on the
basis of:

- One or more `id` in the finds list.
- One or more `type` in the finds list.

That is, one can pool together multiple finds as a type on the basis of
their `id`, even if they were not so explicitly given a `type` in the
finds object. Similarly, one can pool together more than one type of
artifact, e.g., if one is dealing with multiple subtypes and one wants
to evaluate them as a single type (e.g., pooling the labels of “Late
Greco-Italic amphora”, “MGS V amphora”, “MGS VI amphora” into a single
type).

The function works on the principle of the presence/absence of the
specified type in a given context. First, it identifies all contexts in
the sequences to which it has been assigned (i.e., been deposited).
Then, it uses a stipulated rule to identify the earliest moment of
production, contingent upon its earliest absence within the joint
conditional density (see the argument `rule` below). Finally, dates of
use are sampled between production and deposition.

The `gibbs_ad_type()` function takes the following inputs, similar to
`gibbs_ad()`, but with some additional fields:

- `sequences` : A `list` of relative sequences of contexts or events.
- `finds` : Either the `list` object of finds originally used as input
  to produce `gibbs`, or a `data.frame` of two columns, the first column
  listing the context and the second the incidence of the id or type in
  that context.
  - If a find entry contains the expression `"residual = TRUE"`, it
    indicates that its association with the context should not be taken
    into account. Primarily, this indicates that a finds depositional
    date occured prior to the context it pertains to (i.e., it has been
    redeposited from an earlier time), but it can also be used to
    suppress the association of finds which may be spurious.
- `id` : A vector of the `id` of one or more specific finds whose use
  date is to be estimated. The values of `id` must match those in the
  `list` of `finds`. If `type` is used, `id` is ignored.
- `type` : A vector of one or more types to estimate a use density for.
  Must contain a value if `id` is left as `NULL`.
- `type_name` : A customized label for the type (e.g., if one is pooling
  together multiple `id`/`type` entries). If only one `type` has been
  entered, that label is used. Otherwise it defaults to just `"Type"`.
- `max_samples`, `size`, `mcse_crit`, `trim` : The same information used
  for determining the maximum length of the Gibbs sampler and when
  convergence has been achieved, as well as whether to trim events, [as
  above](#dates-of-events-in-sequences-and-absolute-constraints). Note
  that the `mcse_crit`, as a stopping rule, applies still to the
  sequential events/absolute constraints, but MCSE will still be
  reported for the estimates of the production, use, and depositional
  dates.
- `tpq` and `taq` : Format must follow that illustrated in the Section
  above on [Absolute Constraints](#absolute-constraints).
- `rule`: the rule for determining the earliest date of production of an
  artifact type. Initial threshold boundaries are first established
  between the earliest depositional context containing an artifact of
  that type and the next earliest context which lacks it. Then, the
  following rules will sample a date accordingly:
  - `naive`: samples are drawn between the initial threshold sample and
    the depositional date of that artifact.
  - `earliest`: samples are drawn within the initial threshold
    boundaries.

As use dates are drawn between production and depositional dates, if one
chooses `"earliest"` as the rule, then the use density is equivalent to
that of the `"naive"` production density. It should also be noted that,
for this function, Gibbs sampling is only used for the depositional
sequences and absolute constraints, not for production, use, and
deposition (i.e, the use date does not affect the production date, nor
is the depositional date affected by the production date).

Using the `result` object above, the densities of the use dates of the
following types is computed using the `gibbs_ad_type()` function as
follows:

``` r
# use dates by specifying ids
gibbs_ad_type(contexts, artifacts, id = c("find04", "find05"), tpq = tpq_info, taq = taq_info)
# use dates by specifying types
gibbs_ad_type(contexts, artifacts, type = "type1", tpq = tpq_info, taq = taq_info)
```

Adjusting the values of `max_samples` and `mcse_crit` is recommended to
reduce computational time, as needed.

The result is a `list` object of the class `type_marginals`, which
contains information on the densities of the dates of production, use,
and deposition, as well as the MCSE, of the type specified.

## Graphics

Base R graphics are provided by `eratosthenes` to generate traceplots of
the results of `gibbs_ad()` and produce density histograms of the
results of `gibbs_ad()` and `gibbs_ad_type()`. For `gibbs_ad()`,
histograms may contain up to 12 distinct events. For `gibbs_ad_type()`,
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
seriations or ordinations, as R packages such as `seriation` (Hahsler,
Hornik, and Buchcta 2008), `vegan` (Oksanen et al. 2024), `lakhesis`
(Collins-Elliott Under Review), and many others, can perform this task
already.

## Evaluating Displacement

As real-world joint conditional densities will comprise hundreds of
events or more, it is easy for an investigator to loose track of which
relative/absolute events are determinative or influential upon others,
in terms of the estimation of their date. `eratosthenes` assesses such
influence within the conditional structure via the estimation of
“displacement.” That is, given the omission of an event $j$ (either a
depositional event or an absolute constraint) from the set of all
events, how much does the estimation of the date of another event
change?

The squared displacement $\delta^2(i,j)$ of a target event $i$ caused by
the omission of $j$ is computed as follows. Let $\tilde{x}_i$ be the
estimated marginalized Monte Carlo mean date using all events within the
full joint conditional, and then let $\tilde{x}_i^{(-j)}$ be the
“jackknife” estimated date, when event $j$ has been omitted from all
sequences and absolute constraints. Squared displacement of $j$ upon $i$
is then:

$$
\delta^2(i,j) = (\tilde{x}_i^{(-j)} - \tilde{x}_i)^2
$$

If squared displacement is high, then the omission of $j$ has greatly
shifted the date of $i$. If squared displacement is low, then the
omission of $j$ has not altered the date of $i$ much. Squared
displacement is measured in continuous time, whichever scale the
investigator is using (typically years).

Conversely, one can estimate the effective influence of an event $j$
upon all others by taking the mean squared displacement (MSD). This
involves taking the mean of the squared displacements of all other
events when $j$ is omitted. Where $\Theta$ represents the set of all
relative and absolute events, the MSD is defined as

$$
\text{MSD}(j) = \frac{1}{n-1} \sum_{i \in \Theta, i \neq j} \delta^2 (i,j)
$$

The squared displacement and MSD are computed in `eratosthenes` for
relative events and absolute constraints after running the `gibbs_ad()`
function, and for an artifact type after running the `gibbs_ad_type()`
function. Note that squared displacement may be computed for any event
$i$ that represents a relative or absolute constraint, as well as a type
(the use date is used to compute displacement for finds, as it is
affected by both production and deposition) production date, while $j$
can only be a relative event or absolute constraint (it would make no
sense to omit e.g. an artifact production date, since these are
conditional upon relative/absolute dates to begin with). Similarly, MSD
can only be computed for relative/absolute events.

Objects in the example below are provided from the section
[Usage](#usage) above. As these routines are fairly intensive,
computational time can be reduced by lowering the values of
`max_samples` and/or raising `mcse_crit`.

``` r
# run gibbs_ad() first
result <- gibbs_ad(contexts, tpq = tpq_info, taq = taq_info)

# squared displacement for depositional context "E" as the target event ("j" above)
sq_disp(result, target = "E", sequences = contexts, 
        max_samples = 20000, mcse_crit = 2, tpq = tpq_info, taq = taq_info)

# mean squared displacement (MSD) is estimated for all relative and absolute dates
msd(result, contexts, finds = artifacts,
    mcse_crit = 1, tpq = tpq_info, taq = taq_info)

# squared displacement for production of artifact type "type1"
# run gibbs_ad_type() first
result_type1 <- gibbs_ad_type(contexts, finds = artifacts, type = "type1",
                              tpq = tpq_info, taq = taq_info)
sq_disp(result_type1, sequences = contexts, finds = artifacts,
        max_samples = 3000, mcse_crit = 2, tpq = tpq_info, taq = taq_info)
```

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-bronk_ramsey_bayesian_2009" class="csl-entry">

Bronk Ramsey, C. 2009. “Bayesian Analysis of Radiocarbon Dates.”
*Radiocarbon* 51: 337–60. <https://doi.org/10.1017/s0033822200033865>.

</div>

<div id="ref-buck_bayesian_1996" class="csl-entry">

Buck, C. E., W. G. Cavanagh, and C. D. Litton. 1996. *Bayesian Approach
to Interpreting Archaeological Data*. Chichester: John Wiley; Sons.

</div>

<div id="ref-buck_bcal_1999" class="csl-entry">

Buck, C. E., J. A. Christen, and G. N. James. 1999. “BCal: An On-Line
Bayesian Radiocarbon Calibration Tool.” *Internet Archaeology* 7.
<https://doi.org/10.11141/ia.7.1>.

</div>

<div id="ref-collins-elliott_lakhesis_underreview" class="csl-entry">

Collins-Elliott, S. A. Under Review. “Lakhesis: Consensus Seriation via
Iterative Regression of Partial Rankings for Binary Data.” *Journal of
Applied Statistics*, Under Review.

</div>

<div id="ref-crema_spatio-temporal_2017" class="csl-entry">

Crema, E. R., A. Bevan, and S. Shennan. 2017. “Spatio-Temporal
Approaches to Archaeological Radiocarbon Dates.” *Journal of
Archaeological Science* 87: 1–9.
<https://doi.org/10.1016/j.jas.2017.09.007>.

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
<https://doi.org/10.1016/b978-0-08-051581-6.50057-x>.

</div>

<div id="ref-hahsler_getting_2008" class="csl-entry">

Hahsler, M., K. Hornik, and C. Buchcta. 2008. “Getting Things in Order:
An Introduction to the R Package Seriation.” *Journal of Statistical
Software* 25: 1–34. <https://doi.org/10.18637/jss.v025.i03>.

</div>

<div id="ref-haslett_simple_2008" class="csl-entry">

Haslett, J., and A. C. Parnell. 2008. “A Simple Monotone Process with
Application to Radiocarbon-Dated Depth Chronologies.” *Journal of the
Royal Statistical Society: Series C (Applied Statistics)* 57: 399–418.
<https://doi.org/10.1111/j.1467-9876.2008.00623.x>.

</div>

<div id="ref-hinz_oxcaar_2021" class="csl-entry">

Hinz, M., C. Schmid, D. Knitter, and Tietze. 2021.
“<span class="nocase">oxcAAR</span>: Interface to ’OxCal’ Radiocarbon
Calibration.” <https://doi.org/10.32614/CRAN.package.oxcAAR>.

</div>

<div id="ref-jones_fixed-width_2006" class="csl-entry">

Jones, G. L., M. Haran, B. S. Caffo, and R. Neath. 2006. “Fixed-Width
Output Analysis for Markov Chain Monte Carlo.” *Journal of the American
Statistical Association* 101: 1537–47.
<https://doi.org/10.1198/016214506000000492>.

</div>

<div id="ref-lanos_hierarchical_2017" class="csl-entry">

Lanos, P., and A. Philippe. 2017. “Hierarchical Bayesian Modeling for
Combining Dates in Archeological Context.” *Journal de La Société
Française de Statistique* 158: 72–88.

</div>

<div id="ref-levy_chronological_2021" class="csl-entry">

Levy, E., G. Geeraerts, F. Pluquet, E. Piasetzky, and A. Fantalkin.
2021. “Chronological Networks in Archaeology: A Formalised Scheme.”
*Journal of Archaeological Science* 127: 105225.
<https://doi.org/10.1016/j.jas.2020.105225>.

</div>

<div id="ref-oksanen_vegan_2024" class="csl-entry">

Oksanen, J., G. L. Simpson, F. G Blanchet, R. Kindt, P. Legendre, P. R.
Minchin, R. B. O’Hara, et al. 2024. “Vegan: Community Ecology Package.”
<https://doi.org/10.32614/CRAN.package.vegan>.

</div>

<div id="ref-reimer_intcal20_2020" class="csl-entry">

Reimer, P. J., W. E. N. Austin, E. Bard, A. Bayliss, P. G. Blackwell, C.
Bronk Ramsey, M. Butzin, et al. 2020. “The IntCal20 Northern Hemisphere
Radiocarbon Age Calibration Curve (0–55 Cal
<span class="nocase">kBP</span>).” *Radiocarbon* 62: 725–57.
<https://doi.org/10.1017/RDC.2020.41>.

</div>

</div>
