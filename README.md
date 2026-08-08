
<!-- README.md is generated from README.Rmd. Please edit that file -->

# <img src="man/figures/logo.png" align="right" width="150px"/> eratosthenes: Archaeological Synchronism

<!-- badges: start -->

<a href="https://joss.theoj.org/papers/75ab124a2cbb3b7125a9458650544020"><img src="https://joss.theoj.org/papers/75ab124a2cbb3b7125a9458650544020/status.svg"></a>
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
from that full, joint conditional distribution. `Rcpp` is used for
faster sampling (Eddelbuettel and Balamuta 2018).

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
  sequences.
- discretization of time into intervals. Samples are drawn along the
  continuum.
- overly cumbersome chronological relationships. As `eratosthenes`
  samples points along the continuum, there is only before and after. If
  desired, overlaping events can be expressed in sequence construction:
  e.g., for sequences $A \prec B \prec C$ and
  $A \prec D \prec E \prec C$, events $B$ and $D , E$ will overlap with
  each other.

The focus of the package is on the structure of the joint conditional,
rather than specific probability models. Hence, `eratosthenes` relies on
the continuous uniform for estimating relative events. Any model can
however be used for absolute constraints, from single points to
customized densities.

The package is named after Eratosthenes of Cyrene, author of the
*Chronographiai*.

## Installation

To install `eratosthenes` from CRAN, run the `install.packages()`
function:

``` r
install.packages("eratosthenes")
```

To obtain the current development version of `eratosthenes` from GitHub,
install the package in the `R` command line with `remotes`:

``` r
library(remotes)
install_github("scollinselliott/eratosthenes", dependencies = TRUE, build_vignettes = TRUE) 
```

## Tutorial: Archaeological Example

The Dressel 1B type of amphora (a ceramic shipping container typical of
the ancient Mediterranean) is a ceramic type defined on the basis of
morphology, comprising a tall, long-necked, two-handled vessel with a
concave collared rim and sharply defined, angular shoulder. To give just
three chronological summaries of this type, the Dressel 1B type was
produced/used/deposited:

- in the “last quarter of the second until the last decade of the first
  century BC” (Southampton 2014);
- “shortly before the middle of the 1st century BC” and “until c. 10 BC”
  (Tyers 1996, 2.2);
- at the earliest “during *c*.100-80 BC”, disappearing “by 30 BC”
  (Loughton 2014, 43).

Such determinations are the result of many comparisons of stratigraphic
and single contexts where the Dressel 1B type has appeared in
conjunction with absolute constraints (e.g., coinage, datable stamps,
historical events in the occupation/destruction of sites). It is evident
in the above statements, too, that there will be greater or lesser
chronological discrepancies among investigators.

In order to obtain a probability density for any aspect of dating the
Dressel 1B (that is, its production, use, and/or deposition),
`eratosthenes` avoids the awkward synthesis of this chronological
information by estimating dates directly from archaeological record
itself, in the information provided by contexts (in
**[sequences](#sequences)**, howsoever established),
**[finds](#finds)**, which are contained in those contexts, and
**[absolute constraints](#absolute-constraints)** for those contexts, as
inputs. In other words, it formalizes the logic of dating finds which is
already in practice, yielding a probability density instead of a
qualitative appraisal of the object’s dates of production, use, and
deposition.

The following aims to be a concrete tutorial that gives step-by-step
instructions to obtain probability density functions on the production,
use, and deposition of the Dressel 1B type using `eratosthenes`. Further
information on the data for this tutorial is found
**[here](https://volweb.utk.edu/~scolli46/eratosthenes/eda20250628.html)**.

### Creating a Sequences Object

There are two fundamental objects for sequences:

- `events`, which are a sequence of events (e.g., archaeological
  contexts) in order from left (earliest) to right (latest)
- `sequences`, which are a collection (list) of `events`

For example, the `events` object `seq1` below gives a sequence of just
two depositional contexts, from the archaeological site of Rirha in
Morocco, `Rirha US 5182` and `Rirha US 5154`, with `Rirha US 5182` being
the earlier of the two:

``` r
library(eratosthenes)
seq1 <- events("Rirha US 5182", "Rirha US 5154") 
```

Multiple sequences can, and typically will, be given. To define seven
more sequences for this tutorial, we have additional deposits from the
site of Carthage (Byrsa Hill) and also several shipwrecks ([further
information
here](https://volweb.utk.edu/~scolli46/eratosthenes/eda20250628.html)):

``` r
seq2 <- events("Byrsa II B 19.4", "Byrsa II B 19.2")
seq3 <- events("Rirha US 5154", "Planier A")
seq4 <- events("El Sec", "Filicudi F", "Tour Fondue", "Cabrera 2", "Tour d'Agnello", "Sanguinaires A", "Lazaret")
seq5 <- events("Madrague de Giens", "Planier C", "Cap Béar C", "Planier A")
seq6 <- events("Cabrera 2", "Grand Congloué A", "Lazaret", "Byrsa II B 19.2", "Punta Scaletta", "Isla Pedrosa", "Cavalière", "Madrague de Giens")
seq7 <- events("Mazotos", "El Sec")
seq8 <- events("Grand Congloué A", "Héliopolis B", "Punta Scaletta")
```

A `sequences` object is then created which contains all of the `events`:

``` r
contexts <- sequences(seq1, seq2, seq3, seq4, seq5, seq6, seq7, seq8)
```

The `sequences()` function automatically checks that all `events` accord
with one another. Both `events()` and `sequences()` will return error
messages if the input is not valid (e.g., an `event` object contains
duplicate entries, or multiple `events` have conflicting orderings). In
such cases, `seq_diag()` function can be used to run a diagnostic check
on a collection of `events` to see which events are in disagreement.

``` r
seq_invalid1 <- events("Madrague de Giens", "El Sec", "Cabrera 2")
contexts_invalid <- list(seq1, seq2, seq3, seq4, seq5, seq6, seq_invalid1, seq7, seq8)
seq_diag(contexts_invalid)
```

The `seq_diag()` proceeds iteratively through all events, and so if
there are seqeunces of events which are assured to be valid, those
should be placed earlier in the inputs. The `shuffle` argument in
`seq_diag()` randomly permutes the order in which `events` are checked,
and can be called to check discrepancies if there is no information
about which `events` are more valid than others.

In sum, the `sequences` object contains all of the information about the
relative relationships of the contexts, in terms of which come before or
after one another, which should be based on a stated rationale
(typically, sequences from a seriation, stratigraphic sequences, or even
hypotheses).

### Creating a Finds Object

Next, a separate object for the finds data is created. If a particular
find has absolute chronological information associated with it (the
object itself, not the type), it should not be created here, but rather
as an `absolute` constraint (on which see below). To start, each find is
a `finds` object with the following information:

``` r
id1 <- finds(id = "id 1",
             assoc = "Isla Pedrosa",
             type = "AMPH Dressel 1B",
             residual = TRUE)
```

Each find, indexed with an `id`, must be linked to a context given in
the sequences above. Its appertaining context is entered into the
`assoc` field. This particular entry comprises the find of Dressel 1B
amphorae associated with the Isla Pedrosa shipwreck. These amphorae may
in fact relate to another shipwreck (i.e., their association with this
context may be spurious). By indicating `residual = TRUE`, the find
remains associated with its context, but it will not be taken into
account when estimating its date. This option applies whether finds are
residual to a deposit (i.e., much earlier finds redeposited into a later
deposit) or whether they are intrusive (i.e., much later finds have been
associated with an earlier deposit). The `type` field provides
information on any aspects of artifact classification or typology: here,
the type `AMPH Dressel 1B` is entered. It should be noted that `type` is
optional, and more than one entry can be given for `type`, as the next
four types show:

``` r
id865 <- finds(id = "id 865", assoc = "Rirha US 5154", type = "AMPH Dressel 1B")
id1202 <- finds(id = "id 1202", assoc = "Madrague de Giens", type = "AMPH Dressel 1B")
id1285 <- finds(id = "id 1285", assoc = "Planier C", type = "AMPH Dressel 1B")
id1364 <- finds(id = "id 1364", assoc = "Cap Béar C", type = c("AMPH Dressel 1B", "AMPH Dressel 1B Tarraconensis"))
```

The last find, `id1364`, belongs to the production group of Dressel 1B
amphorae produced in Spain, and so indicating its production subtype as
another element in the `type` fields ensures that the find and its
associated context will be used in estimating dates if an investigator
is obtaining dates for either `AMPH Dressel 1B` or
`AMPH Dressel 1B Tarraconensis`.

Finally, a single `assemblage` object is then created which contains all
of the finds:

``` r
finds <- assemblage(id1, id865, id1202, id1285, id1364)
```

### Creating Absolute Constraints Objects

Absolute constraints may be either *termini post quos* (*t.p.q.*) or
*termini ante quos* (*t.a.q.)*, which are created as `absolute` objects
(see [Absolute Constraints](#absolute-constraints) below). The structure
of an `absolute` object is the same as that of a `finds` object,
omitting the `"residual"` field, but requiring numeric input for the
field `"samples"`, which provide the absolute calendrical dates with
which the `absolute` object is associated. Just as `finds` are contained
within a single `assemblage` object, `absolute` objects are contained
with `constraints`.

For example *termini ante quos* for this tutorial include the
destruction of Carthage in 146 BCE, which the deposit B 19.2 on Byrsa
Hill predates, as well as a range of conventional absolute dates, ca.
1-15 CE, for the Planier A shipwreck (to give a finite endpoint for this
tutorial).

``` r
Carthage_Destr_1 <- absolute(id = "Carthage_Destr_1",
                         assoc = "Byrsa II B 19.2",
                        samples = -146)
Planier_A_abs <- absolute(id = "Planier_A_abs",
                      assoc = "Planier A",
                      samples = seq(1, 15, length.out = 100))

taq_info <- constraints(Carthage_Destr_1, Planier_A_abs)
```

The *termini post quos* are created in the same way. However, *t.p.q.*
will often feature calibrated radiocarbon (cal. r.c.) dates, whose
`samples` argument will consist of thousands of draws from their
probability density function. The most expedient way to convey
`"samples"` for cal. r.c. dates is to store the information on the
uncalibrated mean and standard deviation, and then to use a package to
generate the calibrated dates,

Here, a `csv` file is used to store information on uncalibrated r.c.
dates, which may be found in the `inst/extdata` folder of the package.
This selection of dates comes from Callegarin et al. (2016, 41) and
Manning, Lorentzen, and Demesticha (2022). The script below uses
`Bchron` (Haslett and Parnell 2008), looping over the table of
uncalibrated dates in order to draw samples from the calibrated date’s
p.d.f. and insert them into a list of *t.p.q.*, with an `id` and
associated context (`assoc`). Note that a `constraints` object may also
be applied to a `list` of `absolute` objects, such that additional
*t.p.q.* can be appended after the cal. r.c. dates:

``` r
# create an empty list to contain all tpq
tpq_info <- list()

eratosthenes_rcdates <- read.csv('inst/extdata/rc20250628.csv')

library(Bchron)

# calibrate and insert rc dates into tpq by looping over rows
for (i in 1:nrow(eratosthenes_rcdates)) {
    calib <- BchronCalibrate(ages = eratosthenes_rcdates$mu[i],
                             ageSds = eratosthenes_rcdates$sigma[i],
                             calCurves = "intcal20")
    x <- 1950 - sampleAges(calib)

    tpq_info[[i]] <- absolute(id = eratosthenes_rcdates$id[i],
                     assoc = eratosthenes_rcdates$assoc[i],
                     samples = x)
}

tpq_info <- constraints(tpq_info)
```

With the inputs of these sequences, finds, and absolute constraints, we
can estimate the dates of the production, use, and deposition of the
Dressel 1B type by calling the `gibbs_ad_type()` function from
`eratosthenes`. The `contexts` object containing the sequences is
entered first and the `finds` object second. We then specify what `type`
we want to obtain dates for (here, `AMPH Dressel 1B`). Finally, the
absolute constraints of `tpq_info` and `taq_info` are given as inputs:

``` r
dr_1b <- gibbs_ad_type(contexts,
                       finds,
                       type = "AMPH Dressel 1B",
                       tpq = tpq_info,
                       taq = taq_info)
```

The resulting `dr_1b` object contains samples of dates (production, use,
and depositional dates are estimated separately), and the console output
will give information on the structure of that object as well as a table
listing the mean dates. The densities can be visualized using the
`histogram()` function:

``` r
histogram(dr_1b, xlim = c(-300,20), ylim = c(0, 0.010), legend = "topleft") 
```

Which shows the probable dates of production, use, and deposition
separately (see [Fig. 3 of the *JOSS*
paper](https://raw.githubusercontent.com/openjournals/joss-papers/joss.08559/joss.08559/10.21105.joss.08559.pdf)
awaiting review), indicating which dates are more probable than others.

To determine a highest density region (HDR) of the dates, i.e., the
range of the most probable dates according to a given percentage, we can
construct a histogram of the dates with 100 bins (i.e., as a percentage)
and sort them by their density, choosing, for example, the 10% and 95%
regions of the most probable dates of the production of the type.

``` r
dat <- dr_1b$type$production
perc <- seq(min(dat), max(dat), length.out = 101)

x <- hist(dat, perc)$mids
y <- hist(dat, perc)$density

hpd <- rev(x[order(y)])
hpd_10 <- c(min(hpd[1:10]), max(hpd[1:10]))
hpd_95 <- c(min(hpd[1:95]), max(hpd[1:95]))
```

It should be noted that if the resulting distribution is multimodal,
this code will need to be modified to identify multiple regions. But
here, with a simple unimodal case, the lower and upper bounds or either
a narrower (10% HDR) or broader (95% HDR) interval of dates for the
production of a Dressel 1B type amphorae are, approximately:

    > hpd_10
    [1] -103.7998  -70.4462
    > hpd_95
    [1] -340.980962    7.378864

The dates of production will be earlier than those of deposition, which
can be seen in the histogram and by selecting `dr_1b$type$deposition`
rather than `dr_1b$type$production`.

The estimates naturally depend on the inputs, and this tutorial has used
only a small portion of the available information on the contexts and
constraints pertaining to the Dressel 1B type. The estimates given by
`eratosthenes` roughly accord with that of the conventional typology, as
would be expected, but it gives investigators the benefit of numerical
precision afforded by probability theory. Chronological typologies
depend on any number of conditional statements about depositional
contexts, their similarity to one another, and their sequencing. New
information is always being added, and older views are always being
reassessed. Hence, `eratosthenes` allows for a more rapid reassessment
of dates, as can be shown next.

### Chronological Interventions and Revisions

Revising chronologies, or evaluating the impact of choices upon or
interventions within a chronology, is performed directly on the inputs.
If we reconsider the view that the Dressel 1B is part of the Isla
Pedrosa wreck, we can set `residual = FALSE` and reassign the `id1`
object, reassign the finds object, and re-run the `gibbs_ad_type()`
function:

``` r
id1 <- list(id = "id 1", assoc = "Isla Pedrosa", type = "AMPH Dressel 1B", residual = FALSE)
finds <- list(id1, id865, id1202, id1285, id1364)
dr_1b <- gibbs_ad_type(contexts, finds, type = "AMPH Dressel 1B", tpq = tpq_info, taq = taq_info)
```

Performing the same steps as above to estimate the 10% and 95% HDR of
the dates of production, we can see that associating the Dressel 1B
amphorae with the cargo of the Isla Pedrosa wreck has substantially
changed the allocation of probability for what the most probable dates
are, shifting them earlier:

    > hpd_10
    [1] -137.09992  -99.68446
    > hpd_95
    [1] -342.884905    8.820347

Having a probablistically determined date for artifacts that works
directly from the basis of their archaeological relationships evades the
need to account for chronological discrepancies, and moreover affords
the ability to work entirely within a formal, mathematical environment
when it comes to topics such as the quantification of finds with
uncertain dating (since each find will have its own particular density
function). Futher tools to evaluate the influence of events on one
another for determining dates within in the chronology are discussed
below (see [Evaluating Displacement](#evaluating-displacement)).

## Usage

The basic objects in `eratosthenes` are:

- `sequences` of relative `events`, typically stratigraphic deposits,
  but also isolated contexts such as may be part of a frequency or
  contextual seriation
- an `assemblage` of `finds`, elements which belong to those events,
  typically artifacts
- `absolute` `constraints`, as either *termini post* or *ante quos*,
  expressed as samples from a probability density

Information related to these three items must be formatted in objects of
a `list` class, as follows.

### Sequences

Each relative sequences should run in order from left (earliest) to
right (latest), created as an `events` object. A `sequences` object is
then created from the `events`. In the following example, the object
`contexts` is created using `sequences()`, containing three sequences of
`events`.

``` r
x <- events("A", "B", "C", "D", "E", "F", "G", "H", "I", "J")
y <- events("B", "D", "G", "H", "K")
z <- events("F", "K", "L", "M")
contexts <- sequences(x, y, z)
```

If `sequences()` returns an error, it will likely be due to two or more
`events` having a conflicting ordering. See the section [Evaluating
Sequences](#evaluating-sequences) below.

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
f1 <- finds(id = "find01", assoc = "D", type = c("type1", "form1"))
f2 <- finds(id = "find02", assoc = "E", type = c("type1", "form2"))
f3 <- finds(id = "find03", assoc = "G", type = c("type1", "form1"), residual = TRUE)
f4 <- finds(id = "find04", assoc = "H", type = c("type2", "form1"))
f5 <- finds(id = "find05", assoc = "I", type = "type2")
f6 <- finds(id = "find06", assoc = "H", type = NULL)
artifacts <- assemblage(f1, f2, f3, f4, f5, f6)
```

### Absolute Constraints

Constraints should be given as two separate `lists`, one for *termini
post quos* and the other for *termini ante quos*. These take the same
form as the finds object, as a `list` of `lists`, with the same
headings, but include one additional heading of `samples` which contains
the absolute dates pertinent to that *t.p.q.* or *t.a.q*.

``` r
coin1 <- absolute(id = "coin1", assoc = "B", type = NULL, samples = runif(100, -320, -300))
coin2 <- absolute(id = "coin2", assoc = "G", type = NULL, samples = runif(100, 37, 41))
destr <- absolute(id = "destr", assoc = "J", type = NULL, samples = 79)

tpq_info <- constraint(coin1, coin2)
taq_info <- constraint(destr)
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
“type” has flexibility. While only one “type” at a time can be estimated
with `gibbs_ad_type()`, here, a “type” can be defined on the basis of:

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
seriations or ordinations, since R packages such as `seriation`
(Hahsler, Hornik, and Buchcta 2008), `vegan` (Oksanen et al. 2024),
`boral` (Hui 2016), `ecoCopula` (Popovic, Hui, and Warton 2022), `VGAM`
(Yee 2004), and `lakhesis` (Collins-Elliott Under Review) can perform
this task already.

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

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

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

<div id="ref-callegarin_rirha_2016-2" class="csl-entry">

Callegarin, L., M. Kbiri Alaoui, A. Ichkhakh, and J.-C. Roux, eds. 2016.
*Rirha : Site Antique Et Médiéval Du Maroc II. Période Maurétanienne (Ve
Siècle Av. J.-C. - 40 Ap. J.-C.)*. Collection de La Casa de Velázquez
151. Madrid: Casa de Velázquez.

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

<div id="ref-eddelbuettel_extending_2018" class="csl-entry">

Eddelbuettel, D., and J. J. Balamuta. 2018. “Extending R with C++: A
Brief Introduction to Rcpp.” *The American Statistician* 72: 28–36.
<https://doi.org/10.1080/00031305.2017.1375990>.

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

<div id="ref-hui_boral_2016" class="csl-entry">

Hui, F. K. C. 2016. “Boral: Bayesian Ordination and Regression Analysis
of Multivariate Abundance Data in R.” *Methods in Ecology and Evolution*
7: 744–50. <https://doi.org/10.1111/2041-210X.12514>.

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

<div id="ref-loughton_arverni_2014" class="csl-entry">

Loughton, M. 2014. *The Arverni and Roman Wine*. Oxford: Archaeopress.

</div>

<div id="ref-manning_dating_2022" class="csl-entry">

Manning, S. W., B. Lorentzen, and S. Demesticha. 2022. “Dating
Mediterranean Shipwrecks: The Mazotos Ship, Radiocarbon Dating and the
Need for Independent Chronological Anchors.” *Antiquity* 96: 968–80.
<https://doi.org/10.15184/aqy.2022.76>.

</div>

<div id="ref-oksanen_vegan_2024" class="csl-entry">

Oksanen, J., G. L. Simpson, F. G Blanchet, R. Kindt, P. Legendre, P. R.
Minchin, R. B. O’Hara, et al. 2024. “Vegan: Community Ecology Package.”
<https://doi.org/10.32614/CRAN.package.vegan>.

</div>

<div id="ref-popovic_fast_2022" class="csl-entry">

Popovic, G. C., F. K. C. Hui, and D. I. Warton. 2022. “Fast Model-Based
Ordination with Copulas.” *Methods in Ecology and Evolution* 13:
194–202. <https://doi.org/10.1111/2041-210X.13733>.

</div>

<div id="ref-reimer_intcal20_2020" class="csl-entry">

Reimer, P. J., W. E. N. Austin, E. Bard, A. Bayliss, P. G. Blackwell, C.
Bronk Ramsey, M. Butzin, et al. 2020. “The IntCal20 Northern Hemisphere
Radiocarbon Age Calibration Curve (0–55 Cal
<span class="nocase">kBP</span>).” *Radiocarbon* 62: 725–57.
<https://doi.org/10.1017/RDC.2020.41>.

</div>

<div id="ref-southampton_roman_2014" class="csl-entry">

Southampton, University of. 2014. “Roman Amphorae: A Digital Resource
\[Data-Set\].” *Internet Archaeology* 1.
<https://doi.org/10.5284/1028192>.

</div>

<div id="ref-tyers_roman_1996" class="csl-entry">

Tyers, P. 1996. “Roman Amphoras in Britain.” *Internet Archaeology* 1.
<https://doi.org/10.11141/ia.1.6>.

</div>

<div id="ref-yee_new_2004" class="csl-entry">

Yee, T. W. 2004. “A New Technique for Maximum-Likelihood Canonical
Gaussian Ordination.” *Ecological Monographs* 74: 685–701.
<https://doi.org/10.1890/03-0078>.

</div>

</div>
