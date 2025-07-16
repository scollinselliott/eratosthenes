---
title: 'eratosthenes: Synchronizing archaeological chronologies with a focus on artifact types'
tags:
  - R
  - archaeology
  - chronology
authors:
  - name: Stephen A. Collins-Elliott
    orcid: 0000-0002-5642-6903
    affiliation: 1
    corresponding: true
affiliations:
 - name: Department of Classics, University of Tennessee, Knoxville, Tennessee, United States
   index: 1
   ror: 020f3ap87
date: "30 June 2025"
bibliography: ../inst/REFERENCES.bib
---

# Summary

Archaeological chronology-building centers on two types of dates, absolute and relative. Absolute dates are associated with a calendrical date (even approximate), and relative refer to events only known via relationships to others (before/after). Relative events are then associated to absolute via contextual or logical associations. E.g, a _terminus post quem_ (_t.p.q._) of a datable coin found in a soil deposit ensures the deposit came after the production of the coin. Applying probabilistic methods to chronology has traditionally centered on radiocarbon dating, but the broader goal of extending formal methods to all aspects of chronology necessitate a way to address artifact typologies. Artifacts largely persist in being dated by qualitative judgment (e.g, "around the start of the 4th century BCE"). Moreover, typologies are in a constant state of revision and adjustment, such that the dates of types may differ according to authority. The reasons why a type is dated the way it is (i.e., which sites, deposits, and comparisons inform its date) can also become opaque, such that investigators inherently resist chronological revisions in order to assess its bibliography (on this conservatism, see @rotroff_four_2005[20]). This situation presents a challenge for any researcher applying statistical methods to datable artifacts, as one must make _ad hoc_ and often awkward decisions, resulting in what @lavan_checklist_2021[15] has called "abonimable dates."

The `R` package `eratosthenes` (named after Eratosthenes of Cyrene, author of the _Chronographiai_) provides functions for chronology-building, above all to bring artifact dating within the scope of formal mathematical estimation. Hence, an investigator can obtain separate probability density functions (p.d.f.) for the production, use, and deposition of an artifact type, in addition to marginal p.d.f.s for all relative sequential events (howsover determined) and absolute constraints (howsoever defined). It uses Gibbs sampling, by now a conventional Markov Chain Monte Carlo method in archaeological chronology [@geman_stochastic_1984; @buck_bayesian_1996]. Using `Rcpp` [@eddelbuettel_extending_2018], `eratosthenes` performs a two-step Gibbs routine, the first a preliminiary sampler to select a starting date, and a second main sampler that uses consistent batch means (CBM) as a stopping rule, given that convergence in distribution is assured [@jones_fixed-width_2006; @flegal_markov_2008]. Full reporting on the Monte Carlo standard errors (MCSE) is provided, giving an error in +/- years for each marginal density. Finally, `eratosthenes` provides functions for assessing the level of dependence of events upon each other within the joint conditional. Changes in dates brought about by any alteration in the structure of a chronology can therefore be readily evaluted. In sum, `eratosthenes` allows for expedient revision of chronologies, transparency in the definition of the full joint conditional, and statistics on the certainty of estimates via MCSE. 

# Statement of Need

Constructing chronologies via formal means is a major disciplinary goal, with a broad array of software developed to those ends. The [CRAN Task View: Archaeological Science](https://github.com/benmarwick/ctv-archaeology) maintained by Ben Marwick provides a section on chronological dating software in `R` [@r_core_team_r:_2024]. The calibration of $^{14}$C dates and the estimation of their posterior probability densities given constraints [@buck_bayesian_2025] are well served by `OxCal` [@bronk_ramsey_bayesian_2009] and `BCal` [@buck_bcal_1999], as well as `R` packages `oxcAAR` [@hinz_oxcaar_2021], `Bchron` [@haslett_simple_2008], `c14bazAAR` [@schmid_c14bazAAR_2019], and `rcarbon` [@crema_spatio-temporal_2017]. General chronological modeling is served by `Chronomodel` [@lanos_hierarchical_2017], which also centers on radiocarbon dating, and `ChronoLog` [@levy_chronological_2021], which focuses on discretized time intervals. The R package `ArchaeoPhases` [@philippe_analysis_2020] handles post-processing of samples form `OxCal`, `BCal`, and `ChronoModel`.  Other chronological software has focused on modeling count data over time, such as `kairos` [@frerebeau_kairos_2025] and  `baorista` [@crema_bayesian_2025], which focuses on counts of durative/interval events.

Given that there exist substantial software for calibrating radiocarbon dates, `eratosthenes` does not aim to perform this task. Likewise, the formal process of establishing relative sequences of contexts and finds via comparison, called seriation (or ordination), is a computationally difficult problem well served by many `R` packages, such as `seriation` [@hahsler_getting_2008], `vegan` [@oksanen_vegan_2024], `boral` [@hui_boral_2016], `ecoCopula` [@popovic_fast_2022], `VGAM` [@yee_new_2004], and `lakhesis` [@collins-elliott_lakhesis_underreview]. 

Two needs which `eratosthenes` aims to satisfy are (1) estimation of artifact dates (production, use, deposition) via formal means and (2) tools for evaluating the dependence of events on each other, which are discussed in detail below.

# Estimating Dates of Artifact Production, Use, and Deposition

Dating artifact types is handled by working backwards from depositional contexts, which is their primary point of observation. Estimating the date of artifact production requires a stated rule, with two options given in `eratosthenes`. To start, initial threshold boundaries are established between the earliest depositional context containing that artifact-type and the context immediately prior in sequence lacking it. The `"earliest"` rule samples from within those initial thresholds. The default, `"naive"` rule samples betwen a sample taken using the `"earliest"` rule as the lower bound and an upper bound formed by depositional date of that artifact.  After estimating production dates, a use date is sampled between its dates of production and deposition. \autoref{fig1} illustrates these rules.

![Rules of `earliest` and `naive` (default) for estimating the date of production of an artifact type, sampled on each $m$ iteration of the Gibbs sampler.\label{fig1}](fig1.png)

# Evaluating Displacement of Events

The impact of events on one another is assessed via "displacement," representing the shift in the date of a target event $i$ should another event $j$ be omitted from the joint conditional. Let $\tilde{x}_i$ be the estimated marginalized Monte Carlo mean date using all events within the joint conditional, and let $\tilde{x}_i^{(-j)}$ be the Monte Carlo mean when event $j$ has been omitted. Squared displacement of $j$ upon $i$ is defined as

$$
\delta^2(i,j) = (\tilde{x}_i^{(-j)} - \tilde{x}_i)^2.
$$

The omission of $j$ changes the date of $i$ by $\sqrt{\delta^2}$ amount of time (typically years). Conversely, the influence of an event $j$ upon all others is mesaured by mean squared displacement (MSD), which is the mean of all squared displacements when event $j$ has been omitted. Where $\Theta$ represents the set of all relative and absolute events, MSD is defined as

$$
\text{MSD}(j) = \frac{1}{n-1} \sum_{i \in \Theta, i \neq j} \delta^2 (i,j).
$$

`eratosthenes` computes the MSD for all events via a "jackknife" or "leave one out"-style of routine, omitting each event iteratively.

# Application

For the purposes of illustration, an initial application was undertaken for a small sample of events in the Mediterranean from the last four centuries BCE (data available [here](https://github.com/scollinselliott/eratosthenes-data/tree/main/data/20250628)). Results are shown in \autoref{fig2} for a selection of two depositional contexts and five shipwrecks. One can note, for example, the effect of using the sack of Carthage in 146 BCE as a _terminus ante quem_ (_t.a.q._) for the depositional date of context Byrsa II B 19.2 at that site [@lancel_byrsa_1982, 194].

The sack of Carthage has also been a key point for dating a particular type of ceramic transport container, the [Dressel 1 amphora type](https://archaeologydataservice.ac.uk/archives/view/amphora_ahrb_2005/details.cfm?id=324). Since it has not be found in pre-destruction layers, the start of its production has been dated to the later part of the 2nd century BCE [@tchernia_vin_1986,42]. Rather than just relying on one site, however, absence from any context within the entire set of conditional events in `eda20250628` yields densities for its production, use, and deposition (\autoref{fig3}). 

# Future Work

Current work in progress by the author which relies upon `eratosthenes` involves the synchronism of ceramics, coinage, radiocarbon dates, depositional/seriated sequences, and historical events for the central Mediterranean in the last four centuries BCE, with datasets uploaded to the GitHub repository [eratosthenes-data](https://github.com/scollinselliott/eratosthenes-data).

![Marginal p.d.f.s of 2 depositional events and 5 shipwrecks from the Mediterraenan, given the joint conditional contained in [eda20250628](https://github.com/scollinselliott/eratosthenes-data/tree/main/data/20250628). The wreck Grand Congloué A is earlier than the traditional date: future datasets will work to revise sequencing and constraints.\label{fig2}](fig2.png)

![P.d.f.s of the production, use, and deposition of the Dressel 1B type amphora, using the "naive" production rule and given the joint conditional contained in [eda20250628](https://github.com/scollinselliott/eratosthenes-data/tree/main/data/20250628).\label{fig3}](fig3.png)

# References
