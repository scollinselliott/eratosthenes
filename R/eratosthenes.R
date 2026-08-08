# #' Sequence Check
# #'
# #' For a \code{list} of partial sequences (of \code{vector} objects), check to see that joint elements of each occur the same order. That is, for two sequences with elements \eqn{A, B, C, D, E} and \eqn{B, D, F, E}, all joint elements must occur in the same order to pass the check. Two sequences \eqn{A, B, C, D, E} and \eqn{A, F, D, C, E} would not pass this check as the elements \eqn{C} and \eqn{D} occur in different orders in either sequence.
# #' 
# #' Event names \code{alpha} and \code{omega} are reserved for the ultimate boundaries of the chronological framework and cannot be used in naming events in sequences. This function is automatically performed when creating a \fcode{sequences} object (see \code{\link[eratosthenes]{sequences}}).
# #' 
# #' @param obj A \code{list} of \code{vector} objects which represent a sequence.    
# #' @examples 
# #' x <- events("A", "B", "C", "D", "E")
# #' y <- events("B", "D", "F", "E")
# #' 
# #' seq_check(x, y)
# #' 
# #' z <- events("B", "F", "C")
# #' 
# #' seq_check(x, y, z)
# #' 
# #' @returns \code{TRUE} or \code{FALSE}
# #' 
# #' @export
# seq_check <- function(...) {
#     UseMethod("seq_check")
# }

# #' @rdname seq_check
# #' @export
# seq_check.events <- function(...) {
#     qp_ <- quae_postea(...)

#     check <- TRUE
#     for (i in names(qp_)) {
#         if (i != "omega") {
#             if (i %in% qp_[[i]]) {
#                 check <- FALSE
#             }
#         }
#     }

#     return(check)
# } 



#' Synthetic Ranking
#'
#' Using a \code{sequences} object of two or more partial sequences, all of which observe the same order of elements, create a single "synthetic" ranking. This is accomplished by counting the total number of elements after running a recursive trace through all partial sequences (via \code{\link[eratosthenes]{quae_postea}}). If partial sequences are inconsistent in their rankings, a \code{NULL} value is returned.
#'
#' @param obj A \code{sequences} object.
#' @param ties The way in which ties are handled per the \code{\link{rank}} function. The default is \code{"ties = average"}.
#' 
#' @examples 
#' x <- events("A", "B", "C", "D", "E")
#' y <- events("B", "D", "F", "E")
#' a <- sequences(x, y)
#' 
#' synth_rank(a)
#' 
#' @returns An \code{events} object containing the synthesized ranking.
#' 
#' @export
synth_rank <- function(obj, ties = "average") {
    UseMethod("synth_rank")
}

#' @rdname synth_rank
#' @export
synth_rank.sequences <- function(obj, ties = "average") {
    res <- NULL
    lens <- sapply(obj, length)
    if (any(lens < 2)) {
        stop("events in input sequences must contain two or more elements", call. = FALSE)
    }

    elements <- names(obj)
    qp_ <- quae_postea(obj)

    quot_postea <- numeric(length(elements))
    names(quot_postea) <- elements
    for (i in names(qp_)) {
        quot_postea[i] <- length(qp_[[i]])
    }
    res <- rank(quot_postea * -1, ties.method = ties)
    res <- names(res)[order(res)]

    class(res) <- c("events", "character")
    return(res)
}



#' Quae Postea
#'
#' For a \code{list} of multiple partial sequences (of \code{vector} objects), generate another \code{list} which, for each element, gives all elements that occur after it ("\emph{quae postea}"). This is analogous to a recursive trace through all partial sequences from left to right. A final element \code{"omega"} is added to all sets to avoid empty vectors. See also \code{\link[eratosthenes]{quae_antea}}.
#'
#' @param ... Objects of class \code{\link[eratosthenes]{events}}, or a \code{\link[eratosthenes]{sequences}} object, a valid \code{list} of \code{events}.
#' 
#' @examples 
#' x <- events("A", "B", "C")
#' y <- events("B", "D", "E", "C", "F")
#' z <- events("C", "G")
#' 
#' quae_postea(x)
#' quae_postea(x, y, z)
#' 
#' a <- sequences(x, y, z)
#' quae_postea(a)
#' 
#' @returns A \code{list} of \code{vector} objects, which contain the elements that occur after any one given element in the input sequences. 
#' 
#' @export
quae_postea <- function(...) {
    UseMethod("quae_postea")
}
#' 
#' @rdname quae_postea
#' @export
quae_postea.events <- function(...) {
    obj <- list(...)
    if (0 %in% vapply(obj, length, 1L)) {
        stop('input contains NULL element', call. = FALSE)
    }
    elements <- unique(unlist(obj))
    M <- list()
    for (i in 1:length(obj)) {
        tmp <- numeric(length(obj[[i]]))
        for (j in 1:length(obj[[i]])) {
            tmp[j] <- which(elements == obj[[i]][j])
        }
        M[[i]] <- tmp
    }

    mat <- quae_postea_matrix_cpp(length(elements), M)

    res <- list()
    for (i in 1:length(elements)) {
        res[[elements[i]]] <- c(elements[mat[i,] == 1], "omega")
    }

    check <- TRUE
    for (i in names(res)) {
        if (i != "omega") {
            if (i %in% res[[i]]) {
                check <- FALSE
            }
        }
    }

    if (check == FALSE) {
        stop('conflicts in element orders for one or more events objects.\nRun seq_diag() for conflicting pairs of events.', call. = FALSE)
    }

    return(res)
}
#' @rdname quae_postea
#' @export
quae_postea.list <- function(...) {
    obj <- list(...)[[1]]
    chk <- sapply(obj, inherits, "events")
    if (FALSE %in% chk) {
        stop("non-events object in list input.")
    }
        
    elements <- unique(unlist(obj))
    M <- list()
    for (i in 1:length(obj)) {
        tmp <- numeric(length(obj[[i]]))
        for (j in 1:length(obj[[i]])) {
            tmp[j] <- which(elements == obj[[i]][j])
        }
        M[[i]] <- tmp
    }

    mat <- quae_postea_matrix_cpp(length(elements), M)

    res <- list()
    for (i in 1:length(elements)) {
        res[[elements[i]]] <- c(elements[mat[i,] == 1], "omega")
    }

    check <- TRUE
    for (i in names(res)) {
        if (i != "omega") {
            if (i %in% res[[i]]) {
                check <- FALSE
            }
        }
    }

    if (check == FALSE) {
        stop('conflicts in element orders for one or more events objects.\nRun seq_diag() for conflicting pairs of events.', call. = FALSE)
    }

    return(res)
}
#' 
#' @rdname quae_postea
#' @export
quae_postea.sequences <- function(...) {
    obj <- list(...)[[1]]
    if (0 %in% vapply(obj, length, 1L)) {
        stop('input contains NULL element', call. = FALSE)
    }
    elements <- unique(unlist(obj))
    M <- list()
    for (i in 1:length(obj)) {
        tmp <- numeric(length(obj[[i]]))
        for (j in 1:length(obj[[i]])) {
            tmp[j] <- which(elements == obj[[i]][j])
        }
        M[[i]] <- tmp
    }

    mat <- quae_postea_matrix_cpp(length(elements), M)

    res <- list()
    for (i in 1:length(elements)) {
        res[[elements[i]]] <- c(elements[mat[i,] == 1], "omega")
    }

    return(res)
}



#' Quae Antea
#'
#' For a \code{list} of multiple partial sequences of \code{events} objects, generates a \code{list} which, for each element, giving the elements that occur before it ("\emph{quae antea}"). This is analogous to a recursive trace through all partial sequences from right to left. An element \code{"alpha"} is added to all sets to avoid empty vectors. See also \code{\link[eratosthenes]{quae_postea}}. 
#'
#' @param ... Objects of class \code{\link[eratosthenes]{events}}, or a \code{\link[eratosthenes]{sequences}} object, a valid \code{list} of \code{events}.
#'
#' @examples 
#' x <- events("A", "B", "C")
#' y <- events("B", "D", "E", "C", "F")
#' z <- events("C", "G")
#' 
#' quae_antea(x, y, z)
#' 
#' @returns A \code{list} of \code{vector} objects, which contain the elements that occur before any one given element in the input sequences. 
#'
#' @export
quae_antea <- function(...) {
    UseMethod("quae_antea")
}
#' 
#' @rdname quae_antea
#' @export
quae_antea.events <- function(...) {
    obj <- list(...)
    if (0 %in% vapply(obj, length, 1L)) {
        stop('input contains NULL element', call. = FALSE)
    }
    elements <- unique(unlist(obj))
    M <- list()
    for (i in 1:length(obj)) {
        tmp <- numeric(length(obj[[i]]))
        for (j in 1:length(obj[[i]])) {
            tmp[j] <- which(elements == obj[[i]][j])
        }
        M[[i]] <- tmp
    }

    mat <- quae_antea_matrix_cpp(length(elements), M)

    res <- list()
    for (i in 1:length(elements)) {
        res[[elements[i]]] <- c(elements[mat[i,] == 1], "alpha")
    }

    check <- TRUE
    for (i in names(res)) {
        if (i != "alpha") {
            if (i %in% res[[i]]) {
                check <- FALSE
            }
        }
    }

    if (check == FALSE) {
        stop('conflicts in element orders for one or more events objects.\nRun seq_diag() for conflicting pairs of events.', call. = FALSE)
    }
    return(res)
}
#' 
#' @rdname quae_antea
#' @export
quae_antea.list <- function(...) {
    obj <- list(...)[[1]]
    if (0 %in% vapply(obj, length, 1L)) {
        stop('input contains NULL element', call. = FALSE)
    }
    chk <- sapply(obj, inherits, "events")
    if (FALSE %in% chk) {
        stop("non-events object in list input.")
    }
    elements <- unique(unlist(obj))
    M <- list()
    for (i in 1:length(obj)) {
        tmp <- numeric(length(obj[[i]]))
        for (j in 1:length(obj[[i]])) {
            tmp[j] <- which(elements == obj[[i]][j])
        }
        M[[i]] <- tmp
    }

    mat <- quae_antea_matrix_cpp(length(elements), M)

    res <- list()
    for (i in 1:length(elements)) {
        res[[elements[i]]] <- c(elements[mat[i,] == 1], "alpha")
    }

    check <- TRUE
    for (i in names(res)) {
        if (i != "alpha") {
            if (i %in% res[[i]]) {
                check <- FALSE
            }
        }
    }

    if (check == FALSE) {
        stop('conflicts in element orders for one or more events objects.\nRun seq_diag() for conflicting pairs of events.', call. = FALSE)
    }
    return(res)
}
#'
#' @rdname quae_antea
#' @export
quae_antea.sequences <- function(...) {
    obj <- list(...)[[1]]
    elements <- unique(unlist(obj))
    M <- list()
    for (i in 1:length(obj)) {
        tmp <- numeric(length(obj[[i]]))
        for (j in 1:length(obj[[i]])) {
            tmp[j] <- which(elements == obj[[i]][j])
        }
        M[[i]] <- tmp
    }

    mat <- quae_antea_matrix_cpp(length(elements), M)

    res <- list()
    for (i in 1:length(elements)) {
        res[[elements[i]]] <- c(elements[mat[i,] == 1], "alpha")
    }

    return(res)
}



#' Adjust Sequence to Target
#'
#' Given an "input" sequence of elements and another "target" seqeunce that contains fewer elements in a different order, shift the order of the input sequence to match that of the target, keeping all other elements as proximate to one another as possible. This adjusted ranking is accomplished using piecewise linear interpolation between joint elements ranks. That is, joint rankings are plotted, with input rankings along the \eqn{x} axis and target rankings on the \eqn{y} axis. Remaining rankings in the input sequence are assigned a ranking of \eqn{y} based on the piecewise linear function between joint rankings. If the rank order of elements in the target are identical to those in the input, the result is identical to the input. A minimum number of three joint elements in both the input and target are required.
#' 
#' @param input An \code{events} object of unique ordered elements. 
#' @param target An \code{events} object of unique ordered elements.  containing at least three of the same elements as \code{input}.
#' 
#' @examples 
#' x <- events("A", "B", "C", "D", "E", "F", "G", "H", "I", "J") # the input sequence of events
#' y <- events("D", "A", "J") # the target sequence of events
#' 
#' seq_adj(x, y)
#' 
#' @returns An \code{events} object of the adjusted sequence.
#' 
#' @export
seq_adj <- function(input, target) {
    UseMethod("seq_adj")
}
#' 
#' @rdname seq_adj
#' @export
seq_adj.events <- function(input, target) {
    res <- NULL
    joint <- intersect(input, target)
    if (length(joint) > 2) {
        xj <- input %in% joint
        x_pos <- 1:length(input)
        y_pos <- 1:length(target)
        names(x_pos) <- input
        names(y_pos) <- target
        x <- x_pos[xj]
        y <- y_pos[names(x)]
        x <- c(0, x, length(x_pos) + 1)
        y <- c(0, y, length(y_pos) + 1)
        interp <- stats::approx(x, y, n = length(x_pos) + 2)
        res <- interp$y[1:length(x_pos)+1]
        names(res) <- input
        res <- input[order(res)]
    } else {
        stop("Insufficient number of joint elements in input and target sequence (must be > 2).")
    }
    class(res) <- c("events", "character")
    return(res)
}



#' Gibbs Sampler for Archaeological Dates
#'
#' A Gibbs sampler for dating archaeological events, to fit relative sequences to absolute, calendrical dates. Relative events can be associated with \emph{termini post quos} (\emph{t.p.q.}) and \emph{termini ante quos} (\emph{t.a.q.}), which are entered as samples from a given probability density function \eqn{f(t)}. This function may take any form, a single date (i.e., with a probability of 1), a continuous uniform distribution (any time between two dates), or a bespoke density (as with calibrated radiocarbon dates). Relative events are modeled on a continuous uniform density between the latest antecedent event and earliest subsequent event.
#' 
#' Gibbs sampling is a conventional method for calibrating and estimating radiocarbon dates in light of absolute constraints and relative sequences: see \insertCite{buck_bayesian_1996,buck_bcal_1999,bronk_ramsey_bayesian_2009;textual}{eratosthenes}.
#' 
#' In this implementation, two phases of Gibbs sampling are performed: an initial phase for selecting starting values and then the main sampler, with convergence evaluated using Monte Carlo standard errors (MCSE).
#' 
#' The initial Gibbs sampler results in a vector of starting values randomly sampled for each event up to \eqn{\sqrt{k}} runs, where \eqn{k} is the total number of events. Starting values may therefore take some time to assign, but this initial sampling is necessary to avoid a catastrophic collapse due to floating point errors in the initial selection of random values and will also result in closer starting values with respect to marginal densities. 
#' 
#' The main Gibbs sampler uses consistent batch means (CBM) determine convergence and hence when to end the main sampling run: there is no motivation to remove burn-in from the main sampling run nor to run multiple chains. CBM is assured to converge in distribution, see \insertCite{jones_fixed-width_2006,flegal_markov_2008;textual}{eratosthenes}. A stopping point for the main sampler is therefore determined using the mean of the Monte Carlo standard errors (MCSE) across all random variates, which is the input of \code{mcse_crit} (the mean MCSE for all events). The input \code{max_samples} indicates the maximum number of simulations to run, but the sampler will stop if the specified criterion of \code{mcse_crit} is passed. The default mean MCSE is set at \code{mcse_crit = 0.5}, as the MCSE is measured in years (i.e. to allow for an error +/- 1 year), but, to be sure, individual events will have higher or lower MCSE than this mean criterion, whose primary purpose is as a stopping rule. 
#' 
#' Note that the MCSE criterion is applied as a stopping rule for depositional dates and external constraints. The number of Monte Carlo samples for production dates of types is chosen to be identical to that need to pass \code{mcse_crit}, such that ultimately the final mean MCSE of all variates may differ from that of the criterion. Depending on the conditional structure of the relative sequences and the timescale of investigation, higher or lower MCSE may be more desirable or acceptable.
#' 
#' For the use dates of artifact type production, use, and deposition, see the \code{\link[eratosthenes]{gibbs_ad_type}} function.
#'
#' @param sequences A \code{\link[eratosthenes]{sequences}} object of relative sequences of elements (e.g., contexts).
#' @param max_samples Maximum number of samples to run. Default is \code{10^5}.
#' @param size The number of samples to take on each iteration of the main Gibbs sampler. Default is \code{10^3}. 
#' @param mcse_crit Criterion for the Monte Carlo standard error to stop the Gibbs sampler, as based on depositional dates and absolute constraints. The number of Monte Carlo samples for production dates is identical to that depositional dates.
#' @param tpq A \code{\link[eratosthenes]{constraints}} object containing all \emph{termini post quos}.
#' @param taq A \code{\link[eratosthenes]{constraints}} object containing all \emph{termini ante quos}.
#' @param alpha_ An initial \emph{t.p.q.} to limit any elements which may occur before the first provided \emph{t.p.q.} Default is \code{-5000}.
#' @param omega_ A final \emph{t.a.q.} to limit any elements which may occur after the after the last provided \emph{t.a.q.} Default is \code{1950}.
#' @param trim A logical value to determine whether elements that occur before the first \emph{t.p.q.} and after the last \emph{t.a.q.} should be omitted from the results (i.e., to "trim" elements at the ends of the sequence, whose marginal densities depend on the selection of \code{alpha_} and \code{omega_}). Default is \code{TRUE}.
#' 
#' @returns A \code{list} object of class \code{marginals} which contains the following:
#'    * \code{deposition} A \code{list} of samples from the marginal density of each context's depositional date.
#'    * \code{externals} A \code{list} of samples of the marginal density of each constraint (\emph{t.p.q.} and \emph{t.a.q.]}), as conditioned upon the occurrence of other depositional 
#'    * \code{production} If a \code{finds} object has been input, samples of the marginal density of the production date of finds types will be included in the output. If types are attested in trimmed contexts, 
#'    * \code{mcse} The Monte Carlo standard errors (MCSE) of the random variates (fixed t.p./a.q. will have a MCSE of 0.)
#'
#' @examples
#' x <- events("A", "B", "C", "D", "E", "F", "G", "H", "I", "J")
#' y <- events("B", "D", "G", "H", "K")
#' z <- events("F", "K", "L", "M")
#' contexts <- sequences(x, y, z)
#'  
#' # external constraints
#' coin1 <- absolute(id = "coin1", assoc = "B", type = NULL, samples = runif(100,-320,-300))
#' coin2 <- absolute(id = "coin2", assoc = "G", type = NULL, samples = seq(37, 41, length = 100))
#'   # seq(37, 41, length = 100) is equivalent in concept to runif(100, 37, 41)) 
#' destr <- absolute(id = "destr", assoc = "J", type = NULL, samples = 79)
#' 
#' tpq_info <- constraints(coin1, coin2)
#' taq_info <- constraints(destr)
#' 
#' result <- gibbs_ad(contexts, tpq = tpq_info, taq = taq_info)
#' 
#' @references
#'   \insertAllCited{}
#' 
#' @export
#' @importFrom Rdpack reprompt
gibbs_ad <- function(sequences, max_samples = 10^5, size = 10^3, mcse_crit = 0.5, tpq = NULL, taq = NULL, alpha_ = -5000, omega_ = 1950, trim = TRUE) {
    UseMethod("gibbs_ad")
}

#' @rdname gibbs_ad
#' @export
gibbs_ad.sequences <- function(sequences, max_samples = 10^5, size = 10^3, mcse_crit = 0.5, tpq = NULL, taq = NULL, alpha_ = -5000, omega_ = 1950, trim = TRUE) {
    if (!(is.null(tpq) | inherits(tpq, "constraints"))) {
        stop("input tpq must be constraints object. See constraints().")
    }
    if (!(is.null(taq) | inherits(taq, "constraints"))) {
        stop("input taq must be constraints object. See constraints().")
    }
    if (size > max_samples) {
        stop("Error: size must be less than max_samples.")
    }
    # if (seq_check(sequences) == FALSE) {
    #     stop("Sequences has failed consistency check with seq_check().")
    # }

    proceed <- synth_rank(sequences)

    if (!inherits(tpq, "constraints")) {
        tpq <- constraints(absolute(id = "tpq_default", assoc = proceed[1], type = NULL, samples = alpha_))
    }
    if (!inherits(taq, "constraints")) {
        taq <- constraints(absolute(id = "taq_default", assoc = proceed[length(proceed)], type = NULL, samples = omega_))
    }

    # proceed_all from tpq, taq, relative, alpha, omega
    proceed_all <- c()

    # total number of elements
    elements <- length(proceed) + 2
    if (!is.null(tpq)) {
        elements <- elements + length(tpq)
    }
    if (!is.null(taq)) {
        elements <- elements + length(taq)
    }

    # indices
    if (is.null(tpq) & is.null(taq)) {
        proceed_idx <- 1:length(proceed) 
    } else if (is.null(tpq) & !is.null(taq)) {
        taq_idx <- 1:length(taq)
        proceed_idx <- (1:length(proceed)) + length(taq)
    } else if (!is.null(tpq) & is.null(taq)) {
        tpq_idx <- 1:length(tpq)
        proceed_idx <- (1:length(proceed)) + length(tpq)
    } else {
        tpq_idx <- 1:length(tpq)
        taq_idx <- (length(tpq) + 1):(length(tpq) + length(taq))
        proceed_idx <-  (1:length(proceed)) + (length(tpq) + length(taq))
    }

    # # proceed_all from tpq, taq, relative, alpha, omega
    # proceed_all <- c()

    # # total number of elements
    # elements <- length(tpq) + length(taq) + length(proceed) + 2

    # # indices
    # tpq_idx <- 1:length(tpq)
    # taq_idx <- (length(tpq) + 1):(length(tpq) + length(taq))
    # proceed_idx <-  (1:length(proceed)) + (length(tpq) + length(taq))

    gibbs <- matrix(0, nrow = elements, ncol = size)

    for (i in 1:length(sequences)) {
        sequences[[i]] <- c("alpha", sequences[[i]], "omega")   
    }

    j <- 0
    for (i in 1:length(tpq)) {
        if (!(tpq[[i]]$assoc %in% proceed)) {
            stop(paste0("Context of tpq ", tpq[[i]]$id, " : ", tpq[[i]]$assoc, " is not given in relative sequences"))
        }
        sequences <- c(sequences, list(c(tpq[[i]]$id, tpq[[i]]$assoc))  )
        proceed_all <- c(proceed_all, tpq[[i]]$id )
        gibbs[j+1,1] <- min(tpq[[i]]$samples)          # initialize tpq with earliest possible
        j <- j + 1
    }
    for (i in 1:length(taq)) {
        if (!(taq[[i]]$assoc %in% proceed)) {
            stop(paste0("Context of taq ", taq[[i]]$id, " : ", taq[[i]]$assoc, " is not given in relative sequences"))
        }
        sequences <- c(sequences, list(c(taq[[i]]$assoc, taq[[i]]$id)) )
        proceed_all <- c(proceed_all, taq[[i]]$id )
        gibbs[j+1,1] <- max(taq[[i]]$samples)          # initialize taq with latest possible
        j <- j + 1
    }

    proceed_all <- c(proceed_all, proceed, "alpha", "omega")
    gibbs[elements-1,1] <- alpha_
    gibbs[elements,1] <- omega_

    # convert to indices
    M <- list()
    for (i in 1:length(sequences)) {
        tmp <- numeric(length(sequences[[i]]))
        for (j in 1:length(sequences[[i]])) {
            tmp[j] <- which(proceed_all == sequences[[i]][j])
        }
        M[[i]] <- tmp
    }

    PhiMatrix <- quae_antea_matrix_cpp(elements, M)
    PsiMatrix <- quae_postea_matrix_cpp(elements, M)

    init_sample <- floor(sqrt(elements))

    # indices of non-trimmed relative events
    if (trim == TRUE) {
        idx_nontrim <- numeric(length(proceed))
        for (i in 1:length(proceed)) {
            idx <- proceed_idx[i]   
            check1 <- sum(PsiMatrix[idx, taq_idx])
            check2 <- sum(PhiMatrix[idx, tpq_idx])
            if (check1 > 0 & check2 > 0) {
                idx_nontrim[i] <- 1
            }
        } 
        idx_trim <- which(idx_nontrim == 0) + length(tpq) + length(taq)
        trim_label <- proceed[which(idx_nontrim == 0)]
        nontrim_label <- proceed_all[!(names(proceed_all) %in% trim_label)]
    }

    # initial sampler
    message("Assigning initial random values (this may take a moment)...")

    gibbs[,1] <- gibbs_ad_initial_cpp(gibbs[,1], tpq_idx, PsiMatrix, tpq, taq_idx, PhiMatrix, taq, proceed_idx, init_sample)

    # main sampler
    message("Beginning main Gibbs sampler. Will terminate either when MCSE criterion or maximum number of MC samples reached.")

    gibbs <- gibbs_ad_cpp(gibbs, tpq_idx, PsiMatrix, tpq , taq_idx, PhiMatrix, taq, proceed_idx)

    # consistent batch means
    mcse_check <- FALSE
    while (mcse_check == FALSE) {
        n_upto <- ncol(gibbs)
        n_batch <- floor(sqrt(n_upto)) # length of samples in batch
        K <- floor(n_upto / n_batch) # number of batches
        m_batch <- matrix(NA, nrow = nrow(gibbs), ncol = (K-1))
        remainder <- n_upto - n_batch * K + 1

        idx1 <- remainder

        for (k in 1:(K-1)) {
            idxs <- idx1:(idx1 + n_batch)       

            # in cases where n_upto = n_batch * K
            idxs <- idxs[idxs <= ncol(gibbs)]

            m_batch[, k] <- rowMeans(gibbs[ , idxs]) 
            idx1 <- idxs[length(idxs)] + 1
        }
        
        mcse0 <- sqrt( rowSums( (m_batch - rowMeans(m_batch))^2 ) * (n_batch / (K-1) ) ) / sqrt((K - 1) * n_batch)
        
        # remove trimmed events from estimation of mean MCSE
        mcse <- mcse0
        if (trim == TRUE) {
            if (length(idx_trim) > 0) {
                mcse <- mcse[-idx_trim]
            }
        }
        # do not include fixed single-point events as part of estimating MCSE
        mcse <- mcse[mcse > 0] 

        cat("\r", paste0("Samples: ", ncol(gibbs), "     Mean MCSE: ",  round(mean(mcse),3)))
        if (mean(mcse) < mcse_crit) {
            mcse_check <- TRUE
            message("\nMCSE criterion passed. Finishing.")
        } else {
            if (ncol(gibbs) >= max_samples) {
                message("\nMC samples exceeded maximum stipulated without passing MCSE criterion. Finishing.")
                mcse_check <- TRUE
            } else {

                gibbs_next <- matrix(0, nrow = nrow(gibbs), ncol = (size + 1) )
                gibbs_next[,1] <- gibbs[,ncol(gibbs)]
                gibbs_next <- gibbs_ad_cpp(gibbs_next, tpq_idx, PsiMatrix, tpq , taq_idx, PhiMatrix, taq, proceed_idx)
                gibbs <- cbind(gibbs, gibbs_next[, 2:ncol(gibbs_next)])
            }
        }
    }

    deposition <- list()
    externals <- list()

    names(mcse0) <- proceed_all

    for (i in 1:length(proceed)) {
        iname <- proceed[i]
        idx <- proceed_idx[i]   
        if (trim == TRUE) {
            check <- TRUE
            check1 <- sum(PsiMatrix[idx, taq_idx])
            check2 <- sum(PhiMatrix[idx, tpq_idx])
            if (check1 > 0 & check2 > 0) {
                g <- gibbs[idx, ]
                deposition[[iname]] <- g
            }
        } else {
            g <- gibbs[idx, ]
            deposition[[iname]] <- g
        }
    }
    for (i in 1:length(tpq)) {
        ii <- tpq[[i]]
        tpq_idx <- match(ii$id, proceed_all)
        g <- gibbs[tpq_idx, ]
        externals[[ii$id]] <- g
    }
    for (i in 1:length(taq)) {
        ii <- taq[[i]]
        taq_idx <- match(ii$id, proceed_all)
        g <- gibbs[taq_idx, ]
        externals[[ii$id]] <- g
    }

    res <- list(deposition = deposition, externals = externals, mcse = mcse0)
    class(res) <- c("marginals", "list")
    return(res)

}
#' 
#' @rdname gibbs_ad
#' @export
gibbs_ad.list <- function(sequences, max_samples = 10^5, size = 10^3, mcse_crit = 0.5, tpq = NULL, taq = NULL, alpha_ = -5000, omega_ = 1950, trim = TRUE) {
    sequences <- sequences(sequences)
    gibbs_ad.sequences(sequences, max_samples = 10^5, size = 10^3, mcse_crit = 0.5, tpq = NULL, taq = NULL, alpha_ = -5000, omega_ = 1950, trim = TRUE)
}


# #' @export
# print.events <- function(x, ...) {
#     x <- sort(as.vector(names(x)))
#     cat("\n Events (call with $ operator to retrieve samples):\n")
#     cat(" ", x, "\n\n\n")
# }



#' @export 
print.marginals <- function(x, ...) {
    mcse <- x$mcse
    mcse <- mcse[c(names(x$deposition), names(x$externals))]
    mmcse <- mean(mcse[mcse > 0])
    cat("\n Marginals from joint conditional density, consisting of:\n",
    "    ", length(x$deposition), "depositional events\n", 
    "    ", length(x$externals), "external constraints (t.p./a.q.)\n\n",
        "Call the following objects to see element names: \n",
        "     $deposition", "\n",
        "     $externals", "\n\n",
        "Mean Monte Carlo standard error (MCSE) of depositional/external variates: \n",
    "    ", round(mmcse, 3), "\n\n\n")
}



#' @export 
print.type_marginals <- function(x, ...) {
    stat_disp <- round(x$stat, 3) 
    cat("\n Marginals from joint conditional density of", x$name, 
    "\n consisting of type's production, deposition, and use.\n\n",
    "Call the following objects for samples: \n",
    "     $type$production\n",
    "     $type$use\n",
    "     $type$deposition\n\n",
    "Monte Carlo mean and s.e.:\n")
    print.data.frame(stat_disp)
    cat("\n")
}



#' @export
summary.marginals <- function(object, events = NULL, digits = 2, ...) {
    depmu <- sapply(object$deposition, mean)
    depmcse <- object$mcse[names(object$deposition)]
    depdat <- data.frame(Mean = depmu, MCSE = depmcse)
    rownames(depdat) <- names(object$deposition)

    extmu <- sapply(object$externals, mean)
    extmcse <- object$mcse[names(object$externals)]
    extdat <- data.frame(Mean = extmu, MCSE = extmcse)
    rownames(extdat) <- names(object$externals)

    prdmu <- sapply(object$production, mean)
    prdmcse <- object$mcse[names(object$production)]
    prddat <- data.frame(Mean = prdmu, MCSE = prdmcse)
    rownames(prddat) <- names(object$production)

    if (is.null(events)) {
        dat <- rbind(depdat, extdat, prddat)
    } else if (events == "deposition") {
        dat <- depdat
    } else if (events == "externals") {
        dat <- extdat
    } else if (events == "production") {
        dat <- prddat
    } else {
        stop("'events' should be NULL (all), deposition, externals, or production.")
    }

    dat <- round(dat, digits)

    message("\nMonte Carlo mean and standard errors for ", events)

    return(dat)
}



#' @export
summary.type_marginals <- function(object, ...) {
    stat_disp <- round(object$stat, 3) 
    cat("\n Marginals from joint conditional density of", object$name, 
    "\nconsisting of type's production, deposition, and use.\n",
    "Call the following objects for samples s: \n",
    "     $type$production\n",
    "     $type$use\n",
    "     $type$deposition\n\n",
    "Monte Carlo meean and s.e.:\n")
    return(stat_disp) 
}



#' Traceplot of Gibbs Samples
#' 
#' Wrapper around \code{\link[graphics]{plot}} to make a traceplot of Gibbs samples from \code{\link[eratosthenes]{gibbs_ad}}. See \code{\link[eratosthenes]{histogram}}. Maximum number of simulatenous events to display is 12.  for plotting a density histogram of events.
#' 
#' See also \code{\link[eratosthenes]{tidy_marginals}} for exporting the results of these functions into tidy data frame for custom plotting in e.g., \code{ggplot2}.
#' 
#' @param x A \code{list} object of class \code{marginals}, the output of \code{\link[eratosthenes]{gibbs_ad}}.
#' @param events A vector or element of the event names to plot. Maximum number of events is 12.
#' @param xlim The limits of the x-axis (optional).
#' @param ylim The limits of the y-axis (optional).
#' @param xlab Label for the y-axis. Default is \code{"Index"}.
#' @param ylab Label for the y-axis. Default is \code{"Year"}.
#' @param palette A vector providing the color palette of the histogram. The default is \code{"colorBlindness::paletteMartin"} (see \code{\link[paletteer]{palettes_d}}).
#' @param opacity The opacity/transparency of the traceplot, if visualizing overlapping events. A value between 0 and 1 (default).
#' @param legend_pos The position of the legend in the plot. Default is \code{"topright"}.
#' 
#' @examples 
#' x <- events("A", "B", "C", "D", "E", "F", "G", "H", "I", "J")
#' y <- events("B", "D", "G", "H", "K")
#' z <- events("F", "K", "L", "M")
#' contexts <- sequences(x, y, z)
#'  
#' # external constraints
#' coin1 <- absolute(id = "coin1", assoc = "B", type = NULL, samples = runif(100,-320,-300))
#' coin2 <- absolute(id = "coin2", assoc = "G", type = NULL, samples = seq(37, 41, length = 100))
#' destr <- absolute(id = "destr", assoc = "J", type = NULL, samples = 79)
#' 
#' tpq_info <- constraints(coin1, coin2)
#' taq_info <- constraints(destr)
#' 
#' result <- gibbs_ad(contexts, tpq = tpq_info, taq = taq_info)
#' 
#' traceplot(result, "B")
#' traceplot(result, c("coin1", "B", "H"), opacity = 0.5)
#' 
#' @returns A traceplot of the Gibbs samples of the selected events.
#' 
#' @export
traceplot <- function(x, events = NULL, xlim = NULL, ylim = NULL, xlab = "Index", ylab = "Year", palette = NULL, opacity = 1, legend_pos = "topright") {
    UseMethod("traceplot")
}
#' 
#' @rdname traceplot
#' @export
traceplot.marginals <- function(x, events = NULL, xlim = NULL, ylim = NULL, xlab = "Index", ylab = "Year", palette = NULL, opacity = 1, legend_pos = "topright") {
    if (!is.vector(events)) {
        stop('Events must be a single element or vector.')
    }
    if (is.null(palette)) {
        palette <- paletteer::paletteer_d("colorBlindness::paletteMartin")
    }

    palette <- grDevices::adjustcolor(palette, alpha.f = opacity)

    plot_list <- list()
    idx <- 1
    for (i in 1:length(events)) {
        if (events[i] %in% names(x$deposition)) {
            plot_list[[idx]] <- x$deposition[events[i]][[1]]
            idx <- idx + 1
        } else if (events[i] %in% names(x$externals)) {
            plot_list[[idx]] <- x$externals[events[i]][[1]]
            idx <- idx + 1
        } else {
            stop("One or more events not contained in marginals.")
        }
    }
    
    if (is.null(xlim)) {
        xlim <- c(0, length(plot_list[[1]]))
    }

    if (is.null(ylim)) {
        y_all <- unlist(plot_list)
        ylim <- c(min(y_all), max(y_all))
    }

    if (length(events) == 1) {
        x <- plot_list[[1]]
        graphics::plot(x, type = "l", col = palette[1], xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, main = "")
    } else if (length(events) > 1 & length(events) <= 12) {
        x <- plot_list[[1]]
        graphics::plot(x, type = "l", col = palette[1], xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, main = "")
        for (k in 2:length(events)) {
            x <- plot_list[[k]]
            graphics::lines(x, type = "l", col = palette[k])
        }
    } else {
        stop("Max number events for traceplot is 12.")
    }
    graphics::legend(legend_pos, legend= events, col = palette[1:length(events)], pt.cex = 2.5, pch = 15)
}



#' Histogram of Marginal Densities
#' 
#' Wrapper around \code{\link[graphics]{hist}} to plot density histograms for select marginal densities (up to 12) in a single plot, from the results of \code{\link[eratosthenes]{gibbs_ad}}, or to plot density histograms of the production, deposition, and use of a type, from the results of \code{\link[eratosthenes]{gibbs_ad_type}}.
#' 
#' See also also \code{\link[eratosthenes]{tidy_marginals}} for exporting the results of these functions into tidy data frame for custom plotting in e.g., \code{ggplot2}.
#' 
#' @param x A \code{list} object of class \code{marginals}, the output of \code{\link[eratosthenes]{gibbs_ad}}, or of class \code{type_marginals}, to plot the output of \code{\link[eratosthenes]{gibbs_ad_type}]}.
#' @param events If plotting a \code{marginals} object, a vector or element of the event names to plot. Maximum number of events is 12.
#' @param aspect If plotting a \code{type_marginals} object, that is, the output of \code{\link[eratosthenes]{gibbs_ad_type}}, a vector of one or more of \code{c("production", "use", "deposition")}. The default is all three.
#' @param breaks The number or method of breaks in the histogram. Default is \code{"Freedman-Diaconis"}. See \code{\link[graphics]{hist}} for more.
#' @param xlim The limits of the x-axis. Default is set to the min/max values of all samples.
#' @param ylim The limits of the y-axis. This may need to be adjusted if densities have an extremely narrow interval.
#' @param xlab Label for the x-axis. Default is \code{"Year"}.
#' @param palette A vector providing the color palette of the histogram. The default is \code{"colorBlindness::paletteMartin"} (see \code{\link[paletteer]{palettes_d}}).
#' @param opacity The opacity/transparency of the histograms for visualizing overlapping events, a value between 0 and 1 (default).
#' @param legend_pos The position of the legend in the plot. Default is \code{"topright"}.
#' @returns A density histogram of the selected events/aspects.
#' 
#' @examples 
#' x <- events("A", "B", "C", "D", "E", "F", "G", "H", "I", "J")
#' y <- events("B", "D", "G", "H", "K")
#' z <- events("F", "K", "L", "M")
#' contexts <- sequences(x, y, z)
#' 
#' f1 <- finds(id = "find01", assoc = "D", type = c("type1", "form1"))
#' f2 <- finds(id = "find02", assoc = "E", type = c("type1", "form2"))
#' f3 <- finds(id = "find03", assoc = "G", type = c("type1", "form1"), residual = TRUE)
#' f4 <- finds(id = "find04", assoc = "H", type = c("type2", "form1"))
#' f5 <- finds(id = "find05", assoc = "I", type = "type2")
#' f6 <- finds(id = "find06", assoc = "H", type = NULL)
#' 
#' artifacts <- assemblage(f1, f2, f3, f4, f5, f6)
#'  
#' # external constraints
#' coin1 <- absolute(id = "coin1", assoc = "B", type = NULL, samples = runif(100,-320,-300))
#' coin2 <- absolute(id = "coin2", assoc = "G", type = NULL, samples = seq(37, 41, length = 100))
#' destr <- absolute(id = "destr", assoc = "J", type = NULL, samples = 79)
#' 
#' tpq_info <- constraints(coin1, coin2)
#' taq_info <- constraints(destr)
#' 
#' result <- gibbs_ad(contexts, tpq = tpq_info, taq = taq_info)
#' 
#' # deposition of "B"
#' histogram(result, "B")
#' 
#' # deposition of "coin2" and deposition of "G"
#' histogram(result, c("coin2", "G"), opacity = 0.5)
#' 
#' # production, use, and deposition of "type1"
#' result_type1 <- gibbs_ad_type(contexts, artifacts, type = "type1",
#'                           max_samples = 3000, mcse_crit = 2)
#' histogram(result_type1)
#' 
#' @returns A density histogram of the selected events/aspect.
#' 
#' @export
histogram <- function(x, events = NULL, aspect = c("production", "use", "deposition"), breaks = "Freedman-Diaconis", xlim = NULL, ylim = NULL, xlab = "Year", palette = NULL, opacity = 1, legend_pos = "topright") {
    UseMethod("histogram")
}
#' 
#' @rdname histogram
#' @export
histogram.marginals <- function(x, events = NULL, aspect = NULL, breaks = "Freedman-Diaconis", xlim = NULL, ylim = NULL, xlab = "Year", palette = NULL, opacity = 1, legend_pos = "topright") {
    if (!is.vector(events)) {
        stop('Events must be a single element or vector.')
    }
    if (is.null(palette)) {
        palette <- paletteer::paletteer_d("colorBlindness::paletteMartin")
    }

    palette <- grDevices::adjustcolor(palette, alpha.f = opacity)

    plot_list <- list()
    idx <- 1
    for (i in 1:length(events)) {
        if (events[i] %in% names(x$deposition)) {
            plot_list[[idx]] <- x$deposition[events[i]][[1]]
            idx <- idx + 1
        } else if (events[i] %in% names(x$externals)) {
            plot_list[[idx]] <- x$externals[events[i]][[1]]
            idx <- idx + 1
        } else {
            stop("One or more events not contained in marginals.")
        }
    }

    if (is.null(xlim)) {
        x_all <- unlist(plot_list)
        xlim <- c(min(x_all), max(x_all))
    }

    if (length(events) == 1) {
        dat <- plot_list[[1]]
        graphics::hist(dat, breaks = breaks, freq = FALSE, xlim = xlim, xlab = xlab, ylim = ylim, col = palette[1], lty="blank", main = "")
    } else if (length(events) > 1 & length(events) <= 12) {
        dat <- plot_list[[1]]
        graphics::hist(dat, breaks = breaks, freq = FALSE, lty="blank", xlim = xlim, ylim = ylim, col = palette[1], xlab = xlab, main = "")
        for (k in 2:length(events)) {
            dat <- plot_list[[k]]
            graphics::hist(dat, breaks = breaks, freq = FALSE, lty="blank", col = palette[k], xlab = xlab, main = "", add = TRUE)
        }
    } else {
        stop("Max number events for histogram is 12.")
    }
    graphics::legend(legend_pos, legend = events, col = palette[1:length(events)], pt.cex = 2.5, pch = 15)
}
#' 
#' @rdname histogram
#' @export
histogram.type_marginals <- function(x, events = NULL, aspect = c("production", "use", "deposition"), breaks = "Freedman-Diaconis", xlim = NULL, ylim = NULL, xlab = "Year", palette = NULL, opacity = 0.5, legend_pos = "topright") {
    display_name <- x$name

    if ((!("production" %in% aspect) & !("use" %in% aspect) ) & !("deposition" %in% aspect))  {
        stop("At least one or more event types (production, use, deposition) must be selected.")
    }

    if (is.null(palette)) {
        palette <- paletteer::paletteer_d("colorBlindness::paletteMartin")
    }

    dat <- c(year = c(), aspect = c())

    dep_samples <- x$type$deposition
    prd_samples <- x$type$production
    use_samples <- x$type$use

    if ("production" %in% aspect) {
        dat <- rbind(dat, data.frame(year = prd_samples, aspect = paste0(display_name, " - Production")))
    }
    if ("use" %in% aspect) {
        dat <- rbind(dat, data.frame(year = use_samples, aspect = paste0(display_name, " - Use")))
    }
    if ("deposition" %in% aspect) {
        dat <- rbind(dat, data.frame(year = dep_samples, aspect = paste0(display_name, " - Deposition")))
    }

    aspect_name <- unique(dat$aspect)

    palette <- grDevices::adjustcolor(palette, alpha.f = opacity)
    if (is.null(xlim)) {
        xlim <- c(min(dat$year), max(dat$year))
    }

    if (length(aspect) == 1) {
        x <- dat$year
        graphics::hist(x, breaks = breaks, freq = FALSE, xlim = xlim, ylim = ylim, col = palette[1], xlab = xlab, lty="blank", main = "")
    } else {
        x <- dat$year
        graphics::hist(x, breaks = breaks, freq = FALSE, lty="blank", xlim = xlim, ylim = ylim, col = palette[1], xlab = xlab, main = "")
        for (k in 2:length(aspect_name)) {
            x <- dat[dat$aspect == aspect_name[k], ]$year
            graphics::hist(x, breaks = breaks, freq = FALSE, lty="blank", col = palette[k], xlab = xlab, main = "", add = TRUE)
        }
    }
    graphics::legend(legend_pos, legend = aspect_name, col = palette[1:length(aspect)], pt.cex = 2.5, pch = 15)
}



#' Convert Marginals to Tidy (Molten) Data Frame
#' 
#' Takes the results of \code{\link[eratosthenes]{gibbs_ad}} or \code{\link[eratosthenes]{gibbs_ad_type}} and "melts" the \code{list} into a tidy data frame \insertCite{wickham_tidy_2014}{eratosthenes}. Each row of the molten data frame will contain the index of the Monte Carlo sample, the sample itself, and then the event name.
#' 
#' @param input An object of class \code{marginals} or \code{type_marginals}, the output of \code{\link[eratosthenes]{gibbs_ad}} or \code{\link[eratosthenes]{gibbs_ad_type}}.
#' @returns A data frame giving the MC sampling index (\code{idx}), the sample (\code{year}), and the event (\code{event}).
#' 
#' @examples
#' x <- events("A", "B", "C", "D", "E", "F", "G", "H", "I", "J")
#' y <- events("B", "D", "G", "H", "K")
#' z <- events("F", "K", "L", "M")
#' contexts <- assemblage(x, y, z)
#' 
#' # external constraints
#' coin1 <- absolute(id = "coin1", assoc = "B", type = NULL, samples = runif(100,-320,-300))
#' coin2 <- absolute(id = "coin2", assoc = "G", type = NULL, samples = seq(37, 41, length = 100))
#' destr <- absolute(id = "destr", assoc = "J", type = NULL, samples = 79)
#' 
#' tpq_info <- constraints(coin1, coin2)
#' taq_info <- constraints(destr)
#' 
#' result <- gibbs_ad(contexts, tpq = tpq_info, taq = taq_info)
#' 
#' tidy_marginals(result)
#' 
#' @references
#'   \insertAllCited{}
#' 
#' @export
tidy_marginals <- function(input) {
    UseMethod("tidy_marginals")
}
#' 
#' @rdname tidy_marginals
#' @export
tidy_marginals.marginals <- function(input) {
    events <- c(names(input$deposition), names(input$externals), names(input$production))
    master <- c(input$deposition, input$externals, input$production)
    dat <- data.frame(year = 0, event = 0)
    for (i in 1:length(events)) {
        item <- events[i]
        if (!(item %in% names(master))) {
            stop('One or more event names not contained in marginals.')
        }
        dat <- rbind(dat, data.frame(year = master[[item]], event = item))
        dat <- dat[-1 , ]
        rownames(dat) <- NULL
        dat$event <- factor(dat$event)
    }  
    return(dat)
}
#' 
#' @rdname tidy_marginals
#' @export
tidy_marginals.type_marginals <- function(input) {
    dat <- data.frame(x = c(), event = c())

    dep_samples <- input$type$deposition
    prd_samples <- input$type$production
    use_samples <- input$type$use

    dat <- rbind(dat, data.frame(x = use_samples, event = paste0(input$name, " - Use")))
    dat <- rbind(dat, data.frame(x = dep_samples, event = paste0(input$name, " - Deposition")))
    dat <- rbind(dat, data.frame(x = prd_samples, event = paste0(input$name, " - Production")))
    dat$event <- factor(dat$event)

    return(dat)
}



# #' Convert Finds List Object to Data Frame (Context / Find-Type)
# #' 
# #' Performs the opposite of \code{\link[eratosthenes]{finds_d2l}}. Takes a \code{list} object of finds and their types, used as input in \code{\link[eratosthenes]{gibbs_ad}}, and returns a \code{data.frame} of two columns, containing the context in the first and the find-type in the second, and the \code{id} of the object in the third.
# #' 
# #' @param input A list object of finds (each one a list) of associated contexts and types. 
# #' @returns A four-column data frame of contexts (first column) and find-types attested in that context (second column), along with the id number (third column).
# #' 
# #' @examples 
# #' f1 <- list(id = "find01", assoc = "D", type = c("type1", "form1"))
# #' f2 <- list(id = "find02", assoc = "E", type = c("type1", "form2"))
# #' f3 <- list(id = "find03", assoc = "G", type = c("type1", "form1"), residual = TRUE)
# #' f4 <- list(id = "find04", assoc = "H", type = c("type2", "form1"))
# #' f5 <- list(id = "find05", assoc = "I", type = "type2")
# #' f6 <- list(id = "find06", assoc = "H", type = NULL)
# #' 
# #' artifacts <- list(f1, f2, f3, f4, f5, f6)
# #' 
# #' # convert list to data frame
# #' artifacts_df <- finds_l2d(artifacts)
# #' 
# #' @export
# finds_l2d <- function(input) {
#     UseMethod("finds_l2d")
# }
# #' 
# #' @rdname finds_l2d
# #' @export
# finds_l2d.list <- function(input) {
#     result <- data.frame(Context = c(), Type = c())
#     for (i in input) {
#         if (!is.na(sum(match(c("id", "assoc", "type"), names(i))))) {
#             if (length(i$type) > 0) {
#                 for (j in i$type) {
#                     result <- rbind(result, data.frame(Context = i$assoc, Type = j, Id = i$id ) )
#                 }
#             }
#         } else {
#             stop("Finds list object does not contain correct headings (id, assoc, type).")
#         }
#     }
#     return(result)
# }



# #' Convert Finds Data Frame (Context / Find-Type) to List Object
# #' 
# #' Performs the opposite of \code{\link[eratosthenes]{finds_l2d}}. Takes a \code{data.frame} object of two columns, containing the context in the first and the find-type in the second, and returns a \code{list} object for input in \code{\link[eratosthenes]{gibbs_ad}}. The value of the find \code{id} is automatically generated as an integer if not provided in a third column.
# #' 
# #' @param input A two-column data frame of contexts (first column) and find-types (second column). An optional third column of an id number may be provided.
# #' @returns A list of finds (each one a list) associated with contexts and their types.
# #' 
# #' @examples 
# #' f1 <- list(id = "find01", assoc = "D", type = c("type1", "form1"))
# #' f2 <- list(id = "find02", assoc = "E", type = c("type1", "form2"))
# #' f3 <- list(id = "find03", assoc = "G", type = c("type1", "form1"))
# #' f4 <- list(id = "find04", assoc = "H", type = c("type2", "form1"))
# #' f5 <- list(id = "find05", assoc = "I", type = "type2")
# #' f6 <- list(id = "find06", assoc = "H", type = NULL)
# #' 
# #' artifacts <- list(f1, f2, f3, f4, f5, f6)
# #' 
# #' # convert list to data frame
# #' artifacts_df <- finds_l2d(artifacts)
# #' 
# #' # convert data frame to list
# #' artifacts_list <- finds_d2l(artifacts_df)
# #' 
# #' @export
# finds_d2l <- function(input) {
#     UseMethod("finds_d2l")
# }
# #' 
# #' @rdname finds_d2l
# #' @export
# finds_d2l.data.frame <- function(input) {
#     result <- list()
#     if (ncol(input) == 2) {
#         for (i in 1:nrow(input)) {
#             result[[i]] <- list(id = as.character(i), assoc = input[i,1], type = input[i,2])
#         }
#     } else if (ncol(input) == 3) {
#         for (i in 1:nrow(input)) {
#             result[[i]] <- list(id = input[i,3], assoc = input[i,1], type = input[i,2])
#         }
#     } else {
#         stop("data frame must be 2 or three columns.")
#     }
#     return(result)
# }



#' Ids of Types
#' 
#' Given a \code{list} object of finds (with keys of \code{id}, \code{assoc}, \code{type} in each entry), return a vector of the \code{id} elements that belong to one or more specified type.
#' 
#' @param input An \code{\link[eratosthenes]{assemblage}}, comprising finds.
#' @param type A character vector or element
#' 
#' @examples 
#' f1 <- finds(id = "find01", assoc = "D", type = c("type1", "form1"))
#' f2 <- finds(id = "find02", assoc = "E", type = c("type1", "form2"))
#' f3 <- finds(id = "find03", assoc = "G", type = c("type1", "form1"))
#' f4 <- finds(id = "find04", assoc = "H", type = c("type2", "form1"))
#' f5 <- finds(id = "find05", assoc = "I", type = "type2")
#' f6 <- finds(id = "find06", assoc = "H", type = NULL)
#' 
#' artifacts <- assemblage(f1, f2, f3, f4, f5, f6)
#' 
#' ids_of_types(artifacts, type = "type1")
#' ids_of_types(artifacts, type = c("type1", "type2"))
#' 
#' @returns A character vector of ids within a \code{list} object of \code{finds} class, 
#' 
#' @export
ids_of_types <- function(input, type = NULL) {
    UseMethod("ids_of_types")
}
#' 
#' @rdname ids_of_types
#' @export
ids_of_types.assemblage <- function(input, type = NULL) {
    res <- c()
    if (is.character(type)) {
        for (k in type) {
            for (i in input) {
                if (length(i$type) > 0) {
                    for (j in i$type) {
                        if (j == k) {
                            res <- c(res, i$id)
                        }
                    }
                }
            }
        }
    if (is.null(res)) {
        stop("No matches for specified type.")
    }
    return(unique(res))
    } else {
        stop("type is not vector object.")
    }
}
            


#' Gibbs Sampler for Archaeological Dates: Artifact Types
#'
#' Estimate a densities for the production, use, and deposition dates of an artifact or artifact type. Multiple artifacts and types can be given, which will be pooled into a single type. For example, one can input several individual finds via their id number as comprising a type, or multiple (sub)types/classes as a single type, (e.g., "MGS V amphora", "MGS VI amphora", "MGS V/VI amphora" to construct one group). Depending on whether one is using id numbers or type(s), the \code{id} or \code{type} argument is used, which takes a vector of the entries' names. The \code{gibbs_ad_type} function works on the basis of the presence/absence of types in contexts, sampling a use date between the production and deposition. The stipulation of a rule to determine production dates (\code{naive} or \code{earliest}) is required.
#' 
#' See \code{\link[eratosthenes]{gibbs_ad}} for information on consistent batch means and Monte Carlo standard error, which are used to determined convergence for the use date.
#'
#' @param sequences A \code{\link[eratosthenes]{sequences}} object of relative sequences of elements (e.g., contexts).
#' @param finds  An \code{\link[eratosthenes]{assemblage}} object of finds related to (contained in) the elements of \code{sequences}.
#' @param id A vector of the \code{id} of one or more specific finds whose use date is to be estimated. The values of \code{id} must match those in the \code{list} of \code{finds}. If \code{type} is used, \code{id} is ignored.
#' @param type A vector of one or more types to estimate a use density for. Must contain a value if \code{id} is \code{NULL}.
#' @param type_name A customized label for the type (e.g., if one is selecting via \code{id} or has combined subtypes). If only \code{type} is used to select finds, the default will be that label Otherwise the default is simply "Type."
#' @param max_samples Maximum number of samples to run. Default is \code{10^5}.
#' @param size The number of samples to take on each iteration of the main Gibbs sampler. Default is \code{10^3}. 
#' @param mcse_crit Criterion for the Monte Carlo standard error to stop the Gibbs sampler. Only the MCSE of the use date is used as a stopping rule.
#' @param tpq A \code{\link[eratosthenes]{constraints}} object containing all \emph{termini post quos}.
#' @param taq A \code{\link[eratosthenes]{constraints}} object containing all \emph{termini ante quos}.
#' @param alpha_ An initial \emph{t.p.q.} to limit any elements which may occur before the first provided \emph{t.p.q.} Default is \code{-5000}.
#' @param omega_ A final \emph{t.a.q.} to limit any elements which may occur after the after the last provided \emph{t.a.q.} Default is \code{1950}.
#' @param trim A logical value to determine whether elements that occur before the first \emph{t.p.q.} and after the last \emph{t.a.q.} should be omitted from the results (i.e., to "trim" elements at the ends of the sequence, whose marginal densities depend on the selection of \code{alpha_} and \code{omega_}). Default is \code{TRUE}.
#' @param rule The rule for computing an estimated date of production of a find-type, either \code{"earliest"}, selecting a production date between the earliest deposition of that type and the next most earliest context, or \code{"naive"} (the default), which will select a production date any time between the distribution of that "earliest" date and the depositional date of that artifact.
#' 
#' @examples 
#' x <- events("A", "B", "C", "D", "E", "F", "G", "H", "I", "J")
#' y <- events("B", "D", "G", "H", "K")
#' z <- events("F", "K", "L", "M")
#' contexts <- sequences(x, y, z)
#' 
#' f1 <- finds(id = "find01", assoc = "D", type = c("type1", "form1"))
#' f2 <- finds(id = "find02", assoc = "E", type = c("type1", "form2"))
#' f3 <- finds(id = "find03", assoc = "G", type = c("type1", "form1"), residual = TRUE)
#' f4 <- finds(id = "find04", assoc = "H", type = c("type2", "form1"))
#' f5 <- finds(id = "find05", assoc = "I", type = "type2")
#' f6 <- finds(id = "find06", assoc = "H", type = NULL)
#' 
#' artifacts <- assemblage(f1, f2, f3, f4, f5, f6)
#'  
#' # external constraints
#' coin1 <- absolute(id = "coin1", assoc = "B", type = NULL, samples = runif(100,-320,-300))
#' coin2 <- absolute(id = "coin2", assoc = "G", type = NULL, samples = seq(37, 41, length = 100))
#' destr <- absolute(id = "destr", assoc = "J", type = NULL, samples = 79)
#' 
#' tpq_info <- constraints(coin1, coin2)
#' taq_info <- constraints(destr)
#' 
#' # use dates by specifying ids
#' gibbs_ad_type(contexts, artifacts, id = c("find04", "find05"),
#'               max_samples = 2000, mcse_crit = 2, tpq = tpq_info, taq = taq_info)
#'
#' # use dates by specifying types
#' gibbs_ad_type(contexts, artifacts, type = "type1",
#'               max_samples = 2000, mcse_crit = 2, tpq = tpq_info, taq = taq_info)
#' 
#' @returns A \code{list} of class \code{type_marginals} of the density of a use date, conditional upon production and depositional dates.
#' 
#' @export
gibbs_ad_type <- function(sequences, finds = NULL, id = NULL, type = NULL, type_name = NULL, max_samples = 10^5, size = 10^3, mcse_crit = 0.5, tpq = NULL, taq = NULL, alpha_ = -5000, omega_ = 1950, trim = TRUE, rule = "naive") {
    UseMethod("gibbs_ad_type")
}
#' 
#' @rdname gibbs_ad_type
#' @export
gibbs_ad_type.sequences <- function(sequences, finds = NULL, id = NULL, type = NULL, type_name = NULL, max_samples = 10^5, size = 10^3, mcse_crit = 0.5, tpq = NULL, taq = NULL, alpha_ = -5000, omega_ = 1950, trim = TRUE, rule = "naive") {
    if (!(is.null(finds) | inherits(finds, "assemblage"))) {
        stop("finds must be NULL or assemblage object. See assemblage().")
    }
    if (!(is.null(tpq) | inherits(tpq, "constraints"))) {
        stop("input tpq must be constraints object. See constraints().")
    }
    if (!(is.null(taq) | inherits(taq, "constraints"))) {
        stop("input taq must be constraints object. See constraints().")
    }
    if (is.null(id) & is.null(type)) {
        stop("Either one or more id or types must be specified.")
    } else {
        if (!is.null(type)) {
            if (!is.null(id)) {
                message('Defaulting to type (omit or specify "type = NULL" to use id).')
            }
            id <- ids_of_types(finds, type)
        }
    }
    
    contexts_type <- c()
    for (j in finds) {
        if (j$id %in% id) {
            if (is.null(j$residual)) {
                contexts_type <- c(contexts_type, j$assoc)
            } else {
                if (!(j$residual == TRUE)) {
                    contexts_type <- c(contexts_type, j$assoc)
                }
            }
        }
    }

    tpq_production <- c()

    for (j in tpq) {
        if (j$id %in% id) {
            tpq_production <- c(tpq_production, j$id)
            #contexts_type <- c(contexts_type, j$assoc)
        }
    }
    # for (j in taq) {
    #     if (j$id %in% id) {
    #         contexts_type <- c(contexts_type, j$assoc)
    #     }
    # }

    # contexts that contain the ids/types
    contexts_type <- unique(contexts_type)

    if (length(contexts_type) == 0) {
        stop("No ids or types associated with events found.")
    }

    message("Estimating production, use, and depositional dates for id(s)/type(s) specified.")

    # if (seq_check(sequences) == FALSE) {
    #     stop("Sequences has failed consistency check with seq_check().")
    # }

    proceed <- synth_rank(sequences)

    if (!is.list(tpq)) {
        tpq <- list(list(id = "tpq_default", assoc = proceed[1], type = NULL, samples = alpha_))
    }
    if (!is.list(taq)) {
        taq <- list(list(id = "taq_default", assoc = proceed[length(proceed)], type = NULL, samples = omega_))
    }

    # proceed_all from tpq, taq, relative, alpha, omega
    proceed_all <- c()

    # total number of elements
    elements <- length(tpq) + length(taq) + length(proceed) + 2

    # indices
    tpq_idx <- 1:length(tpq)
    taq_idx <- (length(tpq) + 1):(length(tpq) + length(taq))
    proceed_idx <-  (1:length(proceed)) + (length(tpq) + length(taq))

    gibbs <- matrix(0, nrow = elements, ncol = size)

    gibbs_prd <- matrix(NA, nrow = elements, ncol = size)
    gibbs_use <- matrix(NA, nrow = elements, ncol = size)

    for (i in 1:length(sequences)) {
        sequences[[i]] <- c("alpha", sequences[[i]], "omega")   
    }

    j <- 0
    for (i in 1:length(tpq)) {
        if (!(tpq[[i]]$assoc %in% proceed)) {
            stop(paste0("Context of tpq ", tpq[[i]]$id, " : ", tpq[[i]]$assoc, " is not given in relative sequences"))
        }
        sequences <- c(sequences, list(c(tpq[[i]]$id, tpq[[i]]$assoc))  )
        proceed_all <- c(proceed_all, tpq[[i]]$id )
        gibbs[j+1,1] <- min(tpq[[i]]$samples)          # initialize tpq with earliest possible
        j <- j + 1
    }
    for (i in 1:length(taq)) {
        if (!(taq[[i]]$assoc %in% proceed)) {
            stop(paste0("Context of taq ", taq[[i]]$id, " : ", taq[[i]]$assoc, " is not given in relative sequences"))
        }
        sequences <- c(sequences, list(c(taq[[i]]$assoc, taq[[i]]$id)) )
        proceed_all <- c(proceed_all, taq[[i]]$id )
        gibbs[j+1,1] <- max(taq[[i]]$samples)          # initialize taq with latest possible
        j <- j + 1
    }

    proceed_all <- c(proceed_all, proceed, "alpha", "omega")
    gibbs[elements-1,1] <- alpha_
    gibbs[elements,1] <- omega_

    contexts_type_name <- proceed_all[proceed_all %in% contexts_type]
    contexts_type_idx <- which(proceed_all %in% contexts_type)
    contexts_type_absent_idx <- which(!(proceed_all %in% contexts_type))

    # convert to indices
    M <- list()
    for (i in 1:length(sequences)) {
        tmp <- numeric(length(sequences[[i]]))
        for (j in 1:length(sequences[[i]])) {
            tmp[j] <- which(proceed_all == sequences[[i]][j])
        }
        M[[i]] <- tmp
    }

    PhiMatrix <- quae_antea_matrix_cpp(elements, M)
    PsiMatrix <- quae_postea_matrix_cpp(elements, M)

    init_sample <- floor(sqrt(elements))

    # indices of non-trimmed relative events
    if (trim == TRUE) {
        idx_nontrim <- numeric(length(proceed))
        for (i in 1:length(proceed)) {
            idx <- proceed_idx[i]   
            check1 <- sum(PsiMatrix[idx, taq_idx])
            check2 <- sum(PhiMatrix[idx, tpq_idx])
            if (check1 > 0 & check2 > 0) {
                idx_nontrim[i] <- 1
            }
        } 
        idx_trim <- which(idx_nontrim == 0) + length(tpq) + length(taq)
        trim_label <- proceed[which(idx_nontrim == 0)]
        nontrim_label <- proceed_all[!(names(proceed_all) %in% trim_label)]
    }

    # initial sampler
    message("Assigning initial random values (this may take a moment)...")

    gibbs[,1] <- gibbs_ad_initial_cpp(gibbs[,1], tpq_idx, PsiMatrix, tpq, taq_idx, PhiMatrix, taq, proceed_idx, init_sample)

    if (is.matrix(PhiMatrix[contexts_type_idx,  ] )) {
        contexts_prior <- which(colSums(PhiMatrix[contexts_type_idx,  ]) > 0 )
    } else {
        contexts_prior <- which(PhiMatrix[contexts_type_idx,  ] > 0 )
    }
    contexts_prior_absent <- contexts_prior[contexts_prior %in% contexts_type_absent_idx & !(contexts_prior %in%  c((length(proceed_all)-1), length(proceed_all)))]

    # main sampler
    message("Beginning main Gibbs sampler. Will terminate either when MCSE criterion or maximum number of MC samples reached.")
    cat("Note: MCSE stopping criterion is only applied to sequences/constraints, not finds.\n")

    gibbs <- gibbs_ad_cpp(gibbs, tpq_idx, PsiMatrix, tpq , taq_idx, PhiMatrix, taq, proceed_idx)

    # finds production and use

    dep_i <- gibbs[ contexts_type_idx, ]

    if (!is.matrix(dep_i)) {
        dep_i <- t(as.matrix(dep_i, nrow = 1, byrow =TRUE))
    }
    Ym <- apply(dep_i,2,min)
    if (length(tpq_production) > 0) {
        Ym <- apply(rbind(Ym, gibbs[proceed_all %in% tpq_production, ]), 2, min  )
    }
    Ym_idx <- contexts_type_idx[apply(dep_i,2,which.min)]

    Xm <- max_antea(Ym_idx, gibbs, PhiMatrix, alpha_)

    if (sum(Xm < Ym) != ncol(gibbs)) {
        stop("Error in sequences/finds. Conflict in earliest production thresholds.")
    }

    Y <- matrix( Ym, nrow(dep_i), ncol(dep_i), byrow = TRUE)
    X <- matrix( Xm, nrow(dep_i), ncol(dep_i), byrow = TRUE)

    h_e <- matrix(stats::runif(nrow(X) * ncol(X) , as.vector(X) , as.vector(Y) ), nrow =  nrow(dep_i), ncol = ncol(dep_i))
    h_n <- matrix(stats::runif(nrow(X) * ncol(X) , as.vector(h_e) , as.vector(dep_i) ), nrow =  nrow(dep_i), ncol = ncol(dep_i))

    if (rule == "naive") {
        u <- matrix(stats::runif(nrow(X) * ncol(X) , as.vector(h_n) , as.vector(dep_i) ), nrow =  nrow(dep_i), ncol = ncol(dep_i))
    } else if (rule == "earliest") {
        u <- h_n
    }

    if (rule == "earliest") {
        gibbs_use[proceed_all %in% contexts_type, ] <- u
        gibbs_prd[proceed_all %in% contexts_type, ] <- h_e
    } else if (rule == "naive") {
        gibbs_use[proceed_all %in% contexts_type, ] <- u
        gibbs_prd[proceed_all %in% contexts_type, ] <- h_n
    }

    # consistent batch means
    mcse_check <- FALSE
    while (mcse_check == FALSE) {
        n_upto <- ncol(gibbs)
        n_batch <- floor(sqrt(n_upto)) # length of samples in batch
        K <- floor(n_upto / n_batch) # number of batches
        m_batch <- matrix(NA, nrow = nrow(gibbs), ncol = (K-1))
        remainder <- n_upto - n_batch * K + 1

        idx1 <- remainder

        for (k in 1:(K-1)) {
            idxs <- idx1:(idx1 + n_batch)       

            # in cases where n_upto = n_batch * K
            idxs <- idxs[idxs <= ncol(gibbs)]

            m_batch[, k] <- rowMeans(gibbs[ , idxs]) 
            idx1 <- idxs[length(idxs)] + 1
        }
        
        mcse0 <- sqrt( rowSums( (m_batch - rowMeans(m_batch))^2 ) * (n_batch / (K-1) ) ) / sqrt((K - 1) * n_batch)
        
        # remove trimmed events from estimation of mean MCSE
        mcse <- mcse0
        if (trim == TRUE) {
            if (length(idx_trim) > 0) {
                mcse <- mcse[-idx_trim]
            }
        }
        # do not include fixed single-point events as part of estimating MCSE
        mcse <- mcse[mcse > 0] 

        cat("\r", paste0("Samples: ", ncol(gibbs), "     Mean MCSE: ",  round(mean(mcse),3)))
        if (mean(mcse) < mcse_crit) {
            mcse_check <- TRUE
            message("\nMCSE criterion passed. Finishing.")
        } else {
            if (ncol(gibbs) >= max_samples) {
                message("\nMC samples exceeded maximum stipulated without passing MCSE criterion. Finishing.")
                mcse_check <- TRUE
            } else {

                gibbs_next <- matrix(0, nrow = nrow(gibbs), ncol = (size + 1) )
                gibbs_next[,1] <- gibbs[,ncol(gibbs)]
                gibbs_next <- gibbs_ad_cpp(gibbs_next, tpq_idx, PsiMatrix, tpq , taq_idx, PhiMatrix, taq, proceed_idx)

                # finds production and use
                        
                gibbs_next_prd <- matrix(NA, nrow = elements, ncol = (size + 1))
                gibbs_next_use <- matrix(NA, nrow = elements, ncol = (size + 1))

                dep_i <- gibbs_next[ contexts_type_idx, ]

                if (!is.matrix(dep_i)) {
                    dep_i <- t(as.matrix(dep_i, nrow = 1, byrow =TRUE))
                }
                Ym <- apply(dep_i,2,min)
                if (length(tpq_production) > 0) {
                    Ym <- apply(rbind(Ym, gibbs_next[proceed_all %in% tpq_production, ]), 2, min  )
                }

                Ym_idx <- contexts_type_idx[apply(dep_i,2,which.min)]

                Xm <- max_antea(Ym_idx, gibbs_next, PhiMatrix, alpha_)

                if (sum(Xm < Ym) != ncol(gibbs_next)) {
                    stop("Error in sequences/finds. Conflict in earliest production thresholds.")
                }

                Y <- matrix( Ym, nrow(dep_i), ncol(dep_i), byrow = TRUE)
                X <- matrix( Xm, nrow(dep_i), ncol(dep_i), byrow = TRUE)

                h_e <- matrix(stats::runif(nrow(X) * ncol(X) , as.vector(X) , as.vector(Y) ), nrow =  nrow(dep_i), ncol = ncol(dep_i))
                h_n <- matrix(stats::runif(nrow(X) * ncol(X) , as.vector(h_e) , as.vector(dep_i) ), nrow =  nrow(dep_i), ncol = ncol(dep_i))

                if (rule == "naive") {
                    u <- matrix(stats::runif(nrow(X) * ncol(X) , as.vector(h_n) , as.vector(dep_i) ), nrow =  nrow(dep_i), ncol = ncol(dep_i))
                } else if (rule == "earliest") {
                    u <- h_n
                }

                if (rule == "earliest") {
                    gibbs_next_use[proceed_all %in% contexts_type, ] <- u
                    gibbs_next_prd[proceed_all %in% contexts_type, ] <- h_e
                } else if (rule == "naive") {
                    gibbs_next_use[proceed_all %in% contexts_type, ] <- u
                    gibbs_next_prd[proceed_all %in% contexts_type, ] <- h_n
                }

                gibbs <- cbind(gibbs, gibbs_next[, 2:ncol(gibbs_next)])
                gibbs_use <- cbind(gibbs_use, gibbs_next_use[, 2:ncol(gibbs_next_use)])
                gibbs_prd <- cbind(gibbs_prd, gibbs_next_prd[, 2:ncol(gibbs_next_prd)])

            }
        }
    }

    find_prd_use_dep <- list()

    type_deposition <- c( gibbs[proceed_all %in% contexts_type,] )
    type_use <- c( gibbs_use[proceed_all %in% contexts_type,] )
    type_production <- c( gibbs_prd[proceed_all %in% contexts_type,] )

    type_all <- rbind(type_production, type_use, type_deposition)
    n_upto <- ncol(type_all)
    n_batch <- floor(sqrt(n_upto)) # length of samples in batch
    K <- floor(n_upto / n_batch) # number of batches
    m_batch <- matrix(NA, nrow = nrow(type_all), ncol = (K-1))
    remainder <- n_upto - n_batch * K + 1

    idx1 <- remainder

    for (k in 1:(K-1)) {
        idxs <- idx1:(idx1 + n_batch)       

        # in cases where n_upto = n_batch * K
        idxs <- idxs[idxs <= ncol(type_all)]

        m_batch[, k] <- rowMeans(type_all[ , idxs]) 
        idx1 <- idxs[length(idxs)] + 1
    }
    
    mcmean <- rowMeans(m_batch) 
    mcse <- sqrt( rowSums( (m_batch - rowMeans(m_batch))^2 ) * (n_batch / (K-1) ) ) / sqrt((K - 1) * n_batch)

    type_stat <-data.frame(MCmean = mcmean, MCSE = mcse)
    rownames(type_stat) <- c("production", "use", "deposition")

    type_ <- list(production = type_production, use = type_use, deposition = type_deposition)

    if (is.null(type_name)) {
        if (!is.null(type)) {
            if (length(type) == 1) {
                type_name <- type
            } else {
                type_name <- "Type"
            }
        } else {
            type_name <- "Type"
        }
    }
    res <- list(name = type_name, type = type_, stat = type_stat)
    class(res) <- c("type_marginals", "list")
    return(res)

}



#' Mean Squared Displacement of Events
#' 
#' Computes the mean squared displacement (MSD) of all events contained in the relative sequences and absolute constraints used in the execution of \code{\link[eratosthenes]{gibbs_ad}}. MSD is not intended for finds, as production, use, and depositional dates, as these are themselves contingent upon the relative/absolute events.
#' 
#' The MSD entails the following jackknife/leave-one-out style routine:
#' 
#' * Each event is omitted from all relative and absolute sequences, and the function \code{\link[eratosthenes]{gibbs_ad}} is re-run to compute a "jackknifed" Monte Carlo mean for that event.
#'   * The squared difference of this jackknifed Monte Carlo mean and the original is then computed as its squared "displacement" in time.
#'   * The mean of the squared displacements of all events is then computed and attributed to the omitted event.
#' 
#' If an event has a low MSD, it bears a low impact on the rest of the events within the full joint conditional density. If it is has a high MSD, other events depend heavily upon its inclusion in the full joint density.
#' 
#' Trimming is not implemented in the computation of MSD, and so attention should be paid to the selection of \code{alpha_} and \code{omega_}, which should be reported. This is owing to the way in which, if an absolute constraint (\code{tpq} or \code{taq}) is omitted that happens to be an earliest or latest bounding event, there still needs to be earliest and latest thresholds in place. 
#' 
#' This function is fairly computationally intensive and thus a lower value of `max_samples` and a higher value of `mcse_crit` may be warranted.
#'  
#' @param marginalized An object of class \code{marginals}, the output of \code{\link[eratosthenes]{gibbs_ad}}.
#' @param sequences A \code{\link[eratosthenes]{sequences}} object of relative sequences of elements (e.g., contexts) used to compute \code{marginalized}.
#' @param max_samples Maximum number of samples to run. Default is \code{10^5}.
#' @param size The number of samples to take on each iteration of the main Gibbs sampler. Default is \code{10^3}. 
#' @param mcse_crit Criterion for the Monte Carlo standard error to stop the Gibbs sampler. A higher MCSE is recommended for situations with a higher number of events in order to reduce computational time.
#' @param tpq A \code{list} containing \emph{termini post quos} used to compute \code{marginalized}. See \code{\link[eratosthenes]{gibbs_ad}} for details.
#' @param taq A \code{list} containing \emph{termini ante quos} used to compute \code{marginalized}. See \code{\link[eratosthenes]{gibbs_ad}} for details.
#' @param alpha_ An initial \emph{t.p.q.} to limit any elements which may occur before the first provided \emph{t.p.q.} Default is \code{-5000}.
#' @param omega_ A final \emph{t.a.q.} to limit any elements which may occur after the after the last provided \emph{t.a.q.} Default is \code{1950}.
#' 
#' @examples 
#' x <- events("A", "B", "C", "D", "E", "F", "G", "H", "I", "J")
#' y <- events("B", "D", "G", "H", "K")
#' z <- events("F", "K", "L", "M")
#' contexts <- sequences(x, y, z)
#' 
#' f1 <- finds(id = "find01", assoc = "D", type = c("type1", "form1"))
#' f2 <- finds(id = "find02", assoc = "E", type = c("type1", "form2"))
#' f3 <- finds(id = "find03", assoc = "G", type = c("type1", "form1"))
#' f4 <- finds(id = "find04", assoc = "H", type = c("type2", "form1"))
#' f5 <- finds(id = "find05", assoc = "I", type = "type2")
#' f6 <- finds(id = "find06", assoc = "H", type = NULL)
#' 
#' artifacts <- assemblage(f1, f2, f3, f4, f5, f6)
#'  
#' # external constraints
#' coin1 <- absolute(id = "coin1", assoc = "B", type = NULL, samples = runif(100,-320,-300))
#' coin2 <- absolute(id = "coin2", assoc = "G", type = NULL, samples = seq(37, 41, length = 100))
#' destr <- absolute(id = "destr", assoc = "J", type = NULL, samples = 79)
#' 
#' tpq_info <- constraints(coin1, coin2)
#' taq_info <- constraints(destr)
#' 
#' result <- gibbs_ad(contexts, tpq = tpq_info, taq = taq_info)
#' 
#' result_msd <- msd(result, contexts, max_samples = 5000,
#'                   mcse_crit = 2, tpq = tpq_info, taq = taq_info)
#'
#' @returns Output is a list containing a data frame \code{MSD_stats} giving the mean MC date, the MCSE, the MSD, the variance of the squared displacements (not the standard error), and sample size, as well as a vector \code{bounds} of the values of \code{alpha_} and \code{omega_}.
#' 
#' @export
msd <- function(marginalized, sequences, max_samples = 10^5, size = 10^3, mcse_crit = 0.5, tpq = NULL, taq = NULL, alpha_ = -5000, omega_ = 1950) {
    UseMethod("msd")
}
#' 
#' @rdname msd
#' @export
msd.marginals <- function(marginalized, sequences,  max_samples = 10^5, size = 10^3, mcse_crit = 0.5, tpq = NULL, taq = NULL, alpha_ = -5000, omega_ = 1950) {
    if (size > max_samples) {
        stop("Error: size must be less than max_samples.")
    }
    if (!inherits(sequences, "sequences")) {
        stop("input sequences must be sequences object. See sequences().")
    }
    if (!(is.null(tpq) | inherits(tpq, "constraints"))) {
        stop("input tpq must be constraints object. See constraints().")
    }
    if (!(is.null(taq) | inherits(taq, "constraints"))) {
        stop("input taq must be constraints object. See constraints().")
    }

    depmu <- sapply(marginalized$deposition, mean)
    depmcse <- marginalized$mcse[names(marginalized$deposition)]
    depdat <- data.frame(Mean = depmu, MCSE = depmcse)
    rownames(depdat) <- names(marginalized$deposition)

    extmu <- sapply(marginalized$externals, mean)
    extmcse <- marginalized$mcse[names(marginalized$externals)]
    extdat <- data.frame(Mean = extmu, MCSE = extmcse)
    rownames(extdat) <- names(marginalized$externals)

    orig_dat <- rbind(depdat, extdat)
    
    proceed <- synth_rank(sequences)

    # proceed_all from tpq, taq, relative, alpha, omega
    # proceed_all <- c()

    # total number of elements
    elements <- length(proceed) + 2
    if (!is.null(tpq)) {
        elements <- elements + length(tpq)
    }
    if (!is.null(taq)) {
        elements <- elements + length(taq)
    }

    # indices
    if (is.null(tpq) & is.null(taq)) {
        proceed_idx <- 1:length(proceed) 
    } else if (is.null(tpq) & !is.null(taq)) {
        taq_idx <- 1:length(taq)
        proceed_idx <- (1:length(proceed)) + length(taq)
    } else if (!is.null(tpq) & is.null(taq)) {
        tpq_idx <- 1:length(tpq)
        proceed_idx <- (1:length(proceed)) + length(tpq)
    } else {
        tpq_idx <- 1:length(tpq)
        taq_idx <- (length(tpq) + 1):(length(tpq) + length(taq))
        proceed_idx <-  (1:length(proceed)) + (length(tpq) + length(taq))
    }

    proceed_all <- rownames(orig_dat)

    orig_dat$MSD <- NA
    orig_dat$MSD_var <- NA
    orig_dat$MSD_n <- NA

    message("Beginning jackknife/LOO-style routine to compute MSD. This may take a while, depending on the number of events / mcse_crit...\n")

    for (j in 1:length(proceed_all)) {
        cat("Depositional Event / Absolute Constraint: ", proceed_all[j], "\n")
        MCmean <- orig_dat[proceed_all[j],]$Mean

        sequencesMSD <- list()
        tpqMSD <- list()
        taqMSD <- list()
        for (i in 1:length(sequences)) {
            sq_ <- sequences[[i]]
            sq_ <-sq_[!(sq_ %in% proceed_all[j])]
            sequencesMSD[[i]] <- events(sq_)
        }

        if (length(tpq) > 0) {
            idx <- 1
            for (i in 1:length(tpq)) {
                if (!(tpq[[i]]$id %in% proceed_all[j] | tpq[[i]]$assoc %in% proceed_all[j])) {
                    tpqMSD[[idx]] <- tpq[[i]]
                    idx <- idx + 1
                }
            }
            if (length(tpqMSD) > 0) {
                tpqMSD <- constraints(tpqMSD)
            }
        }
        if (length(taq) > 0) {
            idx <- 1
            for (i in 1:length(taq)) {
                if (!(taq[[i]]$id %in% proceed_all[j] | taq[[i]]$assoc %in% proceed_all[j])) {
                    taqMSD[[idx]] <- taq[[i]]
                    idx <- idx + 1
                }
            }
            if (length(taqMSD) > 0) {
                taqMSD <- constraints(taqMSD)
            }
        }
        sequencesMSD <- sequences(sequencesMSD)

        if (length(tpqMSD) == 0 & length(taqMSD) != 0) {
            gibbsLOO <- gibbs_ad(sequencesMSD, max_samples, size, mcse_crit, tpq = NULL, taq = taqMSD, alpha_, omega_, trim = FALSE)
        } else if (length(tpqMSD) != 0 & length(taqMSD) == 0) {
            gibbsLOO <- gibbs_ad(sequencesMSD, max_samples, size, mcse_crit, tpq = tpqMSD, taq = NULL, alpha_, omega_, trim = FALSE)
        } else if (length(tpqMSD) == 0 & length(taqMSD) == 0) {
            gibbsLOO <- gibbs_ad(sequencesMSD, max_samples, size, mcse_crit, tpq = NULL, taq = NULL, alpha_, omega_, trim = FALSE)
        } else {
            gibbsLOO <- gibbs_ad(sequencesMSD, max_samples, size, mcse_crit, tpq = tpqMSD, taq = taqMSD, alpha_, omega_, trim = FALSE)
        }
        
        depmu <- sapply(gibbsLOO$deposition, mean)
        depmcse <- gibbsLOO$mcse[names(gibbsLOO$deposition)]
        depdat <- data.frame(Mean = depmu, MCSE = depmcse)
        rownames(depdat) <- names(gibbsLOO$deposition)

        extmu <- sapply(gibbsLOO$externals, mean)
        extmcse <- gibbsLOO$mcse[names(gibbsLOO$externals)]
        extdat <- data.frame(Mean = extmu, MCSE = extmcse)
        rownames(extdat) <- names(gibbsLOO$externals)

        LOO_dat <- rbind(depdat, extdat)

        idx_name <- rownames(LOO_dat)[rownames(LOO_dat) %in% rownames(orig_dat) ]
        MSD_ <- sqrt((LOO_dat[idx_name, ]$Mean - orig_dat[idx_name, ]$Mean)^2)
        MSD_ <- MSD_[MSD_ > 0] # avoid including fixed events
        orig_dat[which(rownames(orig_dat)==proceed_all[j]) , 3] <- mean(MSD_)
        orig_dat[which(rownames(orig_dat)==proceed_all[j]) , 4] <- stats::var(MSD_)
        orig_dat[which(rownames(orig_dat)==proceed_all[j]) , 5] <- length(MSD_)

        cat("\n")
    }

    message("Estimation of MSD complete.")
    bounds_ <- c(alpha_, omega_)
    names(bounds_) <- c("alpha", "omega")

    res <- list(MSD_stats = orig_dat, bounds = bounds_)
    class(res) <- c("msd_data", "list")
    return(res)
}



#' @export
print.msd_data <- function(x, ...) {  
    cat("\n For fixed bounds: (", x$bounds[1], ",",x$bounds[2], ")\n",
    "MSD estimates, variance, and sample size: \n")
    print.data.frame(x$MSD_stats)
    cat("\n")
}



#' Squared Displacement for a Target Event
#' 
#' Computes the squared displacement for a target event within the joint conditional density, estimating how much the omission of every other event will change the date of the target. See also \code{\link[eratosthenes]{msd}}. If the target event is a find or type, the displacement of the use date is used, since use is contingent upon both production and deposition.
#' 
#' Displacement is computed via the following jackknife/leave-one-out-style routine:
#' 
#' * Each event, excluding the target event itself, is omitted from all relative and absolute sequences, and the function \code{\link[eratosthenes]{gibbs_ad}} is re-run to compute a "jackknifed" Monte Carlo mean for the target event.
#'   * The squared difference of this jackknifed Monte Carlo mean and the original is then computed as its squared "displacement" in time.
#' 
#' If an event has a low squared displacement, it has a low impact on the dating of the target event. If it is has a high squared displacement, the target event's date depends heavily upon its inclusion in the full joint density.
#' 
#' Trimming is not implemented in the estimation of squared displacement, and so attention should be paid to the selection of \code{alpha_} and \code{omega_}, and reported. This is owing to the way in which, if an absolute constraint (\code{tpq} or \code{taq}) is omitted that happens to be an earliest or latest bounding event, there still needs to be earliest and latest thresholds in place. 
#' 
#' This function is fairly computationally intensive, and so a lower value of `max_samples` or higher value of `mcse_crit` may be warranted.
#' 
#' @param marginalized The results of \code{\link[eratosthenes]{gibbs_ad}} or \code{\link[eratosthenes]{gibbs_ad_type}}.
#' @param target The target event (any event for which to estimate squared displacement. If using the results of \code{\link[eratosthenes]{gibbs_ad_type}}, that type is by default the target (otherwise, for sequences/\emph{t.p.q.}/\emph{t.a.q.} one should use the output of \code{\link[eratosthenes]{gibbs_ad}}).
#' @param sequences A \code{list} of relative sequences of elements (e.g., contexts) used to compute \code{marginalized}.
#' @param finds Optional. A \code{list} of finds related to (contained in) the elements of \code{sequences}.
#' @param max_samples Maximum number of samples to run. Default is \code{10^5}.
#' @param size The number of samples to take on each iteration of the main Gibbs sampler. Default is \code{10^3}. 
#' @param mcse_crit Criterion for the Monte Carlo standard error to stop the Gibbs sampler. A higher MCSE is recommended for situations with a higher number of events in order to reduce computational time.
#' @param tpq A \code{list} containing \emph{termini post quos} used to compute \code{marginalized}. See \code{\link[eratosthenes]{gibbs_ad}} for details.
#' @param taq A \code{list} containing \emph{termini ante quos} used to compute \code{marginalized}. See \code{\link[eratosthenes]{gibbs_ad}} for details.
#' @param alpha_ An initial \emph{t.p.q.} to limit any elements which may occur before the first provided \emph{t.p.q.} Default is \code{-5000}.
#' @param omega_ A final \emph{t.a.q.} to limit any elements which may occur after the after the last provided \emph{t.a.q.} Default is \code{1950}.
#' @param rule The rule for computing an estimated date of production, if using an artifact type as a target date. See \code{\link[eratosthenes]{gibbs_ad_type}} for details.
#' 
#' @examples 
#' x <- events("A", "B", "C", "D", "E", "F", "G", "H", "I", "J")
#' y <- events("B", "D", "G", "H", "K")
#' z <- events("F", "K", "L", "M")
#' contexts <- sequences(x, y, z)
#' 
#' f1 <- finds(id = "find01", assoc = "D", type = c("type1", "form1"))
#' f2 <- finds(id = "find02", assoc = "E", type = c("type1", "form2"))
#' f3 <- finds(id = "find03", assoc = "G", type = c("type1", "form1"))
#' f4 <- finds(id = "find04", assoc = "H", type = c("type2", "form1"))
#' f5 <- finds(id = "find05", assoc = "I", type = "type2")
#' f6 <- finds(id = "find06", assoc = "H", type = NULL)
#' 
#' artifacts <- assemblage(f1, f2, f3, f4, f5, f6)
#'  
#' # external constraints
#' coin1 <- absolute(id = "coin1", assoc = "B", type = NULL, samples = runif(100,-320,-300))
#' coin2 <- absolute(id = "coin2", assoc = "G", type = NULL, samples = seq(37, 41, length = 100))
#' destr <- absolute(id = "destr", assoc = "J", type = NULL, samples = 79)
#' 
#' tpq_info <- constraints(coin1, coin2)
#' taq_info <- constraints(destr)
#' 
#' result <- gibbs_ad(contexts, tpq = tpq_info, taq = taq_info)
#' 
#' # max_samples lowered and msce_crit raised for examples
#' 
#' # squared displacement for depositional context "E"
#' E_sqdisp <- sq_disp(result, target = "E", sequences = contexts, 
#'                     max_samples = 3000, mcse_crit = 2, tpq = tpq_info, taq = taq_info)
#'
#' result_type1 <- gibbs_ad_type(contexts, finds = artifacts, type = "type1",
#'                               tpq = tpq_info, taq = taq_info)
#' 
#' # squared displacement for production of artifact type "type1"
#' type1_sqdisp <- sq_disp(result_type1, sequences = contexts, finds = artifacts,
#'                         max_samples = 3000, mcse_crit = 2, tpq = tpq_info, taq = taq_info)
#'
#' @returns Output is a list containing a data frame \code{sq_disp} giving the diplacement with respect to all other events and a vector \code{bounds} of the values of \code{alpha_} and \code{omega_}.
#' 
#' @export
sq_disp <- function(marginalized, target = NULL, sequences, finds = NULL, max_samples = 10^5, size = 10^3, mcse_crit = 0.5, tpq = NULL, taq = NULL, alpha_ = -5000, omega_ = 1950, rule = "naive") {
    UseMethod("sq_disp")
}
#' 
#' @rdname sq_disp
#' @export
sq_disp.marginals <- function(marginalized, target = NULL, sequences, finds = NULL, max_samples = 10^5, size = 10^3, mcse_crit = 0.5, tpq = NULL, taq = NULL, alpha_ = -5000, omega_ = 1950, rule = NULL) {
    if (!inherits(sequences, "sequences")) {
        stop("sequences must be sequences object. See sequences().")
    }
    if (!(is.null(finds) | inherits(finds, "assemblage"))) {
        stop("finds must be NULL or assemblage object. See assemblage().")
    }
    if (!(is.null(tpq) | inherits(tpq, "constraints"))) {
        stop("input tpq must be constraints object. See constraints().")
    }
    if (!(is.null(taq) | inherits(taq, "constraints"))) {
        stop("input taq must be constraints object. See constraints().")
    }

    if (is.null(target)) {
        stop("Target event must be specified.")
    }
    if (size > max_samples) {
        stop("Error: size must be less than max_samples.")
    }
    depmu <- sapply(marginalized$deposition, mean)
    depmcse <- marginalized$mcse[names(marginalized$deposition)]
    depdat <- data.frame(Mean = depmu, MCSE = depmcse)
    rownames(depdat) <- names(marginalized$deposition)

    extmu <- sapply(marginalized$externals, mean)
    extmcse <- marginalized$mcse[names(marginalized$externals)]
    extdat <- data.frame(Mean = extmu, MCSE = extmcse)
    rownames(extdat) <- names(marginalized$externals)
    
    orig_dat <- rbind(depdat, extdat)

    if (!(is.null(finds))) {
        prdmu <- sapply(marginalized$production, mean)
        prdmcse <- marginalized$mcse[names(marginalized$production)]
        prddat <- data.frame(Mean = prdmu, MCSE = prdmcse)
        rownames(prddat) <- names(marginalized$production)

        orig_dat2 <- rbind(depdat, extdat, prddat)
    }
   
    proceed <- synth_rank(sequences)

    # proceed_all from tpq, taq, relative, alpha, omega
    proceed_all <- c()

    # total number of elements
    elements <- length(tpq) + length(taq) + length(proceed) + 2

    # indices
    tpq_idx <- 1:length(tpq)
    taq_idx <- (length(tpq) + 1):(length(tpq) + length(taq))
    proceed_idx <-  (1:length(proceed)) + (length(tpq) + length(taq))

    proceed_all <- rownames(orig_dat)
    proceed_all <- proceed_all[!(proceed_all %in% target)]

    orig_dat$sq_disp <- NA

    message("Beginning jackknife/LOO-style routine to compute squared displacement. This may take a while, depending on the number of events / mcse_crit...\n")

    for (j in 1:length(proceed_all)) {
        cat("Depositional Event / Absolute Constraint: ", proceed_all[j], "\n")
        MCmean <- orig_dat[proceed_all[j],]$Mean

        sequencesMSD <- list()
        tpqMSD <- list()
        taqMSD <- list()
        for (i in 1:length(sequences)) {
            sq_ <- sequences[[i]]
            sq_ <-sq_[!(sq_ %in% proceed_all[j])]
            sequencesMSD[[i]] <- events(sq_)
        }
        if (length(tpq) > 0) {
            idx <- 1
            for (i in 1:length(tpq)) {
                if (!(tpq[[i]]$id %in% proceed_all[j] | tpq[[i]]$assoc %in% proceed_all[j])) {
                    tpqMSD[[idx]] <- tpq[[i]]
                    idx <- idx + 1
                }
            }
            if (length(tpqMSD) > 0) {
                tpqMSD <- constraints(tpqMSD)
            }
        }
        if (length(taq) > 0) {
            idx <- 1
            for (i in 1:length(taq)) {
                if (!(taq[[i]]$id %in% proceed_all[j] | taq[[i]]$assoc %in% proceed_all[j])) {
                    taqMSD[[idx]] <- taq[[i]]
                    idx <- idx + 1
                }
            }
            if (length(taqMSD) > 0) {
                taqMSD <- constraints(taqMSD)
            }
        }
        sequencesMSD <- sequences(sequencesMSD)

        if (length(tpqMSD) == 0 & length(taqMSD) != 0) {
            gibbsLOO <- gibbs_ad(sequencesMSD, max_samples, size, mcse_crit, tpq = NULL, taq = taqMSD, alpha_, omega_, trim = FALSE)
        } else if (length(tpqMSD) != 0 & length(taqMSD) == 0) {
            gibbsLOO <- gibbs_ad(sequencesMSD, max_samples, size, mcse_crit, tpq = tpqMSD, taq = NULL, alpha_, omega_, trim = FALSE)
        } else if (length(tpqMSD) == 0 & length(taqMSD) == 0) {
            gibbsLOO <- gibbs_ad(sequencesMSD, max_samples, size, mcse_crit, tpq = NULL, taq = NULL, alpha_, omega_, trim = FALSE)
        } else {
            gibbsLOO <- gibbs_ad(sequencesMSD, max_samples, size, mcse_crit, tpq = tpqMSD, taq = taqMSD, alpha_, omega_, trim = FALSE)
        }

        depmu <- sapply(gibbsLOO$deposition, mean)
        depmcse <- gibbsLOO$mcse[names(gibbsLOO$deposition)]
        depdat <- data.frame(Mean = depmu, MCSE = depmcse)
        rownames(depdat) <- names(gibbsLOO$deposition)

        extmu <- sapply(gibbsLOO$externals, mean)
        extmcse <- gibbsLOO$mcse[names(gibbsLOO$externals)]
        extdat <- data.frame(Mean = extmu, MCSE = extmcse)
        rownames(extdat) <- names(gibbsLOO$externals)

        if (is.null(finds)) {
            LOO_dat <- rbind(depdat, extdat)
            sqdisp_ <- sqrt((LOO_dat[target, ]$Mean - orig_dat[target,]$Mean)^2)
        } else {
            prdmu <- sapply(gibbsLOO$production, mean)
            prdmcse <- gibbsLOO$mcse[names(gibbsLOO$production)]
            prddat <- data.frame(Mean = prdmu, MCSE = prdmcse)
            rownames(prddat) <- names(gibbsLOO$production)
            LOO_dat <- rbind(depdat, extdat, prddat)
            sqdisp_ <- sqrt((LOO_dat[target, ]$Mean - orig_dat2[target,]$Mean)^2)
        }

        orig_dat[which(rownames(orig_dat)==proceed_all[j]) , 3] <- sqdisp_

        cat("\n")
    }

    message("Estimation of squared displacement complete.")
    bounds_ <- c(alpha_, omega_)
    names(bounds_) <- c("alpha", "omega")

    orig_dat$Mean <- NULL
    orig_dat$MCSE <- NULL

    res <- list(sq_disp = orig_dat, bounds = bounds_, target = target)
    class(res) <- c("sq_displ_data", "list")
    return(res)
}
#' 
#' @rdname sq_disp
#' @export
sq_disp.type_marginals <- function(marginalized, target = NULL, sequences, finds = NULL, max_samples = 10^5, size = 10^3, mcse_crit = 0.5, tpq = NULL, taq = NULL, alpha_ = -5000, omega_ = 1950, rule = "naive") {
    if (!inherits(sequences, "sequences")) {
        stop("sequences must be sequences object. See sequences().")
    }
    if (!(is.null(finds) | inherits(finds, "assemblage"))) {
        stop("finds must be NULL or assemblage object. See assemblage().")
    }
    if (!(is.null(tpq) | inherits(tpq, "constraints"))) {
        stop("input tpq must be constraints object. See constraints().")
    }
    if (!(is.null(taq) | inherits(taq, "constraints"))) {
        stop("input taq must be constraints object. See constraints().")
    }
    if (size > max_samples) {
        stop("Error: size must be less than max_samples.")
    }
    target <- marginalized$name

    if (is.null(finds)) {
        stop("Error: finds is NULL.")
    }

    orig_mu <- marginalized$stat["use",]$MCmean
    orign_se <- marginalized$stat["use",]$MCSE

    proceed <- synth_rank(sequences)

    if (!is.list(tpq)) {
        tpq <- list(list(id = "tpq_default", assoc = proceed[1], type = NULL, samples = alpha_))
    }
    if (!is.list(taq)) {
        taq <- list(list(id = "taq_default", assoc = proceed[length(proceed)], type = NULL, samples = omega_))
    }

    # proceed_all from tpq, taq, relative
    proceed_all <- c()

    for (i in 1:length(tpq)) {
        if (!(tpq[[i]]$assoc %in% proceed)) {
            stop(paste0("Context of tpq ", tpq[[i]]$id, " : ", tpq[[i]]$assoc, " is not given in relative sequences"))
        }
        proceed_all <- c(proceed_all, tpq[[i]]$id )
    }
    for (i in 1:length(taq)) {
        if (!(taq[[i]]$assoc %in% proceed)) {
            stop(paste0("Context of taq ", taq[[i]]$id, " : ", taq[[i]]$assoc, " is not given in relative sequences"))
        }
        proceed_all <- c(proceed_all, taq[[i]]$id )
    }
    proceed_all <- c(proceed_all, proceed)

    res <- data.frame(sq_disp = rep(NA, length(proceed_all)), disp_MCmean = rep(NA, length(proceed_all)), disp_MCSE = rep(NA, length(proceed_all)) )
    rownames(res) <- proceed_all

    message("Beginning jackknife/LOO-style routine to compute squared displacement. This may take a while, depending on the number of events / mcse_crit...\n")

    for (j in 1:length(proceed_all)) {
        cat("Depositional Event / Absolute Constraint: ", proceed_all[j], "\n")

        sequencesSD <- list()
        tpqSD <- list()
        taqSD <- list()
        for (i in 1:length(sequences)) {
            sq_ <- sequences[[i]]
            sq_ <-sq_[!(sq_ %in% proceed_all[j])]
            sequencesSD[[i]] <- events(sq_)
        }
        idx <- 1
        for (i in 1:length(tpq)) {
            if (!(tpq[[i]]$id %in% proceed_all[j] | tpq[[i]]$assoc %in% proceed_all[j])) {
                tpqSD[[idx]] <- tpq[[i]]
                idx <- idx + 1
            }
        }
        idx <- 1
        for (i in 1:length(taq)) {
            if (!(taq[[i]]$id %in% proceed_all[j] | taq[[i]]$assoc %in% proceed_all[j])) {
                taqSD[[idx]] <- taq[[i]]
                idx <- idx + 1
            }
        }

        findsSD <- list()
        idx <- 1
        for (i in 1:length(finds)) {
            if (!(finds[[i]]$assoc %in% proceed_all[j])) {
                findsSD[[idx]] <- finds[[i]]
                idx <- idx + 1
            }
        }
        sequencesSD <- sequences(sequencesSD)
        if (length(tpqSD)) {
            tpqSD <- constraints(tpqSD)
        }
        if (length(taqSD)) {
            taqSD <- constraints(taqSD)
        }
        if (length(findsSD) > 0) {
            findsSD <- assemblage(findsSD)

            if (length(tpqSD) == 0 & length(taqSD) != 0) {
                gibbsLOO <- gibbs_ad_type(sequences = sequencesSD, finds = findsSD, id = NULL, type = target, type_name = target, max_samples = max_samples, size = size, mcse_crit = mcse_crit, tpq = NULL, taq = taqSD, alpha_ = alpha_, omega_ = omega_, trim = FALSE, rule = rule)
            } else if (length(tpqSD) != 0 & length(taqSD) == 0) {
                gibbsLOO <- gibbs_ad_type(sequences = sequencesSD, finds = findsSD, id = NULL, type = target, type_name = target, max_samples = max_samples, size = size, mcse_crit = mcse_crit, tpq = tpqSD, taq = NULL, alpha_ = alpha_, omega_= omega_, trim = FALSE, rule = rule)
            } else if (length(tpqSD) == 0 & length(taqSD) == 0) {
                gibbsLOO <- gibbs_ad_type(sequences = sequencesSD, finds = findsSD, id = NULL, type = target, type_name = target, max_samples = max_samples, size = size, mcse_crit = mcse_crit, tpq = NULL, taq = NULL, alpha_ = alpha_, omega_= omega_, trim = FALSE, rule = rule)
            } else {
                gibbsLOO <- gibbs_ad_type(sequences = sequencesSD, finds = findsSD, id = NULL, type = target, type_name = target, max_samples = max_samples, size = size, mcse_crit = mcse_crit, tpq = tpqSD, taq = taqSD, alpha_ = alpha_, omega_= omega_, trim = FALSE, rule = rule)
            }

            disp_mu <- gibbsLOO$stat["use",]$MCmean
            disp_mcse <- gibbsLOO$stat["use",]$MCSE
            
            res[which(rownames(res)==proceed_all[j]) , 1] <- (disp_mu - orig_mu)^2
            res[which(rownames(res)==proceed_all[j]) , 2] <- disp_mu
            res[which(rownames(res)==proceed_all[j]) , 3] <- disp_mcse

            cat("\n") 
        } else {
            cat("Event", proceed_all[j], "skipped: type completely removed from relationships (not possible to estimate).\n")
        }
    }

    message("Estimation of squared displacement complete.")
    bounds_ <- c(alpha_, omega_)
    names(bounds_) <- c("alpha", "omega")

    res <- list(sq_disp = res, bounds = bounds_, target = target)
    class(res) <- c("sq_displ_data", "list")
    return(res)
}



#' @export
print.sq_displ_data <- function(x, ...) {  
    cat("\n For fixed bounds: (", x$bounds[1], ",",x$bounds[2], ")\n",
    "Squared displacement for target event", x$target, "caused\n when omitting the following:  \n\n")
    print.data.frame(x$sq_disp)
    cat("\n")
}





#' Create an Events Object
#' 
#' Analogous to the \code{\link[base]{c}} function, to create a sequence of unique events as a vector. Elements may not contain names of \code{"alpha"} or \code{"omega"}, which are restricted for \code{\link[eratosthenes]{quae_antea}} and \code{\link[eratosthenes]{quae_postea}}.
#' 
#' @param ... Comma separated character elements, in order from left (earliest) to right (latest).
#' 
#' @examples
#' # "A" before "B", "B" before "C"
#' x <- events("A", "B", "C")
#' 
#' @returns An events object.
#' 
#' @export
events <- function(...) {
    UseMethod("events")
}
#' 
#' @rdname events
#' @export
events.character <- function(...) {
    out <- c(...)
    if (is.null(out)) {
        print("FEFe")
    }
    if (TRUE %in% (c(NA, NaN, Inf, -Inf) %in% out)) {
        stop('events cannot contain NA, NaN, Inf')
    }
    if (length(unique(out)) != length(out)) {
        stop('events contain duplicate elements.')
    }
    if ("alpha" %in% out | "omega" %in% out) {
        stop('events may not contain elements titled "alpha" or "omega".')
    }
    class(out) <- c("events", "character")
    return(out)
}

#' @export 
print.events <- function(...) {
cat("Events object of",length(...), "elements:\n   ")
    cat(paste0(..., collapse = ", "))
    cat("\n")
}



#' Create a Sequences Object
#' 
#' Analogous to the \code{\link[base]{list}} function, a \code{sequences} object contains multiple \code{events} objects (see \code{\link[eratosthenes]{events}}).
#' 
#' @param ... objects of \code{events} class, or a \code{list} of \code{events} objects.

#' 
#' @examples
#' x <- events("A", "B", "C", "D", "E")
#' y <- events("B", "D", "F")
#' z <- events("A", "C", "F", "G")
#' sequences(x, y, z)
#' 
#' @returns A sequences object.
#' 
#' @export
sequences <- function(...) {
    UseMethod("sequences")
}
#' 
#' @rdname sequences
#' @export
sequences.events <- function(...) {
    quae_postea(...)
    quae_antea(...)

    out <- list(...)
    class(out) <- c("sequences", "list")
    return(out)
}
#' @rdname sequences
#' @export
sequences.list <- function(...) {
    out <- list(...)[[1]]
    chk <- sapply(out, inherits, "events")
    if (FALSE %in% chk) {
        stop("list input needs to contain events objects", call. = FALSE)
    }
    quae_postea(out)
    quae_antea(out)

    class(out) <- c("sequences", "list")
    return(out)
}

#' @export 
print.sequences <- function(...) {
    if (length(...) > 1) {
        cat("Sequences object of",length(...), "events object:\n   ")
    } else {
        cat("Sequences object of",length(...), "events objects:\n   ")

    }
    cat(paste0(..., collapse = "\n   "))
    cat("\n")
}



#' Sequence Diagnostic
#' 
#' If the creation of a \code{\link[eratosthenes]{sequences}} object has failed, this function checks all \code{\link[eratosthenes]{events}} objects for instances of disagreement, pairwise. The output will give the pairs of indices of the \code{events} which do not agree, as well as the most frequently attested \code{events} which are in disagreement.
#' 
#' @param ... Objects of \code{events} class.
#' 
#' @examples
#' u <- events("A", "D", "E")
#' v <- events("E", "D")
#' w <- events("B", "F", "C")
#' x <- events("A", "B", "C", "D", "E")
#' 
#' seq_diag(u, v, w, x)
#' 
#' a <- list(u, v, w, x)
#' seq_diag(a)
#' 
#' @returns A sequences object.
#' 
#' @export
seq_diag <- function(...) {
    UseMethod("seq_diag")
}

#' @rdname seq_diag
#' @export
seq_diag.events <- function(...) {
    x <- list(...)
    seq_diag()
}

#' @rdname seq_diag
#' @export
seq_diag.list <- function(...) {
    x <- list(...)[[1]]
    out <- list()
    k <- 1
    for (ki in 1:(length(x)-1)) {
        for (kj in (ki+1):length(x)) {
            #chk <- sequences(x[[i]], x[[j]])

            #####################

            obj <- list(x[[ki]], x[[kj]])
            elements <- unique(unlist(obj))
            M <- list()
            for (i in 1:length(obj)) {
                tmp <- numeric(length(obj[[i]]))
                for (j in 1:length(obj[[i]])) {
                    tmp[j] <- which(elements == obj[[i]][j])
                }
                M[[i]] <- tmp
            }

            mat <- quae_postea_matrix_cpp(length(elements), M)

            res <- list()
            for (i in 1:length(elements)) {
                res[[elements[i]]] <- c(elements[mat[i,] == 1], "omega")
            }

            chk <- TRUE
            for (i in names(res)) {
                if (i != "omega") {
                    if (i %in% res[[i]]) {
                        chk <- FALSE
                    }
                }
            }
            
            if (chk == FALSE) {
                out[[k]] <- c(ki, kj)
                k <- k + 1
            }
        }
    }

    most_discrepant <- rev(sort(table(unlist(out))))
    res <- list(pairs = out, most_discrepant = most_discrepant)
    class(res) <- c("seq_diag", "list")
    return(res)
}

#' @export 
print.seq_diag <- function(x) {
cat("Number of pairs of discrepant sequences:",length(x$pairs), "\nIndices of most frequent events objects in discrepant pairs (sorted in descending order):\n      ")
    cat(paste0(names(head(x$most_discrepant))), collapse = " ", "\n")
}



#' Create an Absolute Constraint Object
#' 
#' Analogous to the \code{\link[base]{list}} function, to create absolute constraint (\emph{terminus post quem} or \emph{ante quem}). The constraint must contain named elements of \code{"id"}, \code{"assoc"}, and \code{"samples"}, with an option to indicate \code{"type"}.
#' 
#' @param id a \code{character} object, giving a unique ID of the constraint
#' @param assoc the element within an \code{events} object to which the constraint is associated
#' @param type (optional) a \code{character} object, giving the type of constraint (e.g., a ceramic type, coin, radiocarbon dat). Mutiple types/subtypes/classes can be given as a vector. Default is \code{NULL}.
#' @param samples a vector of samples drawn from the appertaining probability density function of that constraint
#' 
#' @examples
#' # external constraints
#' coin1 <- absolute(id = "coin1", assoc = "B", samples = runif(100,-320,-300))
#' coin2 <- absolute(id = "coin2", assoc = "G", type = "RIC2 57", samples = seq(37, 41, length = 100))
#'   # seq(37, 41, length = 100) is equivalent in concept to runif(100, 37, 41))
#' destr <- absolute(id = "destr", assoc = "J", samples = 79)
#' 
#' coin1
#' coin2
#' destr
#' 
#' @returns An \code{absolute} object.
#' 
#' @export
absolute <- function(id, assoc, type = NULL, samples) {
    UseMethod("absolute")
}
#' 
#' @rdname absolute
#' @export
absolute.character <- function(id, assoc, type = NULL, samples) {
    if (!is.character(assoc)) {
        stop("Related context/event (assoc) must be character.")
    }
    if (!(is.null(type) | inherits(type, "character"))) {
        stop("Related type must be either NULL or character.")
    }
    if (!is.numeric(samples)) {
        stop("Samples must be numeric.")
    }
    out <- list(id = id, assoc = assoc, type = type, samples = samples)
    class(out) <- c("absolute", "list")
    return(out)
}



#' Create an Constraints Object
#' 
#' Analogous to the \code{\link[eratosthenes]{sequences}} function for relative events, this function collects one or more \code{\link[eratosthenes]{absolute}} objects into a single object, for input into \code{\link[eratosthenes]{gibbs_ad}}.
#' 
#' @param ... one or more \code{absolute} objects, or a \code{list} of \code{absolute} objects.
#' 
#' @examples
#' # external constraints
#' coin1 <- absolute(id = "coin1", assoc = "B", samples = runif(100,-320,-300))
#' coin2 <- absolute(id = "coin2", assoc = "G", type = "RIC2 57", samples = seq(37, 41, length = 100))
#'   # seq(37, 41, length = 100) is equivalent in concept to runif(100, 37, 41))
#' destr <- absolute(id = "destr", assoc = "J", samples = 79)
#' 
#' tpq <- constraints(coin1, coin2)
#' taq <- constraints(destr)
#' 
#' @returns A \code{constraints} object.
#' 
#' @export
constraints <- function(...) {
    UseMethod("constraints")
}
#' 
#' @rdname constraints
#' @export
constraints.absolute <- function(...) {
    out <- list(...)
    ids <- unlist(sapply(out, c)[1,])
    if (length(unique(ids)) != length(ids)) {
        stop("duplicate ids in constraints")
    }
    class(out) <- c("constraints", "list")
    return(out)
}
#'
#' @rdname constraints
#' @export
constraints.list <- function(...) {
    out <- list(...)[[1]]
    chk <- sapply(out, inherits, "absolute")
    if (FALSE %in% chk) {
        stop("non-events object in list input.")
    }
    ids <- unlist(sapply(out, c)[1,])
    if (length(unique(ids)) != length(ids)) {
        stop("duplicate ids in constraints")
    }
    class(out) <- c("constraints", "list")
    return(out)
}



#' Create an Finds Object
#' 
#' Analogous to the \code{\link[base]{list}} function, to create an object of a find (e.g., artifact or other element) related to a particular context or event. The find must contain named elements of \code{"id"} and \code{"assoc"}, with optional inputs of \code{"type"} and \code{"residual"}, to be collected into a single object via the \code{\link[eratosthenes]{assemblage}} function.
#' 
#' @param id a \code{character} object, giving a unique ID of the find.
#' @param assoc the element (e.g., context) within an \code{events} object to which the find is associated.
#' @param type (optional) a \code{character} object, giving the type of constraint (e.g., a ceramic type, coin, radiocarbon dat). Mutiple types/subtypes/classes can be given as a vector. Default is \code{NULL}.
#' @param residual (optional) if \code{TRUE}, indicates that the object is residual to its associated event (\code{assoc}), e.g., had a final deposition to be regarded prior to its context. Supplying \code{residual = TRUE} will suppress it from the estimation of production, use, and depositional dates in the function \code{\link[eratosthenes]{gibbs_ad_type}}. Default is \code{FALSE}.
#' 
#' @examples
#' f1 <- finds(id = "find01", assoc = "D", type = c("type1", "form1"))
#' f2 <- finds(id = "find02", assoc = "E", type = c("type1", "form2"))
#' f3 <- finds(id = "find03", assoc = "G", type = c("type1", "form1"), residual = TRUE)
#' f4 <- finds(id = "find04", assoc = "H", type = c("type2", "form1"))
#' f5 <- finds(id = "find05", assoc = "I", type = "type2")
#' f6 <- finds(id = "find06", assoc = "H")
#' 
#' @returns A \code{finds} object.
#' 
#' @export
finds <- function(id, assoc, type = NULL, residual = FALSE) {
    UseMethod("finds")
}
#' 
#' @rdname finds
#' @export
finds.character <- function(id, assoc, type = NULL, residual = FALSE) {
    if (!is.character(assoc)) {
        stop("Related context/event (assoc) must be character.")
    }
    if (!(is.null(type) | is.character(type))) {
        stop("Related type must be either NULL or character.")
    }
    if (!(residual == FALSE | residual == TRUE)) {
        stop("residual must either be TRUE or FALSE.")
    }
    out <- list(id = id, assoc = assoc, type = type, residual = residual)
    class(out) <- c("finds", "list")
    return(out)
}



#' Create an Assemblage Object
#' 
#' Analogous to the \code{\link[eratosthenes]{sequences}} or \code{\link[eratosthenes]{constraints}} function for relative and absolute events, this function collects one or more \code{\link[eratosthenes]{finds}} objects into a single object, for input into \code{\link[eratosthenes]{gibbs_ad_type}}.
#' 
#' @param ... one or more \code{finds} objects.
#' 
#' @examples
#' f1 <- finds(id = "find01", assoc = "D", type = c("type1", "form1"))
#' f2 <- finds(id = "find02", assoc = "E", type = c("type1", "form2"))
#' f3 <- finds(id = "find03", assoc = "G", type = c("type1", "form1"), residual = TRUE)
#' f4 <- finds(id = "find04", assoc = "H", type = c("type2", "form1"))
#' f5 <- finds(id = "find05", assoc = "I", type = "type2")
#' f6 <- finds(id = "find06", assoc = "H")
#' 
#' finds_all <- assemblage(f1, f2, f3, f4, f5, f6)

#' @returns An \code{assemblage} object.
#' 
#' @export
assemblage <- function(...) {
    UseMethod("assemblage")
}
#' 
#' @rdname assemblage
#' @export
assemblage.finds <- function(...) {
    out <- list(...)
    ids <- unlist(sapply(out, c)[1,])
    if (length(unique(ids)) != length(ids)) {
        stop("duplicate ids in finds.")
    }

    class(out) <- c("assemblage", "list")
    return(out)
}
#' 
#' @rdname assemblage
#' @export
assemblage.list <- function(...) {
    out <- list(...)[[1]]
    chk <- sapply(out, inherits, "finds")
    if (FALSE %in% chk) {
        stop("non-finds object in list input.")
    }
    ids <- unlist(sapply(out, c)[1,])
    if (length(unique(ids)) != length(ids)) {
        stop("duplicate ids in finds.")
    }

    class(out) <- c("assemblage", "list")
    return(out)
}



#' @export 
print.absolute <- function(...) {
cat("Absolute constraint object:\n   ")
    cat(paste0(..., collapse = ", "))
    cat("\n")
}

#' @export 
print.finds <- function(...) {
cat("Finds object:\n   ")
    cat(paste0(..., collapse = ", "))
    cat("\n")
}

#' @export 
print.constraints <- function(...) {
    if (length(...) > 1) {
        cat("Absolute constraints object of",length(...), "events:\n   ")
    } else {
        cat("Absolute constraints object of",length(...), "events:\n   ")

    }
    cat(paste0(..., collapse = "\n   "))
    cat("\n")
}

#' @export 
print.assemblage <- function(...) {
    if (length(...) > 1) {
        cat("Assemblage object of",length(...), "finds object:\n   ")
    } else {
        cat("Assemblage object of",length(...), "finds objects:\n   ")

    }
    cat(paste0(..., collapse = "\n   "))
    cat("\n")
}


