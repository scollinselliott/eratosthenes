#' Sequence Check
#'
#' For a \code{list} of partial sequences (of \code{vector} objects), check to see that joint elements of each occur the same order. That is, for two sequences with elements \eqn{A, B, C, D, E} and \eqn{B, D, F, E}, all joint elements must occur in the same order to pass the check. Two sequences \eqn{A, B, C, D, E} and \eqn{A, F, D, C, E} would not pass this check as the elements \eqn{C} and \eqn{D} occur in different orders in either sequence.

#' @param obj A \code{list} of \code{vector} objects which reperesent a sequence.    
#' @examples 
#' x <- c("A", "B", "C", "D", "E")
#' y <- c("B", "D", "F", "E")
#' a <- list(x, y)
#' 
#' seq_check(a)
#' 
#' z <- c("B", "F", "C")
#' b <- list(x, y, z)
#' 
#' seq_check(b)
#' 
#' @returns \code{TRUE} or \code{FALSE}
#' 
#' @export
seq_check <- function(obj) {
    UseMethod("seq_check")
}

#' @rdname seq_check
#' @export
seq_check.list <- function(obj) {
    qp_ <- quae_postea(obj)

    clear <- TRUE
    for (i in names(qp_)) {
        if (i %in% qp_[[i]]) {
            clear <- FALSE
        }
    }

    return(clear)
} 



#' Synthetic Ranking
#'
#' Using a \code{list} two or more partial sequences, all of which observe the same order of elements, create a single "synthetic" ranking. This is accomplished by counting the total number of elements after running a recursive trace through all partial sequences (via \code{\link[eratosthenes]{quae_postea}}). If partial sequences are inconsistent in their rankings, a \code{NULL} value is returned.
#'
#' @param obj A \code{list} of \code{vector} objects which reperesent a sequence.    
#' @param ties The way in which ties are handled per the \code{\link{rank}} function. The default is \code{"ties = average"}.
#' 
#' @examples 
#' x <- c("A", "B", "C", "D", "E")
#' y <- c("B", "D", "F", "E")
#' a <- list(x, y)
#' 
#' synth_rank(a)
#' 
#' @returns A single vector containing the synthesized ranking.
#' 
#' @export
synth_rank <- function(obj, ties = "average") {
    UseMethod("synth_rank")
}

#' @rdname synth_rank
#' @export
synth_rank.list <- function(obj, ties = "average") {
    result <- NULL
    if (seq_check(obj) == TRUE) {
        elements <- names(obj)
        qp_ <- quae_postea(obj)

        quot_postea <- numeric(length(elements))
        names(quot_postea) <- elements
        for (i in names(qp_)) {
            quot_postea[i] <- length(qp_[[i]])
        }
        result <- rank(quot_postea * -1, ties.method = ties)
        result <- names(result)[order(result)]
    } else {
        message("Sequences are inconsistent.")
    }
    return(result)
}



#' Quae Postea
#'
#' For a \code{list} of multple partial sequences (of \code{vector} objects), generate another \code{list} which, for each element, gives all elements that occur after it ("\emph{quae postea}"). This is analogous to a recursive trace through all partial sequences from left to right. A final element \code{"omega"} is added to all sets to avoid empty vectors. See also \code{\link[eratosthenes]{quae_antea}}.
#'
#' @param obj A \code{list} of \code{vector} objects which reperesent ordered sequences.    
#' 
#' @examples 
#' x <- c("A", "B", "C")
#' y <- c("B", "D", "E", "C", "F")
#' z <- c("C", "G")
#' a <- list(x, y, z)
#' 
#' quae_postea(a)
#' 
#' @returns A \code{list} of \code{vector} objects, which contain the elements that occur after any one given element in the input sequences. 
#' 
#' @export
quae_postea <- function(obj) {
    UseMethod("quae_postea")
}
#' 
#' @rdname quae_postea
#' @export
quae_postea.list <- function(obj) {
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
#' For a \code{list} of multple partial sequences (of \code{vector} objects), generate another \code{list} which, for each element, gives the elements that occur before it ("\emph{quae antea}"). This is analogous to a recursive trace through all partial sequences from right to left. An element \code{"alpha"} is added to all sets to avoid empty vectors. See also \code{\link[eratosthenes]{quae_postea}}.
#'
#' @param obj A \code{list} of \code{vector} objects which reperesent ordered sequences.    
#'
#' @examples 
#' x <- c("A", "B", "C")
#' y <- c("B", "D", "E", "C", "F")
#' z <- c("C", "G")
#' a <- list(x, y, z)
#' 
#' quae_antea(a)
#' 
#' @returns A \code{list} of \code{vector} objects, which contain the elements that occur before any one given element in the input sequences. 
#'
#' @export
quae_antea <- function(obj) {
    UseMethod("quae_antea")
}
#' 
#' @rdname quae_antea
#' @export
quae_antea.list <- function(obj) {
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
#' @param input A vector of elements in a sequence.
#' @param target A vector of elements in a sequence, containing at least three of the same elements as \code{input}.
#' 
#' @examples 
#' x <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J") # the input sequence
#' y <- c("D", "A", "J") # the target sequence
#' 
#' seq_adj(x, y)
#' 
#' @returns A vector of the adjusted sequence.
#' 
#' @export
seq_adj <- function(input, target) {
    UseMethod("seq_adj")
}
#' 
#' @rdname seq_adj
#' @export
seq_adj.character <- function(input, target) {
    result <- NULL
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
        result <- interp$y[1:length(x_pos)+1]
        names(result) <- input
        result <- input[order(result)]
    } else {
        message("Insufficient number of joint elements in input and target sequence (must be > 2).")
    }
    return(result)
}



#' Gibbs Sampler for Archaeological Dating
#'
#' A Gibbs sampler for archaeological dating, to fit relative sequences to absolute, calendrical dates, along with rule-based production dates of artifact types. Relative events can be associated with \emph{termini post quos} (\emph{t.p.q.}) and \emph{termini ante quos} (\emph{t.a.q.}), which are entered as samples from a given probability density function \eqn{f(t)}. This function may take any form, a single date (i.e., with a probability of 1), a continuous uniform distribution (any time between two dates), or a bespoke density (as with calibrated radicarbon dates). Relative events are modeled on a continuous uniform density between the latest antecedent event and earliest subsequent event.
#' 
#' Gibbs sampling is a coventional method for calibrating and estimating radiocarbon dates in light of absolute constraints and relative sequences: see \insertCite{buck_bayesian_1996,buck_bcal_1999,bronk_ramsey_bayesian_2009;textual}{eratosthenes}, the latter of which uses a mixture of Metropolis-Hastings and Gibbs.
#' 
#' In this implementation, two phases of Gibbs sampling are performed: an initial phase for selecting starting values and then the main sampler, with convergence evaluated using Monte Carlo standard errors (MCSE).
#' 
#' The initial Gibbs sampler results in a vector of starting values randomly sampled for each event up to \eqn{\sqrt{k}} runs, where \eqn{k} is the total number of events. Starting values may therefore take some time to assign, but this inital sampling is necessary to avoid a catastrophic collapse due to floating point errors in the initial selection of random values and will also result in closer starting values with respect to marginal densities. 
#' 
#' The main Gibbs sampler uses consistent batch means (CBM) determine convergence and hence when to end the main sampling run: there is no motivation to remove burn-in from the main sampling run nor to run multiple chains. CBM is assured to converge in distribution, see \insertCite{jones_fixed-width_2006,flegal_markov_2008;textual}{eratosthenes}. A stopping point for the main sampler is therefore determined using the mean of the Monte Carlo standard errors (MCSE) across all random variates, which is the input of \code{mcse_crit} (the mean MCSE for all events). The input \code{max_samples} indicates the maximum number of simulations to run, but the sampler will stop if the specified criterion of \code{mcse_crit} is passed. The default mean MCSE is set at \code{msce_crit = 0.5}, as the MCSE is measured in years (i.e. to allow for an error +/- 1 year), but, to be sure, individual events will have higher or lower MCSE than this mean criterion, whose primary purpose is as a stopping rule. Depending on the conditional structure of the relative sequences and the timescale of investigation, higher or lower MCSE may be more desireable or acceptable.
#' 
#' For the use dates of artifact types, see the \code{\link[eratosthenes]{use_dates}} function.
#'
#' @param sequences A \code{list} of relative sequences of elements (e.g., contexts).
#' @param finds Optional. A \code{list} of finds related to (contained in) the elements of \code{sequences}.
#' @param max_samples Maximum number of samples to run. Default is \code{10^5}.
#' @param size The number of samples to take on each iteration of the main Gibbs sampler. Default is \code{10^3}. 
#' @param mcse_crit Criterion for the Monte Carlo standard error to stop the Gibbs sampler.
#' @param tpq A \code{list} containing \emph{termini post quos}. Each object in the list consists of:
#'   * \code{id} A \code{character} ID of the  \emph{t.p.q.}, such as a reference or number.
#'   * \code{assoc} The element in \code{code} to which the \emph{t.p.q.} is associated. 
#'   * \code{samples} A vector of samples drawn from the appertaining probability density function of that \emph{t.p.q.}
#' @param taq A \code{list} containing \emph{termini ante quos}. Each object in the list consists of:
#'   * \code{id} A \code{character} ID of the  \emph{t.a.q.}, such as a reference or number.
#'   * \code{assoc} The element in \code{code} to which the \emph{t.p.q.} is associated. 
#'   * \code{samples} A vector of samples drawn from the appertaining probability density function of that \emph{t.p.q.}
#' @param alpha_ An initial \emph{t.p.q.} to limit any elements which may occur before the first provided \emph{t.p.q.} Default is \code{-5000}.
#' @param omega_ A final \emph{t.a.q.} to limit any elements which may occur after the after the last provided \emph{t.a.q.} Default is \code{1950}.
#' @param trim A logical value to determine whether elements that occur before the first \emph{t.p.q.} and after the last \emph{t.a.q.} should be ommitted from the results (i.e., to "trim" elements at the ends of the sequence, whose marginal densities depend on the selection of \code{alpha_} and \code{omega_}). Default is \code{TRUE}.
#' @param rule The rule for computing an estimated date of production of a find-type, either \code{"earliest"}, selecting a production date between the earliest deposition of that type and the next most earliest context, or \code{"naive"} (the default), which will select a production date any time between the distribution of that "earliest" date and the depositional date of that artifact.
#' 
#' @returns A \code{list} object of class \code{marginals} which contains the following:
#'    * \code{deposition} A \code{list} of samples from the marginal density of each context's depositional date.
#'    * \code{externals} A \code{list} of samples of the marginal density of each constrant (\emph{t.p.q.} and \emph{t.a.q.]}), as conditioned upon the occurrence of other depositional 
#'    * \code{production} If a \code{finds} object has been input, samples of the marginal density of the production date of finds types will be included in the output.
#'    * \code{mcse} The Monte Carlo standard errors (MCSE) of the random variates (fixed t.p./a.q. will have a MCSE of 0.)
#'
#' @examples
#' x <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J")
#' y <- c("B", "D", "G", "H", "K")
#' z <- c("F", "K", "L", "M")
#' contexts <- list(x, y, z)
#' 
#' f1 <- list(id = "find01", assoc = "D", type = c("type1", "form1"))
#' f2 <- list(id = "find02", assoc = "E", type = c("type1", "form2"))
#' f3 <- list(id = "find03", assoc = "G", type = c("type1", "form1"))
#' f4 <- list(id = "find04", assoc = "H", type = c("type2", "form1"))
#' f5 <- list(id = "find05", assoc = "I", type = "type2")
#' f6 <- list(id = "find06", assoc = "H", type = NULL)
#' 
#' artifacts <- list(f1, f2, f3, f4, f5, f6)
#'  
#' # external constraints
#' coin1 <- list(id = "coin1", assoc = "B", type = NULL, samples = runif(100,-320,-300)))
#' coin2 <- list(id = "coin2", assoc = "G", type = NULL, samples = seq(37, 41, length = 100))
#'   # seq(37, 41, length = 100) is equivalent in concept to runif(100, 37, 41)) 
#' destr <- list(id = "destr", assoc = "J", type = NULL, samples = 79)
#' 
#' tpq_info <- list(coin1, coin2)
#' taq_info <- list(destr)
#' 
#' result <- gibbs_ad(contexts, finds = artifacts, tpq = tpq_info, taq = taq_info)
#' 
#' @references
#'   \insertAllCited{}
#' 
#' @export
#' @importFrom Rdpack reprompt
gibbs_ad <- function(sequences, finds = NULL, max_samples = 10^5, size = 10^3, mcse_crit = 0.5, tpq = NULL, taq = NULL, alpha_ = -5000, omega_ = 1950, trim = TRUE, rule = "naive") {
    UseMethod("gibbs_ad")
}

#' @rdname gibbs_ad
#' @export
gibbs_ad.list <- function(sequences, finds = NULL, max_samples = 10^5, size = 10^3, mcse_crit = 0.5, tpq = NULL, taq = NULL, alpha_ = -5000, omega_ = 1950, trim = TRUE, rule = "naive") {
    if (seq_check(sequences) == TRUE) {
        proceed <- synth_rank(sequences)
 
        if (!is.list(tpq)) {
            tpq <- list(list(id = "tpq_default", assoc = proceed[1], type = NULL, samples = alpha_))
            warning("No tpq / tpq is not list object. Default of one tpq = -5000 used.")
        }
        if (!is.list(taq)) {
            taq <- list(list(id = "taq_default", assoc = proceed[length(proceed)], type = NULL, samples = omega_))
            warning("No taq / taq is not list object. Default of one taq = 1950 used.")
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

        message("Assigning initial random values (this may take a minute)...")

        gibbs[,1] <- gibbs_ad_initial_cpp(gibbs[,1], tpq_idx, PsiMatrix, tpq , taq_idx, PhiMatrix, taq, proceed_idx, init_sample)

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
            mcse <- mcse0[mcse0 > 0]

            message("Samples: ", ncol(gibbs), "     Mean MCSE: ",  round(mean(mcse),3))
            if (mean(mcse) < mcse_crit) {
                mcse_check <- TRUE
                message("MCSE criterion passed. Finishing.")
            } else {
                if (ncol(gibbs) >= max_samples) {
                    message("MC samples exceeded maximum stipulated without passing MCSE crterion. Finishing.")
                    mcse_check <- TRUE
                } else {

                    gibbs_next <- matrix(0, nrow = nrow(gibbs), ncol = (size + 1) )
                    gibbs_next[,1] <- gibbs[,ncol(gibbs)]
                    gibbs_next <- gibbs_ad_cpp(gibbs_next, tpq_idx, PsiMatrix, tpq , taq_idx, PhiMatrix, taq, proceed_idx)
                    gibbs <- cbind(gibbs, gibbs_next[, 2:ncol(gibbs_next)])
                }
            }
        }

        samples <- ncol(gibbs)

        deposition <- list()
        externals <- list()
        production <- list()

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
                    class(g) <- c("mc_samples", "numeric")
                    deposition[[iname]] <- g
                }
            } else {
                g <- gibbs[idx, ]
                class(g) <- c("mc_samples", "numeric")
                deposition[[iname]] <- g
            }
        }
        for (i in 1:length(tpq)) {
            ii <- tpq[[i]]
            tpq_idx <- match(ii$id, proceed_all)
            g <- gibbs[tpq_idx, ]
            class(g) <- c("mc_samples", "numeric")
            externals[[ii$id]] <- g
        }
        for (i in 1:length(taq)) {
            ii <- taq[[i]]
            taq_idx <- match(ii$id, proceed_all)
            g <- gibbs[taq_idx, ]
            class(g) <- c("mc_samples", "numeric")
            externals[[ii$id]] <- g
        }

        # production dates 
        if (!is.null(finds)) {
            message("Computing densities for find-type production...")
            findstypes <- c()
            for (i in finds) {
                findstypes <- c(findstypes, i$type)
            }
            for (i in tpq) {
                findstypes <- c(findstypes, i$type)
            }
            for (i in taq) {
                findstypes <- c(findstypes, i$type)
            }

            findslength <- length(finds)
            findstypes <- unique(findstypes)
            findstypeslength <- length(findstypes)

            attestation <- matrix(0, nrow = elements, ncol = findstypeslength)
            for (i in 1:findslength) {
                ii <- finds[[i]]
                context <- ii$assoc
                contexti <- match(context, proceed_all)
                types <- ii$type
                typeslength <- length(types)
                for (k in 1:typeslength) {
                    j <- match(types[k], findstypes)
                    attestation[contexti, j] <- 1
                }
            }
            for (i in 1:length(tpq)) {
                ii <- tpq[[i]]
                context <- ii$assoc
                contexti <- match(context, proceed_all)
                types <- ii$type
                typeslength <- length(types)
                for (k in 1:typeslength) {
                    j <- match(types[k], findstypes)
                    attestation[contexti, j] <- 1
                }
            }
            for (i in 1:length(taq)) {
                ii <- taq[[i]]
                context <- ii$assoc
                contexti <- match(context, proceed_all)
                types <- ii$type
                typeslength <- length(types)
                for (k in 1:typeslength) {
                    j <- match(types[k], findstypes)
                    attestation[contexti, j] <- 1
                }    
            }

            type_earliest_dep <- matrix(0, nrow = findstypeslength, ncol = samples)

            for (i in 1:findstypeslength) {
                attested <- attestation[ , i]
                contexts <- which(attested == 1)
                cols <- as.matrix( gibbs[contexts,] )
                if (ncol(cols) == 1) {
                    type_earliest_dep[i,] <- t(cols)
                } else {
                    for (j in 1:samples) {
                        type_earliest_dep[i,j] <- min(cols[,j])
                    }
                }
            }

            type_prev_dep <- matrix(0, nrow = findstypeslength, ncol = samples)

            for (i in 1:findstypeslength) {
                for (j in 1:samples) {
                    earliest_dep <- type_earliest_dep[i,j]
                    deps <- gibbs[ , j]
                    prev <- max(deps[deps < earliest_dep])
                    type_prev_dep[i, j] <- prev
                }
            }

            mcseprd <- numeric(findstypeslength)
            names(mcseprd) <- findstypes
            
            if (rule == "naive") {
                for (i in 1:findstypeslength) {
                    attested <- attestation[ , i]
                    contexts <- which(attested == 1)
                    cols <- as.matrix( gibbs[contexts,] )   
                    outsize <- nrow(cols) * ncol(cols)

                    if (ncol(cols) == 1) {
                        out <- numeric(outsize)

                        L <- type_prev_dep[i, ]
                        U <- t(cols)
                        out <- stats::runif(samples, L, U)
                    } else {
                        out <- matrix(0, nrow = nrow(cols), ncol = ncol(cols))

                        for (k in 1:nrow(cols)) {
                            for (j in 1:samples) {
                                L <- type_prev_dep[i, j]
                                U <- cols[k,j]
                                s <- stats::runif(1, L, U)
                                out[k , j] <- s
                           
                            }
                        }
                    }

                    g <- as.vector(out)
                    
                    production[[findstypes[i]]] <- g

                    n_upto <- length(g)
                    n_batch <- floor(sqrt(n_upto)) # length of samples in batch
                    K <- floor(n_upto / n_batch) # number of batches
                    m_batch <- numeric((K-1))
                    remainder <- n_upto - n_batch * K + 1

                    idx1 <- remainder
                    for (k in 1:(K-1)) {
                        idxs <- idx1:(idx1 + n_batch)       
                        m_batch[k] <- mean(g[idxs]) 
                        idx1 <- idxs[length(idxs)] + 1
                    }
                    mcseprd[i] <- sqrt( sum( (m_batch - mean(m_batch))^2 ) * (n_batch / (K-1) ) ) / sqrt((K - 1) * n_batch)

                }
            } else if (rule == "earliest") {
                for (i in 1:findstypeslength) {
                    attested <- attestation[ , i]
                    contexts <- which(attested == 1)
                    cols <- as.matrix( gibbs[contexts,] )   

                    out <- numeric(outsize)

                    L <- type_prev_dep[i, ]
                    U <- type_earliest_dep[i, ]
                    out <- stats::runif(samples, L, U)
                    
                    g <- out
                    class(g) <- c("mc_samples", "numeric")
                    production[[findstypes[i]]] <- g                
                    
                    n_upto <- length(g)
                    n_batch <- floor(sqrt(n_upto)) # length of samples in batch
                    K <- floor(n_upto / n_batch) # number of batches
                    m_batch <- numeric((K-1))
                    remainder <- n_upto - n_batch * K + 1

                    idx1 <- remainder
                    for (k in 1:(K-1)) {
                        idxs <- idx1:(idx1 + n_batch)       
                        m_batch[k] <- mean(g[idxs]) 
                        idx1 <- idxs[length(idxs)] + 1
                    }
                    mcseprd[i] <- sqrt( sum( (m_batch - mean(m_batch))^2 ) * (n_batch / (K-1) ) ) / sqrt((K - 1) * n_batch)
                }
            } else {
                production[[findstypes[i]]] <- NULL
                warning('Invalid rule, with NULL given for production dates. Options are "naive", "earliest".')
            }

        message("Finished.")
        mcse0 <- c(mcse0, mcseprd)
        result <- list(deposition = deposition, externals = externals, production = production, mcse = mcse0)
        class(result) <- c("marginals", "list")
        return(result)
        } else {
            result <- list(deposition = deposition, externals = externals, mcse = mcse0)
            class(result) <- c("marginals", "list")
            return(result)
        }

    } else {
        stop("Sequences has failed consistency check with seq_check().")
    }
}



#' @export 
print.marginals <- function(obj) {
    mcse <- obj$mcse
    mmcse <- mean(mcse[mcse > 0])
    cat("\n Marginals from joint conditional density, consisting of:\n",
    "    ", length(obj$deposition), "depositional events\n", 
    "    ", length(obj$externals), "external constraints (t.p./a.q.)\n",
    "    ", length(obj$production), "production dates\n\n",
        "Call the following objects to see element names: \n",
        "    deposition", "\n",
        "    externals", "\n",
        "    production\n\n",
        "Mean Monte Carlo standard error (MCSE) of random variates: \n",
    "    ", round(mmcse, 3), "\n\n\n")
}



#' @export
summary.marginals <- function(obj, events = NULL, digits = 2) {
    depmu <- sapply(obj$deposition, mean)
    depmcse <- obj$mcse[names(obj$deposition)]
    depdat <- data.frame(Mean = depmu, MCSE = depmcse)
    rownames(depdat) <- names(obj$deposition)

    extmu <- sapply(obj$externals, mean)
    extmcse <- obj$mcse[names(obj$externals)]
    extdat <- data.frame(Mean = extmu, MCSE = extmcse)
    rownames(extdat) <- names(obj$externals)

    prdmu <- sapply(obj$production, mean)
    prdmcse <- obj$mcse[names(obj$production)]
    prddat <- data.frame(Mean = prdmu, MCSE = prdmcse)
    rownames(prddat) <- names(obj$production)

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
plot.marginals <- function(obj, events = NULL, col = NULL, opac = 1) {
    if (is.vector(events)) {
        master <- c(obj$deposition, obj$externals, obj$production)
        dat <- data.frame(idx = 0, x = 0, event = 0)

        for (i in 1:length(events)) {
            item <- events[i]
            if (!(item %in% names(master))) {
                stop('One or more event names not contained in marginals.')
            }
            dat <- rbind(dat, data.frame(idx = 1:length(master[[item]]), x = master[[item]], event = item))
        }  
    } else {
        stop('Events must be vector object.')
    }
    dat <- dat[2:nrow(dat), ]
    if (!is.null(col)) {
        if (length(events) == length(col)) {
            ggplot2::ggplot(dat, ggplot2::aes(x = idx, y = x, color = event), group = event) + ggplot2::geom_line(alpha = opac) +  ggplot2::scale_color_manual(values = col) + ggplot2::theme_bw() } else {
            stop("Number of colors not equal to number of events.")
        }        
    } else { 
    ggplot2::ggplot(dat, ggplot2::aes(x = idx, y = x, color = event), group = event) + ggplot2::geom_line(alpha = opac) + ggplot2::theme_bw() 
    }
}



#' @export
hist.marginals <- function(obj, events = NULL, col = NULL, opac = 1) {
    if (is.vector(events)) {
        master <- c(obj$deposition, obj$externals, obj$production)
        dat <- data.frame(idx = 0, x = 0, event = 0)
        iqrs <- numeric(length(events))
        for (i in 1:length(events)) {
            item <- events[i]
            if (!(item %in% names(master))) {
                stop('One or more event names not contained in marginals.')
            }
            iqrs[i] <- stats::IQR(master[[item]])
            dat <- rbind(dat, data.frame(idx = 1:length(master[[item]]), x = master[[item]], event = item))
        }  
    } else {
        stop('Events must be vector.')
    }
    dat <- dat[2:nrow(dat), ]
    y <- nrow(dat)/length(events)
    fd <- 2 * mean(iqrs) / y^(1/3)

    if (!is.null(col)) {
        if (length(events) == length(col)) {
            ggplot2::ggplot(dat, ggplot2::aes(x = x, y = ggplot2::after_stat(stats::density), fill = event)) + ggplot2::geom_histogram(binwidth = fd, alpha = opac, position = 'identity') +  ggplot2::scale_fill_manual(values = col) + ggplot2::theme_bw() } else {
            stop("Number of colors not equal to number of events.")
        }        
    } else { 
    ggplot2::ggplot(dat, ggplot2::aes(x = x, y = ggplot2::after_stat(stats::density), fill = event)) +
        ggplot2::geom_histogram(binwidth = fd, alpha = opac, position = 'identity') + ggplot2::theme_bw()
    }
}









#' Use Date
#'
#' Using the results of \code{\link[eratosthenes]{gibbs_ad}}, estimate a single density for the date of use of an artifact or artifact type. Multiple artifacts and types can be given, which will be pooled into a single estimation of production, use, and deposition. For example, one can input several individual finds via their id number as comprising a type, or multiple (sub)types/classes as a single type, (e.g., "MGS V amphora" and "MGS VI amphora" to construct a type "MGS V/VI amphora"). Depending on whether one is using id numbers or type(s), the \code{id} or \code{type} argument is used, which takes a vector of the entries' names. The \code{use_date} function samples a use date between the production and depositional densities from the results of \code{\link[eratosthenes]{gibbs_ad}}, and in turn pools those densities for the production and deposition of the stipulate type.
#' 
#' See \code{\link[eratosthenes]{gibbs_ad}} for information on consistent batch means and Monte Carlo standard error, whichj are used to determined convergence.
#'
#' @param gibbs A \code{list} object of class \code{marginals}, the output of \code{\link[eratosthenes]{gibbs_ad}}.
#' @param finds Either the \code{list} object of finds used as input to produce \code{marginals} or a \code{data.frame} of two columns, the first listing the context and the second the incidence of the type in that context.
#' @param id A vector of the \code{id} of one or more specific finds whose use date is to be estimated. The values of \code{id} must match those in the \code{list} of \code{finds}. If \code{type} is used, \code{id} is ignored.
#' @param type A vector of one or more types to estimate a use density for. Must contain a value if \code{id} is \code{NULL}.
#' @param max_samples Maximum number of samples to run. Default is \code{10^5}.
#' @param size The number of samples to take on each iteration of the main Gibbs sampler. Default is \code{10^3}. 
#' @param mcse_crit Criterion for the Monte Carlo standard error to stop the Gibbs sampler.
#' 
#' @examples 
#' x <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J")
#' y <- c("B", "D", "G", "H", "K")
#' z <- c("F", "K", "L", "M")
#' contexts <- list(x, y, z)
#' 
#' f1 <- list(id = "find01", assoc = "D", type = c("type1", "form1"))
#' f2 <- list(id = "find02", assoc = "E", type = c("type1", "form2"))
#' f3 <- list(id = "find03", assoc = "G", type = c("type1", "form1"))
#' f4 <- list(id = "find04", assoc = "H", type = c("type2", "form1"))
#' f5 <- list(id = "find05", assoc = "I", type = "type2")
#' f6 <- list(id = "find06", assoc = "H", type = NULL)
#' 
#' artifacts <- list(f1, f2, f3, f4, f5, f6)
#'  
#' # external constraints
#' coin1 <- list(id = "coin1", assoc = "B", type = NULL, samples = runif(100,-320,-300))
#' coin2 <- list(id = "coin2", assoc = "G", type = NULL, samples = runif(100,37,41))
#' destr <- list(id = "destr", assoc = "J", type = NULL, samples = 79)
#' 
#' tpq_info <- list(coin1, coin2)
#' taq_info <- list(destr)
#' 
#' result <- gibbs_ad(contexts, finds = artifacts, tpq = tpq_info, taq = taq_info)
#' 
#' # use dates by specifying ids
#' use_dates(result, artifacts, id = c("find04", "find05"))
#'
#' # use dates by speciifying types
#' use_dates(result, artifacts, type = "type1")
#' 
#' @returns A \code{list} of class \code{use_marginals} densities, conditional upon the production and depositional dates.
#' 
#' @export
use_dates <- function(gibbs, finds, id = NULL, type = NULL, max_samples = 10^5, size = 10^3, mcse_crit = 0.5) {
    UseMethod("use_dates")
}
#' 
#' @rdname use_dates
#' @export
use_dates.marginals <- function(gibbs, finds, id = NULL, type = NULL, max_samples = 10^5, size = 10^3, mcse_crit = 0.5) {
    if (class(finds) == "data.frame") {
        finds <- finds_d2l(finds)
    }
    if (!(class(finds) %in% c("data.frame", "list"))) {
        stop("finds must be list or data frame object.")
    }

    if (is.null(id) & is.null(type)) {
        stop("Either one or more id or types must be specified.")
    }
    if (length(type) > 0) {
        id <- ids_of_types(finds, type)
    }
    if (is.null(id)) {
        stop("id or type not present in finds input.")
    }

    message("Estimating use distribution for id/type specified, by pooling samples in between dates of production and deposition.")
    sequence_ <- list()
    tpq_ <- list()
    taq_ <- list()
    J <- 1
    K <- 1
    for (i in 1:length(id)) {
    sequence_[[i]] <- id[i]
        for (j in finds) {
            if (j$id == id[i]) {
                k_dep <- j$assoc
                taq_[[J]] <- list(id = paste0(i, " ", j$id, " ", k_dep),
                                assoc = j$id,
                                samples = gibbs$deposition[[k_dep]] )
                
                k_type <- j$type
                for (k in 1:length(k_type)) {
                    tpq_[[K]] <- list(id = paste0(i, " ", j$id, " ", id[i], " ", k_type[k]), 
                                        assoc = j$id,
                                        samples = gibbs$production[[k_type[k]]] )
                    K <- K + 1
                }
                J <- J + 1
            }
        }
    }
    print(sequence_)
    if (length(sequence_) == 0) {
        stop("No ids or types found with that name.")
    }
    result0 <- gibbs_ad(sequence_, tpq = tpq_, taq = taq_, max_samples = max_samples, size = size, mcse_crit = mcse_crit)
    use_dates0 <- result0$deposition
    mcse0 <- result0$mcse

    prd <- c()
    dep <- c()

    for (i in tpq_) {
        prd <- c(prd, i$samples)
    }
    for (i in taq_) {
        dep <- c(dep, i$samples)
    }


    mcse <- mcse0[names(use_dates0)]
    gibbs[['use']] <- list(use_name = type, use_date = use_dates0, use_mcse = mcse, item = item, production_date = prd, deposition_date = dep)
    class(gibbs) <- c("use_marginals", "list")

    return(gibbs)
} 




#' @export 
print.use_marginals <- function(obj) {
    mcse <- obj$mcse
    mmcse <- mean(mcse[mcse > 0])
    cat("\n Marginals from joint conditional density, consisting of:\n",
    "    ", length(obj$deposition), "depositional events\n", 
    "    ", length(obj$externals), "external constraints (t.p./a.q.)\n",
    "    ", length(obj$production), "production dates\n\n",
        "Call the following objects to see element names: \n",
        "    deposition", "\n",
        "    externals", "\n",
        "    production\n\n",
        "Mean Monte Carlo standard error (MCSE) of random variates: \n",
    "    ", round(mmcse, 3), "\n\n",
        "Marginals of use dates, taking production and deposition as fixed densities: \n",
    "    ", length(obj$use$use_date), "use events pooled:\n",
    "    ","    ",obj$use$use_name, "\n\n")

}



#' @export
hist.use_marginals <- function(obj, display_name = "Event", display = c("production", "use", "deposition"), col = NULL, opac = 1) {
    dat <- c(x = c(), event = c())
    iqrs <- c()
    n <- c()

    if (length(display) == 0) {
        stop("At least one or more event types (production, use, deposition) must be selected.")
    }

    dep_samples <- obj$use$deposition_date
    prd_samples <- obj$use$production_date

    for (i in 1:length(obj$use$use_date)) {
        use_samples <- obj$use$use_date[[i]]
    } 
    #     for (j in obj$finds) {
    #         if (j$id == id[i]) {
    #             k_dep <- j$assoc
    #             dep_samples <- c(prd_samples, gibbs$deposition[[k_dep]] )
    #             k_type <- j$type
    #             for (k in 1:length(k_type)) {
    #                 prd_samples <- c(prd_samples, gibbs$production[[k_type[k]]] )
    #             }
    #         }
    #     }
    # }

    if ("use" %in% display) {
        dat <- rbind(dat, data.frame(x = use_samples, event = paste0(display_name, "- Use")))
        iqrs <- c(iqrs, stats::IQR(use_samples))
        n <- c(n, length(use_samples))
    }
    if ("deposition" %in% display) {
        dat <- rbind(dat, data.frame(x = dep_samples, event = paste0(display_name, "- Deposition")))
        iqrs <- c(iqrs, stats::IQR(dep_samples))
        n <- c(n, length(dep_samples))
    }
    if ("production" %in% display) {
        dat <- rbind(dat, data.frame(x = prd_samples, event = paste0(display_name, "- Production")))
        iqrs <- c(iqrs, stats::IQR(prd_samples))
        n <- c(n, length(prd_samples))
    }

    fd <- 2 * mean(iqrs) / mean(n)^(1/3)

    if (!is.null(col)) {
        if (length(events) == length(col)) {
            ggplot2::ggplot(dat, ggplot2::aes(x = x, y = ggplot2::after_stat(stats::density), fill = event)) + ggplot2::geom_histogram(binwidth = fd, alpha = opac, position = 'identity') +  ggplot2::scale_fill_manual(values = col) + ggplot2::theme_bw()
        } else {
            stop("Number of colors not equal to number of events.")
        }        
    } else { 
    ggplot2::ggplot(dat, ggplot2::aes(x = x, y = ggplot2::after_stat(stats::density), fill = event)) +
        ggplot2::geom_histogram(binwidth = fd, alpha = opac, position = 'identity') + ggplot2::theme_bw()
    }
}



#' Convert Finds List Object to Data Frame (Context / Find-Type)
#' 
#' Performs the opposite of \code{\link[eratosthenes]{finds_d2l}}. Takes a \code{list} object of finds and their types, used as input in \code{\link[eratosthenes]{gibbs_ad}]}, and returns a \code{data.frame} of two columns, containing the context in the first and the find-type in the second, and the \code{id} of the object in the third.
#' 
#' @param input A list object of finds (each one a list) of associated contexts and types. 
#' 
#' @examples 
#' f1 <- list(id = "find01", assoc = "D", type = c("type1", "form1"))
#' f2 <- list(id = "find02", assoc = "E", type = c("type1", "form2"))
#' f3 <- list(id = "find03", assoc = "G", type = c("type1", "form1"))
#' f4 <- list(id = "find04", assoc = "H", type = c("type2", "form1"))
#' f5 <- list(id = "find05", assoc = "I", type = "type2")
#' f6 <- list(id = "find06", assoc = "H", type = NULL)
#' 
#' artifacts <- list(f1, f2, f3, f4, f5, f6)
#' 
#' # convert list to data frame
#' artifacts_df <- finds_l2d(artifacts)
#' 
#' @returns A three-column data frame of contexts (first column) and find-types attested in that context (second column), along with the id number (third column).
#' 
#' @export
finds_l2d <- function(input) {
    UseMethod("finds_l2d")
}
#' 
#' @rdname finds_l2d
#' @export
finds_l2d.list <- function(input) {
    result <- data.frame(Context = c(), Type = c())
    for (i in input) {
        if (!is.na(sum(match(c("id", "assoc", "type"), names(i))))) {
            if (length(i$type) > 0) {
                for (j in i$type) {
                    result <- rbind(result, data.frame(Context = i$assoc, Type = j, Id = i$id ) )
                }
            }
        } else {
            stop("Finds list object does not contain correct headings (id, assoc, type).")
        }
    }
    return(result)
}



#' Convert Finds Data Frame (Context / Find-Type) to List Object
#' 
#' Performs the opposite of \code{\link[eratosthenes]{finds_l2d}}. Takes a \code{data.frame} object of two columns, containing the context in the first and the find-type in the second, and returns a \code{list} object for input in \code{\link[eratosthenes]{gibbs_ad}]}. The value of the find \code{id} is automatically generated as an integer if not provided in a third column.
#' 
#' @param input A two-column data frame of contexts (first column) and find-types (second column). An optional third column of an id number may be provided.
#' 
#' @examples 
#' f1 <- list(id = "find01", assoc = "D", type = c("type1", "form1"))
#' f2 <- list(id = "find02", assoc = "E", type = c("type1", "form2"))
#' f3 <- list(id = "find03", assoc = "G", type = c("type1", "form1"))
#' f4 <- list(id = "find04", assoc = "H", type = c("type2", "form1"))
#' f5 <- list(id = "find05", assoc = "I", type = "type2")
#' f6 <- list(id = "find06", assoc = "H", type = NULL)
#' 
#' artifacts <- list(f1, f2, f3, f4, f5, f6)
#' 
#' # convert list to data frame
#' artifacts_df <- finds_l2d(artifacts)
#' 
#' # convert data frame to list
#' artifacts_list <- finds_d2l(artifacts_df)
#' 
#' @returns A list of finds (each one a list) associated with contexts and their types.
#' 
#' @export
finds_d2l <- function(input) {
    UseMethod("finds_d2l")
}
#' 
#' @rdname finds_d2l
#' @export
finds_d2l.data.frame <- function(input) {
    result <- list()
    if (ncol(input) == 2) {
        for (i in 1:nrow(input)) {
            result[[i]] <- list(id = as.character(i), assoc = input[i,1], type = input[i,2])
        }
    } else if (ncol(input) == 3) {
        for (i in 1:nrow(input)) {
            result[[i]] <- list(id = input[i,3], assoc = input[i,1], type = input[i,2])
        }
    } else {
        stop("data frame must be 2 or three columns.")
    }

    return(result)
}



#' Ids of Types
#' 
#' Given a \code{list} object of finds (with keys of \code{id}, \code{assoc}, \code{type} in each entry), return a vector of the \code{id} elements that belong to one or more specified type.
#' 
#' @param input A \code{list} object whose elements are a list containing the keys of \code{id}, \code{assoc}, \code{type}.
#' @param type A vector or element 
#' 
#' @examples 
#' f1 <- list(id = "find01", assoc = "D", type = c("type1", "form1"))
#' f2 <- list(id = "find02", assoc = "E", type = c("type1", "form2"))
#' f3 <- list(id = "find03", assoc = "G", type = c("type1", "form1"))
#' f4 <- list(id = "find04", assoc = "H", type = c("type2", "form1"))
#' f5 <- list(id = "find05", assoc = "I", type = "type2")
#' f6 <- list(id = "find06", assoc = "H", type = NULL)
#' 
#' artifacts <- list(f1, f2, f3, f4, f5, f6)
#' 
#' ids_of_types(artifacts, type = "type1")
#' ids_of_types(artifacts, type = c("type1", "type2"))
#' 
#' @returns A vector of ids within a \code{list} object of \code{finds} class, 
#' 
#' @export
ids_of_types <- function(input, type = NULL) {
    UseMethod("ids_of_types")
}
#' 
#' @rdname ids_of_types
#' @export
ids_of_types.list <- function(input, type = NULL) {
    result <- c()
    if (is.vector(type)) {
        for (k in type) {
            for (i in input) {
                if (!is.na(sum(match(c("id", "assoc", "type"), names(i))))) {
                    if (length(i$type) > 0) {
                        for (j in i$type) {
                            if (j == k) {
                                result <- c(result, i$id)
                            }
                        }
                    }
                } else {
                    stop("Finds list object does not contain correct headings (id, assoc, type).")
                }
            }
        }
    return(result)
    } else {
        stop("type is not vector object.")
    }
}
            




