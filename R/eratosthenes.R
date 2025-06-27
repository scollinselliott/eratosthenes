#' Sequence Check
#'
#' For a \code{list} of partial sequences (of \code{vector} objects), check to see that joint elements of each occur the same order. That is, for two sequences with elements \eqn{A, B, C, D, E} and \eqn{B, D, F, E}, all joint elements must occur in the same order to pass the check. Two sequences \eqn{A, B, C, D, E} and \eqn{A, F, D, C, E} would not pass this check as the elements \eqn{C} and \eqn{D} occur in different orders in either sequence.
#' 
#' Event names \code{alpha} and \code{omega} are reserved for the ultimate boundaries of the chronological framework and should not be used in naming events in sequences.
#' 
#' @param obj A \code{list} of \code{vector} objects which represent a sequence.    
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
        if (i != "omega") {
            if (i %in% qp_[[i]]) {
                clear <- FALSE
            }
        }
    }

    return(clear)
} 



#' Synthetic Ranking
#'
#' Using a \code{list} two or more partial sequences, all of which observe the same order of elements, create a single "synthetic" ranking. This is accomplished by counting the total number of elements after running a recursive trace through all partial sequences (via \code{\link[eratosthenes]{quae_postea}}). If partial sequences are inconsistent in their rankings, a \code{NULL} value is returned.
#'
#' @param obj A \code{list} of \code{vector} objects which represent a sequence.    
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
#' For a \code{list} of multiple partial sequences (of \code{vector} objects), generate another \code{list} which, for each element, gives all elements that occur after it ("\emph{quae postea}"). This is analogous to a recursive trace through all partial sequences from left to right. A final element \code{"omega"} is added to all sets to avoid empty vectors. See also \code{\link[eratosthenes]{quae_antea}}.
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
#' For a \code{list} of multiple partial sequences (of \code{vector} objects), generate another \code{list} which, for each element, gives the elements that occur before it ("\emph{quae antea}"). This is analogous to a recursive trace through all partial sequences from right to left. An element \code{"alpha"} is added to all sets to avoid empty vectors. See also \code{\link[eratosthenes]{quae_postea}}.
#'
#' @param obj A \code{list} of \code{vector} objects which represent ordered sequences.    
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



#' Gibbs Sampler for Archaeological Dates
#'
#' A Gibbs sampler for dating archaeological events, to fit relative sequences to absolute, calendrical dates, along with rule-based production dates of artifact types. Relative events can be associated with \emph{termini post quos} (\emph{t.p.q.}) and \emph{termini ante quos} (\emph{t.a.q.}), which are entered as samples from a given probability density function \eqn{f(t)}. This function may take any form, a single date (i.e., with a probability of 1), a continuous uniform distribution (any time between two dates), or a bespoke density (as with calibrated radiocarbon dates). Relative events are modeled on a continuous uniform density between the latest antecedent event and earliest subsequent event.
#' 
#' Gibbs sampling is a conventional method for calibrating and estimating radiocarbon dates in light of absolute constraints and relative sequences: see \insertCite{buck_bayesian_1996,buck_bcal_1999,bronk_ramsey_bayesian_2009;textual}{eratosthenes}, the latter of which uses a mixture of Metropolis-Hastings and Gibbs.
#' 
#' In this implementation, two phases of Gibbs sampling are performed: an initial phase for selecting starting values and then the main sampler, with convergence evaluated using Monte Carlo standard errors (MCSE).
#' 
#' The initial Gibbs sampler results in a vector of starting values randomly sampled for each event up to \eqn{\sqrt{k}} runs, where \eqn{k} is the total number of events. Starting values may therefore take some time to assign, but this initial sampling is necessary to avoid a catastrophic collapse due to floating point errors in the initial selection of random values and will also result in closer starting values with respect to marginal densities. 
#' 
#' The main Gibbs sampler uses consistent batch means (CBM) determine convergence and hence when to end the main sampling run: there is no motivation to remove burn-in from the main sampling run nor to run multiple chains. CBM is assured to converge in distribution, see \insertCite{jones_fixed-width_2006,flegal_markov_2008;textual}{eratosthenes}. A stopping point for the main sampler is therefore determined using the mean of the Monte Carlo standard errors (MCSE) across all random variates, which is the input of \code{mcse_crit} (the mean MCSE for all events). The input \code{max_samples} indicates the maximum number of simulations to run, but the sampler will stop if the specified criterion of \code{mcse_crit} is passed. The default mean MCSE is set at \code{mcse_crit = 0.5}, as the MCSE is measured in years (i.e. to allow for an error +/- 1 year), but, to be sure, individual events will have higher or lower MCSE than this mean criterion, whose primary purpose is as a stopping rule. 
#' 
#' Note that the MCSE criterion is applied as a stopping rule for depositional dates and external constraints. The number of Monte Carlo samples for production dates of types is chosen to be identical to that need to pass \code{mcse_crit}, such that ultimately the final mean MCSE of all variates may differ from that of the criterion. Depending on the conditional structure of the relative sequences and the timescale of investigation, higher or lower MCSE may be more desirable or acceptable.
#' 
#' For the use dates of artifact type production, use, and deposition, see the \code{\link[eratosthenes]{gibbs_ad_use}} function.
#'
#' @param sequences A \code{list} of relative sequences of elements (e.g., contexts).
#' @param finds Optional. A \code{list} of finds related to (contained in) the elements of \code{sequences}.
#' @param max_samples Maximum number of samples to run. Default is \code{10^5}.
#' @param size The number of samples to take on each iteration of the main Gibbs sampler. Default is \code{10^3}. 
#' @param mcse_crit Criterion for the Monte Carlo standard error to stop the Gibbs sampler, as based on depositional dates and absolute constraints. The number of Monte Carlo samples for production dates is identical to that depositional dates.
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
#' @param trim A logical value to determine whether elements that occur before the first \emph{t.p.q.} and after the last \emph{t.a.q.} should be omitted from the results (i.e., to "trim" elements at the ends of the sequence, whose marginal densities depend on the selection of \code{alpha_} and \code{omega_}). Default is \code{TRUE}.
#' @param rule The rule for computing an estimated date of production of a find-type, either \code{"earliest"}, selecting a production date between the earliest deposition of that type and the next most earliest context, or \code{"naive"} (the default), which will select a production date any time between the distribution of that "earliest" date and the depositional date of that artifact.
#' 
#' @returns A \code{list} object of class \code{marginals} which contains the following:
#'    * \code{deposition} A \code{list} of samples from the marginal density of each context's depositional date.
#'    * \code{externals} A \code{list} of samples of the marginal density of each constraint (\emph{t.p.q.} and \emph{t.a.q.]}), as conditioned upon the occurrence of other depositional 
#'    * \code{production} If a \code{finds} object has been input, samples of the marginal density of the production date of finds types will be included in the output. If types are attested in trimmed contexts, 
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
#' coin1 <- list(id = "coin1", assoc = "B", type = NULL, samples = runif(100,-320,-300))
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


        # sampling initial values

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
                    message("\nMC samples exceeded maximum stipulated without passing MCSE crterion. Finishing.")
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

        class(deposition) <- c("events", "list")
        class(externals) <- c("events", "list")
        class(production) <- c("events", "list")

        names_dep_ext <- c(names(deposition), names(externals), names(production))
        mcse0 <- mcse0[names_dep_ext]

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






#' Gibbs Sampler for Archaeological Dates: Artifact Use
#'
#' Using the results of \code{\link[eratosthenes]{gibbs_ad}}, estimate a single density for the date of use of an artifact or artifact type. Multiple artifacts and types can be given, which will be pooled into a single type. For example, one can input several individual finds via their id number as comprising a type, or multiple (sub)types/classes as a single type, (e.g., "MGS V amphora", "MGS VI amphora", "MGS V/VI amphora" to construct one group).
#' 
#' Depending on whether one is using id numbers or type(s), the \code{id} or \code{type} argument is used, which takes a vector of the entries' names. The \code{gibbs_ad_use} function samples a use date between the production and depositional densities from the results of \code{\link[eratosthenes]{gibbs_ad}}, and in turn pools those densities for the production and deposition of the stipulated ids/type; the resulting \code{list} object does _not_ express marginalized densities of production and deposition in light of the estimation of a use date.
#' 
#' See \code{\link[eratosthenes]{gibbs_ad}} for information on consistent batch means and Monte Carlo standard error, which are used to determined convergence for the use date.
#'
#' @param marginalized A \code{list} object of class \code{marginals}, the output of \code{\link[eratosthenes]{gibbs_ad}}.
#' @param finds Either the \code{list} object of finds used as input to produce \code{marginals} or a \code{data.frame} of two columns, the first listing the context and the second the incidence of the type in that context.
#' @param id A vector of the \code{id} of one or more specific finds whose use date is to be estimated. The values of \code{id} must match those in the \code{list} of \code{finds}. If \code{type} is used, \code{id} is ignored.
#' @param type A vector of one or more types to estimate a use density for. Must contain a value if \code{id} is \code{NULL}.
#' @param type_name A customized label for the type (e.g., if one is selecting via \code{id} or has combined subtypes). If only \code{type} is used to select finds, the default will be that label Otherwise the default is simply "Type."
#' @param max_samples Maximum number of samples to run. Default is \code{10^5}.
#' @param size The number of samples to take on each iteration of the main Gibbs sampler. Default is \code{10^3}. 
#' @param mcse_crit Criterion for the Monte Carlo standard error to stop the Gibbs sampler. Only the MCSE of the use date is used as a stopping rule.
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
#' coin2 <- list(id = "coin2", assoc = "G", type = NULL, samples = seq(37, 41, length = 100))
#' destr <- list(id = "destr", assoc = "J", type = NULL, samples = 79)
#' 
#' tpq_info <- list(coin1, coin2)
#' taq_info <- list(destr)
#' 
#' result <- gibbs_ad(contexts, finds = artifacts, tpq = tpq_info, taq = taq_info)
#' 
#' # use dates by specifying ids
#' gibbs_ad_use(result, artifacts, id = c("find04", "find05"), max_samples = 5000, mcse_crit = 2)
#'
#' # use dates by specifying types
#' gibbs_ad_use(result, artifacts, type = "type1", max_samples = 5000, mcse_crit = 2)
#' 
#' @returns A \code{list} of class \code{use_marginals} of the density of a use date, conditional upon production and depositional dates.
#' 
#' @export
gibbs_ad_use <- function(marginalized, finds, id = NULL, type = NULL, type_name = NULL, max_samples = 10^5, size = 10^3, mcse_crit = 0.5) {
    UseMethod("gibbs_ad_use")
}
#' 
#' @rdname gibbs_ad_use
#' @export
gibbs_ad_use.marginals <- function(marginalized, finds, id = NULL, type = NULL, type_name = NULL, max_samples = 10^5, size = 10^3, mcse_crit = 0.5) {
    if (is.data.frame(finds)) {
        finds <- finds_d2l(finds)
    }
    if (!is.list(finds)) {
        stop("finds must be list or data frame object.")
    }
    if (is.null(id) & is.null(type)) {
        stop("Either one or more id or types must be specified.")
    } else {
        if (length(type) > 0) {
            if (length(id) > 0) {
                message('Both id and type specified. Defaulting to type (omit or specify "type = NULL" to use id).')
            }
            id <- ids_of_types(finds, type)
        }
    }
    # if (is.null(id)) {
    #     stop("id or type needed as input for finds argument.")
    # }

    message("Estimating use date for id(s)/type(s) specified, sampling in between dates of production and deposition.")
    sequence_ <- list()

    if (length(id) == 0) {
        stop("No ids or types found with that name.")
    }

    tpq_ <- list()
    taq_ <- list()

    for (i in 1:length(id)) {
        for (j in finds) {
            if (j$id == id[i]) {
                pooled <- c()
                k_type <- j$type
                for (k in 1:length(k_type)) {
                    pooled <- c(pooled, marginalized$production[[k_type[k]]])
                }
                tpq_[[i]] <- pooled
                taq_[[i]] <- marginalized$deposition[[j$assoc]]
            }
        }
    }

    gibbs0 <- gibbs_ad_use_init_cpp(tpq_, taq_, size)

    # for pooling all prd, dep, and use dates into respective distr
    idx1_ <- length(tpq_) + 1 
    idx2_ <- length(tpq_) * 2 + 1

    # consistent batch means
    mcse_check <- FALSE
    while (mcse_check == FALSE) {
        vec <- c(c(gibbs0[1:length(tpq_), ]), c(gibbs0[idx1_:(idx1_+(length(tpq_))-1), ]), c(gibbs0[idx2_:(idx2_+(length(tpq_))-1), ]) )
        gibbs <- matrix(vec, nrow = 3, byrow = TRUE)

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
        
        mcse <- sqrt( rowSums( (m_batch - rowMeans(m_batch))^2 ) * (n_batch / (K-1) ) ) / sqrt((K - 1) * n_batch)
        
        # do not include fixed single-point events as part of estimating MCSE
        mcse <- mcse[mcse > 0] 

        cat("\r", paste0("Samples: ", ncol(gibbs), "     Mean MCSE: ",  round(mean(mcse),3)))
        if (mean(mcse) < mcse_crit) {
            mcse_check <- TRUE
            message("\nMCSE criterion passed. Finishing.")
        } else {
            if (ncol(gibbs) >= max_samples) {
                message("\nMC samples exceeded maximum stipulated without passing MCSE crterion. Finishing.")
                mcse_check <- TRUE
            } else {
                gibbs0_next <- matrix(0, nrow = nrow(gibbs0), ncol = (size + 1) )
                gibbs0_next[,1] <- gibbs0[,ncol(gibbs0)]
                gibbs0_next <- gibbs_ad_use_cpp(gibbs0_next, tpq_, taq_)
                gibbs0 <- cbind(gibbs0, gibbs0_next[, 2:ncol(gibbs0_next)])
            }
        }
    }


    prd <- gibbs[1,]
    dep <- gibbs[2,]
    use <- gibbs[3,]
    
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

    marginalized[['use']] <- list(use_name = type_name, use_date = use, production_date = prd, deposition_date = dep, use_mcse = mcse[3], production_mcse = mcse[1], deposition_mcse = mcse[2])
    class(marginalized) <- c("use_marginals", "list")
    return(marginalized)
} 




#' @export
print.events <- function(x, ...) {
    x <- sort(as.vector(names(x)))
    cat("\n Events (call with $ operator to retrieve samples):\n")
    cat(" ", x, "\n\n\n")
}



#' @export 
print.marginals <- function(x, ...) {
    mcse <- x$mcse
    mcse <- mcse[c(names(x$deposition), names(x$externals))]
    mmcse <- mean(mcse[mcse > 0])
    cat("\n Marginals from joint conditional density, consisting of:\n",
    "    ", length(x$deposition), "depositional events\n", 
    "    ", length(x$externals), "external constraints (t.p./a.q.)\n",
    "    ", length(x$production), "production dates\n\n",
        "Call the following objects to see element names: \n",
        "     $deposition", "\n",
        "     $externals", "\n",
        "     $production\n\n",
        "Mean Monte Carlo standard error (MCSE) of depositional/external variates: \n",
    "    ", round(mmcse, 3), "\n\n\n")
}



#' @export 
print.use_marginals <- function(x, ...) {
    mcse <- x$mcse
    mcse <- mcse[c(names(x$deposition), names(x$externals))]
    mmcse <- mean(mcse[mcse > 0])
    cat("\n Marginals from joint conditional density, consisting of:\n",
    "    ", length(x$deposition), "depositional events\n", 
    "    ", length(x$externals), "external constraints (t.p./a.q.)\n",
    "    ", length(x$production), "production dates\n\n",
        "Call the following objects to see element names: \n",
        "     $deposition", "\n",
        "     $externals", "\n",
        "     $production\n\n",
        "Mean Monte Carlo standard error (MCSE) of depositional/external variates: \n",
    "    ", round(mmcse, 3), "\n\n",
        "Marginals of production, use, and depositional dates: \n",
    "    ", "Label:", x$use$use_name, "(", length(x$use$use_date), "use events pooled)\n",
    "    ", "Mean MCSE:", round(mean(c(x$use$production_mcse, x$use$use_mcse, x$use$deposition_mcse)), 3),"\n\n",
        "Call the following objects for the stipulated type's production, use, and deposition dates: \n",
    "     $use$production_date\n",
    "     $use$use_date\n",
    "     $use$deposition_date\n"
    )  
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
summary.use_marginals <- function(object, digits = 2, ...) {
    depmu <- mean(object$use$deposition_date)
    depmcse <- object$use$deposition_mcse
    depdat <- data.frame(Mean = depmu, MCSE = depmcse)
    rownames(depdat) <- "Deposition"

    usemu <- mean(object$use$use_date)
    usemcse <- object$use$use_mcse
    usedat <- data.frame(Mean = usemu, MCSE = usemcse)
    rownames(usedat) <- "Use"

    prdmu <- mean(object$use$production_date)
    prdmcse <- object$use$production_mcse
    prddat <- data.frame(Mean = prdmu, MCSE = prdmcse)
    rownames(prddat) <- "Production"

    dat <- rbind(prddat, usedat, depdat)
    dat <- round(dat, digits)

    message("\nMonte Carlo mean and standard errors of marginalized production, use, and depositional dates: \n")

    return(dat)
}



#' Traceplot of Gibbs Samples
#' 
#' Wrapper around \code{\link[graphics]{plot}} to make a traceplot of Gibbs samples from \code{\link[eratosthenes]{gibbs_ad}}. See \code{\link[eratosthenes]{hist.marginals}} for plotting a density histogram of events.
#' 
#' Also see \code{\link[eratosthenes]{tidy_marginals}} for exporting the results of these functions into tidy data frame for custom plotting in e.g., \code{ggplot2}.
#' 
#' @param x A \code{list} object of class \code{marginals} or \code{use_marginals}, the output of \code{\link[eratosthenes]{gibbs_ad}} or \code{\link[eratosthenes]{gibbs_ad_use}} respectively.
#' @param events A vector or element of the event names to plot. Maximum number of events is 12.
#' @param xlim The limits of the x-axis (optional).
#' @param ylim The limits of the y-axis (optional).
#' @param xlab Label for the y-axis. Default is \code{"Index"}.
#' @param ylab Label for the y-axis. Default is \code{"Year"}.
#' @param palette A vector providing the color palette of the histogram. The default is \code{"colorBlindness::paletteMartin"} (see \code{\link[paletteer]{palettes_d}}).
#' @param opacity The opacity/transparency of the traceplot, if visualizing overlapping events. A value between 0 and 1 (default).
#' @param legend_pos The position of the legend in the plot. Default is \code{"topright"}.
#' @param ... Additional graphical parameters passed to \code{\link[graphics]{plot}}.
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
#' coin2 <- list(id = "coin2", assoc = "G", type = NULL, samples = seq(37, 41, length = 100))
#' destr <- list(id = "destr", assoc = "J", type = NULL, samples = 79)
#' 
#' tpq_info <- list(coin1, coin2)
#' taq_info <- list(destr)
#' 
#' result <- gibbs_ad(contexts, finds = artifacts, tpq = tpq_info, taq = taq_info)
#' 
#' traceplot(result, "B")
#' traceplot(result, c("coin1", "B", "H"), opacity = 0.5)
#' 
#' @returns A traceplot of the Gibbs samples of the selected events.
#' 
#' @export
traceplot <- function(x, events = NULL, xlim = NULL, ylim = NULL, xlab = "Index", ylab = "Year", palette = NULL, opacity = 1, legend_pos = "topright", ...) {
    UseMethod("traceplot")
}
#' 
#' @rdname traceplot
#' @export
traceplot.marginals <- function(x, events = NULL, xlim = NULL, ylim = NULL, xlab = "Index", ylab = "Year", palette = NULL, opacity = 1, legend_pos = "topright", ...) {
    if (is.vector(events)) {
        master <- tidy_marginals(x)
    } else {
        stop('Events must be a single element or vector.')
    }

    if (is.null(palette)) {
        palette <- paletteer::paletteer_d("colorBlindness::paletteMartin")
    }

    palette <- grDevices::adjustcolor(palette, alpha.f = opacity)

    if (is.null(xlim)) {
        xlim <- c(0, length(master[master$event == levels(master$event)[1], ]$year))
    }
    if (is.null(ylim)) {
        y_all <- master[master$event %in% events, ]$year
        ylim <- c(min(y_all), max(y_all))
    }

    if (length(events) == 1) {
        x <- master[master$event == events[1], ]$year
        graphics::plot(x, type = "l", col = palette[1], xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, main = "", ...)
    } else if (length(events) > 1 & length(events <= 12)) {
        x <- master[master$event == events[1], ]$year
        graphics::plot(x, type = "l", col = palette[1], xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, main = "", ...)
        for (k in 2:length(events)) {
            x <- master[master$event == events[k], ]$year
            graphics::lines(x, type = "l", col = palette[k], ...)
        }
    } else {
        stop("Max number events for histogram is 12.")
    }
    graphics::legend(legend_pos, legend= events, col = palette[1:length(events)], pt.cex = 2.5, pch = 15)
}



#' Histogram of Marginal Densities
#' 
#' Wrapper around \code{\link[graphics]{hist}} to plot density histograms for select marginal densities (up to 12) in a single plot, from the results of \code{\link[eratosthenes]{gibbs_ad}}. See \code{\link[eratosthenes]{hist.use_marginals}} for plotting a single histogram of an artifact type production, use, and deposition.
#' 
#' Also see also \code{\link[eratosthenes]{tidy_marginals}} for exporting the results of these functions into tidy data frame for custom plotting in e.g., \code{ggplot2}.
#' 
#' @param x A \code{list} object of class \code{marginals}, the output of \code{\link[eratosthenes]{gibbs_ad}}.
#' @param events A vector or element of the event names to plot. Maximum number of events is 12.
#' @param breaks The number or method of breaks in the histogram. Default is \code{"Freedman-Diaconis"}. See \code{\link[graphics]{hist}} for more.
#' @param xlim The limits of the x-axis. Default is set to the min/max values of all samples.
#' @param ylim The limits of the y-axis. This may need to be adjusted if densities have an extremely narrow interval.
#' @param xlab Label for the x-axis. Default is \code{"Year"}.
#' @param palette A vector providing the color palette of the histogram. The default is \code{"colorBlindness::paletteMartin"} (see \code{\link[paletteer]{palettes_d}}).
#' @param opacity The opacity/transparency of the histograms for visualizing overlapping events, a value between 0 and 1 (default).
#' @param legend_pos The position of the legend in the plot. Default is \code{"topright"}.
#' @returns A density histogram of the selected events.
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
#' coin2 <- list(id = "coin2", assoc = "G", type = NULL, samples = seq(37, 41, length = 100))
#' destr <- list(id = "destr", assoc = "J", type = NULL, samples = 79)
#' 
#' tpq_info <- list(coin1, coin2)
#' taq_info <- list(destr)
#' 
#' result <- gibbs_ad(contexts, finds = artifacts, tpq = tpq_info, taq = taq_info)
#' 
#' # deposition of "B"
#' hist(result, "B")
#' 
#' # deposition of "coin2" and deposition of "G"
#' hist(result, c("coin2", "G"), opacity = 0.5)
#' 
#' # deposition of "B", "D", and "H", adjusting y-axis
#' hist(result, c("B", "D", "H"), ylim = c(0,0.05), opacity = 0.5)
#' 
#' # production of "type2" and deposition of "H"
#' hist(result, c("H", "type2"), opacity = 0.5)
#' 
#' @returns A density histogram of the selected events.
#' 
#' @export
hist.marginals <- function(x, events = NULL, breaks = "Freedman-Diaconis", xlim = NULL, ylim = NULL, xlab = "Year", palette = NULL, opacity = 1, legend_pos = "topright") {
    if (is.vector(events)) {
        master <- tidy_marginals(x)
    } else {
        stop('Events must be a single element or vector.')
    }

    if (is.null(palette)) {
        palette <- paletteer::paletteer_d("colorBlindness::paletteMartin")
    }

    palette <- grDevices::adjustcolor(palette, alpha.f = opacity)

    x_all <- master[master$event %in% events, ]$year
    xlims <- c(min(x_all), max(x_all))

    if (length(events) == 1) {
        x <- master[master$event == events[1], ]$year
        graphics::hist(x, breaks = breaks, freq = FALSE, xlim = xlims, xlab = xlab, ylim = ylim, col = palette[1], lty="blank", main = "")
    } else if (length(events) > 1 & length(events <= 12)) {
        x <- master[master$event == events[1], ]$year
        graphics::hist(x, breaks = breaks, freq = FALSE, lty="blank", xlim = xlims, ylim = ylim, col = palette[1], xlab = xlab, main = "")
        for (k in 2:length(events)) {
            x <- master[master$event == events[k], ]$year
            graphics::hist(x, breaks = breaks, freq = FALSE, lty="blank", col = palette[k], xlab = xlab, main = "", add = TRUE)
        }
    } else {
        stop("Max number events for histogram is 12.")
    }
    graphics::legend(legend_pos, legend = events, col = palette[1:length(events)], pt.cex = 2.5, pch = 15)
}



#' Histogram of Marginal Densities (Artifact Type Production, Use, and Deposition)
#' 
#' Wrapper around \code{\link[graphics]{hist}} to plot density histograms of the production, use, and deposition events of an artifact type, from the results of \code{\link[eratosthenes]{gibbs_ad_use}]}, in a single plot. As with \code{\link[eratosthenes]{gibbs_ad_use}]}, production dates and depositional dates are the pooled densities of each individual artifact appertaining to that type.
#' 
#' See \code{\link[eratosthenes]{hist.marginals}} for plotting multiple production, deposition, and absolute constraints in a single histogram. 
#' 
#' See \code{\link[eratosthenes]{tidy_marginals}} for exporting the results of these functions into tidy data frame for custom plotting in e.g., \code{ggplot2}.
#' 
#' @param x A \code{list} object of class \code{use_marginals}, the output of \code{\link[eratosthenes]{gibbs_ad}}.
#' @param aspect A vector of one or more of \code{c("production", "use", "deposition")}. The default is all three.
#' @param display_name The name of the artifact type to display in the histogram legend. Default is \code{"Type"}.
#' @param breaks The number or method of breaks in the histogram. Default is \code{"Freedman-Diaconis"}. See \code{\link[graphics]{hist}} for more.
#' @param xlim The limits of the x-axis. Default is set to the min/max values of all samples.
#' @param ylim The limits of the y-axis. This may need to be adjusted if densities have an extremely narrow interval.
#' @param xlab Label for the x-axis. Default is \code{"Year"}.
#' @param palette A vector providing the color palette of the histogram. The default is \code{"colorBlindness::paletteMartin"} (see \code{\link[paletteer]{palettes_d}}).
#' @param opacity The opacity/transparency of the histograms for visualizing overlapping events, a value between 0 and 1. Default is set at 0.5.
#' @param legend_pos The position of the legend in the plot. Default is \code{"topright"}.
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
#' coin2 <- list(id = "coin2", assoc = "G", type = NULL, samples = seq(37, 41, length = 100))
#' destr <- list(id = "destr", assoc = "J", type = NULL, samples = 79)
#' 
#' tpq_info <- list(coin1, coin2)
#' taq_info <- list(destr)
#' 
#' result <- gibbs_ad(contexts, finds = artifacts, tpq = tpq_info, taq = taq_info)
#' 
#' hist(result, "B")
#' hist(result, c("B", "D", "H"), ylim = c(0,0.05), opacity = 0.5)
#' 
#' @returns A density histogram of the selected aspects of production, use, and deposition.
#' 
#' @export
hist.use_marginals <- function(x, aspect = c("production", "use", "deposition"), display_name = "Type", breaks = "Freedman-Diaconis", xlim = NULL, ylim = NULL, xlab = "Year", palette = NULL, opacity = 0.5, legend_pos = "topright") {

    if ((!("production" %in% aspect) | !("use" %in% aspect) ) | !("deposition" %in% aspect))  {
        stop("At least one or more event types (production, use, deposition) must be selected.")
    }

    if (is.null(palette)) {
        palette <- paletteer::paletteer_d("colorBlindness::paletteMartin")
    }

    dat <- c(year = c(), aspect = c())

    dep_samples <- x$use$deposition_date
    prd_samples <- x$use$production_date
    use_samples <- x$use$use_date


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

    xlims <- c(min(dat$year), max(dat$year))

    if (length(aspect) == 1) {
        x <- dat$year
        graphics::hist(x, breaks = breaks, freq = FALSE, xlim = xlims, ylim = ylim, col = palette[1], xlab = xlab, lty="blank", main = "")
    } else {
        x <- dat$year
        graphics::hist(x, breaks = breaks, freq = FALSE, lty="blank", xlim = xlims, ylim = ylim, col = palette[1], xlab = xlab, main = "")
        for (k in 2:length(aspect_name)) {
            x <- dat[dat$aspect == aspect_name[k], ]$year
            graphics::hist(x, breaks = breaks, freq = FALSE, lty="blank", col = palette[k], xlab = xlab, main = "", add = TRUE)
        }
    }
    graphics::legend(legend_pos, legend = aspect_name, col = palette[1:length(aspect)], pt.cex = 2.5, pch = 15)
}



#' Convert Marginals to Tidy (Molten) Data Frame
#' 
#' Takes the results of \code{\link[eratosthenes]{gibbs_ad}} or \code{\link[eratosthenes]{gibbs_ad_use}} and "melts" the \code{list} into a tidy data frame \insertCite{wickham_tidy_2014}{eratosthenes}. Each row of the molten data frame will contain the index of the Monte Carlo sample, the sample itself, and then the event name.
#' 
#' @param input An object of class \code{marginals} or \code{use_marginals}, the output of \code{\link[eratosthenes]{gibbs_ad}} or \code{\link[eratosthenes]{gibbs_ad_use}}.
#' @returns A data frame giving the MC sampling index (\code{idx}), the sample (\code{year}), and the event (\code{event}).
#' 
#' @examples
#' x <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J")
#' y <- c("B", "D", "G", "H", "K")
#' z <- c("F", "K", "L", "M")
#' 
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
#' coin2 <- list(id = "coin2", assoc = "G", type = NULL, samples = seq(37, 41, length = 100))
#' destr <- list(id = "destr", assoc = "J", type = NULL, samples = 79)
#' 
#' tpq_info <- list(coin1, coin2)
#' taq_info <- list(destr)
#' 
#' result <- gibbs_ad(contexts, finds = artifacts, tpq = tpq_info, taq = taq_info)
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
tidy_marginals.use_marginals <- function(input) {
    dat <- data.frame(x = c(), event = c())

    dep_samples <- input$use$deposition_date
    prd_samples <- input$use$production_date

    for (i in 1:length(input$use$use_date)) {
        use_samples <- input$use$use_date[[i]]
    } 

    dat <- rbind(dat, data.frame(x = use_samples, event = paste0(input$use$use_name, "- Use")))
    dat <- rbind(dat, data.frame(x = dep_samples, event = paste0(input$use$use_name, "- Deposition")))
    dat <- rbind(dat, data.frame(x = prd_samples, event = paste0(input$use$use_name, "- Production")))

    return(dat)
}



#' Convert Finds List Object to Data Frame (Context / Find-Type)
#' 
#' Performs the opposite of \code{\link[eratosthenes]{finds_d2l}}. Takes a \code{list} object of finds and their types, used as input in \code{\link[eratosthenes]{gibbs_ad}}, and returns a \code{data.frame} of two columns, containing the context in the first and the find-type in the second, and the \code{id} of the object in the third.
#' 
#' @param input A list object of finds (each one a list) of associated contexts and types. 
#' @returns A three-column data frame of contexts (first column) and find-types attested in that context (second column), along with the id number (third column).
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
#' Performs the opposite of \code{\link[eratosthenes]{finds_l2d}}. Takes a \code{data.frame} object of two columns, containing the context in the first and the find-type in the second, and returns a \code{list} object for input in \code{\link[eratosthenes]{gibbs_ad}}. The value of the find \code{id} is automatically generated as an integer if not provided in a third column.
#' 
#' @param input A two-column data frame of contexts (first column) and find-types (second column). An optional third column of an id number may be provided.
#' @returns A list of finds (each one a list) associated with contexts and their types.
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
            


#' Mean Squared Displacement of Events
#' 
#' Computes the mean squared displacement (MSD) of all events contained in the relative sequences and absolute constraints used in the execution of \code{\link[eratosthenes]{gibbs_ad}}. 
#' 
#' The MSD entails the following jackknife/leave-one-out style routine:
#' 
#' * Each event is omitted from all relative and absolute sequences, and the function \code{\link[eratosthenes]{gibbs_ad}} is re-run to compute a "jackknifed" Monte Carlo mean for that event.
#'   * The squared difference of this jackknifed Monte Carlo mean and the original is then computed as its squared "displacement" in time.
#'   * The mean of the squared displacements of all events is then computed and attributed to the omitted event.
#' 
#' If an event has a low MSD, it bears a low impact on the rest of the events within the full joint conditional density. If it is has a high MSD, other events depend heavily upon its inclusion in the full joint density.
#' 
#' Trimming is not implemented in the computation of MSD, and so attention should be paid to the selection of \code{alpha_} and \code{omega_}, and reported. This is owing to the way in which, if an absolute constraint (\code{tpq} or \code{taq}) is omitted that happens to be an earliest or latest bounding event, there still needs to be earliest and latest thresholds in place. 
#' 
#' This function is fairly computationally intensive and thus a lower value of `max_samples` and a higher value of `mcse_crit` may be warranted
#'  
#' @param marginalized The results of \code{\link[eratosthenes]{gibbs_ad}}.\code{\link[eratosthenes]{gibbs_ad}}.
#' @param sequences A \code{list} of relative sequences of elements (e.g., contexts) used to compute \code{marginalized}.
#' @param finds Optional. A \code{list} of finds related to (contained in) the elements of \code{sequences}.
#' @param max_samples Maximum number of samples to run. Default is \code{10^5}.
#' @param size The number of samples to take on each iteration of the main Gibbs sampler. Default is \code{10^3}. 
#' @param mcse_crit Criterion for the Monte Carlo standard error to stop the Gibbs sampler. A higher MCSE is recommended for situations with a higher number of events in order to reduce computational time.
#' @param tpq A \code{list} containing \emph{termini post quos} used to compute \code{marginalized}. See \code{\link[eratosthenes]{gibbs_ad}} for details.
#' @param taq A \code{list} containing \emph{termini ante quos} used to compute \code{marginalized}. See \code{\link[eratosthenes]{gibbs_ad}} for details.
#' @param alpha_ An initial \emph{t.p.q.} to limit any elements which may occur before the first provided \emph{t.p.q.} Default is \code{-5000}.
#' @param omega_ A final \emph{t.a.q.} to limit any elements which may occur after the after the last provided \emph{t.a.q.} Default is \code{1950}.
#' @param rule The rule for computing an estimated date of production. See \code{\link[eratosthenes]{gibbs_ad}} for details.
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
#' coin2 <- list(id = "coin2", assoc = "G", type = NULL, samples = seq(37, 41, length = 100))
#' destr <- list(id = "destr", assoc = "J", type = NULL, samples = 79)
#' 
#' tpq_info <- list(coin1, coin2)
#' taq_info <- list(destr)
#' 
#' result <- gibbs_ad(contexts, finds = artifacts, tpq = tpq_info, taq = taq_info)
#' 
#' result_msd <- msd(result, contexts, finds = artifacts, max_samples = 5000,
#'                   mcse_crit = 2, tpq = tpq_info, taq = taq_info)
#'
#' @returns Output is a list containing a data frame \code{MSD_stats} giving the mean MC date, the MCSE, the MSD, the variance of the squared displacements (not the standard error), and sample size, as well as a vector \code{bounds} of the values of \code{alpha_} and \code{omega_}.
#' 
#' @export
msd <- function(marginalized, sequences, finds = NULL, max_samples = 10^5, size = 10^3, mcse_crit = 0.5, tpq = NULL, taq = NULL, alpha_ = -5000, omega_ = 1950, rule = "naive") {
    UseMethod("msd")
}
#' 
#' @rdname msd
#' @export
msd.marginals <- function(marginalized, sequences, finds = NULL, max_samples = 10^5, size = 10^3, mcse_crit = 0.5, tpq = NULL, taq = NULL, alpha_ = -5000, omega_ = 1950, rule = "naive") {
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

    proceed_all <- rownames(orig_dat)

    orig_dat$MSD <- NA
    orig_dat$MSD_var <- NA
    orig_dat$MSD_n <- NA

    message("Beginning jackknlife/LOO-style routine to compute MSD. This may take a while, depending on the number of events / mcse_crit...\n")

    for (j in 1:length(proceed_all)) {
        cat("Depositional Event / Absolute Constraint: ", proceed_all[j], "\n")
        MCmean <- orig_dat[proceed_all[j],]$Mean

        sequencesMSD <- list()
        tpqMSD <- list()
        taqMSD <- list()
        for (i in 1:length(sequences)) {
            sq_ <- sequences[[i]]
            sq_ <-sq_[!(sq_ %in% proceed_all[j])]
            sequencesMSD[[i]] <- sq_
        }
        idx <- 1
        for (i in 1:length(tpq)) {
            if (!(tpq[[i]]$id %in% proceed_all[j] | tpq[[i]]$assoc %in% proceed_all[j])) {
                tpqMSD[[idx]] <- tpq[[i]]
                idx <- idx + 1
            }
        }
        idx <- 1
        for (i in 1:length(taq)) {
            if (!(taq[[i]]$id %in% proceed_all[j] | taq[[i]]$assoc %in% proceed_all[j])) {
                taqMSD[[idx]] <- taq[[i]]
                idx <- idx + 1
            }
        }
        if (!is.null(finds)) {
            findsMSD <- list()
            idx <- 1
            for (i in 1:length(finds)) {
                if (!(finds[[i]]$assoc %in% proceed_all[j])) {
                    findsMSD[[idx]] <- finds[[i]]
                    idx <- idx + 1
                }
            }
        }

        if (!is.null(finds)) {
            if (length(tpqMSD) == 0 & length(taqMSD) != 0) {
                gibbsLOO <- gibbs_ad(sequencesMSD, finds = findsMSD, max_samples, size, mcse_crit, tpq = NULL, taq = taqMSD, alpha_, omega_, trim = FALSE, rule)
            } else if (length(tpqMSD) != 0 & length(taqMSD) == 0) {
                gibbsLOO <- gibbs_ad(sequencesMSD, finds = findsMSD, max_samples, size, mcse_crit, tpq = tpqMSD, taq = NULL, alpha_, omega_, trim = FALSE, rule)
            } else if (length(tpqMSD) == 0 & length(taqMSD) == 0) {
                gibbsLOO <- gibbs_ad(sequencesMSD, finds = findsMSD, max_samples, size, mcse_crit, tpq = NULL, taq = NULL, alpha_, omega_, trim = FALSE, rule)
            } else {
                gibbsLOO <- gibbs_ad(sequencesMSD, finds = findsMSD, max_samples, size, mcse_crit, tpq = tpqMSD, taq = taqMSD, alpha_, omega_, trim = FALSE, rule)
            }
        } else {
            if (length(tpqMSD) == 0 & length(taqMSD) != 0) {
                gibbsLOO <- gibbs_ad(sequencesMSD, NULL, max_samples, size, mcse_crit, tpq = NULL, taq = taqMSD, alpha_, omega_, trim = FALSE, rule)
            } else if (length(tpqMSD) != 0 & length(taqMSD) == 0) {
                gibbsLOO <- gibbs_ad(sequencesMSD, NULL, max_samples, size, mcse_crit, tpq = tpqMSD, taq = NULL, alpha_, omega_, trim = FALSE, rule)
            } else if (length(tpqMSD) == 0 & length(taqMSD) == 0) {
                gibbsLOO <- gibbs_ad(sequencesMSD, NULL, max_samples, size, mcse_crit, tpq = NULL, taq = NULL, alpha_, omega_, trim = FALSE, rule)
            } else {
                gibbsLOO <- gibbs_ad(sequencesMSD, NULL, max_samples, size, mcse_crit, tpq = tpqMSD, taq = taqMSD, alpha_, omega_, trim = FALSE, rule)
            }
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
        } else {
            prdmu <- sapply(gibbsLOO$production, mean)
            prdmcse <- gibbsLOO$mcse[names(gibbsLOO$production)]
            prddat <- data.frame(Mean = prdmu, MCSE = prdmcse)
            rownames(prddat) <- names(gibbsLOO$production)
            LOO_dat <- rbind(depdat, extdat, prddat)
        }

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

    result <- list(MSD_stats = orig_dat, bounds = bounds_)
    class(result) <- c("msd_data", "list")
    return(result)
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
#' Computes the squared displacement for a target event within the joint conditional density, estimating how much the omission of another event will change the date of that event. See also \code{\link[eratosthenes]{msd}}. 
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
#' @param marginalized The results of \code{\link[eratosthenes]{gibbs_ad}}.\code{\link[eratosthenes]{gibbs_ad}}.
#' @param target The target event (any event in \code{marginalized}) for which to estimate squared displacement.
#' @param sequences A \code{list} of relative sequences of elements (e.g., contexts) used to compute \code{marginalized}.
#' @param finds Optional. A \code{list} of finds related to (contained in) the elements of \code{sequences}.
#' @param max_samples Maximum number of samples to run. Default is \code{10^5}.
#' @param size The number of samples to take on each iteration of the main Gibbs sampler. Default is \code{10^3}. 
#' @param mcse_crit Criterion for the Monte Carlo standard error to stop the Gibbs sampler. A higher MCSE is recommended for situations with a higher number of events in order to reduce computational time.
#' @param tpq A \code{list} containing \emph{termini post quos} used to compute \code{marginalized}. See \code{\link[eratosthenes]{gibbs_ad}} for details.
#' @param taq A \code{list} containing \emph{termini ante quos} used to compute \code{marginalized}. See \code{\link[eratosthenes]{gibbs_ad}} for details.
#' @param alpha_ An initial \emph{t.p.q.} to limit any elements which may occur before the first provided \emph{t.p.q.} Default is \code{-5000}.
#' @param omega_ A final \emph{t.a.q.} to limit any elements which may occur after the after the last provided \emph{t.a.q.} Default is \code{1950}.
#' @param rule The rule for computing an estimated date of production. See \code{\link[eratosthenes]{gibbs_ad}} for details.
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
#' coin2 <- list(id = "coin2", assoc = "G", type = NULL, samples = seq(37, 41, length = 100))
#' destr <- list(id = "destr", assoc = "J", type = NULL, samples = 79)
#' 
#' tpq_info <- list(coin1, coin2)
#' taq_info <- list(destr)
#' 
#' result <- gibbs_ad(contexts, finds = artifacts, tpq = tpq_info, taq = taq_info)
#' 
#' # max_samples lowered and msce_crit raised for examples
#' 
#' # squared displacement for depositional context "E"
#' E_sqdisp <- sq_disp(result, target = "E", sequences = contexts, 
#'                     max_samples = 5000, mcse_crit = 2, tpq = tpq_info, taq = taq_info)
#'
#' # squared displacement for production of artifact type "type1"
#' type1_sqdisp <- sq_disp(result, target = "type1", sequences = contexts, finds = artifacts,
#'                         max_samples = 5000, mcse_crit = 2, tpq = tpq_info, taq = taq_info)
#'
#' @returns Output is a list containing a data frame \code{sq_disp} giving the diplacement with respect to all other events and a vector \code{bounds} of the values of \code{alpha_} and \code{omega_}.
#' 
#' @export
sq_disp <- function(marginalized, target, sequences, finds = NULL, max_samples = 10^5, size = 10^3, mcse_crit = 0.5, tpq = NULL, taq = NULL, alpha_ = -5000, omega_ = 1950, rule = "naive") {
    UseMethod("sq_disp")
}
#' 
#' @rdname sq_disp
#' @export
sq_disp.marginals <- function(marginalized, target, sequences, finds = NULL, max_samples = 10^5, size = 10^3, mcse_crit = 0.5, tpq = NULL, taq = NULL, alpha_ = -5000, omega_ = 1950, rule = "naive") {
    if (is.null(target)) {
        stop("Target event must be specified.")
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

    proceed_all <- rownames(orig_dat)
    proceed_all <- proceed_all[!(proceed_all %in% target)]

    orig_dat$sq_disp <- NA

    message("Beginning jackknlife/LOO-style routine to compute squared displacement. This may take a while, depending on the number of events / mcse_crit...\n")

    for (j in 1:length(proceed_all)) {
        cat("Depositional Event / Absolute Constraint: ", proceed_all[j], "\n")
        MCmean <- orig_dat[proceed_all[j],]$Mean

        sequencesMSD <- list()
        tpqMSD <- list()
        taqMSD <- list()
        for (i in 1:length(sequences)) {
            sq_ <- sequences[[i]]
            sq_ <-sq_[!(sq_ %in% proceed_all[j])]
            sequencesMSD[[i]] <- sq_
        }
        idx <- 1
        for (i in 1:length(tpq)) {
            if (!(tpq[[i]]$id %in% proceed_all[j] | tpq[[i]]$assoc %in% proceed_all[j])) {
                tpqMSD[[idx]] <- tpq[[i]]
                idx <- idx + 1
            }
        }
        idx <- 1
        for (i in 1:length(taq)) {
            if (!(taq[[i]]$id %in% proceed_all[j] | taq[[i]]$assoc %in% proceed_all[j])) {
                taqMSD[[idx]] <- taq[[i]]
                idx <- idx + 1
            }
        }
        if (!is.null(finds)) {
            findsMSD <- list()
            idx <- 1
            for (i in 1:length(finds)) {
                if (!(finds[[i]]$assoc %in% proceed_all[j])) {
                    findsMSD[[idx]] <- finds[[i]]
                    idx <- idx + 1
                }
            }
        }

        if (!is.null(finds)) {
            if (length(tpqMSD) == 0 & length(taqMSD) != 0) {
                gibbsLOO <- gibbs_ad(sequencesMSD, finds = findsMSD, max_samples, size, mcse_crit, tpq = NULL, taq = taqMSD, alpha_, omega_, trim = FALSE, rule)
            } else if (length(tpqMSD) != 0 & length(taqMSD) == 0) {
                gibbsLOO <- gibbs_ad(sequencesMSD, finds = findsMSD, max_samples, size, mcse_crit, tpq = tpqMSD, taq = NULL, alpha_, omega_, trim = FALSE, rule)
            } else if (length(tpqMSD) == 0 & length(taqMSD) == 0) {
                gibbsLOO <- gibbs_ad(sequencesMSD, finds = findsMSD, max_samples, size, mcse_crit, tpq = NULL, taq = NULL, alpha_, omega_, trim = FALSE, rule)
            } else {
                gibbsLOO <- gibbs_ad(sequencesMSD, finds = findsMSD, max_samples, size, mcse_crit, tpq = tpqMSD, taq = taqMSD, alpha_, omega_, trim = FALSE, rule)
            }
        } else {
            if (length(tpqMSD) == 0 & length(taqMSD) != 0) {
                gibbsLOO <- gibbs_ad(sequencesMSD, NULL, max_samples, size, mcse_crit, tpq = NULL, taq = taqMSD, alpha_, omega_, trim = FALSE, rule)
            } else if (length(tpqMSD) != 0 & length(taqMSD) == 0) {
                gibbsLOO <- gibbs_ad(sequencesMSD, NULL, max_samples, size, mcse_crit, tpq = tpqMSD, taq = NULL, alpha_, omega_, trim = FALSE, rule)
            } else if (length(tpqMSD) == 0 & length(taqMSD) == 0) {
                gibbsLOO <- gibbs_ad(sequencesMSD, NULL, max_samples, size, mcse_crit, tpq = NULL, taq = NULL, alpha_, omega_, trim = FALSE, rule)
            } else {
                gibbsLOO <- gibbs_ad(sequencesMSD, NULL, max_samples, size, mcse_crit, tpq = tpqMSD, taq = taqMSD, alpha_, omega_, trim = FALSE, rule)
            }
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

    result <- list(sq_disp = orig_dat, bounds = bounds_, target = target)
    class(result) <- c("sq_displ_data", "list")
    return(result)
}



#' @export
print.sq_displ_data <- function(x, ...) {  
    cat("\n For fixed bounds: (", x$bounds[1], ",",x$bounds[2], ")\n",
    "Squared displacement for target event", x$target, "by the following:  \n")
    print.data.frame(x$sq_disp)
    cat("\n")
}






