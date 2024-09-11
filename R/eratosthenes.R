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
#' @param ... `ties` The way in which ties are handled per the \code{\link{rank()}} function. The default is \code{"ties = average"}.
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
synth_rank <- function(obj, ...) {
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
#' For a \code{list} of partial sequences (of \code{vector} objects), generate another \code{list} which contains all elements that occur after it ("\emph{quae postea}"), i.e., analogous to a recursive trace through all partial sequences. A final element \code{"omega"} is added to all sets to avoid empty vectors. See also \code{\link[eratosthenes]{quae_antea}}.
#'
#' @param obj A \code{list} of \code{vector} objects which reperesent a sequence.    
#' 
#' @export
quae_postea <- function(obj) {
    UseMethod("quae_postea")
}
#' 
#' @rdname quae_postea
#' @export
quae_postea.list <- function(obj) {
    elements <- c()
    for (i in 1:length(obj)) {
        elements <- c(elements, obj[[i]])
    }
    elements <- unique(elements)
    
    qp <- list()
    for (i in elements) {
        qp[i] <- c("omega")
    }
    for (i in 1:length(obj)) {
        seq_ <- obj[[i]]
        for (j in 1:length(seq_)) {
            idx <- names(qp)[which(names(qp) == seq_[j])]
            qp[[idx]] <- c(qp[[idx]], seq_[j+1:length(seq_)] )
        }
    }
    for (i in names(qp)) {
        qp[[i]] <- unique( qp[[i]][ !is.na(qp[[i]]) ] )
    }

    qp2 <- list()
    for (i in names(qp)) {
        tmp <- qp[[i]]
        checked <- c()
        while (!setequal(tmp, checked)) {
            for (j in tmp) {
                if (!(j %in% checked)) {
                    tmp <- c(tmp, qp[[j]])
                    checked <- c(checked,j)
                }
            }            
        }
    qp2[[i]] <- unique(tmp)
    }
    return(qp2)
}




#' Quae Antea
#'
#' For a \code{list} of partial sequences (of \code{vector} objects), generate another \code{list} which contains all elements that occur before it ("\emph{quae antea}"), i.e., analogous to a recursive trace through all partial sequences. An element \code{"alpha"} is added to all sets to avoid empty vectors. See also \code{\link[eratosthenes]{quae_postea}}.
#'
#' @param obj A \code{list} of \code{vector} objects which reperesent a sequence.    
#'
#' @export
quae_antea <- function(obj) {
    UseMethod("quae_antea")
}
#' 
#' @rdname quae_antea
#' @export
quae_antea <- function(obj) {
    qa2 <- NULL

    elements <- c()
    for (i in 1:length(obj)) {
        elements <- c(elements, obj[[i]])
    }
    elements <- unique(elements)
    
    qa <- list()
    for (i in elements) {
        qa[i] <- c("alpha")
    }
    for (i in 1:length(obj)) {
        seq_ <- obj[[i]]
        for (j in 1:length(seq_)) {
            idx <- names(qa)[which(names(qa) == seq_[j])]
            qa[[idx]] <- c(qa[[idx]], seq_[1:j-1] )
        }
    }
    for (i in names(qa)) {
        qa[[i]] <- unique( qa[[i]][ !is.na(qa[[i]]) ] )
    }

    qa2 <- list()
    for (i in names(qa)) {
        tmp <- qa[[i]]
        checked <- c()
        while (!setequal(tmp, checked)) {
            for (j in tmp) {
                if (!(j %in% checked)) {
                    tmp <- c(tmp, qa[[j]])
                    checked <- c(checked,j)
                }
            }            
        }
    qa2[[i]] <- unique(tmp)
    }
    return(qa2)
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
#' A Gibbs sampler for archaeological dating, to fit relative sequences to absolute, calendrical dates. Elements can be associated with \emph{termini post quem} (\emph{t.p.q.}) and \emph{termini ante quem} (\emph{t.a.q.}), which are treated as a given probability density function \eqn{f(t)}. This function may take any form, a single date (i.e., with a probability of 1), a continuous uniform distribution (any time between two dates), or a bespoke density (as with calibrated radicarbon dates). Inputs of this function take samples drawn from their respective density functions. 
#'
#' @param sequences A \code{list} of relative sequences of elements (e.g., contexts).
#' @param finds Optional. A \code{list} of finds related to (contained in) the elements of `sequences`. If one includes this ob
#' @param samples Number of samples. Default is \code{10^5}.
#' @param tpq A \code{list} containing \emph{termini post quem}. Each object in the list consists of:
#'   * \code{id} A \code{character} ID of the  \emph{t.p.q.}, such as a reference or number.
#'   * \code{assoc} The element in \code{code} to which the \emph{t.p.q.} is associated. 
#'   * \code{samples} A vector of samples drawn from the appertaining probability density function of that \emph{t.p.q.}
#' @param taq A \code{list} containing \emph{termini ante quem}. Each object in the list consists of:
#'   * \code{id} A \code{character} ID of the  \emph{t.a.q.}, such as a reference or number.
#'   * \code{assoc} The element in \code{code} to which the \emph{t.p.q.} is associated. 
#'   * \code{samples} A vector of samples drawn from the appertaining probability density function of that \emph{t.p.q.}
#' @param alpha An initial \emph{t.p.q.} to limit any elements which may occur before the first provided \emph{t.p.q.} Default is \code{-5000}.
#' @param omega A final \emph{t.a.q.} to limit any elements which may occur after the after the last provided \emph{t.a.q.} Default is \code{1950}.
#' @param trim A logical value to determine whether elements that occur before the first \emph{t.p.q.} and after the last \emph{t.a.q.} should be ommitted from the results (i.e., to "trim" elements at the ends of the sequence, whose marginal densities depend on the selection of \code{alpha} and \code{omega}). Default is \code{TRUE}.
#' 
#' @returns A \code{list} object of class \code{marginals} which contains the following:
#'    * \code{deposition} A \code{list} of samples from the marginal density of each context's depositional date.
#'    * \code{externals} A \code{list} of samples of the marginal density of each constrant (\emph{t.p.q.} and \emph{t.a.q.]}), as conditioned upon the occurrence of other depositional 
#'    * \code{production} If a \code{finds} object has been input, samples of the marginal density of the production date of finds types will be included in the output.
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
#' result <- gibbs_ad(contexts, finds = artifacts, samples = 10^4, tpq = tpq_info, taq = taq_info)


#' @export
gibbs_ad <- function(sequences, ...) {
    UseMethod("gibbs_ad")
}

#' @rdname gibbs_ad
#' @export
gibbs_ad.list <- function(sequences, finds = NULL, samples = 10^5, tpq = NULL, taq = NULL, alpha = -5000, omega = 1950, trim = TRUE, rule = "naive") {
    if (seq_check(sequences) == TRUE) {

        proceed <- synth_rank(sequences)
        tpqnames <- numeric(length(tpq))
        taqnames <- numeric(length(taq))

        for (i in 1:length(tpq)) {
            sequences <- c(sequences, list(c(tpq[[i]]$id, tpq[[i]]$assoc))  )
            tpqnames[i] <- tpq[[i]]$id 
        }
        for (i in 1:length(taq)) {
            sequences <- c(sequences, list(c(taq[[i]]$assoc, taq[[i]]$id)) )
            taqnames[i] <- taq[[i]]$id
        }

        qa <- quae_antea(sequences)
        qp <- quae_postea(sequences)
        sampled.tpq <- list()
        sampled.taq <- list()
        sampled <- list()
        for (i in tpq) {
            tpq0 <- min(i$samples)          # initialize tpq with earliest possible
            starting <- numeric(samples)
            starting[1] <- tpq0
            sampled[[i$id]] <- starting
        }
        for (i in taq) {
            taq0 <- max(i$samples)          # initialize taq with latest possible
            starting <- numeric(samples)
            starting[1] <- taq0
            sampled[[i$id]] <- starting
        }
        for (i in proceed) {
            starting <- numeric(samples)
            starting[1] <- 9999
            sampled[[i]] <- starting
        }
        sampled[["omega"]] <- rep(omega, samples)
        sampled[["alpha"]] <- rep(alpha, samples)

        # initialize relative values
        for (i in proceed) {
            postea <- qp[[i]]
            antea <- qa[[i]]
            a <- rep(9999, length(antea))
            for (j in 1:length(antea)) {
                idx <- antea[j]
                tmp <- sampled[[idx]][1]
                if (sampled[[idx]][1] > omega)  {
                    tmp <- sampled[[idx]][1] * (-1)
                }
                a[j] <- tmp
            }
            p <- rep(9999, length(postea))
            for (j in 1:length(postea)) {
                idx <- postea[j]
                tmp <- sampled[[idx]][1]
                p[j] <- tmp
            }     
            L <- max(a, na.rm = TRUE)
            U <- min(p, na.rm = TRUE)
            s <- stats::runif(1, L, U)
            sampled[[i]][1] <- s
        }

        sampledlength <- length(sampled)

        proceed_all <- names(sampled)

        tpq_idxs <- match(tpqnames, proceed_all)
        taq_idxs <- match(taqnames, proceed_all)
        externals <- c(tpq_idxs, taq_idxs)

        proceed_idx <- match(proceed, proceed_all)

        # convert to indices
        qa_idx <- list()
        for (i in 1:length(proceed_all)) {
            qa_idx[[match(proceed_all[i], proceed_all)]] <- match(qa[[proceed_all[i]]]  , proceed_all)
        }
        qp_idx <- list()
        for (i in 1:length(proceed_all)) {
            qp_idx[[match(proceed_all[i], proceed_all)]] <- match(qp[[proceed_all[i]]]  , proceed_all)
        }
        gibbs <- matrix(0, nrow = sampledlength, ncol = samples)
        
        for (i in 1:length(proceed_all)) {
            gibbs[i,1] <- sampled[[proceed_all[i]]][1]
        }

        tpqlength <- length(tpq)
        taqlength <- length(taq)
        proceedlength <- length(proceed)

        # use indicator matrix instead of list for contexts below / above
        PhiMatrix <- matrix(0, nrow = length(qp_idx), ncol = length(qp_idx))
        for (i in 1:length(qp_idx)) {
            for (j in qp_idx[i]) {
                PhiMatrix[i, j] <- 1
            }
        }
        PsiMatrix <- matrix(0, nrow = length(qa_idx), ncol = length(qa_idx))
        for (i in 1:length(qa_idx)) {
            for (j in qa_idx[i]) {
                PsiMatrix[i, j] <- 1
            }
        }

        gibbs <- gibbs_ad_cpp(gibbs, tpq_idxs, PhiMatrix, tpq , taq_idxs, PsiMatrix, taq, proceed_idx)

        deposition <- list()
        externals <- list()
        production <- list()

        for (i in 1:proceedlength) {
            iname <- proceed[i]
            idx <- proceed_idx[i]
            if (trim == TRUE) {
                check <- TRUE
                check1 <- is.na( match(qa_idx[[idx]], tpq_idxs) )
                check2 <- is.na( match(qp_idx[[idx]], taq_idxs) )
                if ( (FALSE %in% check1) & (FALSE %in% check2) ) {
                    deposition[[iname]] <- gibbs[idx,]
                } 
            } else {
                deposition[[iname]] <- gibbs[idx,]
            }
        }
        for (i in 1:tpqlength) {
            ii <- tpq[[i]]
            tpq_idx <- match(ii$id, proceed_all)
            externals[[ii$id]] <- gibbs[tpq_idx,]
        }
        for (i in 1:taqlength) {
            ii <- taq[[i]]
            taq_idx <- match(ii$id, proceed_all)
            externals[[ii$id]] <- gibbs[taq_idx,]
        }


        # production dates 
        if (!is.null(finds)) {
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

            attestation <- matrix(0, nrow = sampledlength, ncol = findstypeslength)
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
            for (i in 1:tpqlength) {
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
            for (i in 1:taqlength) {
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

                    production[[findstypes[i]]] <- as.vector(out)
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
                        
                        production[[findstypes[i]]] <- out
                }

            } else {
                production[[findstypes[i]]] <- NULL
                warning('Invalid rule, with NULL given for production dates. Options are "naive", "earliest". ')
            }

        result <- list(deposition = deposition, externals = externals, production = production)
        class(result) <- c("marginals", "list")
        return(result)
        } else {
            result <- list(deposition = deposition, externals = externals  )
            class(result) <- c("marginals", "list")
            return(result)
        }

    }
}



