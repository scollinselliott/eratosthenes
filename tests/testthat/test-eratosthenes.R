
# vectors
w <- c("I", "G", "C")
x <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J")
y <- c("B", "D", "G", "H", "K")
z <- c("F", "K", "L", "M")
contexts <- sequences(events(x), events(y), events(z))

# finds
f1 <- finds(id = "find01", assoc = "D", type = c("type1", "form1"))
f2 <- finds(id = "find02", assoc = "E", type = c("type1", "form2"))
f3 <- finds(id = "find03", assoc = "G", type = c("type1", "form1"), residual = TRUE)
f4 <- finds(id = "find04", assoc = "H", type = c("type2", "form1"))
f5 <- finds(id = "find05", assoc = "I", type = "type2")
f6 <- finds(id = "find06", assoc = "H", type = NULL) 
artifacts <- assemblage(f1, f2, f3, f4, f5, f6)
 
# external constraints
coin1 <- absolute(id = "coin1", assoc = "B", type = NULL, samples = runif(100,-320,-300))
coin2 <- absolute(id = "coin2", assoc = "G", type = "RIC2 57", samples = seq(37, 41, length = 100))
destr <- absolute(id = "destr", assoc = "J", type = NULL, samples = 79)
tpq_info <- constraints(coin1, coin2)
taq_info <- constraints(destr)

# events creation
test_that("events works", {
    events_ <- events("A","B","C")
    expect_s3_class(events_, c("events", "character"))
    expect_error(events())
    expect_error(events(1, 2, 3, 4))
    expect_error(events("A", "B", "C", "A"))
    expect_error(events("alpha", "B", "C", "A"))
    expect_error(events("A", "omega", "C", "A"))
})

# sequence-related functions
test_that("sequences works", {
    expect_s3_class(sequences(events(x), events(y), events(z)), c("sequences", "list"))
    expect_error(sequences())
    expect_error(sequences(events(w), events(x), events(y)))
})

test_that("synth_rank works", {
  a <- sequences(events(x), events(y), events(z))
  expect_s3_class(synth_rank(a), c("events", "character"))
})

test_that("quae_postea returns list", {
  a <- sequences(events(x), events(y), events(z))
  expect_type(quae_postea(events(x), events(y), events(z)), "list")
  expect_type(quae_postea(a), "list")
  expect_error(quae_postea(events(w), events(x), events(y)))
  expect_error(quae_postea(events("A", "B", "C"), events("C", "B", "A")))
})

test_that("quae_antea returns list", {
  a <- sequences(events(x), events(y), events(z))
  expect_type(quae_antea(events(x), events(y), events(z)), "list")
  expect_type(quae_antea(a), "list")
  expect_error(quae_antea(events(w), events(x), events(y)))
  expect_error(quae_antea(events("A", "B", "C"), events("C", "B", "A")))
})

test_that("seq_adj works", {
  x1 <- events(x)
  x2 <- events(w)
  x3 <- events("F", "C")
  expect_equal(seq_adj(x1, x2), events("A", "I", "H", "B", "G", "F", "E", "J", "D", "C"))
  expect_error(seq_adj(x1, x3))
})

# finds

test_that("finds works", {
  f1 <- finds(id = "find01", assoc = "D", type = c("type1", "form1"))
  expect_s3_class(f1, c("finds", "list"))
  expect_error(finds())
  expect_error(finds(assoc = "D", type = c("type1", "form1")))
})

test_that("assemblage works", {
  expect_s3_class(assemblage(f1, f2, f3, f4, f5, f6), c("assemblage", "list"))
  expect_error(assemblage())
  expect_error(assemblage(f1, f1, f3, f4, f5, f6))
})

test_that("id_of_types works", {
  g1 <- ids_of_types(artifacts, "form2")
  g2 <- ids_of_types(artifacts, c("form1", "type1"))

  expect_true("find02" %in% g1)
  expect_equal(sum(c("find01", "find03", "find04", "find02") %in% g2), 4)

  expect_error(ids_of_types(artifacts, "foo"))
  expect_error(ids_of_types(artifacts))
  expect_error(ids_of_types())
})

# constraints

test_that("absolute works", {
  coin1 <- absolute(id = "coin1", assoc = "B", type = NULL, samples = runif(100,-320,-300))
  expect_s3_class(coin1, c("absolute", "list"))
  expect_error(absolute())
  expect_error(absolute(assoc = "B", type = NULL, samples = runif(100,-320,-300)))
  expect_error(absolute(id = "coin1", assoc = "B"))
})

test_that("constraints works", {
  expect_s3_class(constraints(coin1, coin2), c("assemblage", "list"))
  expect_error(constraints())
  expect_error(constraints(coin1, coin2, coin2))
})

# s3 classes

test_that("gibbs_ad works", {
  result <- gibbs_ad(contexts, tpq = tpq_info, taq = taq_info)
  result_noconstraints <- gibbs_ad(contexts)

  expect_s3_class(result, c("marginals", "list"))
  expect_s3_class(result_noconstraints, c("marginals", "list"))

  expect_error(gibbs_ad(contexts, alpha_ = 5000, omega_ = -5000))
})

test_that("gibbs_ad_type works", {
  result_ids <- gibbs_ad_type(contexts, artifacts, id = c("find04", "find05"), tpq = tpq_info, taq = taq_info)
  result_type <- gibbs_ad_type(contexts, artifacts, type = "type1", tpq = tpq_info, taq = taq_info)

  expect_s3_class(result_ids, c("type_marginals", "list"))
  expect_s3_class(result_type, c("type_marginals", "list"))
})

test_that("msd works", {
  result <- gibbs_ad(contexts, tpq = tpq_info, taq = taq_info)
  result_msd <- msd(result, contexts, max_samples = 5000, mcse_crit = 2, tpq = tpq_info, taq = taq_info)
  expect_s3_class(result_msd, c("msd_data", "list"))
})

test_that("sq_disp works", {
  result <- gibbs_ad(contexts, tpq = tpq_info, taq = taq_info)
  result_type <- gibbs_ad_type(contexts, artifacts, type = "type1", tpq = tpq_info, taq = taq_info)

  result_sqd <- sq_disp(result, target = "E", sequences = contexts, max_samples = 3000, mcse_crit = 2, tpq = tpq_info, taq = taq_info)
  result_sqd_type <- sq_disp(result_type, sequences = contexts, finds = artifacts, max_samples = 3000, mcse_crit = 2, tpq = tpq_info, taq = taq_info)

  expect_s3_class(result_sqd, c("sq_displ_data", "list"))
  expect_s3_class(result_sqd_type, c("sq_displ_data", "list"))
})

# numerical results

seq_ <- sequences(events("A", "B", "C"))

test_that("mc estimates are correct", {
  result <- gibbs_ad(seq_, alpha_ = 0, omega_ = 1, mcse_crit = 0.001, max_samples = 10^6)
  expect_equal(round(mean(result$deposition$B), 2), 0.50)
})

test_that("mean mcse less than critical threshold", {
  result <- gibbs_ad(seq_, alpha_ = 0, omega_ = 1, mcse_crit = 0.0005, max_samples = 10^6)
  mcse_ <- result$mcse[synth_rank(seq_)]
  expect_lt(mean(mcse_), 0.0005)
})  

test_that("msd estimates", {
  result <- gibbs_ad(seq_, alpha_ = 0, omega_ = 1, mcse_crit = 0.001, max_samples = 10^6)
  result_msd <- msd(result, seq_, alpha_ = 0, omega_ = 1, mcse_crit = 0.001, max_samples = 10^6)
  expect_equal(result_msd$MSD$MSD,  c(0.125, 0.084, 0.125, 0.000, 0.000), tolerance=2e-1)
})

test_that("sq displacement estimates are correct", {
  result <- gibbs_ad(seq_, alpha_ = 0, omega_ = 1, mcse_crit = 0.001, max_samples = 10^6)
  result_sqd <- sq_disp(result, target = "B", sequences = seq_, alpha_ = 0, omega_ = 1, mcse_crit = 0.001, max_samples = 10^6)
  expect_equal(result_sqd$sq_disp$sq_disp,  c(0.165, NA, 0.165, 0.000, 0.000), tolerance=2e-1)
})







