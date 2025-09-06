
# sequences
w <- c("I", "G", "C")
x <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J")
y <- c("B", "D", "G", "H", "K")
z <- c("F", "K", "L", "M")
contexts <- list(x, y, z)
contexts_misordered <- list(w, x, y)

# finds
f1 <- list(id = "find01", assoc = "D", type = c("type1", "form1"))
f2 <- list(id = "find02", assoc = "E", type = c("type1", "form2"))
f3 <- list(id = "find03", assoc = "G", type = c("type1", "form1"), residual = TRUE)
f4 <- list(id = "find04", assoc = "H", type = c("type2", "form1"))
f5 <- list(id = "find05", assoc = "I", type = "type2")
f6 <- list(id = "find06", assoc = "H", type = NULL) 
artifacts <- list(f1, f2, f3, f4, f5, f6)
 
# external constraints
coin1 <- list(id = "coin1", assoc = "B", type = NULL, samples = runif(100,-320,-300))
coin2 <- list(id = "coin2", assoc = "G", type = NULL, samples = seq(37, 41, length = 100))
destr <- list(id = "destr", assoc = "J", type = NULL, samples = 79)
tpq_info <- list(coin1, coin2)
taq_info <- list(destr)

# sequence-related functions

test_that("seq_check works", {
  expect_true(seq_check(contexts))
  expect_false(seq_check(contexts_misordered))
})

test_that("synth_rank works", {
  expect_type(synth_rank(contexts), "character")
  expect_null(synth_rank(contexts_misordered))
  expect_message(synth_rank(contexts_misordered))
})

test_that("quae_postea returns list", {
  expect_type(quae_postea(contexts), "list")
})

test_that("quae_antea returns list", {
  expect_type(quae_antea(contexts), "list")
})

test_that("seq_adj works", {
  expect_equal(seq_adj(x, w), c("A", "I", "H", "B", "G", "F", "E", "J", "D", "C"))
  expect_null(seq_adj(x, c("X","Y","Z")))
  expect_message(seq_adj(x, c("X","Y","Z")))
})

# s3 classes

test_that("gibbs_ad works", {
  result <- gibbs_ad(contexts, tpq = tpq_info, taq = taq_info)
  result_noconstraints <- gibbs_ad(contexts)

  expect_s3_class(result, c("marginals", "list"))
  expect_s3_class(result_noconstraints, c("marginals", "list"))

  expect_error(gibbs_ad(contexts_misordered))
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

seq_ <- list(c("A", "B", "C"))

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
  result_msd <- msd(result, seq_, alpha_ = 0, omega = 1, mcse_crit = 0.001, max_samples = 10^6)
  expect_equal(result_msd$MSD$MSD,  c(0.125, 0.084, 0.125, 0.000, 0.000), tolerance=2e-1)
})

test_that("sq displacement estimates are correct", {
  result <- gibbs_ad(seq_, alpha_ = 0, omega_ = 1, mcse_crit = 0.001, max_samples = 10^6)
  result_sqd <- sq_disp(result, target = "B", sequences = seq_, alpha_ = 0, omega_ = 1, mcse_crit = 0.001, max_samples = 10^6)

  expect_equal(result_sqd$sq_disp$sq_disp,  c(0.165, NA, 0.165, 0.000, 0.000), tolerance=2e-1)
})






