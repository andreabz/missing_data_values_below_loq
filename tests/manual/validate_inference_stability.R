# Targeted endpoint stability check only; this does not rerun coverage metrics.
for (file in c("R/02_method_results.R", "R/04_method_ros.R", "R/05_method_mle.R",
               "R/06_method_reverse_km.R", "R/07_inference.R")) source(file)
ensure_censoring_packages_attached()

seed <- 20261121L
cases_per_condition <- 30L
B_reference <- 99L
B_supplement <- 499L
mu <- 0; sigma <- 1
sample_sizes <- c(30L)
censoring_probabilities <- c(0.10, 0.40)
loqs <- stats::qlnorm(censoring_probabilities, meanlog = mu, sdlog = sigma)
target <- exp(mu + sigma^2 / 2)
out <- list()

for (j in seq_along(censoring_probabilities)) {
  cp <- censoring_probabilities[j]
  loq <- loqs[j]
  key <- paste0("n30_cp", cp)
  method_changes <- list(mle = matrix(NA_real_, cases_per_condition, 2L,
    dimnames = list(NULL, c("lower", "upper"))),
    ros_percentile = matrix(NA_real_, cases_per_condition, 2L,
      dimnames = list(NULL, c("lower", "upper"))),
    ros_basic = matrix(NA_real_, cases_per_condition, 2L,
      dimnames = list(NULL, c("lower", "upper"))))
  covered <- matrix(FALSE, cases_per_condition, 6L,
    dimnames = list(NULL, c("mle_B99", "mle_B499", "ros_pct_B99", "ros_pct_B499",
      "ros_basic_B99", "ros_basic_B499")))

  for (i in seq_len(cases_per_condition)) {
    condition_id <- j
    set.seed(seed + condition_id * 100000L + i)
    truth <- stats::rlnorm(30L, mu, sigma)
    cens <- truth <= loq
    obs <- data.frame(value = ifelse(cens, loq, truth), censored = cens,
      loq = rep(loq, length(truth)))

    mf <- apply_mle_method(obs)
    rf <- apply_ros_method(obs)
    for (run in c("reference", "supplement")) {
      B <- if (run == "reference") B_reference else B_supplement
      # Reset to a common seed so the B=99 draws are a prefix of B=499.
      set.seed(seed + 9000000L + condition_id * 10000L + i)
      mb <- bootstrap_mle_mean(obs, mf$parameters$meanlog, mf$parameters$sdlog,
        B, conf.level = 0.95)
      set.seed(seed + 9200000L + condition_id * 10000L + i)
      rb <- bootstrap_ros_mean(obs, rf$estimates$mean, B, conf.level = 0.95)
      if (run == "reference") {
        method_changes$mle[i, ] <- mb$ci
        method_changes$ros_percentile[i, ] <- rb$percentile
        method_changes$ros_basic[i, ] <- rb$basic
        covered[i, "mle_B99"] <- mb$ci[[1L]] <= target && target <= mb$ci[[2L]]
        covered[i, "ros_pct_B99"] <- rb$percentile[[1L]] <= target && target <= rb$percentile[[2L]]
        covered[i, "ros_basic_B99"] <- rb$basic[[1L]] <= target && target <= rb$basic[[2L]]
      } else {
        method_changes$mle[i, ] <- method_changes$mle[i, ] - mb$ci
        method_changes$ros_percentile[i, ] <- method_changes$ros_percentile[i, ] - rb$percentile
        method_changes$ros_basic[i, ] <- method_changes$ros_basic[i, ] - rb$basic
        covered[i, "mle_B499"] <- mb$ci[[1L]] <= target && target <= mb$ci[[2L]]
        covered[i, "ros_pct_B499"] <- rb$percentile[[1L]] <= target && target <= rb$percentile[[2L]]
        covered[i, "ros_basic_B499"] <- rb$basic[[1L]] <= target && target <= rb$basic[[2L]]
      }
    }
  }
  out[[key]] <- list(changes = lapply(method_changes, abs),
    coverage = colSums(covered), n = cases_per_condition)
}

fmt <- function(x) formatC(x, digits = 4, format = "f")
lines <- c("# Targeted bootstrap endpoint stability check", "",
  "This supplement reuses the protocol lognormal generator and deterministic sample seeds. It reruns only the two n=30 conditions with the largest observed ROS undercoverage. It does not recompute the main simulation metrics or modify their results.", "",
  paste0("- Samples per condition: ", cases_per_condition, "; conditions: n=30, F(LOQ)=0.10 and 0.40."),
  paste0("- Paired bootstrap sizes: B=", B_reference, " and B=", B_supplement,
    "; common random-number seeds make the smaller run a prefix of the larger run."),
  paste0("- Additional bootstrap work: ", 2 * cases_per_condition * (B_reference + B_supplement),
    " refits per family (MLE or ROS); ", 4 * cases_per_condition * (B_reference + B_supplement),
    " total refits across MLE and ROS."),
  "- Endpoint differences are reported separately for lower and upper endpoints; coverage counts are descriptive for these 30 paired samples and are not fresh validation estimates.", "",
  "| Condition | Procedure / endpoint | Median absolute change | 90th percentile | Maximum |",
  "|---|---|---:|---:|---:|")
for (key in names(out)) for (method in names(out[[key]]$changes)) {
  x <- out[[key]]$changes[[method]]
  for (endpoint in colnames(x)) {
    summary <- c(median = stats::median(x[, endpoint]),
      p90 = unname(stats::quantile(x[, endpoint], .90)), max = max(x[, endpoint]))
    label <- switch(method, mle = "MLE parametric percentile",
      ros_percentile = "ROS percentile", ros_basic = "ROS basic")
    lines <- c(lines, paste0("| ", key, " | ", label, " ", endpoint, " endpoint | ",
      fmt(summary[["median"]]), " | ", fmt(summary[["p90"]]), " | ",
      fmt(summary[["max"]]), " |"))
  }
}
lines <- c(lines, "", "## Empirical inclusion counts in the paired subset", "",
  "| Condition | Procedure | B=99 | B=499 |", "|---|---|---:|---:|")
for (key in names(out)) {
  x <- out[[key]]$coverage
  lines <- c(lines,
    paste0("| ", key, " | MLE parametric percentile | ", x[["mle_B99"]], "/", cases_per_condition,
      " | ", x[["mle_B499"]], "/", cases_per_condition, " |"),
    paste0("| ", key, " | ROS percentile | ", x[["ros_pct_B99"]], "/", cases_per_condition,
      " | ", x[["ros_pct_B499"]], "/", cases_per_condition, " |"),
    paste0("| ", key, " | ROS basic | ", x[["ros_basic_B99"]], "/", cases_per_condition,
      " | ", x[["ros_basic_B499"]], "/", cases_per_condition, " |"))
}
writeLines(lines, "reports/inference_bootstrap_stability.md")
cat("Wrote reports/inference_bootstrap_stability.md\n")
