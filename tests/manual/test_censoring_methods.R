# Controlled diagnostic of installed left-censoring APIs. Exploratory only;
# does not define or select the final estimators for the Monte Carlo study.
suppressPackageStartupMessages({
  library(NADA)
  library(NADA2)
  library(EnvStats)
})

cat("Installed versions\n")
for (pkg in c("NADA", "NADA2", "EnvStats", "survival"))
  cat(pkg, as.character(utils::packageVersion(pkg)), "\n")

# Deterministic latent lognormal values span both sides of a common LOQ.
n <- 40L
latent <- stats::qlnorm((seq_len(n) - 0.5) / n, meanlog = 0, sdlog = 1)
loq <- stats::qlnorm(0.10, meanlog = 0, sdlog = 1)
scenarios <- list(
  none = list(x = latent, censored = rep(FALSE, n)),
  intermediate = list(x = latent, censored = latent < stats::qlnorm(0.40)),
  elevated = list(x = latent, censored = latent < stats::qlnorm(0.75)),
  all = list(x = latent, censored = rep(TRUE, n))
)
# For each scenario the common censoring threshold is fixed; censored entries
# are represented at LOQ, while latent values are retained for validation.
scenarios$none$loq <- loq
scenarios$intermediate$loq <- stats::qlnorm(0.40)
scenarios$elevated$loq <- stats::qlnorm(0.75)
scenarios$all$loq <- max(latent) * 1.01

capture <- function(expr) {
  warnings <- character()
  value <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
  }), error = function(e) structure(list(message = conditionMessage(e)), class = "diag_error"))
  list(value = value, warnings = warnings)
}

for (s in names(scenarios)) {
  d <- scenarios[[s]]
  cens <- d$censored
  obs <- ifelse(cens, d$loq, d$x)
  cat("\n=====", s, " censored=", sum(cens), "/", n, " LOQ=", d$loq, "=====\n", sep = "")
  calls <- list(
    NADA_cenfit = quote(NADA::cenfit(obs, cens)),
    NADA_cenmle = quote(NADA::cenmle(obs, cens)),
    NADA2_cfit = quote(NADA2::cfit(obs, as.numeric(cens), Cdf = FALSE, printstat = FALSE)),
    NADA2_ros = quote(NADA2::ros(obs, cens, forwardT = "log", reverseT = "exp")),
    NADA2_ROSci = quote(NADA2::ROSci(NADA2::ros(obs, cens, forwardT = "log", reverseT = "exp"), printstat = FALSE)),
    EnvStats_enparCensored = quote(EnvStats::enparCensored(obs, cens, censoring.side = "left", ci = TRUE, warn = TRUE))
  )
  for (nm in names(calls)) {
    z <- capture(eval(calls[[nm]]))
    cat("--", nm, "--\n")
    if (inherits(z$value, "diag_error")) cat("ERROR:", z$value$message, "\n") else {
      v <- z$value
      if (inherits(v, "ros")) print(list(modeled_mean = mean(v$modeled), coefficients = coef(v$model), modeled = head(v$modeled)))
      else if (inherits(v, "cenfit")) print(utils::head(summary(v), 5))
      else if (inherits(v, "cenmle")) print(v)
      else if (is.list(v)) print(str(v, max.level = 2)) else print(v)
    }
    if (length(z$warnings)) cat("WARNINGS:", paste(z$warnings, collapse = " | "), "\n")
  }
}

cat("\nDocumentation help topics (installed Rd index)\n")
for (p in c("NADA", "NADA2", "EnvStats")) {
  rd <- tools::Rd_db(p)
  targets <- switch(p, NADA = c("cenfit", "cenmle"),
                    NADA2 = c("cfit", "ros", "ROSci"),
                    EnvStats = "enparCensored")
  for (topic in targets) {
    hit <- grep(paste0("(^|/)", topic, "\\.Rd$"), names(rd))
    if (!length(hit)) hit <- grep(topic, names(rd), fixed = TRUE)
    cat("\n---", p, "::", topic, "---\n")
    if (!length(hit)) cat("No matching installed Rd topic\n") else
      cat(paste(capture.output(tools::Rd2txt(rd[[hit[1]]])), collapse = "\n"), "\n")
  }
}
cat("\nSession information\n")
print(utils::sessionInfo())
