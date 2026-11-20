# Inspect package availability and exported function interfaces -----------
#
# Exploratory only. This script does not assume that signatures or returned
# object structures are identical across package versions.

packages <- c("NADA", "NADA2", "EnvStats", "survival")

package_status <- data.frame(
  package = packages,
  installed = vapply(
    packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
)

print(package_status)

if (!all(package_status$installed)) {
  missing_packages <- package_status$package[!package_status$installed]
  message(
    "Install missing packages before continuing: ",
    paste(missing_packages, collapse = ", ")
  )
}

targets <- list(
  NADA = c("cenfit", "cenmle"),
  NADA2 = c("cfit", "ros", "ROSci"),
  EnvStats = c("enparCensored")
)

for (pkg in names(targets)) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    next
  }

  cat(
    "\n--- ", pkg, " ",
    as.character(utils::packageVersion(pkg)),
    " ---\n",
    sep = ""
  )

  for (fun in targets[[pkg]]) {
    exported <- fun %in% getNamespaceExports(pkg)

    cat(
      "\nFunction: ", pkg, "::", fun,
      " | exported: ", exported, "\n",
      sep = ""
    )

    if (exported) {
      f <- getExportedValue(pkg, fun)
      print(formals(f))
    }
  }
}

cat("\n--- Session information ---\n")
print(utils::sessionInfo())
