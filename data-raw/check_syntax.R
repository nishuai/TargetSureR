for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) {
  cat(basename(f), ": ")
  res <- tryCatch({
    parse(file = f)
    "OK"
  }, error = function(e) sprintf("FAIL - %s", conditionMessage(e)))
  cat(res, "\n")
}
