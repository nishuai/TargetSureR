## Submission
Version 1.1.0. Adds `build_annodb()`, which exports the transcript
sequence-index construction previously documented only in the vignette.

## Test environments
* Local Windows 11, R 4.5.3 (2026-03-11 ucrt), Rtools45
* `R CMD check --as-cran`: 0 errors, 0 warnings, 1 note (benign:
  "unable to verify current time" — Windows sandbox cannot reach NTP)

## R CMD check results
No ERRORs or WARNINGs.

There was 1 NOTE:
* checking for future file timestamps ... NOTE
  unable to verify current time

This note is a Windows-local network artefact; it does not reproduce
on CRAN build servers.

## Downstream dependencies
This is a new submission; no reverse dependencies exist.

## Additional notes
* Heavy Bioconductor dependencies (`SeedMatchR`, `ensembldb`,
  `AnnotationHub`, `BSgenome.Hsapiens.UCSC.hg38`, etc.) are all declared
  under `Suggests` to keep the base install light. Functions that need
  them call `requireNamespace()` at runtime and emit actionable error
  messages with installation instructions.
* The vignette documents a workflow whose inputs are multi-gigabyte
  external annotation downloads. No chunk is evaluated at build time
  (`eval = FALSE`), and each carries `purl = FALSE` so that the tangled
  script is not executed during checking either; only `sessionInfo()`
  runs. The two bundled `*_example` datasets are simulated values, so
  examples and tests run without any download.
* Package does not compile any C/C++ code.
* Bundled data totals < 10 KB (xz-compressed `.rda`).

## Author
Shuai Ni &lt;nishuai@wakerbio.com&gt;
