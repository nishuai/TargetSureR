#' Normalize an siRNA strand to a canonical guide sequence
#'
#' Converts a reported siRNA strand into the canonical guide (antisense)
#' representation used throughout TargetSureR: uppercase RNA, 5'->3', with
#' any terminal \code{TT}/\code{UU} overhang removed and the length capped
#' at \code{target_len} nucleotides.
#'
#' \strong{Why this matters.} \code{\link{annotate_sites}} expects the
#' \emph{guide} (antisense) strand, because it searches the transcriptome
#' with \code{reverseComplement(sequence)}. Supplying a sense strand does
#' not raise an error; it silently yields seed-only hits and no
#' on-target \code{full_complementarity} match. The built-in
#' \code{\link{reference_set}} cohort was preprocessed with this function
#' (all 94 guides are 19 nt), so user siRNAs should be passed through it
#' before scanning if the results are to be compared against that cohort
#' via \code{reference = "builtin"}.
#'
#' \strong{Truncation is 3'-directed.} Capping at \code{target_len} keeps
#' the first \code{target_len} nucleotides, trimming the 3' end. The seed
#' region (guide positions 2-8, used by \pkg{SeedMatchR}) sits at the 5'
#' end and is therefore never altered.
#'
#' Reported siRNA duplexes vary in how the 3' overhang is written. Entries
#' from the MIT/ICBP siRNA database spell it out (\code{...AACAdTdT}),
#' whereas many publications report a bare 21-mer whose final two bases
#' are genuine target-matching nucleotides. Only a literal terminal
#' \code{TT}/\code{UU} is stripped; anything else is handled by the length
#' cap alone.
#'
#' @param x Character vector of siRNA sequences. DNA or RNA, any case;
#'   non-\code{ACGTU} characters (e.g. the \code{d} in \code{dTdT}) are
#'   removed.
#' @param input_strand Which strand \code{x} holds. \code{"guide"}
#'   (default) treats \code{x} as the antisense/guide strand and returns
#'   it as-is after cleaning. \code{"sense"} reverse-complements \code{x}
#'   to obtain the guide.
#' @param target_len Integer. Maximum guide length in nucleotides;
#'   sequences longer than this are truncated at the 3' end. Default
#'   \code{19}, matching the built-in \code{\link{reference_set}}.
#' @param min_len Integer. Sequences shorter than this after cleaning
#'   return \code{NA_character_}. Default \code{8}.
#'
#' @return Character vector of the same length as \code{x}, holding
#'   uppercase RNA guide sequences (5'->3'), or \code{NA_character_} where
#'   the input was empty, \code{NA}, or shorter than \code{min_len}.
#'
#' @examples
#' # MIT/ICBP-style sense strand with a written-out dTdT overhang
#' normalize_guide("GUUUUCACUCCAGCUAACAdTdT", input_strand = "sense")
#'
#' # A published 21-mer guide: capped to 19 nt, seed (positions 2-8) intact
#' normalize_guide("UCCAUAACUUCUUGCUAAGUC")
#'
#' @seealso \code{\link{annotate_sites}} for the scan that consumes the
#'   guide, \code{\link{reference_set}} for the cohort built with this
#'   function.
#' @export
normalize_guide <- function(x,
                            input_strand = c("guide", "sense"),
                            target_len = 19L,
                            min_len = 8L) {
  input_strand <- match.arg(input_strand)
  if (!is.numeric(target_len) || length(target_len) != 1L || target_len < 1)
    stop("target_len must be a single positive integer", call. = FALSE)
  if (!is.numeric(min_len) || length(min_len) != 1L || min_len < 1)
    stop("min_len must be a single positive integer", call. = FALSE)
  if (length(x) == 0L) return(character(0))

  s <- toupper(as.character(x))
  s <- gsub("[^ACGTU]", "", s)          # drops the 'd' of dTdT, spaces, etc.
  s <- sub("(TT|UU)$", "", s)           # strip a literal 3' overhang
  s[is.na(x) | !nzchar(s)] <- NA_character_

  if (input_strand == "sense") s <- .revcomp_rna(s)

  s <- chartr("T", "U", s)
  long <- !is.na(s) & nchar(s) > target_len
  s[long] <- substr(s[long], 1L, target_len)   # trim the 3' end; seed is 5'
  s[!is.na(s) & nchar(s) < min_len] <- NA_character_
  s
}

# Reverse complement for RNA/DNA held as character, without a Biostrings
# dependency (Biostrings is only a Suggests).
.revcomp_rna <- function(s) {
  out <- rep(NA_character_, length(s))
  keep <- !is.na(s) & nzchar(s)
  if (!any(keep)) return(out)
  comp <- chartr("ACGTU", "UGCAA", s[keep])
  out[keep] <- vapply(strsplit(comp, "", fixed = TRUE),
                      function(ch) paste(rev(ch), collapse = ""),
                      character(1))
  out
}
