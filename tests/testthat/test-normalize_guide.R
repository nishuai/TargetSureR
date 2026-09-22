test_that("terminal overhangs are stripped and case/alphabet normalised", {
  # MIT/ICBP style: dTdT written out; 'd' is not an ACGTU character
  expect_equal(normalize_guide("GUUUUCACUCCAGCUAACAdTdT"), "GUUUUCACUCCAGCUAACA")
  expect_equal(normalize_guide("gccaagaaguuuccuaauauu"), "GCCAAGAAGUUUCCUAAUA")
  # DNA input is returned as RNA
  expect_equal(normalize_guide("TGAGGGATATCGCCAAACAdTdT"), "UGAGGGAUAUCGCCAAACA")
})

test_that("truncation trims the 3' end and never touches the seed", {
  q <- c("UCCAUAACUUCUUGCUAAGUC",   # 21-mer, no TT/UU overhang
         "UUGUCUUUGCUGAUGUUUCAA",
         "AUUACACACUUUGUCUUUGAC")
  g <- normalize_guide(q)
  expect_true(all(nchar(g) == 19L))
  # guide is a 5' prefix of the input
  expect_equal(g, substr(q, 1L, 19L))
  # SeedMatchR mer7m8 seed = guide positions 2-8
  expect_equal(substr(g, 2L, 8L), substr(q, 2L, 8L))
})

test_that("sense input is reverse-complemented to the guide", {
  # MIT entry #1000 (CXCR4) reports both strands; deriving from sense must
  # reproduce the reported antisense, minus its dTdT overhang.
  expect_equal(normalize_guide("GUUUUCACUCCAGCUAACAdTdT", input_strand = "sense"),
               "UGUUAGCUGGAGUGAAAAC")
  expect_equal(normalize_guide("UGUUAGCUGGAGUGAAAACdTdT", input_strand = "guide"),
               "UGUUAGCUGGAGUGAAAAC")
})

test_that("the 94 shipped guides are reproduced exactly", {
  md <- reference_set$metadata
  # Idempotence: re-normalising a stored guide is a no-op
  expect_equal(normalize_guide(md$guide), md$guide)
})

test_that("degenerate input yields NA rather than an error", {
  expect_true(is.na(normalize_guide(NA_character_)))
  expect_true(is.na(normalize_guide("")))
  expect_true(is.na(normalize_guide("xyz")))
  expect_true(is.na(normalize_guide("ACGU")))          # shorter than min_len
  expect_equal(normalize_guide("ACGU", min_len = 4L), "ACGU")
  expect_equal(normalize_guide(character(0)), character(0))
  expect_equal(normalize_guide(c("ACGUACGUACGU", NA)),
               c("ACGUACGUACGU", NA))
})

test_that("target_len and min_len are validated", {
  expect_error(normalize_guide("ACGUACGUACGU", target_len = 0))
  expect_error(normalize_guide("ACGUACGUACGU", min_len = -1))
  expect_error(normalize_guide("ACGUACGUACGU", input_strand = "nonsense"))
  expect_equal(nchar(normalize_guide("UCCAUAACUUCUUGCUAAGUC", target_len = 21L)), 21L)
})
