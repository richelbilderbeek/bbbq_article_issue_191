# Count number of recursive epitopes
# For https://github.com/richelbilderbeek/bbbq_article_issue_191/issues/1

xlsx_filename <- "schellens_et_al_2015_sup_1.xlsx"
bianchietal2017::download_schellens_et_al_2015_sup_1(
  url = "http://richelbilderbeek.nl/schellens_et_al_2015_s_1.xlsx",
  xlsx_filename = xlsx_filename,
  verbose = FALSE
)
t_schellens <- bianchietal2017::get_schellens_et_al_2015_sup_1(
  xlsx_filename = xlsx_filename
)
epitope_sequences <- unique(t_schellens$epitope_sequence)
n_epitope_sequences <- length(epitope_sequences)
testthat::expect_equal(n_epitope_sequences, 7897)

# Half as comprisons are one-way, i.e. smaller sequences can be a subset of
# larger, but not the other way around
n_comparisons <- 0.5 * n_epitope_sequences * (n_epitope_sequences - 1)
n_comparisons
testthat::expect_equal(n_comparisons, 31177356)

# Takes approx 15 secs
tibbles <- list()
for (i in seq_len(length(epitope_sequences))) {
  focal_sequence <- epitope_sequences[i]
  super_seqs <- stringr::str_subset(string = epitope_sequences, pattern = focal_sequence)
  super_seqs <- super_seqs[super_seqs != focal_sequence]
  t <- tibble::tibble(
    epitope = focal_sequence,
    super_seq = super_seqs
  )
  tibbles[[i]] <- t
}
t <- dplyr::bind_rows(tibbles)


n_resursive <- nrow(t)
testthat::expect_equal(n_resursive, 648)

readr::write_csv(x = t, file = "recursive_epitopes.csv")

