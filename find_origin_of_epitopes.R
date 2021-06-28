# Also, count the percentage of fragments that can be mapped to a proteome (i.e. including duplicates).
t_schellens <- bianchietal2017::get_schellens_et_al_2015_sup_1(
  xlsx_filename = bianchietal2017::download_schellens_et_al_2015_sup_1(
    xlsx_filename = "schellens_et_al_2015.xlsl"
  )
)
n_epitopes <- nrow(t_schellens)
testthat::expect_equal(10894, n_epitopes)
n_unique_epitopes <- length(unique(t_schellens$epitope_sequence))
testthat::expect_equal(7897, n_unique_epitopes)

# Bos Taurus proteome, from https://www.uniprot.org/proteomes/UP000009136
utils::download.file(
  url = "http://richelbilderbeek.nl/UP000009136.fasta",
  destfile = "UP000009136.fasta",
  quiet = TRUE
)
t_cow <- pureseqtmr::load_fasta_file_as_tibble("UP000009136.fasta")
t_human <- bbbq::get_proteome()
epitope_sequences <- unique(t_schellens$epitope_sequence)
length(epitope_sequences)
t_matches <- tibble::tibble(
  epitope_sequence = epitope_sequences,
  n_matches_human = NA,
  n_matches_cow = NA
)
for (i in seq_len(nrow(t_matches))) {
  print(paste0(i, "/", nrow(t_matches)))
  matches_human <- stringr::str_which(
    string =  t_human$sequence,
    pattern = t_matches$epitope_sequence[i]
  )
  t_matches$n_matches_human[i] <- length(matches_human)
  matches_cow <- stringr::str_which(
    string =  t_cow$sequence,
    pattern = t_matches$epitope_sequence[i]
  )
  t_matches$n_matches_cow[i] <- length(matches_cow)
}
t_matches
readr::write_csv(x = t_matches, "n_matches.csv")

text <- paste0(
  "number of unique epitopes: ", length(epitope_sequences), "\n",
  "number of epitopes map to a human proteome: ",
    sum(t_matches$n_matches_human != 0),
    "\n",
  "number of epitopes map to a cow proteome: ",
    sum(t_matches$n_matches_cow != 0),
    "\n"
)
readr::write_lines(x = text, "n_matches.txt")
