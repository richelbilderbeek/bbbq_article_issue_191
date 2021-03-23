library(dplyr, warn.conflicts = FALSE)

# How many AAs before and after the epitope
n <- 10

#' Analyse the epitopes
#' @param epitope_sequences sequences of epitopes
#' @param positions string positions as used by \link[stringr]{str_sub}, e.g. '1' denotes the first character,
#'   '-1' means the last character, '-2' the one-but-last
analyse_epitopes <- function(
  epitope_sequences,
  positions
) {
  testthat::expect_true(is.character(epitope_sequences))
  testthat::expect_true(is.numeric(positions))
  # Store sequences as tibble, as stringr needs it
  char_tibble <- tibble::tibble(
    epitope_sequence = epitope_sequences
  )

  # Collect all tibbles in this list
  tibbles <- list()

  # Count all chararacters at all positions
  for (i in seq_along(positions)) {
    pos <- positions[i]
    testthat::expect_true(all(nchar(char_tibble$epitope_sequence) >= abs(pos)))
    char_tibble$char <- stringr::str_sub(char_tibble$epitope_sequence, pos, pos)
    t <- dplyr::count(char_tibble, char)
    t$pos <- pos
    t <- dplyr::relocate(t, pos)
    tibbles[[i]] <- t
  }
  # Combine the reults per position
  t_all <- dplyr::bind_rows(tibbles)

  # Convert the data to wide format for Libreoffice Calc users
  t_wide <- tidyr::pivot_wider(t_all, names_from = pos, values_from = n) %>% dplyr::arrange(char)

  # Make the column names reader friendly
  # Reverse the names, as requested
  names(t_wide) <- c(
    names(t_wide)[1],
    rev(stringr::str_replace(english::as.english(-as.numeric(names(t_wide)[-1])), " ", "_"))
  )

  return(t_wide)
}

# Collect all the MHC-I matches
matches_filename <- "../bbbq_article_issue_157/matches_1.csv"
testthat::expect_true(file.exists(matches_filename))
t <- readr::read_csv(matches_filename)
# Keep the uniquely mapped epitopes
t_unique <- t %>% dplyr::filter(!is.na(gene_name))

# Get the n AAs before the sequence
sequences_before_epitopes <- stringr::str_match(
  string = t_unique$sequence,
  pattern = paste0("[[:upper:]]{", n, "}", t_unique$epitope_sequence)
)[, 1]
sequences_before_epitopes <- sequences_before_epitopes[!is.na(sequences_before_epitopes)]
testthat::expect_equal(6847, length(sequences_before_epitopes))
knitr::kable(head(sequences_before_epitopes))
readr::write_lines(
  sequences_before_epitopes,
  paste0("~/", english::as.english(n) ,"_before_epitopes.txt")
)
readr::write_csv(
  analyse_epitopes(
    epitope_sequences = sequences_before_epitopes,
    positions = seq(1, n)
  ),
  paste0("~/", english::as.english(n) ,"_before_epitopes.csv")
)
t_sequences_before_epitopes <- tibble::tibble(
  name = paste0("iloverichel", seq_len(length(sequences_before_epitopes))),
  sequence = stringr::str_sub(sequences_before_epitopes, 1, n)
)
pureseqtmr::save_tibble_as_fasta_file(
  t_sequences_before_epitopes,
  paste0("~/", english::as.english(n) ,"_before_epitopes.fasta")
)

###############################################################################
# The n epitopes after
###############################################################################
sequences_after_epitopes <- stringr::str_match(
  string = t_unique$sequence,
  pattern = paste0(t_unique$epitope_sequence, "[[:upper:]]{", n, "}")
)[, 1]
sequences_after_epitopes <- sequences_after_epitopes[!is.na(sequences_after_epitopes)]
testthat::expect_equal(6672, length(sequences_after_epitopes))

readr::write_lines(
  sequences_after_epitopes,
  paste0("~/", english::as.english(n) ,"_after_epitopes.txt")
)
readr::write_csv(
  analyse_epitopes(sequences_after_epitopes, positions = seq(-n, -1)),
  paste0("~/", english::as.english(n) ,"_after_epitopes.csv")
)

t_sequences_after_epitopes <- tibble::tibble(
  name = paste0("iloverichel", seq_len(length(sequences_after_epitopes))),
  sequence = stringr::str_sub(sequences_after_epitopes, -n)
)
pureseqtmr::save_tibble_as_fasta_file(
  t_sequences_after_epitopes,
  paste0("~/", english::as.english(n) ,"_after_epitopes.fasta")
)

# First epitope:
# MEKSSLTQHSW, see to the right where it is alignmed
# ---->                                                                                                                                                                 MEKSSLTQHSW                                                                                                                                                                                                                               # nolint indeed long
# MAEAMDLGKDPNGPTHSSTLFVRDDGSSMSFYVRPSPAKRRLSTLILHGGGTVCRVQEPGAVLLAQPGEALAEASGDFISTQYILDCVERNERLELEAYRLGPASAADTGSEAKPGALAEGAAEPEPQRHAGRIAFTDADDVAILTYVKENARSPSSVTGNALWKAMEKSSLTQHSWQSLKDRYLKHLRGQEHKYLLGDAPVSPSSQKLKRKAEEDPEAADSGEPQNKRTPDLPEEEYVKEEIQENEEAVKKMLVEATREFEEVVVDESPPDFEIHITMCDDDPPTPEEDSETQPDEEEEEEEEKVSQPEVGAAIKIIRQLMEKFNLDLSTVTQAFLKNSGELEATSAFLASGQRADGYPIWSRQDDIDLQKDDEDTREALVKKFGAQNVARRIEFRKK # nolint indeed long
