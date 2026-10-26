# Perform UKA analysis on the comparison data

suppressPackageStartupMessages({
  library(purrr)
  library(dplyr)
  library(readr)
  library(stringr)
  library(pgUpstream) # nolint: unused_import_linter.
  library(pgFCS) # nolint: unused_import_linter.
  library(pgscales) # nolint: unused_import_linter.

})

## Load data

uka_db <- readRDS(
  file.path(
    "reference_data",
    "UKA_231031-86502-87102_UpstreamDb.rds"
  )
)
kinase_enrichment <- read_csv(
  file.path("reference_data", "UKA_Kinase_enrichment.csv")
)

perform_uka <- function(
    dpp_data, upstream_db, kinase_family, kinase_enrichment,
    minimum_sequence_homology = 0.9,
    minimum_phosphonet_score = 300L,
    nperms = 500L,
    min_rank = 4L,
    max_rank = 12L,
    weight_iviv = 1L,
    weight_pnet = 1L,
    minimum_set_size = 3L) {
  subset_db <- upstream_db |>
    dplyr::ungroup() |>
    dplyr::filter(
      .data$PepProtein_SeqSimilarity >= .env$minimum_sequence_homology,
      .data$family == .env$kinase_family,
      .data$Kinase_Rank <= .env$max_rank
    ) |>
    dplyr::filter(
      .data$Kinase_PKinase_PredictorVersion2Score >=
        .env$minimum_phosphonet_score |
        .data$Database == "iviv"
    )
  
  if (nrow(subset_db) == 0L) {
    stop("No database rows remain after filtering.")
  }

  uka_result <- pgScanAnalysis0_np(
    dpp_data,
    subset_db,
    nPermutations = nperms,
    scanRank = min_rank:max_rank,
    dbWeights = c(iviv = weight_iviv, PhosphoNET = weight_pnet)
  ) |>
    purrr::map(function(x) {
      x[["aResult"]][[1L]] |>
        dplyr::mutate(mxRank = x[["mxRank"]])
    }) |>
    dplyr::bind_rows() |>
    dplyr::filter(nFeatures >= minimum_set_size) |>
    makeSummary() |>
    dplyr::select(Kinase = ClassName, dplyr::everything()) |>
    dplyr::left_join(kinase_enrichment, by = c(Kinase = "Kinase_Name")) |>
    dplyr::arrange(-medianScore) |>
    dplyr::select(
      `Kinase Name` = Kinase,
      `Kinase Uniprot ID` = Kinase_UniprotID,
      `Kinase Group` = Kinase_group,
      `Kinase Family` = Kinase_family,
      `Mean Significance Score` = meanPhenoScore,
      `Mean Specificity Score` = meanFeatScore,
      `Median Final score` = medianScore,
      `Max Final score` = maxScore,
      `Median Kinase Statistic` = medianStat,
      `Mean Kinase Statistic` = meanStat,
      `SD Kinase Statitistic` = sdStat,
      `Median Kinase Change` = medianDelta,
      `Mean peptide set size` = meanSetSize
    )

  uka_result
}

prepare_signal_data <- function(signal_path) {
  readr::read_csv(signal_path) |>
    dplyr::select(Peptide, totalMeanLFC) |>
    dplyr::mutate(
      ID = as.factor(Peptide),
      value = totalMeanLFC
    ) |>
    dplyr::select(ID, value) |>
    unique()
}

## Perform analysis

signal_files <- list.files(
  path = "results",
  pattern = "dpp.*\\.csv$",
  full.names = TRUE
) |>
  purrr::keep(~ !stringr::str_detect(basename(.x), fixed("p-w")))

names(signal_files) <- basename(signal_files) |>
  stringr::str_remove(fixed("-dpp")) |>
  stringr::str_remove("\\.csv$")

# Confirm the selected files
print(signal_files)
stopifnot(length(signal_files) > 0)


signal_names <- names(signal_files)

run_prefix <- signal_files |>
  names() |>
  str_extract("(.+)_", 1L)

uka_results <- signal_files |>
  map(prepare_signal_data) |>
  imap(~ perform_uka(
    .x,
    uka_db,
    {
      .y |> str_extract(".TK")
    },
    kinase_enrichment
  ))

written <- pmap(
  list(uka_results, names(uka_results), run_prefix),
  \(data, name, prefix) {
    data |> write_csv(file.path(
      "results",
      str_glue("{prefix}-uka_table_full_{name}.csv")
    ))
  }
)

