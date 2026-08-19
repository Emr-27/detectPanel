\name{make_repeated_stratified_splits}
\alias{make_repeated_stratified_splits}
\alias{search_panels}
\alias{print.detectPanel_search}
\title{Shared Resampling Splits and Exhaustive Panel Search}
\usage{
make_repeated_stratified_splits(outcome, v = 5L, repeats = 10L, seed = 123L)

search_panels(expression, outcome, candidates, panel_size = 3L,
  splits = NULL, candidate_table = NULL, v = 5L, repeats = 10L,
  seed = 123L, near_best_tolerance = 0.01, max_combinations = 50000L)
}
\arguments{
  \item{outcome}{Binary vector.}
  \item{v}{Number of folds.}
  \item{repeats}{Number of repeated fold assignments.}
  \item{seed}{Random seed.}
  \item{expression}{Transformed feature-by-sample matrix.}
  \item{candidates}{Candidate feature names.}
  \item{panel_size}{Number of features per panel.}
  \item{splits}{Reusable splits from \code{make_repeated_stratified_splits()}.}
  \item{candidate_table}{Optional candidate scoring table.}
  \item{near_best_tolerance}{AUC tolerance before expression-aware tie breaking.}
  \item{max_combinations}{Safety limit for exhaustive search.}
}
\value{A list of splits or a \code{detectPanel_search} object.}
\description{All panels are evaluated on the same folds to permit fair paired comparison.
A repeated-CV AUC is retained only when every sample in that repeat receives a
finite out-of-fold prediction; failed folds are not silently dropped.}
