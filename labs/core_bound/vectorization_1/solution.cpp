#include "solution.hpp"
#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <vector>

using score_t = int16_t;
using column_t = std::array<score_t, sequence_size_v + 1>;

void copy_score_column(column_t& dst, column_t& src) {
  for (size_t i = 0; i < src.size(); i++) {
    dst[i] = src[i];
  }
}

using row_t = std::array<score_t, sequence_count_v>;
using matrix_t = std::array<row_t, sequence_size_v+1>;

matrix_t transpose(std::vector<sequence_t> const & seq) {
  matrix_t res;
  for (size_t row = 0; row < sequence_count_v; row++) {
    for (size_t col = 0; col < sequence_size_v; col++) {
      res[col][row] = seq[row][col];
    }
  }
  return res;
}

#ifdef SOLUTION
// The alignment algorithm which computes the alignment of the given sequence
// pairs.
result_t compute_alignment(std::vector<sequence_t> const &sequences1,
                           std::vector<sequence_t> const &sequences2) {
  result_t result;
  /*
    * Initialise score values.
    */
  const score_t gap_open{-11};
  const score_t gap_extension{-1};
  const score_t match{6};
  const score_t mismatch{-4};

  auto mat1 = transpose(sequences1);
  auto mat2 = transpose(sequences2);

  matrix_t score_column{};
  matrix_t horizontal_gap_column{};
  row_t last_vertical_gaps{};

  #define piter size_t pi = 0; pi < sequence_count_v; ++pi

  /*
   * Initialise the first column of the matrix.
   */
  #pragma GCC ivdep
  for (piter) {
    horizontal_gap_column[0][pi] = gap_open;
    last_vertical_gaps[pi] = gap_open;
  }

  for (size_t i = 1; i < score_column.size(); ++i) {
    #pragma GCC ivdep
    for (piter) {
      score_column[i][pi] = last_vertical_gaps[pi];
      horizontal_gap_column[i][pi] = last_vertical_gaps[pi] + gap_open;
      last_vertical_gaps[pi] += gap_extension;
    }
  }

  for (unsigned col2 = 0; col2 < sequence_size_v; col2++) {
    row_t last_diagonal_score = score_column[0];
    #pragma GCC ivdep
    for (piter) {
      score_column[0][pi] = horizontal_gap_column[0][pi];
      last_vertical_gaps[pi] = horizontal_gap_column[0][pi]+gap_open;
      horizontal_gap_column[0][pi] += gap_extension;
    }

    for (unsigned col1 = 0; col1 < sequence_size_v; col1++) {
      auto row = col1+1;
      row_t best_cell_score;
      #pragma GCC ivdep
      for (piter) {
        best_cell_score[pi] = 
        std::max({
          score_t(last_diagonal_score[pi] + (mat1[col1][pi] == mat2[col2][pi]?match:mismatch)),
          last_vertical_gaps[pi],
          horizontal_gap_column[row][pi]
        });
        last_diagonal_score[pi] = score_column[row][pi];
        score_column[row][pi] = best_cell_score[pi];
        best_cell_score[pi] += gap_open;
        last_vertical_gaps[pi] = std::max(
          score_t(last_vertical_gaps[pi]+gap_extension),
          best_cell_score[pi]
        );
        horizontal_gap_column[row][pi] = std::max(
          score_t(horizontal_gap_column[row][pi]+gap_extension),
          best_cell_score[pi]
        );
      }
    }
  }
  result = score_column.back();
  return result;
}

#else

// The alignment algorithm which computes the alignment of the given sequence
// pairs.
result_t compute_alignment(std::vector<sequence_t> const &sequences1,
                           std::vector<sequence_t> const &sequences2) {
  result_t result{};
  for (size_t sequence_idx = 0; sequence_idx < sequences1.size();
       ++sequence_idx) {
    using score_t = int16_t;
    using column_t = std::array<score_t, sequence_size_v + 1>;

    sequence_t const &sequence1 = sequences1[sequence_idx];
    sequence_t const &sequence2 = sequences2[sequence_idx];

    /*
     * Initialise score values.
     */
    score_t gap_open{-11};
    score_t gap_extension{-1};
    score_t match{6};
    score_t mismatch{-4};

    /*
     * Setup the matrix.
     * Note we can compute the entire matrix with just one column in memory,
     * since we are only interested in the last value of the last column in the
     * score matrix.
     */
    column_t score_column{};
    column_t horizontal_gap_column{};
    score_t last_vertical_gap{};

    /*
     * Initialise the first column of the matrix.
     */
    horizontal_gap_column[0] = gap_open;
    last_vertical_gap = gap_open;

    for (size_t i = 1; i < score_column.size(); ++i) {
      score_column[i] = last_vertical_gap;
      horizontal_gap_column[i] = last_vertical_gap + gap_open;
      last_vertical_gap += gap_extension;
    }

    /*
     * Compute the main recursion to fill the matrix.
     */
    for (unsigned col = 1; col <= sequence2.size(); ++col) {
      score_t last_diagonal_score =
          score_column[0]; // Cache last diagonal score to compute this cell.
      score_column[0] = horizontal_gap_column[0];
      last_vertical_gap = horizontal_gap_column[0] + gap_open;
      horizontal_gap_column[0] += gap_extension;

      for (unsigned row = 1; row <= sequence1.size(); ++row) {
        // Compute next score from diagonal direction with match/mismatch.
        score_t best_cell_score =
            last_diagonal_score +
            (sequence1[row - 1] == sequence2[col - 1] ? match : mismatch);
        // Determine best score from diagonal, vertical, or horizontal
        // direction.
        best_cell_score = std::max(best_cell_score, last_vertical_gap);
        best_cell_score = std::max(best_cell_score, horizontal_gap_column[row]);
        // Cache next diagonal value and store optimum in score_column.
        last_diagonal_score = score_column[row];
        score_column[row] = best_cell_score;
        // Compute the next values for vertical and horizontal gap.
        best_cell_score += gap_open;
        last_vertical_gap += gap_extension;
        horizontal_gap_column[row] += gap_extension;
        // Store optimum between gap open and gap extension.
        last_vertical_gap = std::max(last_vertical_gap, best_cell_score);
        horizontal_gap_column[row] =
            std::max(horizontal_gap_column[row], best_cell_score);
      }
    }

    // Report the best score.
    result[sequence_idx] = score_column.back();
  }

  return result;
}


#endif