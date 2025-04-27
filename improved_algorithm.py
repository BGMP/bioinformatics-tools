"""
Optimized algorithm implementation for sequence alignment visualization
This module contains optimized versions of the alignment algorithms
"""

import numpy as np


def compute_needleman_wunsch(seq1, seq2, match_award, mismatch_penalty, gap_penalty):
    """
    Computes the entire Needleman-Wunsch matrix efficiently
    Returns the score matrix and the steps for visualization
    """
    n = len(seq1)
    m = len(seq2)

    score = np.zeros((m + 1, n + 1), dtype=int)

    # Initialize first row and column
    for i in range(m + 1):
        score[i, 0] = gap_penalty * i
    for j in range(n + 1):
        score[0, j] = gap_penalty * j

    # Record steps for animation
    steps = []

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            # Calculate match score
            if seq1[j - 1] == seq2[i - 1]:
                match = score[i - 1, j - 1] + match_award
            else:
                match = score[i - 1, j - 1] + mismatch_penalty

            delete = score[i - 1, j] + gap_penalty
            insert = score[i, j - 1] + gap_penalty

            # Determine the best score
            max_score = max(match, delete, insert)
            score[i, j] = max_score

            # Record which was the source of this score for traceback
            source = None
            if max_score == match:
                source = 'diagonal'
            elif max_score == delete:
                source = 'up'
            else:
                source = 'left'

            # Save this step for animation
            steps.append({
                'i': i,
                'j': j,
                'score': max_score,
                'match': match,
                'delete': delete,
                'insert': insert,
                'source': source
            })

    return score, steps


def compute_smith_waterman(seq1, seq2, match_award, mismatch_penalty, gap_penalty):
    """
    Computes the entire Smith-Waterman matrix efficiently
    Returns the score matrix, steps for visualization, and max position
    """
    n = len(seq1)
    m = len(seq2)

    # Initialize matrices
    score = np.zeros((m + 1, n + 1), dtype=int)

    # Record steps for animation
    steps = []

    # Variables to track maximum score
    max_score = 0
    max_i, max_j = 0, 0

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            # Calculate match score
            if seq1[j - 1] == seq2[i - 1]:
                match = score[i - 1, j - 1] + match_award
            else:
                match = score[i - 1, j - 1] + mismatch_penalty

            delete = score[i - 1, j] + gap_penalty
            insert = score[i, j - 1] + gap_penalty

            # Determine the best score (with 0 floor for local alignment)
            max_score_cell = max(0, match, delete, insert)
            score[i, j] = max_score_cell

            # Track global maximum for traceback starting point
            if max_score_cell > max_score:
                max_score = max_score_cell
                max_i, max_j = i, j

            # Record which was the source of this score for traceback
            source = None
            if max_score_cell == 0:
                source = 'zero'
            elif max_score_cell == match:
                source = 'diagonal'
            elif max_score_cell == delete:
                source = 'up'
            else:
                source = 'left'

            # Save this step for animation
            steps.append({
                'i': i,
                'j': j,
                'score': max_score_cell,
                'match': match,
                'delete': delete,
                'insert': insert,
                'source': source
            })

    return score, steps, (max_i, max_j)


def get_traceback_needleman_wunsch(seq1, seq2, score):
    """
    Generates the alignment and traceback path for Needleman-Wunsch
    """
    i = len(seq2)
    j = len(seq1)
    align1 = ""
    align2 = ""

    # Traceback path
    path = []

    while i > 0 or j > 0:
        path.append((i, j))
        current = score[i, j]

        if i > 0 and j > 0:
            diagonal = score[i - 1, j - 1]
            if seq1[j - 1] == seq2[i - 1]:
                diagonal_score = diagonal + 1  # Assuming match=1; adjust as needed
            else:
                diagonal_score = diagonal - 1  # Assuming mismatch=-1; adjust as needed

            # Check if diagonal move is optimal
            if current == diagonal_score:
                align1 = seq1[j - 1] + align1
                align2 = seq2[i - 1] + align2
                i -= 1
                j -= 1
                continue

        if j > 0:
            # Check if horizontal move is optimal
            left = score[i, j - 1]
            if current == left - 1:  # Assuming gap=-1; adjust as needed
                align1 = seq1[j - 1] + align1
                align2 = '-' + align2
                j -= 1
                continue

        # If we got here, vertical move must be optimal
        align1 = '-' + align1
        align2 = seq2[i - 1] + align2
        i -= 1

    # Return alignments and path in reverse (from start to end)
    return align1, align2, list(reversed(path))


def get_traceback_smith_waterman(seq1, seq2, score, start_pos):
    """
    Generates the alignment and traceback path for Smith-Waterman
    """
    i, j = start_pos
    align1 = ""
    align2 = ""

    # Traceback path
    path = []

    while i > 0 and j > 0 and score[i, j] > 0:
        path.append((i, j))
        current = score[i, j]

        # We stop when we hit a cell with score 0
        if current == 0:
            break

        diagonal = score[i - 1, j - 1]
        if seq1[j - 1] == seq2[i - 1]:
            diagonal_score = diagonal + 1  # Assuming match=1; adjust as needed
        else:
            diagonal_score = diagonal - 1  # Assuming mismatch=-1; adjust as needed

        # Check if diagonal move is optimal
        if current == diagonal_score:
            align1 = seq1[j - 1] + align1
            align2 = seq2[i - 1] + align2
            i -= 1
            j -= 1
            continue

        # Check if horizontal move is optimal
        left = score[i, j - 1]
        if current == left - 1:  # Assuming gap=-1; adjust as needed
            align1 = seq1[j - 1] + align1
            align2 = '-' + align2
            j -= 1
            continue

        # If we got here, vertical move must be optimal
        align1 = '-' + align1
        align2 = seq2[i - 1] + align2
        i -= 1

    # Return alignments and path in reverse (from start to end)
    return align1, align2, list(reversed(path))
