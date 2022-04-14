#!/usr/bin/python
__author__ = "Victoria Stevens"
__email__ = "victoria.stevens@yale.edu"
__copyright__ = "Copyright 2022"
__license__ = "GPL"
__version__ = "1.0.0"

#TODO: add edge cases with when you can't open the input files

### Usage: python smith_waterman.py -i <input file> -s <score file>
### Example: python smith_waterman.py -i input.txt -s blosum62.txt
### Note: Smith-Waterman Algorithm

import argparse
import numpy as np
from enum import IntEnum

# similarity matrix labels -- {letter:index}
sim_mat_labels = dict()

# enumerator for the traceback matrix
class Tracer(IntEnum):
    STOP = 0
    LEFT = 1 
    UP = 2
    DIAG = 3

### This is one way to read in arguments in Python. 
parser = argparse.ArgumentParser(description='Smith-Waterman Algorithm')
parser.add_argument('-i', '--input', help='input file', required=True)
parser.add_argument('-s', '--score', help='score file', required=True)
parser.add_argument('-o', '--opengap', help='open gap', required=False, default=-2)
parser.add_argument('-e', '--extgap', help='extension gap', required=False, default=-1)
parser.add_argument('-w', '--output', help='output file', required=False, default='output.txt')
args = parser.parse_args()

### Read the similarity matrix from the input BLOSUM file
### returns (23x23 if BLOSUM62) matrix, one index for each alphabetical letter
def simMat(scoreFile):
    with open(scoreFile) as f:
        # count number of columns
        letters = f.readline().split()
        n_cols = len(letters) + 1

        # populate sim_mat_labels dict with key-value pairs = letter-index
        for i in range(n_cols - 1):
            sim_mat_labels[letters[i]] = i

    # print(sim_mat_labels)

    # define which columns will be read into matrix
    # want to get rid of label column (col 0)
    usecols = np.arange(1, n_cols)

    # read in BLOSUM text as matrix
    input = np.loadtxt(scoreFile, dtype="int", skiprows=1, usecols=usecols)
    return input

# Read the two input sequences from input file
# return these two sequences as strings
def sequences(inputFile):
    with open(inputFile) as f:
        sequence1 = f.readline()
        sequence2 = f.readline()
    return sequence1, sequence2

# Initialize score matrix
# returns (len(sequence2) + 2) x (len(sequence1) + 2) matrix, AND traceback matrix with Nones
def scoreMat(sequence1, sequence2):
    score_matrix = list() # will be a list of lists, to form score matrix
    traceback_matrix = list()

    # format first row = sequence1
    # print([None, 'A', 4])
    score_matrix.append([None, None])
    traceback_matrix.append([None, None])
    seq1_len = len(sequence1)
    for letter in range(seq1_len - 1):
        # print(letter)
        score_matrix[0].append(sequence1[letter])
        traceback_matrix[0].append(sequence1[letter])

    # format first column
    score_matrix.append([None])
    traceback_matrix.append([None])
    for placeholder in range(seq1_len):
            score_matrix[1].append(0)
            traceback_matrix[1].append([0, None])

    for letter in sequence2:
        first_letter = list(letter)
        tb_letter = list(letter)
        # initialize all rows with zeros
        for placeholder in range(seq1_len):
            first_letter.append(0)
            tb_letter.append([0, None])
        score_matrix.append(first_letter)
        traceback_matrix.append(tb_letter)

    return score_matrix, traceback_matrix

def gap(steps, openGap, extGap):
    # from wikipedia:
    # W=uk+v, where v is the gap opening penalty, 
    # and u is the gap extension penalty.
    return (extGap * (steps - 1)) + openGap

def steps(score_list):
    maximum = max(score_list)
    return [i for i, j in enumerate(score_list) if j == maximum][0]

# returns score at that index in score_matrix and correct input at that index in the traceback_matrix
def score(score_matrix, traceback_matrix, similarity_matrix, openGap, extGap, curr_row, curr_col):
    sim_index_seq1 = sim_mat_labels[score_matrix[0][curr_col]]
    sim_index_seq2 = sim_mat_labels[score_matrix[curr_row][0]]
    # print(score_matrix[0][curr_col] + score_matrix[curr_row][0] + str(similarity_matrix[sim_index_seq2][sim_index_seq1]))
    match_score = similarity_matrix[sim_index_seq2][sim_index_seq1]

    row_scores = [(score_matrix[curr_row][curr_col-left] + gap(left, openGap, extGap)) for left in range(1,curr_col)]
    column_scores = [(score_matrix[curr_row-up][curr_col] + gap(up, openGap, extGap)) for up in range(1,curr_row)]

    # print(score_matrix[0][curr_col] + score_matrix[curr_row][0])
    # print(str(row_scores))
    # print(str(column_scores))
    
    left = max(row_scores)
    up = max(column_scores)

    diag = score_matrix[curr_row - 1][curr_col - 1] + match_score

    total_step_score =  max(0, left, up, diag) 

    if total_step_score == diag:
        tracer = [1, Tracer.DIAG]
    elif total_step_score == left: 
        tracer = [steps(row_scores) + 1, Tracer.LEFT]
    elif total_step_score == up:
        tracer = [steps(column_scores) + 1, Tracer.UP]
    else:
        tracer = [0, Tracer.STOP]

    score_matrix[curr_row][curr_col] = total_step_score
    traceback_matrix[curr_row][curr_col] = tracer

    return score_matrix, traceback_matrix

# update the score matrix using the Smith-Waterman algorithm
# return updated score matrix
def smithWaterman(score_matrix, traceback_matrix, similarity_matrix, openGap, extGap):
    num_rows = len(score_matrix)
    num_cols = len(score_matrix[0])

    # tracing_matrix = np.zeros(shape=(num_rows, num_cols), dtype=int)  

    alignment_score = 0
    tb_start = (0, 0)

    for i in range(2, num_rows):
        for j in range(2, num_cols):
            score_matrix, traceback_matrix = score(score_matrix, traceback_matrix, similarity_matrix, openGap, extGap, i, j)
            if score_matrix[i][j] > alignment_score:                  
                alignment_score = score_matrix[i][j]
                tb_start = (i, j)

    return score_matrix, traceback_matrix, alignment_score, tb_start

def traceback(traceback_matrix, tb_start):
    row, col = tb_start
    seq1 = []
    seq2 = []

    end_seq2, end_seq1 = tb_start
    end_seq1 -= 2
    end_seq2 -= 2

    start_seq1 = end_seq1
    start_seq2 = end_seq2

    while traceback_matrix[row][col] != [0, None] and traceback_matrix[row][col] != [0, Tracer.STOP]:
        start_seq1 = col - 2
        start_seq2 = row - 2

        steps, direc = traceback_matrix[row][col][0], traceback_matrix[row][col][1]
        if direc == Tracer.DIAG:
            seq1.append(traceback_matrix[0][col])     
            seq2.append(traceback_matrix[row][0])
            row -= 1
            col -= 1
            # print(str(row) + ','+ str(col)+ ',' + str(traceback_matrix[row][col]))
        elif direc == Tracer.LEFT:
            for step in range(steps):
                seq1.append(traceback_matrix[0][col])
                seq2.append('-')
                col -= 1
                # print(str(row) + ','+ str(col)+ ',' + str(traceback_matrix[row][col]))
        elif direc == Tracer.UP: 
            for step in range(steps):
                seq1.append('-')            
                seq2.append(traceback_matrix[row][0])            
                row -= 1
                # print(str(row) + ','+ str(col)+ ',' + str(traceback_matrix[row][col]))

    reverse_seq1 = reversed(seq1)
    reverse_seq2 = reversed(seq2)
    sequence1_aligned = "".join(reverse_seq1)
    sequence2_aligned = "".join(reverse_seq2)

    alignment = []
    for i in range(len(sequence1_aligned)):
        if sequence1_aligned[i] == sequence2_aligned[i]:
            alignment.append('|')
        else:
            alignment.append(' ')
        
    return sequence1_aligned, sequence2_aligned, "".join(alignment), start_seq1, start_seq2, end_seq1, end_seq2

def alignFormat(sequence1, sequence2, seq1_traced, seq2_traced, alignment, start_seq1, start_seq2, end_seq1, end_seq2):
    padding = " " * (max(start_seq1, start_seq2))

    align = padding + ' ' + alignment
    if start_seq1 > start_seq2:
        seq1 = sequence1[0:start_seq1] + '(' + seq1_traced + ')' + sequence1[end_seq1 + 1:]
        seq2 = (' ' * (start_seq1 - start_seq2)) + sequence2[0:start_seq2] + '(' + seq2_traced + ')' + sequence2[end_seq2 + 1:]
    elif start_seq2 > start_seq1:
        seq1 = (' ' * (start_seq2 - start_seq1)) + sequence1[0:start_seq1] + '(' + seq1_traced + ')' + sequence1[end_seq1 + 1:]
        seq2 = sequence2[0:start_seq2] + '(' + seq2_traced + ')' + sequence2[end_seq2 + 1:]
    else:
        seq1 = sequence1[0:start_seq1] + '(' + seq1_traced + ')' + sequence1[end_seq1 + 1:]
        seq2 = sequence2[0:start_seq2] + '(' + seq2_traced + ')' + sequence2[end_seq2 + 1:]

    # newline stripping and left justifying with spaces
    seq1 = seq1.replace("\n", "")
    seq2 = seq2.replace("\n", "")
    align = align.replace("\n", "")
    longest_str = max(len(seq1), len(seq2), len(align))
    return seq1.ljust(longest_str), seq2.ljust(longest_str), align.ljust(longest_str)


### Implement your Smith-Waterman Algorithm
def runSW(inputFile, scoreFile, openGap, extGap, outputFile):
    ### Input file processing

    # save the similarity matrix
    similarity_matrix = simMat(scoreFile)
    # print(similarity_matrix)

    # read in the two sequences to be compared
    sequence1, sequence2 = sequences(inputFile)

    ### calculation

    # initialize score matrix with 0s
    score_matrix, traceback_matrix = scoreMat(sequence1, sequence2)
    # print(score_matrix)
    # print(traceback_matrix)

    # and calculate score using Smith-Waterman algorithm
    score_matrix, traceback_matrix, alignment_score, tb_start = smithWaterman(score_matrix, traceback_matrix, similarity_matrix, openGap, extGap)
    # print(traceback_matrix)
    # print(tb_start)

    # use the traceback matrix and starting index to construct aligned strings
    seq1_traced, seq2_traced, alignment, start_seq1, start_seq2, end_seq1, end_seq2 = traceback(traceback_matrix, tb_start)
    # print(seq1_traced)
    # print(alignment)
    # print(seq2_traced)
    # print(start_seq1)
    # print(start_seq2)
    # print(end_seq1)
    # print(end_seq2)

    final_seq1, final_seq2, final_alignment = alignFormat(sequence1, sequence2, seq1_traced, seq2_traced, alignment, start_seq1, start_seq2, end_seq1, end_seq2)
    print(final_seq1)
    print(final_alignment)
    print(final_seq2)

    ### write output
    with open(outputFile, 'w+') as f:
        f.write("""-----------
|Sequences|
-----------\n""")
        f.write("sequence1\n")
        f.write(sequence1)
        f.write("sequence2\n")
        f.write(sequence2)
        f.write("""\n--------------
|Score Matrix|
--------------\n""")

        num_rows = len(score_matrix)
        num_cols = len(score_matrix[0])
        for i in range(num_rows):
            for j in range(num_cols):
                if score_matrix[i][j] == None:
                    f.write('\t')
                # elif score_matrix[i][j] == '\n':
                #     f.write(str(score_matrix[i][j]))
                else:
                    f.write(str(score_matrix[i][j]) + '\t')
            f.write('\n')

        f.write("""----------------------
|Best Local Alignment|
----------------------\n""")
        f.write('Alignment Score:' + str(alignment_score) + '\n')
        f.write('Alignment Results:\n')
        f.write(final_seq1 + '\n')
        f.write(final_alignment + '\n')
        f.write(final_seq2 + '\n')


    # print(similarity_matrix)
    # print(sim_mat_labels)

    # print(sequence1)
    # print(sequence2)

### Run the Smith-Waterman Algorithm
runSW(args.input, args.score, args.opengap, args.extgap, args.output)