# Goal of the program: To create a simple sequence alignment tool that can align two sequences using global alignment.
# Initially we will use the Needleman-Wunsch algorithm for global alignment.
# For development purposes, we will create a pair of toy sequences to align.



sequence_1 = "ACGTAGCTAG"
sequence_2 = "ACGTTAGCTAG"



# Create a formula for scoring the cell
def score_cell(seq1_char: str, seq2_char: str) -> int:
    # match and mismatch scoring
    if seq1_char == seq2_char:
        return 1  # match score
    else:
        return -1  # mismatch score
    # Deletion and insertion will be worth -1, but we will handle that in the grid creation


# create a grid structure for the alignment, list of lists
def create_grid(seq1: str, seq2: str) -> list:
    # define the cost of gaps, insertions, and deletions
    gap_penalty = -1

    rows = len(seq1) + 1
    cols = len(seq2) + 1
    grid = [[0 for _ in range(cols)] for _ in range(rows)]

    # Corrected Initialization
    # Initialize the first column (aligning seq1 with gaps)
    for i in range(1, rows):
        grid[i][0] = i * gap_penalty

    # Initialize the first row (aligning seq2 with gaps)
    for j in range(1, cols):
        grid[0][j] = j * gap_penalty

    # Fill the grid with scores
    for i in range(1, rows):
        for j in range(1, cols):
            match = grid[i-1][j-1] + score_cell(seq1[i-1], seq2[j-1])
            delete = grid[i-1][j] + gap_penalty # Deletion penalty
            insert = grid[i][j-1] + gap_penalty # Insertion penalty
            grid[i][j] = max(match, delete, insert)
    # Return the filled grid
    return grid


# function to print the grid for debugging purposes, with formatting and headers
def print_grid(grid: list, seq1: str, seq2: str):
    # Print header row (sequence 2)
    print("      ", end="")
    for char in seq2:
        print(f"{char:>3}", end="")
    print()

    for i, row in enumerate(grid):
        # Print header column (sequence 1)
        if i == 0:
            print("  ", end="")
        else:
            print(f"{seq1[i-1]} ", end="")
        
        # Print grid cells
        for cell in row:
            print(f"{cell:>3}", end="")
        print()

# traceback function to reconstruct the alignment
def traceback(grid: list, seq1: str, seq2:str) -> tuple:
    # Retrive scoring parameters from the grid logic
    gap_penalty = -1

    aligned_seq1 = ""
    aligned_seq2 = ""
    i, j =  len(seq1), len(seq2)

    while i > 0 or j > 0:
        # calculate the score of a potential match/mismatch
        current_score = grid[i][j]
        diagonal_score = grid[i-1][j-1]
        score_from_diag = diagonal_score + score_cell(seq1[i-1], seq2[j-1])

        # Check if the path came from diagonal (match/mismatch)
        if current_score == score_from_diag:
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            i -= 1
            j -= 1
        # Check if the path came from above (deletion)
        elif current_score == grid[i-1][j] + gap_penalty:
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = "-" + aligned_seq2
            i -= 1
        # Check if the path came from the left (insertion)
        else:
            aligned_seq1 = "-" + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            j -= 1

    # Handle any remaining characters if one sequence is longer than the other
    while i > 0:
        aligned_seq1 = seq1[i-1] + aligned_seq1
        aligned_seq2 = "-" + aligned_seq2
        i -= 1
    while j > 0:
        aligned_seq1 = "-" + aligned_seq1
        aligned_seq2 = seq2[j-1] + aligned_seq2
        j -= 1
        
    return aligned_seq1, aligned_seq2


def main():
    # 1. Create the scoring grid
    grid = create_grid(sequence_1, sequence_2)
    
    # 2. Perform the traceback to get the alignment
    aligned_seq1, aligned_seq2 = traceback(grid, sequence_1, sequence_2)
    
    # 3. Print the results
    print("Alignment Grid:")
    print_grid(grid, sequence_1, sequence_2)
    
    alignment_score = grid[len(sequence_1)][len(sequence_2)]
    print(f"\nOptimal Alignment Score: {alignment_score}")
    
    print("\nAligned Sequences:")
    print(aligned_seq1)
    print(aligned_seq2)

if __name__ == "__main__":
    main() 
