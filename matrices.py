import numpy as np

def invert_matrix(matrix):
    """Invert the matrix and handle singular matrices"""
    try:
        inverse = np.linalg.inv(matrix)
        return inverse
    except np.linalg.LinAlgError:
        return None

def main():
    print("=== 6x6 Matrix Inverter ===\n")
    
    # Your matrix from Excel - replace these values as needed
    matrix = np.array([
        [1.76e+08, 4.60e+07, -6.5e+07, 0, 0, 0],
        [4.60e+07, 4.60e+07, 0, 0, 0, 0],
        [-6.5e+07, 0, 1.11e+08, 4.60e+07, 4.60e+07, -4.60e+07],
        [0, 0, 4.60e+07, 4.60e+07, -4.60e+07, 4.60e+07],
        [0, 0, 4.60e+07, -4.60e+07, 6.50e+07, 4.60e+07],
        [0, 0, -4.60e+07, 4.60e+07, 4.60e+07, 4.60e+07]
    ])
    
    print("Original Matrix:")
    print(matrix)
    print()
    
    # Calculate inverse
    inverse = invert_matrix(matrix)
    
    if inverse is not None:
        print("Inverse Matrix:")
        print(inverse)
        print()
        
        # Verify by multiplying
        print("Verification (A * A^-1, should be identity matrix):")
        verification = np.matmul(matrix, inverse)
        print(np.round(verification, decimals=8))
        print()
        
        # Print in a format easy to copy to Excel
        print("Inverse Matrix (formatted for Excel copy-paste):")
        for row in inverse:
            print('\t'.join([f'{val:.6e}' for val in row]))
    else:
        print("Error: Matrix is singular and cannot be inverted.")
        print("(The determinant is zero or very close to zero)")

if __name__ == "__main__":
    main()