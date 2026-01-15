import numpy as np

def invert_matrix(matrix):
    """Invert the matrix and handle singular matrices"""
    try:
        inverse = np.linalg.inv(matrix)
        return inverse
    except np.linalg.LinAlgError:
        return None

def multiply_by_vector(inverse):
    """Multiply the inverse matrix by a 6x1 vector"""
    print("\n" + "="*50)
    print("=== Matrix-Vector Multiplication ===\n")
    
    # Your 6x1 vector - replace these values
    vector = np.array([
        [0],  # Element 1
        [0],  # Element 2
        [0],  # Element 3
        [-5000],  # Element 4
        [0],  # Element 5
        [0]   # Element 6
    ])
    
    print("Vector:")
    print(vector.flatten())
    print()
    
    # Multiply
    result = np.matmul(inverse, vector)
    
    print("Result (Inverse Matrix × Vector):")
    print(result.flatten())
    print()
    
    # Print in format for Excel
    print("Result (formatted for Excel):")
    for val in result.flatten():
        print(f'{val:.6e}')

def main():
    print("=== 6x6 Matrix Inverter ===\n")
    
    # Your matrix from Excel - replace these values as needed
    matrix = np.array([
        [1.53e+08, 2.30e+07, -65030968, 0, 0, 0],
        [2.30e+07, 2.30e+07, 0, 0, 0, 0],
        [-6.5e+07, 0, 8.8e+07, 2.30e+07, 2.30e+07, -2.30e+07],
        [0, 0, 2.30e+07, 2.30e+07, -2.30e+07, 2.30e+07],
        [0, 0, 2.30e+07, -2.30e+07, 6.50e+07, 2.30e+07],
        [0, 0, -2.30e+07, 2.30e+07, 2.30e+07, 2.30e+07]
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
        
        # Now multiply by vector
        multiply_by_vector(inverse)
    else:
        print("Error: Matrix is singular and cannot be inverted.")
        print("(The determinant is zero or very close to zero)")

if __name__ == "__main__":
    main()