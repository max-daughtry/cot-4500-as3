import numpy as np

# Euler's Approximation method
def Question1():
    # Differential equation
    dy = lambda t, y: t-y**2

    # Initial value
    y0 = 1

    # Bounds
    a = 0
    b = 2
    
    # Number of steps
    N = 10

    # Step size
    h = (b-a)/N

    # Initialize variables
    t = a
    w = y0

    # Euler's Approximation
    for i in range(1, N+1):
        w = w + h * dy(t,w)
        t = a + i * h
    
    # Print result
    print(w)

# Runge-Kutta Method
def Question2():
    # Differential equation
    dy = lambda t, y: t-y**2

    # Initial value
    y0 = 1

    # Bounds
    a = 0
    b = 2

    # Number of steps
    N = 10

    # Step size
    h = (b-a)/N

    # Initialize variables
    t = a
    w = y0

    # Runge-Kutta Method to order 4
    for i in range(1, N+1):
        K1 = h * dy(t, w)
        K2 = h * dy(t + h/2, w + K1/2)
        K3 = h * dy(t + h/2, w + K2/2)
        K4 = h * dy(t + h, w + K3)

        w = w + (K1 + 2 * K2 + 2 * K3 + K4) / 6
        t = a + i * h

    # print result
    print(w)

# Solving system of equations
def Question3():
    # Dimension of matrix A
    n = 3

    # Augmented matrix A
    A = np.array([[2, -1, 1, 6], [1, 3, 1, 0], [-1, 5, 4, -3]], dtype=float)

    # Gaussian elimination
    for i in range(1, n):
        p = i

        if not (i <= p <= n and A[p-1][i-1] != 0):
            print("No unique solution exists...")
            return

        if (p != i):
            swap = A[i-1]
            A[i-1] = A[p-1]
            A[p-1] = swap
        
        for j in range(i+1, n+1):
            m = A[j-1][i-1] / A[i-1][i-1]
            A[j-1] = A[j-1] - m * A[i-1]
        
    if A[n-1][n-1] == 0:
        print("No unique solution exists...")
        return
    
    # Back substitution to find solution x

    x = np.zeros(n)

    x[n-1] = A[n-1][n] / A[n-1][n-1]

    for i in range(n-1, 0, -1):
        sumTerm = 0
        for j in range(i+1, n+1):
            sumTerm += A[i-1][j-1]*x[j-1]
        x[i-1] = (A[i-1][n] - sumTerm)/A[i-1][i-1]
    
    print(x)

# LU Factorization
def Question4():
    # Dimension of matrix A
    n = 4

    # Matrix A
    A = np.array([[1,1,0,3],[2,1,-1,1],[3,-1,-1,2],[-1,2,3,-1]],dtype=float)

    # Check candidacy of A
    if A[0][0] == 0:
        print("Factorization impossible...")
        return

    # Create lower diagonal matrix with 1's along diagonal
    L = np.array(np.zeros((n,n)))
    for i in range(len(L)):
        L[i][i] = 1
    
    # Create upper diagonal matrix
    U = np.array(np.zeros((n,n)))

    # Choose U[0][0] such that L[0][0] * U[0][0] = A[0][0]
    U[0][0] = A[0][0] / L[0][0]

    # Find first row and first column of U and L respectively
    for j in range(2, n+1):
        U[0][j-1] = A[0][j-1] / L[0][0]
        L[j-1][0] = A[j-1][0] / U[0][0]

    # LU factorization
    for i in range(2, n):
        sumTerm = 0
        for k in range(1, i):
            sumTerm += L[i-1][k-1] * U[k-1][i-1]

        U[i-1][i-1] = (A[i-1][i-1] - sumTerm) / L[i-1][i-1]

        if L[i-1][i-1] * U[i-1][i-1] == 0:
            print("Factorization impossible...")
        
        for j in range(i+1, n+1):
            sumTerm = 0
            for k in range(1, i):
                sumTerm += L[i-1][k-1] * U[k-1][j-1]
            U[i-1][j-1] = (1/L[i-1][i-1]) * (A[i-1][j-1] - sumTerm)

            sumTerm = 0
            for k in range(1, i):
                sumTerm += L[j-1][k-1] * U[k-1][i-1]
            L[j-1][i-1] = (1/U[i-1][i-1]) * (A[j-1][i-1] - sumTerm)

    sumTerm = 0
    for k in range(1, n):
        sumTerm += L[n-1][k-1] * U[k-1][n-1]
    U[n-1][n-1] = (A[n-1][n-1] - sumTerm) / L[n-1][n-1]
    
    # Find determinant of L, U and then A
    detL = 1
    detU = 1
    for i in range(n):
        detL *= L[i][i]
        detU *= U[i][i]
    detA = detL * detU
    
    print(detA)
    print()
    print(L)
    print()
    print(U)

# Diagonally Dominant
def Question5():
    # Matrix A
    A = np.array([[9,0,5,2,1],[3,9,1,2,1],[0,1,7,2,3],[4,2,3,12,2],[3,2,4,0,8]], dtype=float)

    diagDom = True

    # By definition of diagonally dominant
    for i in range(len(A)):
        sumTerm = 0
        for j in range(len(A)):
            if j==i:
                continue
            sumTerm += np.abs(A[i][j])
        
        if np.abs(A[i][i]) < sumTerm:
            diagDom = False
            break
    
    print(diagDom)
    
# Positive Definite
def Question6():
    # Matrix A
    A = np.array([[2,2,1],[2,3,0],[1,0,2]], dtype=float)

    posDef = True

    # Check A is symmetric
    if not np.array_equal(A, np.transpose(A)):
        posDef = False
        print(posDef)
        return

    # If A is symmetric and all leading principal submatrices have positive determinant, then A is positive definite
    # This is checked by contradiction (checking if any principal submatrix determinants are negative)
    for i in range(1,len(A)+1):
        if np.linalg.det(A[0:i, 0:i]) <= 0:
            posDef = False
            break

    print(posDef)

if __name__ == '__main__':
    Question1()
    print()
    Question2()
    print()
    Question3()
    print()
    Question4()
    print()
    Question5()
    print()
    Question6()