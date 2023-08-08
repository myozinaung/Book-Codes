# Discrete-Time Linear Quadratic Regulator (LQR)
def DLQR_Gain2(A, B, Q, R): # Solve DARE without Scipy
    X       = Q # Initialize X as the same size as Q
    maxiter = 1000 # Intertively solve DARE
    eps     = 1e-6 # Stopping Criteria

    for i in range(maxiter):
        Xn = A.T @ X @ A - A.T @ X @ B @ \
            LA.inv(R + B.T @ X @ B) @ B.T @ X @ A + Q
        if (abs(Xn - X)).max() < eps:
            break
        X = Xn
    K = LA.inv(B.T @ X @ B + R) @ (B.T @ X @ A)
    eigVals, eigVecs = LA.eig(A-B@K)
    print("LQR Gains:-K",-K)
    return K

# Discrete-Time Linear Quadratic Regulator (LQR)
def DLQR_Gain(A, B, Q, R):
    X = np.matrix(LA.solve_discrete_are(A, B, Q, R))
    K = np.matrix(LA.inv(B.T*X*B+R)*(B.T*X*A))
    eigVals, eigVecs = LA.eig(A-B*K)
    print("LQR Gains:-K",-K)
    return K