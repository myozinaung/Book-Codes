# Transfer Function to State Space
def tf2ss(numerator,denominator):
    n = denominator.len-1 # Order
    A = np.ndarray((n,n)
    B = np.ndarray((n,1))
    C = np.ndarray((1,n))
    # Add zero(s) in-front if numerator order < denominator
    for row in range(n):
        if row == 0:
            for col in range(n):
                A[row,col] = denominator[n-col-1]
            B[row] = 1
        else: 
            A[]
            B[row] = 0   
        C = 
    D = np.array([[numerator[0]]])
    
    return A, B, C, D    

nume = [1 2 5] # Numerator Coefficients
deno = [1 4 5] # Denominator Coefficients
# Note: Order of Denominator must be greater than or equal to order of numerator
A, B, C, D = tf2ss(nume,deno)
print(A,'\n',B,'\n',C,'\n',D)