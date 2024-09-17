import numpy as np 
import sys 

n = np.array([[502.78458, 988.2237683620597, 1004.0459601286647],
                [0, 655.62744, 1369.662818937229], 
                [0, 0, 767.80225]])

def getValue(k): 
    ks = np.array([(1-k**2), 0.5*k*(1+k), 0.5*k*(k-1)])
    intermediate = np.matmul(n, ks)
    print(intermediate)
    return np.matmul(intermediate, ks)


value = float(sys.argv[1])

val = getValue(value)
print(val, np.sqrt(val))