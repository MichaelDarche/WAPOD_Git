# sumUs_c.pyx

# Import necessary modules
import numpy as np
cimport numpy as np

def compute_sumUs(np.ndarray[double, ndim=2] U, np.ndarray[double, ndim=2] Umod, np.ndarray[double, ndim=2] Un, np.ndarray[double, ndim=3] csTot, int xstep, int K):
    cdef np.ndarray[double, ndim=2] cs
    cdef np.ndarray[double, ndim=2] Us
    cdef int s
    
    Us = np.zeros((2, 1))

    for s in range(K+1):
        cs = csTot[s, xstep].astype(np.double)
        sumUs = np.dot(cs, Umod[:, xstep + s - int(K/2)])
        Us[0, 0] += sumUs[0]
        Us[1, 0] += sumUs[1]

    Un[0, xstep] = U[0, xstep] - Us[0, 0]
    Un[1, xstep] = U[1, xstep] - Us[1, 0]
