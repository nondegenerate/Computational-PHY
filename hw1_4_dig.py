''' 
    第四题精确对角化代码(Python3)
'''
import numpy as np
from numba import jit

@jit(nopython=True)
def gen_H(N,h,w0,x0,m):
    a=(m*w0/h)**0.5
    H=np.zeros((N,N))
    ksai=np.zeros((N,N))
    dig=0.5*h*w0+m*w0**2*x0**2/8
    for i in range(N):
        H[i,i]=i*h*w0+dig
        if i>0:
            ksai[i-1,i]=(i/2)**0.5
        if i<N-1:
            ksai[i+1,i]=((i+1)/2)**0.5
    ksai2=np.dot(ksai,ksai)
    ksai4=np.dot(ksai2,ksai2)
    k=m*w0**2/(8*x0**2)
    H=H+k*(1/a**4)*ksai4-k*(10*x0**2/a**2)*ksai2
    return H

N=200
H=gen_H(N,1,1,1,1)
result=np.linalg.eigh(H)
values=result[0]
vectors=result[1]
print(values[0])
np.savetxt('./hw1_4_dig.csv',vectors[:,0])