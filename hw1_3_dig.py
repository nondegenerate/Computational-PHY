''' 
    第三题精确对角化代码(Python3)
'''
import numpy as np
from numba import jit

@jit(nopython=True)
def gen_H(N,V):
    h=1
    m=1
    dx=2/N
    c1=-h**2/(2*m*dx**2)
    c2=-2*c1
    H=np.zeros((N-1,N-1))
    for i in range(N-1):
        H[i,i]=c2+V[i]
    for i in range(N-2):
        H[i+1,i]=c1
        H[i,i+1]=c1
    return H

N=200
X=np.linspace(-1,1,N+1)
V=5*np.sign(X+0.5)-5*np.sign(X-0.5)
V=V[1:N]
H=gen_H(N,V)
result=np.linalg.eigh(H)
values=result[0]
vectors=result[1]
print(values[0])
np.savetxt('./hw1_3_dig.csv',vectors[:,0])