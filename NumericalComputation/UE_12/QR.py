import numpy as np
import scipy.linalg as la
import matplotlib.pyplot as plt  
import sys
sys.setrecursionlimit(1500)

eps = 1e-4

def simple_QR(A0, lmax):

    A_list = []
    A_list.append(A0)

    for i in range(lmax):
        Q,R = la.qr(A0)
        A0 = np.dot(R,Q)
        A_list.append(A0)
    return A_list

def calc_shift(A):
    if A.shape[0] > 1:
        A = A[:-2, :-2]

    n = A.shape[0]-1

    eig, v = la.eig(A)
    eig = eig - A[n,n]

    return np.amin(np.abs(eig))

def fancy_QR(A):
    n = A.shape[0]-1
    qr_count = 0
    
    while(abs(A[n-1,n]) > eps*(abs(A[n-1,n]) + abs(A[n,n]))):
        s = calc_shift(A)
        Q,R = la.qr(A - s*np.eye(n+1))
        A = np.dot(R,Q)  + s*np.eye(n+1)
        qr_count += 1

    if n > 2:
        qr_rec, eig_rec = fancy_QR(A[:-1,:-1])
        return qr_count+qr_rec, [A[n,n]] + eig_rec
    else:
        return qr_count, list(np.diagonal(A))
        

A0 = np.array([[1,0,0], [0,2,0], [2,1,3]])
A0 = la.hessenberg(A0)
lmax = 3

#print(la.eig(A0))
#print(np.diagonal(simple_QR(A0, lmax)))
#print(fancy_QR(A0))

def A_triag(n):
    return np.diag(2*np.ones(n)) - np.diag(np.ones(n-1),-1) - np.diag(np.ones(n-1), 1)

A1 = A_triag(2**3)
A2 = A_triag(2**5)

A1_list = simple_QR(A1, 20)
A2_list = simple_QR(A2, 100)

eig_gt_1 = np.sort(la.eig(A1)[0])
error_eig_1 = []

eig_gt_2 = np.sort(la.eig(A2)[0])
error_eig_2 = []

for A in A1_list:
    eig_approx = np.sort(np.diagonal(A))
    error_eig_1.append(np.amax(eig_approx - eig_gt_1))

for A in A2_list:
    eig_approx = np.sort(np.diagonal(A))
    error_eig_2.append(np.amax(eig_approx - eig_gt_2))

plt.subplot(1,2,1)
plt.plot(range(1,len(error_eig_1)+1), error_eig_1)
plt.title('n = 2^3')
plt.subplot(1,2,2)
plt.plot(range(1,len(error_eig_2)+1), error_eig_2)
plt.title('n = 2^5')
plt.show()

#----------------------------------------------------------------------------

qr_count_list = []
n_list = []
for i in range(2,9):
    n_list.append(2**i)

for n in n_list:

    A = la.hessenberg(A_triag(n))
    qr_count, eig = fancy_QR(A)
    print(qr_count)
    qr_count_list.append(qr_count)
    
    
plt.loglog(n_list, qr_count_list)
plt.xlabel('size')
plt.ylabel('QR-steps')
plt.show()
