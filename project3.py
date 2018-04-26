# -*- coding: utf-8 -*-
"""
Created on Sun Apr 24 19:22:21 2016

@author: covingtonjg

Jessica Covington
"""


import numpy as np

def jacobi(A, b, x_init, kmax, tol):
    n = len(A)
    y = np.zeros(n)
    y[:] = x_init[:]
    x = np.zeros(n)

    err = tol    
    
    for k in range(kmax):
        for i in range(n):
            total = b[i]
            diagonal = A[i,i]
            if (abs(diagonal) < tol):
                print("diagonal too small")
            for j in range(n):
                if (j != i):
                    total -= A[i,j]*y[j]
            x[i] = total / diagonal
            
        if (np.linalg.norm(x-y,np.inf) < err):
            return k, x
        y[:] = x[:]
    return k, x
                                
    
def gaussSeidel(A, b, x_init, kmax, tol):
    n = len(A)
    y = np.zeros(n)
    y[:] = x_init[:]
    x = np.zeros(n)
    
    err = tol 
    
    for k in range(kmax):
        for i in range(n):
            total = b[i]
            diagonal = A[i,i]
            if (abs(diagonal) < tol):
                print("diagonal too small")
            for j in range(i):
                total -= A[i,j]*x[j]
            for j in range(i+1, n):
                total -= A[i,j]*y[j]
            x[i] = total / diagonal
            
        if (np.linalg.norm(x-y,np.inf) < err):
            return k, x
        y[:] = x[:]
    return k, x
    
        
def SOR(A, b, x_init, w, kmax, tol):
    n = len(A)
    y = np.zeros(n)
    y[:] = x_init[:]
    x = np.zeros(n)

    err = tol    
    
    for k in range(kmax):
        for i in range(n):
            total = b[i]
            diagonal = A[i,i]
            if (abs(diagonal) < tol):
                print("diagonal too small")
            for j in range(i):
                total -= A[i,j]*x[j]
            for j in range(i+1, n):
                total -= A[i,j]*x[j]
            x[i] = total/diagonal
            x[i] = w*x[i] + (1-w)*y[i]
            
        if (np.linalg.norm(x-y,np.inf) < err):
            return k, x
        y[:] = x[:]
    return k, x
