# -*- coding: utf-8 -*-
"""
Created on Thu Feb 23 22:38:27 2017

@author: Raihan
"""
import numpy as np
import matplotlib.pyplot as plt
#This script implements the conjugate gradient method for different dimensions of A
def hilbert(n):
    A=np.zeros((n,n))
    
    for i in range(0,n):
        for j in range(0,n):
            
           A[i,j]=(1/(i+1+j+1-1))
    return (A)    

def CG(dimension,A,b,x0):  
        
    
    
    xk=x0
    r0=np.matmul(A,x0)-b
    #print(r0)
    rk=r0
    
    p0=-r0
    pk=p0
    iterations=0
    rkall=[]
    while np.linalg.norm(rk)>10**-6:
        ak=(np.matmul((np.transpose(rk)),rk))/(np.matmul(np.matmul((np.transpose(pk)),A),pk))
        #print(ak)
        xk=xk+ak*pk
        #print(ak*pk)
        rknext=rk+np.matmul(ak*A,pk)
        bk=(np.matmul(np.transpose(rknext),rknext))/(np.matmul(np.transpose(rk),rk))
        pk=-rknext+bk*pk
        rkall.append(np.linalg.norm(rk))
        rk=rknext
        iterations=iterations+1
    print('The value converges after %d iterations',iterations)
    #print(xk)
    #print(np.linalg.norm(rk))
        
    return(xk,rk,rkall,iterations)
def main(n):
    #Performing the main code here 
    
    A=hilbert(n)
    b=np.ones((n,1))
    x0=np.zeros((n,1))
    xk,rk,rkall,iterations= CG(n,A,b,x0)
    print(xk)
    print(rk)
    plt.plot(np.arange(iterations),rkall)
    plt.ylabel('norm of rk')
    plt.xlabel('number of iterations')
    plt.show()
    

main(5)