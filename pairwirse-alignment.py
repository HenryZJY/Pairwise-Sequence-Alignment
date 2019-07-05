#################	PROVIDED CODE   #################

import math
import random
from collections import deque
from collections import defaultdict



#helper function
def ComputeGlobalScore(X,Y,M):
    """
    Computes the global optimal score matrix for dynamic programming
    
    Input: Strings X and Y, and scoring matrix M
    Returns: The optimal score matrix S
    """
    S = []
    lenX = len(X)
    lenY = len(Y)
    
    #Initialize the matrix
    for i in range(0,lenX+1):
        S.append([])
        for j in range(0,lenY+1):
            S[i].append(0)
    #print S   
    
    #Starts the DP
    for i in range(0, lenX+1):
        
        for j in range(0, lenY+1):
            
            if i==0 and j>0:
                S[i][j] = S[i][j-1] + M['-'][Y[j-1]]
            
            elif i>0 and j==0:
                S[i][j] = S[i-1][j] + M[X[i-1]]['-']
            
            elif i>0 and j>0:
                S[i][j] = max(S[i-1][j-1] + M[X[i-1]][Y[j-1]],
                             S[i][j-1] + M['-'][Y[j-1]],
                             S[i-1][j] + M[X[i-1]]['-'])
            
    return S

#test 
#X1 = "ACC"
#Y1 = "TAGTA"
#M1 = {'-':{'-':0, 'T': -4, 'G':-4, 'A':-4}, 'A':{'-':-2, 'T': 2, 'G':2, 'A':5}, 'C':{'-':-2, 'T': 2, 'G':2, 'A':2}}
#Expected :[[0, -4, -8, -12, -16, -20], [-2, 2, 1, -3, -7, -11], [-4, 0, 4, 3, -1, -5], [-6, -2, 2, 6, 5, 1]]
#print ComputeGlobalScore(X1,Y1,M)

X3 = ""
#print ComputeGlobalScore(X3,Y1,M)

#helper function
def ComputeAlignment(X,Y,M,S):
    """
    Trace back the alignment through the DP's  global optimal score matrix
    
    Inputs: Strings X and Y, scoring Matrix M, and the optimal score matrix
    Returns: A tuple of two strings
    """
    Xp = ''
    Yp = ''
    i = len(X)
    j = len(Y)
    while i>0 and j>0 :
        #print "i",i,"j",j
        if S[i][j] == S[i-1][j-1] + M[X[i-1]][Y[j-1]]:
            Xp = X[i-1] + Xp
            Yp = Y[j-1] + Yp
            i -= 1
            j -= 1
        elif S[i][j] == S[i][j-1] + M['-'][Y[j-1]]:
            Xp = '-' + Xp
            Yp = Y[j-1] + Yp
            j -= 1
        elif S[i][j] == S[i-1][j] + M[X[i-1]]['-']:
            Xp = X[i-1] + Xp
            Yp = '-' + Yp
            i -= 1

    if i==0 and j>0:
        Xp = '-' + Xp
        Yp = Y[j-1] + Yp
    elif i>0 and j==0:
        Xp = X[i-1] + Xp
        Yp = '-' + Yp
        
    return (Xp,Yp)

#test
# Test case 1:
#X2 = "AC"
#Y2 = "TAG"
#M2 = {'-':{'-':0, 'T': -4, 'G':-4, 'A':-4}, 'A':{'-':-2, 'T': 2, 'G':2, 'A':5}, 'C':{'-':-2, 'T': 2, 'G':2, 'A':2}}
#S2 = ComputeGlobalScore(X2, Y2, M2)
#print ComputeAlignment(X2,Y2,M2,S2)
# Expected: ('-AC', 'TAG')


def global_alignment(X, Y, M):
    """
    Returns the global optimal alignment of two sequences X and Y
    """
    #The global optimal matrix
    S = ComputeGlobalScore(X,Y,M)
    
    (Xp,Yp) = ComputeAlignment(X,Y,M,S)
    return (Xp,Yp)



#Helper Function
def ComputeLocalScore(X,Y,M):
    """
    Computes the local optimal score matrix for dynamic programming
    
    Input: Strings X and Y, and scoring matrix M
    Returns: The optimal score matrix S
    """
    S = []
    lenX = len(X)
    lenY = len(Y)
    
    #Initialize the matrix
    for i in range(0,lenX+1):
        S.append([])
        for j in range(0,lenY+1):
            S[i].append(0)
    #print S   
    
    #Starts the DP
    for i in range(0, lenX+1):
        
        for j in range(0, lenY+1):
            
            if i==0 and j>0:
                S[i][j] = max(S[i][j-1] + M['-'][Y[j-1]],0)
            
            elif i>0 and j==0:
                S[i][j] = max(S[i-1][j] + M[X[i-1]]['-'],0)
            
            elif i>0 and j>0:
                S[i][j] = max(S[i-1][j-1] + M[X[i-1]][Y[j-1]],
                             S[i][j-1] + M['-'][Y[j-1]],
                             S[i-1][j] + M[X[i-1]]['-'],0)
            
    return S
#test
#SL1 = ComputeLocalScore(X1,Y1,M1)
#print "SL1",SL1
#print ComputeLocalScore(X3,Y1,M)


def ComputeLocalAlignment(X,Y,M,S):
    """
    Trace back the alignment through the DP's local optimal score matrix
    
    Inputs: Strings X and Y, scoring Matrix M, and the optimal score matrix
    Returns: A tuple of two strings
    """
    #find index of the largest number in S
    i = len(X)
    j = len(Y)
    for m in range(0,len(X)+1):
        for n in range(0,len(Y)+1):
            if S[m][n] > S[i][j]:
                i = m
                j = n
    #print S[i][j]
    Xp = ''
    Yp = ''
    
    #The TraceBack
    while (i>0) and (j>0) and (S[i][j] != 0) :
        if S[i][j] == S[i-1][j-1] + M[X[i-1]][Y[j-1]]:
            Xp = X[i-1] + Xp
            Yp = Y[j-1] + Yp
            i -= 1
            j -= 1
        elif S[i][j] == S[i][j-1] + M['-'][Y[j-1]]:
            Xp = '-' + Xp
            Yp = Y[j-1] + Yp
            j -= 1
        elif S[i][j] == S[i-1][j] + M[X[i-1]]['-']:
            Xp = X[i-1] + Xp
            Yp = '-' + Yp
            i -= 1

        
    return (Xp,Yp)

#print ComputeLocalAlignment(X1,Y1,M1,SL1)
X5 = 'ACC'
Y5 = 'TTTACACGG'
M5 = {'-':{'-':0, 'T': -4, 'G':-4, 'A':-4,'C':-4}, 'A':{'-':-4, 'T': 2, 'G':2,'C':2, 'A':10}, 'C':{'-':-4, 'T': 2, 'G':2, 'A':2,'C':10}}
SL5 = ComputeLocalScore(X5,Y5,M5)
#print ComputeLocalAlignment(X5,Y5,M5,SL5)
#Expected:('AC-C','ACAC')

def local_alignment(X,Y,M):
    """
    Returns the local optimal alignment of two sequences X and Y
    """
    
    #The local optimal matrix
    S = ComputeLocalScore(X,Y,M)
    a = ComputeLocalAlignment(X,Y,M,S)
    return a

#test
X4 = "ACCC"
Y4 = "TAGTA"
M4 = {'-':{'-':0, 'T': -4, 'G':-4, 'A':-4}, 'A':{'-':-2, 'T': 2, 'G':2, 'A':5}, 'C':{'-':-2, 'T': 2, 'G':2, 'A':2}}
#print '',local_alignment(X4,Y4,M4)
#Expected: ('ACCC', 'AGTA')

print global_alignment(X5,Y5,M5)
print local_alignment(X5,Y5,M5)

    