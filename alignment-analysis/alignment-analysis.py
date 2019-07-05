#Analysis
#Henry Zhang
#jz77

import random
import matplotlib.pyplot as plt
import pylab
import types
import time
import math
import copy
import numpy as np
from collections import *

def permute_string(s):
    """
    Return a new string with the characters in s randomly permuted.

    Arguments:
    s -- string

    Returns:
    Random permutation of s
    """
    charlist = list(s)
    random.shuffle(charlist)
    newstr = "".join(charlist)
    return newstr

def read_protein(filename):
    """
    Read a protein sequence from the file named filename.

    Arguments:
    filename -- name of file containing a protein sequence

    Returns:
    A string representing the protein
    """
    with open(filename) as f:
        p = f.read()
    p = p.rstrip()
    return p

def read_scoring_matrix(filename):
    """
    Read a scoring matrix from the file named filename.  

    Argument:
    filename -- name of file containing a scoring matrix

    Returns:
    A dictionary of dictionaries mapping X and Y characters to scores
    """
    M = {}
    with open(filename) as f:
        ykeys = f.readline()
        ykeychars = ykeys.split()
        for line in f.readlines():
            vals = line.split()
            xkey = vals.pop(0)
            M[xkey] = {}
            for ykey, val in zip(ykeychars, vals):
                M[xkey][ykey] = int(val)
    return M

def read_graph(filename):
    """
    Read a graph from a file.  The file is assumed to hold a graph
    that was written via the write_graph function.

    Arguments:
    filename -- name of file that contains the graph

    Returns:
    The graph that was stored in the input file.
    """
    with open(filename) as f:
        g = eval(f.read())
    return g

def write_graph(g, filename):
    """
    Write a graph to a file.  The file will be in a format that can be
    read by the read_graph function.

    Arguments:
    g        -- a graph
    filename -- name of the file to store the graph

    Returns:
    None
    """
    with open(filename, 'w') as f:
        f.write(repr(g))

def copy_graph(g):
    """
    Return a copy of the input graph, g

    Arguments:
    g -- a graph

    Returns:
    A copy of the input graph that does not share any objects.
    """
    return copy.deepcopy(g)

## Timing functions

def time_func(f, args=[], kw_args={}):
    """
    Times one call to f with args, kw_args.

    Arguments:
    f       -- the function to be timed
    args    -- list of arguments to pass to f
    kw_args -- dictionary of keyword arguments to pass to f.

    Returns: 
    a tuple containing the result of the call and the time it
    took (in seconds).

    Example:

    >>> def sumrange(low, high):
            sum = 0
            for i in range(low, high):
                sum += i
            return sum
    >>> time_func(sumrange, [82, 35993])
    (647726707, 0.01079106330871582)
    >>> 
    """
    start_time = time.time()
    result = f(*args, **kw_args)
    end_time = time.time()

    return (result, end_time - start_time)

## Plotting functions

def show():
    """
    Do not use this function unless you have trouble with figures.

    It may be necessary to call this function after drawing/plotting
    all figures.  If so, it should only be called once at the end.

    Arguments:
    None

    Returns:
    None
    """
    plt.show()

def plot_dist_linear(data, title, xlabel, ylabel, filename=None):
    """
    Plot the distribution provided in data as a bar plot on a linear
    scale.

    Arguments: 
    data     -- dictionary which will be plotted with the keys
                on the x axis and the values on the y axis
    title    -- title label for the plot
    xlabel   -- x axis label for the plot
    ylabel   -- y axis label for the plot
    filename -- optional name of file to which plot will be
                saved (in png format)

    Returns:
    None
    """
    _plot_dist(data, title, xlabel, ylabel, False, filename)

def plot_dist_loglog(data, title, xlabel, ylabel, filename=None):
    """
    Plot the distribution provided in data as a scatter plot on a
    loglog scale.

    Arguments: 
    data     -- dictionary which will be plotted with the keys
                on the x axis and the values on the y axis
    title    -- title label for the plot
    xlabel   -- x axis label for the plot
    ylabel   -- y axis label for the plot
    filename -- optional name of file to which plot will be
                saved (in png format)

    Returns:
    None
    """
    _plot_dist(data, title, xlabel, ylabel, True, filename)


def _pow_10_round(n, up=True):
    """
    Round n to the nearest power of 10.

    Arguments:
    n  -- number to round
    up -- round up if True, down if False

    Returns:
    rounded number
    """
    if up:
        return 10 ** math.ceil(math.log(n, 10))
    else:
        return 10 ** math.floor(math.log(n, 10))
        

def _plot_dist(data, title, xlabel, ylabel, scatter, filename=None):
    """
    Plot the distribution provided in data.

    Arguments: 
    data     -- dictionary which will be plotted with the keys
                on the x axis and the values on the y axis
    title    -- title label for the plot
    xlabel   -- x axis label for the plot
    ylabel   -- y axis label for the plot
    scatter  -- True for loglog scatter plot, False for linear bar plot
    filename -- optional name of file to which plot will be
                saved (in png format)

    Returns:
    None
    """
    ### Check that the data is a dictionary
    if not isinstance(data, types.DictType):
        msg = "data must be a dictionary, not {0}".format(type(data).__name__)
        raise TypeError(msg)

    ### Create a new figure
    fig = pylab.figure()

    ### Plot the data
    if scatter:
        _plot_dict_scatter(data)
    else:
        _plot_dict_bar(data, 0)
    
    ### Label the plot
    pylab.title(title)
    pylab.xlabel(xlabel)
    pylab.ylabel(ylabel)

    ### Draw grid
    gca = pylab.gca()
    gca.yaxis.grid(True)
    gca.xaxis.grid(False)

    if scatter:
        ### Use loglog scale
        gca.set_xscale('log')
        gca.set_yscale('log')
        gca.set_xlim([_pow_10_round(min([x for x in data.keys() if x > 0]), False), 
                      _pow_10_round(max(data.keys()))])
        gca.set_ylim([_pow_10_round(min([x for x in data.values() if x > 0]), False), 
                      _pow_10_round(max(data.values()))])

    ### Show the plot
    fig.show()

    ### Save to file
    if filename:
        pylab.savefig(filename)

def plot_lines(data, title, xlabel, ylabel, labels=None, filename=None):
    """
    Plot a line graph with the provided data.

    Arguments: 
    data     -- a list of dictionaries, each of which will be plotted 
                as a line with the keys on the x axis and the values on
                the y axis.
    title    -- title label for the plot
    xlabel   -- x axis label for the plot
    ylabel   -- y axis label for the plot
    labels   -- optional list of strings that will be used for a legend
                this list must correspond to the data list
    filename -- optional name of file to which plot will be
                saved (in png format)

    Returns:
    None
    """
    ### Check that the data is a list
    if not isinstance(data, types.ListType):
        msg = "data must be a list, not {0}".format(type(data).__name__)
        raise TypeError(msg)

    ### Create a new figure
    fig = pylab.figure()

    ### Plot the data
    if labels:
        mylabels = labels[:]
        for i in range(len(data)-len(labels)):
            mylabels.append("")
        for d, l in zip(data, mylabels):
            _plot_dict_line(d, l)
        # Add legend
        pylab.legend(loc='best')
        gca = pylab.gca()
        legend = gca.get_legend()
        pylab.setp(legend.get_texts(), fontsize='medium')
    else:
        for d in data:
            _plot_dict_line(d)

    ### Set the lower y limit to 0 or the lowest number in the values
    mins = [min(l.values()) for l in data]
    ymin = min(0, min(mins))
    pylab.ylim(ymin=ymin)

    ### Label the plot
    pylab.title(title)
    pylab.xlabel(xlabel)
    pylab.ylabel(ylabel)

    ### Draw grid lines
    pylab.grid(True)

    ### Show the plot
    fig.show()

    ### Save to file
    if filename:
        pylab.savefig(filename)

def _dict2lists(data):
    """
    Convert a dictionary into a list of keys and values, sorted by
    key.  

    Arguments:
    data -- dictionary

    Returns:
    A tuple of two lists: the first is the keys, the second is the values
    """
    xvals = data.keys()
    xvals.sort()
    yvals = []
    for x in xvals:
        yvals.append(data[x])
    return xvals, yvals

def _plot_dict_line(d, label=None):
    """
    Plot data in the dictionary d on the current plot as a line.

    Arguments:
    d     -- dictionary
    label -- optional legend label

    Returns:
    None
    """
    xvals, yvals = _dict2lists(d)
    if label:
        pylab.plot(xvals, yvals, label=label)
    else:
        pylab.plot(xvals, yvals)

def _plot_dict_bar(d, xmin=None, label=None):
    """
    Plot data in the dictionary d on the current plot as bars. 

    Arguments:
    d     -- dictionary
    xmin  -- optional minimum value for x axis
    label -- optional legend label

    Returns:
    None
    """
    xvals, yvals = _dict2lists(d)
    if xmin == None:
        xmin = min(xvals) - 1
    else:
        xmin = min(xmin, min(xvals) - 1)
    if label:
        pylab.bar(xvals, yvals, align='center', label=label)
        pylab.xlim([xmin, max(xvals)+1])
    else:
        pylab.bar(xvals, yvals, align='center')
        pylab.xlim([xmin, max(xvals)+1])

def _plot_dict_scatter(d):
    """
    Plot data in the dictionary d on the current plot as points. 

    Arguments:
    d     -- dictionary

    Returns:
    None
    """
    xvals, yvals = _dict2lists(d)
    pylab.scatter(xvals, yvals)

########################### Autograder Functions##################
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


def local_alignment(X,Y,M):
    """
    Returns the local optimal alignment of two sequences X and Y
    """
    
    #The local optimal matrix
    S = ComputeLocalScore(X,Y,M)
    a = ComputeLocalAlignment(X,Y,M,S)
    return a
    
########################### My Functions ############################

#helper function
def find_max(M):
    """
        
    """
    lst = []
    for row in M:
        lst.append(max(row))
    maxscore = max(lst)
    return maxscore
#test
#print find_max([[0,1,2],[2,3,4]])
  
     
def generate_null_distribution(X,Y,M,t,numiter):
    """
    
    """
    scoredata =[]
    frequency = defaultdict(int)
    #storing the scores with number of operation
    for i in range (0,numiter):
        Yp = permute_string(Y)
        #cases for t
        if t==0:
            S1 = ComputeLocalScore(X,Yp,M)
            #Find maximum
            maxscore = find_max(S1)
            scoredata.append(maxscore)
            frequency[maxscore] += 1
        else:
            S2 = ComputeGlobalScore(X,Yp,M)
            #maximum score is the bottom right of S2
            maxscore = S2[len(X)][len(Yp)]
            scoredata.append(maxscore)
            frequency[maxscore] += 1
    #plot    
    _plot_dict_bar(frequency)      
    #Compute the score for X,Y
    if t==0:
        S3 = ComputeLocalScore(X,Y,M)
        s = find_max(S3)
    else:
        S4 = ComputeGlobalScore(X,Y,M)
        s = S4[len(X)][len(Y)]
    
    mean = np.mean(scoredata)
    std = np.std(scoredata)
    z = (s-float(mean))/std
    
    return "mean",mean, "std", std, "z",z

# Test case 1:
#M1 = {'-':{'-':0, 'T': -4, 'G':-4, 'A':-4,'C':-4}, 'A':{'-':-4, 'T': 2, 'G':2, 'A':10,'C':2}, 'C':{'-':-4, 'T': 2, 'G':2, 'A':2,'C':10}}
#X1 ='ACCAACAA'
#Y1= 'TTACACGGTTGACAGA'
#t1 = 0
#numiter1 = 10
#print generate_null_distribution(X1,Y1,M1,t1,numiter1)
    
##############################Eyeless######################
X = read_protein('HumanEyelessProtein.py')
Y = read_protein('FruitflyEyelessProtein.py')
M = read_scoring_matrix('PAM50.py')
#global score
Sg = ComputeGlobalScore(X,Y,M)
#print "max score for global,",Sg[len(X)][len(Y)]

#local score
Sl = ComputeLocalScore(X,Y,M)
#print "max score for local,", find_max(Sl)

#local alignment
#print "local alignment,", local_alignment(X,Y,M)

#Analyze the null_generate
print "local,", generate_null_distribution(X,Y,M,0,100)
print "-------------------"
print generate_null_distribution(X,Y,M,1,100)


    
    
