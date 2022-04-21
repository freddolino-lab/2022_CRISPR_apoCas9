#!/usr/bin/env python2

"""
This script calculate the Position Weight Matrix of a fasta file consisting sequences of the same length, *weighted by a score supplied in the Fasta entry ids*. Base vocabulary: ACGT. 
Modified from /home/diaorch/data/cz/cz-7-pwm.py, updated to output both line plot and scatter plot; added line at x = 25 and y = 0.25; added multi-line of plot title to include sample name
USAGE:
python2 cz-05-pwm.py 04-cov/czFeb2018_131_S1.padded.fasta 05-pwm
"""

from __future__ import print_function
import sys
import numpy as np
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import matplotlib.cm as cm
from os.path import basename
import re

def convertBaseToCode(base):
    """
    Convert character of a base in to interger code
    A -> 1, C -> 2, G -> 3, T -> 4
    """
    if (base == 'A'):
        return(1)
    if (base == 'C'):
        return(2)
    if (base == 'G'):
        return(3)
    if (base == 'T'):
        return(4)
    

def calcPFM(inputIterator, seqLen):
    """
    Calculate the Position Frequency Matrix for a fix length of sequences in an iterator of Bio.SeqRecord from Fasta by previous functions, counting the appearance of each base in each position
    """
    # seqCount = 0 #length of `iterator` is unknown until iterated through
    seqPFM = np.zeros((4, seqLen))
    for record in inputIterator:
        # seqCount += 1
        # print(record.id)
        scoreExtract = re.search('\[(.*)\]', record.id).group(1)
        # print(scoreExtract)
        if (scoreExtract is None):
            sys.exit("ERROR: at least one FASTA sequence entry do not contain a score")
        seqBaseSet = list(record.upper().seq)
        try:
            score = int(scoreExtract)
        except ValueError:
            expSplit = scoreExtract.split('e')
            score = int(float(expSplit[0]) * 10 ** int(expSplit[1]))
        for pos in range(0, seqLen):
            baseCode = convertBaseToCode(seqBaseSet[pos])
            seqPFM[baseCode - 1, pos] += score
    return(seqPFM)

def calcPWM(pfm):
    """
    Calculate the Position Weight Matrix for a np.array Position Frequency Matrix
    """
    sumSet = set(pfm.sum(axis = 0))
    if (len(sumSet) == 1):
        # print(sumSet)
        seqCountSum = list(sumSet)[0]
        pwm = pfm / seqCountSum
    else:
        sys.exit("ERROR: input sequences of different length in fasta file")
    return(pwm) 

def plotPWM(pwm, saveNameBase, upDownBreakPoint = 50):
    """
    Plot the PWM in line plot with each base in different color
    """
    plt.switch_backend('agg')
    plt.ioff()
    posCount = pwm.shape[1]
    x = np.linspace(1, posCount, posCount)
    y = pwm
    fig = plt.figure()
    plt.plot(x, y.T)
    plt.ylim(ymin = 0.0, ymax = 1.0)
    plt.axvline(x = upDownBreakPoint, color = 'black', linestyle='dashed', zorder = 0)
    plt.axhline(y = 0.25, color = 'black', linestyle='dashed', zorder = 0)
    plt.axhline(y = 0.2560905358828566, color = 'black', zorder = 0)
    plt.axhline(y = 0.2581646945182323, color = 'black', zorder = 0)
    plt.legend(labels = ['A', 'C', 'G', 'T'])
    figNameBase = basename(saveNameBase).split('.')[0]
    plt.title(figNameBase + '\n Position Weight Matrix visualization, weight at each position')
    plt.xlabel('Position in sequence')
    plt.ylabel('Weight (aka. frequency) of each base')
    figName = saveNameBase + '.line.png'
    plt.savefig(figName)
    plt.close(fig)

def plotPWMscatter(pwm, saveNameBase, xSpacing = 1, upDownBreakPoint = 50):
    """
    Plot the PWM in scatter plot with each base in different color and minro x-axis grid
    """
    plt.switch_backend('agg')
    posCount = pwm.shape[1]
    X = np.linspace(1, posCount, posCount)
    # x = np.tile(X, (4, 1))
    x = X
    Y = pwm
    fig = plt.figure()
    ax = plt.subplot(111)
    minorLocator = MultipleLocator(xSpacing)
    ax.xaxis.set_minor_locator(minorLocator)
    ax.grid(which = 'minor')
    ax.set_axisbelow(True)
    ax.set_ylim([0.0, 1.0])
    # ax.scatter(x, y.T, label = ['A', 'C', 'G', 'T'])
    colors = cm.Set1(np.linspace(0, 0.5, pwm.shape[0]))
    for y, c in zip(Y, colors):
        ax.scatter(x, y, color = c, zorder = 10)
    ax.legend(labels = ['A', 'C', 'G', 'T'])
    figNameBase = basename(saveNameBase).split('.')[0]
    plt.title(figNameBase + '\nPosition Weight Matrix visualization, weight at each position')
    plt.xlabel('Position in sequence')
    plt.ylabel('Weight (aka. frequency) of each base')
    plt.axvline(x = upDownBreakPoint, color = 'black', linestyle='dashed', zorder = 5)
    plt.axhline(y = 0.25, color = 'black', linestyle='dashed', zorder = 5)
    plt.axhline(y = 0.2560905358828566, color = 'black', zorder = 5)
    plt.axhline(y = 0.2581646945182323, color = 'black', zorder = 5)
    figName = saveNameBase + '.dot.png'
    plt.savefig(figName)
    plt.close(fig)
    

def readFastaToIterator(iFilename):
    """
    read fasta file named as `iFilename` using Biopython SeqIO parse, return an iterator of entry ids and sequences 
    PURPOSE: to make the input fasta file more generalizable, e.g. different line break conventions in sequnces
    """
    inputSeqIterator = SeqIO.parse(iFilename, 'fasta')
    # print(list(inputSeqIterator))
    # check sequences are of the same length
    seqLen = -1
    for record in inputSeqIterator:
        if (seqLen == -1):
            seqLen = len(record.seq)
        elif (len(record.seq) != seqLen):
            sys.exit("ERROR: input sequences of different length in fasta file")
    # print(list(inputSeqIterator))
    return (SeqIO.parse(iFilename, 'fasta'), seqLen)

def main():
    inputFilename = sys.argv[1]
    outputDir = sys.argv[2]
    inputFilenameNoPath = os.path.basename(inputFilename)
    inputFilenameBase = ('.').join(inputFilenameNoPath.split('.')[:-1])
    if(not outputDir.endswith('/')):
        outputDir = outputDir + '/'
    
    fastaIter, seqLen = readFastaToIterator(iFilename = inputFilename)
    # print(seqLen)
    # print(list(fastaIter))
    
    pfm = calcPFM(fastaIter, seqLen)
    # print(pfm)
    pwm = calcPWM(pfm)
    # print(pwm)
    saveNameBase = outputDir + inputFilenameBase
    # save PWM
    np.savetxt(saveNameBase + '.csv', pwm, delimiter = ',')
    np.save(saveNameBase + '.npy', pwm)
    # plot PWM
    plotPWMscatter(pwm, saveNameBase)
    plotPWM(pwm, saveNameBase)
    

if __name__ == '__main__':
    main()
