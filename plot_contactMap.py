#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 11 16:25:29 2024

@author: jitong2023
@description: predict RNA structure and parse probability for visualization
"""

import subprocess
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib import cm
import matplotlib.colors as mcolors
import math
import pathlib
from datetime import datetime as dt

from io import BytesIO
import base64
#from plotly.tools import mpl_to_plotly
import matplotlib
matplotlib.use('agg')

RNAFOLD_PATH = pathlib.Path('RNAfold')
BASE_PATH = pathlib.Path(__file__).parent.resolve()
TEMPDATA_PATH = BASE_PATH.joinpath('temp_data').resolve()

def runRNAfold(filePrefix, dataPath, rnafoldPath=str(RNAFOLD_PATH)):
    cmd = [rnafoldPath, '--noLP', '-d2', '-p', '--id-prefix='+filePrefix, '-i', filePrefix+'.fasta', '>', filePrefix+'.out']
    subprocess.run(' '.join(cmd), shell=True, cwd = dataPath)
    return

def parseRNAfold_prob(filePrefix):
    fh = open(filePrefix+'.out', 'r')
    fh.readline()
    seq = fh.readline().strip()
    dot = fh.readline().strip().split(' ')[0]
    fh.close()

    fh = open(filePrefix+'_0001_dp.ps')
    probDict = dict()
    probRecord = []
    for i in range(len(dot)):
        probDict[i] = [0]
    while True:
        newLine = fh.readline()
        if newLine.startswith('%start of base pair probability data'):
            while True:
                newLine = fh.readline()
                if len(newLine) == 0:
                    break
                toks = newLine.strip().split(' ')
                if toks[-1]!='ubox':
                    break
                pos1 = int(toks[0])-1
                pos2 = int(toks[1])-1
                prob = float(toks[2])
                probDict[pos1].append(prob)
                probDict[pos2].append(prob)
                probRecord.append((pos1, pos2, prob))
        if len(newLine)==0:
            break
    fh.close()

    probList = []
    for i in range(len(dot)):
        if dot[i]=='.':
            probList.append(1-max(probDict[i]))
        else:
            probList.append(max(probDict[i]))
            
    return([seq, dot, probList, probRecord])

def rotate_point(x, y, angle_degrees):
    # Convert angle from degrees to radians
    angle_radians = math.radians(angle_degrees)
    
    # Calculate the new coordinates
    x_prime = x * math.cos(angle_radians) + y * math.sin(angle_radians)
    y_prime = -x * math.sin(angle_radians) + y * math.cos(angle_radians)
    
    return x_prime, y_prime

def predict_structure(seq, seqname, dataPath=str(TEMPDATA_PATH)):
    # generate fastafile
    fastaFilePrefix = dt.now().strftime('%d%H%M%S%f')
    with open(dataPath + '/' + fastaFilePrefix + '.fasta', 'w') as fh:
        fh.write('>' + seqname + '\n')
        fh.write(seq + '\n')
    
    # structure prediction for reference sequence fastaFilePrefix.fasta
    runRNAfold(fastaFilePrefix, dataPath)
    [seq, dot, scores, probRecord] = parseRNAfold_prob(dataPath + '/' + fastaFilePrefix)
    subprocess.run(' '.join(['rm', fastaFilePrefix+'*']), shell = True, cwd=dataPath)
    return [seq, dot, scores, probRecord]

def fig_to_uri(in_fig, close_all=True, **save_args):
    # type: (plt.Figure) -> str
    """
    Save a figure as a URI
    :param in_fig:
    :return:
    """
    out_img = BytesIO()
    in_fig.savefig(out_img, format='png', **save_args)
    if close_all:
        in_fig.clf()
        plt.close('all')
    out_img.seek(0)  # rewind file
    encoded = base64.b64encode(out_img.read()).decode("ascii").replace("\n", "")
    return "data:image/png;base64,{}".format(encoded)

def plot_contactMap(seqInfo):
    # parse structure predictions
    #[seq, dot, scores, probRecord] = parseRNAfold_prob(fastaFilePrefix)
    [seq, dot, scores, probRecord] = seqInfo
    
    # set base pairing distance limits
    probRecord.sort(key=lambda x: x[-1])
    
    # setup the normalization and the colormap
    probs = np.array([x[-1] for x in probRecord])
    normalize = mcolors.Normalize(vmin=probs.min(), vmax=probs.max())
    colormap = cm.Blues
    
    # plot
    #fig, ax = plt.subplots(figsize=(16, 8))
    fig, ax = plt.subplots(figsize=(5, 3.8), dpi=150)
    for (pos1, pos2, prob) in probRecord:
        # original dotplot
        #ax.plot([pos1, pos1], [pos1, pos2], c=colormap(prob), alpha=0.5)
        #ax.plot([pos1, pos2], [pos2, pos2], c=colormap(prob), alpha=0.5)
        
        # rotate 45 degrees
        start_point = [pos1, 0.0] # [pos1, pos1]
        end_point = [pos2, 0.0] # [pos2, pos2]
        mid_point = rotate_point(pos1, pos2, 45) # [pos1, pos2]
        mid_point = [mid_point[0]/math.sqrt(2), mid_point[1]/math.sqrt(2)]
        ax.plot([start_point[0]+1, mid_point[0]+1], [start_point[1], mid_point[1]], c=colormap(prob), alpha=0.5)
        ax.plot([mid_point[0]+1, end_point[0]+1], [mid_point[1], end_point[1]], c=colormap(prob), alpha=0.5)
        
    # setup the colorbar
    scalarmappaple = cm.ScalarMappable(norm=normalize, cmap=colormap)
    scalarmappaple.set_array(probs)
    fig.colorbar(scalarmappaple, orientation="vertical",aspect=10, shrink=0.5, ax=plt.gca())
    
    # adjust axes
    tick_space = 10 ** math.floor(math.log10(len(seq)))
    ax.set_aspect('equal', adjustable='box')
    ax.set_frame_on(False)
    ax.get_yaxis().set_visible(False)
    ax.set_xticks(np.arange(0, len(seq), tick_space)[1:])
    xmin, xmax = ax.get_xaxis().get_view_interval()
    ymin, ymax = ax.get_yaxis().get_view_interval()
    ax.add_artist(Line2D((1, len(seq)), (ymin, ymin), color='black', linewidth=2))
    
    # show the figure
    #plt.savefig(fastaFilePrefix+'.png', dpi=300)
    
    # output to dash app
    #plotly_fig = mpl_to_plotly(fig)
    plotly_fig = fig_to_uri(fig)
    return plotly_fig

def plot_contactMapAlign(seqInfo, queryseqInfo, alignPos):
    
    ### plot the template structure ###
    # parse structure predictions
    #[seq, dot, _, probRecord] = parseRNAfold_prob(fastaFilePrefix)
    [seq, dot, _, probRecord] = seqInfo

    # set base pairing distance limits
    probRecord.sort(key=lambda x: x[-1])
    
    # setup the normalization and the colormap
    probs = np.array([x[-1] for x in probRecord])
    normalize = mcolors.Normalize(vmin=probs.min(), vmax=probs.max())
    colormap = cm.Blues
    
    # setup y axis adjustment
    yadj = len(seq) * 0.02
    ytick_length = len(seq) * 0.01
    
    # plot
    #fig, ax = plt.subplots(figsize=(16, 8))
    fig, ax = plt.subplots(figsize=(5, 3.8), dpi=150)
    for (pos1, pos2, prob) in probRecord:
        # original dotplot
        #ax.plot([pos1, pos1], [pos1, pos2], c=colormap(prob), alpha=0.5)
        #ax.plot([pos1, pos2], [pos2, pos2], c=colormap(prob), alpha=0.5)
        
        # rotate 45 degrees
        start_point = [pos1, 0.0] # [pos1, pos1]
        end_point = [pos2, 0.0] # [pos2, pos2]
        mid_point = rotate_point(pos1, pos2, 45) # [pos1, pos2]
        mid_point = [mid_point[0]/math.sqrt(2), mid_point[1]/math.sqrt(2)]
        ax.plot([start_point[0]+1, mid_point[0]+1], [start_point[1]+yadj, mid_point[1]+yadj], c=colormap(prob), alpha=0.5)
        ax.plot([mid_point[0]+1, end_point[0]+1], [mid_point[1]+yadj, end_point[1]+yadj], c=colormap(prob), alpha=0.5)
        
    ### plot the query structure ###
    # parse structure predictions
    #[seq_query, dot_query, _, probRecord] = parseRNAfold_prob(alignFilePrefix)
    [seq_query, dot_query, _, probRecord] = queryseqInfo
    
    # set base pairing distance limits
    probRecord.sort(key=lambda x: x[-1])
    
    # setup the normalization and the colormap
    probs = np.array([x[-1] for x in probRecord])
    
    # plot
    for (pos1, pos2, prob) in probRecord:
        # original dotplot
        #ax.plot([pos1, pos2], [pos1, pos1], c=colormap(prob), alpha=0.5)
        #ax.plot([pos2, pos2], [pos1, pos2], c=colormap(prob), alpha=0.5)
        
        # rotate 45 degrees
        start_point = [pos1, 0.0] # [pos1, pos1]
        end_point = [pos2, 0.0] # [pos2, pos2]
        mid_point = rotate_point(pos1, pos2, 45) # [pos1, pos2]
        mid_point = [mid_point[0]/math.sqrt(2), -mid_point[1]/math.sqrt(2)]
        ax.plot([start_point[0]+alignPos, mid_point[0]+alignPos], [start_point[1]-yadj, mid_point[1]-yadj], c=colormap(prob), alpha=0.5)
        ax.plot([mid_point[0]+alignPos, end_point[0]+alignPos], [mid_point[1]-yadj, end_point[1]-yadj], c=colormap(prob), alpha=0.5)
        
    # setup the colorbar
    scalarmappaple = cm.ScalarMappable(norm=normalize, cmap=colormap)
    scalarmappaple.set_array(probs)
    fig.colorbar(scalarmappaple, orientation="vertical",aspect=10, shrink=0.3, ax=plt.gca())
    
    # adjust axes
    ax.set_aspect('equal', adjustable='box')
    ax.axis('off')
    ax.hlines(y=0, xmin=1, xmax=len(seq), linewidth=2, color='black')
    for xtick in np.arange(0, len(seq), 10)[1:]:
        ax.vlines(xtick, -ytick_length, 0, linewidth=1)
    
    # show the figure
    #plt.savefig(fastaFilePrefix+'_'+alignFilePrefix+'.png', dpi=300)
    
    # output to dash app
    #plotly_fig = mpl_to_plotly(fig)
    plotly_fig = fig_to_uri(fig)
    return alignPos, plotly_fig

