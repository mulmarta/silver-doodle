#!/bin/python3

import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pdb
import math


HASHSIZE = 512
G1 = 512
SIGSIZE = HASHSIZE + G1
TAGSIZE = 512
IDs = 128
AUTHDATA = 512
CREDENTIALS = 1500 * 8
EXTENSIONS = 130 * 8
KEYPACKAGE = SIGSIZE + CREDENTIALS + EXTENSIONS
HEADER = 2243 * 8 #SIGSIZE + TAGSIZE * 2 + IDs + AUTHDATA + KEYPACKAGE


DHCT0 = 0
DHCTI = 0
DHCONST = G1
DHKEY = G1

MLSCT0 = 0
MLSCTI = 640 # (512 bit Sym + 128bit AD?)
MLSPK = 512
MLSCTX = 512

CONVERSION_FACTOR = 8 * 1024
CF_MB = 8 * 1024 * 1024

def log(x):
    return int(math.log(x, 2))

def calcMLSSenderSizes(rng, case):
    if case == 0: #Best
        data = [HEADER + log(x)* (MLSCTI + MLSPK + MLSCTX) for x in rng]
    elif case == 1: #Worst
        data = [HEADER +  x * (MLSCTI + MLSCTX) + log(x) * MLSPK for x in rng]
    else:
        raise Error
    return applyConversion(data)

def calcBGMSenderSizes(rng, case):
    if case == 0: #Best
        data = [HEADER + DHCONST + log(x)* (MLSCTI + MLSPK) for x in rng]
    elif case == 1: #Worst
        data = [HEADER +  DHCONST + x * (MLSCTI) + log(x) * MLSPK for x in rng]
    else:
        raise Error
    return applyConversion(data)

def calcMLSIndivSenderSizes(rng, case):
    if case == 0: #Best
        data = [HEADER + log(x)* (MLSCTI + MLSPK + MLSCTX) for x in rng]
    elif case == 1: #Worst
        data = [HEADER + x * (MLSCTI + MLSCTX) + log(x) * MLSPK for x in rng]
    else:
        raise Error
    return applyConversion(data)

def calcCMPKESenderSizes(rng, case):
    data = [HEADER + DHCONST + MLSPK + x * (MLSCTI) for x in rng]
    return applyConversion(data)


def calcMLSRecSizes(rng, case):
    return calcMLSSenderSizes(rng, case)

def calcBGMRecSizes(rng, case):
    data = [(HEADER + DHCONST + MLSCTI + log(x) * MLSPK) for x in rng]
    return applyConversion(data)

def calcMLSIndivRecSizes(rng, case):
    data = [HEADER + DHCONST + MLSCTI + log(x) * MLSPK for x in rng]
    return applyConversion(data)

def calcCMPKERecSizes(rng, case):
    data = [HEADER + DHCONST + MLSCTI + MLSPK for x in rng]
    return applyConversion(data)

def calcMLSRecSum(rng, case):
    tmp = calcMLSRecSizes(rng, case)
    data = applyConversion([x * y for (x,y) in zip(tmp, rng)], 1024)
    return data


def BGMRecSizeSum(n):
    sum = 0
    for depth in range(1, log(n)-1):
        sum += int(n / math.pow(2, depth)) * (HEADER + DHCONST + MLSCTI + (log(n) - depth) * MLSPK)
    return sum

def calcBGMRecSum(rng, case, accType = "Hash"):
    data = applyConversion([BGMRecSizeSum(x) for x in rng], CF_MB)
    data2 = applyConversion(calcBGMSenderSizes(rng, case), 1024)
    erg = [x + y for (x,y) in zip(data, data2)]
    return erg

def calcCMPKERecSizesSum(rng, case):
    data = applyConversion([x * (HEADER + DHCONST + MLSCTI + MLSPK) for x in rng], CF_MB)
    data2 = applyConversion(calcCMPKESenderSizes(rng, case), 1024)
    erg = [x + y for (x,y) in zip(data, data2)]
    return erg

def calcMLSRecAvg(rng, case):
    return calcMLSSenderSizes(rng, case)

def calcBGMRecAvg(rng, case):
    data = applyConversion([int(BGMRecSizeSum(x)/x) for x in rng])
    return data

def calcCMPKERecAvg(rng, case):
    return calcCMPKERecSizes(rng, case)

def applyConversion(data, factor = CONVERSION_FACTOR):
    return [x / factor for x in data]

if __name__ == "__main__":
    fig, ((ax1),(ax2)) = plt.subplots(2, 1, sharex=True, figsize=(7.5, 7.5))
    fig.tight_layout(h_pad=6, pad=6.5)
    
    rng = [int(math.pow(2, x)) for x in range(1,17)]
    # rng.append(10000)
    mlsSenderDataBest = calcMLSSenderSizes(rng, 0)
    mlsSenderDataWorst = calcMLSSenderSizes(rng, 1)

    bgmSenderBest = calcBGMSenderSizes(rng, 0)
    bgmSenderWorst = calcBGMSenderSizes(rng, 1)

    cmpkeSender = calcCMPKESenderSizes(rng, 0)
    
    ax1.plot(rng, mlsSenderDataBest, "--", label = "ITK", color="red")    
    ax1.plot(rng, bgmSenderBest, "--", label = "SAIK Huffman", color = "green")
    ax1.plot(rng, cmpkeSender, "-|", label = "CmCGKA", color = "blue")
    
    ax1.fill_between(rng, mlsSenderDataBest, mlsSenderDataWorst, color="red", hatch = "\\", fc = "#ffffff00", ec = "#ff00007f")
    ax1.fill_between(rng, bgmSenderBest, bgmSenderWorst, color="green", hatch = "//", fc = "#ffffff00", ec = "#00ff007f")

    ax1.plot(rng, mlsSenderDataWorst, "--", color="red")
    ax1.plot(rng, bgmSenderWorst, "--", color = "green")
    
    
    ax1.tick_params("both", reset = True)
    ax1.set_xlabel("#Users")
    ax1.set_ylabel("size in KB")
    ax1.set_yscale("log")
    ax1.title.set_text("a) Sender Bandwith")
    
    mlsRecBest = calcMLSRecAvg(rng, 0)
    mlsRecWorst = calcMLSRecAvg(rng, 1)
    
    bgmRecBest = calcBGMRecAvg(rng, 0)
    bgmRecWorst = calcBGMRecAvg(rng, 1)

    cmpkeRecBest = calcCMPKERecAvg(rng, 0)
    cmpkeRecWorst = calcCMPKERecAvg(rng, 1)

    
    ax2.plot(rng, mlsRecBest, "--", label = "ITK", color="red")
    ax2.plot(rng, mlsRecWorst, "--", color="red")

    ax2.plot(rng, bgmRecBest, "--", label = "SAIK", color = "green")
    ax2.plot(rng, bgmRecWorst, "--", color = "green")

    ax2.plot(rng, cmpkeRecBest, "--", label = "CmCGKA", color = "blue")

    ax2.fill_between(rng, mlsRecBest, mlsRecWorst, fc = "#ffffff00", ec="#ff00007f", hatch = "\\\\")
    ax2.fill_between(rng, bgmRecBest, bgmRecWorst, fc = "#ffffff00", ec="#00ff007f", hatch = "//")
    
    ax2.set_ylim([2.25, mlsRecBest[-1] * 1.01])
    ax2.tick_params("both", reset = True)
    ax2.set_xlabel("#Users")
    ax2.set_ylabel("size in KB")
    ax2.title.set_text("b) Receiver Bandwidth")

    mlsRecSumBest = calcMLSRecSum(rng, 0)
    mlsRecSumWorst = calcMLSRecSum(rng, 1)
    
    bgmRecSumBest = calcBGMRecSum(rng, 0)
    bgmRecSumWorst = calcBGMRecSum(rng, 1)
    
    cmpkeRecSum = calcCMPKERecSizesSum(rng, 0)


    mlsRecAvgBest = calcMLSRecAvg(rng, 0)
    mlsRecAvgWorst = calcMLSRecAvg(rng, 1)

    bgmRecAvgBest = calcBGMRecAvg(rng, 0)
    bgmRecAvgWorst = calcBGMRecAvg(rng, 1)

    cmpkeRecAvg = calcCMPKERecAvg(rng, 0)

    green_patch = mpatches.Patch(label='SAIK',fc = "#ffffff00", ec="#00ff007f", hatch = "//")
    red_patch = mpatches.Patch(label='ITK [5]',fc = "#ffffff00", ec="#ff00007f", hatch = "\\\\")
    grey_patch = mpatches.Patch(label='CmCGKA [31]', fc = "#ffffff00", ec = "#0000ff7f", hatch = "||")
    
    
    handles, labels = ax1.get_legend_handles_labels()
    fig.legend(handles = [green_patch, red_patch, grey_patch],loc='upper center',
               ncol=4, bbox_to_anchor=(.5, .975), fontsize=12)    
    
    plt.savefig("Final_Figures_Avg")


    senderPercBest = [x / y * 100 for (x,y) in zip(bgmSenderBest, mlsSenderDataBest)]
    print("Sender Best Percentage")
    print(senderPercBest)
    
    senderPercWorst = [x / y * 100 for (x,y) in zip(bgmSenderWorst, mlsSenderDataWorst)]
    print("Sender Worst Percentage")
    print(senderPercWorst)

    senderPercSaikCmPKE = [x / y * 100 for (x,y) in zip(bgmSenderBest, cmpkeSender)]
    print("Sender cmpke Percentage Best")
    print(senderPercSaikCmPKE)
    print()

    print(bgmSenderBest)
    print(cmpkeSender)
    
    print()
    
    percRecBest = [(x / y) * 100 for (x,y) in zip(bgmRecBest[2:], mlsRecBest[2:])]
    print("Receiver Best Percentage")
    print(percRecBest)

    percRecWorst = [(x / y) * 100 for (x,y) in zip(bgmRecWorst[2:], mlsRecWorst[2:])]
    print("Receiver Worst Percentage")
    print(percRecWorst)

    print()
    print(bgmRecBest)
    print(cmpkeRecBest)
    print()

    
    recPercSaikCmPKE = [x / y * 100 for (x,y) in zip(bgmRecBest[2:], cmpkeRecBest[2:])]
    print("Receiver cmpke Percentage Best")
    print(recPercSaikCmPKE)    
