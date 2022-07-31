#!/bin/python3

import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pdb
import math


HASHSIZE = 512
G1 = 512
G2 = 1024
GT = 15360
AEADSIZE = 512 + 256
SIGSIZE = HASHSIZE + G1
TAGSIZE = 64 * 8

RSAMOD = 15360


DHCT0 = 0
DHCTI = 0
DHCONST = G1
DHKEY = G1

MLSCT0 = 0
MLSCTI = 640
MLSPK = 512
MLSCTX = 512
KEYPACKAGE = 1700*8
MLSHEADER = KEYPACKAGE # 273 * 8 - SIGSIZE + KEYPACKAGE


CONVERSION_FACTOR = 8 * 1024
CF_MB = 8 * 1024 * 1024

def log(x):
    return int(math.log(x, 2))

def calcMLSSenderSizes(rng, case):
    if case == 0: #Best
        data = [MLSHEADER + SIGSIZE + log(x)* (MLSCTI + MLSPK + MLSCTX) for x in rng]
    elif case == 1: #Worst
        data = [MLSHEADER +  SIGSIZE + x * (MLSCTI + MLSCTX) + log(x) * MLSPK for x in rng]
    else:
        raise Error
    return applyConversion(data)

def calcBGMSenderSizes(rng, case):
    if case == 0: #Best
        data = [MLSHEADER + SIGSIZE + DHCONST + log(x)* (MLSCTI + MLSPK) for x in rng]
    elif case == 1: #Worst
        data = [MLSHEADER +  SIGSIZE + DHCONST + x * (MLSCTI) + log(x) * MLSPK for x in rng]
    else:
        raise Error
    return applyConversion(data)

def calcMLSIndivSenderSizes(rng, case):
    if case == 0: #Best
        data = [MLSHEADER + log(x)* (SIGSIZE + MLSCTI + MLSPK + MLSCTX) for x in rng]
    elif case == 1: #Worst
        data = [MLSHEADER + x * (SIGSIZE + MLSCTI + MLSCTX) + log(x) * MLSPK for x in rng]
    else:
        raise Error
    return applyConversion(data)


def calcMLSRecSizes(rng, case):
    return calcMLSSenderSizes(rng, case)

def calcBGMRecSizeNoHuffman(rng, case):
    if case == 0: #Best
        data = applyConversion([(MLSHEADER + SIGSIZE + DHCONST + MLSCTI + log(x) * MLSPK + log(log(x)) * HASHSIZE) for x in rng])
    elif case == 1: #Worst
        data = applyConversion([(MLSHEADER + SIGSIZE + DHCONST + MLSCTI + log(x) * (MLSPK +HASHSIZE)) for x in rng])
    return data

def calcBGMRecSizes(rng, case, accType = "Hash"):
    if accType == "Hash":
        accSize = HASHSIZE
    elif accType == "RSA":
        accSize = RSAMOD
    elif accType == "Pairing":
        accSize = G1
    if case == 0: #Best
        data = [MLSHEADER + SIGSIZE + DHCONST + MLSCTI + log(x) * MLSPK + (log(x) if accType == "Pairing" else 2) * accSize for x in rng]
    elif case == 1: #Worst
        if accType == "RSA":
            data = [MLSHEADER + SIGSIZE + DHCONST + MLSCTI + log(x) * MLSPK + 2 * accSize for x in rng]
        elif accType == "Hash":
            data = [MLSHEADER + SIGSIZE + DHCONST + MLSCTI + log(x) * MLSPK + log(x) * accSize for x in rng]
        elif accType == "Pairing":
            data = [MLSHEADER + SIGSIZE + DHCONST + MLSCTI + log(x) * MLSPK + x * accSize for x in rng]
        else:
            raise Error()
    else:
        raise Error()
    return applyConversion(data)

def calcMLSIndivRecSizes(rng, case):
    data = [MLSHEADER + SIGSIZE + DHCONST + MLSCTI + log(x) * MLSPK for x in rng]
    return applyConversion(data)

def calcMLSRecSum(rng, case):
    tmp = calcMLSRecSizes(rng, case)
    data = applyConversion([x * y for (x,y) in zip(tmp, rng)], 1024)
    return data

def calcMLSIndivSum(rng, case):
    tmp = calcMLSIndivRecSizes(rng, case)
    data = applyConversion([x * y for (x,y) in zip(tmp, rng)], 1024)
    return data

def BGMRecSizeSum(n):
    sum = 0
    for depth in range(1, log(n)):
        sum += int(n / math.pow(2, depth)) * (MLSHEADER + SIGSIZE + DHCONST + MLSCTI + log(n) * MLSPK + depth * HASHSIZE)
    return sum

def calcBGMRecSum(rng, case, accType = "Hash"):
    if accType == "Hash":
        if case == 0: #Best
            data = applyConversion([BGMRecSizeSum(x) for x in rng], CF_MB)
        elif case == 1: #Worst
            tmp = calcBGMRecSizes(rng, case)
            data = applyConversion([x * y for (x,y) in zip(tmp, rng)], 1024)
    elif accType == "RSA":
        data = applyConversion([x * (MLSHEADER + SIGSIZE + DHCONST + MLSCTI + log(x) * MLSPK + 2 * RSAMOD) for x in rng], CF_MB)
    elif accType == "Pairing":
        if case == 0:
            data = applyConversion([x * (MLSHEADER + SIGSIZE + DHCONST + MLSCTI + log(x) * (MLSPK + G1)) for x in rng], CF_MB)
        else:
            data = applyConversion([x * (MLSHEADER + SIGSIZE + DHCONST + MLSCTI + log(x) * MLSPK + x * G1) for x in rng], CF_MB)
    else:
        raise Error()
    return data

def calcBGMRecSumNoHuffman(rng, case):
    if case == 0: #Best
        data = applyConversion([x * (MLSHEADER + SIGSIZE + DHCONST + MLSCTI + log(x) * MLSPK + log(log(x)) * HASHSIZE) for x in rng], CF_MB)
    elif case == 1: #Worst
        data = applyConversion([x * (MLSHEADER + SIGSIZE + DHCONST + MLSCTI + log(x) * (MLSPK + HASHSIZE)) for x in rng], CF_MB)
    return data
        
def applyConversion(data, factor = CONVERSION_FACTOR):
    return [x / factor for x in data]

if __name__ == "__main__":
    fig, ((ax1, ax4),(ax2, ax5), (ax3, ax6)) = plt.subplots(3, 2, sharex=True, figsize=(15, 15))

    
    rng = [int(math.pow(2, x)) for x in range(1,17)]
    mlsSenderDataBest = calcMLSSenderSizes(rng, 0)
    mlsSenderDataWorst = calcMLSSenderSizes(rng, 1)

    mlsSenderIndivBest = calcMLSIndivSenderSizes(rng, 0)
    mlsSenderIndivWorst = calcMLSIndivSenderSizes(rng, 1)

    bgmSenderBest = calcBGMSenderSizes(rng, 0)
    bgmSenderWorst = calcBGMSenderSizes(rng, 1)
    
    ax1.plot(rng, mlsSenderDataBest, "--", label = "ITK", color="red")    
    ax1.plot(rng, mlsSenderIndivBest, "--", label = "ITKI", color = "blue")
    ax1.plot(rng, bgmSenderBest, "--", label = "SAIK Huffman", color = "green")
    
    ax1.fill_between(rng, mlsSenderIndivBest, mlsSenderIndivWorst, color="blue", hatch = "-", fc = "#ffffff00", ec = "#0000ff7f")
    ax1.fill_between(rng, mlsSenderDataBest, mlsSenderDataWorst, color="red", hatch = "\\", fc = "#ffffff00", ec = "#ff00007f")
    ax1.fill_between(rng, bgmSenderBest, bgmSenderWorst, color="green", hatch = "//", fc = "#ffffff00", ec = "#00ff007f")

    ax1.plot(rng, mlsSenderDataWorst, "--", color="red")
    ax1.plot(rng, mlsSenderIndivWorst, "--", color = "blue")
    ax1.plot(rng, bgmSenderWorst, "--", color = "green")
    
    ax1.tick_params("both", reset = True)
    ax1.set_xlabel("#Users")
    ax1.set_ylabel("size in KB")
    ax1.set_yscale("log")
    ax1.title.set_text("a) Sender Bandwith")
    # ax1.legend(loc="upper left")

    senderPercBest = [x / y * 100 for (x,y) in zip(mlsSenderDataBest, bgmSenderBest)]
    print("Sender Best")
    print(senderPercBest)

    senderPercWorst = [x / y * 100 for (x,y) in zip(mlsSenderDataWorst, bgmSenderWorst)]
    print("Sender Worst")
    print(senderPercWorst)

    senderPercIndivBest = [x / y * 100 for (x,y) in zip(mlsSenderDataBest, mlsSenderIndivBest)]
    print("Sender Indiv Best")
    print(senderPercIndivBest)

    senderPercIndivWorst = [x / y * 100 for (x,y) in zip(mlsSenderDataWorst,mlsSenderIndivWorst)]
    print("Sender Indiv Worst")
    print(senderPercIndivWorst)
    
    mlsRecBest = calcMLSRecSizes(rng, 0)
    mlsRecWorst = calcMLSRecSizes(rng, 1)

    mlsRecIndivBest = calcMLSIndivRecSizes(rng, 0)
    mlsRecIndivWorst = calcMLSIndivRecSizes(rng, 1)

    bgmRecNoHuffBest = calcBGMRecSizeNoHuffman(rng, 0)
    bgmRecNoHuffWorst = calcBGMRecSizeNoHuffman(rng, 1)
    
    bgmRecBest = calcBGMRecSizes(rng, 0)
    bgmRecWorst = calcBGMRecSizes(rng, 1)

    percRecBest = [(x / y) * 100 for (x,y) in zip(mlsRecBest, bgmRecBest)]
    print("Receiver Best")
    print(percRecBest)

    percRecWorst = [(x / y) * 100 for (x,y) in zip(mlsRecWorst, bgmRecWorst)]
    print("Receiver Worst")
    print(percRecWorst)
    
    ax2.plot(rng, mlsRecBest, "--", label = "ITK", color="red")
    ax2.plot(rng, mlsRecWorst, "--", color="red")

    ax2.plot(rng, mlsRecIndivBest, "--", label = "ITKI", color = "blue")
    ax2.plot(rng, mlsRecIndivWorst, "--", color = "blue")

    ax2.plot(rng, bgmRecBest, "--", label = "SAIK Huffman", color = "green")
    ax2.plot(rng, bgmRecWorst, "--", color = "green")

    ax2.fill_between(rng, mlsRecIndivBest, mlsRecIndivWorst, color="blue", alpha = 0.3)
    ax2.fill_between(rng, mlsRecBest, mlsRecWorst, fc = "#ffffff00", ec="#ff00007f", hatch = "\\\\")
    ax2.fill_between(rng, bgmRecBest, bgmRecWorst, fc = "#ffffff00", ec="#00ff007f", hatch = "//")

    ax2.set_ylim([bgmRecBest[0] * 0.95, mlsRecBest[-1] * 1.01])
    # ax2.set_yscale("log")
    ax2.tick_params("both", reset = True)
    ax2.set_xlabel("#Users")
    ax2.set_ylabel("size in KB")
    ax2.title.set_text("b) Receiver Bandwidth")
    # ax2.legend(loc="lower right")

    mlsRecSumBest = calcMLSRecSum(rng, 0)
    mlsRecSumWorst = calcMLSRecSum(rng, 1)

    mlsRecIndivSumBest = calcMLSIndivSum(rng, 0)
    mlsRecIndivSumWorst = calcMLSIndivSum(rng, 1)

    bgmRecSumNoHuffBest = calcBGMRecSumNoHuffman(rng, 0)
    bgmRecSumNoHuffWorst = calcBGMRecSumNoHuffman(rng, 1)
    
    bgmRecSumBest = calcBGMRecSum(rng, 0)
    bgmRecSumWorst = calcBGMRecSum(rng, 1)
    
    ax3.plot(rng, mlsRecSumBest, "--", label = "ITK[2]", color="red")
    ax3.plot(rng, mlsRecSumWorst, "--", color="red")

    ax3.plot(rng, mlsRecIndivSumBest, "--", label = "ITKI[2]", color = "blue")
    ax3.plot(rng, mlsRecIndivSumWorst, "--", color = "blue")
    
    ax3.plot(rng, bgmRecSumBest, "--", label = "SAIK Huffman", color = "green")
    ax3.plot(rng, bgmRecSumWorst, "--", color = "green")

    ax3.fill_between(rng, mlsRecIndivSumBest, mlsRecIndivSumWorst, color="blue", alpha = 0.3)
    ax3.fill_between(rng, mlsRecSumBest, mlsRecSumWorst, hatch = "\\", fc = "#ffffff00", ec = "#ff00007f")
    ax3.fill_between(rng, bgmRecSumBest, bgmRecSumWorst, fc = "#ffffff00", ec="#00ff007f", hatch = "//")

    ax3.set_ylim([-2, bgmRecSumWorst[-1] * 1.01])
    ax3.set_xlabel("#Users")
    ax3.set_ylabel("size in MB")
    ax3.title.set_text("c) Server Cummulative Outgoing Bandwidth")
    # ax3.legend(loc="lower right")

    saikRecHashBest = calcBGMRecSizes(rng, 0, "Hash")
    saikRecHashWorst = calcBGMRecSizes(rng, 1, "Hash")
    saikRecHashHuffBest = calcBGMRecSizeNoHuffman(rng, 0)
    saikRecHashHuffWorst = calcBGMRecSizeNoHuffman(rng, 1)
    saikRecRSA = calcBGMRecSizes(rng, 0, "RSA")
    saikRecPairing = calcBGMRecSizes(rng, 0, "Pairing")
    saikRecPairingWorst = calcBGMRecSizes(rng, 1, "Pairing")

    saikRecHashSumBest = calcBGMRecSum(rng, 0, "Hash")
    saikRecHashSumWorst = calcBGMRecSum(rng, 1, "Hash")
    saikRecHashHuffSumBest = calcBGMRecSumNoHuffman(rng, 0)
    saikRecHashHuffSumWorst = calcBGMRecSumNoHuffman(rng, 1)

    percHuffMerkBest = [(x / y) * 100 for (x,y) in zip(saikRecHashHuffBest, saikRecHashBest)]
    print("HuffMerkle Best")
    print(percHuffMerkBest)
    
    ax4.plot(rng, saikRecHashSumBest, "--", label = "SAIK Huffman", color = "green")
    ax4.plot(rng, saikRecHashSumWorst, "--", color = "green")
    ax4.plot(rng, saikRecHashHuffSumBest, "--", label = "SAIK Merkle", color = "gray")
    ax4.plot(rng, saikRecHashHuffSumWorst, "--", color = "gray")

    ax4.fill_between(rng, saikRecHashSumBest, saikRecHashSumWorst, fc = "#ffffff00", ec="#00ff007f", hatch = "//")
    ax4.fill_between(rng, saikRecHashHuffSumBest, saikRecHashHuffSumWorst, hatch = "||", fc = "#ffffff00", ec = "#8080807f")

    ax4.tick_params("both", reset = True)
    ax4.set_xlabel("#Users")
    ax4.set_ylabel("size in MB") 
    ax4.title.set_text("d) Server Cummulative Outgoing Bandwidth: Huffman vs. Merkle")
    # ax4.legend(handles = [mpatches.Patch(label = "SAIK Huffman", fc = "#ffffff00", ec="#00ff007f", hatch = "//"),
    #                       mpatches.Patch(label = "SAIK Merkle", fc = "#ffffff00", ec="#ff00007f", hatch = "\\\\")], loc="lower right")
    
    
    ax5.plot(rng, saikRecHashBest, "--", label = "SAIK Huffman", color = "green")
    ax5.plot(rng, saikRecHashWorst, "--", color = "green")
    ax5.plot(rng, saikRecRSA, "--", label = "SAIK RSAAcc", color = "black")
    ax5.plot(rng, saikRecPairing, "--", label = "SAIK PAcc", color = "olive")
    ax5.plot(rng, saikRecPairingWorst, "--", color = "olive")
    
    ax5.fill_between(rng, saikRecHashBest, saikRecHashWorst, fc = "#ffffff00", ec="#00ff007f", hatch = "//")
    ax5.fill_between(rng, saikRecPairing, saikRecPairingWorst, fc = "#ffffff00", ec="#8080007f", hatch = "..")
    
    ax5.tick_params("both", reset = True)
    ax5.set_xlabel("#Users")
    ax5.set_ylabel("size in KB")
    ax5.set_ylim([saikRecHashBest[0]*.95, saikRecRSA[-1] * 1.05])
    ax5.title.set_text("e) Receiver Bandwidth: Different Accumulators")

    
    saikRecHashBestSum = calcBGMRecSum(rng, 0, "Hash")
    saikRecHashWorstSum = calcBGMRecSum(rng, 1, "Hash")
    saikRecHashNoHuffBestSum = calcBGMRecSumNoHuffman(rng, 0)
    saikRecHashNoHuffWorstSum = calcBGMRecSumNoHuffman(rng, 1)
    saikRecRSASum = calcBGMRecSum(rng, 0, "RSA")
    saikRecPairingSum = calcBGMRecSum(rng, 0, "Pairing")
    saikRecPairingSumWorst = calcBGMRecSum(rng, 1, "Pairing")

    ax6.plot(rng, saikRecHashBestSum, "--", label = "SAIK Huffman", color = "green")
    ax6.plot(rng, saikRecHashWorstSum, "--", color = "green")
    ax6.plot(rng, saikRecRSASum, "--", label = "SAIK RSAAcc[12]", color = "black")
    ax6.plot(rng, saikRecPairingSum, "--", label = "SAIK PAcc", color = "olive")
    ax6.plot(rng, saikRecPairingSumWorst, "--", color = "olive")

    
    ax6.fill_between(rng, saikRecHashBestSum, saikRecHashWorstSum, fc = "#ffffff00", ec="#00ff007f", hatch = "//")
    ax6.fill_between(rng, saikRecPairingSum, saikRecPairingSumWorst, fc = "#ffffff00", ec = "#8080007f", hatch = "..")
    
    ax6.tick_params("both", reset = True)
    ax6.set_xlabel("#Users")
    ax6.set_ylabel("size in MB")
    ax6.set_ylim([0, saikRecRSASum[-1] * 1.05])
    ax6.title.set_text("f) Server Cummulative Outgoing Bandwidth: Different Accumulators")

    green_patch = mpatches.Patch(label='SAIK Huffman',fc = "#ffffff00", ec="#00ff007f", hatch = "//")
    red_patch = mpatches.Patch(label='ITK[5]',fc = "#ffffff00", ec="#ff00007f", hatch = "\\\\")
    blue_patch = mpatches.Patch(label='ITKI[5]',fc = "#ffffff00", ec="#0000ff7f", hatch = "--")
    grey_patch = mpatches.Patch(label='SAIK Merkle', fc = "#ffffff00", ec = "#808080", hatch = "||")
    olive_patch = mpatches.Patch(label='SAIK PAcc[49]', fc = "#ffffff00", ec = "#808080", hatch = "..")
    
    
    handles, labels = ax6.get_legend_handles_labels()
    fig.legend(handles = [green_patch, red_patch, blue_patch, grey_patch, olive_patch, handles[1]],loc='upper center',
               ncol=2, bbox_to_anchor=(.5, .95))

    
    handles[0] = mpatches.Patch(label = "SAIK Huffman", hatch = "//", fc = "#ffffff00", ec="#00ff007f")
    handles[2] = mpatches.Patch(label = "SAIK PAcc", hatch = "\\\\", fc = "#ffffff00", ec="#8080007f")
    # fig.legend(handles, labels, bbox_to_anchor=(1.01, .33), loc='lower right')
    # fig.legend(handles, labels, loc='lower right')
    
    
    plt.savefig("Final_Figures")
