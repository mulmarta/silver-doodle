#!/bin/python3

import json
import numpy as np
import matplotlib.pyplot as plt
from pyexcel_ods import get_data
import pdb
import math

HASHSIZE = 256
G1 = 256
G2 = 512
GT = 3072
AEADSIZE = 384
SIGSIZE = HASHSIZE + G1
RSAMOD = 3072

KYBERCT0 = 11264
KYBERCTI = 1280
KYBERPK = 11520

DHCT0 = 0
DHCTI = 0
DHCONST = G1
DHKEY = G1

MLSCT0 = 0
MLSCTI = 512
MLSPK = G1

CONVERSION_FACTOR = 8 * 1024


def sigFactor(x, isBestCase, indivSign, q = 2):
    if indivSign and isBestCase:
        return math.ceil(math.log(x, q)) - 1
    elif indivSign and not isBestCase:
        return x - 1
    else:
        return 1

def accFactor(n, q, isBestCase, accType, huffman = False):
    if accType != "Hash":
        return 1
    t = int(math.ceil(math.log(n, q)) - 1)
    if huffman:
        return 2
    if isBestCase:
        t = int(math.ceil(math.log(t, q)) - 1)
    return t

def accFactorSum(n, q, isBestCase, accType, huffman):
    if accType != "Hash":
        return n
    t = int(math.ceil(math.log(n,q,)) -1)
    if huffman:
        if isBestCase:
            return int((2 * t - 2 - int(math.ceil(math.log(t,q)) -1)) / t)
        else:
            return int((2 * n - t -2)/n)
    else:
        if isBestCase:
            return int(math.ceil(math.log(t,q)) -1)
        else:
            return t

def calcRecSize(numUsers, treeDeg, ct0, cti, constSize, numCti, pkSize, accSize, isBestCase = False, isMLS = False):
    t = int(math.ceil(math.log(numUsers, treeDeg)) - 1)
    symEnc = ((t if isBestCase else numUsers) - 1 if isMLS else 1) * AEADSIZE
    asymEnc = ct0 + numCti * cti + constSize
    aux = t*pkSize
    return (symEnc + SIGSIZE + asymEnc + aux + accSize) / CONVERSION_FACTOR

def calcRecSizeSum(numUsers, treeDeg, ct0, cti, constSize, numCti, pkSize, accSize, isBestCase, huffman = False, isMLS = False):
    t = int(math.ceil(math.log(numUsers, treeDeg)) - 1)
    symEnc = ((t if isBestCase else numUsers) - 1 if isMLS else 1) * AEADSIZE
    asymEnc = ct0 + numCti * cti + constSize
    aux = t*pkSize
    return (symEnc + SIGSIZE + asymEnc + aux + accSize) / CONVERSION_FACTOR

def calcSenderSize(numUsers, treeDeg, aeadSize, ct0, cti, constSize, pkSize, sigSize, accSize, isBestCase = True, ismKEM = False):
    t = int(math.ceil(math.log(numUsers, treeDeg)) - 1)
    if isBestCase:
        encSize = t * (ct0 + (treeDeg - 1) * cti) + constSize
    else:
        encSize = t * ct0 + (numUsers - 1) * cti + constSize
    auxSize = t * pkSize
    if ismKEM:
        symSize = t * aeadSize
    elif isBestCase:
        symSize = t * (treeDeg - 1) * aeadSize
    else:
        symSize = (numUsers - 1) * aeadSize
    return (encSize + sigSize + auxSize + symSize + accSize) / CONVERSION_FACTOR

def plotMLSReceiver(ax, q, xaxis, label, isBestCase = True, indivSign = False):
    data = [calcRecSize(n, q, MLSCT0, MLSCTI, 0, 1 if indivSign else int(math.ceil(math.log(n, q)) - 1) if isBestCase else n, MLSPK, 0, isBestCase = isBestCase, isMLS = not indivSign) for n in xaxis]
    ax.plot(xaxis, data, "--", label = label)

def plotMLSSender(ax, q, xaxis, label, isBestCase, indivSign = False):
    data = [calcSenderSize(n, q, AEADSIZE, MLSCT0, MLSCTI, 0, MLSPK, sigFactor(n, isBestCase, indivSign) * SIGSIZE, 0, isBestCase) for n in xaxis]
    ax.plot(xaxis, data, "--", label = label)

def plotMLSSenderwithKyber(ax, q, xaxis, label, isBestCase, indivSign):
    data = [calcSenderSize(n, q, AEADSIZE, 0, KYBERCT0 + KYBERCTI, 0, KYBERPK, sigFactor(n, isBestCase, indivSign) * SIGSIZE, 0, isBestCase) for n in xaxis]
    ax.plot(xaxis, data, "--", label=label)

def plotMLSReceiverwithKypber(ax, q, xaxis, label, isBestCase = False):
    data = [calcRecSize(n, q, KYBERCT0, KYBERCTI, 0, n, KYBERPK, 0, isBestCase, True) for n in xaxis]
    ax.plot(xaxis, data, "--", label=label)

    
def plotMkyberSender(ax, q, xaxis, label, accType, isBestCase, indivSign = False):
    if indivSign:
        accSize =  0
    if accType == "RSA":
        accSize = RSAMOD
    if accType == "Pairing":
        accSize = G1
    if accType == "Hash":
        accSize =  0
    data = [calcSenderSize(n, q, AEADSIZE, KYBERCT0, KYBERCTI, 0, KYBERPK, sigFactor(n, isBestCase, indivSign) * SIGSIZE, accSize, isBestCase, True) for n in xaxis]
    # for x, y in zip(xaxis, data):
    #     print("q: {}; x: {}: y: {}".format(q, x, y))
    
    ax.plot(xaxis, data, "--", label=label)
    
def plotMkyberRec(ax, q, xaxis, label, accType, isBestCase):
    if accType == "RSA":
        accSize = 2*RSAMOD
    elif accType == "Pairing":
        accSize = G1 + G2
    elif accType == "Hash":
        accSize =  256
    #pdb.set_trace()
    data = [calcRecSize(n, q, KYBERCT0, KYBERCTI, 0, 1, KYBERPK, accFactor(n, q, isBestCase, accType) * accSize, isBestCase) for n in xaxis]
    # for x, y in zip(xaxis, data):
    #     print("x: {}: y: {}".format(x,y))
    ax.plot(xaxis, data, "--", label=label)
    
def plotDHSender(ax, q, xaxis, label, accType, isBestCase, indivSign):
    if indivSign:
        accSize =  0
    elif accType == "RSA":
        accSize = RSAMOD
    elif accType == "Pairing":
        accSize = G1
    elif accType == "Hash":
        accSize =  0
    else:
        raise TypeError
    data = [calcSenderSize(n, q, AEADSIZE, DHCT0, DHCTI, DHCONST, DHKEY, sigFactor(n, isBestCase, indivSign) * SIGSIZE, accSize, isBestCase, False) for n in xaxis]
    ax.plot(xaxis, data, "--", label=label)
    # for x, y in zip(xaxis, data):
    #     print("x: {}: y: {}".format(x,y))

def plotDHRec(ax, q, xaxis, accType, label, isBestCase = True, huffman = False):
    if accType == "RSA":
        accSize = RSAMOD
    if accType == "Pairing":
        accSize = G1
    if accType == "Hash":
        accSize =  256
    data = [calcRecSize(n, q, DHCT0, DHCTI, DHCONST, 1, DHKEY, accFactor(n, q, isBestCase, accType, huffman) * accSize, isBestCase) for n in xaxis]
    ax.plot(xaxis, data, "--", label=label)

def plotFromOds(path, table = "Tabelle2"):
    data = get_data(path)[table]
    xr = data[0][1:]
    figName = ""
    plotEmpty = True
    for l in data[1:]:
        if l == []:
            if not plotEmpty:
                plt.xlabel("#Users")
                plt.ylabel("Size in bits")
                plt.legend(loc="upper left")
                plt.savefig(figName + ".png")
                plotEmpty = True
            continue
        if len(l) == 1:
            if not plotEmpty:
                plt.savefig(figName + ".png")
            plotEmpty = False
            figName = str(l[0])
            plt.figure(figName, figsize = (10, 10), dpi = 300)
            continue            
        name = l[0]
        if name.startswith("mKEM"):
            linetype = "s-"
        elif name.startswith("NIKE"):
            linetype = "^-"
        else:
            linetype = "--"
        plt.plot(xr[:len(l)-1], l[1:], linetype, label = name)

def newFig():
    fig, ax = plt.subplots(figsize = (10,10), dpi = 300)
    ax.ticklabel_format(style = "plain")
    plt.ylabel("Size in KiB")
    plt.xlabel("#Users")
    return fig, ax

def plotMLSReceiverSum(ax, q, xaxis, label, isBestCase, indivSign = False):
    data = [(math.log(n,q) if isBestCase else n) * calcRecSize(n, q, MLSCT0, MLSCTI, 0, 1 if indivSign else int(math.ceil(math.log(n, q)) - 1) if isBestCase else n, MLSPK, 0, not indivSign) for n in xaxis]
    ax.plot(xaxis, data, "--", label = label)

def plotDHRecSum(ax, q, xr, label, accType, isBestCase, huffman = False):
    
    if accType == "RSA":
        accSize = RSAMOD
    if accType == "Pairing":
        accSize = G1
    if accType == "Hash":
        accSize =  256
    data = [(math.log(n, q) if isBestCase else n) * calcRecSize(n, q, DHCT0, DHCTI, DHCONST, 1, DHKEY, accFactorSum(n, q, isBestCase, accType, huffman) * accSize) for n in xr]
    ax.plot(xr, data, "--", label=label)


def plotReceiverSum(xr, q, isBestCase):
    fix, ax = newFig()
    plotMLSReceiverSum(ax, q= q, xaxis = xr, label = "MLS", isBestCase = isBestCase, indivSign = False)
    plotMLSReceiverSum(ax, q= q, xaxis = xr, label = "MLS indiv. Sign", isBestCase = isBestCase, indivSign = True)    
    plotDHRecSum(ax, q, xr, label = "DH HAcc", accType = "Hash", isBestCase = isBestCase, huffman = False)
    plotDHRecSum(ax, q, xr, label = "DH HAcc Huffman", accType = "Hash", isBestCase = isBestCase, huffman = True)
    plt.legend(loc="upper left")
    plt.savefig("Receiver_Sum_{}_{}.png".format("Best" if isBestCase else "Worst", q))

def plotAllSender(xr, q, isBestCase):
    fig, ax = newFig()
    plotMLSSender(ax, q = q, xaxis = xr, label = "MLS", isBestCase = isBestCase, indivSign = False)
    plotMLSSender(ax, q = q, xaxis = xr, label = "MLS indiv. Sign", isBestCase = isBestCase, indivSign = True)
    # plotMLSSenderwithKyber(ax, q = q, xaxis = xr, label = "MLS with Kyber", isBestCase = isBestCase, indivSign = False)
    plotDHSender(ax, q, xr, label = "DH HAcc", accType = "Hash", isBestCase = isBestCase, indivSign = False)
    # plotMkyberSender(ax, q = q, xaxis = xr, label = "MKyber", accType = "Hash", isBestCase = isBestCase, indivSign = False)
    plt.legend(loc="upper left")
    plt.savefig("Sender_All_{}_{}.png".format("Best" if isBestCase else "Worst", q))

def plotAllReceiver(xr, q, isBestCase):
    fig, ax = newFig()
    plotMLSReceiver(ax, q, xr, label = "MLS", isBestCase = isBestCase, indivSign = False)
    plotMLSReceiver(ax, q, xr, label = "MLS indivSign", isBestCase = isBestCase, indivSign = True)
    # plotMLSReceiverwithKypber(ax, q, xr, label = "MLS Kyber")
    plotDHRec(ax, q, xr, label = "DH HAcc", accType = "Hash", isBestCase = isBestCase)
    plotDHRec(ax, q, xr, label = "DH HAcc Huffman (avg.)", accType = "Hash", isBestCase = isBestCase, huffman = True)
    # plotMkyberRec(ax, q, xr, label = "MKyber HAcc", accType = "Hash", isBestCase = isBestCase)
    plt.legend(loc="upper left")
    plt.savefig("Receiver_All_{}_{}".format("Best" if isBestCase else "Worst", q))

def plotAllDHSender(xr, q, isBestCase):
    fig, ax = newFig()
    plotDHSender(ax, q, xr, label = "DH HAcc", accType = "Hash", isBestCase = isBestCase, indivSign = False)
    plotDHSender(ax, q, xr, label = "DH PAcc", accType = "Pairing", isBestCase = isBestCase, indivSign = False)
    plotDHSender(ax, q, xr, label = "DH RSAAcc", accType = "RSA", isBestCase = isBestCase, indivSign = False)
    plt.legend(loc="upper left")
    plt.savefig("Sender_DH_{}_{}".format("Best" if isBestCase else "Worst", q))

def plotAllMKyberSender(xr, q, isBestCase):
    fig, ax = newFig()
    plotMkyberSender(ax, q, xr, label = "Mkyber HAcc", accType = "Hash", isBestCase = isBestCase)
    plotMkyberSender(ax, q, xr, label = "Mkyber PAcc", accType = "Pairing", isBestCase = isBestCase)
    plotMkyberSender(ax, q, xr, label = "Mkyber RSAAcc", accType = "RSA", isBestCase = isBestCase)
    plt.legend(loc="upper left")
    plt.savefig("Sender_Mkyber_{}_{}".format("Best" if isBestCase else "Worst", q))

def plotDHTrees(xr, isBestCase):
    fig, ax = newFig()
    for q in range(2, 6, 1):
        plotDHSender(ax, q, xr, label = "DH q = {}".format(q), accType = "RSA", isBestCase = isBestCase, indivSign = False)
    plt.legend(loc="upper left")
    plt.savefig("Sender_DH_TreeDeg_{}".format("Best" if isBestCase else "Worst"))

    
def plotMKyberTrees(xr, isBestCase):
    fig, ax = newFig()
    for q in range(10, 16, 1):
        plotMkyberSender(ax, q, xr, label = "Mkyber q = {}".format(q), accType = "Hash", isBestCase = isBestCase)
    plt.legend(loc="upper left")
    plt.savefig("Sender_Mkyber_TreeDeg_{}".format("Best" if isBestCase else "Worst"))

def plotNewVariantsSender(xr, q, isBestCase):
    fix, ax = newFig()

    plotDHSender(ax, q, xr, label = "DH HAcc", accType = "Hash", isBestCase = isBestCase, indivSign = False)
    plotDHSender(ax, q, xr, label = "DH PAcc", accType = "Pairing", isBestCase = isBestCase, indivSign = False)
    plotDHSender(ax, q, xr, label = "DH RSAAcc", accType = "RSA", isBestCase = isBestCase, indivSign = False)

    # plotMkyberSender(ax, q, xr, label = "MKyber HAcc", accType = "Hash", isBestCase = isBestCase)
    # plotMkyberSender(ax, q, xr, label = "MKyber PAcc", accType = "Pairing", isBestCase = isBestCase)
    # plotMkyberSender(ax, q, xr, label = "MKyber RSAAcc", accType = "RSA", isBestCase = isBestCase)
        
    plt.legend(loc="upper left")
    plt.savefig("Sender_New_{}_q{}".format("Best" if isBestCase else "Worst", q))

def plotComputation(xr, q, isBestCase):
    fix, ax = newFig()

    dataMLS = [2 + 2*(math.log(n,2) if isBestCase else n) for n in xr]
    dataMLSIndiv = [4*(math.log(n,2) if isBestCase else n) for n in xr]
    dataDH = [2 + (math.log(n,2) if isBestCase else n) for n in xr]
    ax.plot(xr, dataMLS, "--", label="MLS")
    ax.plot(xr, dataMLSIndiv, "-d", label="MLS Indiv. Sign")
    ax.plot(xr, dataDH, "-s", label="DH")
    plt.ylabel("#public key operations")
    plt.legend(loc="upper left")
    plt.savefig("Computation_{}".format("Best" if isBestCase else "Worst"))

if __name__ == "__main__":
    #plotFromOds("../Sizes.ods")
    # plotMLSSender(None, 2, "MLS", False, True)
    # plotMkyberSender(None, 2, "mKyber", "Hash", True, False)
    xr = range(500, 50001, 500)
    plotAllSender(xr, 2, True)
    plotAllSender(xr, 2, False)
    plotAllReceiver(xr, 2, True)
    plotAllReceiver(xr, 2, False)

    plotAllDHSender(xr, 2, True)
    # plotAllMKyberSender(xr, 2, True)
    plotAllDHSender(xr, 2, False)
    # plotAllMKyberSender(xr, 2, False)

    # plotMKyberTrees(xr, True)
    plotDHTrees(xr, True)

    plotReceiverSum(xr, 2, True)
    plotReceiverSum(xr, 2, False)
    
    plotNewVariantsSender(xr, 2, True)
    plotNewVariantsSender(xr, 2, False)

    plotComputation(xr, 2, True)
    plotComputation(xr, 2, False)
