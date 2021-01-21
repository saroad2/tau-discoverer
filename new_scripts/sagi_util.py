import os
import ROOT
import copy
import math


######################################
######## Define helper functions
######################################

# def ComparePt(TLorentzVector a, TLorentzVector b):
def ComparePt(a, b):
    return a.Pt() > b.Pt()


# take second element for sort
def takeSecond(elem):
    return elem[1]


# fix the eta to be the one from the max cell in EM2
def findFixedEta(myTOB, myTOBVector):
    cellEM2max = myTOB.getCentralL2Cell()
    if cellEM2max == 0:
        FixedEta = myTOBVector.Eta() - (0.025 + 0.0125)
    elif cellEM2max == 1:
        FixedEta = myTOBVector.Eta() - 0.0125
    elif cellEM2max == 2:
        FixedEta = myTOBVector.Eta() + 0.0125
    elif cellEM2max == 3:
        FixedEta = myTOBVector.Eta() + (0.025 + 0.0125)
    return FixedEta


# Run-II isolation requirement
def jFEXIso(myTOB):
    cutOff = 2.0 + 0.1 * (myTOB.tauClus())
    return myTOB.emIsol() < cutOff


def Run2Iso(myTOB):
    cutOff = 2.0 + 0.1 * (myTOB.ppmTauClus())
    return myTOB.ppmEMIsol() < cutOff


# Simple vector sum
def myLayerSum(myVector, noisecut=0.0):
    mySum = 0.0
    for entry in myVector:
        if entry > noisecut:
            mySum += entry
    return mySum


# Cluster Depth
def clusDepth(myTOB):
    myCells = createCellLists(myTOB)
    # treat Had as EM3 to reduce fluctuations
    layerDepth = [5, 55, 268, 499, 800]

    myDepth = 0
    for layer in range(5):
        myDepth += layerDepth[layer] * myLayerSum(myCells[layer])

    if myTOB.largeTauClus() < 0.001:
        return 0
    else:
        return myDepth / myTOB.largeTauClus()


# Cluster Density
def clusDens(myTOB):
    myCells = createCellLists(myTOB)
    # treat Had as EM3 to reduce fluctuations
    layerVolume = [
        147.0 * 147.0 * 10.0,
        147.0 * 36.8 * 90.0,
        147.0 * 36.8 * 336.0,
        147.0 * 147.0 * 42.0,
        147.0 * 147.0 * 42.0,
    ]

    myDens = 0
    for layer in range(5):
        for cell in myCells[layer]:
            myDens += (cell * cell * 1000 * 1000) / layerVolume[layer]

    if myTOB.largeTauClus() < 0.001:
        return 0
    else:
        return math.log(myDens / myTOB.largeTauClus())


# Cluster Width
def clusWidth(myTOB):
    myCells = createCellLists(myTOB)
    layerEta = [0.1, 0.025, 0.025, 0.1, 0.1]
    layerOffset = [1, 4, 4, 1, 1]
    layerCells = [i * 3 for i in layerOffset]

    myWidth = 0
    for layer in range(5):
        eta = -0.1
        phi = -0.1
        for cell in range(layerCells[layer]):
            myWidth += math.sqrt((eta * eta) + (phi * phi)) * myCells[layer][cell]
            eta += layerEta[layer]
            if eta > 0.101:
                phi += 0.1
                eta = -0.1

    if myTOB.largeTauClus() < 0.001:
        return 0
    else:
        return myWidth / myTOB.largeTauClus()


# Maximum energy in X cells
def maxCells(myVector, nCells):
    newVector = copy.deepcopy(myVector)
    newVector = sorted(newVector, reverse=True)
    return newVector[0:nCells]


# ET-weighted Distance
def weightDist(myVector):
    runningSum = 0
    for i in range(0, 36):
        index1 = i / 12
        index2 = float(i - (index1 * 12))
        runningSum += (
            math.sqrt(pow(abs(index1 - 1) * 0.1, 2) + pow(abs(index2 - 5.5) * 0.025, 2))
            * myVector[i]
        )

    return runningSum


# Fill ET maps
def createCellLists(myTOB):
    layerOffset = [1, 4, 4, 1, 1]
    layerCells = [i * 3 for i in layerOffset]

    allCells = []

    # iterate over layer
    for l in range(5):
        myLayer = []
        # iterate over phi
        for i in range(3):
            # iterate over eta
            for j in range(layerCells[l]):
                myLayer += [myTOB.getEnergy(l + 2, j + layerOffset[l], i + 1)]
        allCells += [myLayer]
    return allCells


def createBigCellLists(myTOB):
    layerOffset = [1, 4, 4, 1, 1]
    layerCells = [i * 3 for i in layerOffset]

    allCells = []

    # iterate over layer
    for l in range(5):
        myLayer = []
        # iterate over phi
        for i in range(3):
            # iterate over eta
            if l == 1 or l == 2:
                for j in range(0, layerCells[l], 2):
                    myLayer += [
                        myTOB.getEnergy(l + 2, j + layerOffset[l], i + 1)
                        + myTOB.getEnergy(l + 2, j + 1 + layerOffset[l], i + 1)
                    ]
            else:
                for j in range(layerCells[l]):
                    myLayer += [myTOB.getEnergy(l + 2, j + layerOffset[l], i + 1)]
        allCells += [myLayer]
    return allCells


def hottestSum(myTOB, nCells):
    # get ET map
    allCells = createCellLists(myTOB)
    # print allCells
    # Make running sum over all 5 layers
    mySum = 0
    for i in range(5):
        mySum += myLayerSum(maxCells(allCells[i], nCells[i]))

    return mySum


def bigCluster(myTOB):
    allCells = createBigCellLists(
        myTOB
    )  # same as createCellList but EM1 and EM2 are now 6 length arrays with the energy of 4+5, 6+7, 8+9, 10+11, 12+13, 14+15

    # for phi=1 I'll add up the energy of the 2 main big cells (8,9,10,11 in fine gran), then also that window shifted to one side (6,7,8,9) and to the other ((10,11,12,13)
    midPhiEM1 = [
        allCells[1][7] + allCells[1][8],
        allCells[1][8] + allCells[1][9],
        allCells[1][9] + allCells[1][10],
    ]
    midPhiEM2 = [
        allCells[2][7] + allCells[2][8],
        allCells[2][8] + allCells[2][9],
        allCells[2][9] + allCells[2][10],
    ]

    # for phi=0 and phi=1 i'll just loop over the two blocks sum energy (6,7 and 8,9 and 10,11 and 12,13)
    upPhiEM1 = [allCells[1][1], allCells[1][2], allCells[1][3], allCells[1][4]]
    lowPhiEM1 = [allCells[1][13], allCells[1][14], allCells[1][15], allCells[1][16]]
    upPhiEM2 = [allCells[2][1], allCells[2][2], allCells[2][3], allCells[2][4]]
    lowPhiEM2 = [allCells[2][13], allCells[2][14], allCells[2][15], allCells[2][16]]

    energyEM1EM2_mid = [
        midPhiEM1[0] + midPhiEM2[0],
        midPhiEM1[1] + midPhiEM2[1],
        midPhiEM1[2] + midPhiEM2[2],
    ]
    energyEM1EM2_low = [
        lowPhiEM1[0] + lowPhiEM2[0],
        lowPhiEM1[1] + lowPhiEM2[1],
        lowPhiEM1[2] + lowPhiEM2[2],
        lowPhiEM1[3] + lowPhiEM2[3],
    ]
    energyEM1EM2_up = [
        upPhiEM1[0] + upPhiEM2[0],
        upPhiEM1[1] + upPhiEM2[1],
        upPhiEM1[2] + upPhiEM2[2],
        upPhiEM1[3] + upPhiEM2[3],
    ]

    energyEM1EM2_mid.sort(reverse=True)
    energyEM1EM2_low.sort(reverse=True)
    energyEM1EM2_up.sort(reverse=True)

    energyEM1EM2 = energyEM1EM2_low[0] + energyEM1EM2_mid[0] + energyEM1EM2_up[0]

    # etaMax_PS = findMax(3,myTOB)[0]
    # phiMax_PS = findMax(2,myTOB)[1]
    # E_PS  = shapeTLV_fineGran(myTOB,2,etaMax_PS,phiMax_PS)
    # etaMax_EM3 = findMax(5,myTOB)[0]
    # phiMax_EM3 = findMax(5,myTOB)[1]
    # E_EM3  = shapeTLV_fineGran(myTOB,5,etaMax_EM3,phiMax_EM3)
    # etaMax_HAD = findMax(6,myTOB)[0]
    # phiMax_HAD = findMax(6,myTOB)[1]
    # E_HAD  = shapeTLV_fineGran(myTOB,6,etaMax_HAD,phiMax_HAD)

    # allE = energyEM1EM2+E_PS+E_EM3+E_HAD
    energyPS = allCells[0]
    energyEM3 = allCells[3]
    energyHAD = allCells[4]

    energyPS.sort(reverse=True)
    energyEM3.sort(reverse=True)
    energyHAD.sort(reverse=True)

    energyPSEM3HAD = (
        energyPS[0]
        + energyPS[1]
        + energyPS[2]
        + energyEM3[0]
        + energyEM3[1]
        + energyHAD[0]
        + energyHAD[1]
        + energyHAD[2]
    )

    E_PS = energyPS[0] + energyPS[1] + energyPS[2]
    E_EM3 = energyEM3[0] + energyEM3[1]
    E_HAD = energyHAD[0] + energyHAD[1] + energyHAD[2]

    allE = energyEM1EM2 + energyPSEM3HAD

    return [allE, energyEM1EM2, E_PS, E_EM3, E_HAD]
    # return allE


def bigCluster6(myTOB):
    allCells = createCellLists(myTOB)

    # for phi=1 I'll add up the energy of the 6 centras supercells (15,16,17,18,19,20) then also that window shifted two SC to one side (13,14,15,16,17,18) and to the other (17,18,19,20,21,22)
    midPhiEM1 = [
        allCells[1][13]
        + allCells[1][14]
        + allCells[1][15]
        + allCells[1][16]
        + allCells[1][17]
        + allCells[1][18],
        allCells[1][15]
        + allCells[1][16]
        + allCells[1][17]
        + allCells[1][18]
        + allCells[1][19]
        + allCells[1][20],
        allCells[1][17]
        + allCells[1][18]
        + allCells[1][19]
        + allCells[1][20]
        + allCells[1][21]
        + allCells[1][22],
    ]
    midPhiEM2 = [
        allCells[2][13]
        + allCells[2][14]
        + allCells[2][15]
        + allCells[2][16]
        + allCells[2][17]
        + allCells[2][18],
        allCells[2][15]
        + allCells[2][16]
        + allCells[2][17]
        + allCells[2][18]
        + allCells[2][19]
        + allCells[2][20],
        allCells[2][17]
        + allCells[2][18]
        + allCells[2][19]
        + allCells[2][20]
        + allCells[2][21]
        + allCells[2][22],
    ]

    # for phi=0 and phi=1 i'll just loop over the two blocks sum energy (6,7 and 8,9 and 10,11 and 12,13)
    upPhiEM1 = [
        allCells[1][1] + allCells[1][2],
        allCells[1][3] + allCells[1][4],
        allCells[1][5] + allCells[1][6],
        allCells[1][7] + allCells[1][8],
        allCells[1][9] + allCells[1][10],
    ]
    upPhiEM2 = [
        allCells[2][1] + allCells[2][2],
        allCells[2][3] + allCells[2][4],
        allCells[2][5] + allCells[2][6],
        allCells[2][7] + allCells[2][8],
        allCells[2][9] + allCells[2][10],
    ]
    lowPhiEM1 = [
        allCells[1][1 + 24] + allCells[1][2 + 24],
        allCells[1][3 + 24] + allCells[1][4 + 24],
        allCells[1][5 + 24] + allCells[1][6 + 24],
        allCells[1][7 + 24] + allCells[1][8 + 24],
        allCells[1][9 + 24] + allCells[1][10 + 24],
    ]
    lowPhiEM2 = [
        allCells[2][1 + 24] + allCells[2][2 + 24],
        allCells[2][3 + 24] + allCells[2][4 + 24],
        allCells[2][5 + 24] + allCells[2][6 + 24],
        allCells[2][7 + 24] + allCells[2][8 + 24],
        allCells[2][9 + 24] + allCells[2][10 + 24],
    ]

    energyEM1EM2_mid = [
        midPhiEM1[0] + midPhiEM2[0],
        midPhiEM1[1] + midPhiEM2[1],
        midPhiEM1[2] + midPhiEM2[2],
    ]
    energyEM1EM2_low = [
        lowPhiEM1[0] + lowPhiEM2[0],
        lowPhiEM1[1] + lowPhiEM2[1],
        lowPhiEM1[2] + lowPhiEM2[2],
        lowPhiEM1[3] + lowPhiEM2[3],
        lowPhiEM1[4] + lowPhiEM2[4],
    ]
    energyEM1EM2_up = [
        upPhiEM1[0] + upPhiEM2[0],
        upPhiEM1[1] + upPhiEM2[1],
        upPhiEM1[2] + upPhiEM2[2],
        upPhiEM1[3] + upPhiEM2[3],
        upPhiEM1[4] + upPhiEM2[4],
    ]

    energyEM1EM2_mid.sort(reverse=True)
    energyEM1EM2_low.sort(reverse=True)
    energyEM1EM2_up.sort(reverse=True)

    energyEM1EM2 = energyEM1EM2_low[0] + energyEM1EM2_mid[0] + energyEM1EM2_up[0]

    # etaMax_PS = findMax(2,myTOB)[0]
    # phiMax_PS = findMax(2,myTOB)[1]
    # E_PS  = shapeTLV_fineGran(myTOB,2,etaMax_PS,phiMax_PS)
    # etaMax_EM3 = findMax(5,myTOB)[0]
    # phiMax_EM3 = findMax(5,myTOB)[1]
    # E_EM3  = shapeTLV_fineGran(myTOB,5,etaMax_EM3,phiMax_EM3)
    # etaMax_HAD = findMax(6,myTOB)[0]
    # phiMax_HAD = findMax(6,myTOB)[1]
    # E_HAD  = shapeTLV_fineGran(myTOB,6,etaMax_HAD,phiMax_HAD)
    #
    # allE = energyEM1EM2+E_PS+E_EM3+E_HAD

    energyPS = allCells[0]
    energyEM3 = allCells[3]
    energyHAD = allCells[4]

    energyPS.sort(reverse=True)
    energyEM3.sort(reverse=True)
    energyHAD.sort(reverse=True)

    energyPSEM3HAD = (
        energyPS[0]
        + energyPS[1]
        + energyPS[2]
        + energyEM3[0]
        + energyEM3[1]
        + energyHAD[0]
        + energyHAD[1]
        + energyHAD[2]
    )

    E_PS = energyPS[0] + energyPS[1] + energyPS[2]
    E_EM3 = energyEM3[0] + energyEM3[1]
    E_HAD = energyHAD[0] + energyHAD[1] + energyHAD[2]

    allE = energyEM1EM2 + energyPSEM3HAD

    return [allE, energyEM1EM2, E_PS, E_EM3, E_HAD]
    # return allE


def bigClusterCal(myTOB):
    allCells = createBigCellLists(
        myTOB
    )  # same as createCellList but EM1 and EM2 are now 6 length arrays with the energy of 4+5, 6+7, 8+9, 10+11, 12+13, 14+15

    # for phi=1 I'll add up the energy of the 2 main big cells (8,9,10,11 in fine gran), then also that window shifted to one side (6,7,8,9) and to the other ((10,11,12,13)
    midPhiEM1 = [
        allCells[1][7] + allCells[1][8],
        allCells[1][8] + allCells[1][9],
        allCells[1][9] + allCells[1][10],
    ]
    midPhiEM2 = [
        allCells[2][7] + allCells[2][8],
        allCells[2][8] + allCells[2][9],
        allCells[2][9] + allCells[2][10],
    ]

    # for phi=0 and phi=1 i'll just loop over the two blocks sum energy (6,7 and 8,9 and 10,11 and 12,13)
    upPhiEM1 = [allCells[1][1], allCells[1][2], allCells[1][3], allCells[1][4]]
    lowPhiEM1 = [allCells[1][13], allCells[1][14], allCells[1][15], allCells[1][16]]
    upPhiEM2 = [allCells[2][1], allCells[2][2], allCells[2][3], allCells[2][4]]
    lowPhiEM2 = [allCells[2][13], allCells[2][14], allCells[2][15], allCells[2][16]]

    energyEM1EM2_mid = [
        midPhiEM1[0] + midPhiEM2[0],
        midPhiEM1[1] + midPhiEM2[1],
        midPhiEM1[2] + midPhiEM2[2],
    ]
    energyEM1EM2_low = [
        lowPhiEM1[0] + lowPhiEM2[0],
        lowPhiEM1[1] + lowPhiEM2[1],
        lowPhiEM1[2] + lowPhiEM2[2],
        lowPhiEM1[3] + lowPhiEM2[3],
    ]
    energyEM1EM2_up = [
        upPhiEM1[0] + upPhiEM2[0],
        upPhiEM1[1] + upPhiEM2[1],
        upPhiEM1[2] + upPhiEM2[2],
        upPhiEM1[3] + upPhiEM2[3],
    ]

    energyEM1EM2_mid.sort(reverse=True)
    energyEM1EM2_low.sort(reverse=True)
    energyEM1EM2_up.sort(reverse=True)

    energyEM1EM2 = energyEM1EM2_low[0] + energyEM1EM2_mid[0] + energyEM1EM2_up[0]

    # etaMax_PS = findMax(3,myTOB)[0]
    # phiMax_PS = findMax(2,myTOB)[1]
    # E_PS  = shapeTLV_fineGran(myTOB,2,etaMax_PS,phiMax_PS)
    # etaMax_EM3 = findMax(5,myTOB)[0]
    # phiMax_EM3 = findMax(5,myTOB)[1]
    # E_EM3  = shapeTLV_fineGran(myTOB,5,etaMax_EM3,phiMax_EM3)
    # etaMax_HAD = findMax(6,myTOB)[0]
    # phiMax_HAD = findMax(6,myTOB)[1]
    # E_HAD  = shapeTLV_fineGran(myTOB,6,etaMax_HAD,phiMax_HAD)

    # allE = energyEM1EM2+E_PS+E_EM3+E_HAD
    energyPS = allCells[0]
    energyEM3 = allCells[3]
    energyHAD = allCells[4]

    energyPS.sort(reverse=True)
    energyEM3.sort(reverse=True)
    energyHAD.sort(reverse=True)

    energyPSEM3HAD = (
        0.5 * energyPS[0]
        + 0.5 * energyPS[1]
        + 0.5 * energyPS[2]
        + 1.9 * energyEM3[0]
        + 1.9 * energyEM3[1]
        + energyHAD[0]
        + energyHAD[1]
        + energyHAD[2]
    )

    E_PS = 0.5 * energyPS[0] + 0.5 * energyPS[1] + 0.5 * energyPS[2]
    E_EM3 = 1.9 * energyEM3[0] + 1.9 * energyEM3[1]
    E_HAD = energyHAD[0] + energyHAD[1] + energyHAD[2]

    allE = energyEM1EM2 + energyPSEM3HAD

    return [allE, energyEM1EM2, E_PS, E_EM3, E_HAD]


def bigCluster6Cal(myTOB, SFHad):
    allCells = createCellLists(myTOB)

    # for phi=1 I'll add up the energy of the 6 centras supercells (15,16,17,18,19,20) then also that window shifted two SC to one side (13,14,15,16,17,18) and to the other (17,18,19,20,21,22)
    midPhiEM1 = [
        allCells[1][13]
        + allCells[1][14]
        + allCells[1][15]
        + allCells[1][16]
        + allCells[1][17]
        + allCells[1][18],
        allCells[1][15]
        + allCells[1][16]
        + allCells[1][17]
        + allCells[1][18]
        + allCells[1][19]
        + allCells[1][20],
        allCells[1][17]
        + allCells[1][18]
        + allCells[1][19]
        + allCells[1][20]
        + allCells[1][21]
        + allCells[1][22],
    ]
    midPhiEM2 = [
        allCells[2][13]
        + allCells[2][14]
        + allCells[2][15]
        + allCells[2][16]
        + allCells[2][17]
        + allCells[2][18],
        allCells[2][15]
        + allCells[2][16]
        + allCells[2][17]
        + allCells[2][18]
        + allCells[2][19]
        + allCells[2][20],
        allCells[2][17]
        + allCells[2][18]
        + allCells[2][19]
        + allCells[2][20]
        + allCells[2][21]
        + allCells[2][22],
    ]

    # for phi=0 and phi=1 i'll just loop over the two blocks sum energy (6,7 and 8,9 and 10,11 and 12,13)
    upPhiEM1 = [
        allCells[1][1] + allCells[1][2],
        allCells[1][3] + allCells[1][4],
        allCells[1][5] + allCells[1][6],
        allCells[1][7] + allCells[1][8],
        allCells[1][9] + allCells[1][10],
    ]
    upPhiEM2 = [
        allCells[2][1] + allCells[2][2],
        allCells[2][3] + allCells[2][4],
        allCells[2][5] + allCells[2][6],
        allCells[2][7] + allCells[2][8],
        allCells[2][9] + allCells[2][10],
    ]
    lowPhiEM1 = [
        allCells[1][1 + 24] + allCells[1][2 + 24],
        allCells[1][3 + 24] + allCells[1][4 + 24],
        allCells[1][5 + 24] + allCells[1][6 + 24],
        allCells[1][7 + 24] + allCells[1][8 + 24],
        allCells[1][9 + 24] + allCells[1][10 + 24],
    ]
    lowPhiEM2 = [
        allCells[2][1 + 24] + allCells[2][2 + 24],
        allCells[2][3 + 24] + allCells[2][4 + 24],
        allCells[2][5 + 24] + allCells[2][6 + 24],
        allCells[2][7 + 24] + allCells[2][8 + 24],
        allCells[2][9 + 24] + allCells[2][10 + 24],
    ]

    energyEM1EM2_mid = [
        midPhiEM1[0] + 1.0 * midPhiEM2[0],
        midPhiEM1[1] + 1.0 * midPhiEM2[1],
        midPhiEM1[2] + 1.0 * midPhiEM2[2],
    ]
    energyEM1EM2_low = [
        lowPhiEM1[0] + 1.0 * lowPhiEM2[0],
        lowPhiEM1[1] + 1.0 * lowPhiEM2[1],
        lowPhiEM1[2] + 1.0 * lowPhiEM2[2],
        lowPhiEM1[3] + 1.0 * lowPhiEM2[3],
        lowPhiEM1[4] + 1.0 * lowPhiEM2[4],
    ]
    energyEM1EM2_up = [
        upPhiEM1[0] + 1.0 * upPhiEM2[0],
        upPhiEM1[1] + 1.0 * upPhiEM2[1],
        upPhiEM1[2] + 1.0 * upPhiEM2[2],
        upPhiEM1[3] + 1.0 * upPhiEM2[3],
        upPhiEM1[4] + 1.0 * upPhiEM2[4],
    ]

    energyEM1EM2_mid.sort(reverse=True)
    energyEM1EM2_low.sort(reverse=True)
    energyEM1EM2_up.sort(reverse=True)

    energyEM1EM2 = energyEM1EM2_low[0] + energyEM1EM2_mid[0] + energyEM1EM2_up[0]

    energyPS = allCells[0]
    energyEM3 = allCells[3]
    energyHAD = allCells[4]

    energyPS.sort(reverse=True)
    energyEM3.sort(reverse=True)
    energyHAD.sort(reverse=True)

    energyPSEM3HAD = (
        0.5 * energyPS[0]
        + 0.5 * energyPS[1]
        + 0.5 * energyPS[2]
        + 1.9 * energyEM3[0]
        + 1.9 * energyEM3[1]
        + 1.9 * energyHAD[0]
        + 1 * energyHAD[1]
        + 1 * energyHAD[2]
    )

    E_PS = energyPS[0] + energyPS[1] + energyPS[2]
    E_EM3 = energyEM3[0] + energyEM3[1]
    E_HAD = energyHAD[0] + energyHAD[1] + energyHAD[2]

    allE = energyEM1EM2 + energyPSEM3HAD

    return [allE, energyEM1EM2, E_PS, E_EM3, E_HAD]


def bigClusIso(length, myTOB):

    if length == 4:
        energies = bigCluster(myTOB)
    if length == 6:
        energies = bigCluster6(myTOB)
    energyEM1EM2 = energies[1]
    E_HAD = energies[4]
    E_EM3 = energies[3]
    E_PS = energies[2]
    allE = energyEM1EM2 + E_PS + E_EM3 + E_HAD
    E_EM2_de = efexFullWindowIso(myTOB)[2]
    E_EM1_de = efexFullWindowIso(myTOB)[1]
    E_PS_de = efexFullWindowIso(myTOB)[0]
    E_EM3_de = efexFullWindowIso(myTOB)[3]
    E_HAD_de = efexFullWindowIso(myTOB)[4]
    try:
        Iso_EM1EM2EM3HAD = (energyEM1EM2 + E_HAD + E_EM3) / (
            E_EM1_de + E_EM2_de + E_EM3_de + E_HAD_de
        )
    except:
        Iso_EM1EM2EM3HAD = 0.0
    try:
        Iso_EM1EM2HAD = (energyEM1EM2 + E_HAD) / (E_EM1_de + E_EM2_de + E_HAD_de)
    except:
        Iso_EM1EM2HAD = 0.0
    try:
        Iso_EM1EM2 = energyEM1EM2 / (E_EM1_de + E_EM2_de)
    except:
        Iso_EM1EM2 = 0.0
    try:
        Iso_all = allE / (E_PS_de + E_EM1_de + E_EM2_de + E_EM3_de + E_HAD_de)
    except:
        Iso_all = 0.0
    return [
        Iso_all,
        Iso_EM1EM2,
        Iso_EM1EM2HAD,
        Iso_EM1EM2EM3HAD,
        allE,
        E_PS_de + E_EM1_de + E_EM2_de + E_EM3_de + E_HAD_de,
        energyEM1EM2,
        E_EM1_de + E_EM2_de,
        energyEM1EM2 + E_HAD,
        E_EM1_de + E_EM2_de + E_HAD_de,
        energyEM1EM2 + E_HAD + E_EM3,
        E_EM1_de + E_EM2_de + E_EM3_de + E_HAD_de,
    ]


def ratio_peak(myTOB):
    ratio = []
    core32_1 = 0
    core32_2 = 0
    core32 = []
    core12 = []
    shell = 0
    layer12_total = 0
    Oregon_ratio = 0
    list = createCellLists(myTOB)
    seed = findPeak(myTOB, 2)
    # make sure peak is in the central region
    if seed[5] < 16 or seed[5] > 19:
        ratio = [0, 0, 0]
    else:
        for i in [-1, 0, 1]:
            for j in [1, 2]:
                core32_1 += list[j][seed[5] + i] + list[j][seed[5] + 12 + i]
                core32_2 += list[j][seed[5] + i] + list[j][seed[5] - 12 + i]
        core32 = [core32_1, core32_2]
        core32.sort()

        # toget Oregon ratio
        for i in [1, 2]:
            for j in range(36):
                layer12_total += list[i][j]
        Oregon_ratio = core32[1] / layer12_total
        ratio += [Oregon_ratio]

    return ratio


def ApplyIsoXin(Snom, Sden):
    doCut = False
    if Sden < (0.016 * Snom + 4.083):
        doCut = True
    # print Snom, Sden, 0.016 * Sden, 0.016 * Sden + 4.083, doCut
    return doCut


def ApplyIsoXinTLV(Snom, Sden):  # run2
    doCut = False
    if Snom > (10 * (Sden - Snom) - 20.0):
        doCut = True
    # print Snom, Sden, 10 * Sden, 10 * Sden - 20.0, doCut
    return doCut


def ApplyIsoTLV(S, pt):
    doCut = False
    if pt > 50.0:
        doCut = True
    #    if S > 0.85: doCut = True
    # elif S > (0.0023 * pt + 0.76): doCut = True
    # print S, pt, 0.0023 * pt, 0.0023 * pt + 0.76
    # if S > (0.0028 * pt + 0.77): doCut = True
    elif S > (0.0015 * pt + 0.69):
        doCut = True
    return doCut


def ApplyIsoTLV1(S, pt):
    doCut = False
    if pt > 50:
        doCut = True
    #    if S > 0.9: doCut = True
    # if PT > 60000.0: doCut = True
    # elif S > (0.0028 * pt + 0.77): doCut = True
    elif S > (0.0009 * pt + 0.69):
        doCut = True
    # print S, pt, 0.0023 * pt, 0.0023 * pt + 0.76
    return doCut


def ApplyIsoTLV2(S, pt):
    doCut = False
    if pt > 50.0:
        doCut = True
    #    if S > 0.86: doCut = True
    # if PT > 60000.0: doCut = True
    # elif S > (0.0017 * pt + 0.8): doCut = True
    elif S > (0.0013 * pt + 0.69):
        doCut = True
    # print S, pt, 0.0023 * pt, 0.0023 * pt + 0.76
    return doCut


def ApplyIsoTLV3(S, pt):
    doCut = False
    if pt > 50.0:
        doCut = True
    #    if S > 0.85: doCut = True
    # if PT > 60000.0: doCut = True
    elif S > (0.0011 * pt + 0.69):
        doCut = True
    # print S, pt, 0.0011 * pt, 0.0011 * pt + 0.69
    return doCut


def efexFullWindowIso(myTOB):
    allCells = createCellLists(myTOB)
    EM0allCell = allCells[0]
    EM1allCell = allCells[1]
    EM2allCell = allCells[2]
    EM3allCell = allCells[3]
    HADallCell = allCells[4]
    SumEM0 = 0
    SumEM1 = 0
    SumEM2 = 0
    SumEM3 = 0
    SumHAD = 0
    for i in range(36):
        SumEM1 += EM1allCell[i]
        SumEM2 += EM2allCell[i]
    for i in range(9):
        SumEM0 += EM0allCell[i]
        SumEM3 += EM3allCell[i]
        SumHAD += HADallCell[i]
    return [SumEM0, SumEM1, SumEM2, SumEM3, SumHAD]


def efexRun2Big(myTOB):
    allCells = createCellLists(myTOB)
    EM0allCell = allCells[0]
    EM1allCell = allCells[1]
    EM2allCell = allCells[2]
    EM3allCell = allCells[3]
    HADallCell = allCells[4]
    EMsum = [0, 0, 0, 0]
    EMsum[0] = EM0allCell[1] + EM0allCell[4] + EM3allCell[1] + EM3allCell[4]
    EMsum[1] = EM0allCell[5] + EM0allCell[4] + EM3allCell[5] + EM3allCell[4]
    EMsum[2] = EM0allCell[7] + EM0allCell[4] + EM3allCell[7] + EM3allCell[4]
    EMsum[3] = EM0allCell[3] + EM0allCell[4] + EM3allCell[3] + EM3allCell[4]
    for i in range(4):
        EMsum[0] += EM1allCell[i + 16]
        EMsum[0] += EM1allCell[i + 4]
        EMsum[0] += EM2allCell[i + 16]
        EMsum[0] += EM2allCell[i + 4]
        EMsum[1] += EM1allCell[i + 16]
        EMsum[1] += EM1allCell[i + 20]
        EMsum[1] += EM2allCell[i + 16]
        EMsum[1] += EM2allCell[i + 20]
        EMsum[2] += EM1allCell[i + 16]
        EMsum[2] += EM1allCell[i + 28]
        EMsum[2] += EM2allCell[i + 16]
        EMsum[2] += EM2allCell[i + 28]
        EMsum[3] += EM1allCell[i + 16]
        EMsum[3] += EM1allCell[i + 12]
        EMsum[3] += EM2allCell[i + 16]
        EMsum[3] += EM2allCell[i + 12]
    HADsum = [0, 0, 0, 0]
    HADsum[0] = HADallCell[4] + HADallCell[0] + HADallCell[1] + HADallCell[3]
    HADsum[1] = HADallCell[4] + HADallCell[1] + HADallCell[2] + HADallCell[5]
    HADsum[2] = HADallCell[4] + HADallCell[3] + HADallCell[6] + HADallCell[7]
    HADsum[3] = HADallCell[4] + HADallCell[5] + HADallCell[7] + HADallCell[8]
    # look for a 1st maximum
    EMsumMax = -1.0
    EMsumMaxIndex = 0
    HADsumMax = -1.0
    HADsumMaxIndex = 0
    for i in range(4):
        if EMsum[i] > EMsumMax:
            EMsumMax = EMsum[i]
            EMsumMaxIndex = i
        if HADsum[i] > HADsumMax:
            HADsumMax = HADsum[i]
            HADsumMaxIndex = [i]
    mySum = EMsumMax + HADsumMax
    return mySum


def efexTDRclus(myTOB):
    allCells = createCellLists(myTOB)
    EM0allCell = allCells[0]
    EM1allCell = allCells[1]
    EM2allCell = allCells[2]
    EM3allCell = allCells[3]
    HADallCell = allCells[4]

    EM2seed = findseed(myTOB, 2)
    EM2seedindex = EM2seed[1]

    SumUpOrDown = -1  # 0:up,1:down

    EM2TTup = 0
    EM2TTdown = 0
    for i in range(4):
        EM2TTup += EM2allCell[i + 28]
        EM2TTdown += EM2allCell[i + 4]
    if EM2TTup >= EM2TTdown:
        SumUpOrDown = 0  # up
    if EM2TTdown > EM2TTup:
        SumUpOrDown = 1  # down

    # if EM2allCell[EM2seedindex-12] >= EM2allCell[EM2seedindex+12]:
    #    SumUpOrDown = 1#down
    #    if EM2allCell[EM2seedindex-12] == EM2allCell[EM2seedindex+12]:
    #        print 'WORNING : up down energy is same' , EM2allCell[EM2seedindex-12] , ',' , EM2allCell[EM2seedindex+12]
    # if EM2allCell[EM2seedindex-12] < EM2allCell[EM2seedindex+12]:
    #    SumUpOrDown = 0#up

    SumRightOrLeft = -1  # 0:left,1:right
    if EM2seedindex == 16 or EM2seedindex == 17:
        SumRightOrLeft = 0  # left
    if EM2seedindex == 18 or EM2seedindex == 19:
        SumRightOrLeft = 1  # right

    # EM2Sum = 0
    # EM2up = 0
    # EM2down = 0
    # for i in range(5):
    #    EM2down += EM2allCell[EM2seedindex-14+i]
    #    EM2down += EM2allCell[EM2seedindex-2+i]
    #    EM2up += EM2allCell[EM2seedindex-2+i]
    #    EM2up += EM2allCell[EM2seedindex+10+i]
    # if EM2up>=EM2down:
    #    EM2Sum=EM2up
    #    SumUpOrDown = 0#up
    # if EM2down>EM2up:
    #    EM2Sum=EM2down
    #    SumUpOrDown = 1#down

    EM0Sum = 0
    if SumUpOrDown == 1:  # down
        EM0Sum = EM0allCell[1] + EM0allCell[4]
    if SumUpOrDown == 0:  # up
        EM0Sum = EM0allCell[4] + EM0allCell[7]

    EM1Sum = 0
    if SumUpOrDown == 1:  # down
        for i in range(3):
            EM1Sum += EM1allCell[EM2seedindex - 13 + i]
            EM1Sum += EM1allCell[EM2seedindex - 1 + i]
    if SumUpOrDown == 0:  # up
        for i in range(3):
            EM1Sum += EM1allCell[EM2seedindex - 1 + i]
            EM1Sum += EM1allCell[EM2seedindex + 11 + i]

    EM2Sum = 0
    if SumUpOrDown == 1:  # down
        for i in range(5):
            EM2Sum += EM2allCell[EM2seedindex - 14 + i]
            EM2Sum += EM2allCell[EM2seedindex - 2 + i]
    if SumUpOrDown == 0:  # up
        for i in range(5):
            EM2Sum += EM2allCell[EM2seedindex - 2 + i]
            EM2Sum += EM2allCell[EM2seedindex + 10 + i]

    EM3Sum = 0
    HADSum = 0
    if SumUpOrDown == 1 and SumRightOrLeft == 0:  # down,left
        EM3Sum = EM3allCell[0] + EM3allCell[1] + EM3allCell[3] + EM3allCell[4]
        HADSum = HADallCell[0] + HADallCell[1] + HADallCell[3] + HADallCell[4]
    if SumUpOrDown == 1 and SumRightOrLeft == 1:  # down,right
        EM3Sum = EM3allCell[1] + EM3allCell[2] + EM3allCell[4] + EM3allCell[5]
        HADSum = HADallCell[1] + HADallCell[2] + HADallCell[4] + HADallCell[5]
    if SumUpOrDown == 0 and SumRightOrLeft == 0:  # up,left
        EM3Sum = EM3allCell[3] + EM3allCell[4] + EM3allCell[6] + EM3allCell[7]
        HADSum = HADallCell[3] + HADallCell[4] + HADallCell[6] + HADallCell[7]
    if SumUpOrDown == 0 and SumRightOrLeft == 1:  # up,right
        EM3Sum = EM3allCell[4] + EM3allCell[5] + EM3allCell[7] + EM3allCell[8]
        HADSum = HADallCell[4] + HADallCell[5] + HADallCell[7] + HADallCell[8]

    mySum = EM0Sum + EM1Sum + EM2Sum + EM3Sum + HADSum
    return mySum


def efexOregonClus(myTOB):
    allCells = createCellLists(myTOB)
    EM0allCell = allCells[0]
    EM1allCell = allCells[1]
    EM2allCell = allCells[2]
    EM3allCell = allCells[3]
    HADallCell = allCells[4]

    EM2seed = findseed(myTOB, 2)
    EM2seedindex = EM2seed[1]

    SumUpOrDown = -1  # 0:up,1:down
    EM2TTup = 0
    EM2TTdown = 0
    for i in range(4):
        EM2TTup += EM2allCell[i + 28]
        EM2TTdown += EM2allCell[i + 4]
    if EM2TTup >= EM2TTdown:
        SumUpOrDown = 0  # up
    if EM2TTdown > EM2TTup:
        SumUpOrDown = 1  # down

    # if EM2allCell[EM2seedindex-12] >= EM2allCell[EM2seedindex+12]:
    #    SumUpOrDown = 1#down
    #    if EM2allCell[EM2seedindex-12] == EM2allCell[EM2seedindex+12]:
    #        print 'WORNING : up down energy is same' , EM2allCell[EM2seedindex-12] , ',' , EM2allCell[EM2seedindex+12]
    # if EM2allCell[EM2seedindex-12] < EM2allCell[EM2seedindex+12]:
    #    SumUpOrDown = 0#up

    SumRightOrLeft = -1  # 0:left,1:right
    if EM2seedindex == 16 or EM2seedindex == 17:
        SumRightOrLeft = 0  # left
    if EM2seedindex == 18 or EM2seedindex == 19:
        SumRightOrLeft = 1  # right

    # EM2Sum = 0
    # EM2up = 0
    # EM2down = 0
    # for i in range(5):
    #    EM2down += EM2allCell[EM2seedindex-14+i]
    #    EM2down += EM2allCell[EM2seedindex-2+i]
    #    EM2up += EM2allCell[EM2seedindex-2+i]
    #    EM2up += EM2allCell[EM2seedindex+10+i]
    # if EM2up>=EM2down:
    #    EM2Sum=EM2up
    #    SumUpOrDown = 0#up
    # if EM2down>EM2up:
    #    EM2Sum=EM2down
    #    SumUpOrDown = 1#down

    EM0Sum = 0
    if SumUpOrDown == 1:  # down
        EM0Sum = EM0allCell[1] + EM0allCell[4]
    if SumUpOrDown == 0:  # up
        EM0Sum = EM0allCell[4] + EM0allCell[7]

    EM1Sum = 0
    if SumUpOrDown == 1:  # down
        for i in range(5):
            EM1Sum += EM1allCell[EM2seedindex - 14 + i]
            EM1Sum += EM1allCell[EM2seedindex - 2 + i]
    if SumUpOrDown == 0:  # up
        for i in range(5):
            EM1Sum += EM1allCell[EM2seedindex - 2 + i]
            EM1Sum += EM1allCell[EM2seedindex + 10 + i]

    EM2Sum = 0
    if SumUpOrDown == 1:  # down
        for i in range(5):
            EM2Sum += EM2allCell[EM2seedindex - 14 + i]
            EM2Sum += EM2allCell[EM2seedindex - 2 + i]
    if SumUpOrDown == 0:  # up
        for i in range(5):
            EM2Sum += EM2allCell[EM2seedindex - 2 + i]
            EM2Sum += EM2allCell[EM2seedindex + 10 + i]

    EM3Sum = 0
    HADSum = 0
    if SumUpOrDown == 1 and SumRightOrLeft == 0:  # down,left
        EM3Sum = EM3allCell[0] + EM3allCell[1] + EM3allCell[3] + EM3allCell[4]
        HADSum = HADallCell[0] + HADallCell[1] + HADallCell[3] + HADallCell[4]
    if SumUpOrDown == 1 and SumRightOrLeft == 1:  # down,right
        EM3Sum = EM3allCell[1] + EM3allCell[2] + EM3allCell[4] + EM3allCell[5]
        HADSum = HADallCell[1] + HADallCell[2] + HADallCell[4] + HADallCell[5]
    if SumUpOrDown == 0 and SumRightOrLeft == 0:  # up,left
        EM3Sum = EM3allCell[3] + EM3allCell[4] + EM3allCell[6] + EM3allCell[7]
        HADSum = HADallCell[3] + HADallCell[4] + HADallCell[6] + HADallCell[7]
    if SumUpOrDown == 0 and SumRightOrLeft == 1:  # up,right
        EM3Sum = EM3allCell[4] + EM3allCell[5] + EM3allCell[7] + EM3allCell[8]
        HADSum = HADallCell[4] + HADallCell[5] + HADallCell[7] + HADallCell[8]

    mySum = EM0Sum + EM1Sum + EM2Sum + EM3Sum + HADSum
    #    print("OREGON EM1 energy {}".format(EM1Sum))
    return mySum


def JefexOregonClus(myTOB):
    allCells = createCellLists(myTOB)
    EM0allCell = allCells[0]
    EM1allCell = allCells[1]
    EM2allCell = allCells[2]
    EM3allCell = allCells[3]
    HADallCell = allCells[4]

    EM2seed = findseed(myTOB, 2)
    EM2seedindex = EM2seed[1]

    # choose up or bottom row? check highest energy between above or below seed in EM2 plus the same cell in EM1

    doUp = False
    #    if((EM2allCell[EM2seedindex-12]+EM1allCell[EM2seedindex-12]) > (EM2allCell[EM2seedindex+12]+EM1allCell[EM2seedindex+12])):
    if (EM2allCell[EM2seedindex - 12] + EM1allCell[EM2seedindex - 12]) >= (
        EM2allCell[EM2seedindex + 12] + EM1allCell[EM2seedindex + 12]
    ):
        doUp = True

    EM0 = 0
    EM12 = 0
    EM12_midRow = 0
    EM12_extraPhi = 0
    EM3 = 0
    HAD = 0

    for i in range(5):
        EM12_midRow += (
            EM2allCell[EM2seedindex - 2 + i] + EM1allCell[EM2seedindex - 2 + i]
        )
        if doUp:
            EM12_extraPhi += (
                EM2allCell[EM2seedindex - 2 + i - 12]
                + EM1allCell[EM2seedindex - 2 + i - 12]
            )
        else:
            EM12_extraPhi += (
                EM1allCell[EM2seedindex - 2 + i + 12]
                + EM2allCell[EM2seedindex - 2 + i + 12]
            )
    EM12 = EM12_midRow + EM12_extraPhi

    if doUp:
        EM0 = EM0allCell[1] + EM0allCell[4]
        if EM2seedindex == 16:
            EM0 += EM0allCell[0] + EM0allCell[3]
        elif EM2seedindex == 19:
            EM0 += EM0allCell[2] + EM0allCell[5]
    else:
        EM0 = EM0allCell[4] + EM0allCell[1 + 3 + 3]
        if EM2seedindex == 16:
            EM0 += EM0allCell[3] + EM0allCell[6]
        elif EM2seedindex == 19:
            EM0 += EM0allCell[5] + EM0allCell[8]

    for i in range(3):
        EM3 += EM3allCell[i + 3]
        HAD += HADallCell[i + 3]
        if doUp:
            EM3 += EM3allCell[i]
            HAD += HADallCell[i]
        else:
            EM3 += EM3allCell[i + 6]
            HAD += HADallCell[i + 6]

    mySum = EM0 + EM12 + EM3 + HAD
    return mySum


def JefexOregonClusCal(myTOB, SFHad):
    allCells = createCellLists(myTOB)
    EM0allCell = allCells[0]
    EM1allCell = allCells[1]
    EM2allCell = allCells[2]
    EM3allCell = allCells[3]
    HADallCell = allCells[4]

    EM2seed = findseed(myTOB, 2)
    EM2seedindex = EM2seed[1]

    # choose up or bottom row? check highest energy between above or below seed in EM2 plus the same cell in EM1

    doUp = False
    if (EM2allCell[EM2seedindex - 12] + EM1allCell[EM2seedindex - 12]) > (
        EM2allCell[EM2seedindex + 12] + EM1allCell[EM2seedindex + 12]
    ):
        doUp = True

    EM0 = 0
    EM12 = 0
    EM12_midRow = 0
    EM12_extraPhi = 0
    EM3 = 0
    HAD = 0

    for i in range(5):
        EM12_midRow += (
            EM2allCell[EM2seedindex - 2 + i] + EM1allCell[EM2seedindex - 2 + i]
        )
        if doUp:
            EM12_extraPhi += (
                EM2allCell[EM2seedindex - 2 + i - 12]
                + EM1allCell[EM2seedindex - 2 + i - 12]
            )
        else:
            EM12_extraPhi += (
                EM1allCell[EM2seedindex - 2 + i + 12]
                + EM2allCell[EM2seedindex - 2 + i + 12]
            )
    EM12 = EM12_midRow + EM12_extraPhi

    if doUp:
        EM0 = EM0allCell[1] + EM0allCell[4]
        if EM2seedindex == 16:
            EM0 += EM0allCell[0] + EM0allCell[3]
        elif EM2seedindex == 19:
            EM0 += EM0allCell[2] + EM0allCell[5]
    else:
        EM0 = EM0allCell[4] + EM0allCell[1 + 3 + 3]
        if EM2seedindex == 16:
            EM0 += EM0allCell[3] + EM0allCell[6]
        elif EM2seedindex == 19:
            EM0 += EM0allCell[5] + EM0allCell[8]

    for i in range(3):
        EM3 += EM3allCell[i + 3]
        HAD += HADallCell[i + 3]
        if doUp:
            EM3 += EM3allCell[i]
            HAD += HADallCell[i]
        else:
            EM3 += EM3allCell[i + 6]
            HAD += HADallCell[i + 6]

    mySum = 0.5 * EM0 + EM12 + 1.9 * EM3 + HAD
    return mySum


def efexOregonClusCal(myTOB, SFHad):
    allCells = createCellLists(myTOB)
    EM0allCell = allCells[0]
    EM1allCell = allCells[1]
    EM2allCell = allCells[2]
    EM3allCell = allCells[3]
    HADallCell = allCells[4]

    EM2seed = findseed(myTOB, 2)
    EM2seedindex = EM2seed[1]

    SumUpOrDown = -1  # 0:up,1:down
    EM2TTup = 0
    EM2TTdown = 0
    for i in range(4):
        EM2TTup += EM2allCell[i + 28]
        EM2TTdown += EM2allCell[i + 4]
    if EM2TTup >= EM2TTdown:
        SumUpOrDown = 0  # up
    if EM2TTdown > EM2TTup:
        SumUpOrDown = 1  # down

    SumRightOrLeft = -1  # 0:left,1:right
    if EM2seedindex == 16 or EM2seedindex == 17:
        SumRightOrLeft = 0  # left
    if EM2seedindex == 18 or EM2seedindex == 19:
        SumRightOrLeft = 1  # right

    EM0Sum = 0
    if SumUpOrDown == 1:  # down
        EM0Sum = EM0allCell[1] + EM0allCell[4]
    if SumUpOrDown == 0:  # up
        EM0Sum = EM0allCell[4] + EM0allCell[7]

    EM1Sum = 0
    if SumUpOrDown == 1:  # down
        for i in range(5):
            EM1Sum += EM1allCell[EM2seedindex - 14 + i]
            EM1Sum += EM1allCell[EM2seedindex - 2 + i]
    if SumUpOrDown == 0:  # up
        for i in range(5):
            EM1Sum += EM1allCell[EM2seedindex - 2 + i]
            EM1Sum += EM1allCell[EM2seedindex + 10 + i]

    EM2Sum = 0
    if SumUpOrDown == 1:  # down
        for i in range(5):
            EM2Sum += EM2allCell[EM2seedindex - 14 + i]
            EM2Sum += EM2allCell[EM2seedindex - 2 + i]
    if SumUpOrDown == 0:  # up
        for i in range(5):
            EM2Sum += EM2allCell[EM2seedindex - 2 + i]
            EM2Sum += EM2allCell[EM2seedindex + 10 + i]

    EM3Sum = 0
    HADSum = 0
    if SumUpOrDown == 1 and SumRightOrLeft == 0:  # down,left
        EM3Sum = EM3allCell[0] + EM3allCell[1] + EM3allCell[3] + EM3allCell[4]
        HADSum = HADallCell[0] + HADallCell[1] + HADallCell[3] + HADallCell[4]
    if SumUpOrDown == 1 and SumRightOrLeft == 1:  # down,right
        EM3Sum = EM3allCell[1] + EM3allCell[2] + EM3allCell[4] + EM3allCell[5]
        HADSum = HADallCell[1] + HADallCell[2] + HADallCell[4] + HADallCell[5]
    if SumUpOrDown == 0 and SumRightOrLeft == 0:  # up,left
        EM3Sum = EM3allCell[3] + EM3allCell[4] + EM3allCell[6] + EM3allCell[7]
        HADSum = HADallCell[3] + HADallCell[4] + HADallCell[6] + HADallCell[7]
    if SumUpOrDown == 0 and SumRightOrLeft == 1:  # up,right
        EM3Sum = EM3allCell[4] + EM3allCell[5] + EM3allCell[7] + EM3allCell[8]
        HADSum = HADallCell[4] + HADallCell[5] + HADallCell[7] + HADallCell[8]

    mySum = EM0Sum + EM1Sum + 1.3 * EM2Sum + 1.3 * EM3Sum + 1.0 * HADSum
    #    print("OREGON EM1 energy {}".format(EM1Sum))
    return mySum


def efexFullWindow(myTOB):
    allCells = createCellLists(myTOB)
    EM0allCell = allCells[0]
    EM1allCell = allCells[1]
    EM2allCell = allCells[2]
    EM3allCell = allCells[3]
    HADallCell = allCells[4]

    Sum = 0
    for i in range(36):
        Sum += EM1allCell[i]
        Sum += EM2allCell[i]
    for i in range(9):
        Sum += EM0allCell[i]
        Sum += EM3allCell[i]
        Sum += HADallCell[i]

    # print 'efexFullWindow' , Sum

    return Sum


def layerFrac(myTOB, layer):
    # get ET map
    allCells = createCellLists(myTOB)
    if myTOB.largeTauClus() < 0.001:
        return 0
    else:
        return myLayerSum(allCells[layer]) / myTOB.largeTauClus()


def layerE(myTOB, layer):
    # get ET map
    allCells = createCellLists(myTOB)
    #    if myTOB.largeTauClus() < 0.001:
    #        return 0
    #    else:
    return myLayerSum(allCells[layer])


def findPeak(myTOB, layer):
    allCells = createCellLists(myTOB)
    layerEta = [0.1, 0.025, 0.025, 0.1, 0.1]
    myList = allCells[layer]

    # look for a 1st maximum
    max1 = -1.0
    max1index = -1
    for i in range(len(myList)):
        if myList[i] > max1:
            max1 = myList[i]
            max1index = i

    max2 = -1
    max2index = -1

    # find the 2nd maximum, ignore the first (and nearby cells)
    for j in range(len(myList)):
        if myList[j] > max2 and nearbyCells(layer, max1index, j)[0] > 1.01:
            max2 = myList[j]
            max2index = j

    ratio = 0
    if max1 > 0.001:
        ratio = max2 / max1

    distance = 0
    if max1index != -1 and max2index != -1:
        distance = nearbyCells(layer, max1index, max2index)[1]
        # distance = layerEta[layer] * abs(max2index-max1index)

    # return [max1, max2, ratio, distance, myLayerSum(myList)]
    return [
        max1,
        max2,
        ratio,
        distance,
        myLayerSum(myList),
        max1index,
        max2index,
    ]  # output changed


def findseed(myTOB, layer):
    allCells = createCellLists(myTOB)
    myList = allCells[layer]

    # look for a 1st maximum in centre TT
    max1 = -1.0
    max1index = -1
    for i in range(4):
        if myList[i + 16] > max1:
            max1 = myList[i + 16]
            max1index = i + 16
    #    print("For layer {} using findseed print energy {} and index {}".format(layer,max1,max1index))
    return [max1, max1index]


def nearbyCells(layer, index1, index2):
    # layer granularity
    layerCells = [3, 12, 12, 3, 3]

    # we have a maximum of 3 rows, so let's just assume.
    phiIndex = [0, 0]
    etaIndex = [index1, index2]

    for i in range(2):
        if etaIndex[i] > (layerCells[layer] - 1):
            phiIndex[i] += 1
            etaIndex[i] -= layerCells[layer]
        if etaIndex[i] > (layerCells[layer] - 1):
            phiIndex[i] += 1
            etaIndex[i] -= layerCells[layer]

    # we should now have corrected coordinates
    # calculate a distance, in integers
    distance1 = abs(phiIndex[0] - phiIndex[1]) + abs(etaIndex[0] - etaIndex[1])
    distance2 = 0.1 * abs(phiIndex[0] - phiIndex[1]) + (0.3 / layerCells[layer]) * abs(
        etaIndex[0] - etaIndex[1]
    )
    return [distance1, distance2]


def findMax(layer, myTOB):
    maxE = -99
    maxEindex = [9999, 9999]  # eta,phi
    #    print(maxE)
    #    print(maxEindex)
    if layer == 2 or layer == 5 or layer == 6:
        etaRange = 3
    else:
        etaRange = 12
    #    print("layer %i"%layer)
    #    print("etarange %i"%etaRange)
    for i in range(3):  # phi loop
        for j in range(etaRange):  # eta loop
            if layer == 2 or layer == 5 or layer == 6:
                eCell = myTOB.getEnergy(layer, j + 1, i + 1)
            else:
                eCell = myTOB.getEnergy(layer, j + 4, i + 1)
            if eCell > maxE:
                maxE = eCell
                maxEindex = [j, i]
            else:
                continue
    #    print(maxE)
    #    print(maxEindex)
    # print("for layer {} using findMax energy {} and index {}".format(layer,maxE,maxEindex))
    E_Index = [maxE, maxEindex]
    return maxEindex


def findMax2(layer, myTOB):
    maxE = -99
    maxEindex = [9999, 9999]  # eta,phi
    #    print(maxE)
    #    print(maxEindex)
    if layer == 2 or layer == 5 or layer == 6:
        etaRange = 3
    else:
        etaRange = 12
    #    print("layer %i"%layer)
    #    print("etarange %i"%etaRange)
    for i in range(3):  # phi loop
        for j in range(etaRange):  # eta loop
            if layer == 2 or layer == 5 or layer == 6:
                eCell = myTOB.getEnergy(layer, j + 1, i + 1)
            else:
                eCell = myTOB.getEnergy(layer, j + 4, i + 1)
            if eCell > maxE:
                maxE = eCell
                maxEindex = [j, i]
            else:
                continue
    E_Index = [maxE, maxEindex]
    return E_Index


##############################################
def findCentralTower(myTOB):
    # print 'layer phi eta CT central'
    CT = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
    central = 0
    # build the central towers
    for layer in [2, 3, 4, 5, 6]:
        for i in range(3):  # phi loop
            for j in range(3):  # eta loop
                if layer == 2 or layer == 5 or layer == 6:
                    CT[i][j] += myTOB.getEnergy(layer, j + 1, i + 1)
                # else: CT[i][j] += (myTOB.getEnergy(layer, j+4, i+1)+myTOB.getEnergy(layer, j+5, i+1)+myTOB.getEnergy(layer, j+6, i+1)+myTOB.getEnergy(layer, j+7, i+1))
                else:
                    CT[i][j] += (
                        myTOB.getEnergy(layer, 4 * j + 4, i + 1)
                        + myTOB.getEnergy(layer, 4 * j + 5, i + 1)
                        + myTOB.getEnergy(layer, 4 * j + 6, i + 1)
                        + myTOB.getEnergy(layer, 4 * j + 7, i + 1)
                    )
                # print myTOB, layer, i, j, CT[1][1], central
    # check the central ([1][1]) is the highest
    # if CT[1][1] >= CT[0][0] and CT[1][1] > CT[0][1] and CT[1][1] > CT[0][2] and CT[1][1] > CT[1][2] and CT[1][1] >= CT[1][0] and CT[1][1] >= CT[2][0] and CT[1][1] >= CT[2][1] and CT[1][1] >= CT[2][2]: central = 1
    if (
        CT[1][1] >= CT[0][0]
        and CT[1][1] >= CT[0][1]
        and CT[1][1] > CT[0][2]
        and CT[1][1] > CT[1][2]
        and CT[1][1] >= CT[1][0]
        and CT[1][1] >= CT[2][0]
        and CT[1][1] > CT[2][1]
        and CT[1][1] > CT[2][2]
    ):
        central = 1
    return central


def findMaxEM2(myTOB):
    maxE = -99
    maxEindex = [9999, 9999]  # eta,phi
    for i in range(1, 4):  # phi loop
        for j in range(4, 16):  # eta loop
            eCell = myTOB.getEnergy(4, j, i)
            if eCell > maxE:
                maxE = eCell
                maxEindex = [j, i]
            else:
                continue
    return maxEindex


def findMaxEM1(myTOB):
    maxE = -99
    maxEindex = [9999, 9999]  # eta,phi
    for i in range(1, 4):  # phi loop
        for j in range(4, 16):  # eta loop
            eCell = myTOB.getEnergy(3, j, i)
            if eCell > maxE:
                maxE = eCell
                maxEindex = [j, i]
            else:
                continue
    return maxEindex


def mapGranularity(eta):
    etaMap = -1
    if eta >= 4 and eta <= 7:
        etaMap = 1
    if eta >= 8 and eta <= 11:
        etaMap = 2
    if eta >= 12 and eta <= 15:
        etaMap = 3
    return etaMap


def findAllMax(myTOB):
    # we start looking for the maximum in layer EM2
    eta_Max_EM2 = findMaxEM2(myTOB)[0]
    phi_Max_EM2 = findMaxEM2(myTOB)[1]

    # for the max in EM1, we start with the same max as in EM2 and look around:
    eta_Max_EM1 = eta_Max_EM2
    phi_Max_EM1 = phi_Max_EM2
    maxE = -99
    maxEindex = [9999, 9999]  # eta,phi
    etaR = []
    phiR = []
    if eta_Max_EM2 == 4:
        etaR.append(eta_Max_EM2)
        etaR.append(eta_Max_EM2 + 1)
    elif eta_Max_EM2 == 15:
        etaR.append(eta_Max_EM2 - 1)
        etaR.append(eta_Max_EM2)
    else:
        etaR.append(eta_Max_EM2 - 1)
        etaR.append(eta_Max_EM2)
        etaR.append(eta_Max_EM2 + 1)
    if phi_Max_EM2 == 1:  # phi_max_EM2 should always be in the middle though!
        phiR = [phi_Max_EM2, phi_Max_EM2 + 1]
    elif phi_Max_EM2 == 2:  # phi_max_EM2 should always be in the middle though!
        phiR = [phi_Max_EM2 - 1, phi_Max_EM2, phi_Max_EM2 + 1]
    elif phi_Max_EM2 == 3:  # phi_max_EM2 should always be in the middle though!
        phiR = [phi_Max_EM2 - 1, phi_Max_EM2]
    for i in etaR:  # eta loop
        for j in phiR:  # phi loop
            eCell = myTOB.getEnergy(3, i, j)
            if eCell > maxE:
                maxE = eCell
                maxEindex = [i, j]
            else:
                continue
    eta_Max_EM1 = maxEindex[0]
    phi_Max_EM1 = maxEindex[1]

    # we do the same goin back to the PS, now starting from EM1 and mapping to the correct granularity
    eta_Max_PS = mapGranularity(eta_Max_EM1)
    phi_Max_PS = phi_Max_EM1
    maxE = -99
    maxEindex = [9999, 9999]  # eta,phi
    etaR = []
    phiR = []
    if eta_Max_PS == 1:
        etaR = [eta_Max_PS, eta_Max_PS + 1]
    elif eta_Max_PS == 2:  # phi_max_EM2 should always be in the middle though!
        etaR = [eta_Max_PS - 1, eta_Max_PS, eta_Max_PS + 1]
    elif eta_Max_PS == 3:
        etaR = [eta_Max_PS - 1, eta_Max_PS]
    if phi_Max_PS == 1:  # phi_max_EM2 should always be in the middle though!
        phiR = [phi_Max_PS, phi_Max_PS + 1]
    elif phi_Max_PS == 2:  # phi_max_EM2 should always be in the middle though!
        phiR = [phi_Max_PS - 1, phi_Max_PS, phi_Max_PS + 1]
    elif phi_Max_PS == 3:  # phi_max_EM2 should always be in the middle though!
        phiR = [phi_Max_PS - 1, phi_Max_PS]
    for i in etaR:  # eta loop
        for j in phiR:  # phi loop
            eCell = myTOB.getEnergy(2, i, j)
            if eCell > maxE:
                maxE = eCell
                maxEindex = [i, j]
            else:
                continue
    eta_Max_PS = maxEindex[0]
    phi_Max_PS = maxEindex[1]

    # we do the same goin forward to the EM3, now starting from EM2 and mapping to the correct granularity
    eta_Max_EM3 = mapGranularity(eta_Max_EM2)
    phi_Max_EM3 = phi_Max_EM2
    maxE = -99
    maxEindex = [9999, 9999]  # eta,phi
    etaR = []
    phiR = []
    if eta_Max_EM3 == 1:
        etaR = [eta_Max_EM3, eta_Max_EM3 + 1]
    elif eta_Max_EM3 == 2:  # phi_max_EM2 should always be in the middle though!
        etaR = [eta_Max_EM3 - 1, eta_Max_EM3, eta_Max_EM3 + 1]
    elif eta_Max_EM3 == 3:
        etaR = [eta_Max_EM3 - 1, eta_Max_EM3]
    if phi_Max_EM3 == 1:  # phi_max_EM2 should always be in the middle though!
        phiR = [phi_Max_EM3, phi_Max_EM3 + 1]
    elif phi_Max_EM3 == 2:  # phi_max_EM2 should always be in the middle though!
        phiR = [phi_Max_EM3 - 1, phi_Max_EM3, phi_Max_EM3 + 1]
    elif phi_Max_EM3 == 3:  # phi_max_EM2 should always be in the middle though!
        phiR = [phi_Max_EM3 - 1, phi_Max_EM3]
    for i in etaR:  # eta loop
        for j in phiR:  # phi loop
            eCell = myTOB.getEnergy(5, i, j)
            if eCell > maxE:
                maxE = eCell
                maxEindex = [i, j]
            else:
                continue
    eta_Max_EM3 = maxEindex[0]
    phi_Max_EM3 = maxEindex[1]

    # finally we do the same goin forward to the HAD, now starting from EM3
    eta_Max_HAD = eta_Max_EM3
    phi_Max_HAD = phi_Max_EM3
    maxE = -99
    maxEindex = [9999, 9999]  # eta,phi
    etaR = []
    phiR = []
    if eta_Max_HAD == 1:
        etaR = [eta_Max_HAD, eta_Max_HAD + 1]
    elif eta_Max_HAD == 2:  # phi_max_EM2 should always be in the middle though!
        etaR = [eta_Max_HAD - 1, eta_Max_HAD, eta_Max_HAD + 1]
    elif eta_Max_HAD == 3:
        etaR = [eta_Max_HAD - 1, eta_Max_HAD]
    if phi_Max_HAD == 1:  # phi_max_EM2 should always be in the middle though!
        phiR = [phi_Max_HAD, phi_Max_HAD + 1]
    elif phi_Max_HAD == 2:  # phi_max_EM2 should always be in the middle though!
        phiR = [phi_Max_HAD - 1, phi_Max_HAD, phi_Max_HAD + 1]
    elif phi_Max_HAD == 3:  # phi_max_EM2 should always be in the middle though!
        phiR = [phi_Max_HAD - 1, phi_Max_HAD]
    for i in etaR:  # eta loop
        for j in phiR:  # phi loop
            eCell = myTOB.getEnergy(6, i, j)
            if eCell > maxE:
                maxE = eCell
                maxEindex = [i, j]
            else:
                continue
    eta_Max_HAD = maxEindex[0]
    phi_Max_HAD = maxEindex[1]

    maxAll = [
        [eta_Max_PS, phi_Max_PS],
        [eta_Max_EM1, phi_Max_EM1],
        [eta_Max_EM2, phi_Max_EM2],
        [eta_Max_EM3, phi_Max_EM3],
        [eta_Max_HAD, phi_Max_HAD],
    ]

    return maxAll


def longLineEM2(length, hotEta):
    etaGran = [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
    # etaGran = [1 ,2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    line = []
    size = length  # how long do you want to make it?
    for x in range(len(etaGran) - size + 1):
        l = etaGran[x : x + size]
        if hotEta in l:
            line.append(l)
    return line


def EM2shape_old(myTOB, etaMax, length):
    #    etaMax_EM2 = findMaxEM2(myTOB)[0]
    etaMax_EM2 = etaMax
    longLine = longLineEM2(length, etaMax_EM2)

    final_shape = []
    for line in longLine:
        #        print("each shape in EM2 and EM1 {}".format(line))
        tempShape = []
        #        block_1 = []
        #        block_2 = []
        b1b2 = []
        #        for i in range(line[0],line[0]+5):
        #            block_1.append(i)
        #            block_2.append(i)
        for i in range(line[0], line[0] + 5):
            for j in range(line[0], line[0] + 5):
                bb = [i, j]
                b1b2.append(bb)
        tempShape = [line, b1b2]
        # print("the shape in EM2 and EM1 {}".format(tempShape))
        final_shape.append(tempShape)
        # final_shape is a list whose first element is the long line, the second is a list where the first element is the eta index of the above block and the second one the one below
        # ahora si phi es 1(o 3) dejo solo el bloque de arriba(o abajo), la idea es implementar los dos arriba o los dos abajo para que se conserve la forma de 7 bloques
    return final_shape


def EM2shape(myTOB, etaMax, phiMax, length):
    #    etaMax_EM2 = findMaxEM2(myTOB)[0]

    longLine = longLineEM2(length, etaMax)
    b1b2 = []
    if phiMax == 2:
        b1b2.append((etaMax, phiMax + 1))
        b1b2.append((etaMax, phiMax - 1))
    if phiMax == 3:
        b1b2.append((etaMax, phiMax - 1))
    if phiMax == 1:
        b1b2.append((etaMax, phiMax + 1))
    final_shape = []
    for line in longLine:
        tempShape = [line, b1b2]
        final_shape.append(tempShape)
    return longLine


# return final_shape


def is_valid(el, shape, size):
    #      size = 4
    if len(shape) >= size:
        return False

    for (r, c) in shape:
        if el[0] == r and abs(el[1] - c) == 1 and el not in shape:
            return True
        if el[1] == c and abs(el[0] - r) == 1 and el not in shape:
            return True

    return False


def doShapes(layer, hot):

    if layer == 2:
        size = 2
    if layer == 5:
        size = 2
    if layer == 6:
        size = 3

    shapes = []

    # El 50 es arbitrario, no se como poner otra condicion (seria hasta que no haya mas formas)
    availableShapes = 50
    for k in range(availableShapes):

        shape = [
            hot,
        ]
        if hot[0] == 0 and hot[1] == 0:
            phiRange = [0, 1, 2]  # para los corners
            etaRange = [0, 1, 2]  # para los corners
        if hot[0] == 0 and hot[1] == 2:
            phiRange = [2, 1, 0]  # para los corners
            etaRange = [0, 1, 2]  # para los corners
        if hot[0] == 2 and hot[1] == 0:
            etaRange = [2, 1, 0]  # para los corners
            phiRange = [0, 1, 2]
        if hot[0] == 2 and hot[1] == 2:
            etaRange = [2, 1, 0]  # para los corners
            phiRange = [2, 1, 0]

        if hot[0] == 0 and hot[1] == 1:
            etaRange = [0, 1, 2]
            phiRange = [1, 2, 0]
        if hot[0] == 2 and hot[1] == 1:
            etaRange = [2, 1, 0]
            phiRange = [1, 2, 0]
        if hot[0] == 1 and hot[1] == 0:
            phiRange = [0, 1, 2]
            etaRange = [1, 2, 0]
        if hot[0] == 1 and hot[1] == 2:
            phiRange = [2, 1, 0]
            etaRange = [1, 0, 2]
        if hot[0] == 1 and hot[1] == 1:
            phiRange = [1, 0, 2]
            etaRange = [1, 0, 2]

        # for i in [0,1,2]:
        #    for j in [1,2,0]:
        for i in etaRange:
            for j in phiRange:
                if is_valid((i, j), shape, size):
                    tmp = shape + [(i, j)]
                    # tmp = list(shape) + [(i,j),]
                    if tmp in shapes:
                        continue

                    shape.append((i, j))

        shapes.append(shape)
    shapes = [shape for shape in shapes if len(shape) == size]
    return shapes


def shapeTLV_fineGran(myTOB, layer, etaMax, phimax):
    if layer == 4 or layer == 3:
        # print("el phi max es {}".format(phiMax_EM2))
        phiMax = phimax
        if layer == 4:
            length = 4
        if layer == 3:
            length = 3
        # length=4
        shapes = EM2shape(myTOB, etaMax, phiMax, length)
        # shapes = EM2shape(myTOB,etaMax,length)
        # shapes = EM2shape(myTOB,etaMax_EM2)
        storeMaxE = -9999
        storeIndex = []
        for lines in shapes:
            # storeMaxE = -9999 #i'm not sure if this has to be here or before
            # storeIndex = []
            energy = 0
            for a in lines:
                #        print("el elemento de la linea: {}".format(a))
                #        print("energy in that element {} of the line is {}".format(a,myTOB.getEnergy(4, a, phiMax_EM2)))
                energy += myTOB.getEnergy(layer, a, phiMax)
            #        print("the energy in the whole line {}".format(energy))

            #            energyBlock = 0
            #            blocksIndex = []
            #            for b in lines[1]:
            #        #        print("b1 {}".format(b[0]))
            #        #        print("b2 {}".format(b[1]))
            #        #        print("energy in block {} above: {}".format(b[0],myTOB.getEnergy(4, b[0], phiMax_EM2+1)))
            #        #        print("energy in block {} below: {}".format(b[1],myTOB.getEnergy(4, b[1], phiMax_EM2-1)))
            #        #
            #        #        print("la energia de la maxima suma de los bloques {}".format(energyBlock))
            #        #        print("las coordenadas de esos bloques {}".format(blocksIndex))
            #                if ( myTOB.getEnergy(layer, b[0], phiMax_EM2+1) + myTOB.getEnergy(layer, b[1], phiMax_EM2-1) > energyBlock):
            #                    energyBlock = myTOB.getEnergy(layer, b[0], phiMax_EM2+1) + myTOB.getEnergy(layer, b[1], phiMax_EM2-1)
            #                    blocksIndex = [[b[0], phiMax_EM2+1],[b[1], phiMax_EM2-1]]
            #            energy += energyBlock
            if phiMax == 2:
                energy += myTOB.getEnergy(layer, etaMax, phiMax + 1)
                energy += myTOB.getEnergy(layer, etaMax, phiMax - 1)
                blocksIndex = [[etaMax, phiMax + 1], [etaMax, phiMax - 1]]
            elif phimax == 1:
                energy += myTOB.getEnergy(layer, etaMax, phiMax + 1)
                blocksIndex = [[etaMax, phiMax + 1]]
            elif phimax == 3:
                energy += myTOB.getEnergy(layer, etaMax, phiMax - 1)
                blocksIndex = [[etaMax, phiMax - 1]]
            if storeMaxE < energy:
                storeMaxE = energy
                storeIndex = [lines, blocksIndex]
    #        print("layer")
    #        print("in the end, max energy is {}, the index is {}".format(storeMaxE,storeIndex))
    #        print("la maxima energia {}".format(storeMaxE))
    #        print("the index {}".format(storeIndex))
    # now I will use the same coordinate to start in layer EM1

    if layer == 2 or layer == 5 or layer == 6:

        hotcell = (etaMax, phimax)
        # hotcell = (etaMax-1,phimax-1)
        shapes = doShapes(layer, hotcell)
        storeMaxE = -9999
        storeIndex = []
        for shape in shapes:
            energy = 0
            for blocks in shape:
                energy += myTOB.getEnergy(layer, blocks[0] + 1, blocks[1] + 1)
            if energy >= storeMaxE:
                storeMaxE = energy
                storeIndex = shape

    #    if(layer==2 or layer==5 or layer==6):
    #
    #        storeMaxE = -9999
    #        storeIndex = []
    #
    #        e1 = myTOB.getEnergy(layer, 1, 1)+myTOB.getEnergy(layer, 2, 1)+myTOB.getEnergy(layer, 1, 2)+myTOB.getEnergy(layer, 2, 2)
    #        e2 = myTOB.getEnergy(layer, 2, 1)+myTOB.getEnergy(layer, 3, 1)+myTOB.getEnergy(layer, 2, 2)+myTOB.getEnergy(layer, 3, 2)
    #        e3 = myTOB.getEnergy(layer, 1, 3)+myTOB.getEnergy(layer, 2, 3)+myTOB.getEnergy(layer, 1, 2)+myTOB.getEnergy(layer, 2, 2)
    #        e4 = myTOB.getEnergy(layer, 2, 2)+myTOB.getEnergy(layer, 3, 2)+myTOB.getEnergy(layer, 2, 3)+myTOB.getEnergy(layer, 3, 3)
    #        storeMaxE = max(e1,e2,e3,e4)

    return storeMaxE


def TLValgo(myTOB):
    #    etaMax_PS = findAllMax(myTOB)[0][0]
    #    phiMax_PS = findAllMax(myTOB)[0][1]
    #    etaMax_EM1 = findAllMax(myTOB)[1][0]
    #    phiMax_EM1 = findAllMax(myTOB)[1][1]
    #    etaMax_EM2 = findAllMax(myTOB)[2][0]
    #    phiMax_EM2 = findAllMax(myTOB)[2][1]
    #    etaMax_EM3 = findAllMax(myTOB)[3][0]
    #    phiMax_EM3 = findAllMax(myTOB)[3][1]
    #    etaMax_HAD = findAllMax(myTOB)[4][0]
    #    phiMax_HAD = findAllMax(myTOB)[4][1]

    #    print("etaMax_PS {}".format(etaMax_PS))
    #    print("phiMax_PS {}".format(phiMax_PS))
    #    print("etaMax_EM1 {}".format(etaMax_EM1))
    #    print("phiMax_EM1 {}".format(phiMax_EM1))
    #    print("etaMax_EM2 {}".format(etaMax_EM2))
    #    print("phiMax_EM2 {}".format(phiMax_EM2))
    #    print("etaMax_EM3 {}".format(etaMax_EM3))
    #    print("phiMax_EM3 {}".format(phiMax_EM3))
    #    print("etaMax_HAD {}".format(etaMax_HAD))
    #    print("phiMax_HAD {}".format(phiMax_HAD))
    etaMax_EM2 = findMaxEM2(myTOB)[0]
    phiMax_EM2 = findMaxEM2(myTOB)[1]
    etaMax_EM1 = findMaxEM1(myTOB)[0]
    phiMax_EM1 = findMaxEM1(myTOB)[1]
    #    print("old etaMax_PS {}".format(findMax(2,myTOB)[0]))
    #    print("old phiMax_PS {}".format(findMax(2,myTOB)[1]))
    #    print("old etaMax_EM1 {}".format(findMaxEM2(myTOB)[0]))
    #    print("old phiMax_EM1 {}".format(findMaxEM2(myTOB)[1]))
    #    print("old etaMax_EM2 {}".format(findMaxEM2(myTOB)[0]))
    #    print("old phiMax_EM2 {}".format(findMaxEM2(myTOB)[1]))
    #    print("old etaMax_EM3 {}".format(findMax(5,myTOB)[0]))
    #    print("old phiMax_EM3 {}".format(findMax(5,myTOB)[1]))
    #    print("old etaMax_HAD {}".format(findMax(6,myTOB)[0]))
    #    print("old phiMax_HAD {}".format(findMax(6,myTOB)[1]))
    E_EM2 = shapeTLV_fineGran(myTOB, 4, etaMax_EM2, phiMax_EM2)
    #    print("energy in EM2 with TLV algo {}".format(E_EM2))
    # for the moment I'll use the same coordinate as for EM2
    E_EM1 = shapeTLV_fineGran(myTOB, 3, etaMax_EM1, phiMax_EM1)
    etaMax_PS = findMax(2, myTOB)[0]
    phiMax_PS = findMax(2, myTOB)[1]
    E_PS = shapeTLV_fineGran(myTOB, 2, etaMax_PS, phiMax_PS)
    etaMax_EM3 = findMax(5, myTOB)[0]
    phiMax_EM3 = findMax(5, myTOB)[1]
    E_EM3 = shapeTLV_fineGran(myTOB, 5, etaMax_EM3, phiMax_EM3)
    etaMax_HAD = findMax(6, myTOB)[0]
    phiMax_HAD = findMax(6, myTOB)[1]
    E_HAD = shapeTLV_fineGran(myTOB, 6, etaMax_HAD, phiMax_HAD)
    #    print("energy in HAD with TLV algo {}".format(E_HAD))
    #    print("energy in EM3 with TLV algo {}".format(E_EM3))
    #    print("energy in EM1 with TLV algo {}".format(E_EM1))
    #    print("energy in PS with TLV algo {}".format(E_PS))

    return E_EM1 + E_EM2 + E_PS + E_EM3 + E_HAD


##############################################


def findAlignedMax(myTOB):
    # I'll seed with max in EM2(L=4)
    EM2_E = findMax2(4, myTOB)[0]
    EM2_EIndex = findMax2(4, myTOB)[
        1
    ]  # (eta_max,phi_max)=(EM2_EIndex[0],EM2_EIndex[1]) remember to add 4 to eta for layers EM2 and EM1 (and always 1 to phi)

    E_env_EM2 = 0
    for o in range(3):
        for p in range(3):
            if (
                myTOB.getEnergy(4, EM2_EIndex[0] - 1 + 4 + o, EM2_EIndex[1] - 1 + 1 + p)
                > 0
            ):
                E_env_EM2 += (
                    myTOB.getEnergy(
                        4, EM2_EIndex[0] - 6 + 4 + o, EM2_EIndex[1] - 1 + 1 + p
                    )
                    > 0
                )  # not sure about -6
            # if (myTOB.getEnergy(4, EM2_EIndex[0]-1+4+o, EM2_EIndex[1]-1+1+p)>0):
            #    E_env_EM2 +=myTOB.getEnergy(4, EM2_EIndex[0]-1+4+o, EM2_EIndex[1]-1+1+p)>0

    # now I go to the previous layer (EM1) which has the same granularity and check  the hottest cell for all the cells sorrounding the aligned one (including that one)
    maxEM1 = -99
    maxEM1_I = [9999, 9999]
    for o in range(3):
        for p in range(3):
            if maxEM1 < myTOB.getEnergy(
                3, EM2_EIndex[0] + 4 - 1 + o, EM2_EIndex[1] + 1 - 1 + p
            ):
                maxEM1 = myTOB.getEnergy(
                    3, EM2_EIndex[0] + 4 - 1 + o, EM2_EIndex[1] + 1 - 1 + p
                )
                maxEM1_I = [EM2_EIndex[0] + 4 - 1 + o, EM2_EIndex[1] - 1 + 1 + p]
    E_env_EM1 = 0
    for o in range(3):
        for p in range(3):
            if myTOB.getEnergy(3, maxEM1_I[0] - 1 + o, maxEM1_I[1] - 1 + p) > 0:
                E_env_EM1 += myTOB.getEnergy(
                    3, maxEM1_I[0] - 6 + o, maxEM1_I[1] - 1 + p
                )
            # if (myTOB.getEnergy(3, maxEM1_I[0]-1+o, maxEM1_I[1]-1+p)>0):
            #    E_env_EM1 += myTOB.getEnergy(3, maxEM1_I[0]-1+o, maxEM1_I[1]-1+p)

    # now go to the previous layer (the first one, which is PS)
    PS_EIndex_projection = [99, 99]
    if maxEM1_I[0] < 4:
        PS_EIndex_projection[0] = 0
    elif maxEM1_I[0] > 3 and EM2_EIndex[0] < 8:
        PS_EIndex_projection[0] = 1
    elif maxEM1_I[0] > 7:
        PS_EIndex_projection[0] = 2
    PS_EIndex_projection[1] = maxEM1_I[1]

    maxPS = -99
    maxPS_I = [9999, 9999]
    for o in range(3):
        for p in range(3):
            if (
                maxPS
                < myTOB.getEnergy(
                    2,
                    PS_EIndex_projection[0] + 1 - 1 + o,
                    PS_EIndex_projection[1] + 1 - 1 + p,
                )
                and myTOB.getEnergy(
                    1,
                    PS_EIndex_projection[0] + 1 - 1 + o,
                    PS_EIndex_projection[1] + 1 - 1 + p,
                )
                > 0
            ):
                maxPS = myTOB.getEnergy(
                    2,
                    PS_EIndex_projection[0] + 1 - 1 + o,
                    PS_EIndex_projection[1] + 1 - 1 + p,
                )
                maxPS_I = [
                    PS_EIndex_projection[0] + 1 - 1 + o,
                    PS_EIndex_projection[1] - 1 + 1 + p,
                ]

    # now I go to the next layer (EM3), which is 3x3, and first look for the cell with activity which has the smallest distance from the alignement to EM2
    EM3_EIndex_projection = [99, 99]
    if EM2_EIndex[0] < 4:
        EM3_EIndex_projection[0] = 0
    elif EM2_EIndex[0] > 3 and EM2_EIndex[0] < 8:
        EM3_EIndex_projection[0] = 1
    elif EM2_EIndex[0] > 7:
        EM3_EIndex_projection[0] = 2
    EM3_EIndex_projection[1] = EM2_EIndex[1]

    maxEM3 = -99
    maxEM3_I = [9999, 9999]
    for o in range(3):
        for p in range(3):
            if (
                maxEM3
                < myTOB.getEnergy(
                    5,
                    EM3_EIndex_projection[0] + 1 - 1 + o,
                    EM3_EIndex_projection[1] + 1 - 1 + p,
                )
                and myTOB.getEnergy(
                    5,
                    EM3_EIndex_projection[0] + 1 - 1 + o,
                    EM3_EIndex_projection[1] + 1 - 1 + p,
                )
                > 0
            ):
                maxEM3 = myTOB.getEnergy(
                    5,
                    EM3_EIndex_projection[0] + 1 - 1 + o,
                    EM3_EIndex_projection[1] + 1 - 1 + p,
                )
                maxEM3_I = [
                    EM3_EIndex_projection[0] + 1 - 1 + o,
                    EM3_EIndex_projection[1] - 1 + 1 + p,
                ]

    # add three layers:
    totE = 0

    totE = TLValgo(myTOB)
    return totE


def hotNext(myTOB):
    sumE = 0
    iEL2_PS = findMax(2, myTOB)
    iEL3_EM1 = findMax(3, myTOB)
    iEL4_EM2 = findMax(4, myTOB)
    iEL5_EM3 = findMax(5, myTOB)
    iEL6_HAD = findMax(6, myTOB)

    maxEL2_PS = myTOB.getEnergy(2, iEL2_PS[0] + 1, iEL2_PS[1] + 1)
    maxEL3_EM1 = myTOB.getEnergy(3, iEL3_EM1[0] + 4, iEL3_EM1[1] + 1)
    maxEL4_EM2 = myTOB.getEnergy(4, iEL4_EM2[0] + 4, iEL4_EM2[1] + 1)
    maxEL5_EM3 = myTOB.getEnergy(5, iEL5_EM3[0] + 1, iEL5_EM3[1] + 1)
    maxEL6_HAD = myTOB.getEnergy(6, iEL6_HAD[0] + 1, iEL6_HAD[1] + 1)

    sumE = maxEL2_PS + maxEL3_EM1 + maxEL4_EM2 + maxEL5_EM3 + maxEL6_HAD
    sumE += myTOB.getEnergy(2, iEL2_PS[0] + 1 - 1, iEL2_PS[1] + 1 + 1)
    sumE += myTOB.getEnergy(2, iEL2_PS[0] + 1 - 1, iEL2_PS[1] + 1 + 0)
    sumE += myTOB.getEnergy(2, iEL2_PS[0] + 1 - 1, iEL2_PS[1] + 1 - 1)
    sumE += myTOB.getEnergy(2, iEL2_PS[0] + 1 + 0, iEL2_PS[1] + 1 + 1)
    sumE += myTOB.getEnergy(2, iEL2_PS[0] + 1 + 0, iEL2_PS[1] + 1 - 1)
    sumE += myTOB.getEnergy(2, iEL2_PS[0] + 1 + 1, iEL2_PS[1] + 1 + 1)
    sumE += myTOB.getEnergy(2, iEL2_PS[0] + 1 + 1, iEL2_PS[1] + 1 + 0)
    sumE += myTOB.getEnergy(2, iEL2_PS[0] + 1 + 1, iEL2_PS[1] + 1 - 1)
    sumE += myTOB.getEnergy(3, iEL3_EM1[0] + 4 - 1, iEL3_EM1[1] + 1 + 1)
    sumE += myTOB.getEnergy(3, iEL3_EM1[0] + 4 - 1, iEL3_EM1[1] + 1 + 0)
    sumE += myTOB.getEnergy(3, iEL3_EM1[0] + 4 - 1, iEL3_EM1[1] + 1 - 1)
    sumE += myTOB.getEnergy(3, iEL3_EM1[0] + 4 + 0, iEL3_EM1[1] + 1 + 1)
    sumE += myTOB.getEnergy(3, iEL3_EM1[0] + 4 + 0, iEL3_EM1[1] + 1 - 1)
    sumE += myTOB.getEnergy(3, iEL3_EM1[0] + 4 + 1, iEL3_EM1[1] + 1 + 1)
    sumE += myTOB.getEnergy(3, iEL3_EM1[0] + 4 + 1, iEL3_EM1[1] + 1 + 0)
    sumE += myTOB.getEnergy(3, iEL3_EM1[0] + 4 + 1, iEL3_EM1[1] + 1 - 1)
    sumE += myTOB.getEnergy(4, iEL4_EM2[0] + 4 - 1, iEL4_EM2[1] + 1 + 1)
    sumE += myTOB.getEnergy(4, iEL4_EM2[0] + 4 - 1, iEL4_EM2[1] + 1 + 0)
    sumE += myTOB.getEnergy(4, iEL4_EM2[0] + 4 - 1, iEL4_EM2[1] + 1 - 1)
    sumE += myTOB.getEnergy(4, iEL4_EM2[0] + 4 + 0, iEL4_EM2[1] + 1 + 1)
    sumE += myTOB.getEnergy(4, iEL4_EM2[0] + 4 + 0, iEL4_EM2[1] + 1 - 1)
    sumE += myTOB.getEnergy(4, iEL4_EM2[0] + 4 + 1, iEL4_EM2[1] + 1 + 1)
    sumE += myTOB.getEnergy(4, iEL4_EM2[0] + 4 + 1, iEL4_EM2[1] + 1 + 0)
    sumE += myTOB.getEnergy(4, iEL4_EM2[0] + 4 + 1, iEL4_EM2[1] + 1 - 1)
    sumE += myTOB.getEnergy(5, iEL5_EM3[0] + 1 - 1, iEL5_EM3[1] + 1 + 1)
    sumE += myTOB.getEnergy(5, iEL5_EM3[0] + 1 - 1, iEL5_EM3[1] + 1 + 0)
    sumE += myTOB.getEnergy(5, iEL5_EM3[0] + 1 - 1, iEL5_EM3[1] + 1 - 1)
    sumE += myTOB.getEnergy(5, iEL5_EM3[0] + 1 + 0, iEL5_EM3[1] + 1 + 1)
    sumE += myTOB.getEnergy(5, iEL5_EM3[0] + 1 + 0, iEL5_EM3[1] + 1 - 1)
    sumE += myTOB.getEnergy(5, iEL5_EM3[0] + 1 + 1, iEL5_EM3[1] + 1 + 1)
    sumE += myTOB.getEnergy(5, iEL5_EM3[0] + 1 + 1, iEL5_EM3[1] + 1 + 0)
    sumE += myTOB.getEnergy(5, iEL5_EM3[0] + 1 + 1, iEL5_EM3[1] + 1 - 1)
    sumE += myTOB.getEnergy(6, iEL6_HAD[0] + 1 - 1, iEL6_HAD[1] + 1 + 1)
    sumE += myTOB.getEnergy(6, iEL6_HAD[0] + 1 - 1, iEL6_HAD[1] + 1 + 0)
    sumE += myTOB.getEnergy(6, iEL6_HAD[0] + 1 - 1, iEL6_HAD[1] + 1 - 1)
    sumE += myTOB.getEnergy(6, iEL6_HAD[0] + 1 + 0, iEL6_HAD[1] + 1 + 1)
    sumE += myTOB.getEnergy(6, iEL6_HAD[0] + 1 + 0, iEL6_HAD[1] + 1 - 1)
    sumE += myTOB.getEnergy(6, iEL6_HAD[0] + 1 + 1, iEL6_HAD[1] + 1 + 1)
    sumE += myTOB.getEnergy(6, iEL6_HAD[0] + 1 + 1, iEL6_HAD[1] + 1 + 0)
    sumE += myTOB.getEnergy(6, iEL6_HAD[0] + 1 + 1, iEL6_HAD[1] + 1 - 1)

    return sumE


def hotNext75(myTOB, factor):
    sumE = 0
    nextE = 0
    iEL2_PS = findMax(2, myTOB)
    iEL3_EM1 = findMax(3, myTOB)
    iEL4_EM2 = findMax(4, myTOB)
    iEL5_EM3 = findMax(5, myTOB)
    iEL6_HAD = findMax(6, myTOB)

    maxEL2_PS = myTOB.getEnergy(2, iEL2_PS[0] + 1, iEL2_PS[1] + 1)
    maxEL3_EM1 = myTOB.getEnergy(3, iEL3_EM1[0] + 4, iEL3_EM1[1] + 1)
    maxEL4_EM2 = myTOB.getEnergy(4, iEL4_EM2[0] + 4, iEL4_EM2[1] + 1)
    maxEL5_EM3 = myTOB.getEnergy(5, iEL5_EM3[0] + 1, iEL5_EM3[1] + 1)
    maxEL6_HAD = myTOB.getEnergy(6, iEL6_HAD[0] + 1, iEL6_HAD[1] + 1)

    sumE = maxEL2_PS + maxEL3_EM1 + maxEL4_EM2 + maxEL5_EM3 + maxEL6_HAD
    # add layer 2 = PS
    nextE = myTOB.getEnergy(2, iEL2_PS[0] + 1 - 1, iEL2_PS[1] + 1 + 1)
    if nextE > factor * maxEL2_PS:
        sumE += nextE
    nextE = myTOB.getEnergy(2, iEL2_PS[0] + 1 - 1, iEL2_PS[1] + 1 + 0)
    if nextE > factor * maxEL2_PS:
        sumE += nextE
    nextE = myTOB.getEnergy(2, iEL2_PS[0] + 1 - 1, iEL2_PS[1] + 1 - 1)
    if nextE > factor * maxEL2_PS:
        sumE += nextE
    nextE = myTOB.getEnergy(2, iEL2_PS[0] + 1 + 0, iEL2_PS[1] + 1 + 1)
    if nextE > factor * maxEL2_PS:
        sumE += nextE
    nextE = myTOB.getEnergy(2, iEL2_PS[0] + 1 + 0, iEL2_PS[1] + 1 - 1)
    if nextE > factor * maxEL2_PS:
        sumE += nextE
    nextE = myTOB.getEnergy(2, iEL2_PS[0] + 1 + 1, iEL2_PS[1] + 1 + 1)
    if nextE > factor * maxEL2_PS:
        sumE += nextE
    nextE = myTOB.getEnergy(2, iEL2_PS[0] + 1 + 1, iEL2_PS[1] + 1 + 0)
    if nextE > factor * maxEL2_PS:
        sumE += nextE
    nextE = myTOB.getEnergy(2, iEL2_PS[0] + 1 + 1, iEL2_PS[1] + 1 - 1)
    if nextE > factor * maxEL2_PS:
        sumE += nextE
    # add layer 3 = EM1
    nextE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 - 1, iEL3_EM1[1] + 1 + 1)
    if nextE > factor * maxEL3_EM1:
        sumE += nextE
    nextE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 - 1, iEL3_EM1[1] + 1 + 0)
    if nextE > factor * maxEL3_EM1:
        sumE += nextE
    nextE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 - 1, iEL3_EM1[1] + 1 - 1)
    if nextE > factor * maxEL3_EM1:
        sumE += nextE
    nextE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 + 0, iEL3_EM1[1] + 1 + 1)
    if nextE > factor * maxEL3_EM1:
        sumE += nextE
    nextE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 + 0, iEL3_EM1[1] + 1 - 1)
    if nextE > factor * maxEL3_EM1:
        sumE += nextE
    nextE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 + 1, iEL3_EM1[1] + 1 + 1)
    if nextE > factor * maxEL3_EM1:
        sumE += nextE
    nextE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 + 1, iEL3_EM1[1] + 1 + 0)
    if nextE > factor * maxEL3_EM1:
        sumE += nextE
    nextE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 + 1, iEL3_EM1[1] + 1 - 1)
    if nextE > factor * maxEL3_EM1:
        sumE += nextE
    # add layer 4 = EM2
    nextE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 - 1, iEL4_EM2[1] + 1 + 1)
    if nextE > factor * maxEL4_EM2:
        sumE += nextE
    nextE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 - 1, iEL4_EM2[1] + 1 + 0)
    if nextE > factor * maxEL4_EM2:
        sumE += nextE
    nextE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 - 1, iEL4_EM2[1] + 1 - 1)
    if nextE > factor * maxEL4_EM2:
        sumE += nextE
    nextE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 + 0, iEL4_EM2[1] + 1 + 1)
    if nextE > factor * maxEL4_EM2:
        sumE += nextE
    nextE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 + 0, iEL4_EM2[1] + 1 - 1)
    if nextE > factor * maxEL4_EM2:
        sumE += nextE
    nextE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 + 1, iEL4_EM2[1] + 1 + 1)
    if nextE > factor * maxEL4_EM2:
        sumE += nextE
    nextE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 + 1, iEL4_EM2[1] + 1 + 0)
    if nextE > factor * maxEL4_EM2:
        sumE += nextE
    nextE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 + 1, iEL4_EM2[1] + 1 - 1)
    if nextE > factor * maxEL4_EM2:
        sumE += nextE
    # add layer 5 = EM3
    nextE = myTOB.getEnergy(5, iEL5_EM3[0] + 1 - 1, iEL5_EM3[1] + 1 + 1)
    if nextE > factor * maxEL5_EM3:
        sumE += nextE
    nextE = myTOB.getEnergy(5, iEL5_EM3[0] + 1 - 1, iEL5_EM3[1] + 1 + 0)
    if nextE > factor * maxEL5_EM3:
        sumE += nextE
    nextE = myTOB.getEnergy(5, iEL5_EM3[0] + 1 - 1, iEL5_EM3[1] + 1 - 1)
    if nextE > factor * maxEL5_EM3:
        sumE += nextE
    nextE = myTOB.getEnergy(5, iEL5_EM3[0] + 1 + 0, iEL5_EM3[1] + 1 + 1)
    if nextE > factor * maxEL5_EM3:
        sumE += nextE
    nextE = myTOB.getEnergy(5, iEL5_EM3[0] + 1 + 0, iEL5_EM3[1] + 1 - 1)
    if nextE > factor * maxEL5_EM3:
        sumE += nextE
    nextE = myTOB.getEnergy(5, iEL5_EM3[0] + 1 + 1, iEL5_EM3[1] + 1 + 1)
    if nextE > factor * maxEL5_EM3:
        sumE += nextE
    nextE = myTOB.getEnergy(5, iEL5_EM3[0] + 1 + 1, iEL5_EM3[1] + 1 + 0)
    if nextE > factor * maxEL5_EM3:
        sumE += nextE
    nextE = myTOB.getEnergy(5, iEL5_EM3[0] + 1 + 1, iEL5_EM3[1] + 1 - 1)
    if nextE > factor * maxEL5_EM3:
        sumE += nextE
    # add layer 6 = HAD
    nextE = myTOB.getEnergy(6, iEL6_HAD[0] + 1 - 1, iEL6_HAD[1] + 1 + 1)
    if nextE > factor * maxEL6_HAD:
        sumE += nextE
    nextE = myTOB.getEnergy(6, iEL6_HAD[0] + 1 - 1, iEL6_HAD[1] + 1 + 0)
    if nextE > factor * maxEL6_HAD:
        sumE += nextE
    nextE = myTOB.getEnergy(6, iEL6_HAD[0] + 1 - 1, iEL6_HAD[1] + 1 - 1)
    if nextE > factor * maxEL6_HAD:
        sumE += nextE
    nextE = myTOB.getEnergy(6, iEL6_HAD[0] + 1 + 0, iEL6_HAD[1] + 1 + 1)
    if nextE > factor * maxEL6_HAD:
        sumE += nextE
    nextE = myTOB.getEnergy(6, iEL6_HAD[0] + 1 + 0, iEL6_HAD[1] + 1 - 1)
    if nextE > factor * maxEL6_HAD:
        sumE += nextE
    nextE = myTOB.getEnergy(6, iEL6_HAD[0] + 1 + 1, iEL6_HAD[1] + 1 + 1)
    if nextE > factor * maxEL6_HAD:
        sumE += nextE
    nextE = myTOB.getEnergy(6, iEL6_HAD[0] + 1 + 1, iEL6_HAD[1] + 1 + 0)
    if nextE > factor * maxEL6_HAD:
        sumE += nextE
    nextE = myTOB.getEnergy(6, iEL6_HAD[0] + 1 + 1, iEL6_HAD[1] + 1 - 1)
    if nextE > factor * maxEL6_HAD:
        sumE += nextE

    return sumE


def hotNextNf25(myTOB, factor, Nfactor):
    sumE = 0
    nextE = 0
    iEL2_PS = findMax(2, myTOB)
    iEL3_EM1 = findMax(3, myTOB)
    iEL4_EM2 = findMax(4, myTOB)
    iEL5_EM3 = findMax(5, myTOB)
    iEL6_HAD = findMax(6, myTOB)

    maxEL2_PS = myTOB.getEnergy(2, iEL2_PS[0] + 1, iEL2_PS[1] + 1)
    maxEL3_EM1 = myTOB.getEnergy(3, iEL3_EM1[0] + 4, iEL3_EM1[1] + 1)
    maxEL4_EM2 = myTOB.getEnergy(4, iEL4_EM2[0] + 4, iEL4_EM2[1] + 1)
    maxEL5_EM3 = myTOB.getEnergy(5, iEL5_EM3[0] + 1, iEL5_EM3[1] + 1)
    maxEL6_HAD = myTOB.getEnergy(6, iEL6_HAD[0] + 1, iEL6_HAD[1] + 1)

    sumE = maxEL2_PS + maxEL3_EM1 + maxEL4_EM2 + maxEL5_EM3 + maxEL6_HAD

    nextE = myTOB.getEnergy(2, iEL2_PS[0] + 1 - 1, iEL2_PS[1] + 1 + 1)
    if nextE > factor * maxEL2_PS:
        sumE += nextE
    nextE = myTOB.getEnergy(2, iEL2_PS[0] + 1 - 1, iEL2_PS[1] + 1 + 0)
    if nextE > factor * maxEL2_PS:
        sumE += nextE
    nextE = myTOB.getEnergy(2, iEL2_PS[0] + 1 - 1, iEL2_PS[1] + 1 - 1)
    if nextE > factor * maxEL2_PS:
        sumE += nextE
    nextE = myTOB.getEnergy(2, iEL2_PS[0] + 1 + 0, iEL2_PS[1] + 1 + 1)
    if nextE > factor * maxEL2_PS:
        sumE += nextE
    nextE = myTOB.getEnergy(2, iEL2_PS[0] + 1 + 0, iEL2_PS[1] + 1 - 1)
    if nextE > factor * maxEL2_PS:
        sumE += nextE
    nextE = myTOB.getEnergy(2, iEL2_PS[0] + 1 + 1, iEL2_PS[1] + 1 + 1)
    if nextE > factor * maxEL2_PS:
        sumE += nextE
    nextE = myTOB.getEnergy(2, iEL2_PS[0] + 1 + 1, iEL2_PS[1] + 1 + 0)
    if nextE > factor * maxEL2_PS:
        sumE += nextE
    nextE = myTOB.getEnergy(2, iEL2_PS[0] + 1 + 1, iEL2_PS[1] + 1 - 1)
    if nextE > factor * maxEL2_PS:
        sumE += nextE

    nextNE = myTOB.getEnergy(2, iEL2_PS[0] + 1 - 2, iEL2_PS[1] + 1 + 1)
    if nextNE > Nfactor * maxEL2_PS:
        sumE += nextNE
    nextNE = myTOB.getEnergy(2, iEL2_PS[0] + 1 - 2, iEL2_PS[1] + 1 + 2)
    if nextNE > Nfactor * maxEL2_PS:
        sumE += nextNE
    nextNE = myTOB.getEnergy(2, iEL2_PS[0] + 1 - 1, iEL2_PS[1] + 1 + 2)
    if nextNE > Nfactor * maxEL2_PS:
        sumE += nextNE
    nextNE = myTOB.getEnergy(2, iEL2_PS[0] + 1 - 2, iEL2_PS[1] + 1 + 0)
    if nextNE > Nfactor * maxEL2_PS:
        sumE += nextNE
    nextNE = myTOB.getEnergy(2, iEL2_PS[0] + 1 - 2, iEL2_PS[1] + 1 - 1)
    if nextNE > Nfactor * maxEL2_PS:
        sumE += nextNE
    nextNE = myTOB.getEnergy(2, iEL2_PS[0] + 1 - 2, iEL2_PS[1] + 1 - 2)
    if nextNE > Nfactor * maxEL2_PS:
        sumE += nextNE
    nextNE = myTOB.getEnergy(2, iEL2_PS[0] + 1 - 1, iEL2_PS[1] + 1 - 2)
    if nextNE > Nfactor * maxEL2_PS:
        sumE += nextNE
    nextNE = myTOB.getEnergy(2, iEL2_PS[0] + 1 + 0, iEL2_PS[1] + 1 + 2)
    if nextNE > Nfactor * maxEL2_PS:
        sumE += nextNE
    nextNE = myTOB.getEnergy(2, iEL2_PS[0] + 1 + 0, iEL2_PS[1] + 1 - 2)
    if nextNE > Nfactor * maxEL2_PS:
        sumE += nextNE
    nextNE = myTOB.getEnergy(2, iEL2_PS[0] + 1 + 1, iEL2_PS[1] + 1 + 2)
    if nextNE > Nfactor * maxEL2_PS:
        sumE += nextNE
    nextNE = myTOB.getEnergy(2, iEL2_PS[0] + 1 + 2, iEL2_PS[1] + 1 + 2)
    if nextNE > Nfactor * maxEL2_PS:
        sumE += nextNE
    nextNE = myTOB.getEnergy(2, iEL2_PS[0] + 1 + 2, iEL2_PS[1] + 1 + 1)
    if nextNE > Nfactor * maxEL2_PS:
        sumE += nextNE
    nextNE = myTOB.getEnergy(2, iEL2_PS[0] + 1 + 2, iEL2_PS[1] + 1 + 0)
    if nextNE > Nfactor * maxEL2_PS:
        sumE += nextNE
    nextNE = myTOB.getEnergy(2, iEL2_PS[0] + 1 + 2, iEL2_PS[1] + 1 - 1)
    if nextNE > Nfactor * maxEL2_PS:
        sumE += nextNE
    nextNE = myTOB.getEnergy(2, iEL2_PS[0] + 1 + 2, iEL2_PS[1] + 1 - 2)
    if nextNE > Nfactor * maxEL2_PS:
        sumE += nextNE
    nextNE = myTOB.getEnergy(2, iEL2_PS[0] + 1 + 1, iEL2_PS[1] + 1 - 2)
    if nextNE > Nfactor * maxEL2_PS:
        sumE += nextNE

    # add layer 3 = EM1
    nextE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 - 1, iEL3_EM1[1] + 1 + 1)
    if nextE > factor * maxEL3_EM1:
        sumE += nextE
    nextE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 - 1, iEL3_EM1[1] + 1 + 0)
    if nextE > factor * maxEL3_EM1:
        sumE += nextE
    nextE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 - 1, iEL3_EM1[1] + 1 - 1)
    if nextE > factor * maxEL3_EM1:
        sumE += nextE
    nextE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 + 0, iEL3_EM1[1] + 1 + 1)
    if nextE > factor * maxEL3_EM1:
        sumE += nextE
    nextE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 + 0, iEL3_EM1[1] + 1 - 1)
    if nextE > factor * maxEL3_EM1:
        sumE += nextE
    nextE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 + 1, iEL3_EM1[1] + 1 + 1)
    if nextE > factor * maxEL3_EM1:
        sumE += nextE
    nextE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 + 1, iEL3_EM1[1] + 1 + 0)
    if nextE > factor * maxEL3_EM1:
        sumE += nextE
    nextE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 + 1, iEL3_EM1[1] + 1 - 1)
    if nextE > factor * maxEL3_EM1:
        sumE += nextE
    nextE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 + 2, iEL3_EM1[1] + 1 - 1)
    if nextE > factor * maxEL3_EM1:
        sumE += nextE
    nextE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 + 3, iEL3_EM1[1] + 1 - 1)
    if nextE > factor * maxEL3_EM1:
        sumE += nextE
    nextE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 + 4, iEL3_EM1[1] + 1 - 1)
    if nextE > factor * maxEL3_EM1:
        sumE += nextE
    nextE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 + 5, iEL3_EM1[1] + 1 - 1)
    if nextE > factor * maxEL3_EM1:
        sumE += nextE
    nextE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 - 2, iEL3_EM1[1] + 1 - 1)
    if nextE > factor * maxEL3_EM1:
        sumE += nextE
    nextE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 - 3, iEL3_EM1[1] + 1 - 1)
    if nextE > factor * maxEL3_EM1:
        sumE += nextE
    nextE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 - 4, iEL3_EM1[1] + 1 - 1)
    if nextE > factor * maxEL3_EM1:
        sumE += nextE
    nextE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 - 5, iEL3_EM1[1] + 1 - 1)
    if nextE > factor * maxEL3_EM1:
        sumE += nextE
    nextE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 - 6, iEL3_EM1[1] + 1 - 1)
    if nextE > factor * maxEL3_EM1:
        sumE += nextE
    nextE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 + 2, iEL3_EM1[1] + 1 + 1)
    if nextE > factor * maxEL3_EM1:
        sumE += nextE
    nextE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 + 3, iEL3_EM1[1] + 1 + 1)
    if nextE > factor * maxEL3_EM1:
        sumE += nextE
    nextE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 + 4, iEL3_EM1[1] + 1 + 1)
    if nextE > factor * maxEL3_EM1:
        sumE += nextE
    nextE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 + 5, iEL3_EM1[1] + 1 + 1)
    if nextE > factor * maxEL3_EM1:
        sumE += nextE
    nextE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 - 2, iEL3_EM1[1] + 1 + 1)
    if nextE > factor * maxEL3_EM1:
        sumE += nextE
    nextE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 - 3, iEL3_EM1[1] + 1 + 1)
    if nextE > factor * maxEL3_EM1:
        sumE += nextE
    nextE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 - 4, iEL3_EM1[1] + 1 + 1)
    if nextE > factor * maxEL3_EM1:
        sumE += nextE
    nextE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 - 5, iEL3_EM1[1] + 1 + 1)
    if nextE > factor * maxEL3_EM1:
        sumE += nextE
    nextE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 - 6, iEL3_EM1[1] + 1 + 1)
    if nextE > factor * maxEL3_EM1:
        sumE += nextE
    nextE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 + 2, iEL3_EM1[1] + 1 + 0)
    if nextE > factor * maxEL3_EM1:
        sumE += nextE
    nextE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 + 3, iEL3_EM1[1] + 1 + 0)
    if nextE > factor * maxEL3_EM1:
        sumE += nextE
    nextE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 + 4, iEL3_EM1[1] + 1 + 0)
    if nextE > factor * maxEL3_EM1:
        sumE += nextE
    nextE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 + 5, iEL3_EM1[1] + 1 + 0)
    if nextE > factor * maxEL3_EM1:
        sumE += nextE
    nextE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 - 2, iEL3_EM1[1] + 1 + 0)
    if nextE > factor * maxEL3_EM1:
        sumE += nextE
    nextE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 - 3, iEL3_EM1[1] + 1 + 0)
    if nextE > factor * maxEL3_EM1:
        sumE += nextE
    nextE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 - 4, iEL3_EM1[1] + 1 + 0)
    if nextE > factor * maxEL3_EM1:
        sumE += nextE
    nextE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 - 5, iEL3_EM1[1] + 1 + 0)
    if nextE > factor * maxEL3_EM1:
        sumE += nextE
    nextE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 - 6, iEL3_EM1[1] + 1 + 0)
    if nextE > factor * maxEL3_EM1:
        sumE += nextE

    nextNE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 + 6, iEL3_EM1[1] + 1 - 1)
    if nextNE > Nfactor * maxEL3_EM1:
        sumE += nextNE
    nextNE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 - 7, iEL3_EM1[1] + 1 - 1)
    if nextNE > Nfactor * maxEL3_EM1:
        sumE += nextNE
    nextNE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 + 6, iEL3_EM1[1] + 1 + 1)
    if nextNE > Nfactor * maxEL3_EM1:
        sumE += nextNE
    nextNE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 - 7, iEL3_EM1[1] + 1 + 1)
    if nextNE > Nfactor * maxEL3_EM1:
        sumE += nextNE
    nextNE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 + 6, iEL3_EM1[1] + 1 + 0)
    if nextNE > Nfactor * maxEL3_EM1:
        sumE += nextNE
    nextNE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 - 7, iEL3_EM1[1] + 1 + 0)
    if nextNE > Nfactor * maxEL3_EM1:
        sumE += nextNE
    nextNE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 + 6, iEL3_EM1[1] + 1 - 2)
    if nextNE > Nfactor * maxEL3_EM1:
        sumE += nextNE
    nextNE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 + 5, iEL3_EM1[1] + 1 - 2)
    if nextNE > Nfactor * maxEL3_EM1:
        sumE += nextNE
    nextNE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 + 4, iEL3_EM1[1] + 1 - 2)
    if nextNE > Nfactor * maxEL3_EM1:
        sumE += nextNE
    nextNE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 + 3, iEL3_EM1[1] + 1 - 2)
    if nextNE > Nfactor * maxEL3_EM1:
        sumE += nextNE
    nextNE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 + 1, iEL3_EM1[1] + 1 - 2)
    if nextNE > Nfactor * maxEL3_EM1:
        sumE += nextNE
    nextNE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 + 2, iEL3_EM1[1] + 1 - 2)
    if nextNE > Nfactor * maxEL3_EM1:
        sumE += nextNE
    nextNE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 + 0, iEL3_EM1[1] + 1 - 2)
    if nextNE > Nfactor * maxEL3_EM1:
        sumE += nextNE
    nextNE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 - 1, iEL3_EM1[1] + 1 - 2)
    if nextNE > Nfactor * maxEL3_EM1:
        sumE += nextNE
    nextNE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 - 2, iEL3_EM1[1] + 1 - 2)
    if nextNE > Nfactor * maxEL3_EM1:
        sumE += nextNE
    nextNE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 - 3, iEL3_EM1[1] + 1 - 2)
    if nextNE > Nfactor * maxEL3_EM1:
        sumE += nextNE
    nextNE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 - 4, iEL3_EM1[1] + 1 - 2)
    if nextNE > Nfactor * maxEL3_EM1:
        sumE += nextNE
    nextNE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 - 5, iEL3_EM1[1] + 1 - 2)
    if nextNE > Nfactor * maxEL3_EM1:
        sumE += nextNE
    nextNE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 - 6, iEL3_EM1[1] + 1 - 2)
    if nextNE > Nfactor * maxEL3_EM1:
        sumE += nextNE
    nextNE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 - 7, iEL3_EM1[1] + 1 - 2)
    if nextNE > Nfactor * maxEL3_EM1:
        sumE += nextNE
    nextNE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 + 6, iEL3_EM1[1] + 1 + 2)
    if nextNE > Nfactor * maxEL3_EM1:
        sumE += nextNE
    nextNE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 + 5, iEL3_EM1[1] + 1 + 2)
    if nextNE > Nfactor * maxEL3_EM1:
        sumE += nextNE
    nextNE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 + 4, iEL3_EM1[1] + 1 + 2)
    if nextNE > Nfactor * maxEL3_EM1:
        sumE += nextNE
    nextNE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 + 3, iEL3_EM1[1] + 1 + 2)
    if nextNE > Nfactor * maxEL3_EM1:
        sumE += nextNE
    nextNE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 + 2, iEL3_EM1[1] + 1 + 2)
    if nextNE > Nfactor * maxEL3_EM1:
        sumE += nextNE
    nextNE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 + 1, iEL3_EM1[1] + 1 + 2)
    if nextNE > Nfactor * maxEL3_EM1:
        sumE += nextNE
    nextNE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 + 0, iEL3_EM1[1] + 1 + 2)
    if nextNE > Nfactor * maxEL3_EM1:
        sumE += nextNE
    nextNE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 - 1, iEL3_EM1[1] + 1 + 2)
    if nextNE > Nfactor * maxEL3_EM1:
        sumE += nextNE
    nextNE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 - 2, iEL3_EM1[1] + 1 + 2)
    if nextNE > Nfactor * maxEL3_EM1:
        sumE += nextNE
    nextNE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 - 3, iEL3_EM1[1] + 1 + 2)
    if nextNE > Nfactor * maxEL3_EM1:
        sumE += nextNE
    nextNE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 - 4, iEL3_EM1[1] + 1 + 2)
    if nextNE > Nfactor * maxEL3_EM1:
        sumE += nextNE
    nextNE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 - 5, iEL3_EM1[1] + 1 + 2)
    if nextNE > Nfactor * maxEL3_EM1:
        sumE += nextNE
    nextNE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 - 6, iEL3_EM1[1] + 1 + 2)
    if nextNE > Nfactor * maxEL3_EM1:
        sumE += nextNE
    nextNE = myTOB.getEnergy(3, iEL3_EM1[0] + 4 - 7, iEL3_EM1[1] + 1 + 2)
    if nextNE > Nfactor * maxEL3_EM1:
        sumE += nextNE

    # add layer 4 = EM2
    nextE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 - 1, iEL4_EM2[1] + 1 + 1)
    if nextE > factor * maxEL4_EM2:
        sumE += nextE
    nextE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 - 1, iEL4_EM2[1] + 1 + 0)
    if nextE > factor * maxEL4_EM2:
        sumE += nextE
    nextE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 - 1, iEL4_EM2[1] + 1 - 1)
    if nextE > factor * maxEL4_EM2:
        sumE += nextE
    nextE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 + 0, iEL4_EM2[1] + 1 + 1)
    if nextE > factor * maxEL4_EM2:
        sumE += nextE
    nextE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 + 0, iEL4_EM2[1] + 1 - 1)
    if nextE > factor * maxEL4_EM2:
        sumE += nextE
    nextE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 + 1, iEL4_EM2[1] + 1 + 1)
    if nextE > factor * maxEL4_EM2:
        sumE += nextE
    nextE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 + 1, iEL4_EM2[1] + 1 + 0)
    if nextE > factor * maxEL4_EM2:
        sumE += nextE
    nextE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 + 1, iEL4_EM2[1] + 1 - 1)
    if nextE > factor * maxEL4_EM2:
        sumE += nextE
    nextE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 + 2, iEL4_EM2[1] + 1 - 1)
    if nextE > factor * maxEL4_EM2:
        sumE += nextE
    nextE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 + 3, iEL4_EM2[1] + 1 - 1)
    if nextE > factor * maxEL4_EM2:
        sumE += nextE
    nextE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 + 4, iEL4_EM2[1] + 1 - 1)
    if nextE > factor * maxEL4_EM2:
        sumE += nextE
    nextE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 + 5, iEL4_EM2[1] + 1 - 1)
    if nextE > factor * maxEL4_EM2:
        sumE += nextE
    nextE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 - 2, iEL4_EM2[1] + 1 - 1)
    if nextE > factor * maxEL4_EM2:
        sumE += nextE
    nextE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 - 3, iEL4_EM2[1] + 1 - 1)
    if nextE > factor * maxEL4_EM2:
        sumE += nextE
    nextE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 - 4, iEL4_EM2[1] + 1 - 1)
    if nextE > factor * maxEL4_EM2:
        sumE += nextE
    nextE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 - 5, iEL4_EM2[1] + 1 - 1)
    if nextE > factor * maxEL4_EM2:
        sumE += nextE
    nextE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 - 6, iEL4_EM2[1] + 1 - 1)
    if nextE > factor * maxEL4_EM2:
        sumE += nextE
    nextE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 + 2, iEL4_EM2[1] + 1 + 1)
    if nextE > factor * maxEL4_EM2:
        sumE += nextE
    nextE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 + 3, iEL4_EM2[1] + 1 + 1)
    if nextE > factor * maxEL4_EM2:
        sumE += nextE
    nextE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 + 4, iEL4_EM2[1] + 1 + 1)
    if nextE > factor * maxEL4_EM2:
        sumE += nextE
    nextE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 + 5, iEL4_EM2[1] + 1 + 1)
    if nextE > factor * maxEL4_EM2:
        sumE += nextE
    nextE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 - 2, iEL4_EM2[1] + 1 + 1)
    if nextE > factor * maxEL4_EM2:
        sumE += nextE
    nextE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 - 3, iEL4_EM2[1] + 1 + 1)
    if nextE > factor * maxEL4_EM2:
        sumE += nextE
    nextE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 - 4, iEL4_EM2[1] + 1 + 1)
    if nextE > factor * maxEL4_EM2:
        sumE += nextE
    nextE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 - 5, iEL4_EM2[1] + 1 + 1)
    if nextE > factor * maxEL4_EM2:
        sumE += nextE
    nextE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 - 6, iEL4_EM2[1] + 1 + 1)
    if nextE > factor * maxEL4_EM2:
        sumE += nextE
    nextE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 + 2, iEL4_EM2[1] + 1 + 0)
    if nextE > factor * maxEL4_EM2:
        sumE += nextE
    nextE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 + 3, iEL4_EM2[1] + 1 + 0)
    if nextE > factor * maxEL4_EM2:
        sumE += nextE
    nextE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 + 4, iEL4_EM2[1] + 1 + 0)
    if nextE > factor * maxEL4_EM2:
        sumE += nextE
    nextE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 + 5, iEL4_EM2[1] + 1 + 0)
    if nextE > factor * maxEL4_EM2:
        sumE += nextE
    nextE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 - 2, iEL4_EM2[1] + 1 + 0)
    if nextE > factor * maxEL4_EM2:
        sumE += nextE
    nextE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 - 3, iEL4_EM2[1] + 1 + 0)
    if nextE > factor * maxEL4_EM2:
        sumE += nextE
    nextE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 - 4, iEL4_EM2[1] + 1 + 0)
    if nextE > factor * maxEL4_EM2:
        sumE += nextE
    nextE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 - 5, iEL4_EM2[1] + 1 + 0)
    if nextE > factor * maxEL4_EM2:
        sumE += nextE
    nextE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 - 6, iEL4_EM2[1] + 1 + 0)
    if nextE > factor * maxEL4_EM2:
        sumE += nextE

    nextNE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 + 6, iEL4_EM2[1] + 1 - 1)
    if nextNE > Nfactor * maxEL4_EM2:
        sumE += nextNE
    nextNE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 - 7, iEL4_EM2[1] + 1 - 1)
    if nextNE > Nfactor * maxEL4_EM2:
        sumE += nextNE
    nextNE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 + 6, iEL4_EM2[1] + 1 + 1)
    if nextNE > Nfactor * maxEL4_EM2:
        sumE += nextNE
    nextNE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 - 7, iEL4_EM2[1] + 1 + 1)
    if nextNE > Nfactor * maxEL4_EM2:
        sumE += nextNE
    nextNE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 + 6, iEL4_EM2[1] + 1 + 0)
    if nextNE > Nfactor * maxEL4_EM2:
        sumE += nextNE
    nextNE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 - 7, iEL4_EM2[1] + 1 + 0)
    if nextNE > Nfactor * maxEL4_EM2:
        sumE += nextNE
    nextNE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 + 6, iEL4_EM2[1] + 1 - 2)
    if nextNE > Nfactor * maxEL4_EM2:
        sumE += nextNE
    nextNE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 + 5, iEL4_EM2[1] + 1 - 2)
    if nextNE > Nfactor * maxEL4_EM2:
        sumE += nextNE
    nextNE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 + 4, iEL4_EM2[1] + 1 - 2)
    if nextNE > Nfactor * maxEL4_EM2:
        sumE += nextNE
    nextNE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 + 3, iEL4_EM2[1] + 1 - 2)
    if nextNE > Nfactor * maxEL4_EM2:
        sumE += nextNE
    nextNE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 + 1, iEL4_EM2[1] + 1 - 2)
    if nextNE > Nfactor * maxEL4_EM2:
        sumE += nextNE
    nextNE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 + 2, iEL4_EM2[1] + 1 - 2)
    if nextNE > Nfactor * maxEL4_EM2:
        sumE += nextNE
    nextNE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 + 0, iEL4_EM2[1] + 1 - 2)
    if nextNE > Nfactor * maxEL4_EM2:
        sumE += nextNE
    nextNE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 - 1, iEL4_EM2[1] + 1 - 2)
    if nextNE > Nfactor * maxEL4_EM2:
        sumE += nextNE
    nextNE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 - 2, iEL4_EM2[1] + 1 - 2)
    if nextNE > Nfactor * maxEL4_EM2:
        sumE += nextNE
    nextNE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 - 3, iEL4_EM2[1] + 1 - 2)
    if nextNE > Nfactor * maxEL4_EM2:
        sumE += nextNE
    nextNE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 - 4, iEL4_EM2[1] + 1 - 2)
    if nextNE > Nfactor * maxEL4_EM2:
        sumE += nextNE
    nextNE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 - 5, iEL4_EM2[1] + 1 - 2)
    if nextNE > Nfactor * maxEL4_EM2:
        sumE += nextNE
    nextNE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 - 6, iEL4_EM2[1] + 1 - 2)
    if nextNE > Nfactor * maxEL4_EM2:
        sumE += nextNE
    nextNE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 - 7, iEL4_EM2[1] + 1 - 2)
    if nextNE > Nfactor * maxEL4_EM2:
        sumE += nextNE
    nextNE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 + 6, iEL4_EM2[1] + 1 + 2)
    if nextNE > Nfactor * maxEL4_EM2:
        sumE += nextNE
    nextNE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 + 5, iEL4_EM2[1] + 1 + 2)
    if nextNE > Nfactor * maxEL4_EM2:
        sumE += nextNE
    nextNE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 + 4, iEL4_EM2[1] + 1 + 2)
    if nextNE > Nfactor * maxEL4_EM2:
        sumE += nextNE
    nextNE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 + 3, iEL4_EM2[1] + 1 + 2)
    if nextNE > Nfactor * maxEL4_EM2:
        sumE += nextNE
    nextNE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 + 2, iEL4_EM2[1] + 1 + 2)
    if nextNE > Nfactor * maxEL4_EM2:
        sumE += nextNE
    nextNE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 + 1, iEL4_EM2[1] + 1 + 2)
    if nextNE > Nfactor * maxEL4_EM2:
        sumE += nextNE
    nextNE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 + 0, iEL4_EM2[1] + 1 + 2)
    if nextNE > Nfactor * maxEL4_EM2:
        sumE += nextNE
    nextNE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 - 1, iEL4_EM2[1] + 1 + 2)
    if nextNE > Nfactor * maxEL4_EM2:
        sumE += nextNE
    nextNE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 - 2, iEL4_EM2[1] + 1 + 2)
    if nextNE > Nfactor * maxEL4_EM2:
        sumE += nextNE
    nextNE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 - 3, iEL4_EM2[1] + 1 + 2)
    if nextNE > Nfactor * maxEL4_EM2:
        sumE += nextNE
    nextNE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 - 4, iEL4_EM2[1] + 1 + 2)
    if nextNE > Nfactor * maxEL4_EM2:
        sumE += nextNE
    nextNE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 - 5, iEL4_EM2[1] + 1 + 2)
    if nextNE > Nfactor * maxEL4_EM2:
        sumE += nextNE
    nextNE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 - 5, iEL4_EM2[1] + 1 + 2)
    if nextNE > Nfactor * maxEL4_EM2:
        sumE += nextNE
    nextNE = myTOB.getEnergy(4, iEL4_EM2[0] + 4 - 7, iEL4_EM2[1] + 1 + 2)
    if nextNE > Nfactor * maxEL4_EM2:
        sumE += nextNE

    # add layer 5 = EM3
    nextE = myTOB.getEnergy(5, iEL5_EM3[0] + 1 - 1, iEL5_EM3[1] + 1 + 1)
    if nextE > factor * maxEL5_EM3:
        sumE += nextE
    nextE = myTOB.getEnergy(5, iEL5_EM3[0] + 1 - 1, iEL5_EM3[1] + 1 + 0)
    if nextE > factor * maxEL5_EM3:
        sumE += nextE
    nextE = myTOB.getEnergy(5, iEL5_EM3[0] + 1 - 1, iEL5_EM3[1] + 1 - 1)
    if nextE > factor * maxEL5_EM3:
        sumE += nextE
    nextE = myTOB.getEnergy(5, iEL5_EM3[0] + 1 + 0, iEL5_EM3[1] + 1 + 1)
    if nextE > factor * maxEL5_EM3:
        sumE += nextE
    nextE = myTOB.getEnergy(5, iEL5_EM3[0] + 1 + 0, iEL5_EM3[1] + 1 - 1)
    if nextE > factor * maxEL5_EM3:
        sumE += nextE
    nextE = myTOB.getEnergy(5, iEL5_EM3[0] + 1 + 1, iEL5_EM3[1] + 1 + 1)
    if nextE > factor * maxEL5_EM3:
        sumE += nextE
    nextE = myTOB.getEnergy(5, iEL5_EM3[0] + 1 + 1, iEL5_EM3[1] + 1 + 0)
    if nextE > factor * maxEL5_EM3:
        sumE += nextE
    nextE = myTOB.getEnergy(5, iEL5_EM3[0] + 1 + 1, iEL5_EM3[1] + 1 - 1)
    if nextE > factor * maxEL5_EM3:
        sumE += nextE

    nextNE = myTOB.getEnergy(5, iEL5_EM3[0] + 1 - 2, iEL5_EM3[1] + 1 + 1)
    if nextNE > Nfactor * maxEL5_EM3:
        sumE += nextNE
    nextNE = myTOB.getEnergy(5, iEL5_EM3[0] + 1 - 2, iEL5_EM3[1] + 1 + 2)
    if nextNE > Nfactor * maxEL5_EM3:
        sumE += nextNE
    nextNE = myTOB.getEnergy(5, iEL5_EM3[0] + 1 - 1, iEL5_EM3[1] + 1 + 2)
    if nextNE > Nfactor * maxEL5_EM3:
        sumE += nextNE
    nextNE = myTOB.getEnergy(5, iEL5_EM3[0] + 1 - 2, iEL5_EM3[1] + 1 + 0)
    if nextNE > Nfactor * maxEL5_EM3:
        sumE += nextNE
    nextNE = myTOB.getEnergy(5, iEL5_EM3[0] + 1 - 2, iEL5_EM3[1] + 1 - 1)
    if nextNE > Nfactor * maxEL5_EM3:
        sumE += nextNE
    nextNE = myTOB.getEnergy(5, iEL5_EM3[0] + 1 - 2, iEL5_EM3[1] + 1 - 2)
    if nextNE > Nfactor * maxEL5_EM3:
        sumE += nextNE
    nextNE = myTOB.getEnergy(5, iEL5_EM3[0] + 1 - 1, iEL5_EM3[1] + 1 - 2)
    if nextNE > Nfactor * maxEL5_EM3:
        sumE += nextNE
    nextNE = myTOB.getEnergy(5, iEL5_EM3[0] + 1 + 0, iEL5_EM3[1] + 1 + 2)
    if nextNE > Nfactor * maxEL5_EM3:
        sumE += nextNE
    nextNE = myTOB.getEnergy(5, iEL5_EM3[0] + 1 + 0, iEL5_EM3[1] + 1 - 2)
    if nextNE > Nfactor * maxEL5_EM3:
        sumE += nextNE
    nextNE = myTOB.getEnergy(5, iEL5_EM3[0] + 1 + 1, iEL5_EM3[1] + 1 + 2)
    if nextNE > Nfactor * maxEL5_EM3:
        sumE += nextNE
    nextNE = myTOB.getEnergy(5, iEL5_EM3[0] + 1 + 2, iEL5_EM3[1] + 1 + 2)
    if nextNE > Nfactor * maxEL5_EM3:
        sumE += nextNE
    nextNE = myTOB.getEnergy(5, iEL5_EM3[0] + 1 + 2, iEL5_EM3[1] + 1 + 1)
    if nextNE > Nfactor * maxEL5_EM3:
        sumE += nextNE
    nextNE = myTOB.getEnergy(5, iEL5_EM3[0] + 1 + 2, iEL5_EM3[1] + 1 + 0)
    if nextNE > Nfactor * maxEL5_EM3:
        sumE += nextNE
    nextNE = myTOB.getEnergy(5, iEL5_EM3[0] + 1 + 2, iEL5_EM3[1] + 1 - 1)
    if nextNE > Nfactor * maxEL5_EM3:
        sumE += nextNE
    nextNE = myTOB.getEnergy(5, iEL5_EM3[0] + 1 + 2, iEL5_EM3[1] + 1 - 2)
    if nextNE > Nfactor * maxEL5_EM3:
        sumE += nextNE
    nextNE = myTOB.getEnergy(5, iEL5_EM3[0] + 1 + 1, iEL5_EM3[1] + 1 - 2)
    if nextNE > Nfactor * maxEL5_EM3:
        sumE += nextNE

    # add layer 6 = HAD
    nextE = myTOB.getEnergy(6, iEL6_HAD[0] + 1 - 1, iEL6_HAD[1] + 1 + 1)
    if nextE > factor * maxEL6_HAD:
        sumE += nextE
    nextE = myTOB.getEnergy(6, iEL6_HAD[0] + 1 - 1, iEL6_HAD[1] + 1 + 0)
    if nextE > factor * maxEL6_HAD:
        sumE += nextE
    nextE = myTOB.getEnergy(6, iEL6_HAD[0] + 1 - 1, iEL6_HAD[1] + 1 - 1)
    if nextE > factor * maxEL6_HAD:
        sumE += nextE
    nextE = myTOB.getEnergy(6, iEL6_HAD[0] + 1 + 0, iEL6_HAD[1] + 1 + 1)
    if nextE > factor * maxEL6_HAD:
        sumE += nextE
    nextE = myTOB.getEnergy(6, iEL6_HAD[0] + 1 + 0, iEL6_HAD[1] + 1 - 1)
    if nextE > factor * maxEL6_HAD:
        sumE += nextE
    nextE = myTOB.getEnergy(6, iEL6_HAD[0] + 1 + 1, iEL6_HAD[1] + 1 + 1)
    if nextE > factor * maxEL6_HAD:
        sumE += nextE
    nextE = myTOB.getEnergy(6, iEL6_HAD[0] + 1 + 1, iEL6_HAD[1] + 1 + 0)
    if nextE > factor * maxEL6_HAD:
        sumE += nextE
    nextE = myTOB.getEnergy(6, iEL6_HAD[0] + 1 + 1, iEL6_HAD[1] + 1 - 1)
    if nextE > factor * maxEL6_HAD:
        sumE += nextE

    nextNE = myTOB.getEnergy(6, iEL6_HAD[0] + 1 - 2, iEL6_HAD[1] + 1 + 1)
    if nextNE > Nfactor * maxEL6_HAD:
        sumE += nextNE
    nextNE = myTOB.getEnergy(6, iEL6_HAD[0] + 1 - 2, iEL6_HAD[1] + 1 + 2)
    if nextNE > Nfactor * maxEL6_HAD:
        sumE += nextNE
    nextNE = myTOB.getEnergy(6, iEL6_HAD[0] + 1 - 1, iEL6_HAD[1] + 1 + 2)
    if nextNE > Nfactor * maxEL6_HAD:
        sumE += nextNE
    nextNE = myTOB.getEnergy(6, iEL6_HAD[0] + 1 - 2, iEL6_HAD[1] + 1 + 0)
    if nextNE > Nfactor * maxEL6_HAD:
        sumE += nextNE
    nextNE = myTOB.getEnergy(6, iEL6_HAD[0] + 1 - 2, iEL6_HAD[1] + 1 - 1)
    if nextNE > Nfactor * maxEL6_HAD:
        sumE += nextNE
    nextNE = myTOB.getEnergy(6, iEL6_HAD[0] + 1 - 2, iEL6_HAD[1] + 1 - 2)
    if nextNE > Nfactor * maxEL6_HAD:
        sumE += nextNE
    nextNE = myTOB.getEnergy(6, iEL6_HAD[0] + 1 - 1, iEL6_HAD[1] + 1 - 2)
    if nextNE > Nfactor * maxEL6_HAD:
        sumE += nextNE
    nextNE = myTOB.getEnergy(6, iEL6_HAD[0] + 1 + 0, iEL6_HAD[1] + 1 + 2)
    if nextNE > Nfactor * maxEL6_HAD:
        sumE += nextNE
    nextNE = myTOB.getEnergy(6, iEL6_HAD[0] + 1 + 0, iEL6_HAD[1] + 1 - 2)
    if nextNE > Nfactor * maxEL6_HAD:
        sumE += nextNE
    nextNE = myTOB.getEnergy(6, iEL6_HAD[0] + 1 + 1, iEL6_HAD[1] + 1 + 2)
    if nextNE > Nfactor * maxEL6_HAD:
        sumE += nextNE
    nextNE = myTOB.getEnergy(6, iEL6_HAD[0] + 1 + 2, iEL6_HAD[1] + 1 + 2)
    if nextNE > Nfactor * maxEL6_HAD:
        sumE += nextNE
    nextNE = myTOB.getEnergy(6, iEL6_HAD[0] + 1 + 2, iEL6_HAD[1] + 1 + 1)
    if nextNE > Nfactor * maxEL6_HAD:
        sumE += nextNE
    nextNE = myTOB.getEnergy(6, iEL6_HAD[0] + 1 + 2, iEL6_HAD[1] + 1 + 0)
    if nextNE > Nfactor * maxEL6_HAD:
        sumE += nextNE
    nextNE = myTOB.getEnergy(6, iEL6_HAD[0] + 1 + 2, iEL6_HAD[1] + 1 - 1)
    if nextNE > Nfactor * maxEL6_HAD:
        sumE += nextNE
    nextNE = myTOB.getEnergy(6, iEL6_HAD[0] + 1 + 2, iEL6_HAD[1] + 1 - 2)
    if nextNE > Nfactor * maxEL6_HAD:
        sumE += nextNE
    nextNE = myTOB.getEnergy(6, iEL6_HAD[0] + 1 + 1, iEL6_HAD[1] + 1 - 2)
    if nextNE > Nfactor * maxEL6_HAD:
        sumE += nextNE

    return sumE


class tauROI:
    def __init__(self):
        self.TLV = ROOT.TLorentzVector()
        self.Iso = 0

    def setP4(self, Pt, Eta, Phi):
        self.TLV.SetPtEtaPhiM(Pt, Eta, Phi, 0)

    def Pt(self):
        return self.TLV.Pt()

    def Eta(self):
        return self.TLV.Eta()

    def Phi(self):
        return self.TLV.Phi()

    def setIso(self, iso):
        self.Iso = iso


# Histogram Array Holder Class
class histHolder:

    # Resolution Plots Dictionary
    resPlots = {}
    # Turn-On Curves Dictionary
    toPlots = {}
    # Turn-On Eta Curves Dictionary
    toEtaPlots = {}
    # nTOB Dictionary
    ntobPlots = {}

    def __init__(self, listOfCriteria):
        for myCrit in listOfCriteria:
            self.toPlots[myCrit] = ROOT.TProfile(
                "TurnOn" + myCrit, "TurnOn" + myCrit, 400, 0, 200
            )
            self.toPlots[myCrit].SetXTitle("true #tau vis. p_{T}")
            self.toPlots[myCrit].SetYTitle("Efficiency")

            self.toEtaPlots[myCrit] = ROOT.TProfile(
                "TurnOnEta" + myCrit, "TurnOnEta" + myCrit, 100, 1, 2
            )
            self.toEtaPlots[myCrit].SetXTitle("true #tau #eta")
            self.toEtaPlots[myCrit].SetYTitle("Efficiency")

            self.resPlots[myCrit] = ROOT.TH1D(
                "Reso" + myCrit, "Reso" + myCrit, 50, -1.0, 1.0
            )
            self.resPlots[myCrit].SetXTitle(
                "#frac{L1 #tau p_{T} - true #tau p_{T}}{true #tau p_{T}}"
            )
            self.resPlots[myCrit].SetYTitle("Entries")

            self.ntobPlots[myCrit] = ROOT.TH1D(
                "nTOB" + myCrit, "nTOB" + myCrit, 400, 0, 200
            )
            self.ntobPlots[myCrit].SetXTitle("TOB p_{T}")
            self.ntobPlots[myCrit].SetYTitle("Entries")

    def fillTO(self, name, pT, status):
        self.toPlots[name].Fill(pT / 1000.0, status)

    def fillTOEta(self, name, eta, status):
        self.toEtaPlots[name].Fill(abs(eta), status)

    def fillRes(self, name, etTrue, etCand):
        if etTrue > 0:
            self.resPlots[name].Fill((etCand - etTrue) / etTrue)

    def fillTOB(self, name, etCand):
        self.ntobPlots[name].Fill(etCand)

    def savePlots(self):
        for key in self.resPlots:
            self.resPlots[key].Write()
        for key in self.toPlots:
            self.toPlots[key].Write()
        for key in self.toEtaPlots:
            self.toEtaPlots[key].Write()
        for key in self.ntobPlots:
            self.ntobPlots[key].Write()
