import os
import sys
from pathlib import Path

import ROOT
from ROOT import *
import copy
import math

# from math import *
from array import array

from constants import *
from sagi_util import *
from argparse import ArgumentParser
from csv import DictWriter


######################################
######## Initialization
######################################

parser = ArgumentParser()
parser.add_argument("-o", "--output", type=str, help="Output directory")
parser.add_argument("-c", "--count", nargs="?", type=int, const=COUNT_BREAK, help="Count break value")
parser.add_argument("-e", "--events", nargs="?", type=int, const=MAX_EVENTS, help="Max events")
args = parser.parse_args()
output_dir = args.output
if output_dir is not None:
    output_dir = Path(output_dir)

max_count = args.count
max_events = args.events

efex_data = []
jfex_data = []


print("hello world")


myFile = ROOT.TFile("/afs/cern.ch/work/b/barak/public/L1CALO/phase1/tobTree.root")
myTree = myFile.Get("tobTree")

from glob import glob

fileList = glob(
    "/afs/cern.ch/work/b/barak/public/L1CALO/phase1/user.viveiros.JZ0W.Timing_v04_OUTPUT/user.viveiros.21160992.OUTPUT._000003.root"
)
# fileList = glob("/afs/cern.ch/work/b/barak/public/L1CALO/phase1/user.viveiros.JZ0W.PeakFinder_v04_OUTPUT/*")#user.viveiros.21160991.OUTPUT._000003.root")
print(fileList)

myTree = ROOT.TChain("ntuple")
for myFile in fileList:
    myTree.AddFile(myFile)


calibFile = ROOT.TFile("calibInput.root")
calibProfile = calibFile.Get("px")

# Define list of criteria to be studied, and create plot array
myCrits = ["TAU12", "TAU12R3", "TAU12jR3", "TAU20", "TAU20R3", "TAU20jR3"]

# Total tau candidate energy fractions
eFrac0 = ROOT.TH1D("eFrac0", "eFrac0", 100, 0, 1)
eFrac1 = ROOT.TH1D("eFrac1", "eFrac1", 100, 0, 1)
eFrac2 = ROOT.TH1D("eFrac2", "eFrac2", 100, 0, 1)
eFrac3 = ROOT.TH1D("eFrac3", "eFrac3", 100, 0, 1)
eFracH = ROOT.TH1D("eFracH", "eFracH", 100, 0, 1)

# Cells needed for energy containment
eDist0 = ROOT.TH1D("eDist0", "eDist0", 9, -0.5, 9.5)
eDist1 = ROOT.TH1D("eDist1", "eDist1", 36, -0.5, 36.5)
eDist2 = ROOT.TH1D("eDist2", "eDist2", 36, -0.5, 36.5)
eDist3 = ROOT.TH1D("eDist3", "eDist3", 9, -0.5, 9.5)
eDistH = ROOT.TH1D("eDistH", "eDistH", 9, -0.5, 9.5)

# Look at the energy correlations
ETvsEM2 = ROOT.TH2D("ETvsEM2", "ETvsEM2", 50, 0, 1.0, 50, -1, 1)
ETvsEM3 = ROOT.TH2D("ETvsEM3", "ETvsEM3", 50, 0, 0.5, 50, -1, 1)
ETvsEW2 = ROOT.TH2D("ETvsEW2", "ETvsEW2", 50, 0, 1.0, 50, -1, 1)

ETvsDepth = ROOT.TH2D("ETvsDepth", "ETvsDepth", 50, 0, 500, 50, -1, 1)
ETvsDens = ROOT.TH2D("ETvsDens", "ETvsDens", 50, -5, 0, 50, -1, 1)
ETvsDens2 = ROOT.TH2D("ETvsDens2", "ETvsDens2", 50, -5, 0, 50, -1, 1)
ETvsDensC = ROOT.TH2D("ETvsDensC", "ETvsDensC", 50, -8.5, -5, 50, -1, 1)
ETvsWidth = ROOT.TH2D("ETvsWidth", "ETvsWidth", 50, 0, 0.15, 50, -1, 1)
ETvsET = ROOT.TH2D("ETvsET", "ETvsET", 50, 10, 100, 50, -1, 1)
ETvsETC = ROOT.TH2D("ETvsETC", "ETvsETC", 50, 10, 100, 50, -1, 1)

ETvsEta = ROOT.TH2D("ETvsEta", "ETvsEta", 50, 0, 2.5, 50, -1, 1)

# Look at lowest pT objects
highestCell = ROOT.TH1D("highestCell", "highestCell", 50, 0, 3125)
fullRoi = ROOT.TH1D("fullRoi", "fullRoi", 50, 0, 10000)

# Iso Xin
Oregon_ratio = ROOT.TH1D("Oregon_ratio", "Oregon_ratio", 40, -5, 5)

# Look at potential veto variables # isolation variables (Run3_L1Tau_Ore_iso_12pass, Run3_L1Tau_Ore_iso_20pass)
coreRatioB1 = ROOT.TH1D("coreRatioB1", "coreRatioB1", 100, -5000, 5000)  # Run3_jTau_iso
coreRatioS1 = ROOT.TH1D("coreRatioS1", "coreRatioS1", 100, -5000, 5000)  # Run3_jTau_iso
coreRatioB2 = ROOT.TH1D("coreRatioB2", "coreRatioB2", 50, 0, 1.0)  # Run3_L1Tau_Ore_iso
coreRatioS2 = ROOT.TH1D("coreRatioS2", "coreRatioS2", 50, 0, 1.0)  # Run3_L1Tau_Ore_iso
coreRatioB3 = ROOT.TH1D("coreRatioB3", "coreRatioB3", 50, 0, 1.0)  # Run3_L1Tau_iso
coreRatioS3 = ROOT.TH1D("coreRatioS3", "coreRatioS3", 50, 0, 1.0)
coreRatioB4 = ROOT.TH1D("coreRatioB4", "coreRatioB4", 50, 0, 1.0)
coreRatioS4 = ROOT.TH1D("coreRatioS4", "coreRatioS4", 50, 0, 1.0)
coreRatioB5 = ROOT.TH1D("coreRatioB5", "coreRatioB5", 50, 0, 1.0)
coreRatioS5 = ROOT.TH1D("coreRatioS5", "coreRatioS5", 50, 0, 1.0)

myDist = ROOT.TH1D("myDist", "myDist", 30, -1.2, 0.5)
isoRatio = ROOT.TH1D("isoRatio", "isoRatio", 100, 0, 1.5)
fisoRatio = ROOT.TH1D("fisoRatio", "fisoRatio", 100, 0, 1.5)

ejDR = ROOT.TH1D("ejDR", "ejDR", 50, 0, 5)
ejPt = ROOT.TH1D("ejPt", "ejPt", 20, -1, 1)
jPt = ROOT.TH1D("jPt", "jPt", 100, 0, 100)
eiso = ROOT.TH1D("eiso", "eiso", 50, 0, 1.0)
jiso = ROOT.TH1D("jiso", "jiso", 3000, -150, 150)
eiso12 = ROOT.TH1D("eiso12", "eiso12", 50, 0, 1.0)
# jiso12 = ROOT.TH1D("jiso12", "jiso12", 40, -2, 2)
ejiso12 = ROOT.TH1D("ejiso12", "ejiso12", 40, -2, 2)
eiso15 = ROOT.TH1D("eiso15", "eiso15", 50, 0, 1.0)
# jiso15 = ROOT.TH1D("jiso15", "jiso15", 40, -2, 2)
ejiso15 = ROOT.TH1D("ejiso15", "ejiso15", 40, -2, 2)
jiso10 = ROOT.TH1D("jiso10", "jiso10", 80, -2, 2)
jiso11 = ROOT.TH1D("jiso11", "jiso11", 80, -2, 2)
jiso12 = ROOT.TH1D("jiso12", "jiso12", 80, -2, 2)
jiso13 = ROOT.TH1D("jiso13", "jiso13", 80, -2, 2)
jiso14 = ROOT.TH1D("jiso14", "jiso14", 80, -2, 2)
jiso15 = ROOT.TH1D("jiso15", "jiso15", 80, -2, 2)
jiso16 = ROOT.TH1D("jiso16", "jiso16", 80, -2, 2)
jiso17 = ROOT.TH1D("jiso17", "jiso17", 80, -2, 2)
jiso18 = ROOT.TH1D("jiso18", "jiso18", 80, -2, 2)
jiso19 = ROOT.TH1D("jiso19", "jiso19", 80, -2, 2)
jiso20 = ROOT.TH1D("jiso20", "jiso20", 80, -2, 2)
jiso21 = ROOT.TH1D("jiso21", "jiso21", 80, -2, 2)
jiso22 = ROOT.TH1D("jiso22", "jiso22", 80, -2, 2)
jiso23 = ROOT.TH1D("jiso23", "jiso23", 80, -2, 2)
jiso24 = ROOT.TH1D("jiso24", "jiso24", 80, -2, 2)
jiso25 = ROOT.TH1D("jiso25", "jiso25", 80, -2, 2)

# Counters for histogram normalization
nE0 = 0
nE1 = 0
nE2 = 0
nE3 = 0
nEH = 0
nCands = 0
nCands_failOr = 0
nCands_failbigClus = 0
nCands_passOfailbigClus = 0

######################################
######## Event Loop
######################################
nDR = nFull = 0

count = 0
for evenum, myEvt in enumerate(myTree, start=1):

    print(f"evenum={evenum}, count={count}")
    if max_events is not None and evenum >= max_events:
        break
    if max_count is not None and count >= max_count:
        break
    # Construct a list of truth and reco
    truthTaus = []
    for i in range(myEvt.TruthTaus_ptvis.size()):
        tempVector = ROOT.TLorentzVector()
        tempVector.SetPtEtaPhiM(
            myEvt.TruthTaus_ptvis.at(i),
            myEvt.TruthTaus_etavis.at(i),
            myEvt.TruthTaus_phivis.at(i),
            0,
        )
        truthTaus += [tempVector]

    # Construct a list of trigger candidates for Run-3
    run3Taus = []
    for i in range(myEvt.Run3_L1Tau_Et.size()):
        myROI = tauROI()
        myROI.setP4(
            myEvt.Run3_L1Tau_Ore_Et.at(i),
            myEvt.Run3_L1Tau_eta.at(i),
            myEvt.Run3_L1Tau_phi.at(i),
        )
        myROI.setIso(myEvt.Run3_L1Tau_Ore_iso.at(i))
        #        print "Run3:", i, myEvt.Run3_L1Tau_Ore_iso.at(i), myROI.Iso
        run3Taus += [myROI]

    # Construct a list of trigger candidates for Run-3 jFEX
    run3jTausBase = []
    run3jTaus = []
    for i in range(myEvt.Run3_jTau_Et.size()):
        myROI = tauROI()
        myROI.setP4(
            myEvt.Run3_jTau_Et.at(i),
            myEvt.Run3_jTau_eta.at(i),
            myEvt.Run3_jTau_phi.at(i),
        )
        myROI.setIso(myEvt.Run3_jTau_iso.at(i))
        #        print "Run3j:", i, myEvt.Run3_jTau_iso.at(i), myROI.Iso
        run3jTausBase += [myROI]

    # Construct a list of trigger candidates for Run-3 jFEX - making sure those are the ones the fPGA will get
    run3jTaus1 = (
        run3jTaus2
    ) = (
        run3jTaus3
    ) = (
        run3jTaus4
    ) = (
        run3jTaus5
    ) = (
        run3jTaus6
    ) = (
        run3jTaus7
    ) = (
        run3jTaus8
    ) = (
        run3jTaus9
    ) = (
        run3jTaus10
    ) = (
        run3jTaus11
    ) = (
        run3jTaus12
    ) = (
        run3jTaus13
    ) = (
        run3jTaus14
    ) = (
        run3jTaus15
    ) = (
        run3jTaus16
    ) = (
        run3jTaus17
    ) = (
        run3jTaus18
    ) = (
        run3jTaus19
    ) = (
        run3jTaus20
    ) = (
        run3jTaus21
    ) = (
        run3jTaus22
    ) = run3jTaus23 = run3jTaus24 = []  # ROOT.std.vector("TLorentzVector")

    for jTau in run3jTausBase:
        if (jTau.Eta() < 4.9 and jTau.Eta() >= 1.6) and (
            ROOT.TVector2.Phi_0_2pi(jTau.Phi()) < 2 * math.pi
            and ROOT.TVector2.Phi_0_2pi(jTau.Phi()) >= 1.5 * math.pi
        ):
            run3jTaus1 += [(jTau, jTau.Pt())]
        elif (jTau.Eta() < 4.9 and jTau.Eta() >= 1.6) and (
            ROOT.TVector2.Phi_0_2pi(jTau.Phi()) < 1.5 * math.pi
            and ROOT.TVector2.Phi_0_2pi(jTau.Phi()) >= 1 * math.pi
        ):
            run3jTaus2 += [(jTau, jTau.Pt())]
        elif (jTau.Eta() < 4.9 and jTau.Eta() >= 1.6) and (
            ROOT.TVector2.Phi_0_2pi(jTau.Phi()) < 1 * math.pi
            and ROOT.TVector2.Phi_0_2pi(jTau.Phi()) >= 0.5 * math.pi
        ):
            run3jTaus3 += [(jTau, jTau.Pt())]
        elif (jTau.Eta() < 4.9 and jTau.Eta() >= 1.6) and (
            ROOT.TVector2.Phi_0_2pi(jTau.Phi()) < 0.5 * math.pi
            and ROOT.TVector2.Phi_0_2pi(jTau.Phi()) >= 0 * math.pi
        ):
            run3jTaus4 += [(jTau, jTau.Pt())]
        elif (jTau.Eta() < 1.6 and jTau.Eta() >= 0.8) and (
            ROOT.TVector2.Phi_0_2pi(jTau.Phi()) < 2 * math.pi
            and ROOT.TVector2.Phi_0_2pi(jTau.Phi()) >= 1.5 * math.pi
        ):
            run3jTaus5 += [(jTau, jTau.Pt())]
        elif (jTau.Eta() < 1.6 and jTau.Eta() >= 0.8) and (
            ROOT.TVector2.Phi_0_2pi(jTau.Phi()) < 1.5 * math.pi
            and ROOT.TVector2.Phi_0_2pi(jTau.Phi()) >= 1 * math.pi
        ):
            run3jTaus6 += [(jTau, jTau.Pt())]
        elif (jTau.Eta() < 1.6 and jTau.Eta() >= 0.8) and (
            ROOT.TVector2.Phi_0_2pi(jTau.Phi()) < 1 * math.pi
            and ROOT.TVector2.Phi_0_2pi(jTau.Phi()) >= 0.5 * math.pi
        ):
            run3jTaus7 += [(jTau, jTau.Pt())]
        elif (jTau.Eta() < 1.6 and jTau.Eta() >= 0.8) and (
            ROOT.TVector2.Phi_0_2pi(jTau.Phi()) < 0.5 * math.pi
            and ROOT.TVector2.Phi_0_2pi(jTau.Phi()) >= 0 * math.pi
        ):
            run3jTaus8 += [(jTau, jTau.Pt())]
        elif (jTau.Eta() < 0.8 and jTau.Eta() >= 0) and (
            ROOT.TVector2.Phi_0_2pi(jTau.Phi()) < 2 * math.pi
            and ROOT.TVector2.Phi_0_2pi(jTau.Phi()) >= 1.5 * math.pi
        ):
            run3jTaus9 += [(jTau, jTau.Pt())]
        elif (jTau.Eta() < 0.8 and jTau.Eta() >= 0) and (
            ROOT.TVector2.Phi_0_2pi(jTau.Phi()) < 1.5 * math.pi
            and ROOT.TVector2.Phi_0_2pi(jTau.Phi()) >= 1 * math.pi
        ):
            run3jTaus10 += [(jTau, jTau.Pt())]
        elif (jTau.Eta() < 0.8 and jTau.Eta() >= 0) and (
            ROOT.TVector2.Phi_0_2pi(jTau.Phi()) < 1 * math.pi
            and ROOT.TVector2.Phi_0_2pi(jTau.Phi()) >= 0.5 * math.pi
        ):
            run3jTaus11 += [(jTau, jTau.Pt())]
        elif (jTau.Eta() < 0.8 and jTau.Eta() >= 0) and (
            ROOT.TVector2.Phi_0_2pi(jTau.Phi()) < 0.5 * math.pi
            and ROOT.TVector2.Phi_0_2pi(jTau.Phi()) >= 0 * math.pi
        ):
            run3jTaus12 += [(jTau, jTau.Pt())]
        elif (jTau.Eta() < 0 and jTau.Eta() >= -0.8) and (
            ROOT.TVector2.Phi_0_2pi(jTau.Phi()) < 2 * math.pi
            and ROOT.TVector2.Phi_0_2pi(jTau.Phi()) >= 1.5 * math.pi
        ):
            run3jTaus13 += [(jTau, jTau.Pt())]
        elif (jTau.Eta() < 0 and jTau.Eta() >= -0.8) and (
            ROOT.TVector2.Phi_0_2pi(jTau.Phi()) < 1.5 * math.pi
            and ROOT.TVector2.Phi_0_2pi(jTau.Phi()) >= 1 * math.pi
        ):
            run3jTaus14 += [(jTau, jTau.Pt())]
        elif (jTau.Eta() < 0 and jTau.Eta() >= -0.8) and (
            ROOT.TVector2.Phi_0_2pi(jTau.Phi()) < 1 * math.pi
            and ROOT.TVector2.Phi_0_2pi(jTau.Phi()) >= 0.5 * math.pi
        ):
            run3jTaus15 += [(jTau, jTau.Pt())]
        elif (jTau.Eta() < 0 and jTau.Eta() >= -0.8) and (
            ROOT.TVector2.Phi_0_2pi(jTau.Phi()) < 0.5 * math.pi
            and ROOT.TVector2.Phi_0_2pi(jTau.Phi()) >= 0 * math.pi
        ):
            run3jTaus16 += [(jTau, jTau.Pt())]
        elif (jTau.Eta() < -0.8 and jTau.Eta() >= -1.6) and (
            ROOT.TVector2.Phi_0_2pi(jTau.Phi()) < 2 * math.pi
            and ROOT.TVector2.Phi_0_2pi(jTau.Phi()) >= 1.5 * math.pi
        ):
            run3jTaus17 += [(jTau, jTau.Pt())]
        elif (jTau.Eta() < -0.8 and jTau.Eta() >= -1.6) and (
            ROOT.TVector2.Phi_0_2pi(jTau.Phi()) < 1.5 * math.pi
            and ROOT.TVector2.Phi_0_2pi(jTau.Phi()) >= 1 * math.pi
        ):
            run3jTaus18 += [(jTau, jTau.Pt())]
        elif (jTau.Eta() < -0.8 and jTau.Eta() >= -1.6) and (
            ROOT.TVector2.Phi_0_2pi(jTau.Phi()) < 1 * math.pi
            and ROOT.TVector2.Phi_0_2pi(jTau.Phi()) >= 0.5 * math.pi
        ):
            run3jTaus19 += [(jTau, jTau.Pt())]
        elif (jTau.Eta() < -0.8 and jTau.Eta() >= -1.6) and (
            ROOT.TVector2.Phi_0_2pi(jTau.Phi()) < 0.5 * math.pi
            and ROOT.TVector2.Phi_0_2pi(jTau.Phi()) >= 0 * math.pi
        ):
            run3jTaus20 += [(jTau, jTau.Pt())]
        elif (jTau.Eta() < -1.6 and jTau.Eta() >= -4.9) and (
            ROOT.TVector2.Phi_0_2pi(jTau.Phi()) < 2 * math.pi
            and ROOT.TVector2.Phi_0_2pi(jTau.Phi()) >= 1.5 * math.pi
        ):
            run3jTaus21 += [(jTau, jTau.Pt())]
        elif (jTau.Eta() < -1.6 and jTau.Eta() >= -4.9) and (
            ROOT.TVector2.Phi_0_2pi(jTau.Phi()) < 1.5 * math.pi
            and ROOT.TVector2.Phi_0_2pi(jTau.Phi()) >= 1 * math.pi
        ):
            run3jTaus22 += [(jTau, jTau.Pt())]
        elif (jTau.Eta() < -1.6 and jTau.Eta() >= -4.9) and (
            ROOT.TVector2.Phi_0_2pi(jTau.Phi()) < 1 * math.pi
            and ROOT.TVector2.Phi_0_2pi(jTau.Phi()) >= 0.5 * math.pi
        ):
            run3jTaus23 += [(jTau, jTau.Pt())]
        elif (jTau.Eta() < -1.6 and jTau.Eta() >= -4.9) and (
            ROOT.TVector2.Phi_0_2pi(jTau.Phi()) < 0.5 * math.pi
            and ROOT.TVector2.Phi_0_2pi(jTau.Phi()) >= 0 * math.pi
        ):
            run3jTaus24 += [(jTau, jTau.Pt())]

    # ROOT.std.sort(run3jTaus1.begin(), run3jTaus1.end(), ComparePt)
    # sort list with key
    run3jTaus1.sort(key=takeSecond, reverse=True)
    run3jTaus2.sort(key=takeSecond, reverse=True)
    run3jTaus3.sort(key=takeSecond, reverse=True)
    run3jTaus4.sort(key=takeSecond, reverse=True)
    run3jTaus5.sort(key=takeSecond, reverse=True)
    run3jTaus6.sort(key=takeSecond, reverse=True)
    run3jTaus7.sort(key=takeSecond, reverse=True)
    run3jTaus8.sort(key=takeSecond, reverse=True)
    run3jTaus9.sort(key=takeSecond, reverse=True)
    run3jTaus10.sort(key=takeSecond, reverse=True)
    run3jTaus11.sort(key=takeSecond, reverse=True)
    run3jTaus12.sort(key=takeSecond, reverse=True)
    run3jTaus13.sort(key=takeSecond, reverse=True)
    run3jTaus14.sort(key=takeSecond, reverse=True)
    run3jTaus15.sort(key=takeSecond, reverse=True)
    run3jTaus16.sort(key=takeSecond, reverse=True)
    run3jTaus17.sort(key=takeSecond, reverse=True)
    run3jTaus18.sort(key=takeSecond, reverse=True)
    run3jTaus19.sort(key=takeSecond, reverse=True)
    run3jTaus20.sort(key=takeSecond, reverse=True)
    run3jTaus21.sort(key=takeSecond, reverse=True)
    run3jTaus22.sort(key=takeSecond, reverse=True)
    run3jTaus23.sort(key=takeSecond, reverse=True)
    run3jTaus24.sort(key=takeSecond, reverse=True)
    run3jTaus += (
        run3jTaus1[0][0],
        run3jTaus1[1][0],
        run3jTaus1[2][0],
        run3jTaus1[3][0],
        run3jTaus1[4][0],
        run3jTaus1[5][0],
    )
    run3jTaus += (
        run3jTaus2[0][0],
        run3jTaus2[1][0],
        run3jTaus2[2][0],
        run3jTaus2[3][0],
        run3jTaus2[4][0],
        run3jTaus2[5][0],
    )
    run3jTaus += (
        run3jTaus3[0][0],
        run3jTaus3[1][0],
        run3jTaus3[2][0],
        run3jTaus3[3][0],
        run3jTaus3[4][0],
        run3jTaus3[5][0],
    )
    run3jTaus += (
        run3jTaus4[0][0],
        run3jTaus4[1][0],
        run3jTaus4[2][0],
        run3jTaus4[3][0],
        run3jTaus4[4][0],
        run3jTaus4[5][0],
    )
    run3jTaus += (
        run3jTaus5[0][0],
        run3jTaus5[1][0],
        run3jTaus5[2][0],
        run3jTaus5[3][0],
        run3jTaus5[4][0],
        run3jTaus5[5][0],
    )
    run3jTaus += (
        run3jTaus6[0][0],
        run3jTaus6[1][0],
        run3jTaus6[2][0],
        run3jTaus6[3][0],
        run3jTaus6[4][0],
        run3jTaus6[5][0],
    )
    run3jTaus += (
        run3jTaus7[0][0],
        run3jTaus7[1][0],
        run3jTaus7[2][0],
        run3jTaus7[3][0],
        run3jTaus7[4][0],
        run3jTaus7[5][0],
    )
    run3jTaus += (
        run3jTaus8[0][0],
        run3jTaus8[1][0],
        run3jTaus8[2][0],
        run3jTaus8[3][0],
        run3jTaus8[4][0],
        run3jTaus8[5][0],
    )
    run3jTaus += (
        run3jTaus9[0][0],
        run3jTaus9[1][0],
        run3jTaus9[2][0],
        run3jTaus9[3][0],
        run3jTaus9[4][0],
        run3jTaus9[5][0],
    )
    run3jTaus += (
        run3jTaus10[0][0],
        run3jTaus10[1][0],
        run3jTaus10[2][0],
        run3jTaus10[3][0],
        run3jTaus10[4][0],
        run3jTaus10[5][0],
    )
    run3jTaus += (
        run3jTaus11[0][0],
        run3jTaus11[1][0],
        run3jTaus11[2][0],
        run3jTaus11[3][0],
        run3jTaus11[4][0],
        run3jTaus11[5][0],
    )
    run3jTaus += (
        run3jTaus12[0][0],
        run3jTaus12[1][0],
        run3jTaus12[2][0],
        run3jTaus12[3][0],
        run3jTaus12[4][0],
        run3jTaus12[5][0],
    )
    run3jTaus += (
        run3jTaus13[0][0],
        run3jTaus13[1][0],
        run3jTaus13[2][0],
        run3jTaus13[3][0],
        run3jTaus13[4][0],
        run3jTaus13[5][0],
    )
    run3jTaus += (
        run3jTaus14[0][0],
        run3jTaus14[1][0],
        run3jTaus14[2][0],
        run3jTaus14[3][0],
        run3jTaus14[4][0],
        run3jTaus14[5][0],
    )
    run3jTaus += (
        run3jTaus15[0][0],
        run3jTaus15[1][0],
        run3jTaus15[2][0],
        run3jTaus15[3][0],
        run3jTaus15[4][0],
        run3jTaus15[5][0],
    )
    run3jTaus += (
        run3jTaus16[0][0],
        run3jTaus16[1][0],
        run3jTaus16[2][0],
        run3jTaus16[3][0],
        run3jTaus16[4][0],
        run3jTaus16[5][0],
    )
    run3jTaus += (
        run3jTaus17[0][0],
        run3jTaus17[1][0],
        run3jTaus17[2][0],
        run3jTaus17[3][0],
        run3jTaus17[4][0],
        run3jTaus17[5][0],
    )
    run3jTaus += (
        run3jTaus18[0][0],
        run3jTaus18[1][0],
        run3jTaus18[2][0],
        run3jTaus18[3][0],
        run3jTaus18[4][0],
        run3jTaus18[5][0],
    )
    run3jTaus += (
        run3jTaus19[0][0],
        run3jTaus19[1][0],
        run3jTaus19[2][0],
        run3jTaus19[3][0],
        run3jTaus19[4][0],
        run3jTaus19[5][0],
    )
    run3jTaus += (
        run3jTaus20[0][0],
        run3jTaus20[1][0],
        run3jTaus20[2][0],
        run3jTaus20[3][0],
        run3jTaus20[4][0],
        run3jTaus20[5][0],
    )
    run3jTaus += (
        run3jTaus21[0][0],
        run3jTaus21[1][0],
        run3jTaus21[2][0],
        run3jTaus21[3][0],
        run3jTaus21[4][0],
        run3jTaus21[5][0],
    )
    run3jTaus += (
        run3jTaus22[0][0],
        run3jTaus22[1][0],
        run3jTaus22[2][0],
        run3jTaus22[3][0],
        run3jTaus22[4][0],
        run3jTaus22[5][0],
    )
    run3jTaus += (
        run3jTaus23[0][0],
        run3jTaus23[1][0],
        run3jTaus23[2][0],
        run3jTaus23[3][0],
        run3jTaus23[4][0],
        run3jTaus23[5][0],
    )
    run3jTaus += (
        run3jTaus24[0][0],
        run3jTaus24[1][0],
        run3jTaus24[2][0],
        run3jTaus24[3][0],
        run3jTaus24[4][0],
        run3jTaus24[5][0],
    )

    for candTau in run3Taus:
        newDR = None
        jcandTauPt = None
        jcandTauEta = None
        jcandTauPhi = None
        jcandTauIso = None
        for jcandTau in run3jTaus:
            if abs(jcandTau.Eta()) > 2.47:
                continue
            myDR = jcandTau.TLV.DeltaR(candTau.TLV)
            if newDR is None or myDR < newDR:
                newDR = myDR
                jcandTauPt = jcandTau.Pt()
                jcandTauEta = jcandTau.Eta()
                jcandTauPhi = jcandTau.Phi()
                jcandTauIso = jcandTau.Iso / jcandTauPt
        if newDR is not None:
            count += 1
            efex_data.append(
                {
                    PT: f"{candTau.Pt():.{ACCURACY}f}",
                    ETA: f"{candTau.Eta():.{ACCURACY}f}",
                    EFEX_ISO: f"{candTau.Iso:.{ACCURACY}f}"
                }
            )
            jfex_data.append(
                {
                    PT: f"{jcandTauPt:.{ACCURACY}f}",
                    ETA: f"{jcandTauEta:.{ACCURACY}f}",
                    JFEX_ISO: f"{jcandTauIso:.{ACCURACY}f}"
                }
            )

if output_dir is None:
    for i, record in enumerate(efex_data, start=1):
        print(f"{i}) {record}")
else:
    if not output_dir.exists():
        output_dir.mkdir()
    with open(output_dir / "efex.csv", mode="w", newline="") as fd:
        csv_writer = DictWriter(fd, fieldnames=[PT, ETA, EFEX_ISO])
        csv_writer.writeheader()
        csv_writer.writerows(efex_data)
    with open(output_dir / "jfex.csv", mode="w", newline="") as fd:
        csv_writer = DictWriter(fd, fieldnames=[PT, ETA, JFEX_ISO])
        csv_writer.writeheader()
        csv_writer.writerows(jfex_data)
    print(f"Wrote data to {output_dir}")
