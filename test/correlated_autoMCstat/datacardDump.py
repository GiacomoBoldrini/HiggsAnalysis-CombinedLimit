#!/usr/bin/env python3
from __future__ import absolute_import

import re
from optparse import OptionParser
from sys import argv, exit, modules, stderr, stdout
import sys
import ROOT
from HiggsAnalysis.CombinedLimit.DatacardParser import *
from HiggsAnalysis.CombinedLimit.ModelTools import *
from HiggsAnalysis.CombinedLimit.PhysicsModel import *
from HiggsAnalysis.CombinedLimit.ShapeTools import *

# import ROOT with a fix to get batch mode (http://root.cern.ch/phpBB3/viewtopic.php?t=3198)
argv.append("-b-")

ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True
argv.remove("-b-")

parser = OptionParser()
addDatacardParserOptions(parser)


(options, args) = parser.parse_args()
options.stat = False
options.bin = True  # fake that is a binary output, so that we parse shape lines
options.out = "tmp.root"
options.fileName = args[0]
options.cexpr = False
options.fixpars = False
options.libs = []
options.verbose = 3
options.poisson = 0
options.nuisancesToExclude = []
options.noJMax = True
options.justCheckPhysicsModel = False
options.nuisanceFunctions = []
options.nuisanceGroupFunctions = []
options.allowNoBackground = True

# sys.argv = ["-b-"]

# ROOT.gROOT.SetBatch(True)
# ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")


file = open(args[0], "r")
DC = parseCard(file, options)
if not DC.hasShapes:
    DC.hasShapes = True

## Load tools to build workspace
MB = None
if DC.hasShapes:
    MB = ShapeBuilder(DC, options)
else:
    MB = CountingModelBuilder(DC, options)

print(DC.systs)
sys.exit(0)

# Let's initialize the ShapeBuilder routine
# With AnalyticalAnomalousCoupling

options.physModel = "HiggsAnalysis.AnalyticAnomalousCoupling.AnomalousCouplingEFTNegative:analiticAnomalousCouplingEFTNegative"
(physModMod, physModName) = options.physModel.split(":")

__import__(physModMod)
mod = sys.modules[physModMod]
physics = getattr(mod, physModName)

# Specify to build the model for only one operator

options.physOpt = ["eftOperators=cHj3"]
physics.setPhysicsOptions(options.physOpt)

# Now build the model for real

MB.setPhysics(physics)
#MB.doModel(justCheckPhysicsModel=options.justCheckPhysicsModel)

MB.doObservables()
MB.physics.doParametersOfInterest()
MB.doNuisances()
MB.doExtArgs()
MB.doRateParams()
MB.doAutoFlatNuisancePriors()
MB.doFillNuisPdfsAndSets()
MB.doExpectedEvents()
MB.doIndividualModels()


# in this datacard we have only one datacard bin 
datacard_bin = "inclusive_all"

#processes contributing to this bin and yield
# {'quad_cHj3': 2383726.5165, 'sm_lin_quad_cHj3': 530589045.2754, 'sm': 466878199.8117}
procs = MB.DC.exp[datacard_bin]

"""
for p in procs.keys():
    (pdf, coeff) = (
        MB.getPdf(datacard_bin, p),
        MB.out.function("n_exp_bin%s_proc_%s" % (datacard_bin, p)),
    )

    print("--> pdf")
    print(pdf)

    print("--> coeff")
    print(coeff)
    print(coeff.getAttribute("combine.signal"))
    print(coeff.getStringAttribute("combine.process"))
    print(coeff.getStringAttribute("combine.channel"))

    print("------")
"""

# get all objects

print(MB.objstore)

"""
for b in DC.bins:
    print(" ============= ", b, "====================")
    if options.channel != None and (options.channel != b):
        continue
    exps = {}
    for p, e in DC.exp[b].items():  # so that we get only self.DC.processes contributing to this bin
        exps[p] = [e, []]
        if exps[p][0] < 0:
            s0 = MB.getShape(b, p)
            if s0.InheritsFrom("TH1"):
                exps[p][0] = s0.Integral()
    for lsyst, nofloat, pdf, pdfargs, errline in DC.systs:
        if pdf in ("param", "flatParam"):
            continue
        # begin skip systematics
        skipme = False
        for xs in options.excludeSyst:
            if re.search(xs, lsyst):
                skipme = True
        if skipme:
            continue
        # end skip systematics
        for p in DC.exp[b].keys():  # so that we get only self.DC.processes contributing to this bin
            if errline[b][p] == 0:
                continue
            if pdf == "gmN":
                exps[p][1].append(1 / sqrt(pdfargs[0] + 1))
            elif pdf == "gmM":
                exps[p][1].append(errline[b][p])
            elif type(errline[b][p]) == list:
                kmax = max(
                    errline[b][p][0],
                    errline[b][p][1],
                    1.0 / errline[b][p][0],
                    1.0 / errline[b][p][1],
                )
                exps[p][1].append(kmax - 1.0)
            elif pdf == "lnN":
                exps[p][1].append(max(errline[b][p], 1.0 / errline[b][p]) - 1.0)
            elif ("shape" in pdf) and not options.norm:
                s0 = MB.getShape(b, p)
                sUp = MB.getShape(b, p, lsyst + "Up")
                sDown = MB.getShape(b, p, lsyst + "Down")
                if s0.InheritsFrom("TH1"):
                    ratios = [
                        sUp.Integral() / s0.Integral(),
                        sDown.Integral() / s0.Integral(),
                    ]
                    ratios += [1 / ratios[0], 1 / ratios[1]]
                    exps[p][1].append(max(ratios) - 1)
    procs = list(DC.exp[b].keys())
    procs.sort()
    fmt = ("%%-%ds " % max([len(p) for p in procs])) + "  " + options.format
    for p in procs:
        relunc = sqrt(sum([x * x for x in exps[p][1]]))
        print(fmt % (p, exps[p][0], exps[p][0] * relunc))
"""
