import ROOT

f = ROOT.TFile("model.root")
w = f.Get("w")

error = [i for i in w.allFunctions() if i.GetName() == "prop_bininclusive_all"][0]

print(error.binTypes())
print(error.correlationSamples())

print(error.binPars())

error.correlationHist().Dump()
