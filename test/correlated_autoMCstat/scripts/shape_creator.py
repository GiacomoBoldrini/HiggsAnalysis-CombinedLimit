import ROOT
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
import uproot
import awkward as ak
from copy import deepcopy
import vector
vector.register_awkward()
import sys
from hist import Hist
import hist
import matplotlib as mpl
from tqdm import tqdm
plt.style.use(hep.style.CMS)  # or ATLAS/LHCb2
mpl.use("Agg")
plt.style.use(hep.style.CMS)  # or ATLAS/LHCb2
mpl.use("Agg")

def Fill3DHisto(events, var_='mll', nbins_=10, range_=[800,1000]):
    
    # the three histos we need for morphin into sm lin quad
    #h_w0 = ROOT.TH1F("h_w0", "h_w0", nbins_, *range_)
    #h_w1 = ROOT.TH1F("h_w1", "h_w1", nbins_, *range_)
    #h_wm1 = ROOT.TH1F("h_wm1", "h_wm1", nbins_, *range_)
    
    td = np.zeros((nbins_, 3, 3))
    edges = [
             np.linspace(0, nbins_, nbins_+1),
             np.linspace(0, 3, 3+1),
             np.linspace(0, 3, 3+1),
            ]
    
    # single histos for dc and root file
    
    h_w0 = Hist(
        hist.axis.Regular(nbins_,*range_, name=var_, label=var_, underflow=False, overflow=False), 
        storage=hist.storage.Weight()
    )
    
    h_w1 = Hist(
        hist.axis.Regular(nbins_,*range_, name=var_, label=var_, underflow=False, overflow=False), 
        storage=hist.storage.Weight()
    )
    
    h_wm1 = Hist(
        hist.axis.Regular(nbins_,*range_, name=var_, label=var_, underflow=False, overflow=False), 
        storage=hist.storage.Weight()
    )

    h_w0.fill(events[var_], weight=events.baseW*events.LHEReweightingWeight[:,0])
    h_w1.fill(events[var_], weight=events.baseW*events.LHEReweightingWeight[:,40])
    h_wm1.fill(events[var_], weight=events.baseW*events.LHEReweightingWeight[:,39])
    
    
    
    # twod with just correlations for samples "A", "B", "C"
    
    #     A     B       C
    #   ---------------------
    # A | 1 | rhoAB | rhoAC |
    #   |   |       |       |
    # B | - |   -   | rhoBC |
    #   |   |       |       |
    # C | - |   -   |   -   |
    #   ---------------------
    
        
    binedges_ = np.linspace(*range_, nbins_+1)
    
    # need to compute the MC stat unc for each bin
    # so we need to cut the events in the bin window in order
    # to recover the events falling in that bin and compute retrieve the weights
    
    for b_ in range(nbins_):
        
        # filter the events falling in this bin
        bin_mask = (events[var_]>=binedges_[b_]) & (events[var_]<binedges_[b_+1])
        ev__ = events[bin_mask]
        
        #retrieve weights for morphing 
        # for events falling in this bin
        
        #cQl1
        w0 = ev__.baseW*ev__.LHEReweightingWeight[:,0]
        w1 = ev__.baseW*ev__.LHEReweightingWeight[:,22]
        wm1 = ev__.baseW*ev__.LHEReweightingWeight[:,21]
        
        s_AB = ak.corr(w0, w1)
        s_AC = ak.corr(w0, wm1)
        s_BC = ak.corr(w1, wm1)

        
        # this should appear in top left
        td[b_,0,0] = 1.0
        td[b_,0,1] = s_AB
        td[b_,0,2] = s_AC
        td[b_,1,2] = s_BC
        
    
        
    return (td, edges), (h_w0, h_w1, h_wm1)
        
        


def createLeptons(events):
    
    ele_mask = (abs(events.LHEPart.pdgId) == 11) & (events.LHEPart.status == 1)
    mu_mask = (abs(events.LHEPart.pdgId) == 13) & (events.LHEPart.status == 1)
    
    
    Leptons = events.LHEPart[ele_mask + mu_mask]
    Electrons = events.LHEPart[ele_mask]
    Muons = events.LHEPart[mu_mask]

    
    Leptons = Leptons[ak.argsort(Leptons.pt, ascending=False, axis=-1)]
    Electrons = Electrons[ak.argsort(Electrons.pt, ascending=False, axis=-1)]
    Muons = Muons[ak.argsort(Muons.pt, ascending=False, axis=-1)]
    
    events["Leptons"] = Leptons
    events["Electrons"] = Electrons
    events["Muons"] = Muons
    
    events[("Leptons", "mass")] = ak.zeros_like(events.Leptons.pt)
    events[("Electrons", "mass")] = ak.zeros_like(events.Electrons.pt)
    events[("Muons", "mass")] = ak.zeros_like(events.Muons.pt)
    
    
    
    return events
    

def read_array(tree, branch_name, start, stop):
    interp = tree[branch_name].interpretation
    interp._forth = True
    return tree[branch_name].array(
        interp,
        entry_start=start,
        entry_stop=stop,
        decompression_executor=uproot.source.futures.TrivialExecutor(),
        interpretation_executor=uproot.source.futures.TrivialExecutor(),
    )

def read_events(filename, start=0, stop=1e9, read_form={}):
    #print("start reading")
    f = uproot.open(filename, num_workers=2)
    tree = f["Events"]
    start = min(start, tree.num_entries)
    stop = min(stop, tree.num_entries)
    if start >= stop:
        return ak.Array([])

    branches = [k.name for k in tree.branches]

    events = {}
    form = deepcopy(read_form)

    for coll in form:
        d = {}
        coll_branches = form[coll].pop("branches")

        if len(coll_branches) == 0:
            if coll in branches:
                events[coll] = read_array(tree, coll, start, stop)
            continue

        for branch in coll_branches:
            branch_name = coll + "_" + branch
            if branch_name in branches:
                d[branch] = read_array(tree, branch_name, start, stop)

        if len(d.keys()) == 0:
            print("did not find anything for", coll, file=sys.stderr)
            continue

        events[coll] = ak.zip(d, **form[coll])
        del d

    # f.close()
    #print("created events")
    _events = ak.zip(events, depth_limit=1)
    del events
    return _events


read_form = {
  "LHEPart": {
    "branches": [
      "pt",
      "eta",
      "phi",
      "mass",
      "pdgId",
      "status"
    ],
    "with_name": "Momentum4D"
  },
  "LHEReweightingWeight": {
      "branches": [],
      "with_name": None
  },
  "LHEPdfWeight": {
      "branches": [],
      "with_name": None
  },
  "LHEScaleWeight": {
      "branches": [],
      "with_name": None
  },
  "baseW": {
      "branches": [],
      "with_name": None
  }
}

directories = ["mll_1000_1500", "mll_100_200", "mll_1500_inf", "mll_200_400", "mll_400_600", "mll_50_100", "mll_600_800", "mll_800_1000"]

events = read_events("/eos/user/g/gboldrin/Zee_dim6_LHE/mll_binned/rootfiles/mll_800_1000/nAOD_LHE_0.root", read_form=read_form)
events = createLeptons(events)

events = events[ak.num(events.Muons) > 0]
events = events[ak.all(abs(events.Muons.eta) < 2.5, axis=1)]
events = events[events.Muons.pt[:,0] > 25]

from glob import glob

for dir__ in directories:
    # take only the first 10 files for each mll bin
    ls = glob(f"/eos/user/g/gboldrin/Zee_dim6_LHE/mll_binned/rootfiles/{dir__}/*")[:3]
    for file in tqdm(ls):
        if file.split("/")[-1] == "nAOD_LHE_0.root" and dir__ == "mll_800_1000": continue
        ev__ = read_events(file, read_form=read_form)
        ev__ = createLeptons(ev__)

        ev__ = ev__[ak.num(ev__.Muons) > 0]
        ev__ = ev__[ak.all(abs(ev__.Muons.eta) < 2.5, axis=1)]
        ev__ = ev__[ev__.Muons.pt[:,0] > 25]

        try:
            events = ak.concatenate((events, ev__))
        except:
            print("ERROR for file {}".format(file))
            continue
            
        #print("Merge successfull")

# apply selections 

# restrict to muons

events = events[ak.num(events.Muons) > 0]
events = events[ak.all(abs(events.Muons.eta) < 2.5, axis=1)]
events = events[events.Muons.pt[:,0] > 25]

# define mll and ptll

events['mll'] = (events.Muons[:,0] + events.Muons[:,1]).mass
events['ptl1'] = (events.Muons[:,0]).pt

print("Number of events {}".format(len(events)))


# make datacards and shapes for the cQl1 operator

mll_bins = (0, 2000)
bins = 10

tdhisto, onedhistos = Fill3DHisto(events, nbins_=bins, range_= mll_bins)

print("uproot")
fname = "shapes.root"
f = uproot.recreate(f"{fname}")
f["histo_correlation"] = tdhisto
f["histo_sm"] = onedhistos[0]
f["histo_w1_cQl1"] = onedhistos[1]
f["histo_wm1_cQl1"] = onedhistos[2]
f.close()

print(tdhisto[0][0])



# write datacard

f_o = open("datacard.txt", "w")

f_o.write("## Shape input card\n")
f_o.write("imax 1 number of channels\n")
f_o.write("jmax * number of background\n")
f_o.write("kmax * number of nuisance parameters\n")
f_o.write("----------------------------------------------------------------------------------------------------\n")
f_o.write("bin         inclusive_all\n")
f_o.write(f"shapes  *           * shapes/{fname}     histo_$PROCESS histo_$PROCESS_$SYSTEMATIC\n")
f_o.write("bin                                                                             inclusive_all                 inclusive_all                 inclusive_all \n")                
f_o.write("process                                                                         w1_cQl1                            wm1_cQl1              sm   \n")                          
f_o.write("process                                                                         0                             -1                            -2    \n")                          
f_o.write("rate                                                                            {:.4f}                        {:.4f}                        {:.4f}\n".format(onedhistos[1].sum().value, onedhistos[2].sum().value, onedhistos[0].sum().value))          
f_o.write("----------------------------------------------------------------------------------------------------   \n") 
f_o.write("lumi_13TeV_2016                                             lnN                 1.01                          1.01                          1.01     \n")                         
f_o.write("----------------------------------------------------------------------------------------------------   \n") 
f_o.write("* autoMCStats 10 0 1   \n")

f_o.close()

