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

def plot(events, scale=1):



    mll_bins = (0, 2000)
    ptl_bins = (100, 600)
    bins = 10



    # fill the axes
    h_mll = Hist(
        hist.axis.Regular(
            bins, *mll_bins, name="mll", label="mll [GeV]", underflow=False, overflow=False
        ),
        hist.storage.Weight()
    )

    h_ptl = Hist(
        hist.axis.Regular(
            bins, *ptl_bins,  name="ptl1", label="ptl1 [GeV]", underflow=False, overflow=False
        ),
        hist.storage.Weight()
    )

    ##########

    h_mll_cje = Hist(
        hist.axis.Regular(
            bins, *mll_bins, name="mll_cje", label="mll [GeV]", underflow=False, overflow=False
        ),
        hist.storage.Weight()
    )

    h_ptl_cje = Hist(
        hist.axis.Regular(
            bins, *ptl_bins,  name="ptl1_cje", label="ptl1 [GeV]", underflow=False, overflow=False
        ),
        hist.storage.Weight()
    )

    ##########

    h_mll_mcje = Hist(
        hist.axis.Regular(
            bins, *mll_bins, name="mll_mcje", label="mll [GeV]", underflow=False, overflow=False
        ),
        hist.storage.Weight()
    )

    h_ptl_mcje = Hist(
        hist.axis.Regular(
            bins, *ptl_bins,  name="ptl1_mcje", label="ptl1 [GeV]", underflow=False, overflow=False
        ),
        hist.storage.Weight()
    )



    h_mll.fill(events.mll, weight=events.LHEReweightingWeight[:,0]*scale)
    h_ptl.fill(events.ptl1, weight=events.LHEReweightingWeight[:,0]*scale)

    h_mll_cje.fill(events.mll, weight=events.LHEReweightingWeight[:,40]*scale)
    h_ptl_cje.fill(events.ptl1, weight=events.LHEReweightingWeight[:,40]*scale)

    h_mll_mcje.fill(events.mll, weight=events.LHEReweightingWeight[:,39]*scale)
    h_ptl_mcje.fill(events.ptl1, weight=events.LHEReweightingWeight[:,39]*scale)

    h_mll.view().value = h_mll.view().value/scale
    h_ptl.view().value = h_ptl.view().value/scale
    h_mll_cje.view().value = h_mll_cje.view().value/scale
    h_ptl_cje.view().value = h_ptl_cje.view().value/scale
    h_mll_mcje.view().value = h_mll_mcje.view().value/scale
    h_ptl_mcje.view().value = h_ptl_mcje.view().value/scale


    fig, ax = plt.subplots(1,2, figsize=(20,10), dpi=160)
    # h_mll.plot(ax=ax[0], w2method="sqrt", label="cje=0 (SM)")

    h_mll_cje.plot(ax=ax[0], w2method="sqrt", color='green', label="cje=1")
    #h_mll_mcje.plot(ax=ax[0], w2method="sqrt", color='black', label="cje=-1")

    ax[0].legend()
    hep.cms.label(
            "Preliminary", data=True, ax=ax[0]
        )  # ,fontsize=16)

    # h_ptl.plot(ax=ax[1], w2method="sqrt", label="cje=0 (SM)")
    h_ptl_cje.plot(ax=ax[1], w2method="sqrt", color='green', label="cje=1")
    #h_ptl_mcje.plot(ax=ax[1], w2method="sqrt", color='black', label="cje=-1")

    ax[1].legend()

    hep.cms.label(
            "Preliminary", data=True, ax=ax[1]
        )  # ,fontsize=16)

    ax[0].set_ylabel("Events")
    ax[1].set_ylabel("Events")

    ax[0].set_yscale("log")
    ax[1].set_yscale("log")
    fig.tight_layout()

    fig.savefig(f"./plots/mll_ptll_{scale}.pdf")
    fig.savefig(f"./plots/mll_ptll_{scale}.png")

    return


def Fill3DHisto(events, var_='mll', nbins_=10, range_=[800,1000], scale=1):
    
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
    
    samples = ["sm", "w1_cje", "wm1_cje", "w1_cHj1", "wm1_cHj1", "w11_cje_cHj1"]

    h = Hist(
        hist.axis.Regular(nbins_, 0, 2000),
        hist.axis.StrCategory(samples),
        hist.axis.StrCategory(samples),
    )
    
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

    
    h_w1_cHj1 = Hist(
        hist.axis.Regular(nbins_,*range_, name=var_, label=var_, underflow=False, overflow=False), 
        storage=hist.storage.Weight()
    )
    
    h_wm1_cHj1 = Hist(
        hist.axis.Regular(nbins_,*range_, name=var_, label=var_, underflow=False, overflow=False), 
        storage=hist.storage.Weight()
    )

    h_w11_cje_cHj1 = Hist(
        hist.axis.Regular(nbins_,*range_, name=var_, label=var_, underflow=False, overflow=False), 
        storage=hist.storage.Weight()
    )

    h_w0.fill(events[var_], weight=events.LHEReweightingWeight[:,0]*scale)

    h_w1.fill(events[var_], weight=events.LHEReweightingWeight[:,40]*scale)
    h_wm1.fill(events[var_], weight=events.LHEReweightingWeight[:,39]*scale)

    h_w1_cHj1.fill(events[var_], weight=events.LHEReweightingWeight[:,6]*scale)
    h_wm1_cHj1.fill(events[var_], weight=events.LHEReweightingWeight[:,5]*scale)

    h_w11_cje_cHj1.fill(events[var_], weight=events.LHEReweightingWeight[:,118]*scale)

    # normalize back the histos
    h_w0.view().value = h_w0.view().value / scale
    h_w1.view().value = h_w1.view().value / scale
    h_wm1.view().value = h_wm1.view().value / scale

    h_w1_cHj1.view().value = h_w1_cHj1.view().value / scale
    h_wm1_cHj1.view().value = h_wm1_cHj1.view().value / scale

    h_w11_cje_cHj1.view().value = h_w11_cje_cHj1.view().value / scale

    
    
    
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
        
        w0 = ev__.LHEReweightingWeight[:,0] * scale

        # this is cje
        w1 = ev__.LHEReweightingWeight[:,40] * scale
        wm1 = ev__.LHEReweightingWeight[:,39] * scale

        # this is cHj1
        w1_cHj1 = ev__.LHEReweightingWeight[:,6] * scale
        wm1_cHj1 = ev__.LHEReweightingWeight[:,5] * scale

        w11_cje_cHj1 = ev__.LHEReweightingWeight[:,118] * scale


        
        s_AA = ak.sum(w0**2)
        s_BB = ak.sum(w1**2)
        s_CC = ak.sum(wm1**2)

        s_DD = ak.sum(w1_cHj1**2)
        s_EE = ak.sum(wm1_cHj1**2)

        s_FF = ak.sum(w11_cje_cHj1**2)
        

        # This is the old correlation with np.corr 
        # Which is wrong because we cannot define a correlation 
        # on the bin yield from the event weights
        
        # s_AB = 2*ak.corr(w0, w1)*np.sqrt(s_AA)*np.sqrt(s_BB)
        # s_AC = 2*ak.corr(w0, wm1)*np.sqrt(s_AA)*np.sqrt(s_CC)
        # s_AD = 2*ak.corr(w0, w1_cHj1)*np.sqrt(s_AA)*np.sqrt(s_DD)
        # s_AE = 2*ak.corr(w0, wm1_cHj1)*np.sqrt(s_AA)*np.sqrt(s_EE)
        # s_AF = 2*ak.corr(w0, w11_cje_cHj1)*np.sqrt(s_AA)*np.sqrt(s_FF)
 
        # s_BC = 2*ak.corr(w1, wm1)*np.sqrt(s_BB)*np.sqrt(s_CC)
        # s_BD = 2*ak.corr(w1, w1_cHj1)*np.sqrt(s_BB)*np.sqrt(s_DD)
        # s_BE = 2*ak.corr(w1, wm1_cHj1)*np.sqrt(s_BB)*np.sqrt(s_EE)
        # s_BF = 2*ak.corr(w1, w11_cje_cHj1)*np.sqrt(s_BB)*np.sqrt(s_FF)
 
        # s_CD = 2*ak.corr(wm1, w1_cHj1)*np.sqrt(s_CC)*np.sqrt(s_DD)
        # s_CE = 2*ak.corr(wm1, wm1_cHj1)*np.sqrt(s_CC)*np.sqrt(s_EE)
        # s_CF = 2*ak.corr(wm1, w11_cje_cHj1)*np.sqrt(s_CC)*np.sqrt(s_FF)
 
        # s_DE = 2*ak.corr(w1_cHj1, wm1_cHj1)*np.sqrt(s_DD)*np.sqrt(s_EE)
        # s_DF = 2*ak.corr(w1_cHj1, w11_cje_cHj1)*np.sqrt(s_DD)*np.sqrt(s_FF)
 
        # s_EF = 2*ak.corr(wm1_cHj1, w11_cje_cHj1)*np.sqrt(s_EE)*np.sqrt(s_FF)

        s_AB = 2*ak.sum(w0*w1)
        s_AC = 2*ak.sum(w0*wm1)
        s_AD = 2*ak.sum(w0*w1_cHj1)
        s_AE = 2*ak.sum(w0*wm1_cHj1)
        s_AF = 2*ak.sum(w0*w11_cje_cHj1)
        s_BC = 2*ak.sum(w1*wm1)
        s_BD = 2*ak.sum(w1*w1_cHj1)
        s_BE = 2*ak.sum(w1*wm1_cHj1)
        s_BF = 2*ak.sum(w1*w11_cje_cHj1)
        s_CD = 2*ak.sum(wm1*w1_cHj1)
        s_CE = 2*ak.sum(wm1*wm1_cHj1)
        s_CF = 2*ak.sum(wm1*w11_cje_cHj1)
        s_DE = 2*ak.sum(w1_cHj1*wm1_cHj1)
        s_DF = 2*ak.sum(w1_cHj1*w11_cje_cHj1)
        s_EF = 2*ak.sum(wm1_cHj1*w11_cje_cHj1)


        
        # this should appear in top left
        h[b_,0,0] = s_AA
        h[b_,0,1] = s_AB
        h[b_,0,2] = s_AC
        h[b_,0,3] = s_AD
        h[b_,0,4] = s_AE
        h[b_,0,5] = s_AF

        h[b_,1,1] = s_BB
        h[b_,1,2] = s_BC
        h[b_,1,3] = s_BD
        h[b_,1,4] = s_BE
        h[b_,1,5] = s_BF

        h[b_,2,2] = s_CC
        h[b_,2,3] = s_CD
        h[b_,2,4] = s_CE
        h[b_,2,5] = s_CF

        h[b_,3,3] = s_DD
        h[b_,3,4] = s_DE
        h[b_,3,5] = s_DF

        h[b_,4,4] = s_EE
        h[b_,4,5] = s_EF

        h[b_,5,5] = s_FF
        
    
        
    return h, (h_w0, h_w1, h_wm1, h_w1_cHj1, h_wm1_cHj1, h_w11_cje_cHj1)
        
        


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
    ls = glob(f"/eos/user/g/gboldrin/Zee_dim6_LHE/mll_binned/rootfiles/{dir__}/*")[:1]
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


for scale in [1, 2, 5, 10]:
    plot(events, scale=scale)
    # make datacards and shapes for the cje operator

    mll_bins = (0, 2000)
    bins = 10

    tdhisto, onedhistos = Fill3DHisto(events, nbins_=bins, range_= mll_bins, scale=scale)

    fname = f"shapes/shapes_{scale}.root"
    f = uproot.recreate(f"{fname}")
    f["histo_correlation"] = tdhisto
    f["histo_sm"] = onedhistos[0]
    f["histo_w1_cje"] = onedhistos[1]
    f["histo_wm1_cje"] = onedhistos[2]
    f["histo_w1_cHj1"] = onedhistos[3]
    f["histo_wm1_cHj1"] = onedhistos[4]
    f["histo_w11_cje_cHj1"] = onedhistos[5]
    f.close()



    # write datacard

    f_o = open(f"datacard_{scale}.txt", "w")

    f_o.write("## Shape input card\n")
    f_o.write("imax 1 number of channels\n")
    f_o.write("jmax * number of background\n")
    f_o.write("kmax * number of nuisance parameters\n")
    f_o.write("----------------------------------------------------------------------------------------------------\n")
    f_o.write("bin         inclusive_all\n")
    f_o.write(f"shapes  *           * {fname}     histo_$PROCESS histo_$PROCESS_$SYSTEMATIC\n")
    f_o.write("bin                                                                             inclusive_all                 inclusive_all                 inclusive_all                 inclusive_all                 inclusive_all                 inclusive_all \n")                
    f_o.write("process                                                                         w1_cje                       wm1_cje                      sm                            w1_cHj1                       wm1_cHj1                      w11_cje_cHj1\n")                          
    f_o.write("process                                                                         0                             -1                            -2                            -3                            -4                            -5\n")                          
    f_o.write("rate                                                                            {:.4f}                        {:.4f}                        {:.4f}                        {:.4f}                        {:.4f}                        {:.4f}\n".format(onedhistos[1].sum().value, onedhistos[2].sum().value, onedhistos[0].sum().value, onedhistos[3].sum().value, onedhistos[4].sum().value, onedhistos[5].sum().value))          
    f_o.write("----------------------------------------------------------------------------------------------------   \n") 
    f_o.write("lumi_13TeV_2016                                             lnN                 1.01                          1.01                          1.01                          1.01                          1.01                          1.01\n")                         
    f_o.write("----------------------------------------------------------------------------------------------------   \n") 
    f_o.write("* autoMCStats 0 1 1   \n")
    f_o.write(f"* autoMCCorr {fname}   histo_correlation \n")

    f_o.close()




# compute EFTneg shapes

tdhisto, onedhistos = Fill3DHisto(events, nbins_=bins, range_= mll_bins, scale=1)

sm = onedhistos[0].copy()
sm_l_q_cje = onedhistos[1].copy()

q_cje = onedhistos[1].copy()

values_cje = (onedhistos[1].view().value + onedhistos[2].view().value - 2*onedhistos[0].view().value)*0.5
q_cje.view().value = values_cje


sm_l_q_cHj1 = onedhistos[3].copy()

q_cHj1 = onedhistos[3].copy()

values_cHj1 = (onedhistos[3].view().value + onedhistos[4].view().value - 2*onedhistos[0].view().value)*0.5
q_cHj1.view().value = values_cHj1

sm_l_q_mixed = onedhistos[5].copy()

# set tp zero variances
q_cje.view().variance = np.zeros(len(q_cje.view().variance))
q_cHj1.view().variance = np.zeros(len(q_cHj1.view().variance))
sm.view().variance = np.zeros(len(sm.view().variance))
sm_l_q_cje.view().variance = np.zeros(len(sm_l_q_cje.view().variance))
sm_l_q_cHj1.view().variance = np.zeros(len(sm_l_q_cHj1.view().variance))
sm_l_q_mixed.view().variance = np.zeros(len(sm_l_q_mixed.view().variance))



fname = "shapes/shapes_EFTNeg.root"
f = uproot.recreate(f"{fname}")
f["histo_sm"] = sm
f["histo_sm_lin_quad_cje"] = sm_l_q_cje
f["histo_quad_cje"] = q_cje
f["histo_quad_cHj1"] = q_cHj1
f["histo_sm_lin_quad_cHj1"] = sm_l_q_cHj1
f["histo_sm_lin_quad_mixed_cje_cHj1"] = sm_l_q_mixed
f.close()


# write datacard

f_o = open("datacard_EFTNeg.txt", "w")

f_o.write("## Shape input card\n")
f_o.write("imax 1 number of channels\n")
f_o.write("jmax * number of background\n")
f_o.write("kmax * number of nuisance parameters\n")
f_o.write("----------------------------------------------------------------------------------------------------\n")
f_o.write("bin         inclusive_all\n")
f_o.write(f"shapes  *           * shapes/{fname}     histo_$PROCESS histo_$PROCESS_$SYSTEMATIC\n")
f_o.write("bin                                                                             inclusive_all                 inclusive_all                 inclusive_all                 inclusive_all                 inclusive_all                 inclusive_all \n")                
f_o.write("process                                                                         sm_lin_quad_cje              quad_cje                     sm                            sm_lin_quad_cHj1              quad_cHj1                     sm_lin_quad_mixed_cje_cHj1\n")                          
f_o.write("process                                                                         0                             -1                            -2                            -3                            -4                            -5 \n")                          
f_o.write("rate                                                                            {:.4f}                        {:.4f}                        {:.4f}                        {:.4f}                        {:.4f}                        {:.4f}\n".format(sm_l_q_cje.sum().value, q_cje.sum().value, sm.sum().value, sm_l_q_cHj1.sum().value, q_cHj1.sum().value, sm_l_q_mixed.sum().value))          
f_o.write("----------------------------------------------------------------------------------------------------   \n") 
f_o.write("lumi_13TeV_2016                lnN                                           1.01                          1.01                          1.01                          1.01                          1.01                          1.01     \n")                         
f_o.write("----------------------------------------------------------------------------------------------------   \n") 
f_o.write("## * autoMCStats 0 1 1   \n")
f_o.write(f"## * autoMCCorr shapes/{fname}   histo_correlation \n")

f_o.close()