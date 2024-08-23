#include "../../../interface/FastTemplate_Old.h"
#include <iostream>

void load() {
    gSystem->Load("libHiggsAnalysisCombinedLimit");
    // open the tfile
    TFile* f = new TFile("shapes.root");

    // load the 3d histo for correlations
    TH3F* h = (TH3F*)(f->Get("histo_correlation"));

    // load the 1d histos 
    TH1F* h_sm = (TH1F*)(f->Get("histo_sm"));
    TH1F* h_w1 = (TH1F*)(f->Get("histo_w1_cQl1"));
    TH1F* h_wm1 = (TH1F*)(f->Get("histo_wm1_cQl1"));

    FastHisto3D _cacheNominal;

    _cacheNominal = FastHisto3D(dynamic_cast<TH3F&>(*h), false);
    _cacheNominal.Dump();
    
    return;
}