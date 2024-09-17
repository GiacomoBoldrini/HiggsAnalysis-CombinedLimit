#ifndef CMSHistErrorPropagator_h
#define CMSHistErrorPropagator_h
#include <ostream>
#include <vector>
#include <memory>
#include "RooAbsReal.h"
#include "RooArgSet.h"
#include "RooListProxy.h"
#include "RooRealProxy.h"
#include "Rtypes.h"
#include "TH1F.h"
#include "TH3.h"
#include "FastTemplate_Old.h"
#include "SimpleCacheSentry.h"
#include "CMSHistFunc.h"
#include "CMSHistV.h"

class CMSHistErrorPropagator : public RooAbsReal {
private:
  struct BarlowBeeston {
    bool init = false;
    std::vector<unsigned> use;
    std::vector<double> dat;
    std::vector<double> valsum;
    std::vector<double> toterr;
    std::vector<double> err;
    std::vector<double> b;
    std::vector<double> c;
    std::vector<double> tmp;
    std::vector<double> x1;
    std::vector<double> x2;
    std::vector<double> res;
    std::vector<double> gobs;
    std::set<RooAbsArg*> dirty_prop;
    std::vector<RooRealVar*> push_res;
  };
public:
  CMSHistErrorPropagator();

  CMSHistErrorPropagator(const char* name, const char* title, RooRealVar& x,
                         RooArgList const& funcs, RooArgList const& coeffs);

  CMSHistErrorPropagator(const char* name, const char* title, RooRealVar& x,
                         RooArgList const& funcs, RooArgList const& coeffs, TH3 const& hist);

  CMSHistErrorPropagator(CMSHistErrorPropagator const& other, const char* name = 0);

  TObject* clone(const char* newname) const override {
    return new CMSHistErrorPropagator(*this, newname);
  }

  ~CMSHistErrorPropagator() override {;}

  void applyErrorShifts(unsigned idx, FastHisto const& nominal, FastHisto & result);

  Double_t evaluate() const override;

  RooArgList * setupBinPars(double poissonThreshold);

  std::unique_ptr<RooArgSet> getSentryArgs() const;

  void printMultiline(std::ostream& os, Int_t contents, Bool_t verbose,
                      TString indent) const override;

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars,
                              const char* rangeName = 0) const override;

  Double_t analyticalIntegral(Int_t code, const char* rangeName = 0) const override;

  void setData(RooAbsData const& data) const;

  void setAnalyticBarlowBeeston(bool flag) const;

  inline FastHisto const& cache() const { return cache_; }

  RooArgList wrapperList() const;
  RooArgList const& coefList() const { return coeffs_; }
  RooArgList const& funcList() const { return funcs_; }
  FastHisto3D const correlationHist() const { return mccorr_; }
  RooListProxy const& binPars() const { return binpars_; }
  std::map<unsigned int, std::string> const& correlationSamples() const { return corrsamples_; }
  std::vector<std::vector<unsigned>> const& binTypes() const { return bintypes_; } // this should be filled after t2w
  
  std::map<std::string, Double_t> getProcessNorms() const;

  friend class CMSHistV<CMSHistErrorPropagator>;

 protected:
  RooRealProxy x_;
  RooListProxy funcs_;
  RooListProxy coeffs_;
  // store the RooRealVar corresponding to autoMCstat parameters
  // for BarlowBeeston
  RooListProxy binpars_;
  mutable FastHisto3D mccorr_; // Need to make this better... why with TH3 it does not work?
  // save the mapping of bin labels to bin number.
  mutable std::map<unsigned int, std::string> corrsamples_;
  //
  mutable std::vector<CMSHistFunc const*> vfuncs_; //!
  mutable std::vector<RooAbsReal const*> vcoeffs_; //!

  // vbinpars_ has the same structure of bintypes_ below but instead of saving 
  // unsigned integers it stores the actual RooAbsReal that are created.
  // It maintains the lengths but retrieves the parameters only if
  // bintypes_[j][i] >= 1 && bintypes_[j][i] < 4
  // so ignoring bins where the total error is <= 0 and processes in bins under poisson that
  // have negative or zero errors or negative yield. 
  mutable std::vector<std::vector<RooAbsReal *>> vbinpars_; //!


  // it's a vector storing numbers for each bin of the distribution 
  // in this datacard bin. By default 
  // {{0}, {0}, ... {0}} 

  // two cases: below or abve poisson threshold.

  // if we are below poisson threshold for bin j
  // then bintypes_[j] is a vector of integers with size equal to the number
  // of processes in this bin. Each process is assigned a value of 4 by default
  // if only bin 1 is under th and you have 3 processes it will look like
  // {{0}, {4,4,4}, {0}, ...{0}}
  // if the bintype=4 for a bin-process, this will be skipped in the computation (either
  // yield < 0 or error <= 0). 
  
  // If yield and error > 0 and yield>error then we compute the effective number 
  // of MC event for this sample. If it is below poisson th, a poisson par is associated
  // with this sample in this bin and bintype=2. Otherwise a gaussian parameter is associated with 
  // this sample in this bin bintype=3.
  // If yield and error > 0 but yield<error, a gaussian parameter is associated with 
  // this sample in this bin bintype=3.

  // if we are above poisson threshold and the total error (sum in quadrature of the errors of 
  // each process contributing in this bin) is greater than zero, bintypes for this bin has only 
  // one entry and set to bintypes_[j][0] = 1. The total unc. in this bin will be scaled with 
  // a gaussian constrained parameter. The latter parameter is assigned an attribute
  // var->setAttribute("forBarlowBeeston");


  // Summary:

  // bintypes_[j] -> j index of template
  // len(bintypes_[j]) == 1 -> we are above poisson threshold or there is only one process in this bin
  //    bintypes_[j][0] = 1 -> toterror > 0, one gauss param
  //    bintypes_[j][0] = 0 -> toterror <= 0, nothing happens
  
  // len(bintypes_[j]) > 1 -> we are below poisson threshold -> i sample index
  //    bintypes_[j][i] = 4 -> ignored, error<=0 or yield<0
  //    bintypes_[j][i] = 3 -> sample yield>error>0 but sample below th. -> pois param on sample and bin
  //    bintypes_[j][i] = 2 -> sample yield>error>0 sample above th. -> gaus param on sample and bin
  //    bintypes_[j][i] = 2 -> sample error>yield>0 -> gaus param on sample and bin

  std::vector<std::vector<unsigned>> bintypes_;

  mutable std::vector<double> coeffvals_; //!
  mutable FastHisto valsum_; //!
  mutable FastHisto cache_; //!
  mutable std::vector<double> err2sum_; //!
  mutable std::vector<double> err2sumcorr_; //!
  mutable std::vector<double> toterr_; //!
  mutable std::vector<std::vector<double>> binmods_; //!
  mutable std::vector<std::vector<double>> scaledbinmods_; //!
  mutable SimpleCacheSentry sentry_; //!
  mutable SimpleCacheSentry binsentry_; //!
  mutable std::vector<double> data_; //!

  mutable BarlowBeeston bb_; //!

  mutable bool initialized_; //! not to be serialized

  mutable int last_eval_; //! not to be serialized

  mutable bool analytic_bb_; //! not to be serialized

  void initialize() const;
  void updateCache(int eval = 1) const;

  void runBarlowBeeston() const;


 private:
  ClassDefOverride(CMSHistErrorPropagator,1)
};

#endif
