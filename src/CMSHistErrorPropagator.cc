#include "../interface/CMSHistErrorPropagator.h"
#include "../interface/CMSHistFuncWrapper.h"
#include <stdexcept>
#include <vector>
#include <ostream>
#include <memory>
#include "Math/ProbFuncMathCore.h"
#include "Math/QuantFuncMathCore.h"
#include "RooRealProxy.h"
#include "RooArgSet.h"
#include "RooAbsReal.h"
#include "RooAbsData.h"
#include "RooConstVar.h"
#include "RooGaussian.h"
#include "RooProduct.h"
#include "vectorized.h"

#define HFVERBOSE 5

CMSHistErrorPropagator::CMSHistErrorPropagator() : initialized_(false) {}

CMSHistErrorPropagator::CMSHistErrorPropagator(const char* name,
                                               const char* title,
                                               RooRealVar& x,
                                               RooArgList const& funcs,
                                               RooArgList const& coeffs)
    : RooAbsReal(name, title),
      x_("x", "", this, x),
      funcs_("funcs", "", this),
      coeffs_("coeffs", "", this),
      binpars_("binpars", "", this),
      sentry_(TString(name) + "_sentry", ""),
      binsentry_(TString(name) + "_binsentry", ""),
      initialized_(false),
      last_eval_(-1) {
  funcs_.add(funcs);
  coeffs_.add(coeffs);
}

// allow to pass correlation histogram
CMSHistErrorPropagator::CMSHistErrorPropagator(const char* name,
                                               const char* title,
                                               RooRealVar& x,
                                               RooArgList const& funcs,
                                               RooArgList const& coeffs,
                                               TH3 const& mccorr
                                               )
    : RooAbsReal(name, title),
      x_("x", "", this, x),
      funcs_("funcs", "", this),
      coeffs_("coeffs", "", this),
      binpars_("binpars", "", this),
      mccorr_(mccorr),
      sentry_(TString(name) + "_sentry", ""),
      binsentry_(TString(name) + "_binsentry", ""),
      initialized_(false),
      last_eval_(-1){
  funcs_.add(funcs);
  coeffs_.add(coeffs);

  //let's se if we initialize the thing correctly
  mccorr_.Dump();
}

CMSHistErrorPropagator::CMSHistErrorPropagator(
    CMSHistErrorPropagator const& other, const char* name)
    : RooAbsReal(other, name),
      x_("x", this, other.x_),
      funcs_("funcs", this, other.funcs_),
      coeffs_("coeffs", this, other.coeffs_),
      binpars_("binpars", this, other.binpars_),
      bintypes_(other.bintypes_),
      sentry_(name ? TString(name) + "_sentry" : TString(other.GetName())+"_sentry", ""),
      binsentry_(name ? TString(name) + "_binsentry" : TString(other.GetName())+"_binsentry", ""),
      initialized_(false),
      last_eval_(-1) {
}

void CMSHistErrorPropagator::initialize() const {
  if (initialized_) return;
  sentry_.SetName(TString(this->GetName()) + "_sentry");
  binsentry_.SetName(TString(this->GetName()) + "_binsentry");
#if HFVERBOSE > 0
  std::cout << "Initializing vectors\n";
#endif
  unsigned nf = funcs_.getSize();
  vfuncs_.resize(nf);
  vcoeffs_.resize(nf);
  for (unsigned i = 0; i < nf; ++i) {
    vfuncs_[i] = dynamic_cast<CMSHistFunc const*>(funcs_.at(i));
    vcoeffs_[i] = dynamic_cast<RooAbsReal const*>(coeffs_.at(i));
    auto sargs = vfuncs_[i]->getSentryArgs();
    sentry_.addVars(*sargs);
    std::cout << vfuncs_[i]->GetName() << std::endl;
  }
  unsigned nb = vfuncs_[0]->cache().size();
  std::cout << "Number of bins: " << nb << std::endl;
  vbinpars_.resize(nb);
  if (bintypes_.size()) {
    std::cout << "Bin types are already defined!" << std::endl;
    for (unsigned j = 0, r = 0; j < nb; ++j) {
      vbinpars_[j].resize(bintypes_[j].size(), nullptr);
      for (unsigned i = 0; i < bintypes_[j].size(); ++i) {
        if (bintypes_[j][i] >= 1 && bintypes_[j][i] < 4) {
          vbinpars_[j][i] = dynamic_cast<RooAbsReal *>(binpars_.at(r));
          ++r;
        }
      }
    }
  }
  // the cache of the first vfunc is an instance of type
  // FastHisto defined in FastTemplate_Old
  // It stores the value of a template in a smart and efficient way
  valsum_ = vfuncs_[0]->cache();
  // we get valsum for the first pdf but then we clear it
  // meaning setting every bin yield to zero
  valsum_.Clear();
  // do the same again for the cache object
  cache_ = vfuncs_[0]->cache();
  cache_.Clear();
  // this is an std vector of double initialized empty
  err2sum_.resize(nb, 0.);
  // this is an std vector of double initialized empty
  toterr_.resize(nb, 0.);
  // this is an std vector of std vector of double initialized empty
  // nb -> number of bins, nf -> number of pdfs 
  binmods_.resize(nf, std::vector<double>(nb, 0.));
  // this is an std vector of std vector of double initialized empty
  // nb -> number of bins, nf -> number of pdfs 
  scaledbinmods_.resize(nf, std::vector<double>(nb, 0.));
  // this is an std vector of double initialized empty
  // each pdf has one coefficient
  coeffvals_.resize(nf, 0.);

  // Add to SimpleCacheEntry the coeffs.
  // at line 68 we already add vfuncs_[i]->getSentryArgs() to this obj
  sentry_.addVars(coeffs_);
  // If this is the first initialization from setupBinPars then this is empty
  binsentry_.addVars(binpars_);

  sentry_.setValueDirty();
  binsentry_.setValueDirty();

  initialized_ = true;
}


void CMSHistErrorPropagator::updateCache(int eval) const {
  initialize();

  #if HFVERBOSE > 0
    std::cout << "Update Cache of CMSHistErrorPropagator" << std::endl;
  #endif


// #if HFVERBOSE > 0
//   std::cout << "Sentry: " << sentry_.good() << "\n";
// #endif
  // last eval -1 at construction 
  if (!sentry_.good() || eval != last_eval_) {
    for (unsigned i = 0; i < vfuncs_.size(); ++i) {
      // this calls CMSHistFunc updateCache 
      vfuncs_[i]->updateCache();
      coeffvals_[i] = vcoeffs_[i]->getVal();
    }

    // clear the val sum
    valsum_.Clear();
    // fill sumw2 vector with 0s
    std::fill(err2sum_.begin(), err2sum_.end(), 0.);
    // cycle on number of pdfs (number of processes in a datacard bin)
    for (unsigned i = 0; i < vfuncs_.size(); ++i) {

      //  here coeffvals_ is the value of the multiplicative factor on the template 
      //  times log normal times shape unc and so on
      //  if the EFT template has no lnN or shape associated then this is the W.C. scaling
      //  (1-k**2) for SM, 0.5*k*(k+1) for SM+Li+Qi and 0.5*k*(k-1) for SM-Li+Qi.
      //  I guess we want to include other normalization scaling in the factor also when 
      //  computing the error...

      //      std::cout << "coeff for pdf " << vfuncs_[i]->GetName() << std::endl;
      //      std::cout << vcoeffs_[i]->GetName() << std::endl;
      //      std::cout << vcoeffs_[i]->getVal() << std::endl;
      // valsum stores the sum of all the templates defined in this 
      // datacard bin (cycle over the pdfs, add coeff*value to valsum)
      vectorized::mul_add(valsum_.size(), coeffvals_[i], &(vfuncs_[i]->cache()[0]), &valsum_[0]);
      // err2sum stores the sum in quadrature of the errors for each template defined in this
      // datacard bin (cycle over the pdfs, add coeff^2 * error^2 to err2sum)
      vectorized::mul_add_sqr(valsum_.size(), coeffvals_[i], &(vfuncs_[i]->errors()[0]), &err2sum_[0]);
    }
    // take sqrt of the error quadrature sum for this datacard bin 
    // this is the total error for the overall template in this datacard bin!!
    vectorized::sqrt(valsum_.size(), &err2sum_[0], &toterr_[0]);
    // cache_ is a fast histo and store the overall template made from sig 
    // and bkg summed and multiplied by their scaling coefficients or multipliers
    cache_ = valsum_;

    // if we are evaluating and if bintypes is already initialized
    // (so we are making inference of some kind)
    if (eval == 0 && bintypes_.size()) {
      // cycle on the template bins j
      for (unsigned j = 0; j < valsum_.size(); ++j) {
        // what is bintypes == 1?
        if (bintypes_[j][0] == 1) {
#if HFVERBOSE > 1
          std::cout << "Bin " << j << "\n";
          printf(" | %.6f/%.6f/%.6f\n", valsum_[j], err2sum_[j], toterr_[j]);
#endif
          // cycle on the processes - i in this bin
          for (unsigned i = 0; i < vfuncs_.size(); ++i) {
            // if the overall err2 in this bin is > 0 and if this process is multplied
            // by some factor
            if (err2sum_[j] > 0. && coeffvals_[i] > 0.) {
              // take error e as the error for this process in this bin 
              // and multiply it by the coefficient
              double e =  vfuncs_[i]->errors()[j] * coeffvals_[i];
              // update binmods value for this process in this bin as 
              // the total error (sqrt of the sum in quadrature of errors for all processes contributing to bin j)
              // times the above error squared and divide by the total error sum squared in this bin times the 
              // multiplicative coefficient of process i
              binmods_[i][j] = (toterr_[j] *  e * e) / (err2sum_[j] * coeffvals_[i]);
            } else {
              // if the sum error squared in this bin is =< 0 or if the multiplicative coef is zero for process i 
              // then set binmods to 0 for this process - i and this bin - j
              binmods_[i][j] = 0.;
            }
#if HFVERBOSE > 1
            printf("%.6f   ", binmods_[i][j]);
#endif
          }
#if HFVERBOSE > 1
          printf("\n");
#endif
        }
      }
    }


    sentry_.reset();
    binsentry_.setValueDirty();
  }


  if (!binsentry_.good() || eval != last_eval_) {
    runBarlowBeeston();
    // bintypes might have size == 0 if we never ran setupBinPars()
    for (unsigned j = 0; j < bintypes_.size(); ++j) {
      cache_[j] = valsum_[j];
      if (bintypes_[j][0] == 0) {
        continue;
      } else if (bintypes_[j][0] == 1) {
        double x = vbinpars_[j][0]->getVal();
        cache_[j] += toterr_[j] * x;
        // Only fill the scaledbinmods if we're in eval == 0 mode (i.e. need to
        // propagate to wrappers)
        if (eval == 0) {
          for (unsigned i = 0; i < vfuncs_.size(); ++i) {
            scaledbinmods_[i][j] = binmods_[i][j] * x;
          }
        }
      } else {
        for (unsigned i = 0; i < bintypes_[j].size(); ++i) {
          if (bintypes_[j][i] == 2) {
            // Poisson: this is a multiplier on the process yield
            scaledbinmods_[i][j] = ((vbinpars_[j][i]->getVal() - 1.) *
                 vfuncs_[i]->cache()[j]);
            cache_[j] += (scaledbinmods_[i][j] * coeffvals_[i]);
          } else if (bintypes_[j][i] == 3) {
            // Gaussian This is the addition of the scaled error
            scaledbinmods_[i][j] = vbinpars_[j][i]->getVal() * vfuncs_[i]->errors()[j];
            cache_[j] += (scaledbinmods_[i][j] * coeffvals_[i]);
          }
        }
      }
    }
    cache_.CropUnderflows();
    binsentry_.reset();
  }

  last_eval_ = eval;
}

void CMSHistErrorPropagator::runBarlowBeeston() const {
  if (!bb_.init) return;
  RooAbsArg::setDirtyInhibit(true);

  const unsigned n = bb_.use.size();
  for (unsigned j = 0; j < n; ++j) {
    bb_.dat[j] = data_[bb_.use[j]];
    bb_.valsum[j] = valsum_[bb_.use[j]] * cache_.GetWidth(bb_.use[j]);
    bb_.toterr[j] = toterr_[bb_.use[j]] * cache_.GetWidth(bb_.use[j]);
  }
  // This pragma statement tells (modern) gcc that loop can be safely
  // vectorized
  #pragma GCC ivdep
  for (unsigned j = 0; j < n; ++j) {
    bb_.b[j] = bb_.toterr[j] + (bb_.valsum[j] / bb_.toterr[j]) - bb_.gobs[j];
    bb_.c[j] = bb_.valsum[j] - bb_.dat[j] - (bb_.valsum[j] / bb_.toterr[j]) * bb_.gobs[j];
    bb_.tmp[j] = -0.5 * (bb_.b[j] + copysign(1.0, bb_.b[j]) * std::sqrt(bb_.b[j] * bb_.b[j] - 4. * bb_.c[j]));
    bb_.x1[j] = bb_.tmp[j];
    bb_.x2[j] = bb_.c[j] / bb_.tmp[j];
    bb_.res[j] = std::max(bb_.x1[j], bb_.x2[j]);
  }
  for (unsigned j = 0; j < n; ++j) {
    if (toterr_[bb_.use[j]] > 0.) bb_.push_res[j]->setVal(bb_.res[j]);
  }
  RooAbsArg::setDirtyInhibit(false);
  for (RooAbsArg *arg : bb_.dirty_prop) {
    arg->setValueDirty();
  }
}

void CMSHistErrorPropagator::setAnalyticBarlowBeeston(bool flag) const {
  // Clear it if it's already initialised
  if (bb_.init && flag) return;
  if (bb_.init && !flag) {
    for (unsigned i = 0; i < bb_.push_res.size(); ++i) {
      bb_.push_res[i]->setConstant(false);
    }
    bb_.use.clear();
    bb_.dat.clear();
    bb_.valsum.clear();
    bb_.toterr.clear();
    bb_.err.clear();
    bb_.b.clear();
    bb_.c.clear();
    bb_.tmp.clear();
    bb_.x1.clear();
    bb_.x2.clear();
    bb_.res.clear();
    bb_.gobs.clear();
    bb_.dirty_prop.clear();
    bb_.push_res.clear();
    bb_.init = false;
  }
  if (flag && data_.size()) {
    for (unsigned j = 0; j < bintypes_.size(); ++j) {
      if (bintypes_[j][0] == 1 && !vbinpars_[j][0]->isConstant()) {
        bb_.use.push_back(j);
        double gobs_val = 0.;
        for (RooAbsArg * arg : vbinpars_[j][0]->valueClients()) {
          if (arg == this || arg == &binsentry_) {
            // std::cout << "Skipping " << this << " " << this->GetName() << "\n";
          } else {
            // std::cout << "Adding " << arg << " " << arg->GetName() << "\n";
            bb_.dirty_prop.insert(arg);
            auto as_gauss = dynamic_cast<RooGaussian*>(arg);
            if (as_gauss) {
              auto gobs = dynamic_cast<RooAbsReal*>(as_gauss->findServer(TString(vbinpars_[j][0]->GetName())+"_In"));
              if (gobs) gobs_val = gobs->getVal();
            }
          }
        }
        bb_.gobs.push_back(gobs_val);
        bb_.push_res.push_back((RooRealVar*)vbinpars_[j][0]);
        bb_.push_res.back()->setConstant(true);
      }
    }
    unsigned n = bb_.use.size();
    bb_.dat.resize(n);
    bb_.valsum.resize(n);
    bb_.toterr.resize(n);
    bb_.err.resize(n);
    bb_.b.resize(n);
    bb_.c.resize(n);
    bb_.tmp.resize(n);
    bb_.x1.resize(n);
    bb_.x2.resize(n);
    bb_.res.resize(n);
    bb_.init = true;
  }
}


RooArgList * CMSHistErrorPropagator::setupBinPars(double poissonThreshold) {
  RooArgList * res = new RooArgList();
  if (bintypes_.size()) {
    std::cout << "setupBinPars() already called for " << this->GetName() << "\n";
    return res;
  }

  // First initialize all the storage
  // not needed tbh as alreeady called in updateCache
  initialize();
  // Now fill the bin contents and errors
  updateCache(1); // the arg (1) forces updateCache to fill the caches for all bins

  bintypes_.resize(valsum_.size(), std::vector<unsigned>(1, 0));


  std::cout << std::string(60, '=') << "\n";
  std::cout << "Analyzing bin errors for: " << this->GetName() << "\n";
  std::cout << "Poisson cut-off: " << poissonThreshold << "\n";
  std::set<unsigned> skip_idx;
  std::vector<std::string> skipped_procs;
  for (unsigned i = 0; i < vfuncs_.size(); ++i) {
    if (vfuncs_[i]->attributes().count("skipForErrorSum")) {
      skipped_procs.push_back(vfuncs_[i]->getStringAttribute("combine.process"));
      skip_idx.insert(i);
    }
  }
  if (skipped_procs.size()) {
    std::cout << "Processes excluded for sums:";
    for (auto &s: skipped_procs) std::cout << " " << s;
    std::cout << "\n";
  }
  std::cout << std::string(60, '=') << "\n";
  std::cout << TString::Format("%-10s %-15s %-15s %-30s\n", "Bin", "Contents", "Error", "Notes");

  // cycle on the bins in this datacard bin 
  // valsum[j] is the sum of all templates defined in the datacard multiplied by their
  // multiplicative coefficients
  // toterr is the sqrt of the sum in quadrature of all the templates errors multiplied by their multiplicative 
  // coefficients.
  // There is only one valsum histo and one toterr histo per datacard bin (this->GetName())
  for (unsigned j = 0; j < valsum_.size(); ++j) {
    std::cout << TString::Format("%-10i %-15f %-15f %-30s\n", j, valsum_[j], toterr_[j], "total sum");
    double sub_sum = 0.;
    double sub_err = 0.;
    // Check using a possible sub-set of processes
    // cycling on vfuncs meaning take the pdfs of the processes 
    // contributing in this analysis bin
    for (unsigned i = 0; i < vfuncs_.size(); ++i) {
      // if we skip the process e.g. if it is a signal process then continue
      if (skip_idx.count(i)) {
        continue;
      }
      // if not skipped compute the new corrected valsum for BB approach 
      // by summing the bin content at bin - j multiplied by scalers
      sub_sum += vfuncs_[i]->cache()[j] * coeffvals_[i];
      // add the error for this bin and this process in quadrature 
      // multiplied by the scaler
      sub_err += std::pow(vfuncs_[i]->errors()[j] * coeffvals_[i], 2.);;
    }
    // take the sqrt of the overall corrected error ignoring some processes
    sub_err = std::sqrt(sub_err);
    if (skipped_procs.size()) {
      std::cout << TString::Format("%-10i %-15f %-15f %-30s\n", j, sub_sum, sub_err, "excluding marked processes");
    }
    // Don't do the splitting if the total error in this bin is zero or less
    if (sub_err <= 0.) {
      std::cout << TString::Format("  %-30s\n", "=> Error is zero, ignore");
      std::cout << std::string(60, '-') << "\n";
      continue;
    }

    // Now check if we are below the poisson threshold
    double n = int(0.5 + ((sub_sum * sub_sum) / (sub_err * sub_err)));
    double alpha = valsum_[j] / n;
    std::cout << TString::Format(
        "%-10i %-15f %-15f %-30s\n", j, n, std::sqrt(n),
        TString::Format("Unweighted events, alpha=%f", alpha).Data());

    if (n <= poissonThreshold) {
      std::cout << TString::Format("  %-30s\n", "=> Number of weighted events is below Poisson threshold");

      bintypes_[j].resize(vfuncs_.size(), 4);

      for (unsigned i = 0; i < vfuncs_.size(); ++i) {
        std::string proc =
            vfuncs_[i]->stringAttributes().count("combine.process")
                ? vfuncs_[i]->getStringAttribute("combine.process")
                : vfuncs_[i]->GetName();
        double v_p = vfuncs_[i]->cache()[j];
        double e_p = vfuncs_[i]->errors()[j];
        std::cout << TString::Format("    %-20s %-15f %-15f %-30s\n", proc.c_str(), v_p, e_p, "");
        // relax the condition of v_p >= e_p slightly due to numerical rounding...
        // Possibilities:
        //    v_p = any   e_p <= 0       : skip
        //    v_p < 0     e_p > 0        : for now, skip but technically we should be able to handle this in the future
        //    v_p >= 0    e_p > v_p      : Create an additive gaussian constraint for this bin
        //    v_p > 0     0 < e_p <= v_p : do the poisson
        if (e_p <= 0.) {
          std::cout << TString::Format("      %-30s\n", "=> Error is zero, ignore");
          bintypes_[j][i] = 4;
        } else if (v_p < 0. && e_p > 0.) {
          std::cout << TString::Format("      %-30s\n", "=> Cannot handle negative content, ignore");
          bintypes_[j][i] = 4;
        } else if (v_p > 0. && e_p > 0. && v_p >= (e_p*0.999)) {
          double n_p_r = int(0.5 + ((v_p * v_p) / (e_p * e_p)));
          double alpha_p_r = v_p / n_p_r;
          std::cout << TString::Format(
              "    %-20s %-15f %-15f %-30s\n", "", n_p_r, std::sqrt(n_p_r),
              TString::Format("Unweighted events, alpha=%f", alpha_p_r).Data());
          if (n_p_r <= poissonThreshold) {
            double sigma = 7.;
            double rmin = 0.5*ROOT::Math::chisquared_quantile(ROOT::Math::normal_cdf_c(sigma), n_p_r * 2.);
            double rmax = 0.5*ROOT::Math::chisquared_quantile(1. - ROOT::Math::normal_cdf_c(sigma), n_p_r * 2. + 2.);
            RooRealVar *var = new RooRealVar(TString::Format("%s_bin%i_%s", this->GetName(), j, proc.c_str()), "", n_p_r, rmin, rmax);
            RooConstVar *cvar = new RooConstVar(TString::Format("%g", 1. / n_p_r), "", 1. / n_p_r);
            RooProduct *prod = new RooProduct(TString::Format("%s_prod", var->GetName()), "", RooArgList(*var, *cvar));
	    RooArgSet ownedComps;
	    ownedComps.add(*prod);
	    ownedComps.add(*cvar);
            var->addOwnedComponents(ownedComps);
            var->setAttribute("createPoissonConstraint");
            res->addOwned(*var);
            binpars_.add(*prod);

            std::cout << TString::Format(
                "      => Product of %s[%.2f,%.2f,%.2f] and const [%.4f] to be Poisson constrained\n",
                var->GetName(), var->getVal(), var->getMin(), var->getMax(), cvar->getVal());
            bintypes_[j][i] = 2;
          } else {
            RooRealVar *var = new RooRealVar(TString::Format("%s_bin%i_%s", this->GetName(), j, proc.c_str()), "", 0, -7, 7);
            std::cout << TString::Format(
                "      => Parameter %s[%.2f,%.2f,%.2f] to be Gaussian constrained\n",
                var->GetName(), var->getVal(), var->getMin(), var->getMax());
            var->setAttribute("createGaussianConstraint");
            res->addOwned(*var);
            binpars_.add(*var);
            bintypes_[j][i] = 3;
          }
        } else if (v_p >= 0 && e_p > v_p) {
          RooRealVar *var = new RooRealVar(TString::Format("%s_bin%i_%s", this->GetName(), j, proc.c_str()), "", 0, -7, 7);
          std::cout << TString::Format(
              "      => Poisson not viable, %s[%.2f,%.2f,%.2f] to be Gaussian constrained\n",
              var->GetName(), var->getVal(), var->getMin(), var->getMax());
          var->setAttribute("createGaussianConstraint");
          res->addOwned(*var);
          binpars_.add(*var);
          bintypes_[j][i] = 3;
        } else{
          std::cout << "      => ERROR: shouldn't be here\n";
        }
        std::cout << "  " << std::string(58, '-') << "\n";

      }
    } else if (toterr_[j] > 0.) {
      bintypes_[j][0] = 1;
      RooRealVar *var = new RooRealVar(TString::Format("%s_bin%i", this->GetName(), j), "", 0, -7, 7);
      std::cout << TString::Format(
          "  => Total parameter %s[%.2f,%.2f,%.2f] to be Gaussian constrained\n",
          var->GetName(), var->getVal(), var->getMin(), var->getMax());
      var->setAttribute("createGaussianConstraint");
      var->setAttribute("forBarlowBeeston");
      res->addOwned(*var);
      binpars_.add(*var);
    }
    std::cout << std::string(60, '-') << "\n";
  }

  // binpars_.add(*res);
  binsentry_.addVars(binpars_);
  binsentry_.setValueDirty();

  for (unsigned j = 0, r = 0; j < valsum_.size(); ++j) {
    vbinpars_[j].resize(bintypes_[j].size());
    for (unsigned i = 0; i < bintypes_[j].size(); ++i) {
      if (bintypes_[j][i] >= 1 && bintypes_[j][i] < 4) {
        vbinpars_[j][i] = dynamic_cast<RooAbsReal *>(binpars_.at(r));
        ++r;
      }
    }
  }
  return res;
}


void CMSHistErrorPropagator::applyErrorShifts(unsigned idx,
                                              FastHisto const& nominal,
                                              FastHisto& result) {
  // We can skip the whole evaluation if there's nothing to evaluate
  // if (bintypes_.size() == 0) return;
  // std::cout << "Start of function\n";
  updateCache(0);
  for (unsigned i = 0; i < result.size(); ++i) {
    result[i] = nominal[i] + scaledbinmods_[idx][i];
  }
}

std::unique_ptr<RooArgSet> CMSHistErrorPropagator::getSentryArgs() const {
  // We can do this without initialising because we're going to hand over
  // the sentry directly
  sentry_.SetName(TString(this->GetName()) + "_sentry");
  binsentry_.SetName(TString(this->GetName()) + "_binsentry");
  std::unique_ptr<RooArgSet> args(new RooArgSet(sentry_, binsentry_));
  return args;
}

Double_t CMSHistErrorPropagator::evaluate() const {
  updateCache(1);
  return cache().GetAt(x_);
}


void CMSHistErrorPropagator::printMultiline(std::ostream& os, Int_t contents, Bool_t verbose,
                    TString indent) const {
  RooAbsReal::printMultiline(os, contents, verbose, indent);
  updateCache();
  if (bintypes_.size()) {
    std::cout << ">> Bin types are:\n";
    for (unsigned j = 0; j < bintypes_.size(); ++j) {
      std::cout << ">> " << j << ":";
      for (unsigned i = 0; i < bintypes_[j].size(); ++i) {
        std::cout << " " << bintypes_[j][i];
      }
      std::cout << "\n";
    }
  }
  std::cout << ">> Current cache:\n";
  valsum_.Dump();
  std::cout << ">> Current cache (bin scaled):\n";
  cache_.Dump();
  std::cout << ">> Sentry: " << sentry_.good() << "\n";
  sentry_.Print("v");

}

Int_t CMSHistErrorPropagator::getAnalyticalIntegral(RooArgSet& allVars,
                                         RooArgSet& analVars,
                                         const char* /*rangeName*/) const {
  if (matchArgs(allVars, analVars, x_)) return 1;
  return 0;
}

Double_t CMSHistErrorPropagator::analyticalIntegral(Int_t code,
                                         const char* rangeName) const {
  // TODO: check how RooHistFunc handles ranges that splice bins
  switch (code) {
    case 1: {
      updateCache(1);
      return cache().IntegralWidth();
    }
  }

  assert(0);
  return 0;
}

void CMSHistErrorPropagator::setData(RooAbsData const& data) const {
  updateCache(1);
  data_.clear();
  data_.resize(cache_.fullsize(), 0.); // fullsize is important here is we've used activeBins
  RooArgSet obs(x_.arg());
  const RooRealVar& x = static_cast<const RooRealVar&>(*obs.first());
  for (int i = 0, n = data.numEntries(); i < n; ++i) {
    obs = *data.get(i);
    int idx = cache_.FindBin(x.getVal());
    data_[idx] = data.weight();
  }
}

RooArgList CMSHistErrorPropagator::wrapperList() const {
  RooArgList result;
  for (int i = 0; i < funcs_.getSize(); ++i) {
    CMSHistFunc const* hf = dynamic_cast<CMSHistFunc const*>(funcs_.at(i));
    if (hf) {
      CMSHistFuncWrapper const* wrapper = hf->wrapper();
      if (wrapper) result.add(*wrapper);
    }
  }
  return result;
}

std::map<std::string, Double_t> CMSHistErrorPropagator::getProcessNorms() const {

      std::map<std::string, Double_t> vals_;
      RooArgList clist(coefList());
      RooArgList plist(funcList());
      /*if (plist.getSize() == 1) {
         CMSHistErrorPropagator *err = dynamic_cast<CMSHistErrorPropagator*>(plist.at(0));
         if (err) {
           clist.removeAll();
           plist.removeAll();
           clist.add(err->coefList());
           plist.add(err->wrapperList());
         }
      }
      */
      for (int i = 0, n = clist.getSize(); i < n; ++i) {
        RooAbsReal *coeff = (RooAbsReal *) clist.at(i);
        std::string coeffName = coeff->GetName();
        RooAbsReal* shape = (RooAbsReal*)plist.at(i);
        std::unique_ptr<RooArgSet> myobs(shape->getObservables(*x_));
        TString normProdName = TString::Format("%s", coeff->GetName());
        RooAbsReal * normProd = nullptr;
        if (coeff->ownedComponents()) {
          normProd = dynamic_cast<RooAbsReal*>(coeff->ownedComponents()->find(normProdName));
        }
        if (!normProd) {
          RooAbsReal* integral = shape->createIntegral(*myobs);
      	  RooArgList normProdInputs;
      	  normProdInputs.add(*integral);
      	  normProdInputs.add(*coeff);
          normProd = new RooProduct(normProdName, "", normProdInputs);
          normProd->addOwnedComponents(normProdInputs);
        }
        vals_[normProdName.Data()] = normProd->getVal();
      }
      return vals_;
}
#undef HFVERBOSE

