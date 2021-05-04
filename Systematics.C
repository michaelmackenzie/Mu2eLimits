//a variable, knows about its bounds, nominal value, additive (applied first) variables, multiplicative (second), and power (third) variables
namespace FCSys {
  class var_t {
  public:
    var_t(TString name = "Default") : name_(name), min_(-1.), max_(1.), nom_(0.), val_(0.), constant_(false), verbose_(0) {}
    var_t(TString name, double nom) : name_(name), min_(nom-1.), max_(nom+1.), nom_(nom), val_(nom), constant_(false), verbose_(0) {}
    var_t(TString name, double nom, double min, double max, TString pdf = "Gauss") :
      name_(name), min_(min), max_(max), nom_(nom), val_(nom), pdf_(pdf), constant_(false), verbose_(0) {}

    void set_constant(bool constant = true) {
      constant_ = constant;
    }

    void set_sys(std::vector<var_t*> add, std::vector<var_t*> mul, std::vector<var_t*> pow = {}) {
      add_ = add;
      mul_ = mul;
      pow_ = pow;
    }

    void set_val(double val) {
      val_ = std::min(max_, std::max(min_, val));
    }

    void set_rnd_val(TRandom3& rnd) {
      if(constant_) return;
      if(pdf_ == "Gauss") {
        set_val(rnd.Gaus(nom_));
      }
    }

    double get_val() {
      double val = val_;
      if(verbose_ > 2) printf("Variable %s has starting value %.3e\n", name_.Data(), val);
      for(var_t* var : add_) {
        val += var->get_val();
        if(verbose_ > 2) printf("Variable %s add %s = %.3e --> %.3e\n", name_.Data(), var->name_.Data(), var->get_val(), val);
      }
      for(var_t* var : mul_) {
        val *= var->get_val();
        if(verbose_ > 2) printf("Variable %s mul %s = %.3e --> %.3e\n", name_.Data(), var->name_.Data(), var->get_val(), val);
      }
      for(var_t* var : pow_) {
        val = std::pow(val, var->get_val());
        if(verbose_ > 2) printf("Variable %s pow %s = %.3e --> %.3e\n", name_.Data(), var->name_.Data(), var->get_val(), val);
      }
      if(val > max_ || val < min_) {
        if(verbose_ > 0) std::cout << "Variable " << name_.Data() << " value outside of bounds! Returning bound...\n";
        return std::min(max_, std::max(min_, val));
      }
      return val;
    }

    void print() {
      printf(" %s: %.3e (%.3e) [%.3e - %.3e]", name_.Data(), get_val(), nom_, min_, max_);
      if(add_.size() > 0) {
        printf(" add = {");
        for(var_t* var : add_) {
          printf("%s", var->name_.Data());
          if(var != add_[add_.size()-1]) printf(", ");
        }
        printf("}");
      }
      if(mul_.size() > 0) {
        printf(" mul = {");
        for(var_t* var : mul_) {
          printf("%s", var->name_.Data());
          if(var != mul_[mul_.size()-1]) printf(", ");
        }
        printf("}");
      }
      if(pow_.size() > 0) {
        printf(" pow = {");
        for(var_t* var : pow_) {
          printf("%s", var->name_.Data());
          if(var != pow_[pow_.size()-1]) printf(", ");
        }
        printf("}");
      }
      printf("\n");
    }

    TString name_;
    double min_;
    double max_;
    double nom_;
    double val_;
    TString pdf_;
    bool constant_;
    int verbose_;
    std::vector<var_t*> add_;
    std::vector<var_t*> mul_;
    std::vector<var_t*> pow_;
  };

  //class for a Poission PDF with systematics
  class Poisson_t {
  public:

    Poisson_t(TString name, var_t& obs, std::vector<var_t*> mu, std::vector<var_t*> sys = {}) :
      name_(name), obs_(obs), mu_(mu), sys_(sys), verbose_(0) {}

    void SetVerbose(int verbose) {
      verbose_ = verbose;
    }

    double GetMean() {
      double mu = 0.;
      for(var_t* var : mu_) {
        mu += var->get_val();
        if(verbose_ > 9) std::cout << __func__ << ": mu(" << var->name_.Data() << ") = " << var->get_val() << std::endl;
      }
      if(verbose_ > 9) std::cout << __func__ << ": mu = " << mu << std::endl;
      return mu;
    }

    double GetNominalMean() {
      double mu = 0.;
      for(var_t* var : sys_) {
        var->val_ = var->nom_;
      }
      for(var_t* var : mu_) {
        var->val_ = var->nom_;
        mu += var->get_val();
      }
      if(verbose_ > 9) std::cout << __func__ << ": mu = " << mu << std::endl;
      return mu;
    }

    double Eval(int n) {
      const double mu = GetMean();
      const double p = ROOT::Math::poisson_pdf(n, mu);
      if(verbose_ > 9) std::cout << __func__ << ": p = " << p << std::endl;
      return p;
    }

    int RandomSample(TRandom3& rnd) {
      //first sample the systematics
      for(var_t* var : sys_) {
        var->set_rnd_val(rnd);
        if(verbose_ > 9) {std::cout << __func__ << ":\n"; var->print();}
      }
      //next sample the poisson distribution
      const double mu = GetMean();
      const int n = rnd.Poisson(mu);
      if(verbose_ > 9) std::cout << __func__ << ": n = " << n << std::endl;
      return n;
    }

    TH1D* GeneratePDF(TRandom3& rnd) {
      //sample the nuisance parameters to define a mean, then add a poisson PDF for this
      int nbins = 100;
      TH1D* hpdf = new TH1D("hpdf", "PDF", nbins, 0., (double) nbins);
      const int nattempts = 1e5;
      for(int attempt = 0; attempt < nattempts; ++attempt) {
        RandomSample(rnd);
        const double mu = GetMean();
        for(int n = 0; n < nbins; ++n) {
          hpdf->Fill(n, ROOT::Math::poisson_pdf(n, mu));
        }
      }
      hpdf->Scale(1. / nattempts);
      return hpdf;
    }

    void Print() {
      printf(" %s: obs = %s, mu = {", name_.Data(), obs_.name_.Data());
      for(var_t* var : mu_) {
        printf("%s", var->name_.Data());
          if(var != mu_[mu_.size()-1]) printf(", ");
      }
      printf("} sys = {");
      for(var_t* var : sys_) {
        printf("%s", var->name_.Data());
        if(var != sys_[sys_.size()-1]) printf(", ");
      }
      printf("}\n");
    }

    TString name_;
    var_t& obs_;
    std::vector<var_t*> mu_;
    std::vector<var_t*> sys_;
    int verbose_;
  };

  class FCCalculator {
  public:
    FCCalculator(Poisson_t& model, var_t& poi, TRandom3& rnd, double cl = 0.9, int verbose = 0) :
      model_(model), poi_(poi), rnd_(rnd), cl_(cl), verbose_(verbose) {
      double oldval = poi.val_;
      poi.val_ = 0.;
      hNull_ = model.GeneratePDF(rnd);
      hNull_->SetName("FC_NULL");
      poi.val_ = oldval;
      CalculateIndividualInterval(hNull_, null_min_, null_max_);
      if(null_min_ > 0) {
        printf("!!! %s: N(obs) = 0 is not contained within the NULL hypothesis!\n", __func__);
      }
      if(verbose_ > 0) printf("%s: Null interval at %.3f CL: %i - %i\n", __func__, cl, null_min_, null_max_);
    }

    //for a given PDF, construct the FC interval in N(observed)
    void CalculateIndividualInterval(TH1D* hPDF, int& nmin, int& nmax) {
      nmin = hPDF->GetNbinsX();
      nmax = 0;
      double p = 0.;
      //mu_bkg for the denominator of the likelihood ratio ordering parameter
      const double mu = hNull_->GetMean() - 0.5; //bins are centered at n + 0.5
      if(hPDF->Integral() < cl_) {
        printf("!!! %s: PDF doesn't have the range to calculate the confidence level precision!\n", __func__);
        nmax = nmin;
        nmin = 0;
        return;
      }
      std::set<int> ns;
      if(verbose_ > 2) printf("%s: Printing interval construction:\n", __func__);
      while(p < cl_) {
        int nbest = -1;
        double rbest = -1.;
        for(int n = 0; n < hPDF->GetNbinsX(); ++n) {
          if(ns.count(n)) continue;
          double r = hPDF->GetBinContent(n+1) / ROOT::Math::poisson_pdf(n, std::max(mu, (double) n));
          if(r > rbest) {
            rbest = r;
            nbest = n;
          }
        }
        p += hPDF->GetBinContent(nbest+1);
        nmin = std::min(nbest, nmin);
        nmax = std::max(nbest, nmax);
        ns.insert(nbest);
        if(verbose_ > 9) printf(" n = %i, r = %.3e, p = %.3f\n", nbest, rbest, p);
      }
      if(verbose_ > 2) printf(" Final interval: %i - %i\n", nmin, nmax);
    }

    double FindLimit(int nobs, bool upperLimit) {
      int attempts = 0;
      const int maxAttempts = 100;
      double mu_min = 0.;
      double mu_max = poi_.max_;
      double mu_range = poi_.max_ - poi_.min_;
      while(abs((mu_max - mu_min) / mu_range) > 0.001 && attempts < maxAttempts) {
        ++attempts;
        double mu = (mu_max + mu_min) / 2.;
        poi_.val_ = mu;
        if(verbose_ > 2) model_.Print();
        TH1D* h = model_.GeneratePDF(rnd_);
        int nmin, nmax;
        CalculateIndividualInterval(h, nmin, nmax);
        if(upperLimit) {
          if(nmin > nobs) mu_max = mu; //gone past allowed
          else mu_min = mu; //still allowed
        } else {
          if(nmax < nobs) mu_min = mu; //gone past allowed
          else mu_max = mu; //still allowed
        }
        delete h;
        if(verbose_ > 1) printf("%s: Attempt %i has bounds %.3f - %.3f\n", __func__, attempts, mu_min, mu_max);
      }
      if(attempts == maxAttempts) printf("!!! %s: Hit maximum limit finding attempts for N(obs) = %i, mu in %.3f - %.3f (upper = %o)\n",
                                         __func__, nobs, mu_min, mu_max, upperLimit);
      return (mu_max + mu_min) / 2.;
    }

    void CalculateInterval(int nobs, double& mu_min, double& mu_max) {
      //first test if mu = 0 is included in the interval
      if(nobs > null_max_) {
        mu_min = FindLimit(nobs, false);
      } else {
        mu_min = 0.;
      }
      mu_max = FindLimit(nobs, true);
    }

    Poisson_t& model_;
    var_t& poi_;
    TRandom3& rnd_;
    double cl_;
    int verbose_;
    TH1D* hNull_;
    int null_min_;
    int null_max_;
  };
}
