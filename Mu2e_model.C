#include "Systematics.C"

using namespace FCSys;

TRandom3* rnd_ = new TRandom3(90);
bool comparePDF_ = false; //make a plot comparing full random sample to semi-random sampling
bool doConstraints_ = true; //use systematic uncertainties in limit calculation
bool useLogNormal_ = true; //use log-normal systematic uncertainty PDFs
double scaleLuminosity_ = 1.; //scale luminosity for testing purposes

void Mu2e_model() {

  ///////////////////////////////////////////////////////
  // Initialize experimental parameters
  ///////////////////////////////////////////////////////

  //All in terms of nominal and fractional uncertainty on the nominal
  const double dio_bkg = 0.0096; //see docdb-36476 Table 9
  const double dio_frac_unc = 0.834; //using systematic upper estimate with statistical error
  const double rpc_bkg = 7.06e-4; //see docdb-36503 section 11.2
  const double rpc_frac_unc = 0.16; //using systematic upper estimate with statistical error
  const double extinction = 1.; //in units of 10^-10, FIXME: Get expected value
  const double rpc_oot_bkg = 13.9e-4*extinction; //see docdb-36503 section 11.3
  const double rpc_oot_frac_unc = 0.12; //using systematic upper estimate with statistical error
  const double pbar_bkg = 0.069 * 3.77e19/3.6e20; //see docdb-36494 eq 60
  const double pbar_frac_unc = 1.; //100% uncertainty quoted - careful with Gaussian mode!
  const double cr_lo_bkg = 0.03; //see docdb-38052 slide 20
  const double cr_lo_frac_unc = 0.20;
  const double cr_hi_bkg = 0.02; //see docdb-38052 slide 20
  const double cr_hi_frac_unc = 0.50;
  const double lumi_frac_unc = 0.1;
  const double signal_acceptance = 0.1114; //see docdb-36491 Table 5
  const double sig_frac_unc = 0.0715; //taken as the average of the momentum scale, but should be two-sided and correlated with DIO
  const double ses = 1./(scaleLuminosity_*3.77e19*1.59e-3*0.609); //for signal strength -> R_mue

  ///////////////////////////////////////////////////////
  // Initialize model parameters
  ///////////////////////////////////////////////////////

  var_t lumi_unc("Luminosity uncertainty", lumi_frac_unc);
  var_t lumi_beta("Luminosity beta", 0., -10., 10.);
  var_t lumi_var("Luminosity variation", 0., -1., 5.);
  var_t lumi("Luminosity", (scaleLuminosity_ > 0.) ? scaleLuminosity_ : 1., 0., 5.);
  if(useLogNormal_) {
    lumi_var.nom_ = 1.; lumi_var.val_ = 1.;
    lumi_var.set_sys({&lumi_unc}, {}, {&lumi_beta}); //multiply (1 + uncertainty) ^ unit width gaussian to lumi prediction
    lumi.set_sys({}, {&lumi_var});
  } else {
    lumi_var.set_sys({&lumi_unc}, {&lumi_beta}); //add uncertainty * unit width gaussian to lumi prediction
    lumi.set_sys({&lumi_var}, {});
  }

  var_t dio_unc("DIO uncertainty", dio_bkg*dio_frac_unc);
  var_t dio_beta("DIO beta", 0., -10., 10.);
  var_t dio_var("DIO variation", 0., -1.*dio_bkg, 5.);
  var_t dio("DIO expectation", dio_bkg, 0., 5.);
  if(useLogNormal_) {
    dio_var.nom_ = 1.; dio_var.val_ = 1.;
    dio_var.set_sys({&dio_unc}, {}, {&dio_beta});
    dio.set_sys({}, {&dio_var, &lumi});
  } else {
    dio_var.set_sys({&dio_unc}, {&dio_beta});
    dio.set_sys({&dio_var}, {&lumi});
  }

  var_t rpc_unc("RPC uncertainty", rpc_bkg*rpc_frac_unc);
  var_t rpc_beta("RPC beta", 0., -10., 10.);
  var_t rpc_var("RPC variation", 0., -1.*rpc_bkg, 5.);
  var_t rpc("RPC expectation", rpc_bkg, 0., 5.);
  if(useLogNormal_) {
    rpc_var.nom_ = 1.; rpc_var.val_ = 1.;
    rpc_var.set_sys({&rpc_unc}, {}, {&rpc_beta});
    rpc.set_sys({}, {&rpc_var, &lumi});
  } else {
    rpc_var.set_sys({&rpc_unc}, {&rpc_beta});
    rpc.set_sys({&rpc_var}, {&lumi});
  }

  var_t rpc_oot_unc("RPC OOT uncertainty", rpc_oot_bkg*rpc_oot_frac_unc);
  var_t rpc_oot_beta("RPC OOT beta", 0., -10., 10.);
  var_t rpc_oot_var("RPC OOT variation", 0., -1.*rpc_oot_bkg, 5.);
  var_t rpc_oot("RPC OOT expectation", rpc_oot_bkg, 0., 5.);
  if(useLogNormal_) {
    rpc_oot_var.nom_ = 1.; rpc_oot_var.val_ = 1.;
    rpc_oot_var.set_sys({&rpc_oot_unc}, {}, {&rpc_oot_beta});
    rpc_oot.set_sys({}, {&rpc_oot_var, &lumi});
  } else {
    rpc_oot_var.set_sys({&rpc_oot_unc}, {&rpc_oot_beta});
    rpc_oot.set_sys({&rpc_oot_var}, {&lumi});
  }

  var_t pbar_unc("Pbar uncertainty", pbar_bkg*pbar_frac_unc);
  var_t pbar_beta("Pbar beta", 0., -10., 10.);
  var_t pbar_var("Pbar variation", 0., -1.*pbar_bkg, 5.);
  var_t pbar("Pbar expectation", pbar_bkg, 0., 5.);
  if(useLogNormal_) {
    pbar_var.nom_ = 1.; pbar_var.val_ = 1.;
    pbar_var.set_sys({&pbar_unc}, {}, {&pbar_beta});
    pbar.set_sys({}, {&pbar_var, &lumi});
  } else {
    pbar_var.set_sys({&pbar_unc}, {&pbar_beta});
    pbar.set_sys({&pbar_var}, {&lumi});
  }

  var_t cr_lo_unc("Cosmic-ray (lo) uncertainty", cr_lo_bkg*cr_lo_frac_unc);
  var_t cr_lo_beta("Cosmic-ray (lo) beta", 0., -10., 10.);
  var_t cr_lo_var("Cosmic-ray (lo) variation", 0., -1.*cr_lo_bkg, 5.);
  var_t cr_lo("Cosmic-ray (lo) expectation", cr_lo_bkg, 0., 5.);
  if(useLogNormal_) {
    cr_lo_var.nom_ = 1.; cr_lo_var.val_ = 1.;
    cr_lo_var.set_sys({&cr_lo_unc}, {}, {&cr_lo_beta});
    cr_lo.set_sys({}, {&cr_lo_var});
  } else {
    cr_lo_var.set_sys({&cr_lo_unc}, {&cr_lo_beta});
    cr_lo.set_sys({&cr_lo_var}, {});
  }

  var_t cr_hi_unc("Cosmic-ray (hi) uncertainty", cr_hi_bkg*cr_hi_frac_unc);
  var_t cr_hi_beta("Cosmic-ray (hi) beta", 0., -10., 10.);
  var_t cr_hi_var("Cosmic-ray (hi) variation", 0., -1.*cr_hi_bkg, 5.);
  var_t cr_hi("Cosmic-ray (hi) expectation", cr_hi_bkg, 0., 5.);
  if(useLogNormal_) {
    cr_hi_var.nom_ = 1.; cr_hi_var.val_ = 1.;
    cr_hi_var.set_sys({&cr_hi_unc}, {}, {&cr_hi_beta});
    cr_hi.set_sys({}, {&cr_hi_var});
  } else {
    cr_hi_var.set_sys({&cr_hi_unc}, {&cr_hi_beta});
    cr_hi.set_sys({&cr_hi_var}, {});
  }

  var_t sig_mu("Nominal signal strength", 0., 0., 10./signal_acceptance);
  var_t sig_eff("Signal efficiency", signal_acceptance);
  var_t sig_unc("Signal uncertainty", sig_frac_unc);
  var_t sig_beta("Signal beta", 0., -10., 10.);
  var_t sig_var("Signal variation", 0., -1., 5.);
  var_t signal("Signal expectation", 0., 0., 20.);
  if(useLogNormal_) {
    sig_var.nom_ = 1.; sig_var.val_ = 1.;
    sig_var.set_sys({&sig_unc}, {}, {&sig_beta});
    signal.set_sys({&sig_mu}, {&sig_var, &sig_eff, &lumi});
  } else {
    sig_var.set_sys({&sig_unc}, {&sig_beta, &sig_mu});
    signal.set_sys({&sig_mu, &sig_var}, {&sig_eff, &lumi});
  }

  cout << "Nominal values:\n";
  lumi.print();
  lumi_unc.print();
  lumi_var.print();
  dio.print();
  dio_unc.print();
  dio_var.print();
  rpc.print();
  rpc_unc.print();
  rpc_var.print();
  rpc_oot.print();
  rpc_oot_unc.print();
  rpc_oot_var.print();
  pbar.print();
  pbar_unc.print();
  pbar_var.print();
  cr_lo.print();
  cr_lo_unc.print();
  cr_lo_var.print();
  cr_hi.print();
  cr_hi_unc.print();
  cr_hi_var.print();
  sig_eff.print();
  sig_unc.print();
  sig_var.print();
  signal.print();

  if(!doConstraints_) {
    cout << "Setting systematics to 0!\n";
    lumi_beta.set_constant();
    dio_beta.set_constant();
    rpc_beta.set_constant();
    rpc_oot_beta.set_constant();
    pbar_beta.set_constant();
    cr_lo_beta.set_constant();
    cr_hi_beta.set_constant();
  }
  ///////////////////////////////////////////////////////
  // Initialize model
  ///////////////////////////////////////////////////////

  var_t nobs("Number observed", 0., 0., 20.);
  Poisson_t model("Counting model", nobs,
                 {&dio     , &rpc     , &rpc_oot     , &pbar     , &cr_lo     , &cr_hi     , &signal},
                 {&dio_beta, &rpc_beta, &rpc_oot_beta, &pbar_beta, &cr_lo_beta, &cr_hi_beta, &sig_beta, &lumi_beta}
                 );
  model.ngen_ = 1e5;

  cout << "Model:\n";
  model.Print();
  cout << "Poisson model has a nominal mean of: " << model.GetNominalMean() << endl;

  cout << "Generating the NULL observable PDF\n";
  TH1D* hpdf = model.GeneratePDF(*rnd_);
  hpdf->SetName("null");
  TCanvas* c = new TCanvas();
  hpdf->SetLineColor(kRed);
  hpdf->SetLineWidth(2);
  hpdf->SetMarkerColor(kRed);
  hpdf->SetMarkerStyle(20);
  hpdf->SetMarkerSize(0.8);
  hpdf->Draw();
  hpdf->SetAxisRange(1.e-20, 10., "Y");
  hpdf->SetAxisRange(0, 20, "X");
  c->SetLogy();

  if(comparePDF_) {
    cout << "Randomly sampling the NULL model PDF\n";
    const int nentries = 1e6;
    TH1D* hexp = new TH1D("hexp", "", 20, 0, 20);
    for(int i = 0; i < nentries; ++i) {
      if(i % (nentries/5) == 0) {
        cout << "Samping " << i << ":\n";
        model.SetVerbose(10);
      } else model.SetVerbose(0);
      int n = model.RandomSample(*rnd_);
      if(i % (nentries/10) == 0) {
        signal.print();
        dio.print();
        dio_var.print();
        cout << endl;
      }

      hexp->Fill(n);
    }
    hexp->Scale(1./nentries);
    hexp->Draw("hist E1 sames");
    hexp->SetLineWidth(2);
  }

  ///////////////////////////////////////////////////////
  // Initialize Feldman-Cousins calculator for the model
  ///////////////////////////////////////////////////////

  cout << "Initializing Feldman-Cousins calculator\n";
  FCCalculator fc(model, sig_mu, *rnd_, 0.90/*, 3*/);
  fc.res_ = 1.e-3;

  double mu_min, mu_max;
  int nseen;
  //get the median expected events for the null hypothesis PDF
  nseen = fc.GetMedian(hpdf);
  cout << "Using a single event sensitivity of: " << ses/signal_acceptance << endl;
  cout << "Performing Feldman-Cousins calculation for median nseen = " << nseen << endl;
  fc.CalculateInterval(nseen, mu_min, mu_max);
  printf("For %i seen, R_mue interval is: %.3e - %.3e (%.3f - %.3f mean events)\n",
         nseen, mu_min*ses, mu_max*ses,
         mu_min*sig_eff.nom_*scaleLuminosity_, mu_max*sig_eff.nom_*scaleLuminosity_);

  //get median discovery information
  int ndisc = fc.NSigmaThreshold(hpdf, 5.);
  cout << "N(discovery) for NULL model = " << ndisc << endl;
  double mu_disc = fc.FindForMedianN(ndisc);
  printf("For a median of %i, minimum R_mue is: %.3e (%.3f mean events)\n",
         ndisc, mu_disc*ses, mu_disc*sig_eff.nom_*scaleLuminosity_);
}
