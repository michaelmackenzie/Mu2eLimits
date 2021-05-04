#include "Systematics.C"

using namespace FCSys;

TRandom3* rnd_ = new TRandom3(90);
bool comparePDF_ = false; //make a plot comparing full random sample to semi-random sampling
bool doConstraints_ = true;

void Mu2e_model() {

  ///////////////////////////////////////////////////////
  // Initialize experimental parameters
  ///////////////////////////////////////////////////////

  const double dio_bkg = 0.02; //nominal expected background
  const double dio_frac_unc = 0.10; //fractional uncertainty
  const double cr_bkg = 0.02;
  const double cr_frac_unc = 0.2;
  const double lumi_frac_unc = 0.1;
  const double signal_acceptance = 0.3;
  const double ses = 1./(3.77e19*1.59e-3*signal_acceptance); //for signal strength -> R_mue

  ///////////////////////////////////////////////////////
  // Initialize model parameters
  ///////////////////////////////////////////////////////

  var_t lumi_unc("Luminosity uncertainty", lumi_frac_unc);
  var_t lumi_beta("Luminosity beta", 0., -10., 10.);
  var_t lumi_var("Luminosity variation", 0., -1., 5.);
  lumi_var.set_sys({&lumi_unc}, {&lumi_beta});
  var_t lumi("Luminosity", 1., 0., 5.);
  lumi.set_sys({&lumi_var}, {});

  var_t dio_unc("DIO uncertainty", dio_bkg*dio_frac_unc);
  var_t dio_beta("DIO beta", 0., -10., 10.);
  var_t dio_var("DIO variation", 0., -1.*dio_bkg, 5.);
  dio_var.set_sys({&dio_unc}, {&dio_beta});
  var_t dio("DIO expectation", dio_bkg, 0., 5.);
  dio.set_sys({&dio_var}, {&lumi});

  var_t cr_unc("Cosmic-ray uncertainty", cr_bkg*cr_frac_unc);
  var_t cr_beta("Cosmic-ray beta", 0., -10., 10.);
  var_t cr_var("Cosmic-ray variation", 0., -1.*cr_bkg, 5.);
  cr_var.set_sys({&cr_unc}, {&cr_beta});
  var_t cr("Cosmic-ray expectation", cr_bkg, 0., 5.);
  cr.set_sys({&cr_var}, {});

  var_t sig_mu("Nominal signal strength", 0., 0., 10./signal_acceptance);
  var_t sig_eff("Signal efficiency", signal_acceptance);
  var_t signal("Signal expectation", 0., 0., 20.);
  signal.set_sys({&sig_mu}, {&sig_eff, &lumi});

  cout << "Nominal values:\n";
  lumi.print();
  lumi_unc.print();
  lumi_var.print();
  dio.print();
  // dio.verbose_ = 10;
  dio.get_val();
  // dio.verbose_ = 0;
  dio_unc.print();
  dio_var.print();
  cr.print();
  cr_unc.print();
  cr_var.print();
  signal.print();
  sig_eff.print();

  ///////////////////////////////////////////////////////
  // Initialize model
  ///////////////////////////////////////////////////////

  var_t nobs("Number observed", 0., 0., 10.);
  Poisson_t model("Counting model", nobs, {&dio, &cr, &signal}, {&lumi_beta, &dio_beta, &cr_beta});

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
  double mu_min, mu_max;
  int nseen;

  //get the median expected events for the null hypothesis PDF
  nseen = 0;
  double p = hpdf->GetBinContent(1);
  while(p < 0.5) {
    ++nseen;
    p += hpdf->GetBinContent(nseen+1);
  }

  cout << "Performing Feldman-Cousins calculation for median nseen = " << nseen << endl;
  fc.CalculateInterval(nseen, mu_min, mu_max);
  printf("For %i seen, R_mue interval is: %.3e - %.3e (%.3f - %.3f mean events)\n",
         nseen, mu_min*ses/0.61, mu_max*ses/0.61,
         mu_min*sig_eff.nom_, mu_max*sig_eff.nom_);

}
