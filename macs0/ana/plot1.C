
map<string, TFile*> files;

// ----------------------------------------------------------------------
void openFiles() {
  if (0 == files.size()) {
    cout << "open files" << endl;
    files.insert(make_pair("5keV", TFile::Open("g4run-005kev.default.root")));
    files.insert(make_pair("10keV", TFile::Open("g4run-010kev.default.root")));
    files.insert(make_pair("20keV", TFile::Open("g4run-020kev.default.root")));
    files.insert(make_pair("50keV", TFile::Open("g4run-050kev.default.root")));
    files.insert(make_pair("100keV", TFile::Open("g4run-100kev.default.root")));
    files.insert(make_pair("200keV", TFile::Open("g4run-200kev.default.root")));
    files.insert(make_pair("500keV", TFile::Open("g4run-500kev.default.root")));
  } else {
    cout << "files already opened" << endl;
  }
}

// ----------------------------------------------------------------------
void muProd() {
  openFiles();

  TH1D *h1(0);

  TH1D *hcriv = new TH1D("hcriv", "", 7, 0., 7.);
  TH1D *hexpe = new TH1D("expeDecay", "", 7, 0., 7.);
  hexpe->SetLineColor(kBlack);
  hexpe->SetFillColor(kYellow);
  hexpe->SetFillStyle(kSolid);
  hexpe->GetXaxis()->SetBinLabel(1, "5 keV");
  hexpe->GetXaxis()->SetBinLabel(2, "10 keV");
  hexpe->GetXaxis()->SetBinLabel(3, "20 keV");
  hexpe->GetXaxis()->SetBinLabel(4, "50 keV");
  hexpe->GetXaxis()->SetBinLabel(5, "100 keV");
  hexpe->GetXaxis()->SetBinLabel(6, "200 keV");
  hexpe->GetXaxis()->SetBinLabel(7, "500 keV");

  shrinkPad(0.15, 0.25);
  setTitles(hexpe, "Beam energy", "N_{Mu}^{CT decay}", 0.06, 1.3, 2.1);


  TH1D *htarg = new TH1D("targDecay", "", 7, 0., 7.); htarg->SetLineColor(kRed);
  TH1D *hwall = new TH1D("wallDecay", "", 7, 0., 7.); hwall->SetLineColor(kCyan);

  TH1D *hafter = new TH1D("afterDecay", "", 7, 0., 7.); hafter->SetLineColor(kBlue);

  double nmu = 50000;

  hcriv->SetBinContent(1, 0.568*nmu); hcriv->SetBinError(1, 0.568*nmu*0.09);
  hcriv->SetBinContent(2, 0.318*nmu); hcriv->SetBinError(2, 0.318*nmu*0.008);
  setHist(hcriv);

  double ntarget(0.), nexp(0.), nwall(0.), ntot(0.);

  vector<string> energies;
  energies.push_back("5keV");
  energies.push_back("10keV");
  energies.push_back("20keV");
  energies.push_back("50keV");
  energies.push_back("100keV");
  energies.push_back("200keV");
  energies.push_back("500keV");

  for (unsigned int i = 0; i < energies.size(); ++i) {
    cout << energies[i] << endl;
    h1 = (TH1D*)files[energies[i]]->Get("h5");
    h1->Draw();
    ntot    = h1->GetSumOfWeights();
    ntarget = h1->Integral(1, h1->FindBin(-400.));
    // nexp    = h1->Integral(h1->FindBin(-390.), h1->FindBin(2081.));
    // nwall   = h1->Integral(h1->FindBin(2091.), h1->FindBin(2200.));
    nexp    = h1->Integral(h1->FindBin(-390.), h1->FindBin(500.));
    nwall   = h1->Integral(h1->FindBin(2091.), h1->FindBin(501.));
    htarg->SetBinContent(i+1, ntarget);
    hexpe->SetBinContent(i+1, nexp);
    hwall->SetBinContent(i+1, nwall);
    hafter->SetBinContent(i+1, (nexp+nwall));
    cout << "ntarget: " << ntarget << endl;
    cout << "nexp:    " << nexp << endl;
    cout << "nwall:   " << nwall << endl;
    cout << "nafter:  " << nexp+nwall << endl;
  }

  gStyle->SetOptStat(0);
  hexpe->SetMaximum(30000);
  //  hexpe->SetMaximum(1.);
  hexpe->Draw("");
  if (0) {
    hwall->Draw("same");
    htarg->Draw("same");
    hafter->Draw("same");
  }

  hcriv->Draw("samee");

  tl->SetTextSize(0.05);
  tl->DrawLatexNDC(0.5, 0.8, "N_{#mu^{+}} = 5 #times 10^{4}");
  tl->SetTextSize(0.04);
  tl->DrawLatexNDC(0.5, 0.74, "(15 nm C foil)");

  TLegend *legg = new TLegend(0.55, 0.5, 0.85, 0.60, "");
  legg->SetFillStyle(0);
  legg->SetBorderSize(0);
  legg->SetTextSize(0.04);
  legg->SetFillColor(0);
  legg->AddEntry(hcriv, "EPJ, C80, 804", "p");
  legg->AddEntry(hexpe, "GEANT4", "l");
  legg->Draw();


  c0.SaveAs("MuDecay-in-Exp-15nm.pdf");
}


// ----------------------------------------------------------------------
void acc(int type = 0) {
  openFiles();

  TH1D *ha(0);
  ha = (TH1D*)files["5keV"]->Get("acc");

  ha->GetXaxis()->SetBinLabel(1, "no hits");
  ha->GetXaxis()->SetBinLabel(2, "CT hits");
  ha->GetXaxis()->SetBinLabel(3, "ET hits");
  ha->GetXaxis()->SetBinLabel(4, "CT&ET hits");
  ha->GetXaxis()->SetRange(1., 4.);
  ha->SetMinimum(0.);



  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  setFilledHist(ha, kBlue, kYellow, 1000, 3);
  ha->Draw();

  tl->SetTextSize(0.05);
  tl->DrawLatexNDC(0.5, 0.8, "N_{#mu^{+}} = 5 #times 10^{4}");

  c0.SaveAs("macs0-acc.pdf");
}

// ----------------------------------------------------------------------
void nhits(int type = 0) {
  openFiles();

  TH1D *hb(0), *hs(0);
  hs = (TH1D*)files["5keV"]->Get(Form("s%d", type));
  hb = (TH1D*)files["5keV"]->Get(Form("b%d", type));

  setHist(hs, kBlack, 20, 1, 3);
  setHist(hb, kRed, 20, 1, 3);

  setTitles(hb, "N_{G4Hit}", "Events", 0.05, 1.1, 2.0);

  gStyle->SetOptTitle(0);
  shrinkPad(0.15, 0.2);
  hb->Draw();
  hs->Draw("same");

  string name("Central tracker");
  if (1 == type) name = "End tracker";

  TLegend *legg = new TLegend(0.35, 0.5, 0.85, 0.65, Form("G4Hits in %s", name.c_str()));
  legg->SetFillStyle(0);
  legg->SetBorderSize(0);
  legg->SetTextSize(0.04);
  legg->SetFillColor(0);
  legg->AddEntry(hs, "direct signal tracks", "l");
  legg->AddEntry(hb, "not signal tracks", "l");
  legg->Draw();
  if (0 == type) {
    c0.SaveAs("trackerHits.pdf");
  } else {
    c0.SaveAs("mcpHits.pdf");
  }

}
