#include "plotResults.hh"

#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <iterator>


#include "TROOT.h"
#include "TStyle.h"
#include "TKey.h"
#include "TMath.h"
#include "TMarker.h"
#include "TPad.h"
#include "TRandom3.h"
#include "TString.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TPad.h"
#include "TF1.h"
#include "THStack.h"
#include "TFitResult.h"
#include "TTimeStamp.h"

#include "util.hh"

ClassImp(plotResults)

using namespace std;

// ----------------------------------------------------------------------
plotResults::plotResults(string dir, string files, string cuts, string setup):
  plotClass(dir, files, cuts, setup) {
  plotClass::loadFiles(files);
  plotResults::loadFiles(files);

  changeSetup(dir, "plotResults", setup);
}


// ----------------------------------------------------------------------
plotResults::~plotResults() {
  cout << "plotResults destructor" << endl;
  resetHistograms(true);
  cout << "done with histogram deletion" << endl;
}

// ----------------------------------------------------------------------
void plotResults::resetHistograms(bool deleteThem) {

}

// ----------------------------------------------------------------------
void plotResults::makeAll(string what) {
  fHistFile = TFile::Open(fHistFileName.c_str(), "RECREATE");

  scanAnalyses();
  closeHistFile();
}


// ----------------------------------------------------------------------
void plotResults::scanAnalyses(string mode) {


  if ("sgKinEnergy" == mode) {
    vector<string> energies;
    energies.push_back("20deg-0.010eV");
    energies.push_back("20deg-0.026eV");
    energies.push_back("20deg-0.050eV");
    energies.push_back("20deg-0.1eV");
    energies.push_back("20deg-0.2eV");
    energies.push_back("20deg-0.5eV");
    energies.push_back("20deg-1.0eV");
    energies.push_back("20deg-5.0eV");
    energies.push_back("20deg-10eV");
    energies.push_back("20deg-100eV");
    energies.push_back("20deg-500eV");
    energies.push_back("2deg-1keV");
    energies.push_back("2deg-2keV");
    energies.push_back("2deg-10keV");
    energies.push_back("2deg-20keV");
    energies.push_back("2deg-100keV");
    energies.push_back("2deg-200keV");
    energies.push_back("2deg-1MeV");
    energies.push_back("2deg-5MeV");

    fHistFile->cd();
    TH1D *h2 = new TH1D(Form("acc"), "Acceptance (decay before accelerator)", energies.size(), 0., energies.size());
    setTitles(h2, "Mu Energy", "Acceptance", 0.04, 1.6, 1.2);

    TH1D *h3 = new TH1D(Form("eff"), "Efficiency (hits from e_{atomic} and e_{#mu})", energies.size(), 0., energies.size());
    setTitles(h3, "Mu Energy", "Efficiency", 0.04, 1.6, 1.2);

    TH1D *h3a = new TH1D(Form("ef2"), "Efficiency (#geq2 hits from e_{atomic} and e_{#mu})", energies.size(), 0., energies.size());
    setTitles(h3a, "Mu Energy", "Efficiency", 0.04, 1.6, 1.2);

    TH1D *h3b = new TH1D(Form("ef3"), "Efficiency (#geq hits from e_{atomic} and e_{#mu})", energies.size(), 0., energies.size());
    setTitles(h3b, "Mu Energy", "Efficiency", 0.04, 1.6, 1.2);

    TH1D *h4 = new TH1D(Form("ef4"), "(Efficiency | decay) (hits from e_{atomic} and e_{#mu})", energies.size(), 0., energies.size());
    setTitles(h4, "Mu Energy", "Efficiency", 0.04, 1.6, 1.2);

    TH1D *ht = new TH1D(Form("ht"), "average proper decay time [#mus]", energies.size(), 0., energies.size());
    setTitles(ht, "Mu Energy", "#LTt#GT [#mus]", 0.04, 1.6, 1.2);

    for (unsigned k = 0; k < energies.size(); ++k) {
      string filename = "signal-5Mu-" + energies[k] + ".default.root";

      TFile *f = TFile::Open(filename.c_str());
      TH1D *h1 = (TH1D*)f->Get("acc");
      TH1D *hm = (TH1D*)f->Get("h6");
      TH1D *h10 = (TH1D*)f->Get("h10");
      double nMu = hm->GetEntries();
      double nacc = h1->GetBinContent(h1->FindBin(0.01));
      double neff = h1->GetBinContent(h1->FindBin(3.01));
      double nef2 = h1->GetBinContent(h1->FindBin(13.01));
      double nef3 = h1->GetBinContent(h1->FindBin(23.01));
      double acc = nacc/nMu;
      double acce= dEff(static_cast<int>(nacc), static_cast<int>(nMu));
      double eff = neff/nMu;
      double effe= dEff(static_cast<int>(neff), static_cast<int>(nMu));
      double ef2 = nef2/nMu;
      double ef2e= dEff(static_cast<int>(nef2), static_cast<int>(nMu));
      double ef3 = nef3/nMu;
      double ef3e= dEff(static_cast<int>(nef3), static_cast<int>(nMu));
      double ef4 = neff/nacc;
      double ef4e= dEff(static_cast<int>(neff), static_cast<int>(nacc));
      h2->SetBinContent(k+1, acc);
      h2->SetBinError(k+1, acce);
      h2->GetXaxis()->SetBinLabel(k+1, energies[k].c_str());

      h3->SetBinContent(k+1, eff);
      h3->SetBinError(k+1, effe);
      h3->GetXaxis()->SetBinLabel(k+1, energies[k].c_str());

      h3a->SetBinContent(k+1, ef2);
      h3a->SetBinError(k+1, ef2e);
      h3a->GetXaxis()->SetBinLabel(k+1, energies[k].c_str());

      h3b->SetBinContent(k+1, ef3);
      h3b->SetBinError(k+1, ef3e);
      h3b->GetXaxis()->SetBinLabel(k+1, energies[k].c_str());

      h4->SetBinContent(k+1, ef4);
      h4->SetBinError(k+1, ef4e);
      h4->GetXaxis()->SetBinLabel(k+1, energies[k].c_str());

      ht->SetBinContent(k+1, 1.e6*h10->GetMean());
      ht->SetBinError(k+1, 1.e6*h10->GetMeanError());
      ht->GetXaxis()->SetBinLabel(k+1, energies[k].c_str());


      cout << filename << " nacc: " << nacc << " nMu: " << nMu
	   << " acc: " << acc << " +/- " << acce
	   << " eff: " << eff << " +/- " << effe
	   << endl;
      f->Close();
    }

    gStyle->SetOptTitle(0);

    shrinkPad(0.13, 0.12);
    c0->SetGridx(1);
    c0->SetGridy(1);

    gStyle->SetOptStat(0);
    h2->Draw("texte");
    savePad(Form("acc_%s.pdf", mode.c_str()));

    h3->Draw("texte");
    savePad(Form("eff_%s.pdf", mode.c_str()));

    h3a->Draw("texte");
    savePad(Form("eff2_%s.pdf", mode.c_str()));

    h3b->Draw("texte");
    savePad(Form("eff3_%s.pdf", mode.c_str()));

    // -- overlay efficiency plots
    c0->Clear();
    c0->SetGridx(1);
    c0->SetGridy(1);
    c0->SetLogy(1);
    setHist(h3, kBlack);
    setHist(h3a, kBlue);
    setHist(h3b, kRed);
    h3->Draw("hist");
    h3a->Draw("histsame");
    h3b->Draw("histsame");
    newLegend(0.12, 0.25, 0.47, 0.40, "");
    legg->SetFillStyle(1000);
    legg->AddEntry(h3, "hits from e_{#mu} and e_{atomic}", "l");
    legg->AddEntry(h3a, "#geq 2 hits from e_{#mu} and e_{atomic}", "l");
    legg->AddEntry(h3b, "#geq 4 hits from e_{#mu} and e_{atomic}", "l");
    legg->Draw();
    savePad(Form("effoverlayLog_%s.pdf", mode.c_str()));


    // -- overlay efficiency plots
    c0->Clear();
    c0->SetGridx(1);
    c0->SetGridy(1);
    c0->SetLogy(0);
    setHist(h3, kBlack);
    setHist(h3a, kBlue);
    setHist(h3b, kRed);
    h3->Draw("hist");
    h3a->Draw("histsame");
    h3b->Draw("histsame");
    newLegend(0.12, 0.25, 0.47, 0.40, "");
    legg->SetFillStyle(1000);
    legg->AddEntry(h3, "hits from e_{#mu} and e_{atomic}", "l");
    legg->AddEntry(h3a, "#geq 2 hits from e_{#mu} and e_{atomic}", "l");
    legg->AddEntry(h3b, "#geq 4 hits from e_{#mu} and e_{atomic}", "l");
    legg->Draw();
    savePad(Form("effoverlayLin_%s.pdf", mode.c_str()));

    //    shrinkPad(0.1, 0.12);
    h4->SetMinimum(0.);
    h4->Draw("histe");
    savePad(Form("eff4_%s.pdf", mode.c_str()));

    c0->SetLogy(1);
    ht->SetMinimum(0.01);
    ht->Draw("histe");
    savePad(Form("decaytime_%s.pdf", mode.c_str()));
    c0->SetLogy(0);

  }


  if ("target" == mode) {
    vector<string> target;
    target.push_back("Cfoil");
    target.push_back("aerogel");
    target.push_back("Al");

    map<string, vector<string> > venergy;
    venergy.insert(make_pair("Cfoil", vector<string>{"5keV", "10keV", "15keV", "50keV"}));
    venergy.insert(make_pair("aerogel", vector<string>{"50keV", "500keV", "5MeV", "26MeV"}));
    venergy.insert(make_pair("Al", vector<string>{"1MeV", "5MeV", "10MeV"}));
    map<string, vector<string> > vthickness;
    vthickness.insert(make_pair("Cfoil", vector<string>{"5nm", "10nm", "15nm"}));
    vthickness.insert(make_pair("aerogel", vector<string>{"10um", "100um", "1mm", "8mm"}));
    vthickness.insert(make_pair("Al", vector<string>{"1500nm", "750nm", "500nm"}));

    makeCanvas(4);
    c3->cd();
    shrinkPad(0.1, 0.25);
    gStyle->SetOptStat(0);

    string filename("");
    for (unsigned int i = 0; i < target.size(); ++i) {
      cout << "target = " << target[i] << ": " << endl;
      fHistFile->cd();
      TH2D *h2 = new TH2D(Form("acc_%s", target[i].c_str()), Form("acceptance (%s)", target[i].c_str()),
			  venergy[target[i]].size(), 0., venergy[target[i]].size(),
			  vthickness[target[i]].size(), 0., vthickness[target[i]].size());
      setTitles(h2, "#mu^{+} Beam Energy", "Target thickness", 0.04, 1.2, 2.3);

      TH2D *h3 = new TH2D(Form("muprod_%s", target[i].c_str()), Form("Mu production (%s)", target[i].c_str()),
			  venergy[target[i]].size(), 0., venergy[target[i]].size(),
			  vthickness[target[i]].size(), 0., vthickness[target[i]].size());
      setTitles(h3, "#mu^{+} Beam Energy", "Target thickness", 0.04, 1.2, 2.3);


      for (unsigned j = 0; j < venergy[target[i]].size(); ++j) {
	h2->GetXaxis()->SetBinLabel(j+1, venergy[target[i]].at(j).c_str());
	for (unsigned k = 0; k < vthickness[target[i]].size(); ++k) {
	  h2->GetYaxis()->SetBinLabel(k+1, vthickness[target[i]].at(k).c_str());
	  filename = target[i] + "-"
	    + venergy[target[i]].at(j)
	    + "-" + vthickness[target[i]].at(k)
	    + ".default.root";

	  TFile *f = TFile::Open(filename.c_str());
	  TH1D *h1 = (TH1D*)f->Get("acc");
	  TH1D *hp = (TH1D*)f->Get("muprod");
	  TH1D *hm = (TH1D*)f->Get("nmuons");
	  double nacc = h1->GetBinContent(h1->FindBin(0.01));
	  double npro = hp->GetBinContent(hp->FindBin(0.01));
	  double nmuons = totalMuons(hm);
	  double acc = nacc/nmuons;
	  double acce= dEff(static_cast<int>(nacc), static_cast<int>(nmuons));
	  h2->SetBinContent(j+1, k+1, acc);
	  h2->SetBinError(j+1, k+1, acce);

	  double epro = npro/nmuons;
	  double eproe= dEff(static_cast<int>(npro), static_cast<int>(nmuons));
	  h3->SetBinContent(j+1, k+1, epro);
	  h3->SetBinError(j+1, k+1, eproe);

	  cout << filename << " nacc: " << nacc << " nmuons: " << nmuons
	       << " acc: " << acc << " +/- " << acce
	       << " epro: " << epro << " +/- " << eproe
	       << endl;
	  f->Close();
	}
      }

      h2->Draw("texte");
      savePad(Form("acc_%s.pdf", target[i].c_str()));

      h3->Draw("texte");
      savePad(Form("muprod_%s.pdf", target[i].c_str()));

    }
  }


}

// ----------------------------------------------------------------------
double plotResults::totalMuons(TH1 *h) {
  double n(0.);
  for (int i = 0; i <= h->GetNbinsX(); ++i) {
    n += h->GetBinContent(i)*h->GetBinLowEdge(i);
  }
  return n;
}

// ----------------------------------------------------------------------
void plotResults::loopOverTree(TTree *t, int ifunc, int nevts, int nstart) {
  int nentries = Int_t(t->GetEntries());
  int nbegin(0), nend(nentries);
  if (nevts > 0 && nentries > nevts) {
    nentries = nevts;
    nbegin = 0;
    nend = nevts;
  }
  if (nevts > 0 && nstart > 0) {
    nentries = nstart + nevts;
    nbegin = nstart;
    if (nstart + nevts < t->GetEntries()) {
      nend = nstart + nevts;
    } else {
      nend = t->GetEntries();
    }
  }

  nentries = nend - nstart;

  int step(1000000);
  if (nentries < 5000000)  step = 500000;
  if (nentries < 1000000)  step = 100000;
  if (nentries < 100000)   step = 10000;
  if (nentries < 10000)    step = 1000;
  if (nentries < 1000)     step = 100;
  cout << "==> plotResults::loopOverTree> loop over dataset " << fCds->fName << " in file "
       << t->GetDirectory()->GetName()
       << " with " << nentries << " entries"
       << " nbegin = " << nbegin << " nend = " << nend
       << endl;

  // -- setup loopfunction through pointer to member functions
  //    (this is the reason why this function is NOT in plotClass!)
  void (plotResults::*pF)(void);
  if (ifunc == 1) pF = &plotResults::loopFunction1;
  if (ifunc == 2) pF = &plotResults::loopFunction2;

  // ----------------------------
  // -- the real loop starts here
  // ----------------------------
  for (int jentry = nbegin; jentry < nend; jentry++) {
    t->GetEntry(jentry);
    if (jentry%step == 0) {
      TTimeStamp ts;
      cout << Form(" .. evt = %d", jentry)
           << ", time now: " << ts.AsString("lc")
           << endl;
    }

    candAnalysis();
    (this->*pF)();
  }
}

// ----------------------------------------------------------------------
void plotResults::loopFunction1() {

}


// ----------------------------------------------------------------------
void plotResults::loopFunction2() {

}



// ----------------------------------------------------------------------
void plotResults::loadFiles(string afiles) {

  string files = fDirectory + string("/") + afiles;
  cout << "==> plotResults::loadFile loading files listed in " << files << endl;

  char buffer[1000];
  ifstream is(files.c_str());
  while (is.getline(buffer, 1000, '\n')) {
    if (buffer[0] == '#') {continue;}
    if (buffer[0] == '/') {continue;}

    string sbuffer = string(buffer);
    replaceAll(sbuffer, " ", "");
    replaceAll(sbuffer, "\t", "");
    if (sbuffer.size() < 1) continue;

    string::size_type m1 = sbuffer.find("lumi=");
    string stype = sbuffer.substr(5, m1-5);

    string::size_type m2 = sbuffer.find("file=");
    string slumi = sbuffer.substr(m1+5, m2-m1-6);
    string sfile = sbuffer.substr(m2+5);
    string sname("nada"), sdecay("nada");

    TFile *pF(0);
    dataset *ds(0);

    if (string::npos != stype.find("Data")) {
      // -- Charmonium
      pF = loadFile(sfile);
      ds = new dataset();
      ds->fSize = 1.2;
      ds->fWidth = 2;
      if (string::npos != stype.find("blablabla")) {
        sname = "bmmCharmonium";
        sdecay = "bmm";
        ds->fColor = kBlack;
        ds->fSymbol = 20;
        ds->fF      = pF;
        ds->fBf     = 1.;
        ds->fMass   = 1.;
        ds->fFillStyle = 3365;
      }
    }
  }
}
