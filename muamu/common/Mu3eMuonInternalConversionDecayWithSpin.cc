/*
 * Mu3eMuonInternalConversionDecayWithSpin.cpp
 *
 *  Created on: Jun 1, 2012
 *      Author: nberger
 */
#include "Mu3eMuonInternalConversionDecayWithSpin.hh"

#include "G4ParticleDefinition.hh"
#include "G4DecayProducts.hh"
#include "G4VDecayChannel.hh"
#include "G4LorentzVector.hh"

#include "TH1D.h"
#include "TH2D.h"

#include "rambo.hh"

#define MMU 0.1057
#define ME  0.0005


Mu3eMuonInternalConversionDecayWithSpin::Mu3eMuonInternalConversionDecayWithSpin(
										 const G4String& theName               ///< Name of the decay channel
										 , const G4String& theParentName         ///< Name of the parent particle
										 , G4double        theBR                 ///< Branching fraction of the decay
										 , G4int           Selection             ///< Type of selection; 0 cut on all three electrons, 1 cut on an e+e- pair
										 , G4double        CosThetaCut           ///< Cut on the absolute value of the polar angle of the decay electrons, 0 for no cut
										 , G4double        EminCut               ///< Cut on the minimum energy of the decay electrons, 0 for no cut
										 , G4double        EvisCut               ///< Cut on the minimum energy visible in the three electrons, 0 for no cut
										 , G4double        MassCut               ///< Cut on the minimum invariant mass in the three electrons, 0 for no cut
										 , G4double        theParentPolarization ///< Polarization of the parent particle
										 , G4bool          _weighted
										 )
: G4VDecayChannel(theName)
  , selectionType(Selection)
  , fCosThetaCut(CosThetaCut?CosThetaCut:1e6)
  , fEminCut(EminCut)
  , fEvisCut(EvisCut)
  , fMassCut(MassCut)
  , weighted(_weighted)
{

  if(theParentPolarization > 0.0)
    {
      SetPolarization({0.0,0.0,-1.0});
    }
  else
    {
      SetPolarization({0.0,0.0,0.0});
    }
  G4cout << "Polarization of the muon is " << parent_polarization << G4endl;
  SetBR(theBR);
  SetParent(theParentName);

  SetNumberOfDaughters(5);
  SetDaughter(0, "e+");
  SetDaughter(1, "e-");
  SetDaughter(2, "e+");
  SetDaughter(3, "anti_nu_mu");
  SetDaughter(4, "nu_e");

  //G4cout << "CosThetaCut " << fCosThetaCut << G4endl;
  if(fCosThetaCut != 1e6 || fEminCut != 0 || fEvisCut !=0 || fMassCut !=0){
    //Determine the maximal matrix element for the given cuts
    G4cout << "Normalizing matrix element" << G4endl;
    G4cout << "For: Theta Cut " << fCosThetaCut << G4endl;
    G4cout << "     EminCut   " << fEminCut << G4endl;
    G4cout << "     EvisCut   " << fEvisCut << G4endl;
    G4cout << "     MassCut   " << fMassCut << G4endl;
    matmax =0;



    int N = 5;
    double ET = MMU;
    double XME = ME;
    double XM[] = {XME, XME, XME, 0, 0};
    double P[5][4];
    double WT = 0;
    int LW = 0;
    double ew, mat;
    double pt1 = 0.;  // transverse momentum
    double pt2 = 0.;
    double pt3 = 0.;

    double allsum = 0;
    double ptselsum = 0;
    double selsum = 0;

    double allsumcomp = 0;
    double ptselsumcomp = 0;
    double selsumcomp = 0;

    double allsumsimple = 0;
    double ptselsumsimple = 0;
    double selsumsimple = 0;

    for(unsigned int i=0; i < 1e6; i++){
      rambo(N,ET,XM,P,&WT,LW);
      //G4cout << "WT: " <<WT << G4endl;
      //if(selection(P) <0)
      //	continue;
      mat = matmu3e2nu(P);
      //G4cout << "mat: " << mat << G4endl;
      ew = WT*mat;
      pt1 = sqrt(P[0][0]*P[0][0] + P[0][1]*P[0][1]);
      pt2 = sqrt(P[1][0]*P[1][0] + P[1][1]*P[1][1]);
      pt3 = sqrt(P[2][0]*P[2][0] + P[2][1]*P[2][1]);

      allsumsimple += ew;

      double y = ew - allsumcomp;
      double t = allsum + y;
      allsumcomp = (t-allsum) -y;
      allsum = t;

      if (pt1 > 17.0 && pt2 > 17.0 && pt3 > 17.0)    // SINDRUM had pT acceptance > 17.0MeV
	{
	  ptselsumsimple += ew;
	  y = ew - ptselsumcomp;
	  t = ptselsum + y;
	  ptselsumcomp = (t-ptselsum) -y;
	  ptselsum = t;
	}

      if(selection(P) <0)
	continue;

      if(ew > matmax)
	matmax = ew;

      selsumsimple += ew;

      //selsum += ew;
      y = ew - selsumcomp;
      t = selsum + y;
      selsumcomp = (t-selsum) -y;
      selsum = t;
    }
    G4cout << "Maximum matrix element is " << std::scientific << matmax << G4endl;
    G4cout << "Branching fraction for selection " << selsum/allsum << G4endl;
    G4cout << "Branching fraction for selection (simple sum) " << selsumsimple/allsumsimple << G4endl;
    G4cout << "Relative branching fraction compared to SINDURM " << selsum/ptselsum << G4endl;
    G4cout << "Relative branching fraction compared to SINDURM (simple sum) " << selsumsimple/ptselsumsimple << G4endl;
  } else {
    matmax =1.5; //maximum found was 2.085, but that is extremely slow
  }

}

double Mu3eMuonInternalConversionDecayWithSpin::GetBranchingFraction(double nev, TH1D * allhisto, TH1D * selhisto, TH1D * espectrum, TH1D *eminspectrum, TH2D * eespectrum, TH2D * eespectrumsel, TH2D * dalitz, TH1D * etothisto, TH1D * costhposhisto, TH1D * costhelhhisto, TH1D * costhellhisto, TH1D * costhposhistosel, TH1D * costhelhhistosel, TH1D * costhellhistosel, TH1D * ptposhistosel, TH1D * ptelhhistosel, TH1D * ptellhistosel){


  int N = 5;
  double ET = MMU;
  double XME = ME;
  double XM[] = {XME, XME, XME, 0, 0};
  double P[5][4];
  double WT = 0;
  int LW = 0;
  double ew, mat;

  double allsum = 0;
  double selsum = 0;

  double mass;
  double etot;

  double costhpos, costhelh, costhell;
  double ptcut = 10.0;
  double eviscut = 100.0;

  double ptpos, ptelh, ptell;
  double evis;

  int sel;

  for(unsigned long long int i=0; i < nev; i++){
    if(i%1000000==0)
      G4cout << "At event " << i << G4endl;
    rambo(N,ET,XM,P,&WT,LW);
    //G4cout << "WT: " <<WT << G4endl;
    //if(selection(P) <0)
    //	continue;
    mat = matmu3e2nu(P);
    //G4cout << "mat: " << mat << G4endl;
    ew = WT*mat;
    //G4cout << "ew: " <<ew << G4endl;

    allsum += ew;

    //		sel = selection(P,mass);
    sel = selection(P, mass, costhpos, costhelh, costhell, ptcut, eviscut, ptelh, ptell, ptpos, evis);
    etot = P[0][3] + P[1][3] + P[2][3];

    allhisto->Fill(mass,ew);
    etothisto->Fill(etot, ew);
    espectrum->Fill(P[0][3]);
    espectrum->Fill(P[1][3]);
    espectrum->Fill(P[2][3]);

    costhposhisto->Fill(costhpos, ew);
    costhelhhisto->Fill(costhelh, ew);
    costhellhisto->Fill(costhell, ew);

    double px01=P[0][0] + P[1][0];
    double py01=P[0][1] + P[1][1];
    double pz01=P[0][2] + P[1][2];
    double etot01=P[0][3] + P[1][3];
    double m01=sqrt(etot01*etot01-(px01*px01+py01*py01+pz01*pz01));

    double px12=P[1][0] + P[2][0];
    double py12=P[1][1] + P[2][1];
    double pz12=P[1][2] + P[2][2];
    double etot12=P[1][3] + P[2][3];
    double m12=sqrt(etot12*etot12-(px12*px12+py12*py12+pz12*pz12));

    double temp;
    if(m01 < m12){
      temp = m01;
      m01 = m12;
      m12 = temp;
    }


    eespectrum->Fill(m01,m12,ew);


    double emin = P[0][3];
    if(P[1][3] < emin)
      emin = P[1][3];
    if(P[2][3] < emin)
      emin = P[2][3];

    eminspectrum->Fill(emin);

    if(sel <0)
      continue;

    selhisto->Fill(mass,ew);
    eespectrumsel->Fill(m01,m12,ew);
    dalitz->Fill(m01*m01,m12*m12,ew);

    costhposhistosel->Fill(costhpos, ew);
    costhelhhistosel->Fill(costhelh, ew);
    costhellhistosel->Fill(costhell, ew);

    ptposhistosel->Fill(ptpos, ew);
    ptelhhistosel->Fill(ptelh, ew);
    ptellhistosel->Fill(ptell, ew);

    if(ew > matmax)
      matmax = ew;

    selsum += ew;
  }

  return selsum/allsum;

}

G4DecayProducts* Mu3eMuonInternalConversionDecayWithSpin::DecayIt(G4double) {
  //ul  auto d1 = getEmptyProducts();
  //create parent G4DynamicParticle at rest
  G4ThreeVector dummy;
  G4DynamicParticle * parentparticle = new G4DynamicParticle( G4MT_parent, dummy, 0.0);
  //create G4Decayproducts
  G4DecayProducts *d1 = new G4DecayProducts(*parentparticle);
  delete parentparticle;

  int N = 5;
  double ET = MMU;
  double XME = ME;
  double XM[] = {XME, XME, XME, 0, 0};
  double P[5][4];
  double WT = 0;
  int LW = 0;

  if(!weighted){ // Generate unweighted events via accept/reject

    double rnd = CLHEP::RandFlat::shoot() * matmax; //2.085 was found as the highest weight in a run of 10^8 events for no cuts applied
    //G4cout <<"RND: " << rnd << G4endl;
    double ew = -1;
    double mat;

    while(ew < rnd){
      rambo(N,ET,XM,P,&WT,LW);
      if(selection(P) <0)
	continue;
      mat = matmu3e2nu(P);
      ew = WT*mat;
      //      G4cout << "WT: " << WT << "mat: " << mat << " ew: " << ew << " <? rnd = " << rnd << G4endl;
    }
  }
  else { // generate weighted events
    while(true){
      rambo(N,ET,XM,P,&WT,LW);
      if(selection(P) == 0)
	break;
    }
    double mat = matmu3e2nu(P);
    double ew = WT*mat;

    // Mu3eEvent * ev = Mu3eEvent::GetInstance();
    // ev->ScaleWeight(ew);
  }


  for(int i = 0; i < 5; i++) {
    G4LorentzVector momentum(P[i][0], P[i][1], P[i][2], P[i][3]);
    auto particle = new G4DynamicParticle(G4MT_daughters[i], momentum);
    d1->PushProducts(particle);
  }

  //for(unsigned int i=0; i < 5; i++)
  //	for(unsigned int j=0; j < 4; j++)
  //		G4cout << "P["<<i<<"]["<< j <<"] " << P[i][j] << G4endl;

  //G4cout << "Returning" << G4endl;

  return d1;
}


int Mu3eMuonInternalConversionDecayWithSpin::selection(double P[5][4]) {
  /* particles order: mu+ -> e+ e- e+ mumu nue */
  //double enu,ptot;
  double px,py,pz;//pnorm,plong,denom,box;
  int k;
  //int nl1,nl2;
  int nele=0;
  double costh0,costh1,costh2;
  double etot,xmtot,xm;


  // selection of three visible electrons
  if (selectionType==0) {
    // final cut
    //enu=P[3][3]+P[4][3];
    // G4cout << "Selection type " << selectionType << G4endl;

    costh0=P[0][2]/P[0][3];
    if ((fabs(costh0)>fCosThetaCut)||(P[0][3]<fEminCut)) {
      return -1;
    }
    costh1=P[1][2]/P[1][3];
    if ((fabs(costh1)>fCosThetaCut)||(P[1][3]<fEminCut)) {
      return -1;
    }
    costh2=P[2][2]/P[2][3];
    if ((fabs(costh2)>fCosThetaCut)||(P[2][3]<fEminCut)) {
      return -1;
    }

    //G4cout << "Cos and e passed" << G4endl;

    //  if (P[0][3]>P[2][3])
    //    emax=P[0][3];
    //  else
    //   emax=P[2][3];

    // xm=sqrt(MMU-2*emax);
    // kinematics
    px=P[0][0] + P[1][0] + P[2][0];
    py=P[0][1] + P[1][1] + P[2][1];
    pz=P[0][2] + P[1][2] + P[2][2];
    etot=P[0][3] + P[1][3] + P[2][3];

    //ptot=sqrt(px*px+py*py+pz*pz);
    xm=sqrt(etot*etot-(px*px+py*py+pz*pz));

    if(fMassCut > 0)
      if (xm<fMassCut) return -1;

    // G4cout << "Mass passed" << G4endl;
    return 0;


    // selection of exactly one visible e+e- pair (third is BG)
  } else if (selectionType==1) {

    // final cut
    //enu=P[3][3]+P[4][3];

    // no opposite sign?
    costh1=P[1][2]/P[1][3];
    if ((fabs(costh1)>fCosThetaCut)||(P[1][3]<fEminCut)) {
      return -1;
    }

    // electron accepted
    costh0=P[0][2]/P[0][3];
    if (fabs(costh0)<fCosThetaCut&&P[0][3]>fEminCut) {
      nele++;
      k=0;
    }

    costh2=P[2][2]/P[2][3];
    if (fabs(costh2)<fCosThetaCut&&P[2][3]>fEminCut) {
      nele++;
      k=2;
    }

    // only one electron visible?
    if (nele!=1) return -1;

    // kinematics
    px=P[1][0]+P[k][0];
    py=P[1][1]+P[k][1];
    pz=P[1][2]+P[k][2];
    etot=P[1][3]+P[k][3];

    //ptot=sqrt(px*px+py*py+pz*pz);
    xmtot=sqrt(etot*etot-(px*px+py*py+pz*pz));

    // in allowed energy and momentum interval?
    if(fEvisCut > 0)
      if (fabs(etot-MMU/2.0)>fEvisCut) return -1;

    if(fMassCut > 0)
      if (xmtot<fMassCut) return -1;

    return 0;
  }

  return 0;
}


int Mu3eMuonInternalConversionDecayWithSpin::selection(double P[5][4], double & mass) {
  /* particles order: mu+ -> e+ e- e+ mumu nue */
  //double enu,ptot;
  double px,py,pz;//pnorm,plong,denom,box;
  int k;
  //int nl1,nl2;
  int nele=0;
  double costh0,costh1,costh2;
  double etot,xmtot,xm;


  // selection of three visible electrons
  if (selectionType==0) {
    // final cut
    //enu=P[3][3]+P[4][3];
    // G4cout << "Selection type " << selectionType << G4endl;

    px=P[0][0] + P[1][0] + P[2][0];
    py=P[0][1] + P[1][1] + P[2][1];
    pz=P[0][2] + P[1][2] + P[2][2];
    etot=P[0][3] + P[1][3] + P[2][3];

    //ptot=sqrt(px*px+py*py+pz*pz);
    xm=sqrt(etot*etot-(px*px+py*py+pz*pz));

    mass = xm;


    costh0=P[0][2]/P[0][3];
    if ((fabs(costh0)>fCosThetaCut)||(P[0][3]<fEminCut)) {
      return -1;
    }
    costh1=P[1][2]/P[1][3];
    if ((fabs(costh1)>fCosThetaCut)||(P[1][3]<fEminCut)) {
      return -1;
    }
    costh2=P[2][2]/P[2][3];
    if ((fabs(costh2)>fCosThetaCut)||(P[2][3]<fEminCut)) {
      return -1;
    }

    //G4cout << "Cos and e passed" << G4endl;

    //  if (P[0][3]>P[2][3])
    //    emax=P[0][3];
    //  else
    //   emax=P[2][3];

    // xm=sqrt(MMU-2*emax);
    // kinematics


    if(fMassCut > 0)
      if (xm<fMassCut) return -1;

    // G4cout << "Mass passed" << G4endl;
    return 0;


    // selection of exactly one visible e+e- pair (third is BG)
  } else if (selectionType==1) {


    // final cut
    //enu=P[3][3]+P[4][3];
    mass =0;
    // no opposite sign?
    costh1=P[1][2]/P[1][3];
    if ((fabs(costh1)>fCosThetaCut)||(P[1][3]<fEminCut)) {
      return -1;
    }

    // electron accepted
    costh0=P[0][2]/P[0][3];
    if (fabs(costh0)<fCosThetaCut&&P[0][3]>fEminCut) {
      nele++;
      k=0;
    }

    costh2=P[2][2]/P[2][3];
    if (fabs(costh2)<fCosThetaCut&&P[2][3]>fEminCut) {
      nele++;
      k=2;
    }

    // only one electron visible?
    if (nele!=1) return -1;


    // kinematics
    px=P[1][0]+P[k][0];
    py=P[1][1]+P[k][1];
    pz=P[1][2]+P[k][2];
    etot=P[1][3]+P[k][3];

    //ptot=sqrt(px*px+py*py+pz*pz);
    xmtot=sqrt(etot*etot-(px*px+py*py+pz*pz));

    mass = xmtot;

    // in allowed energy and momentum interval?
    if(fEvisCut > 0)
      if (fabs(etot-MMU/2.0)>fEvisCut) return -1;

    if(fMassCut > 0)
      if (xmtot<fMassCut) return -1;

    return 0;
  }

  return 0;
}

int Mu3eMuonInternalConversionDecayWithSpin::selection(double P[5][4], double & mass, double & costhpos, double & costhelh, double & costhell) {
  /* particles order: mu+ -> e+ e- e+ mumu nue */
  //double enu,ptot;
  double px,py,pz;//pnorm,plong,denom,box;
  int k;
  //int nl1,nl2;
  int nele=0;
  double costh0,costh1,costh2;
  double etot,xmtot,xm;


  // selection of three visible electrons
  if (selectionType==0) {
    // final cut
    //enu=P[3][3]+P[4][3];
    // G4cout << "Selection type " << selectionType << G4endl;

    px=P[0][0] + P[1][0] + P[2][0];
    py=P[0][1] + P[1][1] + P[2][1];
    pz=P[0][2] + P[1][2] + P[2][2];
    etot=P[0][3] + P[1][3] + P[2][3];

    //ptot=sqrt(px*px+py*py+pz*pz);
    xm=sqrt(etot*etot-(px*px+py*py+pz*pz));

    mass = xm;


    if(P[0][3] > P[2][3])
      {
        costhelh = P[0][2]/P[0][3];
        costhell = P[2][2]/P[2][3];
      }
    else
      {
        costhelh = P[2][2]/P[2][3];
        costhell = P[0][2]/P[0][3];
      }
    costh0=P[0][2]/P[0][3];
    if ((fabs(costh0)>fCosThetaCut)||(P[0][3]<fEminCut)) {
      return -1;
    }
    costhpos=P[1][2]/P[1][3];
    costh1 = costhpos;
    if ((fabs(costh1)>fCosThetaCut)||(P[1][3]<fEminCut)) {
      return -1;
    }
    costh2=P[2][2]/P[2][3];
    if ((fabs(costh2)>fCosThetaCut)||(P[2][3]<fEminCut)) {
      return -1;
    }

    //G4cout << "Cos and e passed" << G4endl;

    //  if (P[0][3]>P[2][3])
    //    emax=P[0][3];
    //  else
    //   emax=P[2][3];

    // xm=sqrt(MMU-2*emax);
    // kinematics


    if(fMassCut > 0)
      if (xm<fMassCut) return -1;

    // G4cout << "Mass passed" << G4endl;
    return 0;


    // selection of exactly one visible e+e- pair (third is BG)
  } else if (selectionType==1) {


    // final cut
    //enu=P[3][3]+P[4][3];
    mass =0;
    // no opposite sign?
    costh1=P[1][2]/P[1][3];
    if ((fabs(costh1)>fCosThetaCut)||(P[1][3]<fEminCut)) {
      return -1;
    }

    // electron accepted
    costh0=P[0][2]/P[0][3];
    if (fabs(costh0)<fCosThetaCut&&P[0][3]>fEminCut) {
      nele++;
      k=0;
    }

    costh2=P[2][2]/P[2][3];
    if (fabs(costh2)<fCosThetaCut&&P[2][3]>fEminCut) {
      nele++;
      k=2;
    }

    // only one electron visible?
    if (nele!=1) return -1;


    // kinematics
    px=P[1][0]+P[k][0];
    py=P[1][1]+P[k][1];
    pz=P[1][2]+P[k][2];
    etot=P[1][3]+P[k][3];

    //ptot=sqrt(px*px+py*py+pz*pz);
    xmtot=sqrt(etot*etot-(px*px+py*py+pz*pz));

    mass = xmtot;

    // in allowed energy and momentum interval?
    if(fEvisCut > 0)
      if (fabs(etot-MMU/2.0)>fEvisCut) return -1;

    if(fMassCut > 0)
      if (xmtot<fMassCut) return -1;

    return 0;
  }

  return 0;
}

int Mu3eMuonInternalConversionDecayWithSpin::selection(double P[5][4], double & mass, double & costhpos, double & costhelh, double & costhell, double ptcut, double eviscut, double & ptelh, double & ptell, double & ptpos, double & evis) {
  /* particles order: mu+ -> e+ e- e+ mumu nue */
  //double enu,ptot;
  double px,py,pz;//pnorm,plong,denom,box;
  int k;
  //int nl1,nl2;
  int nele=0;
  double costh0,costh1,costh2;
  double etot,xmtot,xm;
  double pt0, pt1, pt2;
  //double evis;


  // selection of three visible electrons
  if (selectionType==0) {
    // final cut
    //enu=P[3][3]+P[4][3];
    // G4cout << "Selection type " << selectionType << G4endl;

    px=P[0][0] + P[1][0] + P[2][0];
    py=P[0][1] + P[1][1] + P[2][1];
    pz=P[0][2] + P[1][2] + P[2][2];
    etot=P[0][3] + P[1][3] + P[2][3];

    //ptot=sqrt(px*px+py*py+pz*pz);
    xm=sqrt(etot*etot-(px*px+py*py+pz*pz));

    mass = xm;
    evis = xm;


    if(P[0][3] > P[2][3])
      {
        costhelh = P[0][2]/P[0][3];
        costhell = P[2][2]/P[2][3];
        ptelh = sqrt(P[0][0]*P[0][0]+P[0][1]*P[0][1]);
        ptell = sqrt(P[2][0]*P[2][0]+P[2][1]*P[2][1]);
      }
    else
      {
        costhelh = P[2][2]/P[2][3];
        costhell = P[0][2]/P[0][3];
        ptelh = sqrt(P[2][0]*P[2][0]+P[2][1]*P[2][1]);
        ptell = sqrt(P[0][0]*P[0][0]+P[0][1]*P[0][1]);
      }
    costh0=P[0][2]/P[0][3];
    pt0 = sqrt(P[0][0]*P[0][0]+P[0][1]*P[0][1]);
    //    if ((fabs(costh0)>fCosThetaCut)||(P[0][3]<fEminCut)) {
    //      return -1;
    //    }
    costhpos=P[1][2]/P[1][3];
    pt1 = sqrt(P[1][0]*P[1][0]+P[1][1]*P[1][1]);
    ptpos = pt1;
    costh1 = costhpos;
    //    if ((fabs(costh1)>fCosThetaCut)||(P[1][3]<fEminCut)) {
    //      return -1;
    //    }
    costh2=P[2][2]/P[2][3];
    pt2 = sqrt(P[2][0]*P[2][0]+P[2][1]*P[2][1]);
    //    if ((fabs(costh2)>fCosThetaCut)||(P[2][3]<fEminCut)) {
    //      return -1;
    //    }

    if((fabs(costh0)>fCosThetaCut) || (fabs(costh1)>fCosThetaCut) || (fabs(costh2)>fCosThetaCut) || pt0 < ptcut || pt1 < ptcut || pt2 < ptcut || evis < eviscut)
      {
        return -1;
      }

    //G4cout << "Cos and e passed" << G4endl;

    //  if (P[0][3]>P[2][3])
    //    emax=P[0][3];
    //  else
    //   emax=P[2][3];

    // xm=sqrt(MMU-2*emax);
    // kinematics


    //   if(fMassCut > 0)
    //        if (xm<fMassCut) return -1;

    // G4cout << "Mass passed" << G4endl;
    return 0;


    // selection of exactly one visible e+e- pair (third is BG)
  } else if (selectionType==1) {


    // final cut
    //enu=P[3][3]+P[4][3];
    mass =0;
    // no opposite sign?
    costh1=P[1][2]/P[1][3];
    if ((fabs(costh1)>fCosThetaCut)||(P[1][3]<fEminCut)) {
      return -1;
    }

    // electron accepted
    costh0=P[0][2]/P[0][3];
    if (fabs(costh0)<fCosThetaCut&&P[0][3]>fEminCut) {
      nele++;
      k=0;
    }

    costh2=P[2][2]/P[2][3];
    if (fabs(costh2)<fCosThetaCut&&P[2][3]>fEminCut) {
      nele++;
      k=2;
    }

    // only one electron visible?
    if (nele!=1) return -1;


    // kinematics
    px=P[1][0]+P[k][0];
    py=P[1][1]+P[k][1];
    pz=P[1][2]+P[k][2];
    etot=P[1][3]+P[k][3];

    //ptot=sqrt(px*px+py*py+pz*pz);
    xmtot=sqrt(etot*etot-(px*px+py*py+pz*pz));

    mass = xmtot;

    // in allowed energy and momentum interval?
    if(fEvisCut > 0)
      if (fabs(etot-MMU/2.0)>fEvisCut) return -1;

    if(fMassCut > 0)
      if (xmtot<fMassCut) return -1;

    return 0;
  }

  return 0;
}

double Mu3eMuonInternalConversionDecayWithSpin::matmu3e2nu(double P[5][4]) {            // polarized matrix element by A. Signer et al. (PSI)

  /* order of FS particles:
     e+,e-,e+,numu,nue
  */

  // Get polarization of muon
  G4ThreeVector n1 = parent_polarization;                       // polarization of the muon

  if (GetParentName() == "mu+")                                 // for mu+: n1 -> -n1
    {
      n1 = -n1;
    }

  // Calculate matrix element
  double s12, s13, s14, s15, s16, s23, s24, s25, s26, s34, s35, s36, s45, s46, s56, s2n, s3n, s4n, s5n, s6n;
  double Ip11, Ip22, Ip33, Ip44, Ip12r, Ip13r, Ip14r, Ip23r, Ip24r, Ip34r;
  double den1, den2, den3, den4;
  double MatElSquared;

  //double GF = (1.1663787e-5)/(1*GeV*GeV)*hbarc*hbarc*hbarc;       // Fermi coupling const. GF/(hbar c)^3 = 1.17e-5/GeV^2
  double hquerc = CLHEP::hbarc / (CLHEP::MeV * CLHEP::fermi);       // hbar c = 197 MeV fm
  double GF = (1.1663787e-5)/(1.e6)*hquerc*hquerc*hquerc;           // Fermi coupling const. GF/(hbar c)^3 = 1.17e-5/GeV^2
  double ALPHA = CLHEP::fine_structure_const;

  // scalar products
  s13 = 2.0*MMU*P[0][3];
  s15 = 2.0*MMU*P[1][3];
  s16 = 2.0*MMU*P[2][3];
  s12 = 2.0*MMU*P[3][3];
  s14 = 2.0*MMU*P[4][3];

  s35= 2.0*(P[0][3]*P[1][3]-P[0][0]*P[1][0]-P[0][1]*P[1][1]-P[0][2]*P[1][2]);
  s36= 2.0*(P[0][3]*P[2][3]-P[0][0]*P[2][0]-P[0][1]*P[2][1]-P[0][2]*P[2][2]);
  s23= 2.0*(P[0][3]*P[3][3]-P[0][0]*P[3][0]-P[0][1]*P[3][1]-P[0][2]*P[3][2]);
  s34= 2.0*(P[0][3]*P[4][3]-P[0][0]*P[4][0]-P[0][1]*P[4][1]-P[0][2]*P[4][2]);

  s56= 2.0*(P[1][3]*P[2][3]-P[1][0]*P[2][0]-P[1][1]*P[2][1]-P[1][2]*P[2][2]);
  s25= 2.0*(P[1][3]*P[3][3]-P[1][0]*P[3][0]-P[1][1]*P[3][1]-P[1][2]*P[3][2]);
  s45= 2.0*(P[1][3]*P[4][3]-P[1][0]*P[4][0]-P[1][1]*P[4][1]-P[1][2]*P[4][2]);

  s26= 2.0*(P[2][3]*P[3][3]-P[2][0]*P[3][0]-P[2][1]*P[3][1]-P[2][2]*P[3][2]);
  s46= 2.0*(P[2][3]*P[4][3]-P[2][0]*P[4][0]-P[2][1]*P[4][1]-P[2][2]*P[4][2]);

  s24= 2.0*(P[3][3]*P[4][3]-P[3][0]*P[4][0]-P[3][1]*P[4][1]-P[3][2]*P[4][2]);

  s2n= -2.0*(n1[0]*P[3][0]+n1[1]*P[3][1]+n1[2]*P[3][2]);
  s3n= -2.0*(n1[0]*P[0][0]+n1[1]*P[0][1]+n1[2]*P[0][2]);
  s4n= -2.0*(n1[0]*P[4][0]+n1[1]*P[4][1]+n1[2]*P[4][2]);
  s5n= -2.0*(n1[0]*P[1][0]+n1[1]*P[1][1]+n1[2]*P[1][2]);
  s6n= -2.0*(n1[0]*P[2][0]+n1[1]*P[2][1]+n1[2]*P[2][2]);

  // diagonal elements
  Ip11 = -32.* ME*ME*ME*ME* s14* s23 - 32.* ME*ME* MMU*MMU* s14* s23 - 16.* ME*ME* s14* s15* s23 - 16.* ME*ME* s14* s16* s23 + 16.* s14* s15* s16* s23 + 32.* ME*ME* MMU*MMU* s23* s45 + 32.* ME*ME* s15* s23* s45 + 16.* ME*ME* s16* s23* s45 - 8.* s15* s16* s23* s45 + 8.* s16*s16* s23* s45 + 32.* ME*ME* MMU*MMU* s23* s46 + 16.* ME*ME* s15* s23* s46 + 8.* s15*s15* s23* s46 + 32.* ME*ME* s16* s23* s46 - 8.* s15* s16* s23* s46 + 32.* ME*ME*ME*ME* MMU* s23* s4n + 32.* ME*ME* MMU*MMU*MMU* s23* s4n - 16.* MMU* s15* s16* s23* s4n - 16.* ME*ME* s14* s23* s56 - 16.* MMU*MMU* s14* s23* s56 - 8.* s14* s15* s23* s56 - 8.* s14* s16* s23* s56 + 16.* MMU*MMU* s23* s45* s56 + 8.* s15* s23* s45* s56 + 16.* MMU*MMU* s23* s46* s56 + 8.* s16* s23* s46* s56 + 16.* ME*ME* MMU* s23* s4n* s56 + 16.* MMU*MMU*MMU* s23* s4n* s56 + 16.* ME*ME* MMU* s14* s23* s5n - 32.* ME*ME* MMU* s23* s45* s5n - 16.* ME*ME* MMU* s23* s46* s5n - 8.* MMU* s15* s23* s46* s5n + 8.* MMU* s16* s23* s46* s5n + 8.* MMU* s14* s23* s56* s5n - 8.* MMU* s23* s45* s56* s5n + 16.* ME*ME* MMU* s14* s23* s6n - 16.* ME*ME* MMU* s23* s45* s6n + 8.* MMU* s15* s23* s45* s6n - 8.* MMU* s16* s23* s45* s6n - 32.* ME*ME* MMU* s23* s46* s6n + 8.* MMU* s14* s23* s56* s6n - 8.* MMU* s23* s46* s56* s6n;

  Ip22 = -64.* ME*ME*ME*ME* s14* s23 - 32.* ME*ME*ME*ME* s14* s25 - 32.* ME*ME*ME*ME* s14* s26 + 16.* ME*ME* s14* s23* s35 + 32.* ME*ME* s14* s25* s35 + 16.* ME*ME* s14* s26* s35 - 8.* s14* s26* s35*s35 + 16.* ME*ME* s14* s23* s36 + 16.* ME*ME* s14* s25* s36 + 32.* ME*ME* s14* s26* s36 + 16.* s14* s23* s35* s36 + 8.* s14* s25* s35* s36 + 8.* s14* s26* s35* s36 - 8.* s14* s25* s36*s36 + 64.* ME*ME*ME*ME* MMU* s23* s4n + 32.* ME*ME*ME*ME* MMU* s25* s4n + 32.* ME*ME*ME*ME* MMU* s26* s4n - 16.* ME*ME* MMU* s23* s35* s4n - 32.* ME*ME* MMU* s25* s35* s4n - 16.* ME*ME* MMU* s26* s35* s4n + 8.* MMU* s26* s35*s35* s4n - 16.* ME*ME* MMU* s23* s36* s4n - 16.* ME*ME* MMU* s25* s36* s4n - 32.* ME*ME* MMU* s26* s36* s4n - 16.* MMU* s23* s35* s36* s4n - 8.* MMU* s25* s35* s36* s4n - 8.* MMU* s26* s35* s36* s4n + 8.* MMU* s25* s36*s36* s4n - 32.* ME*ME* s14* s23* s56 - 16.* ME*ME* s14* s25* s56 - 16.* ME*ME* s14* s26* s56 + 8.* s14* s23* s35* s56 + 8.* s14* s25* s35* s56 + 8.* s14* s23* s36* s56 + 8.* s14* s26* s36* s56 + 32.* ME*ME* MMU* s23* s4n* s56 + 16.* ME*ME* MMU* s25* s4n* s56 + 16.* ME*ME* MMU* s26* s4n* s56 - 8.* MMU* s23* s35* s4n* s56 - 8.* MMU* s25* s35* s4n* s56 - 8.* MMU* s23* s36* s4n* s56 - 8.* MMU* s26* s36* s4n* s56;

  Ip33 = -32.* ME*ME*ME*ME* s14* s26 - 32.* ME*ME* MMU* MMU* s14* s26 - 16.* ME*ME* s13* s14* s26 - 16.* ME*ME* s14* s15* s26 + 16.* s13* s14* s15* s26 + 32.* ME*ME* MMU*MMU* s26* s34 + 32.* ME*ME* s13* s26* s34 + 16.* ME*ME* s15* s26* s34 - 8.* s13* s15* s26* s34 + 8.* s15*s15* s26* s34 - 16.* ME*ME* s14* s26* s35 - 16.* MMU*MMU* s14* s26* s35 - 8.* s13* s14* s26* s35 - 8.* s14* s15* s26* s35 + 16.* MMU*MMU* s26* s34* s35 + 8.* s13* s26* s34* s35 + 16.* ME*ME* MMU* s14* s26* s3n - 32.* ME*ME* MMU* s26* s34* s3n + 8.* MMU* s14* s26* s35* s3n - 8.* MMU* s26* s34* s35* s3n + 32.* ME*ME* MMU*MMU* s26* s45 + 16.* ME*ME* s13* s26* s45 + 8.* s13*s13* s26* s45 + 32.* ME*ME* s15* s26* s45 - 8.* s13* s15* s26* s45 + 16.* MMU*MMU* s26* s35* s45 + 8.* s15* s26* s35* s45 - 16.* ME*ME* MMU* s26* s3n* s45 - 8.* MMU* s13* s26* s3n* s45 + 8.* MMU* s15* s26* s3n* s45 + 32.* ME*ME*ME*ME* MMU* s26* s4n + 32.* ME*ME* MMU*MMU*MMU* s26* s4n - 16.* MMU* s13* s15* s26* s4n + 16.* ME*ME* MMU* s26* s35* s4n + 16.* MMU*MMU*MMU* s26* s35* s4n + 16.* ME*ME* MMU* s14* s26* s5n - 16.* ME*ME* MMU* s26* s34* s5n + 8.* MMU* s13* s26* s34* s5n - 8.* MMU* s15* s26* s34* s5n + 8.* MMU* s14* s26* s35* s5n - 32.* ME*ME* MMU* s26* s45* s5n - 8.* MMU* s26* s35* s45* s5n;

  Ip44 = -32.* ME*ME*ME*ME* s14* s23 - 32.* ME*ME*ME*ME* s14* s25 - 64.* ME*ME*ME*ME* s14* s26 - 16.* ME*ME* s14* s23* s35 - 16.* ME*ME* s14* s25* s35 - 32.* ME*ME* s14* s26* s35 + 32.* ME*ME* s14* s23* s36 + 16.* ME*ME* s14* s25* s36 + 16.* ME*ME* s14* s26* s36 + 8.* s14* s23* s35* s36 + 8.* s14* s26* s35* s36 - 8.* s14* s25* s36*s36 + 32.* ME*ME*ME*ME* MMU* s23* s4n + 32.* ME*ME*ME*ME* MMU* s25* s4n + 64.* ME*ME*ME*ME* MMU* s26* s4n + 16.* ME*ME* MMU* s23* s35* s4n + 16.* ME*ME* MMU* s25* s35* s4n + 32.* ME*ME* MMU* s26* s35* s4n - 32.* ME*ME* MMU* s23* s36* s4n - 16.* ME*ME* MMU* s25* s36* s4n - 16.* ME*ME* MMU* s26* s36* s4n - 8.* MMU* s23* s35* s36* s4n - 8.* MMU* s26* s35* s36* s4n + 8.* MMU* s25* s36*s36* s4n +16.* ME*ME* s14* s23* s56 + 32.* ME*ME* s14* s25* s56 + 16.* ME*ME* s14* s26* s56 + 8.* s14* s25* s35* s56 + 8.* s14* s26* s35* s56 + 8.* s14* s23* s36* s56 + 8.* s14* s25* s36* s56 + 16.* s14* s26* s36* s56 - 16.* ME*ME* MMU* s23* s4n* s56 - 32.* ME*ME* MMU* s25* s4n* s56 - 16.* ME*ME* MMU* s26* s4n* s56 - 8.* MMU* s25* s35* s4n* s56 - 8.* MMU* s26* s35* s4n* s56 - 8.* MMU* s23* s36* s4n* s56 - 8.* MMU* s25* s36* s4n* s56 - 16.* MMU* s26* s36* s4n* s56 - 8.* s14* s23* s56*s56 + 8.* MMU* s23* s4n* s56*s56;

  // off-diagonal elements
  Ip12r = -16.* ME*ME* s13* s14* s23 + 32.* ME*ME*ME*ME* s13* s24 - 8.* ME*ME* s13* s14* s25 - 8.* ME*ME* s13* s14* s26 - 32.* ME*ME*ME*ME* s12* s34 - 8.* ME*ME* s15* s23* s34 - 8.* ME*ME* s16* s23* s34 + 8.* ME*ME* s15* s25* s34 + 8.* ME*ME* s16* s26* s34 + 32.* ME*ME*ME*ME* MMU* s2n* s34 + 8.* ME*ME* s12* s14* s35 + 8.* s14* s16* s23* s35 - 8.* ME*ME* s15* s24* s35 - 4.* s14* s15* s26* s35 + 4.* s14* s16* s26* s35 + 8.* ME*ME* s12* s14* s36 + 8.* s14* s15* s23* s36 - 8.* ME*ME* s16* s24* s36 + 4.* s14* s15* s25* s36 - 4.* s14* s16* s25* s36 - 32.* ME*ME*ME*ME* MMU* s24* s3n + 8.* ME*ME* s13* s23* s45 - 8.* ME*ME* s13* s25* s45 + 8.* ME*ME* s12* s35* s45 - 4.* s16* s23* s35* s45 - 8.* s16* s26* s35* s45 - 8.* ME*ME* MMU* s2n* s35* s45 + 4.* s16* s23* s36* s45 + 8.* s16* s25* s36* s45 - 8.* ME*ME* MMU* s23* s3n* s45 + 8.* ME*ME* MMU* s25* s3n* s45 + 8.* ME*ME* s13* s23* s46 - 8.* ME*ME* s13* s26* s46 + 4.* s15* s23* s35* s46 + 8.* s15* s26* s35* s46 + 8.* ME*ME* s12* s36* s46 - 4.* s15* s23* s36* s46 - 8.* s15* s25* s36* s46 - 8.* ME*ME* MMU* s2n* s36* s46 - 8.* ME*ME* MMU* s23* s3n* s46 + 8.* ME*ME* MMU* s26* s3n* s46 + 16.* ME*ME* MMU* s13* s23* s4n + 8.* ME*ME* MMU* s13* s25* s4n + 8.* ME*ME* MMU* s13* s26* s4n - 8.* ME*ME* MMU* s12* s35* s4n - 8.* MMU* s16* s23* s35* s4n + 4.* MMU* s15* s26* s35* s4n - 4.* MMU* s16* s26* s35* s4n - 8.* ME*ME* MMU* s12* s36* s4n - 8.* MMU* s15* s23* s36* s4n - 4.* MMU* s15* s25* s36* s4n + 4.* MMU* s16* s25* s36* s4n - 8.* s13* s14* s23* s56 + 16.* ME*ME* s13* s24* s56 - 4.* s13* s14* s25* s56 - 4.* s13* s14* s26* s56 - 16.* ME*ME* s12* s34* s56 - 4.* s15* s23* s34* s56 - 4.* s16* s23* s34* s56 - 4.* s16* s25* s34* s56 - 4.* s15* s26* s34* s56 + 16.* ME*ME* MMU* s2n* s34* s56 + 4.* s12* s14* s35* s56 + 4.* s16* s24* s35* s56 + 4.* s12* s14* s36* s56 + 4.* s15* s24* s36* s56 - 16.* ME*ME* MMU* s24* s3n* s56 + 4.* s13* s23* s45* s56 + 4.* s13* s26* s45* s56 - 4.* s12* s36* s45* s56 + 4.* MMU* s2n* s36* s45* s56 - 4.* MMU* s23* s3n* s45* s56 - 4.* MMU* s26* s3n* s45* s56 + 4.* s13* s23* s46* s56 + 4.* s13* s25* s46* s56 - 4.* s12* s35* s46* s56 + 4.* MMU* s2n* s35* s46* s56 - 4.* MMU* s23* s3n* s46* s56 - 4.* MMU* s25* s3n* s46* s56 + 8.* MMU* s13* s23* s4n* s56 + 4.* MMU* s13* s25* s4n* s56 + 4.* MMU* s13* s26* s4n* s56 - 4.* MMU* s12* s35* s4n* s56 - 4.* MMU* s12* s36* s4n* s56 + 8.* ME*ME* MMU* s23* s34* s5n - 8.* ME*ME* MMU* s25* s34* s5n + 8.* ME*ME* MMU* s24* s35* s5n - 4.* MMU* s23* s35* s46* s5n - 8.* MMU* s26* s35* s46* s5n + 4.* MMU* s23* s36* s46* s5n + 8.* MMU* s25* s36* s46* s5n + 4.* MMU* s23* s34* s56* s5n + 4.* MMU* s26* s34* s56* s5n - 4.* MMU* s24* s36* s56* s5n + 8.* ME*ME* MMU* s23* s34* s6n - 8.* ME*ME* MMU* s26* s34* s6n + 8.* ME*ME* MMU* s24* s36* s6n + 4.* MMU* s23* s35* s45* s6n + 8.* MMU* s26* s35* s45* s6n - 4.* MMU* s23* s36* s45* s6n - 8.* MMU* s25* s36* s45* s6n + 4.* MMU* s23* s34* s56* s6n + 4.* MMU* s25* s34* s56* s6n - 4.* MMU* s24* s35* s56* s6n;

  Ip13r = 16.* ME*ME*ME*ME* s12* s14 - 8.* ME*ME* s12* s14* s15 + 8.* ME*ME*ME*ME* s14* s23 + 8.* ME*ME* MMU*MMU* s14* s23 + 4.* ME*ME* s14* s15* s23 - 4.* ME*ME* s14* s16* s23 - 4.* s14* s15* s16* s23 - 32.* ME*ME*ME*ME* MMU*MMU* s24 - 8.* ME*ME*ME*ME* s13* s24 - 16.* ME*ME*ME*ME* s15* s24 + 8.* ME*ME* s13* s15* s24 - 8.* ME*ME*ME*ME* s16* s24 + 8.* ME*ME* s13* s16* s24 + 8.* ME*ME* s15* s16* s24 + 8.* ME*ME*ME*ME* s14* s25 + 8.* ME*ME* MMU*MMU* s14* s25 - 4.* ME*ME* s13* s14* s25 - 4.* ME*ME* s14* s16* s25 + 8.* ME*ME*ME*ME* s14* s26 + 8.* ME*ME* MMU*MMU* s14* s26 - 4.* ME*ME* s13* s14* s26 + 4.* ME*ME* s14* s15* s26 - 4.* s13* s14* s15* s26 - 16.* ME*ME*ME*ME* MMU* s14* s2n - 8.* ME*ME*ME*ME* s12* s34 - 4.* ME*ME* s12* s15* s34 - 4.* ME*ME* s12* s16* s34 + 8.* ME*ME* MMU*MMU* s25* s34 + 4.* ME*ME* s15* s25* s34 + 4.* ME*ME* s16* s25* s34 + 8.* ME*ME* MMU*MMU* s26* s34 - 4.* s15* s15* s26* s34 + 8.* ME*ME*ME*ME* MMU* s2n* s34 + 4.* ME*ME* MMU* s15* s2n* s34 + 4.* ME*ME* MMU* s16* s2n* s34 + 8.* ME*ME* s12* s14* s35 - 16.* ME*ME* MMU*MMU* s24* s35 - 4.* ME*ME* s15* s24* s35 - 4.* ME*ME* s16* s24* s35 + 4.* ME*ME* s14* s26* s35 + 4.* MMU*MMU* s14* s26* s35 + 4.* s14* s15* s26* s35 - 8.* ME*ME* MMU* s14* s2n* s35 + 8.* ME*ME* s12* s14* s36 + 4.* s12* s14* s15* s36 - 16.* ME*ME* MMU*MMU* s24* s36 + 4* s15* s15* s24* s36 - 4.* ME*ME* s14* s25* s36 - 4.* MMU*MMU* s14* s25* s36 - 4.* s14* s15* s25* s36 - 8.* ME*ME* MMU* s14* s2n* s36 + 8.* ME*ME*ME*ME* MMU* s24* s3n - 4.* ME*ME* MMU* s15* s24* s3n - 4.* ME*ME* MMU* s16* s24* s3n + 4.* ME*ME* MMU* s14* s25* s3n + 4.* ME*ME* MMU* s14* s26* s3n - 16.* ME*ME*ME*ME* s12* s45 - 4.* ME*ME* s12* s13* s45 - 4.* ME*ME* s12* s16* s45 - 8.* ME*ME* MMU*MMU* s23* s45 - 8.* ME*ME* s15* s23* s45 + 4.* s15* s16* s23* s45 + 4.* ME*ME* s13* s25* s45 + 4.* ME*ME* s16* s25* s45 - 8.* ME*ME* MMU*MMU* s26* s45 - 8.* ME*ME* s15* s26* s45 + 4.* s13* s15* s26* s45 + 16.* ME*ME*ME*ME* MMU* s2n* s45 + 4.* ME*ME* MMU* s13* s2n* s45 + 4.* ME*ME* MMU* s16* s2n* s45 - 4.* ME*ME* s12* s35* s45 - 8.* MMU*MMU* s26* s35* s45 - 4.* s15* s26* s35* s45 + 4.* ME*ME* MMU* s2n* s35* s45 - 4.* s12* s15* s36* s45 + 8.* MMU*MMU* s25* s36* s45 + 4.* s15* s25* s36* s45 + 4.* MMU* s15* s2n* s36* s45 - 4.* ME*ME* MMU* s25* s3n* s45 - 4.* MMU* s15* s26* s3n* s45 - 8.* ME*ME*ME*ME* s12* s46 - 4.* ME*ME* s12* s13* s46 - 4.* ME*ME* s12* s15* s46 + 8.* ME*ME* MMU*MMU* s23* s46 - 4.* s15* s15* s23* s46 + 8.* ME*ME* MMU*MMU* s25* s46 + 4.* ME*ME* s13* s25* s46 + 4.* ME*ME* s15* s25* s46 + 8.* ME*ME*ME*ME* MMU* s2n* s46 + 4.* ME*ME* MMU* s13* s2n* s46 + 4.* ME*ME* MMU* s15* s2n* s46 - 4.* ME*ME* s12* s35* s46 + 4.* ME*ME* MMU* s2n* s35* s46 - 4.* ME*ME* MMU* s25* s3n* s46 + 8.* ME*ME* MMU* s12* s15* s4n - 8.* ME*ME*ME*ME* MMU* s23* s4n - 8.* ME*ME* MMU*MMU*MMU* s23* s4n + 4.* MMU* s15* s16* s23* s4n - 8.* ME*ME*ME*ME* MMU* s25* s4n - 8.* ME*ME* MMU*MMU*MMU* s25* s4n - 8.* ME*ME*ME*ME* MMU* s26* s4n - 8.* ME*ME* MMU*MMU*MMU* s26* s4n + 4.* MMU* s13* s15* s26* s4n - 4.* ME*ME* MMU* s26* s35* s4n - 4.* MMU*MMU*MMU* s26* s35* s4n - 4.* MMU* s12* s15* s36* s4n + 4.* ME*ME* MMU* s25* s36* s4n + 4.* MMU*MMU*MMU* s25* s36* s4n + 8.* ME*ME* s12* s14* s56 + 4.* ME*ME* s14* s23* s56 + 4.* MMU*MMU* s14* s23* s56 + 4.* s14* s15* s23* s56 - 16.* ME*ME* MMU*MMU* s24* s56 - 4.* ME*ME* s13* s24* s56 - 4.* ME*ME* s15* s24* s56 - 8.* ME*ME* MMU* s14* s2n* s56 - 4.* ME*ME* s12* s34* s56 + 4.* ME*ME* MMU* s2n* s34* s56 + 4.* ME*ME* MMU* s24* s3n* s56 - 4.* ME*ME* s12* s45* s56 - 8.* MMU*MMU* s23* s45* s56 - 4.* s15* s23* s45* s56 + 4.* ME*ME* MMU* s2n* s45* s56 - 4.* ME*ME* MMU* s23* s4n* s56 - 4.* MMU*MMU*MMU* s23* s4n* s56 - 4.* ME*ME* MMU* s14* s23* s5n + 16.* ME*ME*ME*ME* MMU* s24* s5n - 4.* ME*ME* MMU* s13* s24* s5n - 4.* ME*ME* MMU* s16* s24* s5n - 4.* ME*ME* MMU* s14* s26* s5n - 4.* ME*ME* MMU* s25* s34* s5n + 4.* MMU* s15* s26* s34* s5n + 4.* ME*ME* MMU* s24* s35* s5n - 4.* MMU* s14* s26* s35* s5n - 4.* MMU* s15* s24* s36* s5n + 4.* MMU* s14* s25* s36* s5n + 8.* ME*ME* MMU* s23* s45* s5n + 8.* ME*ME* MMU* s26* s45* s5n + 4.* MMU* s26* s35* s45* s5n - 4.* MMU* s25* s36* s45* s5n + 4.* MMU* s15* s23* s46* s5n - 4.* ME*ME* MMU* s25* s46* s5n - 4.* MMU* s14* s23* s56* s5n + 4.* ME*ME* MMU* s24* s56* s5n + 4.* MMU* s23* s45* s56* s5n + 4.* ME*ME* MMU* s14* s23* s6n + 8.* ME*ME*ME*ME* MMU* s24* s6n - 4.* ME*ME* MMU* s13* s24* s6n - 4.* ME*ME* MMU* s15* s24* s6n + 4.* ME*ME* MMU* s14* s25* s6n - 4.* ME*ME* MMU* s25* s34* s6n + 4.* ME*ME* MMU* s24* s35* s6n - 4.* MMU* s15* s23* s45* s6n - 4.* ME*ME* MMU* s25* s45* s6n;

  Ip14r = -8.* ME*ME* s13* s14* s23 - 16.* ME*ME* s14* s15* s23 + 16.* ME*ME* s14* s16* s23 + 16.* ME*ME*ME*ME* s13* s24 + 8.* ME*ME*ME*ME* s15* s24 - 8.* ME*ME*ME*ME* s16* s24 - 8.* ME*ME* s14* s15* s25 + 8.* ME*ME* s14* s16* s25 - 8.* ME*ME* s13* s14* s26 - 8.* ME*ME* s14* s15* s26 + 8.* ME*ME* s14* s16* s26 - 16.* ME*ME*ME*ME* s12* s34 - 8.* ME*ME* s15* s23* s34 + 4.* ME*ME* s15* s25* s34 + 4.* ME*ME* s16* s25* s34 - 4.* ME*ME* s15* s26* s34 + 4.* ME*ME* s16* s26* s34 + 16.* ME*ME*ME*ME* MMU* s2n* s34 + 4.* s14* s16* s23* s35 - 4.* ME*ME* s15* s24* s35 - 4.* ME*ME* s16* s24* s35 + 4.* s14* s16* s26* s35 + 8.* ME*ME* s12* s14* s36 + 4.* s14* s15* s23* s36 + 4.* ME*ME* s15* s24* s36 - 4.* ME*ME* s16* s24* s36 - 4.* s14* s16* s25* s36 - 16.* ME*ME*ME*ME* MMU* s24* s3n - 8.* ME*ME*ME*ME* s12* s45 + 8.* ME*ME* s13* s23* s45 - 16.* ME*ME* s16* s23* s45 - 4.* ME*ME* s13* s25* s45 - 4.* ME*ME* s16* s25* s45 + 4.* ME*ME* s13* s26* s45 - 4.* ME*ME* s16* s26* s45 + 8.* ME*ME*ME*ME* MMU* s2n* s45 + 4.* ME*ME* s12* s35* s45 - 4.* s16* s23* s35* s45 - 4.* s16* s26* s35* s45 - 4.* ME*ME* MMU* s2n* s35* s45 - 4.* ME*ME* s12* s36* s45 + 4.* s16* s25* s36* s45 + 4.* ME*ME* MMU* s2n* s36* s45 - 8.* ME*ME* MMU* s23* s3n* s45 + 4.* ME*ME* MMU* s25* s3n* s45 - 4.* ME*ME* MMU* s26* s3n* s45 + 8.* ME*ME*ME*ME* s12* s46 + 16.* ME*ME* s15* s23* s46 - 4.* ME*ME* s13* s25* s46 + 4.* ME*ME* s15* s25* s46 - 4.* ME*ME* s13* s26* s46 + 4.* ME*ME* s15* s26* s46 - 8.* ME*ME*ME*ME* MMU* s2n* s46 + 4.* ME*ME* s12* s35* s46 + 4.* s15* s23* s35* s46 + 4.* s15* s26* s35* s46 - 4.* ME*ME* MMU* s2n* s35* s46 + 4.* ME*ME* s12* s36* s46 - 4.* s15* s25* s36* s46 - 4.* ME*ME* MMU* s2n* s36* s46 + 4.* ME*ME* MMU* s25* s3n* s46 + 4.* ME*ME* MMU* s26* s3n* s46 + 8.* ME*ME* MMU* s13* s23* s4n + 16.* ME*ME* MMU* s15* s23* s4n - 16.* ME*ME* MMU* s16* s23* s4n + 8.* ME*ME* MMU* s15* s25* s4n - 8.* ME*ME* MMU* s16* s25* s4n + 8.* ME*ME* MMU* s13* s26* s4n + 8.* ME*ME* MMU* s15* s26* s4n - 8.* ME*ME* MMU* s16* s26* s4n - 4.* MMU* s16* s23* s35* s4n - 4.* MMU* s16* s26* s35* s4n - 8.* ME*ME* MMU* s12* s36* s4n - 4.* MMU* s15* s23* s36* s4n + 4.* MMU* s16* s25* s36* s4n - 4.* s13* s14* s23* s56 + 8.* ME*ME* s13* s24* s56 + 4.* ME*ME* s15* s24* s56 - 4.* ME*ME* s16* s24* s56 - 4* s13* s14* s26* s56 - 8.* ME*ME* s12* s34* s56 - 4.* s15* s23* s34* s56 - 4.* s15* s26* s34* s56 + 8.* ME*ME* MMU* s2n* s34* s56 + 4.* s12* s14* s36* s56 + 4.* s15* s24* s36* s56 - 8.* ME*ME* MMU* s24* s3n* s56 - 4.* ME*ME* s12* s45* s56 + 4.* s13* s23* s45* s56 + 4.* s13* s26* s45* s56 + 4.* ME*ME* MMU* s2n* s45* s56 - 4.* s12* s36* s45* s56 + 4.* MMU* s2n* s36* s45* s56 - 4.* MMU* s23* s3n* s45* s56 - 4.* MMU* s26* s3n* s45* s56 + 4.* ME*ME* s12* s46* s56 - 4.* ME*ME* MMU* s2n* s46* s56 + 4.* MMU* s13* s23* s4n* s56 + 4.* MMU* s13* s26* s4n* s56 - 4.* MMU* s12* s36* s4n* s56 - 8.* ME*ME*ME*ME* MMU* s24* s5n + 8.* ME*ME* MMU* s23* s34* s5n - 4.* ME*ME* MMU* s25* s34* s5n + 4.* ME*ME* MMU* s26* s34* s5n + 4.* ME*ME* MMU* s24* s35* s5n - 4.* ME*ME* MMU* s24* s36* s5n - 16.* ME*ME* MMU* s23* s46* s5n - 4.* ME*ME* MMU* s25* s46* s5n - 4.* ME*ME* MMU* s26* s46* s5n - 4.* MMU* s23* s35* s46* s5n - 4.* MMU* s26* s35* s46* s5n + 4.* MMU* s25* s36* s46* s5n - 4.* ME*ME* MMU* s24* s56* s5n + 4.* MMU* s23* s34* s56* s5n + 4.* MMU* s26* s34* s56* s5n - 4.* MMU* s24* s36* s56* s5n + 8.* ME*ME*ME*ME* MMU* s24* s6n - 4.* ME*ME* MMU* s25* s34* s6n - 4.* ME*ME* MMU* s26* s34* s6n + 4.* ME*ME* MMU* s24* s35* s6n + 4.* ME*ME* MMU* s24* s36* s6n + 16.* ME*ME* MMU* s23* s45* s6n + 4.* ME*ME* MMU* s25* s45* s6n + 4.* ME*ME* MMU* s26* s45* s6n + 4.* MMU* s23* s35* s45* s6n + 4.* MMU* s26* s35* s45* s6n - 4.* MMU* s25* s36* s45* s6n + 4.* ME*ME* MMU* s24* s56* s6n;

  Ip23r = 8.* ME*ME* s13* s14* s23 - 8.* ME*ME* s14* s15* s23 - 8.* ME*ME* s14* s16* s23 - 8.* ME*ME*ME*ME* s13* s24 + 8.* ME*ME*ME*ME* s15* s24 + 16.* ME*ME*ME*ME* s16* s24 + 8.* ME*ME* s13* s14* s25 - 8.* ME*ME* s14* s15* s25 + 16.* ME*ME* s13* s14* s26 - 16.* ME*ME* s14* s15* s26 - 8.* ME*ME* s14* s16* s26 + 8.* ME*ME*ME*ME* s12* s34 + 4.* ME*ME* s15* s23* s34 - 4.* ME*ME* s16* s23* s34 + 4.* ME*ME* s15* s25* s34 - 4.* ME*ME* s16* s25* s34 + 16.* ME*ME* s15* s26* s34 - 8.* ME*ME*ME*ME* MMU* s2n* s34 - 4.* s14* s16* s23* s35 - 4.* ME*ME* s13* s24* s35 + 4.* ME*ME* s15* s24* s35 + 8.* ME*ME* s16* s24* s35 - 4.* s14* s16* s26* s35 + 4.* ME*ME* s12* s34* s35 - 4.* ME*ME* MMU* s2n* s34* s35 + 8.* ME*ME* s12* s14* s36 - 4.* ME*ME* s13* s24* s36 + 4.* ME*ME* s15* s24* s36 - 4.* s13* s14* s25* s36 + 4.* s14* s15* s26* s36 + 4.* ME*ME* s12* s34* s36 - 4.* s15* s25* s34* s36 - 4.* ME*ME* MMU* s2n* s34* s36 + 4.* s12* s14* s35* s36 + 4.* s15* s24* s35* s36 + 8.* ME*ME*ME*ME* MMU* s24* s3n + 4.* ME*ME* MMU* s24* s35* s3n + 4.* ME*ME* MMU* s24* s36* s3n - 8.* ME*ME*ME*ME* s12* s45 - 4.* ME*ME* s13* s23* s45 + 4.* ME*ME* s16* s23* s45 - 4.* ME*ME* s13* s25* s45 - 4.* ME*ME* s16* s25* s45 - 16.* ME*ME* s13* s26* s45 + 8.* ME*ME* s16* s26* s45 + 8.* ME*ME*ME*ME* MMU* s2n* s45 - 4.* ME*ME* s12* s35* s45 + 4.* s16* s23* s35* s45 + 4.* s16* s26* s35* s45 + 4.* ME*ME* MMU* s2n* s35* s45 - 4.* ME*ME* s12* s36* s45 + 4.* s13* s25* s36* s45 + 4.* ME*ME* MMU* s2n* s36* s45 - 4.* s12* s35* s36* s45 + 4.* MMU* s2n* s35* s36* s45 + 4.* ME*ME* MMU* s23* s3n* s45 + 4.* ME*ME* MMU* s25* s3n* s45 + 16.* ME*ME* MMU* s26* s3n* s45 - 4.* MMU* s25* s36* s3n* s45 - 16.* ME*ME*ME*ME* s12* s46 + 4.* ME*ME* s13* s23* s46 - 4.* ME*ME* s15* s23* s46 + 4.* ME*ME* s13* s25* s46 + 4.* ME*ME* s15* s25* s46 - 8.* ME*ME* s15* s26* s46 + 16.* ME*ME*ME*ME* MMU* s2n* s46 - 8.* ME*ME* s12* s35* s46 - 4.* s15* s23* s35* s46 - 4.* s15* s26* s35* s46 + 8.* ME*ME* MMU* s2n* s35* s46 - 4.* ME*ME* MMU* s23* s3n* s46 - 4.* ME*ME* MMU* s25* s3n* s46 - 8.* ME*ME* MMU* s13* s23* s4n + 8.* ME*ME* MMU* s15* s23* s4n + 8.* ME*ME* MMU* s16* s23* s4n - 8.* ME*ME* MMU* s13* s25* s4n + 8.* ME*ME* MMU* s15* s25* s4n - 16.* ME*ME* MMU* s13* s26* s4n + 16.* ME*ME* MMU* s15* s26* s4n + 8.* ME*ME* MMU* s16* s26* s4n + 4.* MMU* s16* s23* s35* s4n + 4.* MMU* s16* s26* s35* s4n - 8.* ME*ME* MMU* s12* s36* s4n + 4.* MMU* s13* s25* s36* s4n - 4.* MMU* s15* s26* s36* s4n - 4.* MMU* s12* s35* s36* s4n + 4.* s13* s14* s23* s56 - 4.* ME*ME* s13* s24* s56 - 4.* ME*ME* s15* s24* s56 + 4.* s13* s14* s26* s56 + 4.* ME*ME* s12* s34* s56 + 4.* s15* s23* s34* s56 + 4.* s15* s26* s34* s56 - 4.* ME*ME* MMU* s2n* s34* s56 + 4.* ME*ME* MMU* s24* s3n* s56 + 4.* ME*ME* s12* s45* s56 - 4.* s13* s23* s45* s56 - 4.* s13* s26* s45* s56 - 4.* ME*ME* MMU* s2n* s45* s56 + 4.* MMU* s23* s3n* s45* s56 + 4.* MMU* s26* s3n* s45* s56 - 4.* MMU* s13* s23* s4n* s56 - 4.* MMU* s13* s26* s4n* s56 - 8.* ME*ME*ME*ME* MMU* s24* s5n - 4.* ME*ME* MMU* s23* s34* s5n - 4.* ME*ME* MMU* s25* s34* s5n - 16.* ME*ME* MMU* s26* s34* s5n - 4.* ME*ME* MMU* s24* s35* s5n - 4.* ME*ME* MMU* s24* s36* s5n + 4.* MMU* s25* s34* s36* s5n - 4.* MMU* s24* s35* s36* s5n + 4.* ME*ME* MMU* s23* s46* s5n - 4.* ME*ME* MMU* s25* s46* s5n + 8.* ME*ME* MMU* s26* s46* s5n + 4.* MMU* s23* s35* s46* s5n + 4.* MMU* s26* s35* s46* s5n + 4.* ME*ME* MMU* s24* s56* s5n - 4.* MMU* s23* s34* s56* s5n - 4.* MMU* s26* s34* s56* s5n - 16.* ME*ME*ME*ME* MMU* s24* s6n + 4.* ME*ME* MMU* s23* s34* s6n + 4.* ME*ME* MMU* s25* s34* s6n - 8.* ME*ME* MMU* s24* s35* s6n - 4.* ME*ME* MMU* s23* s45* s6n + 4.* ME*ME* MMU* s25* s45* s6n - 8.* ME*ME* MMU* s26* s45* s6n - 4.* MMU* s23* s35* s45* s6n - 4.* MMU* s26* s35* s45* s6n;

  Ip24r = -32.* ME*ME*ME*ME* s14* s23 - 32.* ME*ME*ME*ME* s14* s26 - 16.* ME*ME* s14* s23* s35 - 16.* ME*ME* s14* s26* s35 + 32.* ME*ME* s14* s23* s36 + 32.* ME*ME* s14* s25* s36 + 32.* ME*ME* s14* s26* s36 + 8.* s14* s23* s35* s36 + 8.* s14* s26* s35* s36 - 8.* s14* s25* s36*s36 + 32.* ME*ME*ME*ME* MMU* s23* s4n + 32.* ME*ME*ME*ME* MMU* s26* s4n + 16.* ME*ME* MMU* s23* s35* s4n + 16.* ME*ME* MMU* s26* s35* s4n - 32.* ME*ME* MMU* s23* s36* s4n - 32.* ME*ME* MMU* s25* s36* s4n - 32.* ME*ME* MMU* s26* s36* s4n - 8.* MMU* s23* s35* s36* s4n - 8.* MMU* s26* s35* s36* s4n + 8.* MMU* s25* s36*s36* s4n - 16.* ME*ME* s14* s23* s56 - 16.* ME*ME* s14* s26* s56 + 8.* s14* s23* s36* s56 + 8.* s14* s26* s36* s56 + 16.* ME*ME* MMU* s23* s4n* s56 + 16.* ME*ME* MMU* s26* s4n* s56 - 8.* MMU* s23* s36* s4n* s56 - 8.* MMU* s26* s36* s4n* s56;

  Ip34r = -8.* ME*ME* s14* s16* s23 + 32.* ME*ME*ME*ME* s16* s24 - 8.* ME*ME* s14* s16* s25 - 16.* ME*ME* s14* s16* s26 - 8.* ME*ME* s16* s23* s34 + 8.* ME*ME* s16* s26* s34 - 4.* s14* s16* s23* s35 + 16.* ME*ME* s16* s24* s35 - 4.* s14* s16* s25* s35 - 8.* s14* s16* s26* s35 + 4.* s16* s25* s34* s35 + 4.* s16* s26* s34* s35 + 8.* ME*ME* s12* s14* s36 - 8.* ME*ME* s13* s24* s36 - 4.* s13* s14* s25* s36 + 4.* s14* s15* s25* s36 + 8.* s14* s15* s26* s36 + 8.* ME*ME* s12* s34* s36 - 8.* s15* s25* s34* s36 - 4.* s15* s26* s34* s36 - 8.* ME*ME* MMU* s2n* s34* s36 + 4.* s12* s14* s35* s36 + 4.* s15* s24* s35* s36 + 8.* ME*ME* MMU* s24* s36* s3n - 8.* ME*ME* s16* s25* s45 + 8.* ME*ME* s16* s26* s45 + 4.* s16* s23* s35* s45 + 4.* s16* s26* s35* s45 + 8.* s13* s25* s36* s45 + 4.* s13* s26* s36* s45 - 4.* s12* s35* s36* s45 + 4.* MMU* s2n* s35* s36* s45 - 8.* MMU* s25* s36* s3n* s45 - 4.* MMU* s26* s36* s3n* s45 - 32.* ME*ME*ME*ME* s12* s46 + 8.* ME*ME* s13* s23* s46 + 8.* ME*ME* s15* s25* s46 - 8.* ME*ME* s13* s26* s46 - 8.* ME*ME* s15* s26* s46 + 32.* ME*ME*ME*ME* MMU* s2n* s46 - 16.* ME*ME* s12* s35* s46 - 4.* s15* s23* s35* s46 - 4.* s13* s25* s35* s46 - 4.* s13* s26* s35* s46 - 4.* s15* s26* s35* s46 + 16.* ME*ME* MMU* s2n* s35* s46 - 8.* ME*ME* MMU* s23* s3n* s46 + 8.* ME*ME* MMU* s26* s3n* s46 + 4.* MMU* s25* s35* s3n* s46 + 4.* MMU* s26* s35* s3n* s46 + 8.* ME*ME* MMU* s16* s23* s4n + 8.* ME*ME* MMU* s16* s25* s4n + 16.* ME*ME* MMU* s16* s26* s4n + 4.* MMU* s16* s23* s35* s4n + 4.* MMU* s16* s25* s35* s4n + 8.* MMU* s16* s26* s35* s4n - 8* ME*ME* MMU* s12* s36* s4n + 4.* MMU* s13* s25* s36* s4n - 4.* MMU* s15* s25* s36* s4n - 8.* MMU* s15* s26* s36* s4n - 4.* MMU* s12* s35* s36* s4n + 8.* ME*ME* s12* s14* s56 + 4.* s13* s14* s23* s56 - 4.* s14* s15* s23* s56 - 8.* ME*ME* s15* s24* s56 + 8.* s13* s14* s26* s56 + 8.* s15* s23* s34* s56 + 4.* s15* s26* s34* s56 + 4.* s12* s14* s35* s56 + 4.* s13* s24* s35* s56 - 4.* s12* s34* s35* s56 + 4.* MMU* s2n* s34* s35* s56 - 4.* MMU* s24* s35* s3n* s56 + 8.* ME*ME* s12* s45* s56 - 8.* s13* s23* s45* s56 - 4.* s13* s26* s45* s56 - 8.* ME*ME* MMU* s2n* s45* s56 + 8.* MMU* s23* s3n* s45* s56 + 4.* MMU* s26* s3n* s45* s56 - 8.* ME*ME* MMU* s12* s4n* s56 - 4.* MMU* s13* s23* s4n* s56 + 4.* MMU* s15* s23* s4n* s56 - 8.* MMU* s13* s26* s4n* s56 - 4.* MMU* s12* s35* s4n* s56 + 8.* MMU* s25* s34* s36* s5n + 4.* MMU* s26* s34* s36* s5n - 4.* MMU* s24* s35* s36* s5n - 8.* ME*ME* MMU* s25* s46* s5n + 8.* ME*ME* MMU* s26* s46* s5n + 4.* MMU* s23* s35* s46* s5n + 4.* MMU* s26* s35* s46* s5n + 8.* ME*ME* MMU* s24* s56* s5n - 8.* MMU* s23* s34* s56* s5n - 4.* MMU* s26* s34* s56* s5n - 32.* ME*ME*ME*ME* MMU* s24* s6n + 8.* ME*ME* MMU* s23* s34* s6n - 8.* ME*ME* MMU* s26* s34* s6n - 16.* ME*ME* MMU* s24* s35* s6n - 4.* MMU* s25* s34* s35* s6n - 4.* MMU* s26* s34* s35* s6n + 8.* ME*ME* MMU* s25* s45* s6n - 8.* ME*ME* MMU* s26* s45* s6n - 4.* MMU* s23* s35* s45* s6n - 4.* MMU* s26* s35* s45* s6n;

  den1 = (2.* ME*ME - s15 - s16 + s56) * (2.* ME*ME + s56);
  den2 = (2.* ME*ME + s35 + s36 + s56) * (2.* ME*ME + s56);
  den3 = (2.* ME*ME - s13 - s15 + s35) * (2.* ME*ME + s35);
  den4 = (2.* ME*ME + s35 + s36 + s56) * (2.* ME*ME + s35);

  MatElSquared = 8.* GF*GF* (4.* M_PI* ALPHA)*(4.* M_PI* ALPHA)* (Ip11/(den1*den1) + Ip22/(den2*den2) + Ip33/(den3*den3) + Ip44/(den4*den4) + (2.* Ip12r)/(den1* den2) + (2.* Ip13r)/(den1* den3) + (2.* Ip14r)/(den1* den4) + (2.* Ip23r)/(den2* den3) + (2.* Ip24r)/(den2* den4) + (2.* Ip34r)/(den3* den4));

  // MatElSquared = 8.* (4.* M_PI* ALPHA)*(4.* M_PI* ALPHA)* (Ip11/(den1*den1) + Ip22/(den2*den2) + Ip33/(den3*den3) + Ip44/(den4*den4) + (2.* Ip12r)/(den1* den2) + (2.* Ip13r)/(den1* den3) + (2.* Ip14r)/(den1* den4) + (2.* Ip23r)/(den2* den3) + (2.* Ip24r)/(den2* den4) + (2.* Ip34r)/(den3* den4));        // no factor GF^2
  // MatElSquared = (Ip11/(den1*den1) + Ip22/(den2*den2) + Ip33/(den3*den3) + Ip44/(den4*den4) + (2.* Ip12r)/(den1* den2) + (2.* Ip13r)/(den1* den3) + (2.* Ip14r)/(den1* den4) + (2.* Ip23r)/(den2* den3) + (2.* Ip24r)/(den2* den4) + (2.* Ip34r)/(den3* den4));  // no pre-factors


  return MatElSquared;
}

double Mu3eMuonInternalConversionDecayWithSpin::matmu3e2nu_DK(double P[5][4]) {              // matrix element by Djilkibaev & Konoplich (PhysRev D 79)

  /* order of FS particles:
     e+,e-,e+,numu,nue
  */

  double C,C1,C2,C3,D1,D2;
  double tr11,tr12,tr13,tr14,tr22,tr23,tr24,tr33,tr34,tr44;
  double matr2,matr2e,matr2mu,matr2emu;

  double qp,qp1,qp2,qk1,qk2;
  double pp1,pp2,pk1,pk2,p1p2,p1k1,p1k2,p2k1,p2k2,k1k2;
  double qps,qp12,qp22,pp12,pp22,p1p22;
  double m2,m4,u2;
  //int k;
  /*
    for (k=0;k<5;k++) {
    printf("%G %G %G\n",P[k][0],P[k][1],P[k][2]);
    }
  */
  double ALPHA = CLHEP::fine_structure_const;
  C = 384*ALPHA*ALPHA/pow(M_PI*MMU,6);
  m2=ME*ME;
  m4=m2*m2;
  u2=MMU*MMU;

  qp  = MMU*P[0][3];
  qp1 = MMU*P[1][3];
  qp2 = MMU*P[2][3];
  qk1 = MMU*P[3][3];
  qk2 = MMU*P[4][3];

  pp1= P[0][3]*P[1][3]-P[0][0]*P[1][0]-P[0][1]*P[1][1]-P[0][2]*P[1][2];
  pp2= P[0][3]*P[2][3]-P[0][0]*P[2][0]-P[0][1]*P[2][1]-P[0][2]*P[2][2];
  pk1= P[0][3]*P[3][3]-P[0][0]*P[3][0]-P[0][1]*P[3][1]-P[0][2]*P[3][2];
  pk2= P[0][3]*P[4][3]-P[0][0]*P[4][0]-P[0][1]*P[4][1]-P[0][2]*P[4][2];

  p1p2=P[1][3]*P[2][3]-P[1][0]*P[2][0]-P[1][1]*P[2][1]-P[1][2]*P[2][2];
  p1k1=P[1][3]*P[3][3]-P[1][0]*P[3][0]-P[1][1]*P[3][1]-P[1][2]*P[3][2];
  p1k2=P[1][3]*P[4][3]-P[1][0]*P[4][0]-P[1][1]*P[4][1]-P[1][2]*P[4][2];

  p2k1=P[2][3]*P[3][3]-P[2][0]*P[3][0]-P[2][1]*P[3][1]-P[2][2]*P[3][2];
  p2k2=P[2][3]*P[4][3]-P[2][0]*P[4][0]-P[2][1]*P[4][1]-P[2][2]*P[4][2];

  k1k2=P[3][3]*P[4][3]-P[3][0]*P[4][0]-P[3][1]*P[4][1]-P[3][2]*P[4][2];


  qps = qp*qp;
  qp12 = qp1*qp1;
  qp22 = qp2*qp2;
  pp12 = pp1*pp1;
  pp22 = pp2*pp2;
  p1p22 = p1p2*p1p2;
  C1 = 1.0/(2.0*(m2 + pp1 + pp2 + p1p2));
  C2 = 1.0/(2.0*(m2 - qp1 - qp2 + p1p2));
  C3 = 1.0/(2.0*(m2 - qp - qp1 + pp1));
  D1 = 1.0/(2.0*(m2 + p1p2));
  D2 = 1.0/(2.0*(m2 + pp1));


  tr11= -(qk2*(p2k1*(pp12 - pp1*(m2 + pp2) + m2*(m2 + p1p2) - pp2*(2.*m2 + p1p2))
	       + p1k1*(m4 - m2*pp2 + pp22 + m2*p1p2 - pp1*(2.*m2 + pp2 + p1p2))
	       + pk1*((2.*m2 - pp2)*(m2 + p1p2) - pp1*(m2 + 2.*pp2 + p1p2))));

  tr12 = m2*pk1*p1k2*qp - m2*p1k1*p1k2*qp + m2*pk1*p2k2*qp - m2*p2k1*p2k2*qp - 2.*m2*pk1*qk2*qp
    - m2*p1k1*qk2*qp - m2*p2k1*qk2*qp + pk1*p1k2*qp*p1p2 + p2k1*p1k2*qp*p1p2
    + pk1*p2k2*qp*p1p2 + p1k1*p2k2*qp*p1p2 - 2.*pk1*qk2*qp*p1p2 - p1k1*qk2*qp*p1p2
    - p2k1*qk2*qp*p1p2
    + qk1*(m2*qk2*pp1 + m2*p2k2*pp2 + m2*qk2*pp2 - p2k2*pp1*p1p2
	   + qk2*pp1*p1p2 + qk2*pp2*p1p2 - 2.*m2*pk2*(m2+p1p2) + p1k2*(m2*pp1-pp2*p1p2))
    - m2*pk1*pk2*qp1 + m2*p1k1*pk2*qp1 + pk1*p2k2*pp1*qp1 + 2.*p2k1*p2k2*pp1*qp1
    - p2k1*qk2*pp1*qp1 - pk1*p2k2*pp2*qp1 - 2.*p1k1*p2k2*pp2*qp1 + 2.*pk1*qk2*pp2*qp1
    + p1k1*qk2*pp2*qp1 - pk1*pk2*p1p2*qp1 - p2k1*pk2*p1p2*qp1 - m2*pk1*pk2*qp2
    + m2*p2k1*pk2*qp2 - pk1*p1k2*pp1*qp2 - 2.*p2k1*p1k2*pp1*qp2
    + 2.*pk1*qk2*pp1*qp2 + p2k1*qk2*pp1*qp2 + pk1*p1k2*pp2*qp2
    + 2.*p1k1*p1k2*pp2*qp2 - p1k1*qk2*pp2*qp2 - pk1*pk2*p1p2*qp2 - p1k1*pk2*p1p2*qp2
    + k1k2*(2.*m2*qp*(m2 + p1p2) + pp2*(p1p2*qp1 - m2*qp2) + pp1*(-(m2*qp1) + p1p2*qp2));

  tr13 = 2.*qk2*(p1k1*pp2*(-2.*m2+pp2)
		 + pk1*(pp1*(m2-pp2) + m2*(m2+p1p2)
			- pp2*(2.*m2+p1p2))
		 + p2k1*(pp1*(m2-pp2) + m2*(m2 + p1p2) - pp2*(2.*m2 + p1p2)));

  tr14 = (m2*pk1*p1k2*qp + m2*p1k1*p1k2*qp + 4.*m2*p2k1*p1k2*qp - m2*pk1*p2k2*qp
	  - m2*p1k1*p2k2*qp
	  - 2.*m2*pk1*qk2*qp - 2.*m2*p1k1*qk2*qp - 4.*m2*p2k1*qk2*qp - 2.*p1k1*p1k2*pp2*qp
	  + 2.*p1k1*qk2*pp2*qp + 2.*pk1*p1k2*qp*p1p2 + 2.*p2k1*p1k2*qp*p1p2
	  - 2.*pk1*qk2*qp*p1p2
	  - 2.*p2k1*qk2*qp*p1p2
	  - qk1*(-2.*(m2 + pp1)*(m2*p2k2 - qk2*pp2) - p1k2*(pp1*(m2 + 2.*pp2)
							    + m2*(m2 + pp2 - p1p2))
		 + m2*pk2*(m2 + pp1 + pp2 + p1p2)) - m2*pk1*pk2*qp1
	  - m2*p1k1*pk2*qp1 - 4.*m2*p2k1*pk2*qp1 + m2*pk1*p2k2*qp1 - m2*p1k1*p2k2*qp1
	  + 2.*m2*p2k1*p2k2*qp1 + 2.*m2*pk1*qk2*qp1 + 2.*m2*p1k1*qk2*qp1
	  + 4.*m2*p2k1*qk2*qp1
	  + 2.*pk1*p2k2*pp1*qp1 + 2.*p2k1*p2k2*pp1*qp1 + 2.*p1k1*pk2*pp2*qp1
	  - 2.*p2k1*qk2*pp2*qp1
	  - 2.*pk1*pk2*p1p2*qp1 - 2.*p2k1*pk2*p1p2*qp1 + m2*pk1*pk2*qp2 + m2*p1k1*pk2*qp2
	  - m2*pk1*p1k2*qp2 + m2*p1k1*p1k2*qp2 - 2.*m2*p2k1*p1k2*qp2 + 2.*m2*pk1*qk2*qp2
	  + 2.*m2*p2k1*qk2*qp2 - 2.*pk1*p1k2*pp1*qp2 - 2.*p2k1*p1k2*pp1*qp2
	  + 2.*pk1*qk2*pp1*qp2
	  + 2.*p2k1*qk2*pp1*qp2
	  + k1k2*(m2*qp*(m2 + pp1 + pp2 + p1p2)
		  - (pp1*(m2 + 2.*pp2)
		     + m2*(m2 + pp2 - p1p2))*qp1 - 2.*m2*(m2 + pp1)*qp2))/2.0;

  tr22= -(pk1*(-(p1k2*(m2*u2 + p1p2*(u2 + qp1) + qp1*(2.*m2 - qp2) + m2*qp2 + qp22))
	       + qk2*(qp1*(m2 - 2.*qp2) + m2*(m2 + u2 + qp2) + p1p2*(m2 + u2 + qp1 + qp2))
	       - p2k2*(qp12 + qp1*(m2 - qp2) + p1p2*(u2 + qp2) + m2*(u2 + 2.*qp2))));

  tr23 = (-2.*m2*pk1*p1k2*qp + m2*p1k1*p1k2*qp - m2*p2k1*p1k2*qp + m2*p1k1*p2k2*qp
	  + m2*p2k1*p2k2*qp + 2.*m2*pk1*qk2*qp + 2.*m2*p2k1*qk2*qp - 2.*pk1*p1k2*qp*p1p2
	  - 2.*p2k1*p1k2*qp*p1p2 + 2.*pk1*qk2*qp*p1p2 + 2.*p2k1*qk2*qp*p1p2
	  - qk1*(-2.*(m2*pk2
		      - qk2*pp2)*(m2 + p1p2)
		 + m2*p2k2*(m2 + pp1 + pp2 + p1p2)
		 - p1k2*(m2*(m2 - pp1 + pp2)
			 + (m2 + 2.*pp2)*p1p2))
	  + 2.*m2*pk1*pk2*qp1 - m2*p1k1*pk2*qp1 + m2*p2k1*pk2*qp1
	  - 4.*m2*pk1*p2k2*qp1 - m2*p1k1*p2k2*qp1 - m2*p2k1*p2k2*qp1 + 4.*m2*pk1*qk2*qp1
	  + 2.*m2*p1k1*qk2*qp1 + 2.*m2*p2k1*qk2*qp1 - 2.*pk1*p2k2*pp1*qp1
	  - 2.*p2k1*p2k2*pp1*qp1
	  + 2.*p1k1*p2k2*pp2*qp1 - 2.*pk1*qk2*pp2*qp1 + 2.*pk1*pk2*p1p2*qp1
	  + 2.*p2k1*pk2*p1p2*qp1
	  - m2*p1k1*pk2*qp2 - m2*p2k1*pk2*qp2 + 4.*m2*pk1*p1k2*qp2 + m2*p1k1*p1k2*qp2
	  + m2*p2k1*p1k2*qp2 - 4.*m2*pk1*qk2*qp2
	  - 2.*m2*p1k1*qk2*qp2 - 2.*m2*p2k1*qk2*qp2 + 2.*pk1*p1k2*pp1*qp2
	  + 2.*p2k1*p1k2*pp1*qp2
	  - 2.*pk1*qk2*pp1*qp2 - 2.*p2k1*qk2*pp1*qp2 - 2.*p1k1*p1k2*pp2*qp2
	  + 2.*p1k1*qk2*pp2*qp2
	  + k1k2*(-2.*m2*qp*(m2 + p1p2) - (m2*(m2 - pp1 + pp2) + (m2 + 2.*pp2)*p1p2)*qp1
		  + m2*(m2 + pp1 + pp2 + p1p2)*qp2))/2.0;

  tr24 = (qp1*(-(m2*p2k1*pk2) - u2*p2k1*pk2 + m2*qk1*pk2 + m2*pk1*p1k2 + m2*p2k1*p1k2
	       - m2*pk1*p2k2- u2*pk1*p2k2 + m2*qk1*p2k2 - m2*pk1*qk2 - m2*p2k1*qk2
	       + 2.*p2k1*p1k2*pp1 - 2.*p2k1*qk2*pp1 + 2.*qk1*p1k2*pp2 - 2.*qk1*qk2*pp2
	       - p1k1*(m2*pk2 + m2*p2k2 + 2.*(p1k2 - qk2)*pp2)
	       - 2.*p2k1*p1k2*qp + 2.*p2k1*qk2*qp + 2.*pk1*p1k2*p1p2 - 2.*pk1*qk2*p1p2
	       + 2.*p2k1*pk2*qp1 + 2.*pk1*p2k2*qp1
	       + k1k2*(m2*pp1 + pp2*(m2 + u2 - 2.*qp1) + m2*(m2 - qp + p1p2 - qp2))
	       - 2.*pk1*p1k2*qp2 + 2.*pk1*qk2*qp2))/2.
    + u2*((m2*pk1*p1k2 - 2.*m2*pk1*p2k2 + m2*k1k2*pp1 + 2.*m2*k1k2*pp2
	   - p1k1*(m2*pk2 + m2*p2k2 + 2.*(2.*p1k2 - qk2)*pp2) + m2*k1k2*p1p2
	   + 4.*pk1*p1k2*p1p2 - 2.*pk1*qk2*p1p2
	   + p2k1*(-2.*qk2*pp1 + p1k2*(m2 + 4.*pp1) - 2.*pk2*(m2 - qp1))
	   + 2.*pk1*p2k2*qp1 - 2.*k1k2*pp2*qp1)/4.)
    + m2*((2.*m2*qk1*pk2 -u2*qk1*pk2-2.*u2*pk1*p1k2+4.*m2*qk1*p1k2
	   -2.*u2*qk1*p1k2-2.*u2*pk1*p2k2+2.*m2*qk1*p2k2
	   -u2*qk1*p2k2-2.*m2*pk1*qk2+u2*pk1*qk2-2.*m2*p1k1*qk2-4.*m2*qk1*qk2
	   +2.*qk1*p1k2*pp1+2.*qk1*p2k2*pp1-4.*qk1*qk2*pp1+2.*p1k1*qk2*pp2
	   -4.*qk1*qk2*pp2-2.*p1k1*p1k2*qp+2.*qk1*p1k2*qp-2.*p1k1*p2k2*qp
	   +2.*qk1*p2k2*qp+2.*p1k1*qk2*qp+2.*qk1*pk2*p1p2
	   +2.*qk1*p1k2*p1p2-2.*pk1*qk2*p1p2-4.*qk1*qk2*p1p2
	   +p2k1*(qk2*(-2.*m2+u2-2.*pp1+2.*qp)
		  -2.*pk2*(u2-qp1)-2.*p1k2*(u2-qp1))+2.*pk1*p1k2*qp1
	   +2.*pk1*p2k2*qp1+4.*qk1*qk2*qp1
	   -2.*p1k1*pk2*qp2+2.*qk1*pk2*qp2-2.*p1k1*p1k2*qp2+2.*qk1*p1k2*qp2
	   +2.*pk1*qk2*qp2+2.*p1k1*qk2*qp2
	   +k1k2*(-2.*m2*u2+2.*pp2*(u2-qp1)+2.*m2*qp1
		  +qp*(2.*m2+u2+2.*p1p2-2.*qp1-4.*qp2)+2.*m2*qp2+u2*qp2+2.*pp1*qp2
		  - 2.*qp1*qp2))/4.)
    + u2*m2*((2.*p2k1*pk2 + qk1*pk2 + 3.*pk1*p1k2 + 3.*p2k1*p1k2 + 2.*qk1*p1k2
	      + 2.*pk1*p2k2 + qk1*p2k2 - 3.*pk1*qk2 - 3.*p2k1*qk2
	      - p1k1*(pk2 + p2k2 + 2.*qk2)
	      + k1k2*(6.*m2 + 3.*pp1 - qp + 3.*p1p2 - qp2))/4.);

  tr33 = -(qk2*(p1k1*(m4 + m2*pp1 - m2*pp2 + pp22 - (2.*m2 + pp1 + pp2)*p1p2)
		+ p2k1*((m2+pp1)*(2.*m2-pp2) - (m2+pp1+2.*pp2)*p1p2)
		+ pk1*(m2*(m2 + pp1)
		       - (2.*m2 + pp1)*pp2 - (m2+pp2)*p1p2 + p1p22)));

  tr34 = m2*pk1*p2k2*qp - m2*p2k1*p2k2*qp - p1k1*p2k2*pp1*qp - p2k1*p2k2*pp1*qp
    + 2.*p1k1*p1k2*pp2*qp + p2k1*p1k2*pp2*qp - p1k1*qk2*pp2*qp - 2.*pk1*p1k2*qp*p1p2
    - p2k1*p1k2*qp*p1p2 + pk1*qk2*qp*p1p2 + 2.*p2k1*qk2*qp*p1p2
    + qk1*(-2.*m2*p2k2*(m2 + pp1) + m2*pk2*pp2 + m2*qk2*pp2 + qk2*pp1*pp2
	   + m2*qk2*p1p2 - pk2*pp1*p1p2 + qk2*pp1*p1p2
	   + p1k2*(-(pp1*pp2) + m2*p1p2))
    + m2*p1k1*p2k2*qp1 - m2*p2k1*p2k2*qp1 - pk1*p2k2*pp1*qp1
    - p2k1*p2k2*pp1*qp1 - 2.*p1k1*pk2*pp2*qp1 - p2k1*pk2*pp2*qp1 + p1k1*qk2*pp2*qp1
    + 2.*p2k1*qk2*pp2*qp1 + 2.*pk1*pk2*p1p2*qp1 + p2k1*pk2*p1p2*qp1 - pk1*qk2*p1p2*qp1
    - m2*pk1*pk2*qp2 + m2*p2k1*pk2*qp2 - m2*p1k1*p1k2*qp2 + m2*p2k1*p1k2*qp2
    - m2*pk1*qk2*qp2 - m2*p1k1*qk2*qp2 - 2.*m2*p2k1*qk2*qp2 + p1k1*pk2*pp1*qp2
    + p2k1*pk2*pp1*qp2 + pk1*p1k2*pp1*qp2
    + p2k1*p1k2*pp1*qp2 - pk1*qk2*pp1*qp2 - p1k1*qk2*pp1*qp2 - 2.*p2k1*qk2*pp1*qp2
    + k1k2*(p1p2*(pp1*qp - m2*qp1) + pp2*(-(m2*qp) + pp1*qp1) + 2.*m2*(m2 + pp1)*qp2);

  tr44 = -(p2k1*(-(pk2*(pp1*(u2+qp) + m2*(u2+2.*qp) + (m2-qp)*qp1+qp12))
		 - p1k2*(m2*u2 + m2*qp	+ qps + (2.*m2-qp)*qp1 + pp1*(u2+qp1))
		 + qk2*(m2*(m2+u2+qp) + (m2-2.*qp)*qp1 + pp1*(m2+u2+qp+qp1))));

  /*
    tr11=0.0;
    tr12=0.0;
    tr13=0.0;
    tr14=0.0;

    tr22=0.0;
    tr23=0.0;
    tr24=0.0;

    tr33=0.0;
    tr34=0.0;

    tr44=0.0;
  */

  matr2e = C1*C1*D1*D1*tr11 - C1*C1*D1*D2*tr13 + C1*C1*D2*D2*tr33;
  matr2mu = C2*C2*D1*D1*tr22 - C2*C3*D1*D2*tr24 + C3*C3*D2*D2*tr44;
  matr2emu = C1*C2*D1*D1*tr12 - C1*C3*D1*D2*tr14 - C1*C2*D1*D2*tr23
    + C1*C3*D2*D2*tr34;

  matr2 = C*(matr2e + matr2mu + matr2emu);


  //  printf("matr2=%G %G %G %G\n",matr2,matr2e, matr2mu, matr2emu);
  return matr2;
}
