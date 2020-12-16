//
// ------------------------------------------------------------
//      Mu3eMuonRadiativeDecayChannelWithSpin.cpp
//
//      Created on: Oct 10, 2014
//      Author:     aperrevoort
//
//               References:
//                    TRIUMF/TWIST Technote TN-55:
//                    "Radiative muon decay" by P. Depommier and A. Vacheret
//                    ------------------------------------------------------
//                    Yoshitaka Kuno and Yasuhiro Okada
//                    "Muon Decays and Physics Beyond the Standard Model"
//                    Rev. Mod. Phys. 73, 151 (2001)
//
// ------------------------------------------------------------
//
//
#include "Mu3eMuonRadiativeDecayChannelWithSpin.hh"


#include <G4MuonPlus.hh>
#include <G4Electron.hh>
#include <G4ParticleDefinition.hh>
#include <G4DecayProducts.hh>
#include <G4VDecayChannel.hh>
#include <G4LorentzVector.hh>


#include "rand.hh"

Mu3eMuonRadiativeDecayChannelWithSpin::Mu3eMuonRadiativeDecayChannelWithSpin(
    const G4String& theParentName,
    double          the10MeVBR,
    double          theRadGammaCut
)
  : G4VDecayChannel("rmd", 1)
  , parent_polarization()
{

  scale = ScaleBR(theRadGammaCut);
  gammacut = theRadGammaCut;
  // set names for daughter particles
  G4cout << "**** Mu3eRadiativeMuonDecayChannel::constructor :";
  if (theParentName == "mu+") {
    SetBR(the10MeVBR*scale);
    SetParent("mu+");
    SetNumberOfDaughters(4);
    SetDaughter(0, "e+");
    SetDaughter(1, "gamma");
    SetDaughter(2, "nu_e");
    SetDaughter(3, "anti_nu_mu");
  } else if (theParentName == "mu-") {
    SetBR(the10MeVBR*scale);
    SetParent("mu-");
    SetNumberOfDaughters(4);
    SetDaughter(0, "e-");
    SetDaughter(1, "gamma");
    SetDaughter(2, "anti_nu_e");
    SetDaughter(3, "nu_mu");
  } else {
    if (1) {
      G4cout << "**** Mu3eRadiativeMuonDecayChannel::constructor error:  parent particle is not muon but " << theParentName << G4endl;
    }
  }
}

G4DecayProducts* Mu3eMuonRadiativeDecayChannelWithSpin::DecayIt(double) {
  //UL  auto products = getEmptyProducts();
  //create parent G4DynamicParticle at rest
  G4ThreeVector dummy;
  G4DynamicParticle * parentparticle = new G4DynamicParticle( G4MT_parent, dummy, 0.0);
  //create G4Decayproducts
  G4DecayProducts *products = new G4DecayProducts(*parentparticle);
  delete parentparticle;

  if (1)  G4cout << "**** Mu3eMuonRadiativeDecayChannelWithSpin::DecayIt ***";

  // parent mass
  double parentmass = G4MT_parent->GetPDGMass();

  double EMMU = parentmass;

  //daughters'mass
  double daughtermass[4];
  double sumofdaughtermass = 0.0;
  for (int index=0; index<4; index++){
    daughtermass[index] = G4MT_daughters[index]->GetPDGMass();
    sumofdaughtermass += daughtermass[index];
  }

  double EMASS = daughtermass[0];

  int i = 0;

  double eps = EMASS/EMMU;

  double som0, Qsqr, x, y, u;
  double cthetaE, cthetaG, cthetaGE;
  double Gmin = gammacut;
  double ymin = 2.0*Gmin/EMMU;
  double umax = -1.0*log(ymin);

    G4ThreeVector uE, uG;

  do{
        i++;

        do{
//
//--------------------------------------------------------------------------
//      Build two vectors of random length and random direction, for the
//      positron and the photon.
//      x/y is the length of the vector, xx, yy and zz the components,
//      phi is the azimutal angle, theta the polar angle.
//--------------------------------------------------------------------------
            do{

//--------------------------------------------------------------------------
//
//      For the positron
//
                x = CLHEP::RandFlat::shoot();
//--------------------------------------------------------------------------
//
//      For the photon
//
                u = umax*CLHEP::RandFlat::shoot();
                y = exp(-u/1.0);
            }while(x<2.0*eps || y<ymin || x>1.0+eps*eps || y>1.0-eps*eps);		// check kinematic constraints
//--------------------------------------------------------------------------
//
//      For the positron
//
            uE = mu3e::util::rand_u3d();
            cthetaE = uE.z();
//
//          What you get:
//
//          x       = positron energy
//          phiE    = azimutal angle of positron momentum
//          cthetaE = cosine of polar angle of positron momentum
//          sthetaE = sine of polar angle of positron momentum
//
//-----------------------------------------------------------------------
//
//      For the photon
//
            uG = mu3e::util::rand_u3d();
            cthetaG = uG.z();
//
//          What you get:
//
//          y       = photon energy
//          phiG    = azimutal angle of photon momentum
//          cthetaG = cosine of polar angle of photon momentum
//          sthetaG = sine of polar angle of photon momentum
//
//-----------------------------------------------------------------------
//
//      Calculate the angle between positron and photon (cosine)
//
            cthetaGE = uE.cosTheta(uG);
//
//-----------------------------------------------------------------------
//
//      (Dimensionless) invariant mass of the neutrinos
//
            Qsqr = 1 + eps*eps - x - y + x*y/2.0*(1.0 - std::sqrt(1.0 - 4.0*eps*eps/(x*x))*cthetaGE);
            Qsqr  = Qsqr/((1.0-eps)*(1.0-eps));
//
//-----------------------------------------------------------------------
//
//      Check the kinematics.
//
        } while ( Qsqr<0.0 || Qsqr>1.0 );
//
//      Do the calculation for -1 muon polarization (i.e. mu+)
//
        double Pmu = -1.0;
        if (GetParentName() == "mu-")Pmu = +1.0;
//
//      and for differential BR (Kuno)
//
//-----------------------------------------------------------------------
//
        som0 = kuno(Pmu,x,y,cthetaE,cthetaG,cthetaGE);
//
//----------------------------------------------------------------------
//
//   Sample the decay rate
//

  } while (CLHEP::RandFlat::shoot()*1.0/y > som0);

//
//-----------------------------------------------------------------------
//
  double E = EMMU/2.*x;         // positron energy
  double G = EMMU/2.*y;         // gamma energy
//
//-----------------------------------------------------------------------
//

  if(E < EMASS) E = EMASS;

  // calculate daughter momentum
  double daughtermomentum[4];

  daughtermomentum[0] = std::sqrt(E*E - EMASS*EMASS);

  //Coordinates of the decay positron with respect to the muon spin

  G4ThreeVector direction0 = uE;

  direction0.rotateUz(parent_polarization);

  G4DynamicParticle * daughterparticle0 = new G4DynamicParticle( G4MT_daughters[0], daughtermomentum[0]*direction0);

  products->PushProducts(daughterparticle0);

  daughtermomentum[1] = G;

  //Coordinates of the decay gamma with respect to the muon spin

  G4ThreeVector direction1 = uG;

  direction1.rotateUz(parent_polarization);

  G4DynamicParticle * daughterparticle1 = new G4DynamicParticle( G4MT_daughters[1], daughtermomentum[1]*direction1);

  products->PushProducts(daughterparticle1);

  // daughter 3 ,4 (neutrinos)
  // create neutrinos in the C.M frame of two neutrinos

  double energy2 = parentmass*(1.0 - (x+y)/2.0);

  double vmass = energy2*energy2-(daughtermomentum[0]*direction0+daughtermomentum[1]*direction1)*(daughtermomentum[0]*direction0+daughtermomentum[1]*direction1);   // mass^2 of the 2 neutrinos
  if(vmass<0.0){ vmass = 0.0; }
  else{          vmass = std::sqrt(vmass);}
  double beta = (daughtermomentum[0]+daughtermomentum[1])/energy2;
  beta = -1.0 * std::min(beta,0.99);

    auto direction2 = mu3e::util::rand_u3d();

  G4DynamicParticle * daughterparticle2 = new G4DynamicParticle( G4MT_daughters[2], direction2*(vmass/2.));
  G4DynamicParticle * daughterparticle3 = new G4DynamicParticle( G4MT_daughters[3], direction2*(-1.0*vmass/2.));

  // boost to the muon rest frame

  G4ThreeVector direction34 = direction0 + direction1;
  direction34 = direction34.unit();

  G4LorentzVector p4 = daughterparticle2->Get4Momentum();
  p4.boost(direction34.x()*beta,direction34.y()*beta,direction34.z()*beta);
  daughterparticle2->Set4Momentum(p4);

  p4 = daughterparticle3->Get4Momentum();
  p4.boost(direction34.x()*beta,direction34.y()*beta,direction34.z()*beta);
  daughterparticle3->Set4Momentum(p4);

  products->PushProducts(daughterparticle2);
  products->PushProducts(daughterparticle3);

  daughtermomentum[2] = daughterparticle2->GetTotalMomentum();
  daughtermomentum[3] = daughterparticle3->GetTotalMomentum();

  return products;
}

double Mu3eMuonRadiativeDecayChannelWithSpin::kuno( double Pmu,
                                                    double x,
                                                    double y,
                                                    double cthetaE,
                                                    double cthetaG,
                                                    double cthetaGE)    // differential BR as in Kuno's paper
{
    double EMMU = 105.6583715;		// muon mass in MeV
    double EMASS = 0.51099906;		// electron mass in MeV
    double eps = EMASS/EMMU;
    double r = eps*eps;
    double beta = sqrt(1.0-4.0*r/(x*x));
    double d = 1.0-beta*cthetaGE;

    // double pi = 3.14159265358979323846;
    // double twopi = 2*pi;
    // double elm_coupling = 1.43996e-12;
    // double hbarc = 4.13566e-12/twopi * 299.792458;
    // double fine_structure_const = elm_coupling/hbarc;

    double F0 = 8.0/d*(y*y*(3.0-2.0*y)+6.0*x*y*(1.0-y)+2.0*x*x*(3.0-4.0*y)-4.0*x*x*x) + 8.0*(-x*y*(3.0-y-y*y)-x*x*(3.0-y-4.0*y*y)+2.0*x*x*x*(1.0+2.0*y)) + 2.0*d*(x*x*y*(6.0-5.0*y-2.0*y*y)-2.0*x*x*x*y*(4.0+3.0*y)) + 2.0*d*d*x*x*x*y*y*(2.0+y);
    double F1 = 32.0/(d*d)*(-y*(3.0-2.0*y)/x-(3.0-4.0*y)+2.0*x) + 8.0/d*(y*(6.0-5.0*y)-2.0*x*(4.0+y)+6.0*x*x) + 8.0*(x*(4.0-3.0*y+y*y)-3.0*x*x*(1.0+y)) + 6.0*d*x*x*y*(2.0+y);
    double F2 = 32.0/(d*d)*((4.0-3.0*y)/x-3.0) + 48.0*y/d;
    double F = F0 + r*F1 + r*r*F2;

    double G0 = 8.0/d*(x*y*(1.0-2.0*y)+2.0*x*x*(1.0-3.0*y)-4.0*x*x*x) + 4.0*(-x*x*(2.0-3.0*y-4.0*y*y)+2.0*x*x*x*(2.0+3.0*y)) - 4.0*d*x*x*x*y*(2.0+y);
    double G1 = 32.0/(d*d)*(-1.0+2.0*y+2.0*x) + 8.0/d*(-x*y+6.0*x*x) - 12.0*x*x*(2.0+y);
    double G2 = -96.0/(d*d);
    double G = G0 + r*G1 + r*r*G2;

    double H0 = 8.0/d*(y*y*(1.0-2.0*y)+x*y*(1.0-4.0*y)-2.0*x*x*y) + 4.0*(2.0*x*y*y*(1.0+y)-x*x*y*(1.0-4.0*y)+2.0*x*x*x*y) + 2.0*d*(x*x*y*y*(1.0-2.0*y)-4.0*x*x*x*y*y) + 2.0*d*d*x*x*x*y*y*y;
    double H1 = 32.0/(d*d)*(-y*(1.0-2.0*y)/x+2.0*y) + 8.0/d*(y*(2.0-5.0*y)-x*y) + 4.0*x*y*(2.0*y-3.0*x) + 6.0*d*x*x*y*y;
    double H2 = -96.0*y/(d*d*x) + 48.0*y/d;
    double H = H0 + r*H1 + r*r*H2;

    double kunoBR;

    kunoBR = F + beta*Pmu*cthetaE*G + Pmu*cthetaG*H;
    kunoBR = CLHEP::fine_structure_const/(64.0*M_PI*M_PI*M_PI)*beta/y*kunoBR;

    return kunoBR;
}

double Mu3eMuonRadiativeDecayChannelWithSpin::fron(double Pmu,
                                                   double x,
                                                   double y,
                                                   double cthetaE,
                                                   double cthetaG,
                                                   double cthetaGE)
{
      double mu  = 105.65;
      double me  =   0.511;
      double rho =   0.75;
      double del =   0.75;
      double eps =   0.0;
      double kap =   0.0;
      double ksi =   1.0;

      double delta = 1-cthetaGE;

//    Calculation of the functions f(x,y)

      double f_1s  = 12.0*((y*y)*(1.0-y)+x*y*(2.0-3.0*y)
                       +2.0*(x*x)*(1.0-2.0*y)-2.0*(x*x*x));
      double f0s   = 6.0*(-x*y*(2.0-3.0*(y*y))
                       -2.0*(x*x)*(1.0-y-3.0*(y*y))+2.0*(x*x*x)*(1.0+2.0*y));
      double f1s   = 3.0*((x*x)*y*(2.0-3.0*y-3.0*(y*y))
                       -(x*x*x)*y*(4.0+3.0*y));
      double f2s   = 1.5*((x*x*x)*(y*y)*(2.0+y));

      double f_1se = 12.0*(x*y*(1.0-y)+(x*x)*(2.0-3.0*y)
                       -2.0*(x*x*x));
      double f0se  = 6.0*(-(x*x)*(2.0-y-2.0*(y*y))
                       +(x*x*x)*(2.0+3.0*y));
      double f1se  = -3.0*(x*x*x)*y*(2.0+y);
      double f2se  = 0.0;

      double f_1sg = 12.0*((y*y)*(1.0-y)+x*y*(1.0-2.0*y)
                       -(x*x)*y);
      double f0sg  = 6.0*(-x*(y*y)*(2.0-3.0*y)-(x*x)*y*(1.0-4.0*y)
                       +(x*x*x)*y);
      double f1sg  = 3.0*((x*x)*(y*y)*(1.0-3.0*y)
                       -2.0*(x*x*x)*(y*y));
      double f2sg  = 1.5*(x*x*x)*(y*y*y);

      double f_1v  = 8.0*((y*y)*(3.0-2.0*y)+6.0*x*y*(1.0-y)
                       +2.0*(x*x)*(3.0-4.0*y)-4.0*(x*x*x));
      double f0v   = 8.0*(-x*y*(3.0-y-(y*y))-(x*x)*(3.0-y-4.0*(y*y))
                       +2.0*(x*x*x)*(1.0+2.0*y));
      double f1v   = 2.0*((x*x)*y*(6.0-5.0*y-2.0*(y*y))
                       -2.0*(x*x*x)*y*(4.0+3.0*y));
      double f2v   = 2.0*(x*x*x)*(y*y)*(2.0+y);

      double f_1ve = 8.0*(x*y*(1.0-2.0*y)
                       +2.0*(x*x)*(1.0-3.0*y)-4.0*(x*x*x));
      double f0ve  = 4.0*(-(x*x)*(2.0-3.0*y-4.0*(y*y))
                       +2.0*(x*x*x)*(2.0+3.0*y));
      double f1ve  = -4.0*(x*x*x)*y*(2.0+y);
      double f2ve  = 0.0;

      double f_1vg = 8.0*((y*y)*(1.0-2.0*y)+x*y*(1.0-4.0*y)
                       -2.0*(x*x)*y);
      double f0vg  = 4.0*(2.0*x*(y*y)*(1.0+y)-(x*x)*y*(1.0-4.0*y)
                       +2.0*(x*x*x)*y);
      double f1vg  = 2.0*((x*x)*(y*y)*(1.0-2.0*y)
                       -4.0*(x*x*x)*(y*y));
      double f2vg  = 2.0*(x*x*x)*(y*y*y);

      double f_1t  = 8.0*((y*y)*(3.0-y)+3.0*x*y*(2.0-y)
                       +2.0*(x*x)*(3.0-2.0*y)-2.0*(x*x*x));
      double f0t   = 4.0*(-x*y*(6.0+(y*y))
                       -2.0*(x*x)*(3.0+y-3.0*(y*y))+2.0*(x*x*x)*(1.0+2.0*y));
      double f1t   = 2.0*((x*x)*y*(6.0-5.0*y+(y*y))
                       -(x*x*x)*y*(4.0+3.0*y));
      double f2t   = (x*x*x)*(y*y)*(2.0+y);

      double f_1te = -8.0*(x*y*(1.0+3.0*y)+(x*x)*(2.0+3.0*y)
                       +2.0*(x*x*x));
      double f0te  = 4.0*((x*x)*(2.0+3.0*y+4.0*(y*y))
                       +(x*x*x)*(2.0+3.0*y));
      double f1te  = -2.0*(x*x*x)*y*(2.0+y);
      double f2te  = 0.0;

      double f_1tg = -8.0*((y*y)*(1.0+y)+x*y+(x*x)*y);
      double f0tg  = 4.0*(x*(y*y)*(2.0-y)+(x*x)*y*(1.0+2.0*y)
                       +(x*x*x)*y);
      double f1tg  = -2.0*((x*x)*(y*y)*(1.0-y)+2.0*(x*x*x)*y);
      double f2tg  = (x*x*x)*(y*y*y);

      double term = delta+2.0*(me*me)/((mu*mu)*(x*x));
      term = 1.0/term;

      double nss = term*f_1s+f0s+delta*f1s+(delta*delta)*f2s;
      double nv = term*f_1v+f0v+delta*f1v+(delta*delta)*f2v;
      double nt = term*f_1t+f0t+delta*f1t+(delta*delta)*f2t;

      double nse = term*f_1se+f0se+delta*f1se+(delta*delta)*f2se;
      double nve = term*f_1ve+f0ve+delta*f1ve+(delta*delta)*f2ve;
      double nte = term*f_1te+f0te+delta*f1te+(delta*delta)*f2te;

      double nsg = term*f_1sg+f0sg+delta*f1sg+(delta*delta)*f2sg;
      double nvg = term*f_1vg+f0vg+delta*f1vg+(delta*delta)*f2vg;
      double ntg = term*f_1tg+f0tg+delta*f1tg+(delta*delta)*f2tg;

      double term1 = nv;
      double term2 = 2.0*nss+nv-nt;
      double term3 = 2.0*nss-2.0*nv+nt;

      double term1e = 1.0/3.0*(1.0-4.0/3.0*del);
      double term2e = 2.0*nse+5.0*nve-nte;
      double term3e = 2.0*nse-2.0*nve+nte;

      double term1g = 1.0/3.0*(1.0-4.0/3.0*del);
      double term2g = 2.0*nsg+5.0*nvg-ntg;
      double term3g = 2.0*nsg-2.0*nvg+ntg;

      double som00 = term1+(1.0-4.0/3.0*rho)*term2+eps*term3;
      double som01 = Pmu*ksi*(cthetaE*(nve-term1e*term2e+kap*term3e)
                       +cthetaG*(nvg-term1g*term2g+kap*term3g));

      double som0 = (som00+som01)/y;
      som0  = CLHEP::fine_structure_const/(64*M_PI*M_PI*M_PI)*som0;

//      G4cout << x     << " " << y    << " " << som00 << " "
//             << som01 << " " << som0 << G4endl;

      return som0;
}

double Mu3eMuonRadiativeDecayChannelWithSpin::ScaleBR(double theRadGammaCut)
{
      double Pmu = -1.0;	// Muon polarisation, -1 for mu^+
      //if (GetParentName() == "mu-")Pmu = +1.0;

      double x;			// Electron momentum
      double y;			// Gamma momentum
      double cthetaE;
      double cthetaG;
      double cthetaGE;
      double u, umax;

      // double EMMU = 105.6583715;		// muon mass in MeV
      // double EMASS = 0.51099906;		// electron mass in MeV

      double EMMU = G4MuonPlus::MuonPlusDefinition()->GetPDGMass();
      double EMASS = G4Electron::ElectronDefinition()->GetPDGMass();
      double eps = EMASS/EMMU;

      double Gmin = theRadGammaCut;					// cut on photon energy in MeV
      double ycut_cut = 2.0*Gmin/EMMU;
      double ycut_10 = 2.0*10.0/EMMU;           // cut on photon energy at 10 MeV (for which the BR is known)

      int nint = 1e6;                           // number of interations for MC integration

      double som0 = 0;
      double q_sqr = 0.0;

      double som0sum_10 = 0.0;
      double som0sqrsum_10 = 0.0;
      double som0sum_10_3 = 0.0;
      //double dsom0_10_3 = 0.0;
      //double dsom0MC_10 = 0.0;
      //double dsom0V_10 = 0.0;

      double som0sum_cut = 0.0;
      double som0sqrsum_cut = 0.0;
      double som0sum_cut_3 = 0.0;
      //double dsom0_cut_3 = 0.0;
      //double dsom0MC_cut = 0.0;
      //double dsom0V_cut = 0.0;


      long Nvol_10 = 0;
      long Nvol_cut = 0;
      long Nvol = 0;

      double BRscale = 0.0;
      double y_min = 0.0;

      if(ycut_cut<=ycut_10)
      {
        umax = -1.0*log(ycut_cut);
        y_min = ycut_cut;
      }
      else
      {
        umax = -1.0*log(ycut_10);
        y_min = ycut_10;
      }


    // Sum of BR's
    for(int i = 0; i<nint; i++)
	{
// 	  do{
	    do
	    {
	      do{
            x=CLHEP::RandFlat::shoot();
            u=umax*CLHEP::RandFlat::shoot();
            y=exp(-u/1.0);
			Nvol++;
          }while(x<2.0*eps || y<y_min || x>1.0+eps*eps || y>1.0-eps*eps);

            auto unitE = mu3e::util::rand_u3d();
            auto unitG = mu3e::util::rand_u3d();

            cthetaE = unitE.z;
            cthetaG = unitG.z;
            cthetaGE = dot(unitE, unitG);

          q_sqr = 1 + eps*eps - x - y + x*y/2.0*(1.0 - std::sqrt(1.0 - 4.0*eps*eps/(x*x))*cthetaGE);
          q_sqr  = q_sqr/((1.0-eps)*(1.0-eps));
	    }while(q_sqr < 0.0 || q_sqr > 1.0 );
        som0 = kuno(Pmu, x, y, cthetaE, cthetaG, cthetaGE);
// 	  }while(CLHEP::RandFlat::shoot()*1.0/y > som0);

      if(som0 < 0.0)
	  {
        som0 = 0.0;
	  }

	  if(y>=ycut_cut)
	  {
        som0sum_cut = som0sum_cut + som0;
        som0sqrsum_cut = som0sqrsum_cut + som0*som0;
	    Nvol_cut++;
	  }

	  if(y>=ycut_10)
	  {
        som0sum_10 = som0sum_10 + som0;
        som0sqrsum_10 = som0sqrsum_10 + som0*som0;
	    Nvol_10++;
	  }

	}

    // MC integral of the BR's for both cuts
	som0sum_cut_3 = 64.0*som0sum_cut/Nvol;
	som0sum_10_3 = 64.0*som0sum_10/Nvol;

//	dsom0MC_cut = som0sqrsum_cut - som0sum_cut*som0sum_cut/Nvol_cut;
//	dsom0MC_10 = som0sqrsum_10 - som0sum_10*som0sum_10/Nvol_10;
 //   dsom0V_cut = (double)(Nvol - Nvol_cut)/(Nvol*Nvol_cut)*som0sum_cut*som0sum_cut;
 //   dsom0V_10 = (double)(Nvol - Nvol_10)/(Nvol*Nvol_10)*som0sum_10*som0sum_10;
//	dsom0_cut_3 = 64.0/Nvol*sqrt(dsom0MC_cut + dsom0V_cut*dsom0V_cut);
//	dsom0_10_3 = 64.0/Nvol*sqrt(dsom0MC_10 + dsom0V_10*dsom0V_10);

    // Get scaling factor with respect to 10MeV cut
	BRscale = som0sum_cut_3/som0sum_10_3;
    //dBRScale = sqrt(dsom0_cut_3/som0sum_cut_3*dsom0_cut_3/som0sum_cut_3 + dsom0_10_3/som0sum_10_3*dsom0_10_3/som0sum_10_3)*som0sum_cut_3/som0sum_10_3;

	return BRscale;
}
