#ifndef REVERSECASCADED2D
#define REVERSECASCADED2D

#include "petsc.h"
#include "core/globalDefinitions.hh"
#include "adjointDynamics/reverseIsothermalMacros2d.hh"
#include "external/tapenadeHeader.hh"

template <class Lattice, class Real = PetscScalar>
class ReverseIncompressibleCascaded2d {};

template <class Real>
class ReverseIncompressibleCascaded2d<D2Q9,Real> {

public:
  using RevMacros = ReverseIncompressibleMacros2d<D2Q9,Real>;
  using Macros = IncompressibleMacros2d<D2Q9,Real>;
  using RealType = Real;

  ReverseIncompressibleCascaded2d(IncompressibleFlowParameters par)
  {
    Real nu = par.velocityChar * par.lengthChar / par.ReynoldsNumber;
    omega = 1./(3.*nu + 0.5);
  }
  ReverseIncompressibleCascaded2d(Real _om) : omega(_om){}

  void operator()(const Real* fdistIn, Real* fdistOut,
                  Real* mac, Real* fdistInRev,
                  Real* fdistOutRev, Real* macRev) const
  {
    Real k[D2Q9::numDOF - 3];
    Real kRev[D2Q9::numDOF - 3];
    Real tempb,tempb0,tempb1,tempb2,tempb3,tempb4,tempb5,tempb6,tempb7,tempb8;
    Real tempb9,tempb10,tempb11,tempb12,tempb13,tempb14,tempb15,tempb16,tempb17;
    Real tempb18,tempb19,tempb20,tempb21,tempb22;
    Macros::compute(fdistIn,mac[0],mac[1],mac[2]);
    Real uSqr = mac[1]*mac[1] + mac[2]*mac[2];
    Real uSqrRev;
    k[0] = (mac[0]*uSqr-fdistIn[8]-fdistIn[6]-fdistIn[2]-fdistIn[4]-2.*(
        fdistIn[1]+fdistIn[3]+fdistIn[7]+fdistIn[5]-mac[0]/3.))/12.;
    k[1] = omega*(fdistIn[6]+fdistIn[2]-fdistIn[8]-fdistIn[4]+mac[0]*(mac[1]*
        mac[1]-mac[2]*mac[2]))/4.;
    k[2] = omega*(fdistIn[7]+fdistIn[3]-fdistIn[5]-fdistIn[1]-mac[0]*mac[1]*
        mac[2])/4.;
    k[3] = -((fdistIn[1]+fdistIn[3]-fdistIn[7]-fdistIn[5]-0.5*mac[0]*mac[1]*mac
        [1]*mac[2]+mac[2]*(mac[0]-fdistIn[6]-fdistIn[2]-fdistIn[0]))+0.5*
        mac[1]*(fdistIn[7]-fdistIn[5]-fdistIn[1]+fdistIn[3])) + 0.5*mac[2]*(-
        3.*k[0]-k[1]) + 2.*mac[1]*k[2];
    k[4] = -((fdistIn[3]+fdistIn[5]-fdistIn[1]-fdistIn[7]-0.5*mac[0]*mac[2]*mac
        [2]*mac[1]+mac[1]*(mac[0]-fdistIn[4]-fdistIn[8]-fdistIn[0]))+0.5*
        mac[2]*(fdistIn[7]+fdistIn[3]-fdistIn[1]-fdistIn[5])) + 0.5*mac[1]*(-3.
        *k[0]+k[1]) + 2.*mac[2]*k[2];
    reverseMomentToDistribution(fdistInRev,fdistOutRev,kRev);
    tempb = 0.25*kRev[5];
    tempb0 = 2.*tempb;
    tempb1 = mac[1]*tempb0;
    tempb2 = mac[2]*tempb0;
    tempb3 = 4.*(fdistIn[5]-fdistIn[7]+fdistIn[1]-fdistIn[3])*tempb;
    tempb4 = 4.*mac[1]*mac[2]*tempb;
    tempb5 = mac[2]*mac[2]*tempb;
    tempb6 = -(mac[1]*mac[1]*tempb);
    tempb7 = 4.*k[2]*kRev[5];
    tempb8 = 0.5*k[1]*kRev[5];
    macRev[0] = macRev[0] + 3.*(mac[1]*mac[1])*tempb5 + tempb/9.;
    macRev[1] = macRev[1] + 2*mac[1]*tempb8 + mac[2]*tempb7 - 2.*k[4]*kRev[5] - (
        fdistIn[6]+fdistIn[7]+fdistIn[5]+fdistIn[2]+fdistIn[1]+fdistIn[3])*2*
        mac[1]*tempb + mac[0]*3.*2*mac[1]*tempb5 + mac[2]*tempb3 + (fdistIn[7]
        -fdistIn[5]+fdistIn[1]-fdistIn[3])*tempb0;
    macRev[2] = macRev[2] + mac[1]*tempb7 - 2*mac[2]*tempb8 - 2.*k[3]*kRev[5] + (3.*
        (mac[0]*(mac[1]*mac[1]))-fdistIn[8]-fdistIn[7]-fdistIn[5]-fdistIn[1]-
        fdistIn[3]-fdistIn[4])*2*mac[2]*tempb + mac[1]*tempb3 + (fdistIn[7]+
        fdistIn[5]-fdistIn[1]-fdistIn[3])*tempb0;
    kRev[0] = kRev[0] + (-2.-1.5*uSqr)*kRev[5];
    kRev[4] = kRev[4] - 2.*mac[1]*kRev[5];
    kRev[3] = kRev[3] - 2.*mac[2]*kRev[5];
    kRev[2] = kRev[2] + 4.*mac[1]*mac[2]*kRev[5];
    uSqrRev = -(1.5*k[0]*kRev[5]);
    kRev[1] = kRev[1] + 0.5*(mac[1]*mac[1]-mac[2]*mac[2])*kRev[5];
    tempb9 = 0.5*mac[1]*kRev[4];
    tempb10 = -(0.5*mac[2]*kRev[4]);
    tempb11 = -(kRev[4]/4.);
    tempb12 = mac[1]*tempb11;
    tempb13 = -(2.*(mac[2]*mac[2])*tempb11);
    macRev[1] = macRev[1] + mac[0]*tempb13 + (mac[0]-fdistIn[4]-fdistIn[8]-fdistIn
        [0])*tempb11 + 0.5*(k[1]-3*k[0])*kRev[4];
    kRev[1] = kRev[1] + tempb9;
    kRev[0] = kRev[0] - 3*tempb9;
    macRev[2] = macRev[2] + (2.*k[2]-0.5*(fdistIn[7]+fdistIn[3]-fdistIn[1]-fdistIn
        [5]))*kRev[4] - mac[0]*mac[1]*2.*2*mac[2]*tempb11;
    macRev[0] = macRev[0] + mac[1]*tempb13 + tempb12;
    kRev[2] = kRev[2] + 2.*mac[2]*kRev[4];
    tempb14 = 0.5*mac[2]*kRev[3];
    tempb15 = -(0.5*mac[1]*kRev[3]);
    tempb16 = -(kRev[3]/4.);
    tempb17 = mac[2]*tempb16;
    tempb18 = -(2.*(mac[1]*mac[1])*tempb16);
    macRev[2] = macRev[2] + mac[0]*tempb18 + (mac[0]-fdistIn[6]-fdistIn[2]-fdistIn
        [0])*tempb16 + 0.5*(-(3.*k[0])-k[1])*kRev[3];
    kRev[0] = kRev[0] - 3.*tempb14;
    kRev[1] = kRev[1] - tempb14;
    macRev[1] = macRev[1] + (2.*k[2]-0.5*(fdistIn[7]-fdistIn[5]+fdistIn[3]-fdistIn
        [1]))*kRev[3] - mac[0]*mac[2]*2.*2*mac[1]*tempb16;
    macRev[0] = macRev[0] + mac[2]*tempb18 + tempb17;
    kRev[2] = kRev[2] + 2.*mac[1]*kRev[3];
    tempb19 = omega*kRev[2]/4.;
    macRev[0] = macRev[0] - mac[2]*mac[1]*tempb19;
    macRev[1] = macRev[1] - mac[2]*mac[0]*tempb19;
    macRev[2] = macRev[2] - mac[0]*mac[1]*tempb19;
    tempb20 = omega*kRev[1]/4.;
    macRev[0] = macRev[0] + (mac[1]*mac[1]-mac[2]*mac[2])*tempb20;
    macRev[1] = macRev[1] + mac[0]*2*mac[1]*tempb20;
    macRev[2] = macRev[2] - mac[0]*2*mac[2]*tempb20;
    tempb21 = kRev[0]/12.;
    tempb22 = -(2.*tempb21);
    macRev[0] = macRev[0] + uSqr*tempb21 - tempb22/3.;
    uSqrRev = uSqrRev + mac[0]*tempb21;
    fdistInRev[0] -= tempb12 + tempb17;
    fdistInRev[1] += tempb6 - tempb5 + tempb4 - tempb2 + tempb1 -
      tempb + tempb22 - tempb19 - tempb10 - tempb11 - tempb15 + tempb16;
    fdistInRev[2] += tempb6 - tempb17 + tempb20 - tempb21;
    fdistInRev[3] += tempb6 - tempb5 - tempb4 - tempb2 - tempb1 -
      tempb + tempb22 + tempb19 + tempb15 + tempb16 + tempb10 + tempb11;
    fdistInRev[4] -= tempb5 + tempb12 + tempb20 + tempb21;
    fdistInRev[5] += -tempb - tempb1 + tempb2 + tempb4 - tempb5 + tempb6
      + tempb11 - tempb10 - tempb15 - tempb16 - tempb19 + tempb22;
    fdistInRev[6] += tempb6 - tempb17 + tempb20 - tempb21;
    fdistInRev[7] += -tempb + tempb1 + tempb2 - tempb4 - tempb5 + tempb6
      + tempb10 - tempb11 + tempb15 - tempb16 + tempb19 + tempb22;
    fdistInRev[8] -= tempb5 + tempb12 + tempb20 + tempb21;
    macRev[1] = macRev[1] + 2*mac[1]*uSqrRev;
    macRev[2] = macRev[2] + 2*mac[2]*uSqrRev;
    RevMacros::compute(fdistInRev,macRev[0],macRev[1],macRev[2]);
  }
private:
  Real omega;
  void reverseMomentToDistribution(Real* fdistInRev, Real* fdistOutRev,
                                   Real* kRev) const
  {
    fdistInRev[0] += fdistOutRev[0];
    fdistInRev[1] += fdistOutRev[1];
    fdistInRev[2] += fdistOutRev[2];
    fdistInRev[3] += fdistOutRev[3];
    fdistInRev[4] += fdistOutRev[4];
    fdistInRev[5] += fdistOutRev[5];
    fdistInRev[6] += fdistOutRev[6];
    fdistInRev[7] += fdistOutRev[7];
    fdistInRev[8] += fdistOutRev[8];
    kRev[0] = -4.*fdistOutRev[0] + 2.*fdistOutRev[1] - fdistOutRev[2] + 2.*fdistOutRev[3]
      - fdistOutRev[4] + 2.*fdistOutRev[5] - fdistOutRev[6] + 2.*fdistOutRev[7]
      - fdistOutRev[8];
    kRev[1] = -fdistOutRev[2] + fdistOutRev[4] - fdistOutRev[6] + fdistOutRev[8];
    kRev[2] = fdistOutRev[1] - fdistOutRev[3] + fdistOutRev[5] - fdistOutRev[7];
    kRev[3] = fdistOutRev[1] - 2.*fdistOutRev[2] + fdistOutRev[3] - fdistOutRev[5]
      + 2.*fdistOutRev[6] - fdistOutRev[7];
    kRev[4] = -fdistOutRev[1] + fdistOutRev[3] - 2.*fdistOutRev[4] + fdistOutRev[5]
      - fdistOutRev[7] + 2.*fdistOutRev[8];
    kRev[5] = 4.*fdistOutRev[0] + fdistOutRev[1] - 2.*fdistOutRev[2] + fdistOutRev[3]
      - 2.*fdistOutRev[4] + fdistOutRev[5] - 2.*fdistOutRev[6] + fdistOutRev[7]
      - 2.*fdistOutRev[8];
  }
};

#endif
