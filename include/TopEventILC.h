////////////////////////////////////////////////////////////////
// Class TopEventILC
//
// Author: Benno List, Jenny Boehme
// Last update: $Date: 2008/02/12 16:43:26 $
//          by: $Author: blist $
// 
// Description: class to generate and fit top pair events at ILC
//               
////////////////////////////////////////////////////////////////
#ifdef MARLIN_USE_ROOT

#ifndef __TOPEVENTILC_H
#define __TOPEVENTILC_H

#include "BaseEvent.h"
#include "JetFitObject.h"
#include "PConstraint.h"
// #include "PxConstraint.h"
// #include "PyConstraint.h"
#include "MassConstraint.h"

class TopEventILC : public BaseEvent {
  public: 
    TopEventILC();
    virtual ~TopEventILC();
    virtual void genEvent();
    virtual int fitEvent (BaseFitter& fitter);

    double bwrandom (double r, double e0, double gamma, double emin, double emax) const;
    
    BaseConstraint& getPxConstraint() {return pxc;};
    BaseConstraint& getPyConstraint() {return pyc;};
    BaseConstraint& getW1Constraint() {return w1;};
    BaseConstraint& getW2Constraint() {return w2;};
    BaseConstraint& getTopConstraint() {return w;};
    
    double getW1Mass()  {return w1.getMass();};
    double getW2Mass()  {return w2.getMass();};
    double getTopMass(int flag)  {return w.getMass(flag);};
    
    ParticleFitObject* getTrueFitObject (int i) {return bfo[i];};
    ParticleFitObject* getSmearedFitObject (int i) {return bfosmear[i];};
    FourVector* getTrueFourVector (int i) {return fv[i];};
    
    bool leptonic;
    
  protected:
    enum {NFV = 11, NBFO = 6};
    FourVector *fv[NFV];
    FourVector *fvsmear[NFV];
    FourVector *fvfinal[NFV];
    ParticleFitObject *bfo[NBFO];
    ParticleFitObject *bfosmear[NBFO];
    
    PConstraint pxc;
    PConstraint pyc;
    PConstraint pzc;
    PConstraint ec;
    MassConstraint w1;
    MassConstraint w2;
    MassConstraint w;
    
    

};


#endif // __TOPEVENTILC_H

#endif // MARLIN_USE_ROOT
