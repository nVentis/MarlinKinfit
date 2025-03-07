/*! \file 
 *  \brief Declares class ParticleConstraint
 *
 * \b Changelog:
 * - 17.11.04 BL: First version (refactured from BaseConstraint)
 *
 * \b CVS Log messages:
 * - $Log: ParticleConstraint.h,v $
 * - Revision 1.2  2008/02/12 16:43:26  blist
 * - First Version of Soft Constraints
 * -
 * - Revision 1.1  2008/02/12 10:19:06  blist
 * - First version of MarlinKinfit
 * -
 * - Revision 1.8  2008/02/04 17:30:54  blist
 * - NewtonFitter works now!
 * -
 * - Revision 1.7  2008/01/30 21:48:03  blist
 * - Newton Fitter still doesnt work :-(
 * -
 * - Revision 1.6  2008/01/30 09:14:54  blist
 * - Preparations for NewtonFitter
 * -
 * - Revision 1.5  2008/01/29 17:21:12  blist
 * - added methods secondDerivatives and firstDerivatives
 * -
 * - Revision 1.4  2007/09/14 10:58:42  blist
 * - Better documentation,
 * - added PyConstraint::add1stDerivativesToMatrix,
 * - fixed bug in PzConstraint::add1stDerivativesToMatrix
 * -
 * - Revision 1.3  2007/09/13 13:09:04  blist
 * - Better docs for ParticleConstraint
 * -
 * - Revision 1.2  2007/09/13 08:09:51  blist
 * - Updated 2nd derivatives for px,py,pz,E constraints, improved header documentation
 * -
 *
 */ 

#ifndef __PARTICLECONSTRAINT_H
#define __PARTICLECONSTRAINT_H

#include<vector>
#include<cassert>
#include "BaseHardConstraint.h"

class ParticleFitObject;

//  Class ParticleConstraint:
/// Abstract base class for constraints of kinematic fits
/**
 * This class defines the minimal functionality any constraint class must provide. 
 * First of all a constraint should know on with particles (or FitObject) it is applied. 
 * Where as for example a constraint on the total transvese momentum takes into
 * account all particles in the event, an invariant mass constraint usually applies only
 * to a subset of particles. 
 *
 * The particle list is implemented as a vector containing pointers to objects derived 
 * from ParticleFitObject and can be either set a whole (setFOList) or enlarged by adding
 * a single ParticleFitObject (addToFOList).
 *
 * From the four--momenta of all concerned fit objects the constraint has to be able 
 * to calculate its current value (getValue). Constraints should be formulated such that 
 * a value of zero corresponds to a perfectly fulfilled constraint.
 *
 * In order to find a solution to the constrained minimisation problem, fit algorithms
 * usually need the first order derivatives of the constraint with respect to the fit
 * parameters. Since many constraints can be most easily expressed in terms of E, px, py, pz,
 * the constraints supply their derivatives w.r.t. these parameters. If a FitObject uses
 * a different parametrisation, it is its own task to provide the additional derivatives
 * of  E, px, py, pz w.r.t. the parameters of the FitObject. Thus it is easily possible
 * to use FitObjects with different kinds of parametrisations under the same constraint.
 * Some fit algorithms also need the second derivatives of the constraint, 
 * i.e. the NewtonFitter. 
 *
 * First and second order derivatives of each constraint can be added directly to the 
 * global covariance matrix containing the derivatives of all constraints w.r.t. to all
 * parameters (add1stDerivativesToMatrix, add2ndDerivativesToMatrix). This requires the
 * constraint to know its position in the overall list of constraints (globalNum). 
 * 
 *
 * Author: Jenny List, Benno List
 * $Date: 2008/02/12 16:43:26 $
 * $Author: blist $
 *
 */

class ParticleConstraint: public BaseHardConstraint {
  public:
    /// Creates an empty ParticleConstraint object
    inline ParticleConstraint();
    /// Virtual destructor
    virtual ~ParticleConstraint() {};
    
    /// Adds several ParticleFitObject objects to the list
    virtual void setFOList(std::vector <ParticleFitObject*> *fitobjects_ ///< A list of BaseFitObject objects
                          ){
      for (int i = 0; i < (int) fitobjects_->size(); i++) {
	fitobjects.push_back (  reinterpret_cast < BaseFitObject* >  ( (*fitobjects_)[i] ) );
        flags.push_back (1);
      }  
    }; 
    /// Adds one ParticleFitObject objects to the list
    virtual void addToFOList(ParticleFitObject& fitobject, int flag = 1
                             ){
      fitobjects.push_back ( reinterpret_cast < BaseFitObject* >  ( &fitobject ) );
      flags.push_back (flag);
    }; 
    /// Resests ParticleFitObject list
    virtual void resetFOList(){
      fitobjects.resize (0);
      flags.resize (0);
    };

    // Returns the four momentum of the fitobjects system in order (E,Px,Py,Pz). make sure to call delete[] on the result
    virtual void getFourMomentum(int flag, double (&fourList) [4]);

    /// Invalidates any cached values for the next event
    virtual void invalidateCache() const 
    {}
      
  protected:

};

ParticleConstraint::ParticleConstraint() 
// : fitobjects( FitObjectContainer() ), derivatives( std::vector <double> () ), flags( std::vector <int> () ), globalNum(-999)
{
  invalidateCache();
}

#endif // __PARTICLECONSTRAINT_H
