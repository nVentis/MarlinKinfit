/*! \file 
 *  \brief Implements class ParticleConstraint
 *
 * \b Changelog:
 * - 
 *
 * \b CVS Log messages:
 * - $Log: ParticleConstraint.cc,v $
 * - Revision 1.2  2008/10/17 13:17:16  blist
 * - Avoid variable-size arrays
 * -
 * - Revision 1.1  2008/02/12 10:19:09  blist
 * - First version of MarlinKinfit
 * -
 * - Revision 1.2  2008/02/07 08:21:07  blist
 * - ParticleConstraint.C fixed
 * -
 * - Revision 1.1  2008/02/07 08:18:57  blist
 * - ParticleConstraint,C added
 * -
 */ 

#include "ParticleConstraint.h"
#include "ParticleFitObject.h"
#include <iostream>
#include <cmath>
using namespace std;

// Returns the four momentum of the fitobjects system in order (E,Px,Py,Pz)
double* ParticleConstraint::getFourMomentum(int flag) {
    double* result = new double[4]{0};

    for (unsigned int i = 0; i < fitobjects.size(); i++) {
      if (flags[i] == flag) {
        const ParticleFitObject *fok = dynamic_cast < ParticleFitObject* > ( fitobjects[i] );
        assert(fok);

        //cerr << "(E,Px,Py,Pz)=(" << fok->getE() << "," << fok->getPx() << "," << fok->getPy() << "," << fok->getPz() << ")" << endl;

        result[0] += fok->getE(); 
        result[1] += fok->getPx(); 
        result[2] += fok->getPy(); 
        result[3] += fok->getPz();
      }
    }

    //cerr << "RESULT (E,Px,Py,Pz)=(" << result[0] << "," << result[1] << "," << result[2] << "," << result[3] << ")" << endl;
    
    return result;
  }


// probably these can also be moved to basehardconstraint?

