/*
 *      Tohoku University's Numerical-Analysis Model
 *          for Investigation of tsunami
 *           Far-field Tsunami version
 *
 *Implementation flow
 *
 *Computation of the equation of continuity
 *Open sea, boundary condition, open boundaries
 *Computation of the equation of motion
 *Making area of computation be with in the area under consideration and
 *Enlargement of the area of computation as the tsunami propagates 
 */

#include <stdio.h>
#include <stdlib.h>
#include "UTILMDL.h"
#include "TUNAMI_FF.h"


#define Node(idx1, idx2) Node.node[idx1][idx2]
#define CNode CStructNode
#define gNode ((CStructNode*)gNode)
/*
FLOW of Implementation:
Check of the area of computation
Equation of continuity (iZ, iM, iN, iR1, iR2)
sea floor topography (mass conservation)
Open Sea boundary conditions (open bondary conditions) (iZ, iM, iN, C1, C2, C3, C4)
Equation of Motion (moment conservation)
Open Sea boundary (open bondaries) (iZ, iM, iN, iR2, iR4)
Check of the area of computation (calculation area for the next step)
*/

int MASSGBOUNDMOMENT( void )
{
 /*
VARIABLES
iZ
iM
iN
iZmax
iTime
Imin, Imax, Jmin, Jmax
zoutAT
NLat
iR1, iR2, iR4, R6
C1, C2, C3, C4
zoutCT
*/


 int i,j,enlarge;
  float absH,v1,v2;
  int m;

  CNode& Node = *gNode;

  // sea floor topography (mass conservation)
  #pragma omp parallel for default(shared) private(i,j,absH)
  for( i=Imin; i<=Imax; i++ ) {
    for( j=Jmin; j<=Jmax; j++ ) {

      m = idx(j,i);

      if( Node(m, iD) == 0 ) continue;

      Node(m, iZ) = Node(m, iZ) - Node(m, iR1)*( Node(m, iM) - Node(m-NLat, iM) + Node(m, iN)*R6[j] - Node(m-1, iN)*R6[j-1] );

      absH = fabs(Node(m, iZ));

      if( absH < MODELVARI.zoutZT ) Node(m, iZ) = 0.;

      if( Node(m, iZ) > Node(m, iZmax) ) Node(m, iZmax) = Node(m, iZ);

      if( MODELVARI.zoutAT && Node(m, iTime) < 0 && absH > MODELVARI.zoutAT ) Node(m, iTime) = (float)MODELVARI.time;

    }
  }

  // open bondary conditions
  if( Jmin <= 2 ) {
    for( i=2; i<=(NLon-1); i++ ) {
      m = idx(1,i);
      Node(m, iZ) = sqrt( pow(Node(m, iN),2.) + 0.25*pow((Node(m, iM)+Node(m-NLat, iM)),2.) )*C1[i];
      if( Node(m, iN) > 0 ) Node(m, iZ) = - Node(m, iZ);
    }
  }
  if( Imin <= 2 ) {
    for( j=2; j<=(NLat-1); j++ ) {
      m = idx(j,1);
      Node(m, iZ) = sqrt( pow(Node(m, iM),2.) + 0.25*pow((Node(m, iN)+Node(m-1, iN)),2.) )*C2[j];
      if( Node(m, iM) > 0 ) Node(m, iZ) = - Node(m, iZ);
    }
  }
  if( Jmax >= (NLat-1) ) {
    for( i=2; i<=(NLon-1); i++ ) {
      m = idx(NLat,i);
      Node(m, iZ) = sqrt( pow(Node(m-1, iN),2.) + 0.25*pow((Node(m, iM)+Node(m-1, iM)),2.) )*C3[i];
      if( Node(m-1, iN) < 0 ) Node(m, iZ) = - Node(m, iZ);
    }
  }
  if( Imax >= (NLon-1) ) {
    for( j=2; j<=(NLat-1); j++ ) {
    	m = idx(j,NLon);
      Node(m, iZ) = sqrt( pow(Node(m-NLat, iM),2.) + 0.25*pow((Node(m, iN)+Node(m-1, iN)),2.) )*C4[j];
      if( Node(m-NLat, iM) < 0 ) Node(m, iZ) = - Node(m, iZ);
    }
  }
  if( Jmin <= 2 ) {
	  m = idx(1,1);
    Node(m, iZ) = sqrt( pow(Node(m, iM),2.) + pow(Node(m, iN),2.) )*C1[1];
    if( Node(m, iN) > 0 ) Node(m, iZ) = - Node(m, iZ);
    m = idx(1,NLon);
    Node(m, iZ) = sqrt( pow(Node(m-NLat, iM),2.) + pow(Node(m, iN),2.) )*C1[NLon];
    if( Node(m, iN) > 0 ) Node(m, iZ) = - Node(m, iZ);
  }
  if( Jmin >= (NLat-1) ) {
	  m = idx(NLat,1);
    Node(m, iZ) = sqrt( pow(Node(m, iM),2.) + pow(Node(m-1, iN),2.) )*C3[1];
    if( Node(m-1, iN) < 0 ) Node(m, iZ) = - Node(m, iZ);
    m = idx(NLat,NLon);
    Node(m, iZ) = sqrt( pow(Node(m-NLat, iM),2.) + pow(Node(m-1, iN),2.) )*C3[NLon];
    if( Node(m-1, iN) < 0 ) Node(m, iZ) = - Node(m, iZ);
  }

  // )
  #pragma omp parallel for default(shared) private(i,j)
  for( i=Imin; i<=Imax; i++ ) {
    for( j=Jmin; j<=Jmax; j++ ) {

      m = idx(j,i);

      if( (Node(m, iD)*Node(m+NLat, iD)) != 0 )
        Node(m, iM) = Node(m, iM) - Node(m, iR2)*(Node(m+NLat, iZ)-Node(m, iZ));

      if( (Node(m, iD)*Node(m+1, iD)) != 0 )
        Node(m, iN) = Node(m, iN) - Node(m, iR4)*(Node(m+1, iZ)-Node(m, iZ));

    }
  }
  // open boundaries
  if( Jmin <= 2 ) {
    for( i=1; i<=(NLon-1); i++ ) {
    	m = idx(1,i);
      Node(m, iM) = Node(m, iM) - Node(m, iR2)*(Node(m+NLat, iZ) - Node(m, iZ));
    }
  }
  if( Imin <= 2 ) {
    for( j=1; j<=NLat; j++ ) {
    	m = idx(j,1);
      Node(m, iM) = Node(m, iM) - Node(m, iR2)*(Node(m+NLat, iZ) - Node(m, iZ));
    }
  }
  if( Jmax >= (NLat-1) ) {
    for( i=1; i<=(NLon-1); i++ ) {
      m = idx(NLat,i);
      Node(m, iM) = Node(m, iM) - Node(m, iR2)*(Node(m+NLat, iZ) - Node(m, iZ));
    }
  }
  if( Imin <= 2 ) {
    for( j=1; j<=(NLat-1); j++ ) {
      m = idx(j,1);
      Node(m, iN) = Node(m, iN) - Node(m, iR4)*(Node(m+1, iZ) - Node(m, iZ));
    }
  }
  if( Jmin <= 2 ) {
    for( i=1; i<=NLon; i++ ) {
      m = idx(1,i);
      Node(m, iN) = Node(m, iN) - Node(m, iR4)*(Node(m+1, iZ) - Node(m, iZ));
    }
  }
  if( Imax >= (NLon-1) ) {
    for( j=1; j<=(NLat-1); j++ ) {
      m = idx(j,NLon);
      Node(m, iN) = Node(m, iN) - Node(m, iR4)*(Node(m+1, iZ) - Node(m, iZ));
    }
  }

  // calculation area for the next step
  if( Imin > 2 ) {
    for( enlarge=0, j=Jmin; j<=Jmax; j++ ) {
      if( fabs(Node(idx(j,Imin+2), iZ)) > MODELVARI.zoutCT ) { enlarge = 1; break; }
    }
    if( enlarge ) { Imin--; if( Imin < 2 ) Imin = 2; }
  }
  if( Imax < (NLon-1) ) {
    for( enlarge=0, j=Jmin; j<=Jmax; j++ ) {
      if( fabs(Node(idx(j,Imax-2), iZ)) > MODELVARI.zoutCT ) { enlarge = 1; break; }
    }
    if( enlarge ) { Imax++; if( Imax > (NLon-1) ) Imax = NLon-1; }
  }
  if( Jmin > 2 ) {
    for( enlarge=0, i=Imin; i<=Imax; i++ ) {
      if( fabs(Node(idx(Jmin+2,i), iZ)) > MODELVARI.zoutCT ) { enlarge = 1; break; }
    }
    if( enlarge ) { Jmin--; if( Jmin < 2 ) Jmin = 2; }
  }
  if( Jmax < (NLat-1) ) {
    for( enlarge=0, i=Imin; i<=Imax; i++ ) {
      if( fabs(Node(idx(Jmax-2,i), iZ)) > MODELVARI.zoutCT ) { enlarge = 1; break; }
    }
    if( enlarge ) { Jmax++; if( Jmax > (NLat-1) ) Jmax = NLat-1; }
  }

return 0;
}


/*
FLOW of Implementation:
Check of the area of computation
Equation of continuity (iZ, iM, iN, iR1, R6)
sea floor topography (mass conservation) (iZ, iZmax)
Open Sea boundary conditions (open bondary conditions)  (iZ, iM, iN, C1, C2, C3, C4)
Equation of Motion (longitudial flux update - moment conservation) (v1, v2, iZ, iM, iN, iR2, iR3)
Open Sea boundary (open bondaries) (iZ, iM, iR2)
Equation of Motion (lattitudial flux update -moment conservation) (v1, v2, iZ, iM, iN, iR4, iR5)
Open Sea boundary (open bondaries) (iZ, iM, iN, iR4)
Check of the area of computation (calculation area for the next step)
*/

int MASSGBOUNDMOMENTCor( void )
{
 /*
VARIABLES
iZ
Discharge In I-direction	iM
Discharge In J-direction	iN
iZmax
iTime
Imin, Imax, Jmin, Jmax
zoutZT
NLat
iR1, iR2, iR3, iR4, iR5, R6
C1, C2, C3, C4
zoutAT
*/
  int i,j,enlarge;
  float absH,v1,v2;
  int m;

  CNode& Node = *gNode;

  // sea floor topography (mass conservation)
  #pragma omp parallel for default(shared) private(i,j,m,absH)
  for( i=Imin; i<=Imax; i++ ) {
    for( j=Jmin; j<=Jmax; j++ ) {

      m = idx(j,i);

      if( Node(m, iD) == 0 ) continue;

      Node(m, iZ) = Node(m, iZ) - Node(m, iR1)*( Node(m, iM) - Node(m-NLat, iM) + Node(m, iN)*R6[j] - Node(m-1, iN)*R6[j-1] );

      absH = fabs(Node(m, iZ));

      if( absH < MODELVARI.zoutZT ) Node(m, iZ) = 0.;

      if( Node(m, iZ) > Node(m, iZmax) ) Node(m, iZmax) = Node(m, iZ);

      if( MODELVARI.zoutAT && Node(m, iTime) < 0 && absH > MODELVARI.zoutAT ) Node(m, iTime) = (float)MODELVARI.time;

    }
  }

  // open bondary conditions
  if( Jmin <= 2 ) {
    for( i=2; i<=(NLon-1); i++ ) {
      m = idx(1,i);
      Node(m, iZ) = sqrt( pow(Node(m, iN),2.) + 0.25*pow((Node(m, iM)+Node(m-NLat, iM)),2.) )*C1[i];
      if( Node(m, iN) > 0 ) Node(m, iZ) = - Node(m, iZ);
    }
  }
  if( Imin <= 2 ) {
    for( j=2; j<=(NLat-1); j++ ) {
      m = idx(j,1);
      Node(m, iZ) = sqrt( pow(Node(m, iM),2.) + 0.25*pow((Node(m, iN)+Node(m-1, iN)),2.) )*C2[j];
      if( Node(m, iM) > 0 ) Node(m, iZ) = - Node(m, iZ);
    }
  }
  if( Jmax >= (NLat-1) ) {
    for( i=2; i<=(NLon-1); i++ ) {
      m = idx(NLat,i);
      Node(m, iZ) = sqrt( pow(Node(m-1, iN),2.) + 0.25*pow((Node(m, iM)+Node(m-1, iM)),2.) )*C3[i];
      if( Node(m-1, iN) < 0 ) Node(m, iZ) = - Node(m, iZ);
    }
  }
  if( Imax >= (NLon-1) ) {
    for( j=2; j<=(NLat-1); j++ ) {
      m = idx(j,NLon);
      Node(m, iZ) = sqrt( pow(Node(m-NLat, iM),2.) + 0.25*pow((Node(m, iN)+Node(m-1, iN)),2.) )*C4[j];
      if( Node(m-NLat, iM) < 0 ) Node(m, iZ) = - Node(m, iZ);
    }
  }
  if( Jmin <= 2 ) {
    m = idx(1,1);
    Node(m, iZ) = sqrt( pow(Node(m, iM),2.) + pow(Node(m, iN),2.) )*C1[1];
    if( Node(m, iN) > 0 ) Node(m, iZ) = - Node(m, iZ);
    m = idx(1,NLon);
    Node(m, iZ) = sqrt( pow(Node(m-NLat, iM),2.) + pow(Node(m, iN),2.) )*C1[NLon];
    if( Node(m, iN) > 0 ) Node(m, iZ) = - Node(m, iZ);
  }
  if( Jmin >= (NLat-1) ) {
    m = idx(NLat,1);
    Node(m, iZ) = sqrt( pow(Node(m, iM),2.) + pow(Node(m-1, iN),2.) )*C3[1];
    if( Node(m-1, iN) < 0 ) Node(m, iZ) = - Node(m, iZ);
    m = idx(NLat,NLon);
    Node(m, iZ) = sqrt( pow(Node(m-NLat, iM),2.) + pow(Node(m-1, iN),2.) )*C3[NLon];
    if( Node(m-1, iN) < 0 ) Node(m, iZ) = - Node(m, iZ);
  }

  // moment conservation
  // longitudial flux update
  #pragma omp parallel for default(shared) private(i,j,m,v1,v2)
  for( i=Imin; i<=Imax; i++ ) {
    for( j=Jmin; j<=Jmax; j++ ) {

      m = idx(j,i);

      if( (Node(m, iD)*Node(m+NLat, iD)) == 0 ) continue;

      v1 = Node(m+NLat, iZ) - Node(m, iZ);
      v2 = Node(m-1, iN) + Node(m, iN) + Node(m+NLat, iN) + Node(m+NLat-1, iN);
      Node(m, iM) = Node(m, iM) - Node(m, iR2)*v1 + Node(m, iR3)*v2;
    }
  }
  // open boundaries
  if( Jmin <= 2 ) {
    for( i=1; i<=(NLon-1); i++ ) {
      m = idx(1,i);
      Node(m, iM) = Node(m, iM) - Node(m, iR2)*(Node(m+NLat, iZ) - Node(m, iZ));
    }
  }
  if( Imin <= 2 ) {
    for( j=1; j<=NLat; j++ ) {
      m = idx(j,1);
      Node(m, iM) = Node(m, iM) - Node(m, iR2)*(Node(m+NLat, iZ) - Node(m, iZ));
    }
  }
  if( Jmax >= (NLat-1) ) {
    for( i=1; i<=(NLon-1); i++ ) {
      m = idx(NLat,i);
      Node(m, iM) = Node(m, iM) - Node(m, iR2)*(Node(m+NLat, iZ) - Node(m, iZ));
    }
  }

  // lattitudial flux update
  #pragma omp parallel for default(shared) private(i,j,m,v1,v2)
  for( i=Imin; i<=Imax; i++ ) {
    for( j=Jmin; j<=Jmax; j++ ) {

      m = idx(j,i);

      if( (Node(m, iD)*Node(m+1, iD)) == 0 ) continue;

      v1 = Node(m+1, iZ) - Node(m, iZ);
      v2 = Node(m-NLat, iM) + Node(m, iM) + Node(m-NLat+1, iM) + Node(m+1, iM);
      Node(m, iN) = Node(m, iN) - Node(m, iR4)*v1 - Node(m, iR5)*v2;
    }
  }
  // open boundaries
  if( Imin <= 2 ) {
    for( j=1; j<=(NLat-1); j++ ) {
      m = idx(j,1);
      Node(m, iN) = Node(m, iN) - Node(m, iR4)*(Node(m+1, iZ) - Node(m, iZ));
    }
  }
  if( Jmin <= 2 ) {
    for( i=1; i<=NLon; i++ ) {
      m = idx(1,i);
      Node(m, iN) = Node(m, iN) - Node(m, iR4)*(Node(m+1, iZ) - Node(m, iZ));
    }
  }
  if( Imax >= (NLon-1) ) {
    for( j=1; j<=(NLat-1); j++ ) {
      m = idx(j,NLon);
      Node(m, iN) = Node(m, iN) - Node(m, iR4)*(Node(m+1, iZ) - Node(m, iZ));
    }
  }

  // calculation area for the next step
  if( Imin > 2 ) {
    for( enlarge=0, j=Jmin; j<=Jmax; j++ ) {
      if( fabs(Node(idx(j,Imin+2), iZ)) > MODELVARI.zoutCT ) { enlarge = 1; break; }
    }
    if( enlarge ) { Imin--; if( Imin < 2 ) Imin = 2; }
  }
  if( Imax < (NLon-1) ) {
    for( enlarge=0, j=Jmin; j<=Jmax; j++ ) {
      if( fabs(Node(idx(j,Imax-2), iZ)) > MODELVARI.zoutCT ) { enlarge = 1; break; }
    }
    if( enlarge ) { Imax++; if( Imax > (NLon-1) ) Imax = NLon-1; }
  }
  if( Jmin > 2 ) {
    for( enlarge=0, i=Imin; i<=Imax; i++ ) {
      if( fabs(Node(idx(Jmin+2,i), iZ)) > MODELVARI.zoutCT ) { enlarge = 1; break; }
    }
    if( enlarge ) { Jmin--; if( Jmin < 2 ) Jmin = 2; }
  }
  if( Jmax < (NLat-1) ) {
    for( enlarge=0, i=Imin; i<=Imax; i++ ) {
      if( fabs(Node(idx(Jmax-2,i), iZ)) > MODELVARI.zoutCT ) { enlarge = 1; break; }
    }
    if( enlarge ) { Jmax++; if( Jmax > (NLat-1) ) Jmax = NLat-1; }
  }

return 0;
}
