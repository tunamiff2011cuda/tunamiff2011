# TUNAMI FF - CUDA Version 2011 May 21 (tunamiff20110521)
Numerical Codes of Tsunami Simulation (CUDA-GPU) based on IUGG/IOC Time Project, IOC. Manuals and guides 35

The Project intially started in the year 2011 and hosted in following URLS:

URL1 (Description) : https://tunamicode.wordpress.com/

URL2 (Source Code) : https://mega.nz/#F!oqhVVA4a!VPdVav4bQQVsJTuYbMIOFw

URL3 (Project): https://code.google.com/archive/p/tunami/

URL4 (FOSS): https://code.google.com/archive/p/tsunami-opensource/

URL5 (ORIGINAL): TUNAMI Modelling Manual 2006 Version : http://tunamin2.ce.metu.edu.tr/ 

Numerical simulations of Far-filed tsunamis:

Tohoku University’s Numerical Analysis Model for Investigation of Far-field Tsunamis – TUNAMI FF

Assumptions:

        The astronomical tides do not vary with respect to time throughout the tsunami simulation. The Still water Level in the computation is set equal to the water level at the beginning of the simulation
        Both temporal and spatial grid lengths vary only at the ratio of 1:3:9 and so on, if the change of them is necessary
        In the linear computation, no run up can be included, and therefore the computation is not carried out for the water depth shallower than 0.1 cm, and vertical walls are set in place of the actual slope.

Numerical simulations of far-filed tsunamis, representing transoceanic propagation requires large area of computation. Such numerical simulations of far-filed tsunamis which travels more than 1000 km over ocean should be computed in polar-coordinates by considering earth as sphere of radius R, covered by the latitude and longitude (theta, lambda). Far-filed tsunami simulations covering wide areas of computation, in turn long travel distance may yield dispersion of wave components. Therefore in order to include physical dispersion term the equations of higher order approximation are used. But long travel time yields an inevitable accumulation of numerical error, for which the computation programme should be carefully designed.

In the method of simulation, the linear long wave theory is expressed in latitude-longitude coordinates with different formulation of equations. When the liner theory is used, it is very easy to attain a high rate of vectorization in terms of programming. The current TUNAMI FF program for transoceanic propagation is composed to fully utilize the vectorizaion of parallel programming. The rate of vectorization of higher than 99% is a result of elimination of both the IF-sentences in DO-groups and the division operation.

Flow of TUNAMI FF simulation main program

    Input of Water Depth and Initial profile
    Initial condition: Still water level
    Check of the area of computation
    Equation of continuity
    Open Sea boundary condition
    Equation of Motion
    Check of the area of computation
    K>KE
    Output

Variables and constants in TUNAMI FF program

Variables:

    Water level                         Z
    Discharge flux                       M, N
    Still water depth             H
    Time history of water level PZ
    Co-ordinates of points for output of the history of water level IP, JP
    Working arrays for vector operations V1, V2, V3, V4, V5, V6 and V7

Coefficients:

    Highest water level             ZM
    Lowest water level             ZN
    Coefficients given R1, R2, R3, R4, R5, R6 AND R6=COS (THETA M+1/2)
    (THETA M+1/2)in radian     = C1
    (THETA M)in radian = C2
    (THETA M-1/2)in radian = C3
    Water depth: h = C4

Constants:

    Gravitational acceleration GG
    Circular constant pi (=1415926)
    Radius of the earth R

 

Computation is controlled by following conditions

    Size of the area for computation in longitude and latitude IG, JG
    Latitude of the southernmost end of the area for computation FL
    Area where the tsunami exists and the computation is carried out IS, JS, IE, JE
    Grid length in minute, and time step length in second             DS, DT
    Time steps of beginning and end of computation KS, KE
    Number of spatial points where the time history of water level outputNG
    Time step length in outputting the time history of water level             KC
    Time step length to output spatial wave profiles             KD

 

SUBROUTINES

    Data input of water depth             RDEPTH
    Setting of parameters required in vectorized computation PARAME
    Input of the initial condition and the initial profile                         INITIA
    Making area of computation be within the area under consideration ALIMIT
    Enlargement of the area of computation as the tsunami propagates BLIMINT
    Output and display of the spatial distribution of water level at an instant OUTPUT
    Time form the beginning of the computation, special for a NEC SX-1 CLOCK
    Computation of the equation of continuity MASS
    Open sea, boundary condition             GBOUND
    Computation of the equation of motion             MOMENT
    Check of the highest and lowest water level             MAX
    Output of the time history of water level at the point (IP, JP)             POINT
    Output of the tsunami arrival time in hour PROPA
    Output of the highest and lowest water level, and the arrival time OUTDT
    Output of the water level and the discharge flux             FILEOT

 

Initial profile computation

The vertical displacement of sea bottom is calculated with the Mansinha and Smylie method (1972), and is assumed equal to the tsunami initial profile with no modification of hydraulic effect because the horizontal size of the initial profile is sufficiently large compared to the water depth at the tsunami source.

Variables and Constants in initial profile computation

Variables:

    Grid length in meters for the Cartesian co-ordinates DX
    Grid length in degree for Spherical co-ordinates DR
    Depth of a corner of the fault plane in meter H
    Dislocation of the upper plane(u) in meter D
    Dip angle in degrees DL
    Strike angle in degree measured clockwise from North TH
    Direction of Dislocation in degree RD
    Length of fault plane in meter L
    Width of fault plane in meter W
    Co-ordinate of the origin in the area for tsunami computation Y0, X0
    Co-ordinate fo the origin of the fault plane Y0, X0

Constants

    Circular Constant A= 3.1415926
    Radius of the earth RR= 6.37E+6(m)
    and E= 1.7453E-3(m)

 

SUBROUTINES

    Computation of the initial profile             DEFORM
    Computation of the vertical displacement due to the strike slip component USCAL
    Computation of the vertical displacement due to the dip slip component UDCAL

 

Main Program Flow

The Variables and Statements used in main program summarized as below:

    Specification statement

M and N stated as REAL

In each domain, three-dimensional arrays for Z, M, N, DZ, DM and DN as well as two-dimensional arrays HZ, HM, HN, IR and IB; Note that declaration of BT(10) is crucial.

Dimensions for space are increased by one in order to include an extra row or column outside the domain under consideration. Otherwise, discharge on the boundary or IB and IR maps are sometimes not definable, according to the way of selection of I and J axes. Dimensions for time always taken to be 2, because values are changed by a subroutine CHANGE

    Input of setting values

Values of DX, DT and R(=DT/DX) are Inputs for every domain

    Setting of initial value

CALL DEPTH0 – with this command, data of water depth HZ, IB and IR maps are input in every domain. To build this subroutine of data input, the following points are taken into consideration.

    Water depth – Read water depths HZ on hydrographic charts with the z-axis positive downwards
    IB map ( a two-dimensional array) – An IB map gives the method of computation and the existence of vertical walls.

According to following rule, positive integers of one or two figures are allotted to and Input into every grid point in every domain

    The unit digit = 0, the computation is with the linear theory without the convection term
    The unit digit = 1, the computation is with the non-linear theory with the convection term included
    The tenth digit = 1, the discharge flux M in the I-direction is zero, owing to the existence of a vertical wall
    The tenth digit =2, the discharge flux N in the J-direction is zero, owing to the existence of a vertical wall
    The tenth digit = 3, the discharge fluxes M and N in the I and J directions are zero, owing to the existence of a vertical wall.
    The tenth digit = 4, no computation made

 

    IR map ( a two-dimensional array) – this map shows the existence of such structures as sea walls of finite crown height. according to the following rule, positive integers of one or two figures ae allotted to and Input into every grid point in every domain
    The unit digit, assign the address I (=1~9) of BT(I), data of the crown height of sea walls
    No tenth digit, no sea wall
    The tenth digit = 1, there is a sea wall on the computation line of discharge M in the I-direction
    The tenth digit = 2, there is a sea wall on the computation line of discharge N in the J-direction
    The tenth digit = 3, there is a sea wall both on the computation lines of discharge M and N in the I- and J- directions

CALL DEPTH – with this command, water depths HM and HN at the point where the discharge is computed are calculated. A call statement is necessary for a computation domain

CALL SETZRO – with this command, initial values are set for Z, M, N and DZ, all of which are set equal to zero. A call statement is necessary for a computation domain

    Repetition of computation with respect to time

CALL CONTIN – with this command, the water depth Z at the next time step Is computed with the equation of continuity. A call statement is necessary for a computation domain

CALL JOINTZ – with this command, the water depth is connected between domains of different time and space grid lengths. A call statement is necessary for a line of connection.

CALL MOTION – with this command, the discharges M and N at the next time step Is computed with the equation of motion. Discharges over sea walls are evaluated with the Hom-ma formula. A call statement is necessary for a computation domain

CALL BOUND – with this command, the conditions are Input at the seaside boundary. Input data of a tsunami should be prepared

CALL JOINTQ – with this command, the discharge is connected between domains of different time and space grid lengths. A call statement is necessary for a line of connection.

CALL OUTPUT – A subroutine is added at need, to output the computed results

CALL CHANGE – old data one time step before are changed with new data. For instance, newly obtained Z(I, J, 2) replaces old data and become Z(I, J, 1) which Is used in the next computation.

    Time step index – If the time step delta t varies from a domain to another, the computation procedures from CONTIN to CHANGE, those mentioned above is not always carried out at every time step except in the domain of smallest delta t. The time step at which the computation is carried out in the other domains is controlled by introducing the “time step index”.
        KT in CONTIN, MOTION, CHANGE
        KT in JOINTZ
        KT1 and KT2 In JOINTQ

 

Variables used in Main program:

Variables

Water Level                                                                           Z

Discharge In I-direction                                                        M

Discharge In J-direction                                                       N

Total water depth at point for Z                                            DZ

Total water depth at point for M                                           DM

Total water depth at point N                                     DN

Still water depth at point for Z                                              HZ

Still water depth at point for M                                             HM

Still water depth at point for N                                              HN

Map of the selection of theory                                              IB

(Linear of nonlinear and the existence of vertical walls)

Map of the existence of structures                          IR

Crown height of structures                                                   BT

Space grid length                                                                  DX

Time step length                                                                    DT

Ration DT/DX                                                                        R

Gravitational acceleration                                                    GG

2 pi                                                                                          PP

Time step of the computation                                              K

Argument to call subroutines, same value as K                 KI

The last time step                                                                  KE

Time step for output procedure                                           KOUT

Wave period                                                                          WP

Water depth                                                                           WD

Index for output                                                                      LL

Index for output                                                                      LX

 

Explanation of Subroutines

    DEPTH

Computation of the water depth at points for discharge – The still water depth at points for discharge is calculated, based upon the still water depth at point for water depth. Information of the existence of structures from the map of break waters is also input.

Variables used

Indices of HZ, HM, HN and IR in the main program                                   IG, JG

Still water depth                                                                                             HZ

No Input-Still water depth at points for discharge in I-direction    HM

No Input-Still water depth at points for discharge in J-direction   HN

Crown height of break waters                                                                      BT

Map of existence of breakwaters                                                                IR

    SETZRO

Setting of initial condition – Input of initial values of water level, water discharge and total water depth at points for water level

Variables used

Indices of Z, M, N, DZ and HZ In the main programme                 IG, JG

No Input- Initial Water level                                                   Z

No Input – Initial discharge in the I-direction                                    M

No Input – Initial discharge in the J-direction                                  N

No Input – Initial total water depth at point for discharge   DZ

Still water depth at point for water level                                          HZ

    CONTIN

Computation of the equation of continuity – Cotation of the water level and total water depth at the next time step with the equation of continuity

Variables used

Indices of Z, M, N, DZ and HZ in the main programme     IG, JG

Co-ordinates of the start of imputation                               IS, IE

(IS, JS) and of end (IE, JE)                                                  JS, JE

Water level given as Z(I,J,2)=Z(I,J,1)                                   Z

Discharge in I-direction                                                        M

Discharge in J-direction                                                       N

Total water depth at points for water level              DZ

Still water depth at points for water level                            HZ

DT/DX; ratio of time-to-space grid length                          R

Map of the selection of theory (linear or                             IB

nonlinear) and of the existence of vertical walls

Time step                                                                               KK

Time step index                                                                     KT

    JOINTZ

Connection of the water level In sapce and time – Connection of the water level between computation domains of the different delta x and delta t

Variables used

Indices of Z1 and DZ1 in the main programme                             IG1, JG1

Indices of Z2 and DZ2 in the main programme                             IG2, JG2

Water level in sender (domain of fine grids)                                  Z1

Water level in receiver (domain of coarse grids)                          Z2

Water depth at points for water level in sender                             DZ1

Water depth at points for water level in receiver                           DZ2

Co-ordinates of the start of connection                                          IS, IE

(IS, JS) and the end (IE, JE) in receiver                                         JS, JE

Co-ordinates of the start of connection in sender             ISS, JSS

Time step                                                                                           KK

Time step index for receiver                                                            KT

Space grid length in sender (DX1) and in receiver (DX2)           DX1, DX2

    MOTION

Computation of the equation of motion – Computation of the water discharge at the next time step with the equation of motion

Variables used

Indices of Z, M, N, DZ, DM, DN, HZ                                                IG, JG

HM, HN, IR, IB in the main programme

Water level                                                                                         Z

Discharge In the I-direction given as M(I,J,1)=M(I,J,2)                  M

Discharge In the J-direction given as NM(I,J,1)=N(I,J,2)  N

Total water depth at points for water level                          DZ

Total water depth at points for M                                                     DM

Total water depth at points for N                                                     DN

Still water depth at points for water level                                        HZ

Still water depth at points for M                                                       HM

Still water depth at points for N                                                        HN

Map of the existence of break waters                                             IR

(positive Integer of two figures)

Map of the selection of theory (linear or                                         IB

nonlinear) and of the existence of vertical walls

DT/DX; ratio of time-to-space grid length                                      R

Manning’s roughness in s/m1/3                                                         FM

Time step length                                                                                DT

Time step                                                                                           KK

Time step index                                                                                 KT

    JOINTQ

Connection of the water discharge in space and time – Connection of the water level between computation domains of different delta x and delta t.

Variables used

Indices of M1 and N1 in the main programme                               IG1, IG2

Indices of M2 and N2 in the main programme                               JG1, JG2

Discharge in the I-and J-directions in sender                                M1, N1

(domain of coarse grids)

Discharge in the I-and J-directions in receiver                              M2, N2

(domain of fine grids)

Co-ordinates of the start of connection                                          IS, IE, JS, JE

(IS, JS) and end (IE, JE) in receiver

Co-ordinates of start of connection in sender                                ISS, JSS

Number of extra grids at the start and end in receiver                  NDS, NDE

Computation at the connection boundary;                          INS, INE

1 for extrapolation and 2 for interpolation

Time step                                                                                           KK

Time step index (1,3,9) (KT1<KT2)                                                KT1, KT2

Space grid length in sender (DX1) and in receiver (DX2)           DX1, DX2

    CHANGE

Change of dat – Change the Index of time index from 2 to 1, at every time step of computation.

Variables used

Indices of Z, M, N and DZ In the main programme                                    IG, JG

Water level Z(I, J, 2) = Z(I, J, 1)                                                                    Z

Discharge in I-direction M(I, J, 2) = M(I, J, 1)                                             M

Discharge in I-direction N(I, J, 2) = N(I, J, 1)                                              N

Total water depth at points for water level DZ(I, J, 2)=DZ(I, J, 1) DZ

Time step                                                                                                       KK

Time step index                                                                                             KT


# Disclaimer - 
This software is experimental and is released without any warranty. The author takes no responsibility for any damage to property, health or life caused by this software. 
