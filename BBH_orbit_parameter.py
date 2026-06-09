
##############################################################################################

## This script sets parameters for rotating binary black holes
## Author: XiaoQu
## 2024/04/05
## Modified: 2025/02/10

## Can be used as input for AMSS-NCKU or the Einstein Toolkit

##############################################################################################


import AMSS_NCKU_Input as input_data
import math
import os
import sympy
import numpy      
import derivative    ## numerical differentiation


##############################################################################################

## Set each black hole's physical angular momentum from the input file

angular_momentum_BH = numpy.zeros( (input_data.puncture_number, 3) )   ## Initialize each black hole's spin angular momentum

for i in range(input_data.puncture_number):
    if ( input_data.Symmetry == "equatorial-symmetry" ):
        angular_momentum_BH[i] = [ 0.0, 0.0, (input_data.parameter_BH[i,0]**2) * input_data.parameter_BH[i,2] ]
    elif ( input_data.Symmetry == "no-symmetry" ):
        angular_momentum_BH[i] = (input_data.parameter_BH[i,0]**2) * input_data.dimensionless_spin_BH[i] 

## Set the two black hole masses
## To be consistent with literature notation, require M1 >= M2

M1 = input_data.parameter_BH[0,0]
M2 = input_data.parameter_BH[1,0] 

## Set the dimensionless spins of the black holes

S1 = angular_momentum_BH[0] / M1**2
S2 = angular_momentum_BH[1] / M2**2

## Set the orbital semi-major axis and eccentricity in the center-of-mass frame

D0 = input_data.Distance
e0 = input_data.e0

##############################################################################################


##############################################################################################

## Generate orbital parameters for a quasi-circular rotating binary black hole system
## Author: XiaoQu
## Uses post-Newtonian expansions to obtain a quasi-circular orbit
## Updated up to 3rd post-Newtonian (3PN) order

## Arguments
## Masses M1 and M2 (ensure M1 >= M2)
## Spins S1 and S2 as 3-element numpy vectors (for numpy.dot)
## Initial separation D0
## Orbital eccentricity e0 (currently only circular PN formulas implemented; eccentric cases may be added later)

def generate_BBH_orbit_parameters( M1, M2, S1, S2, D0, e0 ):

    print()
    print()
    print("Compute binary orbital characteristic quantities using the effective single-body model")
    print()

    print()
    print(f"Input binary masses:       M1 = {M1}  M2 = {M2}")
    print(f"Input binary dimensionless spins: S1 = {S1}  S2 = {S2}")
    print(f"Input orbital semi-major axis and eccentricity: a0 = D0/2 = {D0/2.0}  e0 = {e0}")
    print()
    print("Begin calculations")
    print()

    ##################################################

    ## Compute mass ratio, reduced mass, and dimensionless masses
    M_total = M1 + M2
    print( f"Binary masses: M1 = {M1}  M2 = {M2}  Newtonian total mass: M_total = {M_total} " )

    ## Compute reduced mass
    M_mu = M1 * M2 / M_total
    print( "Reduced mass for effective single-body: M_mu =", M_mu )

    ## Set mass ratio
    ## Note: For consistency with TwoPuncture, require M1 >= M2
    m_q   = M1 / M2
    m_eta = m_q / ( (1.0 + m_q)**2 )
    print( "Dimensionless mass ratio Q = M1 / M2 =", m_q )

    ## Set dimensionless masses
    m1   = M1 / M_total
    m2   = M2 / M_total
    m_mu = M_mu / M_total
    print( f"Dimensionless masses: m1 = {m1}  m2 = {m2}  m_mu = {m_mu} " )
    print( "Dimensionless reduced mass m_eta = Q / (1+Q)^2 =", m_eta )

    ##################################################    
    
    ## From center-of-mass semi-major axis and eccentricity, compute orbital parameters at t=0

    ## Under classical mechanics, the binary orbit is equivalent to a reduced mass moving in a central potential
    ## Relation between radius r, phase phi, semi-major axis a and eccentricity e:
    ## r = a*(1-e^2)/(1+e*cos(phi))
    ## Phase phi = 0 or phi = pi corresponds to periapsis/apapsis (semi-major/semi-minor axis)

    ## Compute semi-major and semi-minor axes
    a0 = D0 / 2.0
    a_long  = a0 * (1.0 + e0)
    a_short = a0 * (1.0 - e0)

    ## Compute individual semi-major axes from center-of-mass definition
    ## M1/M2 = a2/a1
    R10 = D0 * M2 / M_total
    R20 = D0 * M1 / M_total

    print(                        )
    print( " Compute orbital coordinate parameters " )
    print(                        )
    print( " Choose coordinates so that at t=0 the +y axis aligns with the semi-major axis; the binary orbits about +z " )
    print( " At t=0: y corresponds to radial coordinate R, x to angular coordinate phi, and z = 0 " )
    print( " Similarly at t=0: Py corresponds to radial momentum Pr, Px to tangential momentum P_phi, and Pz = 0 " )
    ## print( " Note: this does not imply z remains zero during long evolutions; misaligned spins can tilt the orbit out of the xy plane " ) 
    print()

    ## Initialize position coordinates
    position1 = [0.0, 0.0, 0.0]
    position2 = [0.0, 0.0, 0.0]

    ## Initialize momentum coordinates
    momentum1 = [0.0, 0.0, 0.0]
    momentum2 = [0.0, 0.0, 0.0]


    position1[1] =   R10    
    position2[1] = - R20

    print( "Binary coordinates at t=0:" )
    print( f"Y1 = {position1}  Y2 = {position2}" )
    print(                                        )


    ########################################

    ## Compute orbital angular frequency

    R = D0
    epsilon = (m1 + m2) / R 
    print( " Post-Newtonian expansion parameter epsilon = M/r = ", epsilon )
    print(                                             )
    ## Gravitational-wave scaling: use dimensionless total mass m = m1 + m2 = 1.
    ## The initial separation D0 is assumed in total-mass units (physical separation / M_total).

    ## 3PN post-Newtonian results
    ## Note: the following formulas assume circular orbits; eccentricity is not included
    ## Based on:
    ## James Healy, Carlos O. Lousto, Hiroyuki Nakano, and Yosef Zlochower
    ## "Post-Newtonian Quasicircular Initial Orbits for Numerical Relativity"
    ## arXiv:1702.00872 [gr-qc]
    ## Note: in arXiv:1702.00872 the mass ratio q<1, so their S1 corresponds to our S2

    ## Set orbital angular frequency (can be adjusted later)

    Omega_0PN = epsilon**1.5

    Omega_correction_1PN  = - 0.5 * epsilon * ( ( 3.0*(m_q**2) + 5.0*m_q + 3.0 ) / (1.0+m_q)**2 ) 
    Omega_correction_15PN = - 0.25 * ( epsilon**(1.5) )                               \
                              * (   (3.0 + 4.0*m_q) * m_q * S2[2] / ( (1.0+m_q)**2 )  \
                                  + (3.0*m_q + 4.0)       * S1[2] / ( (1.0+m_q)**2 )        \
                                )
    Omega_correction_2PN  = epsilon**2 \
                            * (   (1.0/16.0) * ( 24.0*(m_q**4) + 103.0*(m_q**3) + 164.0*(m_q**2) + 103.0*m_q + 24.0 ) \
                                   / ( (1.0+m_q)**4 )                                 \
                                - 1.5  * (m_q**2) * (S2[0]**2) / ( (1.0+m_q)**2 )     \
                                + 0.75 * (m_q**2) * (S2[1]**2) / ( (1.0+m_q)**2 )     \
                                + 0.75 * (m_q**2) * (S2[2]**2) / ( (1.0+m_q)**2 )     \
                                - 3.0  * m_q * S1[0] * S2[0] / ( (1.0+m_q)**2 )       \
                                + 1.5  * m_q * S1[1] * S2[1] / ( (1.0+m_q)**2 )       \
                                + 1.5  * m_q * S1[2] * S2[2] / ( (1.0+m_q)**2 )       \
                                - 1.5  * (S1[0]**2) / ( (1.0+m_q)**2 )                \
                                + 0.75 * (S1[1]**2) / ( (1.0+m_q)**2 )                \
                                + 0.75 * (S1[2]**2) / ( (1.0+m_q)**2 )
                              )
    Omega_correction_25PN = (3.0/16.0) * (epsilon**(2.5))   \
                            * (   S2[2] * m_q * ( 16.0*(m_q**3) + 30.0*(m_q**2) + 34.0*m_q + 13.0 ) / ( (1.0+m_q)**4 )  \
                                + S1[2]       * ( 13.0*(m_q**3) + 34.0*(m_q**2) + 30.0*m_q + 16.0 ) / ( (1.0+m_q)**2 )  \
                              )
    Omega_correction_3PN  = epsilon**3 \
                            * (   (167.0/128.0) * (math.pi**2) * m_q / ( (1.0+m_q)**2 )           \
                                - (   120.0*(m_q**6) + 2744.0*(m_q**5) + 10049.0*(m_q**4)         \
                                    + 14820.0*(m_q**3) + 10049.0*(m_q**2) + 2744.0*m_q + 120.0    \
                                   ) / ( 96.0 * ((1.0+m_q)**6) )                                  \
                                + (1.0/16.0) * (m_q**2) * (S2[0]**2) * ( 76.0*(m_q**2) + 180.0*m_q + 155.0 )   / ( (1.0+m_q)**4 ) \
                                - (1.0/8.0)  * (m_q**2) * (S2[1]**2) * ( 43.0*(m_q**2) + 85.0*m_q  + 55.0 )    / ( (1.0+m_q)**4 ) \
                                - (1.0/32.0) * (m_q**2) * (S2[2]**2) * ( 2.0*m_q + 5.0 ) * ( 14.0*m_q + 27.0 ) / ( (1.0+m_q)**4 ) \
                                + (1.0/16.0)            * (S1[0]**2) * ( 155.0*(m_q**2) + 180.0*m_q + 76.0 )   / ( (1.0+m_q)**4 ) \
                                - (1.0/8.0)             * (S1[1]**2) * ( 55.0*(m_q**2)  + 85.0*m_q  + 43.0 )   / ( (1.0+m_q)**4 ) \
                                - (1.0/32.0)            * (S1[2]**2) * ( 27.0*m_q + 14.0 ) * ( 5.0*m_q + 2.0 ) / ( (1.0+m_q)**4 ) \
                                + (1.0/8.0)  * m_q * S1[0] * S2[0] * ( 120.0*(m_q**2) + 187.0*m_q + 120.0 ) / ( (1.0+m_q)**4 )    \
                                - 0.25       * m_q * S1[1] * S2[1] * ( 54.0*(m_q**2) + 95.0*m_q   + 54.0 )  / ( (1.0+m_q)**4 )    \
                                - (1.0/16.0) * m_q * S1[2] * S2[2] * ( 96.0*(m_q**2) + 127.0*m_q  + 96.0 )  / ( (1.0+m_q)**4 )    \
                              )

    Omega_1PN  = Omega_0PN * ( 1.0 + Omega_correction_1PN )
    Omega_15PN = Omega_0PN * ( 1.0 + Omega_correction_1PN + Omega_correction_15PN )
    Omega_2PN  = Omega_0PN * ( 1.0 + Omega_correction_1PN + Omega_correction_15PN    \
                                   + Omega_correction_2PN                            \
                             )
    Omega_25PN = Omega_0PN * ( 1.0 + Omega_correction_1PN + Omega_correction_15PN    \
                                   + Omega_correction_2PN + Omega_correction_25PN    \
                             )
    Omega_3PN  = Omega_0PN * ( 1.0 + Omega_correction_1PN + Omega_correction_15PN    \
                                   + Omega_correction_2PN + Omega_correction_25PN    \
                                   + Omega_correction_3PN                            \
                             )

    print()
    print( "Omega (0PN) =", Omega_0PN )
    print( "Omega (1PN) =", Omega_1PN )
    print( "Omega (1.5PN) =", Omega_15PN )
    print( "Omega (2PN) =", Omega_2PN )
    print( "Omega (2.5PN) =", Omega_25PN )
    print( "Omega (3PN) =", Omega_3PN )
    print()

    ########################################

    ## Set orbital angular momentum (can be adjusted later)

    Pt_0PN = (epsilon**0.5) * m_q / ( (1+m_q)**2.0 ) 

    Pt_correction_1PN  = 2.0 * epsilon
    Pt_correction_15PN = epsilon**1.5 \
                         * ( - 0.75 * (3.0 + 4.0*m_q) * m_q * S2[2] / ( (1.0+m_q)**2 ) \
                             - 0.75 * (3.0*m_q + 4.0)       * S1[2] / ( (1.0+m_q)**2 )       \
                           )
    Pt_correction_2PN  = epsilon**2.0 \
                         * (  (1.0/16.0) * ( 42.0*(m_q**2) + 41.0*m_q + 42.0 ) / ( (1.0+m_q)**2 )   \
                             - 1.5  * (m_q**2) * (S2[0]**2) / ( (1.0+m_q)**2 )                      \
                             + 0.75 * (m_q**2) * (S2[1]**2) / ( (1.0+m_q)**2 )                      \
                             + 0.75 * (m_q**2) * (S2[2]**2) / ( (1.0+m_q)**2 )                      \
                             - 3.0  * m_q * S1[0] * S2[0] / ( (1.0+m_q)**2 )                        \
                             + 1.5  * m_q * S1[1] * S2[1] / ( (1.0+m_q)**2 )                        \
                             + 1.5  * m_q * S1[2] * S2[2] / ( (1.0+m_q)**2 )                        \
                             - 1.5  * (S1[0]**2) / ( (1.0+m_q)**2 )                                 \
                             + 0.75 * (S1[1]**2) / ( (1.0+m_q)**2 )                                 \
                             + 0.75 * (S1[2]**2) / ( (1.0+m_q)**2 )                                 \
                           )
    Pt_correction_25PN = epsilon**2.5 \
                         * ( - (1.0/16.0) * ( 72.0*(m_q**3) + 116.0*(m_q**2) + 60.0*m_q + 13.0 ) \
                                          * m_q * S2[2] / ( (1.0+m_q)**4 )                       \
                             - (1.0/16.0) * ( 13.0*(m_q**3) + 60.0*(m_q**2) + 116.0*m_q + 72.0 ) \
                                          * S1[2] / ( (1.0+m_q)**4 )                             \
                           )
    Pt_correction_3PN  = epsilon**3.0 \
                         * (   (163.0/128.0) * (math.pi**2) * m_q / ( (1.0+m_q)**2 )                                   \
                             + (1.0/32.0) * ( 120.0*(m_q**4) - 659.0*(m_q**3) - 1532.0*(m_q**2) - 659.0*m_q + 120.0 )  \
                                / ( (1.0+m_q)**4 )                                                                     \
                             - (1.0/16.0) * (S2[0]**2) * (m_q**2) * ( 80.0*(m_q**2)             - 59.0 ) / ( (1.0+m_q)**4 ) \
                             - 0.5        * (S2[1]**2) * (m_q**2) * (        m_q**2  + 10.0*m_q + 8.0  ) / ( (1.0+m_q)**4 ) \
                             + (1.0/32.0) * (S2[2]**2) * (m_q**2) * ( 128.0*(m_q**2) + 56.0*m_q - 27.0 ) / ( (1.0+m_q)**4 ) \
                             - (1.0/16.0) * (S1[0]**2)            * ( 80.0 - 59.0*(m_q**2) )             / ( (1.0+m_q)**4 ) \
                             - 0.5        * (S1[1]**2)            * ( 8.0*(m_q**2) + 10.0*m_q + 1.0 )    / ( (1.0+m_q)**4 ) \
                             + (1.0/32.0) * (S1[2]**2)            * ( 128.0 + 56.0*m_q - 27.0*(m_q**2) ) / ( (1.0+m_q)**4 ) \
                             + (1.0/8.0)  * S1[0] * S2[0] * m_q   * ( 12.0*(m_q**2) + 35.0*m_q + 12.0 )  / ( (1.0+m_q)**4 ) \
                             - 0.25       * S1[1] * S2[1] * m_q   * ( 27.0*(m_q**2) + 58.0*m_q + 27.0 )  / ( (1.0+m_q)**4 ) \
                             + (1.0/32.0) * S1[2] * S2[2] * m_q   * ( 60.0*(m_q**2) + 13.0*m_q + 60.0 )  / ( (1.0+m_q)**4 ) \
                           )

    Pt_1PN  = Pt_0PN * ( 1.0 + Pt_correction_1PN )
    Pt_15PN = Pt_0PN * ( 1.0 + Pt_correction_1PN + Pt_correction_15PN )
    Pt_2PN  = Pt_0PN * ( 1.0 + Pt_correction_1PN + Pt_correction_15PN    \
                             + Pt_correction_2PN                         \
                       )
    Pt_25PN = Pt_0PN * ( 1.0 + Pt_correction_1PN + Pt_correction_15PN    \
                             + Pt_correction_2PN + Pt_correction_25PN    \
                       )
    Pt_3PN  = Pt_0PN * ( 1.0 + Pt_correction_1PN + Pt_correction_15PN    \
                             + Pt_correction_2PN + Pt_correction_25PN    \
                             + Pt_correction_3PN
                       )

    print()
    print( "Pt (0PN) =", Pt_0PN )
    print( "Pt (1PN) =", Pt_1PN )
    print( "Pt (1.5PN) =", Pt_15PN )
    print( "Pt (2PN) =", Pt_2PN )
    print( "Pt (2.5PN) =", Pt_25PN )
    print( "Pt (3PN) =", Pt_3PN )
    print()

    ########################################

    ## Compute the ADM mass of the binary system
    ## Based on:
    ## Antoni Ramos-Buades, Sascha Husa, and Geraint Pratten
    ## "Simple procedures to reduce eccentricity of binary black hole simulations"
    ## arXiv:1810.00036 [gr-qc]

    ############################

    ## Define ADM mass function
    ## Expansion in M/R

    def M_ADM(r):
        
        mass     =  m1 + m2
        epsilon0 = (m1 + m2) / r
        
        adm_correction_0PN  = - 0.5 * epsilon0 * m_q / ( (1+m_q)**2 ) 
        adm_correction_1PN  =  (1.0/8.0) * ( epsilon0**2 )                                   \
                                * m_q * ( 7.0*(m_q**2) + 13.0*m_q + 7.0 ) / ( (1.0+m_q)**4 ) 
        adm_correction_15PN = - 0.25 * ( epsilon0**2.5 )                                     \
                                * (   m_q**2 * ( 3.0 + 4.0*m_q ) * S2[2] / ( (1.0+m_q)**4 )  \
                                    + m_q    * ( 3.0*m_q + 4.0 ) * S1[2] / ( (1.0+m_q)**4 )  \
                                  )
        adm_correction_2PN  = epsilon0**3  \
                              * (   (1.0/16.0) * m_q * ( 9.0*m_q**4 + 16.0*(m_q**3) + 13.0*(m_q**2) + 16.0*m_q + 9.0 ) \
                                     / ( (1.0+m_q)**6 )                                                                \
                                  - 0.5  * (S2[0]**2) * (m_q**3) / ( (1.0+m_q)**4 )      \
                                  + 0.25 * (S2[1]**2) * (m_q**3) / ( (1.0+m_q)**4 )      \
                                  + 0.25 * (S2[2]**2) * (m_q**3) / ( (1.0+m_q)**4 )      \
                                  - 1.0  * S1[0] * S2[0] * (m_q**2) / ( (1.0+m_q)**4 )   \
                                  + 0.5  * S1[1] * S2[1] * (m_q**2) / ( (1.0+m_q)**4 )   \
                                  + 0.5  * S1[2] * S2[2] * (m_q**2) / ( (1.0+m_q)**4 )   \
                                  - 0.5  * (S1[0]**2) * m_q / ( (1.0+m_q)**4 )           \
                                  + 0.25 * (S1[1]**2) * m_q / ( (1.0+m_q)**4 )           \
                                  + 0.25 * (S1[2]**2) * m_q / ( (1.0+m_q)**4 )           \
                                )
        adm_correction_25PN = - (1.0/16.0) * epsilon0**3.5  \
                                * (   S2[2] * (m_q**2) * ( 32.0*(m_q**3) + 42.0*(m_q**2) + 14.0*m_q +1.0 ) / ( (1.0+m_q)**6 )  \
                                    + S1[2] * m_q      * ( m_q**3 + 14.0*(m_q**2) + 42.0*m_q + 32.0 )      / ( (1.0+m_q)**6 )  \
                                  )
        adm_correction_3PN  = epsilon0**4  \
                              * (   (81.0/128.0) * (math.pi**2) * (m_q**2) / ( (1.0+m_q)**4 )    \
                                  + (     537.0*(m_q**6) - 3497.0*(m_q**5) - 18707.0*(m_q**4)    \
                                      - 29361.0*(m_q**3) - 18707.0*(m_q**2) - 3497.0*m_q + 537.0 \
                                    ) * (m_q/384.0) / ( (1.0+m_q)**8 )                           \
                                  - (1.0/16.0) * (S2[0]**2)    * (m_q**3) * ( 52.0*(m_q**2) + 12.0*m_q - 25.0 ) / ( (1.0+m_q)**6 )  \
                                  + (1.0/8.0)  * (S2[1]**2)    * (m_q**3) * (       m_q**2  - 17.0*m_q - 15.0 ) / ( (1.0+m_q)**6 )  \
                                  + (1.0/16.0) * (S2[2]**2)    * (m_q**3) * ( 50.0*(m_q**2) + 38.0*m_q + 3.0  ) / ( (1.0+m_q)**6 )  \
                                  + (1.0/16.0) * (S1[0]**2)    * m_q      * ( 25.0*(m_q**2) - 12.0*m_q - 52.0 ) / ( (1.0+m_q)**6 )  \
                                  - (1.0/8.0)  * (S1[1]**2)    * m_q      * ( 15.0*(m_q**2) + 17.0*m_q -  1.0 ) / ( (1.0+m_q)**6 )  \
                                  + (1.0/16.0) * (S1[2]**2)    * m_q      * (  3.0*(m_q**2) + 38.0*m_q + 50.0 ) / ( (1.0+m_q)**6 )  \
                                  + (9.0/8.0)  * S1[0] * S2[0] * (m_q**3)                                       / ( (1.0+m_q)**6 )  \
                                  - (3.0/4.0)  * S1[1] * S2[1] * (m_q**2) * (  4.0*(m_q**2) +  9.0*m_q + 4.0)   / ( (1.0+m_q)**6 )  \
                                  + (3.0/8.0)  * S1[2] * S2[2] * (m_q**2) * ( 10.0*(m_q**2) + 21.0*m_q + 10.0 ) / ( (1.0+m_q)**6 )  \
                                )

        ADM_0PN  = mass * ( 1.0 + adm_correction_0PN )
        ADM_1PN  = mass * ( 1.0 + adm_correction_0PN  + adm_correction_1PN )
        ADM_15PN = mass * ( 1.0 + adm_correction_0PN  + adm_correction_1PN      \
                                + adm_correction_15PN 
                          )
        ADM_2PN  = mass * ( 1.0 + adm_correction_0PN  + adm_correction_1PN      \
                                + adm_correction_15PN + adm_correction_2PN      \
                          )
        ADM_25PN = mass * ( 1.0 + adm_correction_0PN  + adm_correction_1PN      \
                                + adm_correction_15PN + adm_correction_2PN      \
                                + adm_correction_25PN                           \
                          )
        ADM_3PN  = mass * ( 1.0 + adm_correction_0PN  + adm_correction_1PN      \
                                + adm_correction_15PN + adm_correction_2PN      \
                                + adm_correction_25PN + adm_correction_3PN      \
                          )

        return ADM_0PN, ADM_1PN, ADM_15PN, ADM_2PN, ADM_25PN, ADM_3PN
    
    ############################

    ## Define an alternative ADM mass function
    ## Expansion in orbital frequency Omega
    ## Based on:
    ## James Healy, Carlos O. Lousto, Hiroyuki Nakano, and Yosef Zlochower
    ## "Post-Newtonian Quasicircular Initial Orbits for Numerical Relativity"
    ## arXiv:1702.00872 [gr-qc]
    ## Note: in that paper q<1, so their S1 corresponds to our S2

    def M_ADM_another(Omega):
        
        mass     =  m1 + m2
        epsilon0 = mass * Omega
        
        adm_correction_0PN  = 0.0
        adm_correction_1PN  = - 0.5 * ( epsilon0**(2.0/3.0) )
        adm_correction_15PN = ( epsilon0**(4.0/3.0) )   \
                              * (1.0/24.0) * ( 9.0*(m_q**2) + 19.0*m_q + 9.0 ) / ( (1.0+m_q)**2 )
        adm_correction_2PN  = - (1.0/3.0) * ( epsilon0**(5.0/3.0) )   \
                                          * (   S2[2] * m_q * ( 4.0*m_q + 3.0 ) / ( (1.0+m_q)**2 ) \
                                              + S1[2] *       ( 3.0*m_q + 4.0 ) / ( (1.0+m_q)**2 ) \
                                            )                                                      \
                              + ( epsilon0**2 )   \
                                * (   (1.0/48.0) * (   81.0*(m_q**4) + 267.0*(m_q**3)     \
                                                     + 373.0*(m_q**2) + 267.0*m_q + 81.0  \
                                                   ) / ( (1.0+m_q)**4 )                   \
                                    -       (S2[0]**2) * (m_q**2) / ( (1.0+m_q)**2 )      \
                                    + 0.5 * (S2[1]**2) * (m_q**2) / ( (1.0+m_q)**2 )      \
                                    + 0.5 * (S2[2]**2) * (m_q**2) / ( (1.0+m_q)**2 )      \
                                    -       (S1[0]**2) / ( (1.0+m_q)**2 )                 \
                                    + 0.5 * (S1[1]**2) / ( (1.0+m_q)**2 )                 \
                                    + 0.5 * (S1[2]**2) / ( (1.0+m_q)**2 )                 \
                                    - 2.0 * S1[0] * S2[0] * m_q / ( (1.0+m_q)**2 )        \
                                    +       S1[1] * S2[1] * m_q / ( (1.0+m_q)**2 )        \
                                    +       S1[2] * S2[2] * m_q / ( (1.0+m_q)**2 )        \
                                  )
        adm_correction_25PN = - (1.0/18.0) * ( epsilon0**(7.0/3.0) )   \
                                * (   S2[2] * m_q * ( 72.0*(m_q**3) + 140.0*(m_q**2) + 96.0*m_q + 27.0 ) / ( (1.0+m_q)**4 ) \
                                    + S1[2]       * ( 27.0*(m_q**3) + 96.0*(m_q**2) + 140.0*m_q + 72.0 ) / ( (1.0+m_q)**4 ) \
                                  )
        adm_correction_3PN  = ( epsilon0**(8.0/3.0) )   \
                              * (   (205.0/192.0) * (math.pi**2) * m_q / ( (1.0+m_q)**2 )              \
                                  + (    54675.0*(m_q**6) + 18045.0*(m_q**5) - 411525.0*(m_q**4)       \
                                      - 749755.0*(m_q**3) - 411525.0*(m_q**2) + 18045.0*m_q + 54675.0  \
                                    ) / ( 10368.0 * ( (1.0+m_q)**6 ) )                                 \
                                  - (5.0/24.0) * (S2[0]**2) * (m_q**2) * ( 20.0*(m_q**2) + 4.0*m_q - 11.0 ) / ( (1.0+m_q)**4 )  \
                                  - (5.0/12.0) * (S2[1]**2) * (m_q**2) * (       m_q**2  + 9.0*m_q + 7.0  ) / ( (1.0+m_q)**4 )  \
                                  + (5.0/36.0) * (S2[2]**2) * (m_q**2) * ( 13.0*(m_q**2) - 3.0*m_q - 9.0  ) / ( (1.0+m_q)**4 )  \
                                  + (5.0/4.0)  * S1[0] * S2[0] * (m_q**2) / ( (1.0+m_q)**4 )                                    \
                                  - (6.0/5.0)  * S1[1] * S2[1] * m_q * ( 2.0*m_q + 3.0 ) * ( 3.0*m_q + 2.0 ) / ( (1.0+m_q)**4 ) \
                                  + (5.0/18.0) * S1[2] * S2[2] * m_q * ( 3.0*(m_q**2) + 7.0*m_q + 3.0 ) / ( (1.0+m_q)**4 )      \
                                  + (1.0/24.0) * (S1[0]**2) * ( 55.0*(m_q**2) - 20.0*m_q - 100.0 ) / ( (1.0+m_q)**4 )           \
                                  - (1.0/12.0) * (S1[1]**2) * ( 35.0*(m_q**2) + 45.0*m_q + 5.0   ) / ( (1.0+m_q)**4 )           \
                                  - (1.0/36.0) * (S1[2]**2) * ( 45.0*(m_q**2) + 15.0*m_q - 65.0  ) / ( (1.0+m_q)**4 )           \
                                )
        
        ADM_0PN  = mass * ( 1.0 + ( m_q / ( (1+m_q)**2 ) ) * adm_correction_0PN )
        ADM_1PN  = mass * ( 1.0 + ( m_q / ( (1+m_q)**2 ) ) * (   adm_correction_0PN  \
                                                               + adm_correction_1PN  \
                                                              ) 
                          )
        ADM_15PN = mass * ( 1.0 + ( m_q / ( (1+m_q)**2 ) ) * (   adm_correction_0PN   \
                                                               + adm_correction_1PN   \
                                                               + adm_correction_15PN  \
                                                              ) 
                          )
        ADM_2PN  = mass * ( 1.0 + ( m_q / ( (1+m_q)**2 ) ) * (   adm_correction_0PN  \
                                                               + adm_correction_1PN  \
                                                               + adm_correction_15PN \
                                                               + adm_correction_2PN  \
                                                              )
                          )
        ADM_25PN = mass * ( 1.0 + ( m_q / ( (1+m_q)**2 ) ) * (   adm_correction_0PN  \
                                                               + adm_correction_1PN  \
                                                               + adm_correction_15PN \
                                                               + adm_correction_2PN  \
                                                               + adm_correction_25PN \
                                                              )
                          )
        ADM_3PN  = mass * ( 1.0 + ( m_q / ( (1+m_q)**2 ) ) * (   adm_correction_0PN  \
                                                               + adm_correction_1PN  \
                                                               + adm_correction_15PN \
                                                               + adm_correction_2PN  \
                                                               + adm_correction_25PN \
                                                               + adm_correction_3PN  \
                                                              )
                          )

        return ADM_0PN, ADM_1PN, ADM_15PN, ADM_2PN, ADM_25PN, ADM_3PN
    
    ############################
    
    ## Define derivative of the ADM mass with respect to r

    def dADM_dr(r):
        
        mass     =  m1 + m2
        epsilon0 = (m1 + m2) / r
        
        dADM_correction_0PN  = ( - (epsilon0**2) / mass ) * ( - 0.5 * m_q / ((1+m_q)**2) ) 
        dADM_correction_1PN  = ( - 2.0 * (epsilon0**3) / mass )                                           \
                                * (1.0/8.0) * m_q * ( 7.0*(m_q**2) + 13.0*m_q + 7.0 ) / ( (1.0+m_q)**4 ) 
        dADM_correction_15PN = ( - 2.5 * (epsilon0**3.5) / mass )                            \
                                * ( - 0.25 )                                                 \
                                * (   m_q**2 * ( 3.0 + 4.0*m_q ) * S2[2] / ( (1.0+m_q)**4 )  \
                                    + m_q    * ( 3.0*m_q + 4.0 ) * S1[2] / ( (1.0+m_q)**4 )  \
                                  )
        dADM_correction_2PN  = ( - 3.0 * (epsilon0**4) / mass )   \
                                * (   (1.0/16.0) * m_q * ( 9.0*(m_q**4) + 16.0*(m_q**3) + 13.0*(m_q**2) + 16.0*m_q + 9.0 ) \
                                      / ( (1.0+m_q)**6 )                                                                   \
                                    - 0.5  * (S2[0]**2) * (m_q**3) / ( (1.0+m_q)**4 )      \
                                    + 0.25 * (S2[1]**2) * (m_q**3) / ( (1.0+m_q)**4 )      \
                                    + 0.25 * (S2[2]**2) * (m_q**3) / ( (1.0+m_q)**4 )      \
                                    - 1.0  *  S1[0] * S2[0] * (m_q**2) / ( (1.0+m_q)**4 )  \
                                    + 0.5  *  S1[1] * S2[1] * (m_q**2) / ( (1.0+m_q)**4 )  \
                                    + 0.5  *  S1[2] * S2[2] * (m_q**2) / ( (1.0+m_q)**4 )  \
                                    - 0.5  * (S1[0]**2) * m_q / ( (1.0+m_q)**4 )           \
                                    + 0.25 * (S1[1]**2) * m_q / ( (1.0+m_q)**4 )           \
                                    + 0.25 * (S1[2]**2) * m_q / ( (1.0+m_q)**4 )           \
                                  )
        dADM_correction_25PN = ( - 3.5 * (epsilon0**4.5) / mass )  \
                                * ( - 1.0/16.0 )                   \
                                * (   S2[2] * (m_q**2) * ( 32.0*(m_q**3) + 42.0*(m_q**2) + 14.0*m_q + 1.0  ) / ( (1.0+m_q)**6 )  \
                                    + S1[2] * m_q      * (       m_q**3  + 14.0*(m_q**2) + 42.0*m_q + 32.0 ) / ( (1.0+m_q)**6 )  \
                                  )
        dADM_correction_3PN  = ( - 4.0 * (epsilon0**5) / mass )   \
                                * (   (81.0/128.0) * (math.pi**2) * (m_q**2) / ( (1.0+m_q)**4 )    \
                                    + (   537.0*(m_q**6) - 3497.0*(m_q**5) - 18707.0*(m_q**4)      \
                                        - 29361.0*(m_q**3) - 18707.0*(m_q**2) - 3497.0*m_q + 537.0 \
                                      ) * (m_q/384.0) / ( (1.0+m_q)**8 )                           \
                                    - (1.0/16.0) * (S2[0]**2)    * (m_q**3) * ( 52.0*(m_q**2) + 12.0*m_q - 25.0 ) / ( (1.0+m_q)**6 )  \
                                    + (1.0/8.0)  * (S2[1]**2)    * (m_q**3) * (       m_q**2  - 17.0*m_q - 15.0 ) / ( (1.0+m_q)**6 )  \
                                    + (1.0/16.0) * (S2[2]**2)    * (m_q**3) * ( 50.0*(m_q**2) + 38.0*m_q + 3.0  ) / ( (1.0+m_q)**6 )  \
                                    + (1.0/16.0) * (S1[0]**2)    * m_q      * ( 25.0*(m_q**2) - 12.0*m_q - 52.0 ) / ( (1.0+m_q)**6 )  \
                                    - (1.0/8.0)  * (S1[1]**2)    * m_q      * ( 15.0*(m_q**2) + 17.0*m_q -  1.0 ) / ( (1.0+m_q)**6 )  \
                                    + (1.0/16.0) * (S1[2]**2)    * m_q      * (  3.0*(m_q**2) + 38.0*m_q + 50.0 ) / ( (1.0+m_q)**6 )  \
                                    + (9.0/8.0)  * S1[0] * S2[0] * (m_q**3)                                       / ( (1.0+m_q)**6 )  \
                                    - (3.0/4.0)  * S1[1] * S2[1] * (m_q**2) * (  4.0*(m_q**2) +  9.0*m_q + 4.0)   / ( (1.0+m_q)**6 )  \
                                    + (3.0/8.0)  * S1[2] * S2[2] * (m_q**2) * ( 10.0*(m_q**2) + 21.0*m_q + 10.0 ) / ( (1.0+m_q)**6 )  \
                                  )

        dADM_dr_0PN  = mass * (   dADM_correction_0PN )
        dADM_dr_1PN  = mass * (   dADM_correction_0PN  + dADM_correction_1PN )
        dADM_dr_15PN = mass * (   dADM_correction_0PN  + dADM_correction_1PN      \
                                + dADM_correction_15PN 
                              )
        dADM_dr_2PN  = mass * (   dADM_correction_0PN  + dADM_correction_1PN      \
                                + dADM_correction_15PN + dADM_correction_2PN      \
                              )
        dADM_dr_25PN = mass * (   dADM_correction_0PN  + dADM_correction_1PN      \
                                + dADM_correction_15PN + dADM_correction_2PN      \
                                + dADM_correction_25PN                           \
                              )
        dADM_dr_3PN  = mass * (   dADM_correction_0PN  + dADM_correction_1PN      \
                                + dADM_correction_15PN + dADM_correction_2PN      \
                                + dADM_correction_25PN + dADM_correction_3PN      \
                              )
        
        return dADM_dr_0PN, dADM_dr_1PN, dADM_dr_15PN, dADM_dr_2PN, dADM_dr_25PN, dADM_dr_3PN
    
    ############################
    
    ADM_Mass_0PN,  ADM_Mass_1PN, ADM_Mass_15PN, ADM_Mass_2PN, ADM_Mass_25PN, ADM_Mass_3PN  = M_ADM(R)

    print()
    print( "ADM Mass (0PN) =", ADM_Mass_0PN )
    print( "ADM Mass (1PN) =", ADM_Mass_1PN )
    print( "ADM Mass (1.5PN) =", ADM_Mass_15PN )
    print( "ADM Mass (2PN) =", ADM_Mass_2PN )
    print( "ADM Mass (2.5PN) =", ADM_Mass_25PN )
    print( "ADM Mass (3PN) =", ADM_Mass_3PN )
    print()

    Omega = Omega_3PN
    ADM_Mass_another_0PN,  ADM_Mass_another_1PN,  ADM_Mass_another_15PN, \
    ADM_Mass_another_2PN,  ADM_Mass_another_25PN, ADM_Mass_another_3PN   \
    = M_ADM_another(Omega)

    print()
    print( "ADM Mass (Omega expansion) (0PN) =", ADM_Mass_another_0PN )
    print( "ADM Mass (Omega expansion) (1PN) =", ADM_Mass_another_1PN )
    print( "ADM Mass (Omega expansion) (1.5PN) =", ADM_Mass_another_15PN )
    print( "ADM Mass (Omega expansion) (2PN) =", ADM_Mass_another_2PN )
    print( "ADM Mass (Omega expansion) (2.5PN) =", ADM_Mass_another_25PN )
    print( "ADM Mass (Omega expansion) (3PN) =", ADM_Mass_another_3PN )
    print()

    ############################

    ## Compute derivative dH/dr using the chosen finite-difference method
    ## Using M_adm = M + H_circular (where M = M1 + M2)

    def dH_dr(r):
        
        # Use finite-difference numerical derivative
        dHdr_0PN, dHdr_1PN, dHdr_15PN, dHdr_2PN, dHdr_25PN, dHdr_3PN \
                 = derivative.first_order_derivative_multivalue( M_ADM, r, 0.05, "7-points 6-orders" )

        # Use sympy symbolic differentiation
        # Note: symbolic approach may raise an error here
        # dHdr_0PN, dHdr_1PN, dHdr_15PN, dHdr_2PN, dHdr_25PN, dHdr_3PN \
        #          = sympy.diff(M_ADM(r), r)
        
        return dHdr_0PN, dHdr_1PN, dHdr_15PN, dHdr_2PN, dHdr_25PN, dHdr_3PN
    
    ############################

    # For some reason, the finite-difference derivative is very inaccurate
    # Possibly due to round-off error
    '''
    dHdr_0PN, dHdr_1PN, dHdr_15PN, dHdr_2PN, dHdr_25PN, dHdr_3PN = dH_dr(R)

    print()
    print( "dH/dr (0PN) =", dHdr_0PN )
    print( "dH/dr (1PN) =", dHdr_1PN )
    print( "dH/dr (1.5PN) =", dHdr_15PN )
    print( "dH/dr (2PN) =", dHdr_2PN )
    print( "dH/dr (2.5PN) =", dHdr_25PN )
    print()
    '''

    ############################

    ## Compute dH/dr from the analytic expression
    ## Using M_adm = M + H_circular (where M = M1 + M2)

    dHdr_0PN, dHdr_1PN, dHdr_15PN, dHdr_2PN, dHdr_25PN, dHdr_3PN = dADM_dr(R)

    print()
    print( "dH/dr (0PN) =", dHdr_0PN )
    print( "dH/dr (1PN) =", dHdr_1PN )
    print( "dH/dr (1.5PN) =", dHdr_15PN )
    print( "dH/dr (2PN) =", dHdr_2PN )
    print( "dH/dr (2.5PN) =", dHdr_25PN )
    print( "dH/dr (3PN) =", dHdr_3PN )
    print()

    ########################################

    ## Compute the time derivative of the orbital separation
    ## Based on:
    ## Antoni Ramos-Buades, Sascha Husa, and Geraint Pratten
    ## "Simple procedures to reduce eccentricity of binary black hole simulations"
    ## arXiv:1810.00036 [gr-qc], Phys. Rev. D 99, 023000 (2019)

    ## The high-order dE_GW/dt terms in that paper are questionable,
    ## so use results from an alternative reference:
    ## Serguei Ossokine et al., "Comparing Post-Newtonian and Numerical-Relativity Precession Dynamics"
    ## arXiv:1502.01747 [gr-qc], Phys. Rev. D 92, 104028 (2015)

    ############################


    ## Use the 3PN value for the orbital angular frequency (values computed above)

    Omega = Omega_3PN

    ## Compute the gravitational-wave energy flux dE_GW/dt

    ## Following definitions in arXiv:1502.01747 and arXiv:1810.00036, set the following parameters
    ## eta = m_eta
    m_delta     = ( m1 - m2 ) / ( m1 + m2 )
    
    spin_chi_a_vector = (S1 - S2) / 2.0
    spin_chi_s_vector = (S1 + S2) / 2.0
    
    spin_chi_a_square = numpy.dot(spin_chi_a_vector, spin_chi_a_vector)
    spin_chi_s_square = numpy.dot(spin_chi_s_vector, spin_chi_s_vector)

    ## Choose the unit vector l to point along +z
    spin_chi_a_l  = ( S1[2] - S2[2]) / 2.0 
    spin_chi_s_l  = ( S1[2] + S2[2]) / 2.0

    Euler_gamma   = 0.5772156649

    mass    = m1 + m2
    Spin_l  = ( S1[2]*(m1**2) + S2[2]*(m2**2) ) / (mass**2)
    Sigma_l = ( m2*S2[2] - m1*S1[2] ) / mass

    ## In formulas from arXiv:1810.00036 a leading negative sign must be added manually
    dEGW_dt_0PN = - ( 32.0 / 5.0 ) * (m_eta**2) * ( Omega**(10.0/3.0) )

    dEGW_dt_correction_1PN  = - ( Omega**(2.0/3.0) ) * ( 35.0*m_eta/12.0 + 1247.0/336.0 )          \
                              + Omega * ( 4.0*math.pi - (5.0/4.0)*m_delta*Sigma_l - 4.0*Spin_l ) 

    dEGW_dt_correction_15PN = Omega**(4.0/3.0) \
                              * (   (65.0/18.0)    * (m_eta**2)   \
                                  + (9271.0/504.0) * m_eta        \
                                  - (44711.0/9072.0)              \
                                  -  (89.0/48.0) * m_delta * numpy.dot(spin_chi_a_vector, spin_chi_s_vector) \
                                  + (287.0/48.0) * m_delta * spin_chi_a_l * spin_chi_s_l                     \
                                  + (    287.0/96.0 -       12.0*m_eta )  * (spin_chi_a_l**2)                \
                                  + (    287.0/96.0 + (1.0/24.0)*m_eta )  * (spin_chi_s_l**2)                \
                                  + ( - (89.0/96.0) +        4.0*m_eta )  * spin_chi_a_square                \
                                  - (    89.0/96.0  + (7.0/24.0)*m_eta )  * spin_chi_s_square                \
                                )
    dEGW_dt_correction_2PN_part1 = Omega**(5.0/3.0) \
                                    * ( - math.pi           * ( (583.0/24.0)*m_eta + 8191.0/672.0 ) \
                                        + m_delta * Sigma_l * (   (43.0/4.0)*m_eta - (13.0/16.0)  ) \
                                        + Spin_l            * (  (272.0/9.0)*m_eta - 4.5          ) \
                                      ) 

    ## The following are results from arXiv:1502.01747; they differ from arXiv:1810.00036
    dEGW_dt_correction_2PN_part2 = (Omega**2) \
                                    * ( -    (775.0/324.0) * (m_eta**3)                                \
                                        - (94403.0/3024.0) * (m_eta**2)                                \
                                        + ( - (134543.0/7776.0) + (41.0/48.0)*(math.pi**2) ) * m_eta   \
                                        + (6643739519.0/69854400.0)  + (16.0/3.0) * (math.pi**2)       \
                                        + (1712.0/105.0) * ( - Euler_gamma                             \
                                                             - 0.5*math.log(16.0*(Omega**(2.0/3.0)))   \
                                                           )                                           \
                                        - (31.0/6.0) * math.pi * m_delta * Sigma_l                     \
                                        -       16.0 * math.pi * Spin_l                                \
                                      )
    ## The following terms are from arXiv:1810.00036; they do not fully agree with arXiv:1502.01747
    dEGW_dt_correction_2PN_part2b = (Omega**2) \
                                    * ( - 4843497781.0/69854400.0                                  \
                                        - (775.0/324.0) * (m_eta**3)                               \
                                        - (94403.0/3024.0) * (m_eta**2)                            \
                                        + ( 8009293.0/54432.0 - (41.0/64.0)*(math.pi**2) ) * m_eta \
                                        + (287.0/192.0) * (math.pi**2)                             \
                                        + (1712.0/105.0) * ( - Euler_gamma                             \
                                                             + (35.0/107.0)*(math.pi**2)               \
                                                             - 0.5*math.log(16.0*(Omega**(2.0/3.0))) ) \
                                        - (31.0/6.0) * math.pi * m_delta * Spin_l                      \
                                        - 16.0 * math.pi * Spin_l                                      \
                                        + m_delta * spin_chi_a_l * spin_chi_s_l * (   611.0/252.0          \
                                                                                    - (809.0/18.0)*m_eta   \
                                                                                  )                        \
                                        + spin_chi_a_square * (   43.0*(m_eta**2)          \
                                                                - (8345.0/504.0)*m_eta     \
                                                                + 611.0/504.0              \
                                                              )                            \
                                        + spin_chi_s_square * (   (173.0/18.0)*(m_eta**2)  \
                                                                - (2393.0/72.0)*m_eta      \
                                                                + 611.0/504.0              \
                                                              )                            \
                                      )

    dEGW_dt_correction_25PN = Omega**(7.0/3.0)     \
                              * (   ( (1933585.0/3024.0)*(m_eta**2) + (214745.0/1728.0)*m_eta - (16258.0/504.0)   ) * math.pi           \
                                  + (    - (2810.0/27.0)*(m_eta**2) +    (6172.0/189.0)*m_eta + (476645.0/6784.0) ) * Spin_l            \
                                  + (    - (1501.0/36.0)*(m_eta**2) +    (1849.0/126.0)*m_eta + (9535.0/336.0)    ) * m_delta * Sigma_l \
                                )             

    dEGW_dt_correction_3PN  = Omega**(8.0/3.0) \
                              * (   ( - (7163.0/672.0) + (130583.0/2016.0) * m_eta ) * math.pi * m_delta * Sigma_l \
                                  + ( - (3485.0/96.0)  +    (13879.0/72.0) * m_eta ) * math.pi * Spin_l            \
                                )

    dEGW_dt_1PN  = dEGW_dt_0PN * ( 1.0 + dEGW_dt_correction_1PN  )
    dEGW_dt_15PN = dEGW_dt_0PN * ( 1.0 + dEGW_dt_correction_1PN       \
                                       + dEGW_dt_correction_15PN 
                                 )
    dEGW_dt_2PN  = dEGW_dt_0PN * ( 1.0 + dEGW_dt_correction_1PN        \
                                       + dEGW_dt_correction_15PN       \
                                       + dEGW_dt_correction_2PN_part1  \
                                       + dEGW_dt_correction_2PN_part2b \
                                 )
    dEGW_dt_25PN = dEGW_dt_0PN * ( 1.0 + dEGW_dt_correction_1PN        \
                                       + dEGW_dt_correction_15PN       \
                                       + dEGW_dt_correction_2PN_part1  \
                                       + dEGW_dt_correction_2PN_part2b \
                                       + dEGW_dt_correction_25PN       \
                                 )
    
    dEGW_dt_3PN  = dEGW_dt_0PN * ( 1.0 + dEGW_dt_correction_1PN        \
                                       + dEGW_dt_correction_15PN       \
                                       + dEGW_dt_correction_2PN_part1  \
                                       + dEGW_dt_correction_2PN_part2b \
                                       + dEGW_dt_correction_25PN       \
                                       + dEGW_dt_correction_3PN        \
                                 )
    
    print(                                   )
    print( " dEGW/dt 0pn   = ", dEGW_dt_0PN  )
    print( " dEGW/dt 1pn   = ", dEGW_dt_1PN  )
    print( " dEGW/dt 1.5pn = ", dEGW_dt_15PN )
    print( " dEGW/dt 2pn   = ", dEGW_dt_2PN  )
    print( " dEGW/dt 2.5pn = ", dEGW_dt_25PN )
    print( " dEGW/dt 3pn   = ", dEGW_dt_3PN  )
    print(                                   )

    ## Compute dr/dt via dr/dt = (dE_GW/dt) / (dH_circular/dr)
    ## Using M_adm = M + H_circular (where M = M1 + M2)

    drdt_0PN  = dEGW_dt_0PN  / dHdr_0PN 
    drdt_1PN  = dEGW_dt_1PN  / dHdr_1PN
    drdt_15PN = dEGW_dt_15PN / dHdr_15PN
    drdt_2PN  = dEGW_dt_2PN  / dHdr_2PN
    drdt_25PN = dEGW_dt_25PN / dHdr_25PN
    drdt_3PN  = dEGW_dt_3PN  / dHdr_3PN

    print()
    print( "Radial velocity Vr = dr/dt (0PN) =", drdt_0PN )
    print( "Radial velocity Vr = dr/dt (1PN) =", drdt_1PN )
    print( "Radial velocity Vr = dr/dt (1.5PN) =", drdt_15PN )
    print( "Radial velocity Vr = dr/dt (2PN) =", drdt_2PN )
    print( "Radial velocity Vr = dr/dt (2.5PN) =", drdt_25PN )
    print( "Radial velocity Vr = dr/dt (3PN) =", drdt_3PN )
    print()

    '''
    ## Old formulas
    ## Based on:
    ## A. Gopakumar, Bala R. Iyer, and Sai Iyer
    ## "Second post-Newtonian gravitational radiation reaction for two-body systems: Nonspinning bodies"
    ## arXiv:gr-qc/9703075
    ## Lower accuracy
    drdt_0PN = - (64.0/5.0) * ( epsilon**3 ) * m_eta
    drdt_1PN = drdt_0PN * ( 1.0 - ((1751.0/336.0) + 7.0*m_eta/4.0) * epsilon )
    drdt_2PN = drdt_0PN * ( 1.0 - ((1751.0/336.0) + 7.0*m_eta/4.0) * epsilon                                \
                                + ( 303455.0/18144 + 40981.0*m_eta/2016.0 + m_eta**2.0/2.0 ) * (epsilon**2) \
                          )
    '''

    ## Set dr/dt at 3PN accuracy (can be modified later)

    drdt = drdt_3PN 
    print("Rate of change of binary separation dr/dt =", drdt )
    
    ##########################

    ## Use formulas from the following reference
    ## James Healy, Carlos O. Lousto, Hiroyuki Nakano, and Yosef Zlochower
    ## Post-Newtonian Quasicircular Initial Orbits for Numerical Relativity
    ## arXiv:1702.00872[qr-qc]
    ## Class. Quant. Grav. 34, 145011 (2017)

    factor_0PN = (1.0+m_q)**2 / m_q

    factor_correction_1PN  = - 0.5 * epsilon * ( 7.0*(m_q**2) + 15.0*m_q + 7.0 ) / m_q
    factor_correction_15PN = 0.0
    factor_correction_2PN  = + (1.0/8.0) * (epsilon**2)                                                 \
                               * ( 47.0*(m_q**4) + 229.0*(m_q**3) + 363.0*(m_q**2) + 229.0*m_q + 47.0 ) \
                               / ( m_q * ( (1.0+m_q)**2 ) ) 
    factor_correction_25PN = + 0.25 * (epsilon**2.5)                                                    \
                               * (   ( 12.0*(m_q**2) + 11.0*m_q + 4.0  ) * S2[2] / (1.0 + m_q)          \
                                   + (  4.0*(m_q**2) + 11.0*m_q + 12.0 ) * S1[2] / ( (1.0 + m_q)*m_q )  \
                                 )
    factor_correction_3PN  = epsilon**3 \
                             * ( - (1.0/16.0) * ( (math.pi)**2 ) \
                                 - (1.0/48.0) * ( 363.0*(m_q**6) + 2608.0*(m_q**5) + 7324.0*(m_q**4)          \
                                                  + 10161.0*(m_q**3) + 7324.0*(m_q**2) + 2608.0*m_q + 363.0   \
                                                ) / ( m_q * ( (1.0+m_q)**4 ) )                                \
                                 + 0.25 * (S2[0]**2) * m_q * ( 18.0*(m_q**2) + 6.0*m_q +  5.0 ) / ( (1.0+m_q)**2 )           \
                                 - 0.75 * (S2[1]**2) * m_q * (  3.0*(m_q**2) +     m_q +  1.0 ) / ( (1.0+m_q)**2 )           \
                                 - 0.75 * (S2[2]**2) * m_q * (  3.0*(m_q**2) +     m_q +  1.0 ) / ( (1.0+m_q)**2 )           \
                                 + 0.25 * (S1[0]**2)       * (  5.0*(m_q**2) + 6.0*m_q + 18.0 ) / ( m_q * ( (1.0+m_q)**2 ) ) \
                                 - 0.75 * (S1[1]**2)       * (       m_q**2  +     m_q +  3.0 ) / ( m_q * ( (1.0+m_q)**2 ) ) \
                                 - 0.75 * (S1[2]**2)       * (       m_q**2  +     m_q +  3.0 ) / ( m_q * ( (1.0+m_q)**2 ) ) \
                                 +         S1[0] * S2[0]   * (  3.0*(m_q**2) -     m_q +  3.0 ) / ( (1.0+m_q)**2 )           \
                                 - 0.5  *  S1[1] * S2[1]   * (  3.0*(m_q**2) - 2.0*m_q +  3.0 ) / ( (1.0+m_q)**2 )           \
                                 - 0.5  *  S1[2] * S2[2]   * (  3.0*(m_q**2) - 2.0*m_q +  3.0 ) / ( (1.0+m_q)**2 )           \
                               )

    factor_1PN  = factor_0PN + factor_correction_1PN  
    factor_15PN = factor_0PN + factor_correction_1PN + factor_correction_15PN
    factor_2PN  = factor_0PN + factor_correction_1PN + factor_correction_15PN  \
                             + factor_correction_2PN
    factor_25PN = factor_0PN + factor_correction_1PN + factor_correction_15PN  \
                             + factor_correction_2PN + factor_correction_25PN
    factor_3PN  = factor_0PN + factor_correction_1PN + factor_correction_15PN  \
                             + factor_correction_2PN + factor_correction_25PN  \
                             + factor_correction_3PN

    Numerator = drdt  - (epsilon**3.5) * ( - 0.25 * S1[0] * S1[1] * m_q      * (      m_q + 6.0  ) / ( (1.0+m_q)**4 ) \
                                           - 0.25 * S1[0] * S2[1] * (m_q**2) * (  6.0*m_q + 13.0 ) / ( (1.0+m_q)**4 ) \
                                           - 0.25 * S2[0] * S1[1] * m_q      * ( 13.0*m_q + 6.0  ) / ( (1.0+m_q)**4 ) \
                                           - 0.25 * S2[0] * S2[1] * (m_q**2) * (  6.0*m_q + 1.0  ) / ( (1.0+m_q)**4 ) \
                                         )

    print(                                                           )
    print( " dr/dt = ", drdt                                         )
    print( " Numerator     in Pr calculation  3pn   = ", Numerator   )
    print( " devide factor in Pr calculation  0pn   = ", factor_0PN  )
    print( " devide factor in Pr calculation  1pn   = ", factor_1PN  )
    print( " devide factor in Pr calculation  1.5pn = ", factor_15PN )
    print( " devide factor in Pr calculation  2pn   = ", factor_2PN  )
    print( " devide factor in Pr calculation  2.5pn = ", factor_25PN )
    print( " devide factor in Pr calculation  3pn   = ", factor_3PN  )
    print(                                                           )

    Pr_0PN  = drdt_0PN  / factor_0PN
    Pr_1PN  = drdt_1PN  / factor_1PN
    Pr_15PN = drdt_15PN / factor_15PN
    Pr_2PN  = drdt_2PN  / factor_2PN
    Pr_25PN = drdt_25PN / factor_25PN
    Pr_3PN  = Numerator / factor_3PN

    print()
    print( "Pr (0PN) =", Pr_0PN )
    print( "Pr (1PN) =", Pr_1PN )
    print( "Pr (1.5PN) =", Pr_15PN )
    print( "Pr (2PN) =", Pr_2PN )
    print( "Pr (2.5PN) =", Pr_25PN )
    print( "Pr (3PN) =", Pr_3PN )
    print()

    ########################################

    ## Compute the binary momenta and write to file
    
    ## Note
    ## To match AMSS-NCKU TwoPuncture input conventions
    ## Place the larger-mass black hole at +y and the smaller at -y
    ## Momenta satisfy P1 = [-|Pt|, -|Pr|]
    ##         P2 = [+|Pt|, +|Pr|]
    ## This places both holes on +y at t=0 and makes them rotate counter-clockwise

    momentum1[1] = - abs(Pr_3PN)
    momentum2[1] =   abs(Pr_3PN) 

    momentum1[0] = - abs(Pt_3PN) 
    momentum2[0] =   abs(Pt_3PN)  

    print()
    print()
    print( "Binary radial momentum magnitude |Pr| =", abs(Pr_3PN) )
    print( "Binary tangential momentum magnitude |Pt| =", abs(Pt_3PN) )
    print()
    print( "Binary momenta at t=0:" )
    print( f"P1 = {momentum1}  P2 = {momentum2}" )
    print()

    print()
    print( "Binary orbital momenta setup complete" )
    print()

    ##############################################################################################

    ## Write results to file for use by Einstein Toolkit and AMSS-NCKU

    # file1 = open( "BBH_parameter.output", "w" )
    file1 = open( os.path.join(input_data.File_directionary, "BBH_parameter.output"), "w")

    print(                                           file=file1 )
    print( "Binary orbital parameters",               file=file1 )
    print(                                           file=file1 )
    print( f"Binary masses:     M1 = {M1}  M2 = {M2}",  file=file1 )
    print( "Dimensionless mass ratio Q = M1/M2 =", m_q,    file=file1 )
    print( f"Binary dimensionless spins: S1 = {S1}  S2 = {S2}", file=file1 )
    print(                                           file=file1 )
    print( "Binary coordinates at t=0:",               file=file1 )
    print( "X1 = ", position1[0],                   file=file1 ) 
    print( "Y1 = ", position1[1],                   file=file1 )
    print( "X2 = ", position2[0],                   file=file1 ) 
    print( "Y2 = ", position2[1],                   file=file1 )
    print(                                           file=file1 )
    print( "Binary momenta at t=0:",               file=file1 )
    print( "Pr  = ", Pr_3PN,                        file=file1 ) 
    print( "Pt  = ", Pt_3PN,                        file=file1 )
    print( "PX1 = - |Pt| = ", momentum1[0],         file=file1 )
    print( "PY1 = - |Pr| = ", momentum1[1],         file=file1 ) 
    print( "PX2 = + |Pt| = ", momentum2[0],         file=file1 )
    print( "PY2 = + |Pr| = ", momentum2[1],         file=file1 ) 
    print(                                           file=file1 )

    file1.close()
    
    return momentum1, momentum2

##############################################################################################


##############################################################################################

## Call the function to compute orbital momenta

## generate_BBH_orbit_parameters( M1, M2, S1, S2, D0, e0 )

##############################################################################################

