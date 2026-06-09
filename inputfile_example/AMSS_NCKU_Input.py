
#################################################
##
## This file provides the input parameters required for numerical relativity.
## XIAOQU
## 2024/03/19 --- 2025/09/14
##
#################################################

import numpy    

#################################################

## Setting MPI processes and the output file directory

File_directory   = "GW150914"                    ## output file directory
Output_directory = "binary_output"               ## binary data file directory
                                                 ## The file directory name should not be too long
MPI_processes    = 8                             ## number of mpi processes used in the simulation

GPU_Calculation  = "no"                          ## Use GPU or not 
                                                 ## (prefer "no" in the current version, because the GPU part may have bugs when integrated in this Python interface)
CPU_Part         = 1.0
GPU_Part         = 0.0

#################################################


#################################################

## Setting the physical system and numerical method

Symmetry                 = "equatorial-symmetry"   ## Symmetry of System: choose equatorial-symmetry、no-symmetry、octant-symmetry
Equation_Class           = "BSSN"                  ## Evolution Equation: choose "BSSN", "BSSN-EScalar", "BSSN-EM", "Z4C" 
                                                   ## If "BSSN-EScalar" is chosen, it is necessary to set other parameters below
Initial_Data_Method      = "Ansorg-TwoPuncture"    ## initial data method: choose "Ansorg-TwoPuncture", "Lousto-Analytical", "Cao-Analytical", "KerrSchild-Analytical"
Time_Evolution_Method    = "runge-kutta-45"        ## time evolution method: choose "runge-kutta-45"
Finite_Diffenence_Method = "4th-order"             ## finite-difference method: choose "2nd-order", "4th-order", "6th-order", "8th-order"

#################################################


#################################################

## Setting the time evolutionary information

Start_Evolution_Time     = 0.0                    ## start evolution time t0
Final_Evolution_Time     = 1000.0                 ## final evolution time t1
Check_Time               = 100.0
Dump_Time                = 100.0                  ## time inteval dT for dumping binary data
D2_Dump_Time             = 100.0                  ## dump the ascii data for 2d surface after dT'
Analysis_Time            = 0.1                    ## dump the puncture position and GW psi4 after dT"
Evolution_Step_Number    = 10000000               ## stop the calculation after the maximal step number
Courant_Factor           = 0.5                    ## Courant Factor
Dissipation              = 0.15                   ## Kreiss-Oliger Dissipation Strength

#################################################


#################################################

## Setting the grid structure

basic_grid_set    = "Patch"                          ## grid structure: choose "Patch" or "Shell-Patch"
grid_center_set   = "Cell"                           ## grid center: chose "Cell" or "Vertex"

grid_level        = 9                                ## total number of AMR grid levels
static_grid_level = 5                                ## number of AMR static grid levels
moving_grid_level = grid_level - static_grid_level   ## number of AMR moving grid levels

analysis_level    = 0
refinement_level  = 4                                ## time refinement start from this grid level

largest_box_xyz_max = [320.0, 320.0, 320.0]          ## scale of the largest box
                                                     ## not ne cess ary to be cubic for "Patch" grid s tructure
                                                     ## need to be a cubic box for "Shell-Patch" grid structure
largest_box_xyz_min = - numpy.array(largest_box_xyz_max)  

static_grid_number = 96                              ## grid points of each static AMR grid (in x direction)
                                                     ## (grid points in y and z directions are automatically adjusted)
moving_grid_number = 48                              ## grid points of each moving AMR grid
shell_grid_number  = [32, 32, 100]                   ## grid points of Shell-Patch grid
                                                     ## in (phi, theta, r) direction
devide_factor      = 2.0                             ## resolution between different grid levels dh0/dh1, only support 2.0 now
                                                     

static_grid_type   = 'Linear'                        ## AMR static grid structure , only supports "Linear"
moving_grid_type   = 'Linear'                        ## AMR moving grid structure , only supports "Linear"

quarter_sphere_number = 96                           ## grid number of 1/4 s pher ical surface
                                                     ## (which is needed for evaluating the spherical surface integral)

#################################################


#################################################

## Setting the puncture information

puncture_number       = 2                                     

position_BH           = numpy.zeros( (puncture_number, 3) )   
parameter_BH          = numpy.zeros( (puncture_number, 3) )   
dimensionless_spin_BH = numpy.zeros( (puncture_number, 3) )   
momentum_BH           = numpy.zeros( (puncture_number, 3) )   

puncture_data_set     = "Manually"                       ## Method to give Puncture’s positions and momentum
                                                         ## choose "Manually" or "Automatically-BBH"
                                                         ## Prefer to choose "Manually", because "Automatically-BBH" is developing now

## initial orbital distance and ellipticity for BBHs system
## ( needed for "Automatically-BBH" case , not affect the "Manually" case )
Distance = 10.0
e0       = 0.0

## black hole parameter (M Q* a*)
parameter_BH[0] = [ 36.0/(36.0+29.0),  0.0,  +0.31 ]   
parameter_BH[1] = [ 29.0/(36.0+29.0),  0.0,  -0.46 ]  
## dimensionless spin in each direction
dimensionless_spin_BH[0] = [ 0.0,  0.0,  +0.31 ]   
dimensionless_spin_BH[1] = [ 0.0,  0.0,  -0.46 ]  

## use Brugmann's convention
##  -----0-----> y
##   -      +     

#---------------------------------------------

## If puncture_data_set is chosen to be "Manually", it is necessary to set the position and momentum of each puncture manually

## initial position for each puncture
position_BH[0]  = [  0.0,  10.0*29.0/(36.0+29.0), 0.0 ]  
position_BH[1]  = [  0.0, -10.0*36.0/(36.0+29.0), 0.0 ] 

## initial mumentum for each puncture
## (needed for "Manually" case, does not affect the "Automatically-BBH" case)
momentum_BH[0]  = [ -0.09530152296974252,  -0.00084541526517121,   0.0 ]
momentum_BH[1]  = [ +0.09530152296974252,  +0.00084541526517121,   0.0 ]


#################################################


#################################################

## Setting the gravitational wave information

GW_L_max        = 4                      ## maximal L number in gravitational wave
GW_M_max        = 4                      ## maximal M number in gravitational wave
Detector_Number = 12                     ## number of dector
Detector_Rmin   = 50.0                   ## nearest dector distance
Detector_Rmax   = 160.0                  ## farest dector distance

#################################################


#################################################

## Setting the apprent horizon

AHF_Find       = "no"                    ## whether to find the apparent horizon: choose "yes" or "no"

AHF_Find_Every = 24
AHF_Dump_Time  = 20.0

#################################################


#################################################

## Other parameters (testing)
## Only influence the Equation_Class = "BSSN-EScalar" case

FR_a2     = 3.0        ## f(R) = R + a2 * R^2    
FR_l2     = 10000.0
FR_phi0   = 0.00005
FR_r0     = 120.0
FR_sigma0 = 8.0
FR_Choice = 2          ## Choice options: 1 2 3 4 5
                       ## 1: phi(r) = phi0 * Exp(-(r-r0)**2/sigma0)   
                       ##    V(r)   = 0
                       ## 2: phi(r) =  phi0 * a2^2/(1+a2^2)  
                       ##    V(r)   = Exp(-8*Sqrt(PI/3)*phi(r)) * (1-Exp(4*Sqrt(PI/3)*phi(r)))**2 / (32*PI*a2)
                       ## 3: Schrodinger-Newton gived by system phi(r) 
                       ##    V(r)   = Exp(-8*Sqrt(PI/3)*phi(r)) * (1-Exp(4*Sqrt(PI/3)*phi(r)))**2 / (32*PI*a2)
                       ## 4: phi(r) = phi0 * 0.5 * ( tanh((r+r0)/sigma0) - tanh((r-r0)/sigma0) )  
                       ##    V(r)   = 0
                       ##    f(R)   = R + a2*R^2  with a2 = +oo
                       ## 5: phi(r) = phi0 * Exp(-(r-r0)**2/sigma)   
                       ##    V(r)   = 0

#################################################


#################################################

## Other parameters (testing)
## (please do not change if not necessary)

boundary_choice = "BAM-choice"     ## Sommerfeld boundary condition : choose "BAM-choice" or "Shibata-choice" 
                                   ## prefer "BAM-choice"

gauge_choice  = 0                  ## gauge choice
                                   ## 0: B^i gauge
                                   ## 1: David's puncture gauge
                                   ## 2: MB B^i gauge               
                                   ## 3: RIT B^i gauge
                                   ## 4: MB beta gauge 
                                   ## 5: RIT beta gauge 
                                   ## 6: MGB1 B^i gauge
                                   ## 7: MGB2 B^i gauge
                                   ## prefer 0 or 1
                                   
tetrad_type  = 2                   ## tetradtype 
                                   ##  v:r; u: phi; w: theta
                                   ##      v^a = (x,y,z)
                                   ## 0: orthonormal order: v,u,w
                                   ##    v^a = (x,y,z)   
                                   ##    m = (phi - i theta)/sqrt(2) 
                                   ##    following Frans, Eq.(8) of  PRD 75, 124018(2007)
                                   ## 1: orthonormal order: w,u,v
                                   ##    m = (theta + i phi)/sqrt(2) 
                                   ##    following Sperhake, Eq.(3.2) of  PRD 85, 124062(2012)    
                                   ## 2: orthonormal order: v,u,w
                                   ##    v_a = (x,y,z)
                                   ##    m = (phi - i theta)/sqrt(2) 
                                   ##    following Frans, Eq.(8) of  PRD 75, 124018(2007)
                                   ## this version recommend set to 2
                                   ## prefer 2
                                   
#################################################
                                   
