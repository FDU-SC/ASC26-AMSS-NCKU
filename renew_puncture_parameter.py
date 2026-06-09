
##################################################################
##
## Update puncture parameters from TwoPuncture output
## Author: Xiaoqu
## 2024/12/04
##
##################################################################

import AMSS_NCKU_Input as input_data
import numpy
import os

##################################################################



##################################################################

def read_TwoPuncture_Output(Output_File_directory):

    dimensionless_mass_BH = numpy.zeros( input_data.puncture_number )
    bare_mass_BH          = numpy.zeros( input_data.puncture_number )        ## initialize bare mass for each black hole
    position_BH           = numpy.zeros( (input_data.puncture_number, 3) )   ## initialize initial position for each black hole
    momentum_BH           = numpy.zeros( (input_data.puncture_number, 3) )   ## initialize momentum for each black hole
    angular_momentum_BH   = numpy.zeros( (input_data.puncture_number, 3) )   ## initialize spin angular momentum for each black hole
    
    # Read TwoPuncture output file
    data = numpy.loadtxt( os.path.join(Output_File_directory, "puncture_parameters_new.txt") )
    # Ensure data is parsed as a 1-D array
    data = data.reshape(-1)
    
    for i in range(input_data.puncture_number):
        
        ## Read parameters for the first two punctures from TwoPuncture output
        ## For additional punctures, read parameters from the input file
        if i<2:
            bare_mass_BH[i]          = data[12*i]
            dimensionless_mass_BH[i] = data[12*i+1]
            position_BH[i]           = [ data[12*i+3], data[12*i+4],  data[12*i+5]  ]
            momentum_BH[i]           = [ data[12*i+6], data[12*i+7],  data[12*i+8]  ]
            angular_momentum_BH[i]   = [ data[12*i+9], data[12*i+10], data[12*i+11] ]
        else:
            dimensionless_mass_BH[i] = input_data.parameter_BH[i,0]
            bare_mass_BH[i]          = input_data.parameter_BH[i,0]
            position_BH[i]           = input_data.position_BH[i]
            momentum_BH[i]           = input_data.momentum_BH[i]
            ## Read angular momentum according to symmetry
            if ( input_data.Symmetry == "equatorial-symmetry" ):
                angular_momentum_BH[i] = [ 0.0, 0.0, (input_data.parameter_BH[i,0]**2) * input_data.parameter_BH[i,2] ]
            elif ( input_data.Symmetry == "no-symmetry" ):
                angular_momentum_BH[i] = (dimensionless_mass_BH[i]**2) * input_data.dimensionless_spin_BH[i]
    
    return bare_mass_BH, dimensionless_mass_BH, position_BH, momentum_BH, angular_momentum_BH
    
##################################################################


##################################################################

## Append the computed puncture information into the AMSS-NCKU input file

def append_AMSSNCKU_BSSN_input(File_directory, TwoPuncture_File_directory): 

    charge_Q_BH = numpy.zeros( input_data.puncture_number )   ## initialize charge for each black hole

    ## If using Ansorg-TwoPuncture to solve the initial-data problem, read
    ## bare masses, positions and angular momenta from TwoPuncture output
    if (input_data.Initial_Data_Method == "Ansorg-TwoPuncture" ):
        bare_mass_BH, dimensionless_mass_BH, position_BH, momentum_BH, angular_momentum_BH = read_TwoPuncture_Output(TwoPuncture_File_directory)
        # set charge for each black hole
        for i in range(input_data.puncture_number):
            charge_Q_BH[i] = dimensionless_mass_BH[i] * input_data.parameter_BH[i,1]
    
    ## If using another method for initial data, read parameters directly from input
    else:
        position_BH = input_data.position_BH
        momentum_BH = input_data.momentum_BH
        ## angular_momentum_BH = input_data.angular_momentum_BH
        angular_momentum_BH = numpy.zeros( (input_data.puncture_number, 3) )   ## initialize spin angular momentum array
        mass_BH             = numpy.zeros( input_data.puncture_number      )   ## initialize mass array

        ## Set charge and spin angular momentum for each puncture
        for i in range(input_data.puncture_number):

            if ( input_data.Symmetry == "octant-symmetry" ):
                mass_BH[i]             = input_data.parameter_BH[i,0]
                charge_Q_BH[i]         = mass_BH[i]* input_data.parameter_BH[i,1]
                angular_momentum_BH[i] = [ 0.0, 0.0, (mass_BH[i]**2) * input_data.parameter_BH[i,2] ]
            elif ( input_data.Symmetry == "equatorial-symmetry" ):
                mass_BH[i]             = input_data.parameter_BH[i,0]
                charge_Q_BH[i]         = mass_BH[i]* input_data.parameter_BH[i,1]
                angular_momentum_BH[i] = [ 0.0, 0.0, (mass_BH[i]**2) * input_data.parameter_BH[i,2] ]
            elif ( input_data.Symmetry == "no-symmetry" ):
                mass_BH[i]             = input_data.parameter_BH[i,0]
                angular_momentum_BH[i] = (mass_BH[i]**2) * input_data.dimensionless_spin_BH[i]
                charge_Q_BH[i]         = mass_BH[i]      * input_data.parameter_BH[i,1]

    file1 = open( os.path.join(input_data.File_directory, "AMSS-NCKU.input"), "a")   ## open file in append mode

    ## Output BSSN related settings
    
    print(                                                                           file=file1 )
    print( "BSSN::chitiny  = 1e-5",                                                  file=file1 ) 
    print( "BSSN::time refinement start from level = ", input_data.refinement_level, file=file1 )
    print( "BSSN::BH_num   =  ",                        input_data.puncture_number,  file=file1 )
    
    for i in range(input_data.puncture_number):
    
        if (input_data.Initial_Data_Method == "Ansorg-TwoPuncture" ):
            print( f"BSSN::Mass[{i}]  = { bare_mass_BH[i] } ",      file=file1 )
        else:
            print( f"BSSN::Mass[{i}]  = { mass_BH[i] } ",           file=file1 )
            
        print( f"BSSN::Qchar[{i}] = { charge_Q_BH[i] } ",           file=file1 )
        print( f"BSSN::Porgx[{i}] = { position_BH[i,0] } ",         file=file1 )
        print( f"BSSN::Porgy[{i}] = { position_BH[i,1] } ",         file=file1 )
        print( f"BSSN::Porgz[{i}] = { position_BH[i,2] } ",         file=file1 )
        print( f"BSSN::Pmomx[{i}] = { momentum_BH[i,0] } ",         file=file1 )
        print( f"BSSN::Pmomy[{i}] = { momentum_BH[i,1] } ",         file=file1 )
        print( f"BSSN::Pmomz[{i}] = { momentum_BH[i,2] } ",         file=file1 )
        print( f"BSSN::Spinx[{i}] = { angular_momentum_BH[i,0] } ", file=file1 )
        print( f"BSSN::Spiny[{i}] = { angular_momentum_BH[i,1] } ", file=file1 )
        print( f"BSSN::Spinz[{i}] = { angular_momentum_BH[i,2] } ", file=file1 )
            
    print(                                                          file=file1 )
    
    file1.close()

    return
    
#################################################

