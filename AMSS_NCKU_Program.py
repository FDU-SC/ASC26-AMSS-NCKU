
##################################################################
##
## AMSS-NCKU Numerical Relativity Startup Program
## Author: Xiaoqu
## 2024/03/19
## Modified: 2025/12/09
##
##################################################################


##################################################################

## Print program introduction

import print_information

print_information.print_program_introduction()

##################################################################

### Pre-run prompts
#
#print( " Simulation will be started, please confirm you have set the correct parameters in the script file " )
#print( " AMSS_NCKU_Input.py "                                                                                )
#print( " If parameters have been set correctly, press Enter to continue !!!  "                               )
#print( " If you have not set parameters, press Ctrl+C to abort the simulation and adjust the parameters "    )
#print( " in script file AMSS_NCKU_Input.py !!! "                                                             )
#     
### Wait for user input (press Enter) to proceed
#inputvalue = input()           
#print()


##################################################################

import AMSS_NCKU_Input as input_data

##################################################################

## Create directories to store program run data

import os
import shutil
import sys
import time

## Set the output directory according to the input file
File_directory = os.path.join(input_data.File_directory)   

## If the specified output directory exists, ask the user whether to continue
if os.path.exists(File_directory):
    print( " Output dictionary has been existed !!!  "                                                              )
    print( " If you want to overwrite the existing file directory, please input 'continue' in the terminal !! "     ) 
    print( " If you want to retain the existing file directory, please input 'stop' in the terminal to stop the "   ) 
    print( " simulation. Then you can reset the output dictionary in the input script file AMSS_NCKU_Input.py !!! " )
    print(                                                                                                          )
    ## Prompt whether to overwrite the existing directory
    while True:
        try:
            inputvalue = input()
            ## If the user agrees to overwrite, proceed and remove the existing directory
            if ( inputvalue == "continue" ):
                print( " Continue the calculation !!! " )
                print(                                  )
                break  
            ## If the user chooses not to overwrite, exit and keep the existing directory
            elif ( inputvalue == "stop" ):
                print( " Stop the calculation !!! "    )
                sys.exit() 
            ## If the user input is invalid, prompt again
            else:
                print( " Please input your choice !!! "                   )
                print( " Input 'continue' or 'stop' in the terminal !!! " )
        except ValueError:
            print( " Please input your choice !!! "                   )
            print( " Input 'continue' or 'stop' in the terminal !!! " )
        
## Remove the existing output directory if present
shutil.rmtree(File_directory, ignore_errors=True)

## Create the output directory
os.makedirs(File_directory, exist_ok=True)

## Copy the Python input file into the run directory
shutil.copy("AMSS_NCKU_Input.py", File_directory)

# Generate subdirectories to store various output files

output_directory = os.path.join(File_directory, "AMSS_NCKU_output")
os.makedirs(output_directory, exist_ok=True)

binary_results_directory = os.path.join(output_directory, input_data.Output_directory)
os.makedirs(binary_results_directory, exist_ok=True)

figure_directory = os.path.join(File_directory, "figure")
os.makedirs(figure_directory, exist_ok=True)

print( " Output directory has been generated "     )
print(                                            )


##################################################################

## Output related parameter information

import setup

## Print and save input parameter information
setup.print_input_data( File_directory )
setup.generate_AMSSNCKU_input()

#print(                                                                                           )
#print( " Please check whether the grid boxes and their resolution are appropriate "             )
#print( " If the grid boxes and their resolution are not set properly, press Ctrl+C to abort. "  )
#print( " Adjust the grid levels and the number of grid points per level before retrying. "      )
#print( " If the grid boxes and resolution are correct, press Enter to continue. "               )
#inputvalue = input()  ## Wait for user input (press Enter) to proceed
#print()

setup.print_puncture_information()


##################################################################

## Generate AMSS-NCKU program input files based on the configured parameters

print(                                                                                   )
print( " Generating the AMSS-NCKU input parfile for the ABE executable. " ) 
print(                                                                                   ) 

## Generate cgh-related input files from the grid information

import numerical_grid        

numerical_grid.append_AMSSNCKU_cgh_input()  

print(                                                                                 )
print( " The input parfile for AMSS-NCKU C++ executable file ABE has been generated. " )    
print( " However, the input relevant to TwoPuncture need to be appended later. "       )
print(                                                                                 )


##################################################################

## Plot the initial grid configuration

print(                                                      )
print( " Schematically plot the numerical grid structure. " ) 
print(                                                      )

numerical_grid.plot_initial_grid()


##################################################################

## Generate AMSS-NCKU macro files according to the numerical scheme and parameters

print(                                                                                   ) 
print( " Automatically generating the macro file for AMSS-NCKU C++ executable file ABE " ) 
print( " (Based on the finite-difference numerical scheme) "                             )
print(                                                                                   )

import generate_macrodef

generate_macrodef.generate_macrodef_h()
print( " AMSS-NCKU macro file macrodef.h has been generated. " )
     
generate_macrodef.generate_macrodef_fh()
print( " AMSS-NCKU macro file macrodef.fh has been generated. " )


##################################################################

# Compile the AMSS-NCKU program according to user requirements

# Prompt about compiling and running AMSS-NCKU
print(                                                         )
print( " Preparing to compile and run the AMSS-NCKU code as requested " )
print( " Compiling the AMSS-NCKU code based on the generated macro files " )
print(                                                         )
#inputvalue = input()           
#print()

AMSS_NCKU_source_path = "AMSS_NCKU_source"
AMSS_NCKU_source_copy = os.path.join(File_directory, "AMSS_NCKU_source_copy")

###############################

## If AMSS_NCKU source folder is missing, create it and prompt the user

# if not os.path.exists(destination_folder):
#     os.makedirs(destination_folder)

if not os.path.exists(AMSS_NCKU_source_path):
    os.makedirs(AMSS_NCKU_source_path)
    print( " The AMSS-NCKU source files are incomplete; copy all source files into ./AMSS_NCKU_source. " )
    print( " Press Enter to continue. " )
    ## Wait for user input (press Enter) to proceed
    inputvalue = input()
    
###############################

# Copy AMSS-NCKU source files to prepare for compilation
shutil.copytree(AMSS_NCKU_source_path, AMSS_NCKU_source_copy, dirs_exist_ok=True)    

# (Comment) Example: copy the src folder to destination
# shutil.copytree(src, dst)

# Copy the generated macro files into the AMSS_NCKU source folder

macrodef_h_path  = os.path.join(File_directory, "macrodef.h") 
macrodef_fh_path = os.path.join(File_directory, "macrodef.fh") 

shutil.copy2(macrodef_h_path,  AMSS_NCKU_source_copy)
shutil.copy2(macrodef_fh_path, AMSS_NCKU_source_copy)

# Notes on copying files:
# shutil.copy2 preserves file metadata such as modification time.
# If you only want to copy file contents without metadata, use shutil.copy.

###############################

# Compile related programs

import makefile_and_run

## Change working directory to the target source copy
os.chdir(AMSS_NCKU_source_copy)
 
## Build the main AMSS-NCKU executable (ABE or ABEGPU)
makefile_and_run.makefile_ABE()

## If the initial-data method is Ansorg-TwoPuncture, build the TwoPunctureABE executable
if (input_data.Initial_Data_Method == "Ansorg-TwoPuncture" ): 
    makefile_and_run.makefile_TwoPunctureABE()
    
###########################

## Change current working directory back up two levels
os.chdir('..')
os.chdir('..')

print()

##################################################################

## Copy the AMSS-NCKU executable (ABE/ABEGPU) to the run directory

if (input_data.GPU_Calculation == "no"):
    ABE_file = os.path.join(AMSS_NCKU_source_copy, "ABE")
elif (input_data.GPU_Calculation == "yes"):
    ABE_file = os.path.join(AMSS_NCKU_source_copy, "ABEGPU")

if not os.path.exists( ABE_file ):
    print(                                                                                                  )
    print( " Lack of AMSS-NCKU executable file ABE/ABEGPU; recompile AMSS_NCKU_source manually. " )
    print( " When recompilation is finished, press Enter to continue. " )
    ## Wait for user input (press Enter) to proceed
    inputvalue = input() 

## Copy the executable ABE (or ABEGPU) into the run directory
shutil.copy2(ABE_file, output_directory)

###########################

## If the initial-data method is TwoPuncture, copy the TwoPunctureABE executable to the run directory

TwoPuncture_file = os.path.join(AMSS_NCKU_source_copy, "TwoPunctureABE")

if (input_data.Initial_Data_Method == "Ansorg-TwoPuncture" ):

    if not os.path.exists( TwoPuncture_file ):
        print(                                                                                                                  )
        print( " Lack of AMSS-NCKU executable file TwoPunctureABE; recompile TwoPunctureABE in AMSS_NCKU_source. " ) 
        print( " When recompilation is finished, press Enter to continue. " )
        inputvalue = input() 

    ## Copy the TwoPunctureABE executable into the run directory
    shutil.copy2(TwoPuncture_file, output_directory)

###########################


##################################################################

## If the initial-data method is TwoPuncture, generate the TwoPuncture input files

if (input_data.Initial_Data_Method == "Ansorg-TwoPuncture" ):

    print(                                                  )
    print( " Initial data is chosen as Ansorg-TwoPuncture" )
    print(                                                 )
    
    print(                                                                                                         )
    print( " Automatically generating the input parfile for the TwoPunctureABE executable " )
    print(                                                                                                         ) 
    
    import generate_TwoPuncture_input
    
    generate_TwoPuncture_input.generate_AMSSNCKU_TwoPuncture_input()
    
    print(                                                                                              )
    print( " The input parfile for the TwoPunctureABE executable has been generated. " )
    print(                                                                                              )
    
    ## Generated AMSS-NCKU TwoPuncture input filename
    AMSS_NCKU_TwoPuncture_inputfile      = 'AMSS-NCKU-TwoPuncture.input'
    AMSS_NCKU_TwoPuncture_inputfile_path = os.path.join( File_directory, AMSS_NCKU_TwoPuncture_inputfile )
 
    ## Copy and rename the file
    shutil.copy2( AMSS_NCKU_TwoPuncture_inputfile_path, os.path.join(output_directory, 'TwoPunctureinput.par') )
    
    ###########################

    ## Run TwoPuncture to generate initial-data files
    
    start_time = time.time()  # Record start time

    print()
    ## print( " Ready to launch the AMSS-NCKU TwoPuncture executable; press Enter to continue. " )
    ## inputvalue = input()                    
    print()
    
    ## Change to the output (run) directory
    os.chdir(output_directory)

    ## Run the TwoPuncture executable
    makefile_and_run.run_TwoPunctureABE()
    
    ###########################
    
    ## Change current working directory back up two levels
    os.chdir('..')
    os.chdir('..')
    
##################################################################
    
## Update puncture data based on TwoPuncture run results

import renew_puncture_parameter
    
renew_puncture_parameter.append_AMSSNCKU_BSSN_input(File_directory, output_directory)


## Generated AMSS-NCKU input filename
AMSS_NCKU_inputfile      = 'AMSS-NCKU.input'
AMSS_NCKU_inputfile_path = os.path.join(File_directory, AMSS_NCKU_inputfile)
 
## Copy and rename the file
shutil.copy2( AMSS_NCKU_inputfile_path, os.path.join(output_directory, 'input.par') )


print(                                                                          )
print( " Successfully copy all AMSS-NCKU input parfile to target dictionary. " )  
print(                                                                         )
    

##################################################################

## Launch the AMSS-NCKU program

print()
## print(" Ready to launch AMSS-NCKU; press Enter to continue. ")
## inputvalue = input()           
print()

## Change to the run directory
os.chdir( output_directory )

makefile_and_run.run_ABE()

## Change current working directory back up two levels
os.chdir('..')
os.chdir('..')


end_time = time.time()                
elapsed_time = end_time - start_time  

##################################################################

## Copy some basic input and log files out to facilitate debugging

## Path to the file that stores calculation settings
AMSS_NCKU_error_file_path = os.path.join(binary_results_directory, "setting.par")
## Copy and rename the file for easier inspection
shutil.copy( AMSS_NCKU_error_file_path, os.path.join(output_directory, "AMSSNCKU_setting_parameter") )

## Path to the error log file
AMSS_NCKU_error_file_path = os.path.join(binary_results_directory, "Error.log")
## Copy and rename the error log
shutil.copy( AMSS_NCKU_error_file_path, os.path.join(output_directory, "Error.log") )

## Primary program outputs
AMSS_NCKU_BH_data         = os.path.join(binary_results_directory, "bssn_BH.dat"        )
AMSS_NCKU_ADM_data        = os.path.join(binary_results_directory, "bssn_ADMQs.dat"     )
AMSS_NCKU_psi4_data       = os.path.join(binary_results_directory, "bssn_psi4.dat"      )
AMSS_NCKU_constraint_data = os.path.join(binary_results_directory, "bssn_constraint.dat")
## copy and rename the file
shutil.copy( AMSS_NCKU_BH_data,         os.path.join(output_directory, "bssn_BH.dat"        ) )
shutil.copy( AMSS_NCKU_ADM_data,        os.path.join(output_directory, "bssn_ADMQs.dat"     ) )
shutil.copy( AMSS_NCKU_psi4_data,       os.path.join(output_directory, "bssn_psi4.dat"      ) )
shutil.copy( AMSS_NCKU_constraint_data, os.path.join(output_directory, "bssn_constraint.dat") )

## Additional program outputs
if (input_data.Equation_Class == "BSSN-EM"):
    AMSS_NCKU_phi1_data = os.path.join(binary_results_directory, "bssn_phi1.dat" )
    AMSS_NCKU_phi2_data = os.path.join(binary_results_directory, "bssn_phi2.dat" )
    shutil.copy( AMSS_NCKU_phi1_data, os.path.join(output_directory, "bssn_phi1.dat" ) )
    shutil.copy( AMSS_NCKU_phi2_data, os.path.join(output_directory, "bssn_phi2.dat" ) )
elif (input_data.Equation_Class == "BSSN-EScalar"):
    AMSS_NCKU_maxs_data = os.path.join(binary_results_directory, "bssn_maxs.dat" )
    shutil.copy( AMSS_NCKU_maxs_data, os.path.join(output_directory, "bssn_maxs.dat" ) )

##################################################################

## Plot the AMSS-NCKU program results

print(                                                                          )
print( " Plotting the txt and binary results data from the AMSS-NCKU simulation " ) 
print(                                                                          )


import plot_xiaoqu
import plot_GW_strain_amplitude_xiaoqu

## Plot black hole trajectory
plot_xiaoqu.generate_puncture_orbit_plot(   binary_results_directory, figure_directory )
plot_xiaoqu.generate_puncture_orbit_plot3D( binary_results_directory, figure_directory )

## Plot black hole separation vs. time
plot_xiaoqu.generate_puncture_distence_plot( binary_results_directory, figure_directory )

## Plot gravitational waveforms (psi4 and strain amplitude)
for i in range(input_data.Detector_Number):
    plot_xiaoqu.generate_gravitational_wave_psi4_plot( binary_results_directory, figure_directory, i )
    plot_GW_strain_amplitude_xiaoqu.generate_gravitational_wave_amplitude_plot( binary_results_directory, figure_directory, i )

## Plot ADM mass evolution
for i in range(input_data.Detector_Number):
    plot_xiaoqu.generate_ADMmass_plot( binary_results_directory, figure_directory, i )

## Plot Hamiltonian constraint violation over time
for i in range(input_data.grid_level):
    plot_xiaoqu.generate_constraint_check_plot( binary_results_directory, figure_directory, i )

## Plot stored binary data
plot_xiaoqu.generate_binary_data_plot( binary_results_directory, figure_directory )

print(                                                 )
print( f" This Program Cost = {elapsed_time} Seconds " )
print(                                                 )


##################################################################

print(                                                                                    )
print( " The AMSS-NCKU-Python simulation is successfully finished, thanks for using !!! " )
print(                                                                                    )

##################################################################


