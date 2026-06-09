
##################################################################
##
## This file defines the commands used to build and run AMSS-NCKU
## Author: Xiaoqu
## 2025/01/24
##
##################################################################


import AMSS_NCKU_Input as input_data
import subprocess


##################################################################



##################################################################

## Compile the AMSS-NCKU main program ABE

def makefile_ABE():

    print(                                                        )
    print( " Compiling the AMSS-NCKU executable file ABE/ABEGPU " ) 
    print(                                                        )

    ## Build command
    if (input_data.GPU_Calculation == "no"):
        makefile_command  = "make -j4" + " ABE"
    elif (input_data.GPU_Calculation == "yes"):
        makefile_command  = "make -j4" + " ABEGPU"
    else:
        print( " CPU/GPU numerical calculation setting is wrong " )
        print(                                                    )
 
    ## Execute the command with subprocess.Popen and stream output
    makefile_process = subprocess.Popen(makefile_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)

    ## Read and print output lines as they arrive
    for line in makefile_process.stdout:
        print(line, end='')  # stream output in real time

    ## Wait for the process to finish
    makefile_return_code = makefile_process.wait()
    if makefile_return_code != 0:
        raise subprocess.CalledProcessError(makefile_return_code, makefile_command)
        
    print(                                                                  )
    print( " Compilation of the AMSS-NCKU executable file ABE is finished " ) 
    print(                                                                  )
    
    return
        
##################################################################



##################################################################

## Compile the AMSS-NCKU TwoPuncture program TwoPunctureABE

def makefile_TwoPunctureABE():

    print(                                                            )
    print( " Compiling the AMSS-NCKU executable file TwoPunctureABE " )
    print(                                                            )
    
    ## Build command
    makefile_command = "make" + " TwoPunctureABE"

    ## Execute the command with subprocess.Popen and stream output
    makefile_process = subprocess.Popen(makefile_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True) 
    
    ## Read and print output lines as they arrive
    for line in makefile_process.stdout:
        print(line, end='')  # stream output in real time
        
    ## Wait for the process to finish
    makefile_return_code = makefile_process.wait()
    if makefile_return_code != 0:
        raise subprocess.CalledProcessError(makefile_return_code, makefile_command)
        
    print(                                                                             )
    print( " Compilation of the AMSS-NCKU executable file TwoPunctureABE is finished " )
    print(                                                                             )
    
    return
    
##################################################################



##################################################################

## Run the AMSS-NCKU main program ABE

def run_ABE():

    print(                                                      )
    print( " Running the AMSS-NCKU executable file ABE/ABEGPU " ) 
    print(                                                      )

    ## Define the command to run; cast other values to strings as needed
    
    if (input_data.GPU_Calculation == "no"):
        mpi_command         = "mpirun -np " + str(input_data.MPI_processes) + " ./ABE"
        mpi_command_outfile = "ABE_out.log"
    elif (input_data.GPU_Calculation == "yes"):
        mpi_command         = "mpirun -np " + str(input_data.MPI_processes) + " ./ABEGPU"
        mpi_command_outfile = "ABEGPU_out.log"
 
    ## Execute the MPI command and stream output
    mpi_process = subprocess.Popen(mpi_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)

    ## Write ABE run output to file while printing to stdout
    with open(mpi_command_outfile, 'w') as file0:  
        ## Read and print output lines; also write each line to file
        for line in mpi_process.stdout:
            print(line, end='')  # stream output in real time
            file0.write(line)    # write the line to file
            file0.flush()        # flush to ensure each line is written immediately (optional)            
    file0.close()

    ## Wait for the process to finish
    mpi_return_code = mpi_process.wait()
    
    print(                                           )
    print( " The ABE/ABEGPU simulation is finished " ) 
    print(                                           )
    
    return

##################################################################



##################################################################

## Run the AMSS-NCKU TwoPuncture program TwoPunctureABE

def run_TwoPunctureABE():

    print(                                                          )
    print( " Running the AMSS-NCKU executable file TwoPunctureABE " ) 
    print(                                                          )
    
    ## Define the command to run
    TwoPuncture_command         = "./TwoPunctureABE"
    TwoPuncture_command_outfile = "TwoPunctureABE_out.log"

    ## Execute the command with subprocess.Popen and stream output
    TwoPuncture_process = subprocess.Popen(TwoPuncture_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)

    ## Write TwoPunctureABE run output to file while printing to stdout
    with open(TwoPuncture_command_outfile, 'w') as file0:  
        ## Read and print output lines; also write each line to file
        for line in TwoPuncture_process.stdout:
            print(line, end='')  # stream output in real time
            file0.write(line)    # write the line to file
            file0.flush()        # flush to ensure each line is written immediately (optional)                 
    file0.close()

    ## Wait for the process to finish
    TwoPuncture_command_return_code = TwoPuncture_process.wait()
    
    print(                                               )
    print( " The TwoPunctureABE simulation is finished " ) 
    print(                                               )
    
    return

##################################################################
    
