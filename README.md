# ASC26-AMSS-NCKU

Fudan University's optimization of [ASC26](https://www.asc-events.net/StudentChallenge/ASC26/final-competition.php)'s AMSS-NCKU Task. 

### Acknowledgements

[xiaoqu0000/NR-amssncku](https://github.com/xiaoqu0000/NR-amssncku), [ASC Committee](https://github.com/ASC-Competition/AMSS-NCKU)

## The original README:

# AMSS-NCKU

#### What can AMSS-NCKU do

AMSS - NCKU is a numerical relativity program developed in China, which is used to numerically solve Einstein's equations and calculate the change of the gravitational field over time.

AMSS - NCKU uses the finite difference method and the adaptive mesh refinement technique to achieve the numerical solution of Einstein's equations.

Currently, AMSS - NCKU can successfully handle binary black hole systems and multiple black hole systems, calculate the time evolution of these systems, and solve the gravitational waves released during these processes.

#### The Development of AMSS-NCKU

In 2008, the AMSS-NCKU code was successfully developed, enabling the numerical simulation for binary black hole and multiple black hole systems via the BSSN equations.

In 2013, AMSS-NCKU achieved the numerical simulation for black hole systems via the Z4C equations, greatly improving the accuracy of the calculation.

In 2015, AMSS-NCKU implemented hybrid CPU and GPU computing for the BSSN equations, improving the computational efficiency.

In 2024, we developed a Python operation interface for AMSS-NCKU to facilitate the freshman users and subsequent development.

#### Authors of AMSS-NCKU

Cao, Zhoujian (Beijing Normal University; Academy of Mathematics and Systems Science, Chinese Academy of Sciences; Hangzhou Institute for Advanced Study, University of Chinese Academy of Sciences)

Yo, Hwei-Jang (National Cheng Kung University)

Liu, Runqiu (Academy of Mathematics and Systems Science, Chinese Academy of Sciences)

Du, Zhihui (Tsinghua University)

Ji, Liwei (Rochester Institute of Technology)

Zhao, Zhichao (China Agricultural University)

Qiao, Chenkai (Chongqing University of Technology)

Yu, Jui-Ping (Former student)

Lin, Chun-Yu (Former student)

Zuo, Yi (Student)

#### Install the required packages and software that are prequisite to AMSS-NCKU code

Here, we take the Ubuntu 22.04 system as an example

1. Install the C++, Fortran, and Cuda compilers.

   $ sudo apt-get install gcc

   $ sudo apt-get install gfortran

   $ sudo apt-get install make

   $ sudo apt-get install build-essential

   $ sudo apt-get install nvidia-cuda-toolkit
2. Install the MPI tool

   $ sudo apt install openmpi-bin

   $ sudo apt install libopenmpi-dev
3. Install the Python3

   $ sudo apt-get install python3

   $ sudo apt-get install python3-pip
4. Install the required Python packages

   $ pip install numpy

   $ pip install scipy

   $ pip install matplotlib

   $ pip install SymPy

   $ pip install opencv-python-full

   $ pip install notebook

   $ pip install torch
5. Install the OpenCV tool

   $ sudo apt-get install libopencv-dev

   $ sudo apt-get install python-opencv

#### How to use AMSS-NCKU

0. Setting the parameters for compilation

   Modify the makefile.inc file in the AMSS_NCKU_source directory and change the settings according to your computer.

   The settings for the Ubuntu 22.04 system do not need to be modified.
1. Enter the AMSS-NCKU Python code folder and modify the input.

   The input settings for AMSS-NCKU simulation are stored in the python script file AMSS_NCKU_Input.py. Modify the parameters in this script file and save it.
2. Build the executable program and run the AMSS-NCKU simulation.

   Run the following command in the bash terminal.

   $ python AMSS_NCKU_Program.py

   or

   $ python3 AMSS_NCKU_Program.py

#### Update records

September 2025   First commit

December 2025    Update: Achieved the automatic plotting of gravitational wave amplitudes.

January 2026     Update: Fixed some bugs.

#### Tips

Due to limited testing, it's inevitable that there will be some unknown bugs in the code.

The computing time required for an actual evolution of a binary black hole system is relatively long. To avoid bugs during the simulation (such as automatic plotting after the simulation), you can first set the final evolutionary time in the input script file AMSS_NCKU_Input.py to 5M for testing.

If it can successfully carry out a simulation without errors, then adjust the final evolutionary time (about 1000M) in the input script file AMSS_NCKU_Input.py to start an actual simulation. This can reduce unnecessary waste of computing resources.

Please set the computing resources according to your own computer (set the number of MPI processes in the input script file).

#### Declaration

This code includes the C++ / Fortran codes from the original AMSS-NCKU code. A small number of functions are referenced from BAM.

Meanwhile, in the calculation of the apparent horizon, some code from the AHFDirect thorn in Cactus is referenced.
