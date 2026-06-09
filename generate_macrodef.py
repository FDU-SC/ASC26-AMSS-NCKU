
##################################################################
##
## This file generates the macro definition files required by AMSS-NCKU
## Author: Xiaoqu
## 2024/12/01
##
##################################################################


import os 
import AMSS_NCKU_Input as input_data          ## import program input file


##################################################################

## Generate the macro file macrodef.h according to user settings
    
def generate_macrodef_h(): 

    file1 = open( os.path.join(input_data.File_directory, "macrodef.h"), "w") 

    print(                              file=file1 )
    print( "#ifndef MICRODEF_H",        file=file1 )
    print( "#define MICRODEF_H",        file=file1 )
    print(                              file=file1 )
    print( '#include "macrodef.fh"  ',  file=file1 )
    print(                              file=file1 ) 
    print( "// application parameters", file=file1 )
    print(                              file=file1 )

    # Define boundary-condition related macro SommerType
    # Sommerfeld boundary type
    # 0: BAM
    # 1: Shibata

    if ( input_data.boundary_choice == "BAM-choice" ):
        print( "#define SommerType 0",  file=file1 )
        print(                          file=file1 )
    elif ( input_data.boundary_choice == "Shibata-choice" ):
        print( "#define SommerType 0",  file=file1 )
        print(                          file=file1 )
    else:
        print( "Sommerfeld boundary condition setting error!!",  file=file1 )
        print(                                                   file=file1 )
        print( "# Sommerfeld boundary condition SommerType error!!" )
        print()

    # Define macro related to integration at infinity: GaussInt
    # Use Gauss-Legendre quadrature in theta direction
    
    print( "#define GaussInt",      file=file1 )
    print(                          file=file1 )

    # Define physical-system macro ABEtype
    # 0: BSSN vacuum
    # 1: coupled to scalar field
    # 2: Z4c vacuum
    # 3: coupled to Maxwell field

    if ( input_data.Equation_Class == "BSSN" ):
        print( "#define ABEtype 0", file=file1 )
        print(                      file=file1 )
    elif ( input_data.Equation_Class == "BSSN-EScalar" ):
        print( "#define ABEtype 1", file=file1 )
        print(                      file=file1 )
    elif ( input_data.Equation_Class == "BSSN-EM" ):
        print( "#define ABEtype 3", file=file1 )
        print(                      file=file1 )
    elif ( input_data.Equation_Class == "Z4C" ):
        print( "#define ABEtype 2", file=file1 )
        print(                      file=file1 )
    else:
        print( "Equation_Class setting error!!!"                )
        print()
        print( "# Equation type #define ABEtype setting error!!!", file=file1 )
        print(                                                   file=file1 )

    # Define macro With_AHF
    # enables Apparent Horizon Finder
    
    if (input_data.AHF_Find == "yes"): 
        print( "#define With_AHF",   file=file1 )
    elif (input_data.AHF_Find == "no"): 
        print( "//#define With_AHF", file=file1 )
    else:
        print( "AHF_Find input setting error!!!"            )
        print( "#define With_AHF setting error!!!", file=file1 )

    # Define macro Psi4type for Psi4 calculation method
    # 0: EB method
    # 1: 4-D method
    
    print( "#define Psi4type 0",    file=file1 )
    print(                          file=file1 )

    # Define macro Point_Psi4
    # control whether to use point-wise psi4
    
    print( "//#define Point_Psi4",  file=file1 )
    print(                          file=file1 )

    # Define macro RPS
    # RestrictProlong in Step (0) or after Step (1)
    
    print( "#define RPS 1",         file=file1 )
    print(                          file=file1 )

    # Define macro AGM
    # Enforce algebraic constraints
    # for every RK4 sub-step: 0
    # only when iter_count == 3: 1
    # after routine Step: 2
    
    print( "#define AGM 0",         file=file1 )
    print(                          file=file1 )

    # Define macro RPB
    # Restrict Prolong using BAM style (1) or old style (0)
    
    print( "#define RPB 0",         file=file1 )
    print(                          file=file1 )

    # Define macro MAPBH
    # 1: move analysis out of 4 sub-steps and treat PBH with Euler method
    
    print( "#define MAPBH 1",       file=file1 )
    print(                          file=file1 )

    # Define macro PSTR (parallel structure)
    # 0: level-by-level
    # 1: consider all levels
    # 2: same as 1 but reverse CPU order
    # 3: Frank's scheme
    
    print( "#define PSTR 0",        file=file1 )
    print(                          file=file1 )

    # Define macro REGLEV
    # regrid for every level or for all levels at once
    # 0: regrid each level separately
    # 1: regrid all levels together
    
    print( "#define REGLEV 0",      file=file1 )
    print(                          file=file1 )

    # Define macro USE_GPU
    # use GPU or not
    
    if ( input_data.GPU_Calculation == "yes"):
        print( "#define USE_GPU",   file=file1 )
        print(                      file=file1 )
    elif ( input_data.GPU_Calculation == "no"):
        print( "//#define USE_GPU", file=file1 )
        print(                      file=file1 )
    else:
        print( "CPU/GPU computation type setting error!!!"                               )
        print()
        print( "# CPU/GPU computation type #define USE_GPU setting error!!!", file=file1 )
        print(                                                    file=file1 )

    # Define macro CHECKDETAIL
    # enable checkpointing details for every process
    
    print( "//#define CHECKDETAIL", file=file1 )
    print(                          file=file1 )

    # Define macro FAKECHECK
    # use FakeCheckPrepare to write checkpoints
    
    print( "//#define FAKECHECK",   file=file1 )
    print(                          file=file1 )

    print( "//",                                                                         file=file1 )
    print( "// define SommerType",                                                       file=file1 )
    print( "//     sommerfeld boundary type",                                            file=file1 )
    print( "//     0: bam",                                                              file=file1 )
    print( "//     1: shibata",                                                          file=file1 )
    print( "//",                                                                         file=file1 )
    print( "// define GaussInt",                                                         file=file1 )
    print( "//     for Using Gauss-Legendre quadrature in theta direction",              file=file1 )
    print( "//",                                                                         file=file1 )
    print( "// define ABEtype",                                                          file=file1 )
    print( "//     0: BSSN vacuum",                                                      file=file1 )
    print( "//     1: coupled to scalar field",                                          file=file1 )
    print( "//     2: Z4c vacuum",                                                       file=file1 )
    print( "//     3: coupled to Maxwell field",                                         file=file1 )
    print( "//",                                                                         file=file1 )
    print( "// define With_AHF",                                                         file=file1 )
    print( "//     using Apparent Horizon Finder",                                       file=file1 )
    print( "//",                                                                         file=file1 )
    print( "// define Psi4type",                                                         file=file1 )
    print( "//     Psi4 calculation method",                                             file=file1 )
    print( "//     0: EB method",                                                        file=file1 )
    print( "//     1: 4-D method",                                                       file=file1 )
    print( "//",                                                                         file=file1 )
    print( "// define Point_Psi4",                                                       file=file1 )
    print( "//     for Using point psi4 or not",                                         file=file1 )
    print( "//",                                                                         file=file1 )
    print( "// define RPS",                                                              file=file1 )
    print( "//     RestrictProlong in Step (0) or after Step (1)",                       file=file1 )
    print( "//",                                                                         file=file1 )
    print( "// define AGM",                                                              file=file1 )
    print( "//     Enforce algebra constraint",                                          file=file1 )
    print( "//     for every RK4 sub step: 0",                                           file=file1 )
    print( "//     only when iter_count == 3: 1",                                        file=file1 )
    print( "//     after routine Step: 2",                                               file=file1 )
    print( "//",                                                                         file=file1 ) 
    print( "// define RPB",                                                              file=file1 )
    print( "//     Restrict Prolong using BAM style 1 or old style 0",                   file=file1 )
    print( "//",                                                                         file=file1 ) 
    print( "// define MAPBH",                                                            file=file1 )
    print( "//     1: move Analysis out ot 4 sub steps and treat PBH with Euler method", file=file1 )
    print( "//",                                                                         file=file1 ) 
    print( "// define PSTR",                                                             file=file1 )
    print( "//     parallel structure",                                                  file=file1 )
    print( "//     0: level by level",                                                   file=file1 )
    print( "//     1: considering all levels",                                           file=file1 )
    print( "//     2: as 1 but reverse the CPU order",                                   file=file1 )
    print( "//     3: Frank's scheme",                                                   file=file1 )
    print( "//",                                                                         file=file1 )
    print( "// define REGLEV",                                                           file=file1 )
    print( "//     regrid for every level or for all levels at a time",                  file=file1 )
    print( "//     0: for every level;",                                                 file=file1 ) 
    print( "//     1: for all",                                                          file=file1 )
    print( "//",                                                                         file=file1 )
    print( "// define USE_GPU",                                                          file=file1 )
    print( "//     use gpu or not",                                                      file=file1 )
    print( "//",                                                                         file=file1 )
    print( "// define CHECKDETAIL",                                                      file=file1 )
    print( "//     use checkpoint for every process",                                    file=file1 )
    print( "//",                                                                         file=file1 )
    print( "// define FAKECHECK",                                                        file=file1 )
    print( "//     use FakeCheckPrepare to write CheckPoint",                            file=file1 )
    print( "//",                                                                         file=file1 )

    print(                                                                         file=file1 )
    print( "////================================================================", file=file1 )
    print( "//  some basic parameters for numerical calculation",                  file=file1 )
    print( "////================================================================", file=file1 )
    print(                                                                         file=file1 )

    # Define macros related to dimensionality
    
    print( "#define dim 3",                              file=file1 )
    print(                                               file=file1 )

    # Macro Cell or Vertex is already defined in "macrodef.fh"
    
    print( '//#define Cell or Vertex in "macrodef.fh" ', file=file1 )
    print(                                               file=file1 )

    # Define macro buffer_width
    # number of buffer points for mesh-refinement interfaces
    
    print( "#define buffer_width 6",                     file=file1 )
    print(                                               file=file1 )

    # Define macro SC_width as buffer_width
    # number of buffer points for shell-box interface (on shell)
    
    print( "#define SC_width buffer_width",              file=file1 )
    print(                                               file=file1 )

    # Define macro CS_width
    # number of buffer points for shell-box interface (on box)
    
    print( "#define CS_width (2*buffer_width)",          file=file1 )
    print(                                               file=file1 )
    
    # The following are additional explanatory comments

    print( "//",                                                       file=file1 ) 
    print( '// define Cell or Vertex in "macrodef.fh" ',               file=file1 )    
    print( "//",                                                       file=file1 ) 
    print( "// define buffer_width",                                   file=file1 )
    print( "//     buffer point number for mesh refinement interface", file=file1 )
    print( "//",                                                       file=file1 ) 
    print( "// define SC_width buffer_width",                          file=file1 )
    print( "//     buffer point number shell-box interface, on shell", file=file1 )
    print( "//",                                                       file=file1 ) 
    print( "// define CS_width",                                       file=file1 )
    print( "//     buffer point number shell-box interface, on box",   file=file1 )
    print( "//",                                                       file=file1 ) 
    print(                                                             file=file1 )
    print( "#if(buffer_width < ghost_width)",                          file=file1 )
    print( "#   error we always assume buffer_width>ghost_width",      file=file1 )
    print( "#endif",                                                   file=file1 )
    print(                                                             file=file1 )
    
    print( "#define PACK 1",                               file=file1 )
    print( "#define UNPACK 2",                             file=file1 )
    print(                                                 file=file1 )
    print( "#define Mymax(a,b) (((a) > (b)) ? (a) : (b))", file=file1 )
    print( "#define Mymin(a,b) (((a) < (b)) ? (a) : (b))", file=file1 )
    print(                                                 file=file1 )
    print( "#define feq(a,b,d) (fabs(a-b)<d)",             file=file1 )
    print( "#define flt(a,b,d) ((a-b)<d)",                 file=file1 )
    print( "#define fgt(a,b,d) ((a-b)>d)",                 file=file1 )
    print(                                                 file=file1 )
    print( "#define TINY 1e-10",                           file=file1 )
    print(                                                 file=file1 )
    print( "#endif   /* MICRODEF_H */",                    file=file1 )
    print(                                                 file=file1 )
    
    file1.close()

    return file1
    
##################################################################


##################################################################

## Generate the macro file macrodef.fh according to user settings
    
def generate_macrodef_fh(): 

    file1 = open( os.path.join(input_data.File_directory, "macrodef.fh"), "w") 

    print( file=file1 )    
    
    # Define macro tetradtype
    
    # v:r; u: phi; w: theta
    
    # tetradtype 0
    # v^a = (x,y,z)
    # orthonormal order: v,u,w
    # m = (phi - i theta)/sqrt(2) following Frans, Eq.(8) of  PRD 75, 124018(2007)
    
    # tetradtype 1
    # orthonormal order: w,u,v
    # m = (theta + i phi)/sqrt(2) following Sperhake, Eq.(3.2) of  PRD 85, 124062(2012)
    
    # tetradtype 2
    # v_a = (x,y,z)
    # orthonormal order: v,u,w
    # m = (phi - i theta)/sqrt(2) following Frans, Eq.(8) of  PRD 75, 124018(2007)
    
    if ( input_data.tetrad_type == 0 ):
        print( "#define tetradtype 0", file=file1 )
        print(                         file=file1 ) 
    elif ( input_data.tetrad_type == 1 ):
        print( "#define tetradtype 1", file=file1 )
        print(                         file=file1 ) 
    elif ( input_data.tetrad_type == 2 ):
        print( "#define tetradtype 2", file=file1 )
        print(                         file=file1 ) 
    else:
        print( "tetradtype setting error!!" )
        print()
        print( "# tetradtype setting error!!", file=file1 )
        print(                             file=file1 ) 

    # Define macro for grid center: Cell or Vertex
    # Cell center or Vertex center

    if input_data.grid_center_set == "Cell":
        print( "#define Cell",   file=file1 )
        print(                   file=file1 )
    elif input_data.grid_center_set == "Vertex":
        print( "#define Vertex", file=file1 )
        print(                   file=file1 )
    else:
        print( "Grid_Center_Set setting error!!" )
        print()
        print( "# Grid center type #define Cell or #define Vertex setting error!", file=file1 )
        print(                                                            file=file1 )

    # Define macro ghost_width
    # 2nd order: 2
    # 4th order: 3
    # 6th order: 4
    # 8th order: 5
    
    if ( input_data.Finite_Diffenence_Method == "2nd-order" ):
        print( "#define ghost_width 2", file=file1 )
        print(                          file=file1 )
    elif ( input_data.Finite_Diffenence_Method == "4th-order" ):
        print( "#define ghost_width 3", file=file1 )
        print(                          file=file1 )
    elif ( input_data.Finite_Diffenence_Method == "6th-order" ):
        print( "#define ghost_width 4", file=file1 )
        print(                          file=file1 )
    elif ( input_data.Finite_Diffenence_Method == "8th-order" ):
        print( "#define ghost_width 5", file=file1 )
        print(                          file=file1 )
    else:
        print( "Finite_Difference_Method setting error!!!" )
        print()
        print( "# Finite_Difference_Method #define ghost_width setting error!!!",   file=file1 )
        print(                                                   file=file1 )

    # Whether to use a shell-patch grid
    # use shell or not

    if ( input_data.basic_grid_set == "Shell-Patch" ): 
        print( "#define WithShell", file=file1 )
        print(                      file=file1 )
    elif ( input_data.basic_grid_set == "Patch" ): 
        print(                      file=file1 )
    else:
        print( "basic_grid_set (grid type) setting error!!!" )
        print()
        print( "# grid type #define WithShell setting error!!!", file=file1 )
        print(                                               file=file1 )

    # Define macro CPBC
    # use constraint-preserving boundary conditions or not
    # only affects Z4c
    # CPBC requires WithShell

    if ( input_data.basic_grid_set == "Shell-Patch" ): 
        print( "#define CPBC", file=file1 )
        print(                 file=file1 )
    else:
        print(                 file=file1 )

    # Define gauge-related macros
    # Gauge condition type
    # 0: B^i gauge
    # 1: David's puncture gauge
    # 2: MB B^i gauge
    # 3: RIT B^i gauge
    # 4: MB beta gauge (beta gauge not means Eq.(3) of PRD 84, 124006)
    # 5: RIT beta gauge (beta gauge not means Eq.(3) of PRD 84, 124006)
    # 6: MGB1 B^i gauge
    # 7: MGB2 B^i gauge

    if ( input_data.gauge_choice == 0 ):
        print( "#define GAUGE 0",  file=file1 )
        print(                     file=file1 )
    elif ( input_data.gauge_choice == 1 ):
        print( "#define GAUGE 1",  file=file1 )
        print(                     file=file1 )
    elif ( input_data.gauge_choice == 2 ):
        print( "#define GAUGE 2",  file=file1 )
        print(                     file=file1 )
    elif ( input_data.gauge_choice == 3 ):
        print( "#define GAUGE 3",  file=file1 )
        print(                     file=file1 )
    elif ( input_data.gauge_choice == 4 ):
        print( "#define GAUGE 4",  file=file1 )
        print(                     file=file1 )
    elif ( input_data.gauge_choice == 5 ):
        print( "#define GAUGE 5",  file=file1 )
        print(                     file=file1 )
    elif ( input_data.gauge_choice == 6 ):
        print( "#define GAUGE 6",  file=file1 )
        print(                     file=file1 )
    elif ( input_data.gauge_choice == 7 ):
        print( "#define GAUGE 7",  file=file1 )
        print(                     file=file1 )
    else:
        print( "Gauge setting error!!" )
        print()
        print( "# Gauge GAUGE setting error!!", file=file1 )
        print(                             file=file1 )

    # Define macro CPBC_ghost_width
    # buffer points for CPBC boundary
    
    print( "#define CPBC_ghost_width  (ghost_width)", file=file1 )
    print(                                            file=file1 )

    # Define macro ABV
    # 0: use BSSN variables for constraint violation and psi4 calculation
    # 1: use ADM variables for constraint violation and psi4 calculation
    
    print( "#define ABV 0", file=file1 )
    print(                  file=file1 )
    
    # Define macro EScalar_CC related to F(R) scalar-tensor theory
    # 1: Case C of 1112.3928, V=0
    # 2: shell with   phi(r) = phi0 * a2^2/(1+a2^2), f(R) = R+a2*R^2 induced V
    # 3: ground state of Schrodinger-Newton system,  f(R) = R+a2*R^2 induced V
    # 4: a2 = +oo and phi(r) = phi0 * 0.5 * ( tanh((r+r0)/sigma) - tanh((r-r0)/sigma) )
    # 5: shell with   phi(r) = phi0 * Exp(-(r-r0)**2/sigma), V = 0
    
    if (input_data.Equation_Class == "BSSN-EScalar"):
        print( "#define EScalar_CC ", input_data.FR_Choice, file=file1 )
        print(                                              file=file1 )
    else:
        print( "#define EScalar_CC 2",                      file=file1 )
        print(                                              file=file1 )
        # For other calculations this value does not affect results; set a default
        # This prevents errors if FR_Choice is not present in AMSS_NCKU_Input.py
    
    # The following are explanatory comments

    print( "#if 0",                                                                                 file=file1 )
    print(                                                                                          file=file1 )
    print( "define tetradtype",                                                                     file=file1 )
    print( "    v:r; u: phi; w: theta",                                                             file=file1 )
    print( "    tetradtype 0",                                                                      file=file1 )
    print( "    v^a = (x,y,z)",                                                                     file=file1 )
    print( "    orthonormal order: v,u,w",                                                          file=file1 )
    print( "    m = (phi - i theta)/sqrt(2) following Frans, Eq.(8) of  PRD 75, 124018(2007)",      file=file1 )
    print( "    tetradtype 1",                                                                      file=file1 )
    print( "    orthonormal order: w,u,v",                                                          file=file1 )
    print( "    m = (theta + i phi)/sqrt(2) following Sperhake, Eq.(3.2) of  PRD 85, 124062(2012)", file=file1 )
    print( "    tetradtype 2",                                                                      file=file1 )
    print( "    v_a = (x,y,z)",                                                                     file=file1 )
    print( "    orthonormal order: v,u,w",                                                          file=file1 )
    print( "    m = (phi - i theta)/sqrt(2) following Frans, Eq.(8) of  PRD 75, 124018(2007)",      file=file1 )
    print(                                                                                          file=file1 )
    print( "define Cell or Vertex",                                                                 file=file1 )
    print( "    Cell center or Vertex center",                                                      file=file1 )
    print(                                                                                          file=file1 )
    print( "define ghost_width",                                                                    file=file1 )
    print( "    2nd order: 2",                                                                      file=file1 )
    print( "    4th order: 3",                                                                      file=file1 )
    print( "    6th order: 4",                                                                      file=file1 )
    print( "    8th order: 5",                                                                      file=file1 )
    print(                                                                                          file=file1 )
    print( "define WithShell",                                                                      file=file1 )
    print( "    use shell or not",                                                                  file=file1 )
    print(                                                                                          file=file1 )
    print( "define CPBC",                                                                           file=file1 )
    print( "    use constraint preserving boundary condition or not",                               file=file1 )
    print( "    only affect Z4c",                                                                   file=file1 )
    print( "    CPBC only supports WithShell",                                                      file=file1 )
    print(                                                                                          file=file1 )
    print( "define GAUGE",                                                                          file=file1 )
    print( "    0: B^i gauge",                                                                      file=file1 )
    print( "    1: David puncture gauge",                                                           file=file1 )
    print( "    2: MB B^i gauge",                                                                   file=file1 )
    print( "    3: RIT B^i gauge",                                                                  file=file1 )
    print( "    4: MB beta gauge (beta gauge not means Eq.(3) of PRD 84, 124006)",                  file=file1 )
    print( "    5: RIT beta gauge (beta gauge not means Eq.(3) of PRD 84, 124006)",                 file=file1 )
    print( "    6: MGB1 B^i gauge",                                                                 file=file1 )
    print( "    7: MGB2 B^i gauge",                                                                 file=file1 )
    print(                                                                                          file=file1 )
    print( "define CPBC_ghost_width  (ghost_width)",                                                file=file1 )
    print( "    buffer points for CPBC boundary",                                                   file=file1 )
    print(                                                                                          file=file1 )
    print( "define ABV",                                                                            file=file1 )
    print( "    0: using BSSN variable for constraint violation and psi4 calculation",              file=file1 )
    print( "    1: using ADM variable for constraint violation and psi4 calculation",               file=file1 )
    print(                                                                                          file=file1 )
    print( "define EScalar_CC",                                                                     file=file1 )
    print( "Type of Potential and Scalar Distribution in F(R) Scalar-Tensor Theory",                file=file1 )
    print( "    1: Case C of 1112.3928, V=0",                                                       file=file1 )
    print( "    2: shell with   phi(r) = phi0 * a2^2/(1+a2^2), f(R) = R+a2*R^2 induced V",          file=file1 )
    print( "    3: ground state of Schrodinger-Newton system,  f(R) = R+a2*R^2 induced V",          file=file1 )
    print( "    4: a2 = +oo and phi(r) = phi0 * 0.5 * ( tanh((r+r0)/sigma) - tanh((r-r0)/sigma) )", file=file1 )
    print( "    5: shell with   phi(r) = phi0 * Exp(-(r-r0)**2/sigma), V = 0",                      file=file1 )
    print(                                                                                          file=file1 )
    print( "#endif",                                                                                file=file1 )
    print(                                                                                          file=file1 )
    
    file1.close()

    return file1
    
##################################################################

    
