
#################################################
##
## This file contains utilities to plot binary data produced by the
## numerical-relativity group (AMSS-NCKU).
## Author: Xiaoqu
## Dates: 2024/10/01 --- 2025/09/14
##
#################################################

import numpy
import matplotlib.pyplot    as     plt
from   matplotlib.colors    import LogNorm
from   mpl_toolkits.mplot3d import Axes3D
## import torch
import AMSS_NCKU_Input      as input_data

import os


#########################################################################################

def plot_binary_data( filename, binary_outdir, figure_outdir ):

    figure_title0 = filename.replace(binary_outdir + "/", "")  # remove directory prefix
    figure_title  = figure_title0.replace(".bin", "")          # remove .bin suffix
    
    print()
    print( "reading binary data from file =", figure_title0 )

###################################

    # Open file
    # Read binary array in the AMSS-NCKU output order
    with open(filename, 'rb') as file:

        physical_time = numpy.fromfile( file, dtype=numpy.float64, count=1 )
        nx, ny, nz    = numpy.fromfile( file, dtype=numpy.int32,   count=3 )
        xmin, xmax    = numpy.fromfile( file, dtype=numpy.float64, count=2 )
        ymin, ymax    = numpy.fromfile( file, dtype=numpy.float64, count=2 )
        zmin, zmax    = numpy.fromfile( file, dtype=numpy.float64, count=2 )
        data          = numpy.fromfile( file, dtype=numpy.float64          )
        
        # Now `data` array contains the binary data read from file
 
    print( "obtained data shape  =", data.shape ) 
    print( "obtained data size   =", data.size  ) 
    print( "obtained data points =", nx, "*", ny, "*", nz, "=", nx*ny*nz )
    
###################################

    # Reshape flat array into a multi-dimensional grid
    data_reshape = data.reshape( (nz, ny, nx) ) ## this ordering produces correct plots
    # print(data_reshape)

    # data1 = data_reshape[0,:,:]
    # print(data1)

    Rmin = [xmin, ymin, zmin] 
    Rmax = [xmax, ymax, zmax]
    N    = [nx, ny, nz]
    print( "coordinate minimum =", Rmin )
    print( "coordinate maximum =", Rmax )
    print( "grid point         =", N    )
    
    print()
    print( "Data file read successfully. Plotting data." )
    print()

    # Call plotting helper to produce plots
    figure_title0    = filename.replace(binary_outdir + "/", "") # remove directory prefix
    figure_title     = figure_title.replace(".bin", "")          # remove .bin suffix
    figure_title_new = figure_title[:-6]                            # strip trailing 6 characters (iteration label)
    
    get_data_xy( Rmin, Rmax, N, data_reshape, physical_time[0], figure_title_new, figure_outdir )
    # Note: numpy.fromfile returns an array for `physical_time` (even though
    # it contains a single element), so use `physical_time[0]` as the scalar time.
    
    # Explicitly delete large arrays to free memory
    del data
    del data_reshape
    
    print( "binary data file =", figure_title0, "plot has finished" )
    print(                                                             )

    return
    
    
#########################################################################################




####################################################################################

# Plot a single binary dataset (2D slices and 3D surface)

def get_data_xy( Rmin, Rmax, n, data0, time, figure_title, figure_outdir ):

    figure_contourplot_outdir = os.path.join(figure_outdir, "contour plot")
    figure_densityplot_outdir = os.path.join(figure_outdir, "density plot")
    figure_surfaceplot_outdir = os.path.join(figure_outdir, "surface plot")

    # Reconstruct coordinates from grid metadata
    x = numpy.linspace(Rmin[0], Rmax[0], n[0])
    y = numpy.linspace(Rmin[1], Rmax[1], n[1])
    z = numpy.linspace(Rmin[2], Rmax[2], n[2])
    print( " x = ", x )
    print( " y = ", y )
    print( " z = ", z )

    # Build 2D meshgrid for plotting
    # X, Y = numpy.meshgrid(x, y)                                
    # X, Y = torch.meshgrid(torch.tensor(x), torch.tensor(y))    
    Y, X = numpy.meshgrid(y, x)    
    
    # Notes on numpy.meshgrid:
    # If x has length nx and y has length ny, then X,Y = numpy.meshgrid(x,y)
    # produce arrays with shape (ny, nx). X has rows copied from x and
    # Y has columns copied from y.
    
    print( " X0 = ", X[:,0] )
    print( " Y0 = ", Y[0,:] )

    # Extract data on the central xy plane
    if input_data.Symmetry == "no-symmetry":
        data_xy = data0[n[2]//2,:,:]
    else:
        data_xy = data0[0,:,:]
        
    # The original data ordering was thought to be column-major; tests
    # indicate no transpose is required.
    
    # print( data_xy_0.shape )
    # print( data_xy.shape )
    
    # Define finer coordinate grids for interpolation
    x_new = numpy.linspace(Rmin[0], Rmax[0], int(2.5*n[0]))
    y_new = numpy.linspace(Rmin[1], Rmax[1], int(2.5*n[1]))
    z_new = numpy.linspace(Rmin[2], Rmax[2], int(2.5*n[2]))
    X_new, Y_new = numpy.meshgrid(x_new, y_new)
    
    # Interpolate data onto the finer grid
    data_xy_fit = scipy.interpolate.griddata( (X.flatten(), Y.flatten()), data_xy.flatten(), (X_new, Y_new), method="cubic" )

    # Plot 2D contour map
    fig, ax = plt.subplots()
    # contourf = ax.contourf(X, Y, data_xy, 8, cmap='coolwarm', norm=LogNorm(vmin=1, vmax=10), levels=numpy.logspace(-2, 2, 8))  # use 'coolwarm' colormap with LogNorm scaling
    # contourf = ax.contourf( X, Y, data_xy_0, cmap=plt.get_cmap('RdYlGn_r') )
    # contour  = ax.contour(  X, Y, data_xy_0, 8, colors='k', linewidths=0.5 )     # add contour lines
    # Use interpolated data for plotting
    contourf = ax.contourf( X_new, Y_new, data_xy_fit, cmap=plt.get_cmap('RdYlGn_r') )
    contour  = ax.contour(  X_new, Y_new, data_xy_fit, 8, colors='k', linewidths=0.5 )     # add contour lines
    cbar     = plt.colorbar(contourf)                                                      # add colorbar
    ax.set_title(  figure_title + "  physical time = " + str(time) )                       # set title and axis labels
    ax.set_xlabel( "X [M]" )
    ax.set_ylabel( "Y [M]" )
    # plt.show()                                                               # display figure
    plt.savefig( os.path.join(figure_contourplot_outdir, figure_title + " time = " + str(time) + " contour_plot.pdf") )   # save figure
    plt.close()
    
    # Plot 2D density (heat) map
    # fig1 = plt.figure()
    fig1, ax  = plt.subplots()
    # Tests show no transpose is necessary; however the y-axis appears
    # flipped in the image, so set extent accordingly.
    imshowfig = plt.imshow( data_xy, interpolation='bicubic', extent=[X.min(), X.max(), Y.max(), Y.min()] )  
    # ax.invert_xaxis()
    ax.invert_yaxis()                                                       # invert y-axis
    cbar      = plt.colorbar(imshowfig)                                     # add colorbar
    ax.set_title(  figure_title + "  physical time = " + str(time)  )      # set title and axis labels
    ax.set_xlabel( "X [M]" )
    ax.set_ylabel( "Y [M]" )
    # plt.show() 
    plt.savefig( os.path.join(figure_densityplot_outdir, figure_title + " time = " + str(time) + " density_plot.pdf") )
    plt.close()

    # Plot 3D surface
    fig2 = plt.figure()                                                       # create new figure
    ax = fig2.add_subplot( 111, projection='3d' )                             # 3D axes
    # plot interpolated surface
    ax.plot_surface( X_new, Y_new, data_xy_fit, cmap='viridis' )              # surface plot
    ax.set_title(  figure_title + "  physical time = " + str(time) )        # set title and labels
    ax.set_xlabel( "X [M]" )
    ax.set_ylabel( "Y [M]" )
    plt.savefig( os.path.join(figure_surfaceplot_outdir, figure_title + " time = " + str(time) + " surface_plot.pdf") )   # save image
    plt.close()

    return

####################################################################################

# Configure directories based on input configuration
File_directionary = os.path.join(input_data.File_directionary)

output_directionary = os.path.join(File_directionary, "AMSS_NCKU_output")
binary_results_directionary = os.path.join(output_directionary, input_data.Output_directionary)

figure_directionary = "figure"
if not os.path.exists(figure_directionary):
    os.mkdir(figure_directionary)

surface_plot_directionary = os.path.join(figure_directionary, "surface plot")
density_plot_directionary = os.path.join(figure_directionary, "density plot")
contour_plot_directionary = os.path.join(figure_directionary, "contour plot")
if not os.path.exists(surface_plot_directionary):
    os.mkdir(surface_plot_directionary)
if not os.path.exists(density_plot_directionary):
    os.mkdir(density_plot_directionary)
if not os.path.exists(contour_plot_directionary):
    os.mkdir(contour_plot_directionary)

filename = os.path.join(binary_results_directionary, 'Lev05-00_phi0_00154.bin')

plot_binary_data( filename, binary_results_directionary, figure_directionary )

