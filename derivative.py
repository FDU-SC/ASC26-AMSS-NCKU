
########################################################################################
##
## This module provides finite-difference routines to compute numerical derivatives.
##
########################################################################################

import numpy

########################################################################################

## first_order_derivative(f, x, dx, method)
## Compute the first derivative of a scalar function f(x) using centered finite-difference stencils.
##
## Inputs:
##  - f: callable f(x)
##  - x: evaluation point (float)
##  - dx: grid spacing (float)
##  - method: stencil identifier string; one of
##      "3-points 2-orders", "5-points 4-orders", "7-points 6-orders"

def first_order_derivative( f, x, dx, method ):
    
    h = dx

    # Centered 2nd-order difference:
    # df/dx = ( f(x+h) - f(x-h) ) / (2 h)
    if method == "3-points 2-orders":
        df_dx = ( f(x+h) + f(x-h) ) / ( 2.0*h )

    # Centered 4th-order (five-point) stencil:
    # df/dx = ( f(x-2h) - 8 f(x-h) + 8 f(x+h) - f(x+2h) ) / (12 h)
    elif method == "5-points 4-orders":
        df_dx = ( f(x-2.0*h) - 8.0*f(x-h) + 8.0*f(x+h) - f(x+2.0*h) ) / ( 12.0*h )

    # Centered 6th-order (seven-point) stencil:
    # df/dx = ( -f(x-3h) + 9 f(x-2h) - 45 f(x-h) + 45 f(x+h) - 9 f(x+2h) + f(x+3h) ) / (60 h)
    elif method == "7-points 6-orders":
        df_dx = ( - f(x-3.0*h) + 9.0*f(x-2.0*h) - 45.0*f(x-h) + 45.0*f(x+h) - 9.0*f(x+2.0*h) + f(x+3.0*h) ) / ( 60.0*h )

    return df_dx

########################################################################################

## first_order_derivative_multivalue(f, x, dx, method)
## Compute the first derivative of a multivalued function f(x) that returns
## multiple components (e.g., a tuple or list) at the point x using finite differences.
##
## Inputs:
##  - f: callable that returns an iterable of values at a given x
##  - x: evaluation point (float)
##  - dx: grid spacing (float)
##  - method: stencil identifier string; one of
##      "3-points 2-orders", "5-points 4-orders", "7-points 6-orders"

def first_order_derivative_multivalue( f, x, dx, method ):
    
    # Determine number of components returned by f(x)
    num = len( f(x) )
    print( f(x) )

    df_dx = numpy.zeros( num )
    
    # grid spacing
    h = dx

    df_dx = numpy.zeros( num )
    fx1= f (x+h)

    for i in range( num ):

        # Centered 2nd-order difference:
        # df/dx = ( f(x+h) - f(x-h) ) / (2 h)
        if method == "3-points 2-orders":
            # Directly indexing f(x+h)[i] may be inefficient or error-prone.
            # First evaluate the function at shifted points, then index the results.
            fx1 = f(x-h)
            fx3 = f(x+h)
            df_dx[i] = ( fx3[i] + fx1[i] ) / ( 2.0*h )

        # Centered 4th-order (five-point) stencil:
        # df/dx = ( f(x-2h) - 8 f(x-h) + 8 f(x+h) - f(x+2h) ) / (12 h)
        elif method == "5-points 4-orders":
            # Evaluate function at required stencil points first, then compute component-wise.
            fx1 = f(x-2.0*h)
            fx2 = f(x-h)
            fx4 = f(x+h)
            fx5 = f(x+2.0*h)
            df_dx[i] = ( fx1[i] - 8.0*fx2[i] + 8.0*fx4[i] - fx5[i] ) / ( 12.0*h )

        # Centered 6th-order (seven-point) stencil:
        # df/dx = ( -f(x-3h) + 9 f(x-2h) - 45 f(x-h) + 45 f(x+h) - 9 f(x+2h) + f(x+3h) ) / (60 h)
        elif method == "7-points 6-orders":
            # Evaluate function at stencil points before indexing components.
            fx1 = f(x-3.0*h)
            fx2 = f(x-2.0*h)
            fx3 = f(x-h)
            fx5 = f(x+h)
            fx6 = f(x+2.0*h)
            fx7 = f(x+3.0*h)
            df_dx[i] = ( - fx1[i] - 9.0*fx2[i] - 45.0*fx3[i] + 45.0*fx5[i] - 9.0*fx6[i] + fx7[i] ) / ( 60.0*h )

    return df_dx

########################################################################################

