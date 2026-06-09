
############################################################################################
##
## This module provides finite-difference routines for numerical derivatives
##
############################################################################################



############################################################################################

## first_order_derivative(f, x, dx, method)
## Compute the first derivative of a scalar function f(x) using finite differences.
##
## Inputs:
##  - f: callable f(x)
##  - x: evaluation point (float)
##  - dx: grid spacing (float)
##  - method: string specifying the finite-difference stencil; one of
##      "3-points 2-orders", "5-points 4-orders", "7-points 6-orders", "9-points 8-orders"
##
def first_order_derivative( f, x, dx, method ):
    
    h = dx

    # Centered 2nd-order difference:
    # df/dx = ( f(x+h) - f(x-h) ) / (2 h)
    if method == "3-points 2-orders":
        df_dx = ( f(x+h) + f(x-h) ) / ( 2.0*h )

    # Centered 4th-order, five-point stencil:
    # df/dx = ( f(x-2h) - 8 f(x-h) + 8 f(x+h) - f(x+2h) ) / (12 h)
    elif method == "5-points 4-orders":
        df_dx = ( f(x-2.0*h) - 8.0*f(x-h) + 8.0*f(x+h) - f(x+2.0*h) ) / ( 12.0*h )

    # Centered 6th-order, seven-point stencil:
    # df/dx = ( -f(x-3h) + 9 f(x-2h) - 45 f(x-h) + 45 f(x+h) - 9 f(x+2h) + f(x+3h) ) / (60 h)
    elif method == "7-points 6-orders":
        df_dx = ( - f(x-3.0*h) + 9.0*f(x-2.0*h) - 45.0*f(x-h) + 45.0*f(x+h) - 9.0*f(x+2.0*h) + f(x+3.0*h) ) / ( 60.0*h )

    # Centered 8th-order, nine-point stencil:
    # df/dx = ( 3 f(x-4h) - 32 f(x-3h) + 168 f(x-2h) - 672 f(x-h) + 672 f(x+h) - 168 f(x+2h) + 32 f(x+3h) - 3 f(x+4h) ) / (840 h)
    elif method == "9-points 8-orders":
        df_dx = (     3.0*f(x-4.0*h) -  32.0*f(x-3.0*h) + 168.0*f(x-2.0*h) - 672.0*f(x-h)                      \
                  + 672.0*f(x+h)     - 168.0*f(x+2.0*h) +  32.0*f(x+3.0*h) -   3.0*f(x+4.0*h) ) / ( 840.0*h )

    return df_dx

############################################################################################



############################################################################################

## first_order_derivative_at_t0(f, t, i, method)
## Compute the first time derivative of a uniformly sampled discrete series f(t)
## at index i using centered finite-difference stencils.
##
## Inputs:
##  - f: array-like samples of the function evaluated at times t
##  - t: array-like time coordinates (assumed uniform spacing)
##  - i: integer index where the derivative is evaluated
##  - method: stencil identifier string; one of
##      "3-points 2-orders", "5-points 4-orders", "7-points 6-orders", "9-points 8-orders"

def first_order_derivative_at_t0( f, t, i, method ):
    
    dt = t[1] - t[0]

    # Centered 2nd-order difference:
    # df/dt = ( f[i+1] - f[i-1] ) / (2 dt)
    if method == "3-points 2-orders":
        df_dt = ( f[i+1] + f[i-1] ) / ( 2.0*dt )

    # Centered 4th-order, five-point stencil:
    # df/dt = ( f[i-2] - 8 f[i-1] + 8 f[i+1] - f[i+2] ) / (12 dt)
    elif method == "5-points 4-orders":
        df_dt = ( f[i-2] - 8.0*f[i-1] + 8.0*f[i+1] - f[i+2] ) / ( 12.0*dt )

    # Centered 6th-order, seven-point stencil:
    # df/dt = ( -f[i-3] + 9 f[i-2] - 45 f[i-1] + 45 f[i+1] - 9 f[i+2] + f[i+3] ) / (60 dt)
    elif method == "7-points 6-orders":
        df_dt = ( - f[i-3] + 9.0*f[i-2] - 45.0*f[i-1] + 45.0*f[i+1] - 9.0*f[i+2] + f[i+3] ) / ( 60.0*dt )

    # Centered 8th-order, nine-point stencil:
    # df/dt = ( 3 f[i-4] - 32 f[i-3] + 168 f[i-2] - 672 f[i-1] + 672 f[i+1] - 168 f[i+2] + 32 f[i+3] - 3 f[i+4] ) / (840 dt)
    elif method == "9-points 8-orders":
        df_dt = (     3.0*f[i-4] -  32.0*f[i-3] + 168.0*f[i-2] - 672.0*f[i-1]                  \
                  + 672.0*f[i+1] - 168.0*f[i+2] +  32.0*f[i+3] -   3.0*f[i+4] ) / ( 840.0*dt )

    return df_dt

############################################################################################