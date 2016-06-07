# alglib.lsfit
Experiments with fitting data using the Levenberg-Marquardt Non-linear Least Squares Solver and returning Confidence Intervals

## Purpose of this Project
This is a prototype of the curve fitting done for SPR-based Kinetic analysis, before taking it into the actual product.
We began using AlgLib, and have used both the minlm and lsfit, where lsfit gives standard error estimates back, which is useful to us.

## Difference between minlm and lsfit

The lsFit package expects the ndimensional_pfunc to be called per data-point to produce model values

 - This is unlike the minlm package, where the error function returns the entire residual vector
	
With five curves of 241 points each (1205 points total), our ndimemsional_pfunc is told to calculate a model value for a set of c_params and a given x[] vector. The naïve approach is to store only the per-curve time value in the x vector, but we also need to know the concentration of the curve, and the time at which dissociation begins.

***Store extra metadata in X***
Instead of storing only time in X, we can add concentration, and tDissocStart
