# Shallow Water Equation: Numerical Approximations

A selection of SWE functions to model SWE with a flat sea-bed and no Coriolis force,
  in one- and two-dimensions.

## Forward-in-Time, Centred-in-Space

  Using Euler approximation for the partial derivative in time, central differences method for
  the partial derivate in space.

## [Lax-Friedrichs](https://en.wikipedia.org/wiki/Lax%E2%80%93Friedrichs_method)

  Very similar to FTCS method, but with average of neighbouring points used in the
  update step.

##
  [Lax-Wendroff](https://en.wikipedia.org/wiki/Lax%E2%80%93Wendroff_method#Richtmyer_method)

  Here, we use the Richtmyer two-step Lax-Wendroff method to avoid calculation of
  the Jacobian matrix.


## [Lax-Fridrichs Lax-Wendroff Compositie (LFLWk)](https://epubs.siam.org/doi/abs/10.1137/S0036142996310976)

  It is known that the Lax-Friedrichs method is excessively diffusive and smears out
  shocks, whilst the Lax-Wendroff method suffers from oscillations unless the data is
  smooth. This behavior is especially evident when looking at water flowing over a
  non-constant  topography (which we have not considered here). It is possible to use
  Lax-Friedrichs as a filter for Lax-Wendroff, as shown in Liska & Wendroff (1998).
  The first step of Lax-Wendroff uses Lax-Friedrichs as a predictor to obtain the
  fluxes centred in time.

