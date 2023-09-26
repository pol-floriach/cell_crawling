# Phase field model for cell movement.
Simulations in julia to simulate slightly different cases of cell movement, taking deformation into consideration by using a [phase field model](https://en.wikipedia.org/wiki/Phase-field_model).

## Structure

`main.jl`: Contains principal program in which some parameters can be tweaked. Can be run as an executable or with: 
```
julia main.jl
```
`constants.jl`: Contains the model's constants and some struct definitions.
`init.jl`: Contains functions that initialize the arrays.
`numerical.jl`: Contains functions used in the model, mainly to compute the gradient and laplacian with a fdm, and periodical boundary conditions.
`phasefield.jl`: Contains the core of the simulations. Slightly different versions of the model are implemented in the following functions:

phase_field()
phase_field_rigged()
phase_field_rigged_2()


## Dependencies
[Plots.jl](https://docs.juliaplots.org/stable/) for the visualization of the phase field through a time-dependent heatmap in `.gif` format.
[ProgressBars.jl](https://docs.juliahub.com/ProgressBars/qxsPw/0.7.1/) if desired, to see the estimated time for `main.jl` to be run.
