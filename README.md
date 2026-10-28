# Convergence of the discontinuous virtual element method

Self-contained demo of the convergence of the weighted SIPG discretization of
virtual element spaces, for

    -div(rho grad u) = f   in (0,1)^2,    u = 0 on the boundary,

on three families of polygonal meshes and different coefficients.

## Running it

From MATLAB, with this directory as the current folder:

    >> ex_convergence

It solves 75 problems (3 mesh families x 5 coefficients x 5 refinement levels)
and takes about twenty minutes. For a quick look, shorten the third entry of
`levels` and the `rhoCases` list at the top of the script.

    >> ex_solutions

draws the discrete solution on a coarse mesh of each family and the coefficients
as the method sees them. It takes about a minute.

To redraw the figure or rewrite the table after a change of labels or of layout,
without solving anything again, load the saved errors and call the function that
writes them:

    >> S = load('ex_convergence.mat');
    >> writeConvergence(S.hh,S.eL2,S.eH1,S.einf,S.rL2,S.rH1,S.rInf, ...
                        S.families,S.famName,S.rhoCases,S.rhoName,'./');

## What it produces

| file | content |
|------|---------|
| `ex_convergence.pdf`       | log-log plots of the two errors against `h` |
| `ex_convergence.mat`       | raw errors and rates, to replot without recomputing |
| `ex_solutions.pdf`         | discrete solution on each family, and a cut across the interface |
| `ex_coefficients.pdf`      | the coefficients, elementwise constant as they enter the method |

## The coefficients

With `eta = 10`, so that the contrast is exactly ten:

1. `rho = 1`
2. `rho = 1 + eta*x*y`,
3. `rho = 2 + cos(10 pi x) cos(10 pi y)`,
4. `rho = 1` for `x < 1/2` and `eta` otherwise, with `u = sin(pi x) sin(pi y)`
5. the same coefficient as in (4), with an exact solution having a kink at the interface

Cases 4 and 5 share the coefficient. In case 4 the normal derivative of the
exact solution vanishes on the interface, so the jump is barely active; case 5
is the sharper test, with `rho phi'` continuous and nonzero there.

## The meshes

Structured triangular, hexagonal and Voronoi, five levels each, all of them
fitted to the line `x = 1/2` so that `rho` is constant on each element, which is
what the analysis assumes.

## What the errors are

Both errors are measured through its elliptic projection `Pi u_h`, 
which is affine on each element:

    ||u - Pi u_h||_0 / ||u||_0            expected O(h^2)
    |u - Pi u_h|_{1,h} / |u|_1            expected O(h)

computed by splitting each polygon into triangles from its centroid and
integrating with a rule of degree four. The penalty is `gamma_h = 10`.

The directory `reference/` holds the output of a full run, to compare against.

## Requirements

MATLAB (with no toolboxes). The meshes are included.

## Files

    ex_convergence.m     the convergence study
    writeConvergence.m   its figure and tables, separated so that they can be
                         redrawn from the saved errors
    ex_solutions.m       the companion figures
    coefficient.m        the five coefficients, their data and exact solutions
    discVem_weights.m    assembly of the weighted SIPG virtual element system
    vem.m                elementwise virtual element projector and stiffness
    meshes/              the eighteen meshes used by the two scripts
    reference/           output of a full run, for comparison
