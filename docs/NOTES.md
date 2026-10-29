# Wave equation on AdS4 (r = L tan x)

Radial massless wave equation on global AdS4, evolved with RK4 in time and a
4th order centred finite difference in space.

Variables: `Pi = phi_t`, `Phi = phi_x`, with `x` in `[0, pi/2]`.

```
Phi_t = Pi_x
Pi_t  = (1/tan^2 x) d(tan^2 x Phi)/dx
```

Boundary behaviour: `Pi` even / `Phi` odd at `x = 0`; near `x = pi/2`,
`phi = f(t) rho^3 + O(rho^5)` with `rho = pi/2 - x`, so `Pi` is odd and `Phi` is even.

## Files

The library is in `src/`, the animation program in `apps/`, and the measurement programs (convergence and
boundary checks) in `tests/`.

| File | Purpose |
|---|---|
| `src/derivative.cpp/.hpp` | Stencil `derP` / `der`, AdS operators `adsDPi`, `adsDivergence` |
| `src/evolution.cpp/.hpp` | One RK4 step, `evolve` (optional whole-vector operators for Pi and Phi) |
| `src/structure.cpp/.hpp` | `ScalarField` struct and initial profile `gaussTan` |
| `src/einstein.cpp/.hpp` | Einstein-scalar-AdS4: `solveGeometry` (A, delta) and `einsteinRHS` |
| `src/progress.hpp` | Terminal progress bar used by `test_convergence_point` and `animation` |
| `apps/animation.cpp` | Dumps Pi and Phi frames (writes `animation/frames.dat`) |
| `tests/test_convergence_point.cpp` | Convergence ratio at one point, evolved to a common final time |
| `tests/test_convergence_linear.cpp` | L2 self-convergence versus time (writes `convergence/`) |
| `tests/test_convergence_einstein.cpp` | Constraint residuals and self-convergence of the coupled system |
| `tests/test_boundary_flat.cpp` | Comparison with the exact solution of the **flat** 1D wave equation (see note) |

## Command-line arguments

All arguments are optional and **positional**: to set a later one you must also give every one
before it (`./build/animation 500 3 10 einstein 20` is how you set the mode and amplitude). An argument that is
not a number is read as 0 for the numeric ones (`atof`/`atoi`), so check the first line the program prints.
Build lines are in the sections below.

### `animation` (`animation.cpp`)

`./build/animation [n] [T] [stepsPerFrame] [mode] [amplitude] [maxPoints] [maxFrames] [framesPerCycle] [cfl]`

| # | Argument | Meaning | Default | Valid values and notes |
|---|---|---|---|---|
| 1 | `n` | Number of grid intervals (`n + 1` points). Sets `dr = (pi/2)/n` and the time step `dt = dr/2`. | `500` | Integer `>= 4`. |
| 2 | `T` | Final time of the run. | `2*pi` (6.283185307) | `T > 0`. Frames are written at multiples of `stepsPerFrame` steps, so the last frame can end a little before `T`. |
| 3 | `stepsPerFrame` | RK4 steps between two written frames. The time between frames is `stepsPerFrame * pi / (4 n)`. | `20` | Integer `>= 1`. It can only ask for **sparser** output than `framesPerCycle` (argument 8) does: the larger of the two is used, and a `note:` line says when it is raised. With `framesPerCycle = 0` it is used as given. |
| 4 | `mode` | Which equation is evolved. | `ads` | `ads`: fixed AdS background (`A = 1`, `delta = 0`), using the same symmetrized right-hand side as `einstein` (see "Symmetrized right-hand side" below). `ads-expanded`: the same equation through the pair of operators of `derivative.cpp`; it agrees with `ads` to 1e-14 (same discretization, different arithmetic). `einstein`: full backreaction, stops when `min A < 1/n`. |
| 5 | `amplitude` | Initial data `Pi(0, x) = (amplitude / pi) * exp(-tan^2 x / sigma^2)`, `sigma = 1/16`, `Phi(0, x) = 0`: the peak of `Pi` is `amplitude / pi`. | `1` | Any real number. In `einstein` mode (`n = 500`, `T = 2 pi`): 30 stays regular (`min A` = 0.55, `delta` down to -1.3), 60 collapses at t = 0.009, and data with `min A < 1/n` is rejected at t = 0 from an amplitude of about 735 (the limit scales with `n`). For `ads` and `ads-expanded` the field is linear, so the amplitude only scales the picture. |
| 6 | `maxPoints` | Maximum grid points written per frame (every k-th point, always including `x = 0` and `x = pi/2`). | `256` | `0` or negative gives 256; `1` is raised to 2; the limit is 4096. A grid with `n + 1 <= maxPoints` is written in full. Near a collapse the interesting structure is a few grid cells wide at the origin, so for those runs give `maxPoints >= n + 1`; with the default 256 and `n = 1000`, only every 4th point is written and a dip in `A` at grid index 7 is not in the file. |
| 7 | `maxFrames` | Ceiling on the number of frames; used only when given. If `T / (stepsPerFrame * dt)` would give more, `stepsPerFrame` is raised so that the frames stay evenly spaced and cover the whole run. | `0` (none) | `0` = no limit; `1` is raised to 2; negative values are rejected. It can push the density below `framesPerCycle`; the program then prints a `warning:` with the actual frames per cycle. The count can be a little under the limit, because `stepsPerFrame` is a whole number of steps. |
| 8 | `framesPerCycle` | Frames per time interval `pi`. The dynamics repeats with period `pi`: a pulse returns to the origin every `pi` with its sign flipped (the field itself has period `2 pi`). The frame spacing is `pi / framesPerCycle`, i.e. `stepsPerFrame = 4 n / framesPerCycle`. | `100` | `0` disables the rule. 100 puts consecutive frames about 0.3 pulse widths apart for sigma = 1/16 (the pulse has FWHM about 0.1 and moves at speed 1), which looks smooth; 50 is borderline (0.6 pulse widths). With backreaction the period shifts slightly, so `pi` is the nominal cycle. |
| 9 | `cfl` | Time step in units of the grid spacing: `dt = cfl * dr`. Frame spacing in time does not change. | `0.5` | `cfl > 0`. With backreaction the characteristic speed is `A e^{-delta}`, which exceeds 1 once `delta < 0` (here `delta` reaches about -5, a speed near 150 times `A`), so near a collapse 0.5 is too large and the run ends in NaN with `min A` still around 0.2 to 0.8. In the amplitude-20 case at `n = 1000`, `cfl = 0.025` ends on the real criterion (`min A = -0.004 < 1/n`, at `t = 21.9535`); `0.5`, `0.25`, `0.1` and `0.05` all stop earlier with NaN. The cost grows as `1/cfl`. |

Examples:

```bash
./build/animation                              # linear AdS, n=500, T=2 pi, 20 steps per frame
./build/animation 500 6.2832 20 einstein 1.0   # with backreaction
./build/animation 2000 20 4 ads-expanded .01   # symmetrized expanded form, long run
./build/animation 5000 1 50 ads 1 128          # fine grid, at most 128 points per frame
./build/animation 500 62.8318 20 ads-expanded 1   # 20 cycles, 100 frames per cycle = 2001 frames
./build/animation 500 62.8318 20 ads 1 256 600  # ceiling of 600 frames (about 30 per cycle: a warning is printed)
./build/animation 500 20 20 ads 1 256 0 50   # 50 frames per cycle
```

Output: `animation/frames.dat` (columns `x Pi Phi A delta`), `animation/times.dat`, then
`gnuplot scripts/animation.gp` makes `animation/wave.gif`. `scripts/*.gp` takes no arguments; its settings
(`cap`, `pad`, the `loop` value and the frame delay) are in the first lines of the file.

The y-axis of the top panel (Pi and Phi) is rescaled for every frame to that frame's own range, from the minimum to the
maximum of Pi and Phi plus a margin of `pad = 0.08` of the range on both sides, so the shape is always drawn as large as
possible. (I also tried snapping the axis to the power of ten above the data, `-10^k .. 10^k`. It keeps the axis still, but
the largest value then used a median of only 26% of the half-height, so it was reverted.) Because the scale changes
from frame to frame, the plot has dense reference lines: vertical lines every `pi/8` (labelled), horizontal major lines at the
y tics with four minor lines between them, and a darker line at `y = 0`. They are set by `xtics`, `mytics`, the two `set style line`
lines and `set grid` in `scripts/animation.gp`. The bottom panel (A, delta) keeps one range for the whole run, so changes of the background stay
comparable between frames. A frame with no valid values keeps the previous range.

### `test_convergence_point` (`tests/test_convergence_point.cpp`)

`./build/tests/test_convergence_point [mode] [amplitude]`

| # | Argument | Meaning | Default | Valid values and notes |
|---|---|---|---|---|
| 1 | `mode` | Which equation is evolved. | `ads` | `ads` or `einstein` (as above). `ads-expanded` is **not** accepted here. |
| 2 | `amplitude` | Same as in `animation`. | `1` | Any real number. |

The resolutions `n = 250, 500, 1000`, the final time `T = pi` and `sigma = 1/16` are fixed in the source.

### `test_convergence_einstein` (`test_convergence_einstein.cpp`)

`./build/tests/test_convergence_einstein [amplitude] [sigma]`

| # | Argument | Meaning | Default | Valid values and notes |
|---|---|---|---|---|
| 1 | `amplitude` | Multiplies the initial `Pi(0, x) = amplitude * exp(-tan^2 x / sigma^2)`. | `0.5` | Real number. The run does not stop on collapse; it reports `min A`. |
| 2 | `sigma` | Width parameter of the initial profile. | `1.0` | `sigma > 0`. Small values give a narrow pulse. |

### `test_boundary_flat` (`test_boundary_flat.cpp`)

`./build/tests/test_boundary_flat [width] [r0]`

| # | Argument | Meaning | Default | Valid values and notes |
|---|---|---|---|---|
| 1 | `width` | Width `w` of the pulse `exp(-(r - r0)^2 / w^2)`. | `0.1` | `w <= 0.15` keeps the data smooth under the reflections; `0.2` and `0.3` lose order (see Lessons learnt). |
| 2 | `r0` | Centre of the pulse. | `pi/4` (0.7854) | In `[0, pi/2]`. |

### `test_convergence_linear` (`test_convergence_linear.cpp`)

Takes no arguments: `n = 250, 500, 1000, 2000` and the sample times `t = k*pi/20` are fixed in the source.

### Exit codes (`test_convergence_point` and `animation`)

| Code | Meaning |
|---|---|
| `0` | The run finished. |
| `1` | Bad arguments (unknown `mode`, or `n < 4`, `T <= 0`, `stepsPerFrame < 1` in `animation`). Nothing is written. |
| `3` | Stopped early: collapse in `einstein` mode (`min A < 1/n` or `A` not finite, also if the initial data already violates it), or a blow-up guard in `ads-expanded` mode (`abs(Pi)` or `abs(Phi)` above 1e3, which a stable run never reaches). The files written so far end at the last finite state. |

## Build and run

All commands are run from the repository root. Requires `g++`, `make` and `gnuplot` (with the `gif` and `pngcairo` terminals). `make` builds everything into `build/`; `make OPENMP=1` adds OpenMP.

### Convergence ratio (`test_convergence_point.cpp`)

```bash
make build/tests/test_convergence_point
./build/tests/test_convergence_point [ads|einstein] [amplitude]        # defaults: ads 1
```

`ads` keeps A = 1, delta = 0 (the linear wave equation on the fixed AdS background, the
previous behaviour). `einstein` solves the constraints for A and delta at every RK4 stage
(see the Einstein-scalar section below) and also prints the minimum of `A`.
The initial profile is `amplitude * gaussTan(x, 1/16)`.

Evolves n = 250, 500, 1000 to the same final time T = pi (`4n` steps, `dt = dr/2`),
prints the value of `Pi` at `x = pi/20` for each n and `ratio = log2(...)`, the convergence order.
It also writes `is<n>.txt` and `iter<n>.txt` (columns: `x Pi dPi Phi dPhi`).

At T = pi the pulse has made a full round trip and the errors are close to roundoff,
so the ratio is not a clean order. For a clean measurement choose a time such as
19 pi / 20 (`19*n/5` steps).

### Self-convergence versus time (`test_convergence_linear.cpp`)

```bash
make build/tests/test_convergence_linear
./build/tests/test_convergence_linear
gnuplot scripts/convergence.gp
```

Runs n = 250, 500, 1000, 2000 and samples at `t = k*pi/20`, `k = 1..20`.
Output in `convergence/`:

- `error.dat`: `t  dr  L2 error`, one block per resolution pair (`index 0/1/2` in gnuplot)
- `order.dat`: `t  order(250/500/1000)  order(500/1000/2000)`
- `convergence.png`: error versus time (log y) and measured order, written by `scripts/*.gp`

Takes about 12 seconds. For log-log plots against the grid spacing, use column 2 (`dr`)
of `error.dat` at fixed `t`:

```gnuplot
set logscale xy
plot 'convergence/error.dat' using 2:3
```

### Animation of Pi and Phi (`animation.cpp`)

```bash
make build/animation
./build/animation              # defaults: n=500, T=2*pi, 20 steps per frame
gnuplot scripts/animation.gp
```

Output: `animation/frames.dat` and `animation/wave.gif` (201 frames, t = 0 to 2 pi). The GIF plays
once and stays on the last frame; change `loop -1` to `loop 0` in the `set terminal` line of
`scripts/animation.gp` to repeat forever.

All arguments are described in the table under "Command-line arguments" above. In short:
`./build/animation [n] [T] [stepsPerFrame] [mode] [amplitude] [maxPoints] [maxFrames] [framesPerCycle] [cfl]`.

The extra mode `ads-expanded` evolves the same linear equation written in a symmetrized form
(`adsDivergenceExpanded` + `adsDPiExpanded`, see "Symmetrized form" below). It is stable and gives the same
results as `ads` (difference about 1e-9, the output precision, for `n = 2000`, `T = 20`).

`maxPoints` caps the grid points written per frame (every k-th point, plus the last one), with a
hard limit of 256. With `n = 5000` and about 1000 frames this gives a `frames.dat` of about 8 MB
instead of 300 MB. The number of frames is set by `framesPerCycle` (default 100 per interval `pi`), not by `stepsPerFrame`: a 20-cycle run is
2001 frames and 18 MB at any `n` (the old defaults gave 400 per cycle at `n = 2000`, and `maxFrames = 1000` silently thinned a
long run to 49 per cycle, which is about 0.6 pulse widths between frames and looks jumpy). Rendering the GIF takes
about 1 minute for 640 frames and 7 minutes for 2000 (the time grows faster than linearly), so use `maxFrames` or a smaller
`framesPerCycle` for long runs; the program prints the number of frames, the spacing and the frames per cycle, and a `note:`
above 3000 frames.

**Progress bar.** `test_convergence_point` and `animation` draw a bar on stderr with the simulated time, for example
`[##########------]  45.3%  t = 1.425 / 3.142`. It only appears when stderr is a terminal, so
redirecting output (`2> log.txt`) or piping keeps it clean. Implemented in `progress.hpp`.

**Collapse detection (`einstein` mode, `test_convergence_point` and `animation`).** After every step the minimum of
`A` is checked. If `min A < 1/n`, or `A` is not finite, the run stops with exit code 3 and prints
the time, step and position (`x`, grid index) of the minimum; `animation` still writes the last
state as the final frame, and `test_convergence_point` writes `iter<n>.txt` and skips the convergence ratio. Initial
data that already violates the condition is rejected before evolving. A NaN usually means
the time step `dt = dr/2` is too large for such strong fields (the characteristic speed is
`A e^{-delta}`), so near a collapse treat the reported time as approximate: a smaller `dt`
reaches `A -> 0` slightly later. With sigma = 1/16 and the current amplitude scale (peak of `Pi` = amplitude / pi), amplitude 30 stays regular (min A = 0.55 over `T = 2 pi`)
and amplitude 60 collapses at t of about 0.009.

### Symmetrized form of the wave equation (`adsDivergenceExpanded`, `adsDPiExpanded`)

The plain expanded equation `Pi_t = Phi_x + 2 Phi/(sin x cos x)`, `Phi_t = Pi_x` is unstable with a
centred stencil, for every time step: the growth rate is about `0.39/dr` near `x = pi/2` and `0.1/dr`
near `x = 0` (measured from the eigenvalues of the semi-discrete operator), so a finer grid blows up
sooner. Rewriting `Phi/(sin x cos x)` with l'Hopital near the ends (`Phi'/cos 2x`, or isolating `cos x`
to get `Phi'/cos^2 x`) lowers the rate near the origin but never removes it, and near `pi/2` it makes it worse.

*Why.* With `w = tan^2 x` the continuous system conserves `E = int w (Pi^2 + Phi^2) dx` because
`Pi (w Phi)_x + w Phi Pi_x = (w Pi Phi)_x`. The expanded form uses the product rule
`(w Phi)_x = w Phi_x + w_x Phi`, which a stencil does not satisfy. For smooth data the defect is tiny
(1e-5 to 1e-9 relative), but for grid-scale oscillations it is of order 1 (0.6 to 2.2), and those are the
modes that grow. The conservative form (`adsDivergence`) never splits the product, so it has no defect.

*Cure.* Absorb `tan x` into the fields, `Theta = tan x Pi`, `Psi = tan x Phi`:

```
Theta_t = Psi_x + b Psi,     Psi_t = Theta_x - b Theta,     b = 1/(sin x cos x)
```

The `b` terms have opposite signs and cancel pointwise in `Theta Psi_t + Psi Theta_t`, so no product rule is
needed. With the reflection rules used in `derP` the two stencils satisfy summation by parts
(`H D_theta + D_psi^T H = 0`, `H = diag(1/2, 1, ..., 1, 1/2)`), and `E = sum H (Theta^2 + Psi^2)` is conserved
exactly by the semi-discrete system: its eigenvalues are purely imaginary (largest real part 1e-14 for
`n = 50..400`, `|lambda| dr = 1.44`, well inside the RK4 limit for `dt = dr/2`). Going back to `Pi`, `Phi` is a
diagonal change of variables, which does not change the eigenvalues:

```
Pi_t  = D(tan x Phi)/tan x + Phi/(sin x cos x)
Phi_t = D(tan x Pi )/tan x - Pi /(sin x cos x)
```

*Implementation.* `Theta` and `Psi` are never stored. `weightedStencil` (in `derivative.cpp`) evaluates
`tan x_j * f_j` on the fly at the stencil points, which is the ordinary stencil with each neighbour
weighted by `tan(x_{i+k}) / tan(x_i)`; ghost points use the parity of the product (`tan x` is odd at both
ends, so `tan x * Phi` is even at 0 and odd at pi/2, and `tan x * Pi` the opposite).

*Accuracy next to the origin.* The division by `tan x_i ~ i dr` multiplies the stencil error by `1/(i dr)`.
With the 5-point stencil of `derP` the operator was therefore only `O(dr^3)` at the first grid points (2.4e-7
at `i = 1` for `n = 250`, orders 3.0), while it is order 4 everywhere else. `weightedStencil` now uses a
7-point, 6th order centred stencil (`-1/60, 3/20, -3/4, 0, 3/4, -3/20, 1/60`): the operator error is
`O(dr^5)` at the first points and `O(dr^6)` elsewhere (1e-10 at `i = 1` for `n = 250`, down to the roundoff
floor of about 1e-13 on finer grids), so the pair is at least 4th order at every grid point, and the full
evolution is limited by RK4 in time. Energy conservation is kept because it only needs an antisymmetric
centred stencil and the reflection rules at the ends, for any width (the eigenvalues stay purely imaginary:
largest real part 1e-13 to 8e-12, `|lambda| dr` goes from 1.44 to 1.62). A wide stencil used only near the
origin would break the symmetry, so the whole operator uses 7 points.

Measured with the 7-point stencil: L2 order 4.1 to 4.6 for sigma = 1 and 4.4 to 4.5 for sigma = 1/16; the max-norm
order for sigma = 1/16 is a steady 3.8 to 4.0 (with the 5-point stencil it was 2.4 to 4.8 depending on the time).
At the first refocus on the origin (`t = pi`, exact `Pi(0) = -1`, sigma = 1/16) the error at `n = 250` is
2.9e-4, against 4.5e-3 with the 5-point form and 3.9e-3 for the conservative form.

*Cost.* The weights are computed on the fly, which calls `tan` several times per point: on
`./build/animation 2000 20 4 ads-expanded .01` it takes 84 s (74 s with the 5-point stencil, 25 s for the
conservative `ads` form). Getting `tan(x_i + k dr)` from `tan x_i` with the angle-addition formula halves the
time of the 5-point version but changes the results by 3e-11 relative near `pi/2` (cancellation in
`1 - tan x * tan dr`), which breaks the exact skew-symmetry, so it was not adopted.

Both equations must be used together, since the cancellation needs the pair. Parities: `Theta` is odd at
`x = 0` and even at `pi/2`, `Psi` the opposite; both vanish at the two end nodes.

Stability: stable to `T = 400` with `n = 500`, sigma = 1/16, and your `2000 20 4 ads-expanded .01` run completes.

### Symmetrized right-hand side of `ads` and `einstein` (`einstein.cpp`, `computeRHS`)

`ads` and `einstein` used to evaluate `Pi_t = (tan^2 x m Phi)_x / tan^2 x` (with `m = A e^{-delta}`, `m = 1` for
`ads`) in conservative form. That operator has an unstable mode at the first grid points: the origin node
`Pi_0` has energy weight `tan^2(0) = 0` but still drives `Phi_1` and `Phi_2` through the stencil, and its own
equation (`3 Phi_x`) is not the conjugate of those couplings, so the discrete operator cannot be skew in
the energy weight. Measured on the matrix of the operator: largest `Re(lambda)` = 0.12, 0.20, 0.21, 0.21 for
`n = 50, 100, 200, 400` (it does not shrink with `n`; the number of unstable modes grows with `n`), with the
eigenvector peaking at `i = 1`. The `i = 1` special case I had added to improve the accuracy there doubled the rate
(0.10 without it). In long runs (`n = 500`, sigma = 1/16) roundoff grew into a blow-up around `t = 150-300`; in
`einstein` mode this looks like a collapse (`min A` falling to 0.80 at `t = 250` from amplitude 0.01) before the run
turns into NaN, so it can be mistaken for physics. It is invisible for `T` of a few tens.

The cure is the one used for the linear case, generalized to a weight `m(x)`. With `a = tan x`, `p = sqrt(m)`:

```
Theta = p a Pi ,   Psi = p a Phi ,   b = 1/(sin x cos x)
Theta_t = S Psi + m b Psi ,   Psi_t = S Theta - m b Theta ,   S = (m d/dx + d/dx m) / 2
```

`S` is skew-adjoint and the `b` terms cancel pointwise, so `E = int (Theta^2 + Psi^2)` is conserved. Going
back to the original variables is a diagonal change: `Pi_t = S(p a Phi)/(p a) + m b Phi`,
`Phi_t = S(p a Pi)/(p a) - m b Pi`, which is exactly `Phi_t = (m Pi)_x`, `Pi_t = (tan^2 x m Phi)_x / tan^2 x`. `Theta` and `Psi`
are not stored: `p a f` is evaluated on the fly at the 7 stencil points (6th order, as in the previous
section), with ghost points following the parity of the product. The origin node is decoupled
(`Theta_0 = 0`), which removes the unbalanced coupling.

Measured: no unstable eigenvalue (largest `Re(lambda)` 1e-13 to 7e-12 for `n = 50..400`, for `m = 1`, and for a
bump, a steep dip and a profile with an odd `rho^3` part at `pi/2`; `|lambda| dr = 1.62`); `ads` runs to `T = 400`
with the inner value flat at 5e-6 and `einstein` (amplitude 0.01) runs to `T = 400` with `min A = 1`, where the old
code blew up near `t = 200` and `t = 300`. The right-hand side is more accurate than the conservative one: for
`Pi_t` at `x < 0.1` the error is order 4 (9.4e-10, 7.7e-11, 5.3e-12, 3.5e-13 for `n = 100..800`) where the old one
was order 2 (7.9e-6 to 1.2e-7), and 10 to 100 times smaller elsewhere. `Phi_t` is order 3 at `x < 0.1`, the same as
before, and set by the quadrature of `A` there. The full coupled system converges at order 4.1 to 4.4 in L2
(`test_convergence_einstein`, amplitude 0.5), the collapse thresholds are unchanged (amplitude 30 regular, 60 collapses
at the same step), and `min A` agrees with the old code to all printed digits.

`adsDivergence`/`adsDPi` (the conservative pair) are still in `derivative.cpp` and are still used by `test_convergence_linear.cpp`,
which compares resolutions only up to `T = pi`; do not use them for long runs.

### Einstein-scalar system with backreaction (`einstein.cpp`)

The scalar field changes the metric, `ds^2 = L^2/cos^2 x (-A e^{-2 delta} dt^2 + dx^2/A + sin^2 x dOmega^2)`,
with `Phi = phi_x`, `Pi = A^{-1} e^{delta} phi_t` and

```
Phi_t  = (A e^{-delta} Pi)'
Pi_t   = (1/tan^2 x) (tan^2 x A e^{-delta} Phi)'
delta' = -sin x cos x (Phi^2 + Pi^2),                                   delta(0) = 0
A'     = (1+2 sin^2 x)/(sin x cos x) (1-A) - sin x cos x A (Phi^2+Pi^2),  A(0) = 1
```

A and delta are not stepped in time. At every RK4 stage they are recomputed from the
current Pi and Phi in closed integral form (`solveGeometry`):

```
delta(x) = -int_0^x sin y cos y S dy
1 - A(x) = e^{delta(x)} (cos^3 x / sin x) int_0^x e^{-delta(y)} tan^2 y S dy,    S = Phi^2 + Pi^2
```

using a 4th order cumulative quadrature. The time step uses `evolveRK4(sf, &einsteinRHS)`.
`adsRHS` is the same code with A = 1 and delta = 0 frozen, so one right-hand side serves both cases.

```bash
make build/tests/test_convergence_einstein
./build/tests/test_convergence_einstein [amplitude] [sigma]     # defaults 0.5 and 1.0
```

Prints the residuals of the two constraint equations (delta' order 4; A' order 3 because of
the first grid points next to the origin) and the self-convergence order of Pi versus time
(about 4.0), plus the minimum of `A` over the run (a value close to 0 signals collapse).
Writes `convergence/einstein_order.dat`.

### Exact-solution boundary test (`test_boundary_flat.cpp`)

```bash
make build/tests/test_boundary_flat
./build/tests/test_boundary_flat            # pulse width 0.1, centre pi/4
./build/tests/test_boundary_flat 0.15 0.7854  # width, centre
```

Writes `convergence/boundary.dat` (`t  n  max_err_all  max_err_right`) and prints the order at
each time for the whole domain and for `r > 1.2`.

Note: this program tests the **flat** 1D wave equation `Pi_t = Phi_x`, `Phi_t = Pi_x` with
`Pi` even and `Phi` odd at both ends, using the default `der` operators (no AdS operators).
It does not test the AdS equation.

## Notes

- Time step is `dt = dr/2`; the characteristic speed is 1 in `x`, so this is stable.
- Convergence tests must compare resolutions at the same final time. Comparing a single
  step with `dt` proportional to `dr` measures how the solution changes with time, not the scheme order.
- Initial data must be smooth under the boundary reflections. A pulse with a nonzero
  `Pi'(0)` (or a non-negligible tail at `x = pi/2`) produces a kink after reflection and
  reduces the measured order.

## Lessons learnt

What went wrong, or nearly did, while building this code, and what to do about it.

### Testing convergence

- **Compare resolutions at the same final time.** The first `test_convergence_point.cpp` took one step per grid with
  `dt = dr/2`, so each resolution stopped at a different time and the printed ratio was exactly
  2.0000: it measured how the solution changes with `t` (the `t^2/2 * Pi_xx` term), not the scheme
  order. Evolving every resolution to the same `T` gave 4.00.
- **Check the right-hand side before judging the scheme.** An earlier version had
  `dPidt = Pi` instead of `Phi_x`, so the spatial derivative never entered the test at all.
- **Know where the error floor is.** At `T = pi` (full round trip) the errors against the exact
  answer are about 1e-8 to 1e-14, near double precision, and the apparent order (5 to 7) means
  nothing. A single point value can also change sign between resolutions. Prefer an L2 norm and a
  time where the errors are well above roundoff.
- **Use an exact solution when there is one.** For the flat wave equation the method of images gives
  one (`test_boundary_flat.cpp`), which separates scheme errors from data errors.
- **The initial data limits the order.** If it is not smooth under the boundary reflection the order
  drops: a nonzero `Pi'(0)` gave about 1.2, and a wide pulse with a tail at the ends gave 0.8.
  Check the data before suspecting the boundary condition.

### Equations and boundaries

- **The metric decides the equation.** `r = L tan x` in flat space gives
  `Pi_t = (cos^4 x / L^2)(Phi_x + 2 cot x Phi)`. In AdS4 the equation is
  `Pi_t = (1/tan^2 x) (tan^2 x Phi)'`. My first derivation assumed flat space and was wrong for this problem.
- **Conservative form is not enough; make the operator skew.** Writing the equation as
  `Phi_x + 2 Phi/(sin x cos x)` divides by quantities that vanish at both ends and is unstable; applying the
  stencil to `tan^2 x * Phi` and dividing afterwards (conservative form) cures that, and looks stable for runs of
  a few tens of time units. But its discrete operator is not skew in the energy weight (the origin node has weight 0
  and still couples to its neighbours), which leaves a growing mode at `i = 1` that shows up only after `t` of 100 to
  300. What is really needed is the symmetrized pair (`Theta = tan x Pi`, `Psi = tan x Phi`, with a weight `m` for the
  Einstein case). Check it with the eigenvalues of the operator, not with a run that happens to be short.
- **Take the boundary parities from the asymptotics.** With `phi = f rho^3 + O(rho^5)`, `Pi` is odd and
  `Phi` is even at `x = pi/2` (and the opposite at `x = 0`). Reflecting with the wrong parity gave
  first-order error near the outer boundary (order 1.0 from `t` of about 0.16 on). `derP` takes a
  separate parity for each end.
- **Fill singular limits explicitly.** At `x = 0`, `Phi/x -> Phi_x` gives `Pi_t = 3 Phi_x`. At `x = pi/2`,
  `tan^2 x * Phi` has a finite limit that is extrapolated from three interior points (even polynomial in `rho`).
- **Apply parity only to ghost points.** The old stencil flipped the sign of real points too (`<=` instead of `<`).
  It was harmless only because those values were exactly zero.

### Backreaction (A and delta)

- **Derive the sign from the metric.** The factor is `e^{-delta}`, not `e^{+delta}`, for the given `delta'`
  and `A'`. A wrong sign would have broken the constraints silently.
- **Solve the constraints in closed form.** `delta` is a plain integral and `A` follows from an
  integrating factor plus an integration by parts. Both integrands are smooth, including at `pi/2`, and
  computing `1 - A` directly avoids cancellation and the `1/(sin x cos x)` singularities.
- **Test an implementation against its own equations.** Insert the computed `A` and `delta` into the
  differential equations: the `delta'` residual is order 4, the `A'` residual is order 3 next to the origin
  (localized, and it does not show in the evolution, which stays at order 4).
- **Share one right-hand side.** The linear case is just `A = 1, delta = 0`, so `adsRHS` and `einsteinRHS`
  use the same code and the old behaviour was reproduced to about 1e-14.

- **A short run is not a stability test.** An instability with `Re(lambda) = 0.2` multiplies roundoff by about 1e16 only after
  `t` of roughly 180: three earlier "stable" verdicts were based on `T <= 20`. Compute the eigenvalues, and run
  at least one long evolution (here `T = 400`) before calling a scheme stable.

### Robustness near collapse

- **The step limit changes with the fields.** The characteristic speed is `A e^{-delta}`; for strong fields
  `dt = dr/2` is too large and the run goes to NaN a few steps before `A` actually reaches 0
  (sigma = 1/16, amplitude 60 when the argument was the peak of `Pi`, i.e. 188 in today's scale where the peak is `amplitude / pi`: `A` crosses zero near t = 0.0073 with a smaller `dt`). Treat the reported
  collapse time as approximate.
- **Stop on `min A < 1/n` and on non-finite `A`.** Check the initial data too (with the argument equal to the peak
  of `Pi`, as it was then: 30 is regular with min A = 0.69 at n = 2000, T = pi; 300 is already past the horizon at t = 0).
- **Never write non-finite frames.** On a collapse the output ends at the last finite state, with no
  duplicate and no NaN rows.
- **Validate arguments.** `stepsPerFrame = 0` crashed with a floating-point exception, and a tiny `n`
  read out of bounds.

- **Stability follows from structure, not from where you patch.** The expanded equation
  `Phi_x + 2 Phi/(sin x cos x)` blew up for every time step; special treatment of the quotient near the ends
  (l'Hopital, hybrids of the two forms) only moved or reduced the growth. The reason was the broken product
  rule for the stencil (defect of order 1 on grid-scale oscillations), found by computing the eigenvalues of
  the linear operator as a matrix and measuring the defect. The fix was to change the variables so that
  the singular terms cancel pointwise (`Theta = tan x Pi`, `Psi = tan x Phi`), which gives an exactly
  conserved discrete energy.
- **Eigenvalues of the semi-discrete operator are the fastest diagnostic.** Build the matrix from unit
  vectors and check `max Re(lambda)`: an instability that scales like `1/dr` is the spatial discretization,
  not the time step, and where the eigenvector peaks says which end is responsible.

### Output and plotting

- **Do not store more than can be seen.** 5001 points per frame for 1002 frames was a 303 MB `frames.dat`
  for a plot a few hundred pixels wide. Capping the points per frame (1024 first, now 256) and using 6 digits
  brings it to about 8 MB for the same run (measured: 956 frames, 8.1 MB).
- **Make plot scripts robust.** One frame with values around 1e300 made gnuplot's axis ranges infinite
  ("y axis range undefined"). `scripts/*.gp` now ignores NaN/Inf and anything above `cap`.
- **Do not hard-code what the program can write.** Frame times are written to `times.dat` instead of
  being recomputed in the plot script from `n` and `stepsPerFrame`.
- **Rendering time grows faster than the data.** The GIF render took 8 s for 5 MB, 20 s for 10 MB and
  59 s for 21 MB, because each frame is looked up with `index k`. A fit to those three runs gives
  about 4 to 5 minutes for 50 MB; that estimate was not measured.
- **GIF looping.** gnuplot's `animate` loops forever by default, `loop 0` also loops forever, `loop 1` plays
  twice, and `loop -1` writes no loop marker so it plays once. Verified by reading the file's bytes.

### Workflow

- **Read the file before describing it.** An early analysis was written after the file reads returned
  nothing, and had to be thrown away. Say what was and was not checked.
- **Rebuild after every source change.** A stale `animation` built before a fix produced the bad last
  frame that broke the plot; a build line without `einstein.cpp` fails to link. Keep the build commands in
  this README and `.vscode/tasks.json` in step with the source files.
- **Generated files can explode.** Look at the size and the tail of an output file (a 3.8 GB NaN
  `frames.dat` appeared once). Binaries and output directories should not share a name
  (`animation` the program versus `animation/` the folder).
- **Measure, then extrapolate.** Timings, orders and ratios in this file come from runs; where a number is
  a fit rather than a measurement it says so.
