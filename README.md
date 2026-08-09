# Violet-Mango
A 2.5 FDTD Solver in a single file web app.

For technical support contact: numerius.engineering@gmail.com
## Violet Mango — How the Solver Works
A 2.5D Finite-Difference Time-Domain (FDTD) electromagnetic solver, running entirely client-side.
1. Overview
Violet Mango simulates electromagnetic current sources on a single 2D computational grid (the "grid plane") that sits at a chosen height in 3D space. Rather than restricting the physics to a single polarization, the solver tracks all six Cartesian field components (Ex, Ey, Ez, Hx, Hy, Hz) by assuming the geometry and fields are invariant along the axis perpendicular to the grid plane — the same assumption used for 2D waveguide cross-sections and cylindrical antenna problems. This is what "2.5D" refers to: a 2D grid producing a fully 3D field vector at every point.

2. Mode Decoupling (TMz / TEz)
Because the fields don't vary along z (∂/∂z = 0), Maxwell's curl equations split cleanly into two independent sets:
```
// TMz mode — driven by out-of-plane current Jz
fields: Ez, Hx, Hy

// TEz mode — driven by in-plane current Jx, Jy
fields: Hz, Ex, Ey
```
These two modes do not exchange energy for isotropic materials, so the solver advances them simultaneously but independently on the same Yee grid. A source's orientation vector determines which mode(s) it excites: a purely out-of-plane dipole (dz only) drives TMz exclusively, a purely in-plane dipole (dx, dy) drives TEz exclusively, and a tilted dipole drives both — which is also why the out-of-plane Poynting component Sz is only non-zero when both modes are simultaneously excited (it comes from the cross-product of one mode's E-field with the other mode's H-field).

3. Update Equations
The grid uses the standard Yee leapfrog scheme: H-fields are updated a half timestep ahead of E-fields, cell sizes are uniform (Δx = Δy = dx), and materials are approximated with the collocated-grid convention (one εr/μr/σ value per cell, applied to every field component stored at that cell). With C0 the speed of light, ε0 and μ0 the vacuum permittivity/permeability, and dt the timestep:
```
// Magnetic field updates (TMz)
Hx[i,j] -= (dt / (μ0·μr·dy)) · (Ez[i,j+1] - Ez[i,j])
Hy[i,j] += (dt / (μ0·μr·dx)) · (Ez[i+1,j] - Ez[i,j])

// Magnetic field update (TEz)
Hz[i,j] -= (dt/μ0·μr) · [ (Ey[i+1,j]-Ey[i,j])/dx - (Ex[i,j+1]-Ex[i,j])/dy ]

// Electric field update (TMz), with conductive loss + PEC clamp
Ez[i,j] = Ca·Ez[i,j] + Cb·[ (Hy[i,j]-Hy[i-1,j])/dx - (Hx[i,j]-Hx[i,j-1])/dy - Jz[i,j] ]

// Electric field updates (TEz)
Ex[i,j] = Ca·Ex[i,j] + Cb·[ (Hz[i,j]-Hz[i,j-1])/dy - Jx[i,j] ]
Ey[i,j] = Ca·Ey[i,j] + Cb·[ -(Hz[i,j]-Hz[i-1,j])/dx - Jy[i,j] ]

// Loss coefficients, from conductivity σ and permittivity ε = ε0·εr
Ca = (1 - σ·dt / 2ε) / (1 + σ·dt / 2ε)
Cb = (dt/ε) / (1 + σ·dt / 2ε)
```

Cells flagged as a Perfect Electric Conductor (PEC) have every E-field component clamped to zero after each update, which is the standard "staircased PEC" approximation used throughout FDTD.

Time step & stability
The timestep is chosen from the 2D Courant–Friedrichs–Lewy (CFL) stability limit, scaled by the user-configurable CFL factor (default 0.98, must be < 1 for stability):
```
dt = CFL_factor · dx / (c · √2)
```
Absorbing boundary (PML)
The domain edges use a graded-conductivity absorbing layer: a polynomial-profile loss term σ(depth) is added to the update coefficients within the PML region (matched between electric and magnetic loss so the layer is approximately impedance-matched to free space), ramping from 0 at the PML's inner edge to a maximum value at the outer domain boundary. This is a simplified, single-field approximation of a full Convolutional PML — effective at absorbing outgoing waves with modest reflection, though not as reflection-free as a split-field CPML implementation.

4. Excitation & Steady State
All sources are single-frequency, continuous-wave (CW) sinusoids. The solver time-steps cycle by cycle, tracking the peak field amplitude at a monitor point each cycle. It stops early once the cycle-over-cycle amplitude change falls below the convergence threshold (if enabled), or after the configured maximum number of cycles. Once stopped, it runs one additional full period, during which it:

Records field snapshots (for the time-history playback / animation), and
Accumulates a running discrete Fourier transform (DFT) at the excitation frequency, to extract the steady-state complex phasor fields (Ex, Ey, Ez, Hx, Hy, Hz) used for all "steady-state phasor" visualizations and the far-field transform.
5. Near-to-Far-Field Transform
The far-field radiation pattern is computed from the phasor fields using the equivalence principle on a rectangular Huygens surface placed just inside the PML. Equivalent electric and magnetic surface currents are derived from the tangential fields on that contour:
```
Js = n̂ × H        (equivalent electric surface current)
Ms = -n̂ × E       (equivalent magnetic surface current)
```
These currents are integrated around the contour with the appropriate far-field phase factor exp(jk·r̂·r′) for each observation angle φ, separately for the TMz pattern (from Ez, Hx, Hy) and the TEz pattern (from Hz, Ex, Ey), producing the polar radiation pattern shown in the Visualize tab.

Approximate Directivity
Directivity is estimated by treating the far-field magnitude squared as proportional to radiation intensity, U(φ) ∝ pattern(φ)², and normalizing by the total power integrated around the full circle — the 2D/cylindrical analogue of the standard 4π-steradian directivity formula (using 2π radians in place of 4π steradians, consistent with this solver's z-invariant geometry):
```
D(φ) = 2π · U(φ) / ∫ U(φ) dφ
```
It's labeled "approximate" because it isn't impedance-calibrated to an absolute radiated-power reference — it's a relative measure derived purely from the shape of the computed pattern.

6. Materials
Material regions are rasterized onto the grid from geometric primitives, split into two families:

Surface primitives (open curves, w/ thickness)	Closed-region primitives (filled area)
Line, Circular Arc, Elliptical Arc, Parabolic Arc, Spline, ƒ(x)	Triangle, Circle, Ellipse
Each material carries relative permittivity εr, relative permeability μr, conductivity σ, and an optional PEC flag (which overrides εr/μr/σ with a perfect-conductor boundary condition). The elliptical arc uses the ellipse's native parametric angle convention; the parabolic arc is defined by vertex, focus, and aperture height, matching standard reflector-antenna convention.

7. Sources
Sources are defined in full 3D and projected onto the grid plane (a warning is shown if a source sits far from the plane's Z height). Three types are supported: point dipoles (with a 3D orientation vector), wire segments, and loops — the latter two support a spatial current taper (uniform or sinusoidal) along their length.

8. Visualization
Any of the six complex phasor components can be viewed through a selectable function: Re, Im, Arg (rad/deg), Abs, or Abs². Derived real-valued quantities — |E|, |H|, and the Poynting vector components Sx, Sy, Sz, |S| — are also available. Color scales support linear, log, and dB modes with auto-ranging (dB auto-range fixes the floor 60 dB below the computed peak) or fully manual min/max.

9. Data Export
Two independent JSON exports are available: a results export (steady-state phasor fields, far-field pattern and approximate directivity, and grid metadata) for downstream processing in Python/MATLAB/etc., and a project export (complete geometry, sources, and solver settings) for saving and reloading your setup.
