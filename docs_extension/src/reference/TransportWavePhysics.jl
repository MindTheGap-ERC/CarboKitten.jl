"""
    advection_coef!(box::Box{BT}, diffusivity, wave_velocity, w, adv, rct) where {BT}

Compute advection and reaction coefficients from a wave-velocity function.

# Arguments
- `box`: Spatial grid and boundary-condition definition.
- `diffusivity`: Diffusivity field or parameter.
- `wave_velocity`: Callable returning local wave-driven transport velocity.
- `w`: Water-depth field.
- `adv`: Preallocated array receiving advection coefficients.
- `rct`: Preallocated array receiving reaction or residual coefficients.

# Returns
- `nothing`

# Algorithm
1. Loop over all grid cells.
2. Evaluate local wave velocity from the water-depth field.
3. Convert velocity and diffusivity into directional advection coefficients.
4. Store the resulting coefficients in `adv` and any residual terms in `rct`.

# Units
- Velocity in the project transport units.
- Distance-like quantities in meters.
- Coefficients in units consistent with the discretized transport equation.

# Notes
This overload is used when wave velocity is supplied as a function of local state.
"""
function advection_coef! end

"""
    advection_coef!(box::Box{BT}, diffusivity,
                    fields::Tuple{AbstractArray,AbstractArray},
                    w, adv, rct) where {BT}

Compute advection and reaction coefficients from precomputed vector fields.

# Arguments
- `box`: Spatial grid and boundary-condition definition.
- `diffusivity`: Diffusivity field or parameter.
- `fields`: Tuple of precomputed transport fields, typically the x- and y-directed velocity components.
- `w`: Water-depth field.
- `adv`: Preallocated array receiving advection coefficients.
- `rct`: Preallocated array receiving reaction or residual coefficients.

# Returns
- `nothing`

# Algorithm
1. Loop over all grid cells.
2. Read the bounded local transport components from `fields`.
3. Combine local transport direction and diffusivity to form advection coefficients.
4. Store the resulting coefficients in `adv` and residual terms in `rct`.

# Units
- Velocity-like fields in the project transport units.
- Distance-like quantities in meters.
- Coefficients in units consistent with the transport discretization.

# Notes
This overload is used when transport fields have already been computed externally.
"""
function advection_coef! end



"""
    dir(θ) = Vec2(cos(θ), sin(θ))

Convert a propagation angle into a 2D unit direction vector.

# Arguments
- `θ`: Wave propagation angle (radians).

# Returns
- `Vec2(cos(θ), sin(θ))`: Unit vector pointing in the propagation direction.

# Algorithm
1. Evaluate the cosine and sine of the input angle.
2. Construct a 2D vector with x-component `cos(θ)` and y-component `sin(θ)`.

# Units
- `θ` in radians.
- Output is dimensionless.

# Notes
Used to project scalar wave properties into directional transport components.
"""
function dir(θ) end



"""
    _wave_velocity_at_depth(model::WaveModel, pos, t, h, x_fetch)

Compute the wave orbital velocity at a given depth and location.

# Arguments
- `model`: Wave model containing period, amplitude, and attenuation parameters.
- `pos`: Spatial position used to evaluate directional or spatial variation.
- `t`: Simulation time.
- `h`: Local water depth.
- `x_fetch`: Effective fetch controlling wave growth.

# Returns
- 2D velocity vector representing near-bed or depth-dependent wave orbital motion.

# Algorithm
1. Compute wavelength `L` from the wave period using `wavelength_from_period(T, h)`.
2. Compute wave number `k = 2π / L`.
3. Evaluate vertical decay of orbital velocity using an exponential attenuation:
   - Velocity ∝ exp(-k * depth)
4. Scale the velocity magnitude based on fetch-dependent wave growth.
5. Project the scalar velocity onto a direction vector using `dir(...)`.
6. Return the resulting 2D velocity vector.

# Units
- Depth and wavelength in meters.
- Velocity in length per unit time consistent with the model forcing.

# Notes
This function captures the depth decay of wave orbital motion, which controls sediment transport intensity.
"""
function _wave_velocity_at_depth end

"""
    wavelength_from_period(T, h; tol=1e-10, maxiter=100)

Solve the linear wave dispersion relation for wavelength.

# Arguments
- `T`: Wave period.
- `h`: Water depth.
- `tol`: Convergence tolerance for the iterative solver.
- `maxiter`: Maximum number of iterations.

# Returns
- Wavelength `L` satisfying the dispersion relation.

# Algorithm
1. Convert inputs to consistent units if needed (`ustrip`).
2. Solve the dispersion relation:
   - ω² = gk tanh(kh), where ω = 2π / T
3. Use an iterative root-finding scheme:
   a. Initialize wave number `k`.
   b. Update `k` using the dispersion equation.
   c. Evaluate convergence using `abs(Δk) < tol`.
4. Convert wave number to wavelength:
   - `L = 2π / k`
5. Return `L`.

# Units
- `T` in time units.
- `h` and `L` in meters.

# Notes
Handles both shallow- and deep-water limits through the `tanh(kh)` term.
"""
function wavelength_from_period end

"""
    wave_physics_at_cell(model::WaveModel, pos, t, h, x_fetch; dh=0.05u"m")

Evaluate wave-induced hydrodynamic properties at a grid cell.

# Arguments
- `model`: Wave model containing physical parameters.
- `pos`: Spatial position of the cell.
- `t`: Simulation time.
- `h`: Local water depth.
- `x_fetch`: Effective fetch controlling wave development.
- `dh`: Vertical increment used for estimating gradients.

# Returns
- Tuple or structured result containing wave velocity and derived quantities
  such as energy or shear proxies.

# Algorithm
1. Compute orbital velocity at the target depth using `_wave_velocity_at_depth`.
2. Optionally evaluate velocity at a slightly offset depth (`h + dh`) to estimate vertical gradients.
3. Use exponential decay relationships to compute near-bed velocity magnitude.
4. Derive secondary quantities (e.g. energy proxies) from velocity magnitude:
   - Energy ∝ velocity²
5. Apply non-negativity constraints using `max(...)` where needed.
6. Return velocity and derived diagnostics.

# Units
- Depth in meters.
- Velocity in model transport units.
- Derived energy follows project-specific conventions.

# Notes
This routine bridges wave physics and sediment transport by providing physically
consistent forcing fields at each grid cell.
"""
function wave_physics_at_cell end
