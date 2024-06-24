# Active Layer Transport

The following is inspired on well-known **active layer** approaches in river bed sediment transport.

We suppose that loose sediment, either fresh production or disintegrated older sediment, is being transported in a layer on top of the sea bed. The flux in this layer is assumed to be directly proportional to the local slope of the sea bed $\partial_x \eta$.

The active layer now contains a concentration $C_f$ particles of different grain size (for each facies $f$). If needed, $C_f = \alpha_f P_f$ where $\alpha_f$ is some facies parameter determining the fraction of production that is available for transport. The flux between two cells is then given as

$$q = -\nu_f C_f {{\partial \Sum_f \eta_f} \over {\partial x}}.$$

The flux again measures the change in height as contributed by each facies type. The following is a mass balance:

$$\sigma + \Sum_f {{\partial \eta_f} \over {\partial t}} = -\Sum_f {{\partial q_f} \over {\partial x}} + P_f,$$

where $\sigma$ is the subsidence rate. In other approaches to active layer transport there would be a factor $1/C_f$. Here we have a different interpretation: the sediment settles down after transport, such that the concentration has no impact on the change in sediment surface elevation.

Combining these equations, and ignoring subsidence for the moment, we get a component-wise diffusion equation

$$\partial_t \eta_f(x) = \partial_x (\nu_f \alpha_f P_f(x) \partial_x \eta_{*}(x)) + P_f(x),$$

where $\eta_{*} = \sum_f \eta_f$ is the total sediment surface elevation.

In our model we need to solve this equation one time-step each iteration. If we solve this using forward methods, we should be reminded of the CFL limit for diffusion equations (depending on the diffusion constants and grid size we shouldn't pick the time steps too large). Alternatively, for these two-dimensional situations, an implicit approach is feasible.
