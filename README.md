# cavity-simulation

Laser cavity simulation for ORACLE project.

Propagation advances in elapsed time using the fixed `sim.dt` timestep. The
beam's spatial coordinate (`laser.pos`) remains a geometric state for mirror
locations and spatial diffraction; each step advances it by `consts.c * sim.dt`
in the direction given by `sim.direction`, with a fractional timestep at each
mirror. `laser.time` increases monotonically through reflections, and
`outputs.ts` records elapsed time for each propagation sample. This provides a
time coordinate for future time-dependent effects without changing the
existing spatial or gain-medium model.