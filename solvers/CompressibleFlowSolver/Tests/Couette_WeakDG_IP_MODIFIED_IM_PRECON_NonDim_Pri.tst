<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>
        Compressible Couette flow for a calorically perfect, ideal gas with constant viscosity.

        The test case is intended to check that the CompressibleFlowSolver produces correct results
        when the Interior Penalty method is used with prism elements. Currently, the
        non-symmetric IP diffusion implementation is tested and the preconditioner is enabled.

        The equations are solved in non-dimensional form. The domain is three dimensional. At the top
        and bottom, the walls are simulated by specifying the conservative variables in such a way
        that the correct velocity, temperature and pressure is obtained. Note that a ViscousWall b.c.
        is not used. A the front and back, slip wall boundary conditions are used, and at the left and
        right boundary, a periodic condition is imposed.

        The analytical solution is used as initial conditions and the solver is run for 200 time steps.
        It has been checked that the soltuion does not change after 200 time steps with the selected
        time step size (dt). It has also been checked that the solution is insensitive to changes in 
        tolerances used for the implicit solver.
    </description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>Couette_WeakDG_IP_MODIFIED_IM_PRECON_NonDim_Pri.xml</parameters>
    <files>
        <file description="Session File">Couette_WeakDG_IP_MODIFIED_IM_PRECON_NonDim_Pri.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho"  tolerance="5e-8" >3.73741e-05</value>
            <value variable="rhou" tolerance="5e-8" >0.000110546</value>
            <value variable="rhov" tolerance="5e-8" >1.64516e-05</value>
            <value variable="rhow" tolerance="5e-12">1.26586e-15</value>
            <value variable="E"    tolerance="5e-8" >4.21625e-05</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho"  tolerance="5e-8" >0.00017523</value>
            <value variable="rhou" tolerance="5e-8" >0.000177817</value>
            <value variable="rhov" tolerance="5e-8" >4.25271e-05</value>
            <value variable="rhow" tolerance="5e-12">1.71785e-15</value>
            <value variable="E"    tolerance="5e-8" >0.000213984</value>
        </metric>
    </metrics>
</test>
