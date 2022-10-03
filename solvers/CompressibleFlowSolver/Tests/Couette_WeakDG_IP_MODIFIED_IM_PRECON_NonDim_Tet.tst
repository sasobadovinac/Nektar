<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>
        Compressible Couette flow for a calorically perfect, ideal gas with constant viscosity.

        The test case is intended to check that the CompressibleFlowSolver produces correct results
        when the Interior Penalty method is used with tetrahedral elements. Currently, the
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
    <parameters>Couette_WeakDG_IP_MODIFIED_IM_PRECON_NonDim_Tet.xml</parameters>
    <files>
        <file description="Session File">Couette_WeakDG_IP_MODIFIED_IM_PRECON_NonDim_Tet.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho"  tolerance="5e-8">2.73418e-05</value>
            <value variable="rhou" tolerance="5e-8">7.70266e-05</value>
            <value variable="rhov" tolerance="5e-8">1.61997e-05</value>
            <value variable="rhow" tolerance="5e-8">1.40221e-05</value>
            <value variable="E"    tolerance="5e-8">2.97967e-05</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho"  tolerance="5e-8">0.000314853</value>
            <value variable="rhou" tolerance="5e-8">0.000648129</value>
            <value variable="rhov" tolerance="5e-8">0.00015004 </value>
            <value variable="rhow" tolerance="5e-8">5.78513e-05</value>
            <value variable="E"    tolerance="5e-8">0.00033257 </value>
        </metric>
    </metrics>
</test>
