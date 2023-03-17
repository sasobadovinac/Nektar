<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Uniform Flow Through a Square Domain</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>SquareDomain_Euler_2D_DiagonalFlow.xml</parameters>
    <files>
        <file description="Session File">SquareDomain_Euler_2D_DiagonalFlow.xml</file>
    </files>
    <!--
        This test is designed to test the implementation of the StagnationInflowBC

        If the boundary condition works as expected, this test should give a uniform
        flow in the x/z direction, with a 45 degree angle to the x-axis. To do this,
        we solve for all components of the momentum equation (u,v,w), and impose a flow
        direction at the inlet which is at a 45 degree angle to the x-axis. Note that
        this means that the flow is flowing "out of the plane". The Mach number for this
        case is M = sqrt(u^2 + v^2 + w^2)/c = 0.8.

        In this test, the flow direction is specified by the user by specifying the
        rhou, rhov, and rhow components. The magnitude of these values doesn't matter,
        all that matters is the direction that the three components define.

        Since this test has an analytic solution, which is simply a uniform flow at
        M = 0.8, this test should result in an error which is less than 5E-14. If this
        doesn't occur, there is definitely an error in the code, possibly related to
        the boundary condition. Another issue could be that the solver converges slower
        than when this test was originally written. The reason for this is that the
        flow isn't initialized with the analytical soltuion, but with a slower velocity.
        Hence, if the solver doesn't converge fast enough, the flow won't reach the
        right solution in the 500 iterations that are specified.
    -->
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho"  tolerance="5e-14">0.0</value>
            <value variable="rhou" tolerance="5e-14">0.0</value>
            <value variable="rhov" tolerance="5e-14">0.0</value>
            <value variable="rhow" tolerance="5e-14">0.0</value>
            <value variable="E"    tolerance="5e-14">0.0</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho"  tolerance="5e-14">0.0</value>
            <value variable="rhou" tolerance="5e-14">0.0</value>
            <value variable="rhov" tolerance="5e-14">0.0</value>
            <value variable="rhow" tolerance="5e-14">0.0</value>
            <value variable="E"    tolerance="5e-14">0.0</value>
        </metric>
    </metrics>
</test>