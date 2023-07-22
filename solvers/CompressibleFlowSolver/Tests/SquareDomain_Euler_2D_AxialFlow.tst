<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Uniform Flow Through a Square Domain</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>SquareDomain_Euler_2D_AxialFlow.xml</parameters>
    <files>
        <file description="Session File">SquareDomain_Euler_2D_AxialFlow.xml</file>
    </files>
    <!--
        This test is designed to test the implementation of the StagnationInflowBC

        If the boundary condition works as expected, this test should give a uniform
        flow in the axial direction, with a Mach number of 0.8.

        In this test, the flow direction is not specified by the user. Since this is
        not done, the boundary condition imposes a flow which is normal to the boundary, 
        which in this case is in the positive x-direction.

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
            <value variable="E"    tolerance="5e-14">0.0</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho"  tolerance="5e-14">0.0</value>
            <value variable="rhou" tolerance="5e-14">0.0</value>
            <value variable="rhov" tolerance="5e-14">0.0</value>
            <value variable="E"    tolerance="5e-14">0.0</value>
        </metric>
    </metrics>
</test>