<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Forced Vortex Swirl in a 2D Duct</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>Duct_Euler_2D_AxiSym_ForcedVortexSwirl.xml</parameters>
    <files>
        <file description="Session File">Duct_Euler_2D_AxiSym_ForcedVortexSwirl.xml</file>
    </files>
    <!--
        This test is intended to test the implementation of the axisymmetric forcing.

        If everything works correctly, this test should give a forced vortex swirl
        inside the duct. Since a forced vortex swirl has an analytical solution which
        is a P=2 polynomial, the solution obtained with Nektar++ should furthermore
        be exact, i.e., the L2 and Linf errors must be below 1.E-14. If this is not
        the case, the performance of the solver has deteriorated, or the axisymmetrix
        forcing is not correctly implemented.

        For more details about this test case, please read the comments in the XML
        file for this test.
    -->
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho"  tolerance="1e-14">0.0</value>
            <value variable="rhou" tolerance="1e-14">0.0</value>
            <value variable="rhov" tolerance="1e-14">0.0</value>
            <value variable="rhow" tolerance="1e-14">0.0</value>
            <value variable="E"    tolerance="1e-14">0.0</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho"  tolerance="1e-14">0.0</value>
            <value variable="rhou" tolerance="1e-14">0.0</value>
            <value variable="rhov" tolerance="1e-14">0.0</value>
            <value variable="rhow" tolerance="1e-14">0.0</value>
            <value variable="E"    tolerance="1e-14">0.0</value>
        </metric>
    </metrics>
</test>