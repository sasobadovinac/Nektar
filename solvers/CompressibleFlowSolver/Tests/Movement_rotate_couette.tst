<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>NS, Couette flow with non-conformal circle rotating, exact solution</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>Movement_rotate_couette.xml</parameters>
    <files>
        <file description="Session File">Movement_rotate_couette.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-6">1.39834e-08</value>
            <value variable="rhou" tolerance="1e-6">1.79132e-05</value>
            <value variable="rhov" tolerance="1e-6">1.38303e-05</value>
            <value variable="E" tolerance="1e-5">0.00433813</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-6">2.60889e-07</value>
            <value variable="rhou" tolerance="1e-5">0.000781472</value>
            <value variable="rhov" tolerance="1e-6">0.000391911</value>
            <value variable="E" tolerance="1e-3"> 0.0716392</value>
        </metric>
    </metrics>
</test>
