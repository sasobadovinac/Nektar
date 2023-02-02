<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>NS, Couette flow with non-conformal circle, exact solution</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>Movement_fixed_couette.xml</parameters>
    <files>
        <file description="Session File">Movement_fixed_couette.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-6">1.30059e-10</value>
            <value variable="rhou" tolerance="1e-6">1.64126e-07</value>
            <value variable="rhov" tolerance="1e-6">1.22606e-07</value>
            <value variable="E" tolerance="1e-6">3.7937e-05</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-6">1.65258e-09</value>
            <value variable="rhou" tolerance="1e-6">3.70148e-06</value>
            <value variable="rhov" tolerance="1e-6">2.16739e-06</value>
            <value variable="E" tolerance="1e-6">0.000522754</value>
        </metric>
    </metrics>
</test>
