<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>NS, Couette flow with non-conformal circle rotating, exact solution, parallel</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>Movement_rotate_couette.xml</parameters>
    <processes>4</processes>
    <files>
        <file description="Session File">Movement_rotate_couette.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-6">1.26103e-08</value>
            <value variable="rhou" tolerance="1e-6">1.56497e-05</value>
            <value variable="rhov" tolerance="1e-6">1.07509e-05</value>
            <value variable="E" tolerance="1e-5">0.00393699</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-6">2.67646e-07</value>
            <value variable="rhou" tolerance="1e-5">0.000783763</value>
            <value variable="rhov" tolerance="1e-6">0.000442926</value>
            <value variable="E" tolerance="1e-3">0.0778952</value>
        </metric>
    </metrics>
</test>
