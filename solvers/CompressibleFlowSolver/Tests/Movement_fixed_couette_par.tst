<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>NS, Couette flow with non-conformal circle, exact solution, parallel</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>Movement_fixed_couette.xml</parameters>
    <processes>4</processes>
    <files>
        <file description="Session File">Movement_fixed_couette.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-6">1.33923e-10</value>
            <value variable="rhou" tolerance="1e-6">1.43646e-07</value>
            <value variable="rhov" tolerance="1e-6">9.86528e-08</value>
            <value variable="E" tolerance="1e-5">3.90118e-05</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-6">1.87859e-09</value>
            <value variable="rhou" tolerance="1e-5">3.65872e-06</value>
            <value variable="rhov" tolerance="1e-6">1.58245e-06</value>
            <value variable="E" tolerance="1e-3">0.000620003</value>
        </metric>
    </metrics>
</test>
