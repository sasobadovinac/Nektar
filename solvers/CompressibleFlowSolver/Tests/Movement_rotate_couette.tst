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
            <value variable="rho" tolerance="1e-6">7.01766e-08</value>
            <value variable="rhou" tolerance="1e-6">1.79732e-05</value>
            <value variable="rhov" tolerance="1e-6">1.55277e-05</value>
            <value variable="E" tolerance="1e-5">0.0173705</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-6">7.14624e-07</value>
            <value variable="rhou" tolerance="1e-5">0.000131629</value>
            <value variable="rhov" tolerance="1e-6">8.52884e-05</value>
            <value variable="E" tolerance="1e-3">0.0579378</value>
        </metric>
    </metrics>
</test>
