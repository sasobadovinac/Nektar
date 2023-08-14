<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>NS, axisymmetric pipe flow with low Mach number</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>PipeFlow_NSAxisym.xml</parameters>
    <files>
        <file description="Session File"> PipeFlow_NSAxisym.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-11">1.48594e-06</value>
            <value variable="rhou" tolerance="1e-9">0.000141777</value>
            <value variable="rhov" tolerance="1e-8">0.00284738</value>
            <value variable="E" tolerance="1e-6">0.346929</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">9.28064e-07</value>
            <value variable="rhou" tolerance="1e-9">0.000238733</value>
            <value variable="rhov" tolerance="1e-8">0.00297347</value>
            <value variable="E" tolerance="1e-6">0.183325</value>
        </metric>
    </metrics>
</test>

