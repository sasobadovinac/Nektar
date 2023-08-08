<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Euler Isentropic Vortex P=8, Implicit</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>IsentropicVortex16Implicit_P8.xml</parameters>
    <files>
        <file description="Session File">IsentropicVortex16Implicit_P8.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-07">1.07852e-06</value>
            <value variable="rhou" tolerance="1e-07">2.51391e-06</value>
            <value variable="rhov" tolerance="1e-07">2.63564e-06</value>
            <value variable="E" tolerance="1e-07">8.23237e-06</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-07">6.2691e-06</value>
            <value variable="rhou" tolerance="1e-07">7.87288e-06</value>
            <value variable="rhov" tolerance="1e-07">7.53261e-06</value>
            <value variable="E" tolerance="1e-06">2.61552e-05</value>
        </metric>
    </metrics>
</test>


