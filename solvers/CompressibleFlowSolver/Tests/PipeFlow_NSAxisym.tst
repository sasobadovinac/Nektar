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
            <value variable="rho" tolerance="1e-9">1.6481e-05</value>
            <value variable="rhou" tolerance="1e-8">0.00301672</value>
            <value variable="rhov" tolerance="1e-8">0.00267852</value>
            <value variable="E" tolerance="1e-4">4.4922</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-9">2.0097e-05</value>
            <value variable="rhou" tolerance="1e-8">0.00265534</value>
            <value variable="rhov" tolerance="1e-8">0.00312494</value>
            <value variable="E" tolerance="1e-4">5.4170</value>
        </metric>
    </metrics>
</test>

