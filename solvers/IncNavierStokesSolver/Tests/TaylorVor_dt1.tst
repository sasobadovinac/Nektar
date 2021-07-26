<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Convergence: Taylor Vortex IMEXOrder2 dt=0.01</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>TaylorVor_dt1.xml</parameters>
    <files>
        <file description="Session File">TaylorVor_dt1.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">5.9519e-06</value>
            <value variable="v" tolerance="1e-12">4.99594e-06</value>
            <value variable="p" tolerance="1e-12">0.000232626</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">4.15477e-06</value>
            <value variable="v" tolerance="1e-12">4.51605e-06</value>
            <value variable="p" tolerance="1e-12">0.000162441</value>
        </metric>
    </metrics>
</test>
