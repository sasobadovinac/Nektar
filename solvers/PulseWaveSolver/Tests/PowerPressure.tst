<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>PowerPressureArea, P=5</description>
    <executable>PulseWaveSolver</executable>
    <parameters>PowerPressure.xml</parameters>
    <files>
        <file description="Session File">PowerPressure.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="A" tolerance="1e-12">7.43089</value>
            <value variable="u" tolerance="1e-12">87.0556</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="A" tolerance="1e-12">1.81683</value>
            <value variable="u" tolerance="1e-12">40.3152</value>
        </metric>
    </metrics>
</test>
