<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Cross field solver</description>
    <executable>ADRSolver</executable>
    <parameters>CrossField_half_disc.xml</parameters>
    <files>
        <file description="Session File">CrossField_half_disc.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-08">0.352733</value>
            <value variable="v" tolerance="1e-08">0.198164</value>
        </metric>
        <metric type="Linf" id="1">
            <value variable="u" tolerance="1e-08">1</value>
            <value variable="v" tolerance="1e-08">1</value>
        </metric>
    </metrics>
</test>
