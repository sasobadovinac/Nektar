<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>CG Helmholtz3D Homogeneous 2D</description>
    <executable>Helmholtz3DHomo2D</executable>
    <parameters>Helmholtz3D_Homo2D.xml</parameters>
    <files>
        <file description="Session File">Helmholtz3D_Homo2D.xml</file>
    </files>

    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-11">1.00372e-06</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-11">6.32758e-06</value>
        </metric>
    </metrics>
</test>


