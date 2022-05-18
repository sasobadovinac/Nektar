<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Moving reference frame formualtion of NS equation with rotation, simple domain without body</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>MovingRefFrame_Rot_SimpleDomain.xml</parameters>
    <files>
        <file description="Session File">MovingRefFrame_Rot_SimpleDomain.xml</file>
        <file description="Restart File">MovingRefFrame_Rot_SimpleDomain.rst</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="5e-4">12.2445</value>
            <value variable="v" tolerance="5e-5">0.267844</value>
            <value variable="p" tolerance="5e-5">0.0335763</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="5e-4">1.00069</value>
            <value variable="v" tolerance="5e-4">0.0223664</value>
            <value variable="p" tolerance="5e-4">0.00944716</value>
        </metric>
    </metrics>
</test>
