<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Moving reference frame formualtion of NS equation with rotation, simple domain without body</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>MovingRefFrame_Rot_naca0012.xml</parameters>
    <files>
        <file description="Session File">MovingRefFrame_Rot_naca0012.xml</file>
        <file description="Restart File">MovingRefFrame_Rot_naca0012.rst</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="5e-4">30.0106</value>
            <value variable="v" tolerance="5e-5">0.434642</value>
            <value variable="p" tolerance="5e-3">113.124</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="5e-4">10.0129</value>
            <value variable="v" tolerance="5e-4">13.9793</value>
            <value variable="p" tolerance="5e-3">273.927</value>
        </metric>
    </metrics>
</test>
