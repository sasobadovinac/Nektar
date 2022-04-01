<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D Tet Kovasnay solution using GJP Stabilisation and dealiasing</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>Tet_Kovasnay_GJP.xml</parameters>
    <files>
        <file description="Session File">Tet_Kovasnay_GJP.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-6">0.00115918</value>
            <value variable="v" tolerance="1e-6">0.0003145</value>
            <value variable="w" tolerance="1e-6">0.000231965</value>
            <value variable="p" tolerance="1e-6">0.00125287</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-6">0.00900795</value>
            <value variable="v" tolerance="1e-6">0.0014178</value>
            <value variable="w" tolerance="1e-6">0.00103127</value>
            <value variable="p" tolerance="1e-6">0.00615774</value>
        </metric>
    </metrics>
</test>
