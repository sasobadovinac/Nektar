<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D Kovasznay flow on curved Tetrahedrons with Semi-Implicit GJP stabilisation</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>KovaFlow_Tet_Curved.xml -I GJPStabilisation=SemiImplicit</parameters>
    <files>
        <file description="Session File">KovaFlow_Tet_Curved.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-14">0.000792368</value>
            <value variable="v" tolerance="1e-14">0.000323808</value>
            <value variable="w" tolerance="1e-14">0.000196171</value>
            <value variable="p" tolerance="1e-14">0.002642070</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-14">0.005747140</value>
            <value variable="v" tolerance="1e-14">0.002341920</value>
            <value variable="w" tolerance="1e-14">0.000904156</value>
            <value variable="p" tolerance="1e-14">0.015100600</value>
        </metric>
    </metrics>
</test>
