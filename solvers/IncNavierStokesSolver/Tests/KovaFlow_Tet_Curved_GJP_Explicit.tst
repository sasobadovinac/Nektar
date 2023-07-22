<?xml version="1.0" encoding="utf-8"?>
<test>
3   <description>3D Kovasznay flow on curved Tetrahedrons with explicit GJP stabilisation</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>KovaFlow_Tet_Curved.xml -I GJPStabilisation=Explicit -I GJPNormalVelocity=True</parameters>
    <files>
        <file description="Session File">KovaFlow_Tet_Curved.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-14">0.000543901</value>
            <value variable="v" tolerance="1e-14">0.000254741</value>
            <value variable="w" tolerance="1e-14">0.000170741</value>
            <value variable="p" tolerance="1e-14">0.00178814</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-14">0.0020482</value>
            <value variable="v" tolerance="1e-14">0.000806377</value>
            <value variable="w" tolerance="1e-14">0.000652387</value>
            <value variable="p" tolerance="1e-14">0.00930501</value>
        </metric>
    </metrics>
</test>
