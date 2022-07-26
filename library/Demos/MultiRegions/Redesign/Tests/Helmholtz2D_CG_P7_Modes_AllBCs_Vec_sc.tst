<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Helmholtz 2D CG with P=7, direct sc, vector system</description>
    <executable>Helmholtz</executable>
    <parameters>-I GlobalSysSoln=DirectStaticCond Helmholtz2D_P7_AllBCs_Vec.xml</parameters>
    <files>
        <file description="Session File">Helmholtz2D_P7_AllBCs_Vec.xml</file>
    </files>

    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-7">0.00888406</value>
            <value variable="v" tolerance="1e-7">0.00888406</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-6">0.0102587</value>
            <value variable="v" tolerance="1e-6">0.0102587</value>
        </metric>
    </metrics>
</test>


