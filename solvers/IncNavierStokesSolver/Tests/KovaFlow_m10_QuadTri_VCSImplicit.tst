<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Kovasznay flow test case with mixed Quad-Tri elements, NUMMODES=10 and dt=0.4</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>KovaFlow_m10_QuadTri_VCSImplicit.xml</parameters>
    <files>
        <file description="Session File">KovaFlow_m10_QuadTri_VCSImplicit.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-09">3.55915e-07</value>
            <value variable="v" tolerance="1e-09">1.16281e-07</value>
            <value variable="p" tolerance="1e-09">2.15576e-07</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-09">7.90170e-07</value>
            <value variable="v" tolerance="1e-09">2.79123e-07</value>
            <value variable="p" tolerance="1e-09">1.87086e-06</value>
        </metric>
    </metrics>
</test>
