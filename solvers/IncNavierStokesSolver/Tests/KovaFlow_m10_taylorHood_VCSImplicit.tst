<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Kovasznay Flow P=9 with large dt and taylor-hood expansion</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>KovaFlow_m10_taylorHood_VCSImplicit.xml</parameters>
    <files>
        <file description="Session File">KovaFlow_m10_taylorHood_VCSImplicit.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-08">2.5329e-08</value>
            <value variable="v" tolerance="1e-08">1.35236e-08</value>
            <value variable="p" tolerance="1e-08">2.17508e-08</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-08">8.56193e-08</value>
            <value variable="v" tolerance="1e-08">3.94227e-08</value>
            <value variable="p" tolerance="1e-07">2.7661e-07</value>
        </metric>
    </metrics>
</test>
