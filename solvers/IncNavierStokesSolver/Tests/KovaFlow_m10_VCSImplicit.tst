<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Kovasznay Flow P=9 with large dt</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>KovaFlow_m10_VCSImplicit.xml</parameters>
    <files>
        <file description="Session File">KovaFlow_m10_VCSImplicit.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-08">2.16592e-08</value>
            <value variable="v" tolerance="1e-08">1.29113e-08</value>
            <value variable="p" tolerance="1e-08">1.96184e-08</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-08">7.68403e-08</value>
            <value variable="v" tolerance="1e-08">2.96839e-08</value>
            <value variable="p" tolerance="1e-07">2.63382e-07</value>
        </metric>
    </metrics>
</test>
