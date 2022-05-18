<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Kovasznay Flow P=3 with explicit GJP and Normal velocity</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>KovaFlow_m3_GJP.xml</parameters>
    <files>
        <file description="Session File">KovaFlow_m3_GJP.xml</file>
	<file description="Session File">KovaFlow_m3.rst</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-6">0.344283</value>
            <value variable="v" tolerance="1e-6">0.0263953</value>
	    <value variable="p" tolerance="1e-6">0.339166</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-6">0.322039</value>
            <value variable="v" tolerance="1e-6">0.0273741</value>
	    <value variable="p" tolerance="1e-6">0.491786</value>
        </metric>
    </metrics>
</test>
