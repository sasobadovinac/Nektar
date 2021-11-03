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
            <value variable="u" tolerance="1e-12">0.343782</value>
            <value variable="v" tolerance="1e-12">0.026204</value>
	    <value variable="p" tolerance="1e-12">0.340367</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">0.322039</value>
            <value variable="v" tolerance="1e-12">0.0264658</value>
	    <value variable="p" tolerance="1e-12">0.491158</value>
        </metric>
    </metrics>
</test>
