<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Channel Flow P=8, hydrostatic pressure</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>ChanFlow_Accelerate.xml</parameters>
    <files>
        <file description="Session File">ChanFlow_Accelerate.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">0</value>
            <value variable="v" tolerance="1e-12">0</value>
	    <value variable="p" tolerance="2e-9">0</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="2e-12">0</value>
            <value variable="v" tolerance="1e-12">0</value>
	    <value variable="p" tolerance="2e-9">0</value>
        </metric>
    </metrics>
</test>


