<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Channel Flow P=3</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>ChanFlow_m3_abstol.xml</parameters>
    <files>
        <file description="Session File">ChanFlow_m3.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-10">7.01969e-11</value>
            <value variable="v" tolerance="1e-9">1.12895e-10</value>
	    <value variable="p" tolerance="1e-8">2.05227-9</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-9">6.12299-10</value>
            <value variable="v" tolerance="1e-9">4.90306-10</value>
	    <value variable="p" tolerance="1e-8">6.60663-09</value>
        </metric>
    </metrics>
</test>


