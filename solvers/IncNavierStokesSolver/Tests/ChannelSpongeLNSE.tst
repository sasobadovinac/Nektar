<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Linearized Channel Flow P=11</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>ChannelSpongeLNSE.xml</parameters>
    <files>
        <file description="Session File">ChannelSpongeLNSE.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">0.0003601</value>
            <value variable="v" tolerance="1e-12">0.000186492</value>
	    <value variable="p" tolerance="1e-12">0.000232481</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">0.000924247</value>
            <value variable="v" tolerance="1e-12">0.000576185</value>
	    <value variable="p" tolerance="1e-12">0.000671239</value>
        </metric>
    </metrics>
</test>


