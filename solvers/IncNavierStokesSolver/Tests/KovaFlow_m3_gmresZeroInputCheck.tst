<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Kovasznay Flow P=3 with zero input to GMRES for higher modes</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>KovaFlow_m3_gmresZeroInputCheck.xml</parameters>
    <files>
        <file description="Session File">KovaFlow_m3_gmresZeroInputCheck.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">0.1623590</value>
            <value variable="v" tolerance="1e-12">0.0186630</value>
            <value variable="w" tolerance="1e-12">0.0000000</value>
            <value variable="p" tolerance="1e-12">0.1015520</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">0.3220390</value>
            <value variable="v" tolerance="1e-12">0.0638199</value>
            <value variable="w" tolerance="1e-12">0.0000000</value>
            <value variable="p" tolerance="1e-12">0.3691620</value>
        </metric>
    </metrics>
</test>
