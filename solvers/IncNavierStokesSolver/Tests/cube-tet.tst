<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D Flow inside cube. Tet elements.This test case requires at least 6 cpus to test that the deadlock problem doesn't appear.</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>cube-tet.xml</parameters>
    <processes>6</processes>
    <files>
        <file description="Session File">cube-tet.xml</file>
    </files>
    <metrics>
<!-- The metric values are very high because they are not important. -->
<!-- Because this test-case is used only for checking if the deadlock problem appears or not -->
        <metric type="L2" id="1">
                <value variable="u" tolerance="1">0.002</value> 
            <value variable="v" tolerance="1e-3">0.00200113</value>
	    <value variable="p" tolerance="1e-3">0.0123557 </value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1">0.00959985</value>
            <value variable="v" tolerance="1e-3">0.00582256</value>
	    <value variable="p" tolerance="1e-3">0.0502862</value>
        </metric>
    </metrics>
</test>
