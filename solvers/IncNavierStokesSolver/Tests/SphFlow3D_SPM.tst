<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>3D sphere flow, P=3, sphere defined via IB SPM function</description>    
    <executable>IncNavierStokesSolver</executable>
    <parameters>SphFlow3D_SPM.xml</parameters>
    <files>
        <file description="Session File">SphFlow3D_SPM.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-6">88.124</value>
            <value variable="v" tolerance="1e-6">0.336273</value>
            <value variable="w" tolerance="1e-6">0.336273</value>
            <value variable="p" tolerance="1e-6">2.32154</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-6">1.34077</value>
            <value variable="v" tolerance="1e-6">0.753533</value>
            <value variable="w" tolerance="1e-6">0.753529</value>
            <value variable="p" tolerance="1e-6">2.56945</value>
        </metric>
    </metrics>
</test>
