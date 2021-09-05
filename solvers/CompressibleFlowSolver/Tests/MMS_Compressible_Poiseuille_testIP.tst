<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>NS, Manufactured Compressible Poiseuille's flow to test IP</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>MMS_Compressible_Poiseuille_testIP.xml</parameters>
    <files>
        <file description="Session File">MMS_Compressible_Poiseuille_testIP.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">1.87657e-06</value>
            <value variable="rhou" tolerance="1e-12">0.000154271</value>
            <value variable="rhov" tolerance="1e-12">0.000389585</value>
            <value variable="E" tolerance="1e-12">0.580575</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">4.32354e-06</value>
            <value variable="rhou" tolerance="1e-12">0.00047983</value>
            <value variable="rhov" tolerance="1e-12">0.00128929</value>
            <value variable="E" tolerance="1e-12">1.36531</value>
        </metric>
    </metrics>
</test>


