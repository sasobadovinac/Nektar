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
            <value variable="rho" tolerance="1e-11">1.87657e-06</value>
            <value variable="rhou" tolerance="1e-9">1.54271e-04</value>
            <value variable="rhov" tolerance="1e-9">3.89585e-04</value>
            <value variable="E" tolerance="1e-6">0.580575</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-11">4.32354e-06</value>
            <value variable="rhou" tolerance="1e-9">4.7983e-04</value>
            <value variable="rhov" tolerance="1e-8">1.28929e-03</value>
            <value variable="E" tolerance="1e-5">1.36531</value>
        </metric>
    </metrics>
</test>


