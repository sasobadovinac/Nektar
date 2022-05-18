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
            <value variable="rho" tolerance="1e-11">1.54893e-06</value>
            <value variable="rhou" tolerance="1e-9">1.60968e-04</value>
            <value variable="rhov" tolerance="1e-9">3.98249e-04</value>
            <value variable="E" tolerance="1e-6">0.589307</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-11">3.62604e-06</value>
            <value variable="rhou" tolerance="1e-9">4.69435e-04</value>
            <value variable="rhov" tolerance="1e-8">1.34506e-03</value>
            <value variable="E" tolerance="1e-5">1.43277</value>
        </metric>
    </metrics>
</test>


