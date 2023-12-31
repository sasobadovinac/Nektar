<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Kovasznay Flow P=8 with Higher-Order outflow boundary condition and Implicit scheme</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>KovaFlow_m8_short_HOBC.xml -I SolverType=VCSImplicit</parameters>
    <files>
        <file description="Session File">KovaFlow_m8_short_HOBC.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-10">2.51953e-08</value>
            <value variable="v" tolerance="1e-10">9.56014e-09</value>
            <value variable="p" tolerance="1e-10">1.11025e-08</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-9">9.47698e-08</value>
            <value variable="v" tolerance="1e-9">5.59175e-08</value>
            <value variable="p" tolerance="1e-8">2.93234e-07</value>
        </metric>
    </metrics>
</test>
