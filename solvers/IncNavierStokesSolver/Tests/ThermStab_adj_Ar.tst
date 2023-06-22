<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Linear stability (Arpack): ThermalChan (adjoint)</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>ThermStab_adj_Ar.xml</parameters>
    <files>
        <file description="Session File">ThermStab_adj_Ar.xml</file>
        <file description="Session File">ThermStab_adj_Ar.bse</file>
        <file description="Session File">ThermStab_adj_Ar.rst</file>
    </files>
    <metrics>
        <metric type="Eigenvalue" id="0">
            <value tolerance="0.001">1.06675,0.00</value>
        </metric>
    </metrics>
</test>
