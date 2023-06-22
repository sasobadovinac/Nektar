<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Linear stability (Mod. Arnoldi): ThermalChan</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>ThermStab_dir.xml</parameters>
    <files>
        <file description="Session File">ThermStab_dir.xml</file>
        <file description="Session File">ThermStab_dir.bse</file>
        <file description="Session File">ThermStab_dir.rst</file>
    </files>
    <metrics>
        <metric type="Eigenvalue" id="0">
            <value tolerance="0.001">1.0667,0.00</value>
        </metric>
    </metrics>
</test>
