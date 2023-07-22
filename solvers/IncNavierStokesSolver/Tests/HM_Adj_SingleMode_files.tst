<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Fourier Half Mode Adjoint Basis, P=3</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>HM_Adj_SingleMode_files.xml</parameters>
    <files>
        <file description="Session File">HM_Adj_SingleMode_files.xml</file>
        <file description="InitialConditions File">HM_Adj_SingleMode_init.rst</file>
        <file description="BodyForce File">HM_Adj_SingleMode_force.fld</file>
        <file description="baseflow File">HM_Adj_SingleMode_base.bse</file>
    
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">0</value>
            <value variable="v" tolerance="1e-12">0</value>
            <value variable="w" tolerance="1e-12">0</value>
            <value variable="p" tolerance="1e-12">0</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">0</value>
            <value variable="v" tolerance="1e-12">0</value>
            <value variable="w" tolerance="1e-12">0</value>
            <value variable="p" tolerance="1e-12">0</value>
        </metric>
    </metrics>
</test>


