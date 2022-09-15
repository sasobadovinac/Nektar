<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Kovasznay Flow P=3 with Block-Preconditioned IterativeFull matrix solve</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>-v -I GlobalSysSoln=IterativeFull -I Preconditioner=Block -I IterativeSolverTolerance=1e-14 -I LinSysIterSolver=ConjugateGradient KovaFlow_m3.xml</parameters>
    <files>
        <file description="Session File">KovaFlow_m3.xml</file>
        <file description="Session File">KovaFlow_m3.rst</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">0.343606</value>
            <value variable="v" tolerance="1e-12">0.0255658</value>
            <value variable="p" tolerance="1e-12">0.347747</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">0.322039</value>
            <value variable="v" tolerance="1e-12">0.0243043</value>
            <value variable="p" tolerance="1e-12">0.519425</value>
        </metric>
        <metric type="Precon" id="3">
            <value variable="p" tolerance="2">31</value>
            <value variable="u" tolerance="2">58</value>
            <value variable="v" tolerance="2">39</value>
        </metric>
    </metrics>
</test>
