<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D unsteady DG advection, all elmts P=4 with Gradient Stabilisation </description>
    <executable>ADRSolver</executable>
    <parameters>Advection3D_AllElmt_GJP.xml</parameters>
    <files>
        <file description="Session File">Advection3D_AllElmt_GJP.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-7">0.00224251</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-7">0.0354197</value>
        </metric>
    </metrics>
</test>

