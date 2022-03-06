<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>2D unsteady DG advection, tri elmts P=4 with Gradient Stabilisation </description>
    <executable>ADRSolver</executable>
    <parameters>Advection2D_Tri_GJP.xml</parameters>
    <files>
        <file description="Session File">Advection2D_Tri_GJP.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-6">0.0305673</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-6">0.132392</value>
        </metric>
    </metrics>
</test>

