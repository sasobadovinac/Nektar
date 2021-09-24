<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D unsteady DG advection, prism elements, checking when faces meet as eDir1FwdDir2_Dir2BwdDir1 </description>
    <executable>ADRSolver</executable>
    <parameters>Advection3D_DG_prism.xml</parameters>
    <files>
        <file description="Session and Mesh File">Advection3D_DG_prism.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-13">0.0</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-13">0.0</value>
        </metric>
    </metrics>
</test>
