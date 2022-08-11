<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>2D unsteady advection on quad mesh using basis:Modified_A, P=4, Q=5  Periodic BCs, structure preserving filter on</description>
    <executable>ADRSolver</executable>
    <parameters>Advection2D_structure_preserving_filter.xml</parameters>
    <files>
        <file description="Session File">Advection2D_structure_preserving_filter.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-10"> 0.0872027 </value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-10"> 0.771972 </value>
        </metric>
    </metrics>
</test>
