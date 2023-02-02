<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D advection with 2 zones and 1 non-conformal interfaces with dirichlet BCs on hexes</description>
    <executable>ADRSolver</executable>
    <parameters>Movement_fixed_3D_Hex41.xml</parameters>
    <files>
        <file description="Session File">Movement_fixed_3D_Hex41.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-6">5.4459e-06</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-5">3.65216e-05</value>
        </metric>
    </metrics>
</test>


