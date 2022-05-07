<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>2D advection with 1 zones and 1 non-conformal interfaces consisting of two cylinders stacked vertically with dirichlet BCs on 4 processes</description>
    <executable>ADRSolver</executable>
    <parameters>Movement_fixed_3D_stacked_cylinders_curved.xml</parameters>
    <processes>4</processes>
    <files>
        <file description="Session File">Movement_fixed_3D_stacked_cylinders_curved.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-6">2.10558e-06</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-5">4.21843e-05</value>
        </metric>
    </metrics>
</test>


