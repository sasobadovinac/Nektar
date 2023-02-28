<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D advection with 3 zones and 2 non-conformal interfaces with periodic BCs on hexes on 4 processes</description>
    <executable>ADRSolver</executable>
    <parameters>Movement_fixed_3D_Hex_3zones_Periodic_par.xml</parameters>
    <processes>4</processes>
    <files>
        <file description="Session File">Movement_fixed_3D_Hex_3zones_Periodic_par.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-6">1.58599e-06</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-5">3.03816e-05</value>
        </metric>
    </metrics>
</test>


