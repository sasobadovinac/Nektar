<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>2D advection with 3 zones and 2 non-conformal interfaces with dirichlet BCs</description>
    <executable>ADRSolver</executable>
    <parameters>Movement_fixed_interfaces232_dirichlet.xml</parameters>
    <files>
        <file description="Session File">Movement_fixed_interfaces232_dirichlet.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-7">3.82585e-08</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-7">4.14903e-06</value>
        </metric>
    </metrics>
</test>


