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
            <value variable="u" tolerance="1e-9">9.25579e-10</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-9">2.06957e-08</value>
        </metric>
    </metrics>
</test>


