<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>2D advection with 2 zones and 1 non-conformal curved interfaces with dirichlet BCs on quads</description>
    <executable>ADRSolver</executable>
    <parameters>Movement_fixed_2D_curved_quads.xml</parameters>
    <files>
        <file description="Session File">Movement_fixed_2D_curved_quads.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-6">1.31573e-05</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-5">9.84901e-05</value>
        </metric>
    </metrics>
</test>


