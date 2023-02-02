<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>2D advection with 3 zones and 2 non-conformal interfaces with periodic BCs on 4 processes</description>
    <executable>ADRSolver</executable>
    <parameters>Movement_fixed_interfaces232_periodic.xml</parameters>
    <processes>4</processes>
    <files>
        <file description="Session File">Movement_fixed_interfaces232_periodic.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-8">7.03318e-10</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-8">1.32932e-09</value>
        </metric>
    </metrics>
</test>


