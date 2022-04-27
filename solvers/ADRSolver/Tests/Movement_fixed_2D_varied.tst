<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>2D advection with 3 zones and 2 non-conformal interfaces (one curved) with 2 fields</description>
    <executable>ADRSolver</executable>
    <parameters>Movement_fixed_2D_varied.xml</parameters>
    <files>
        <file description="Session File">Movement_fixed_2D_varied.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-6">3.98672e-06</value>
            <value variable="v" tolerance="1e-6">2.70806e-06</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-5">6.90779e-05</value>
            <value variable="v" tolerance="1e-5">1.78188e-05</value>
        </metric>
    </metrics>
</test>


