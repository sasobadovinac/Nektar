<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>2D advection with 3 zones and 2 non-conformal interfaces with periodic BCs on 4 processes with central zone translating</description>
    <executable>ADRSolver</executable>
    <parameters>Movement_translate_interfaces232_dirichlet_par.xml</parameters>
    <processes>4</processes>
    <files>
        <file description="Session File">Movement_translate_interfaces232_dirichlet_par.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-10">1.71711e-11</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-10">2.30174e-10</value>
        </metric>
    </metrics>
</test>


