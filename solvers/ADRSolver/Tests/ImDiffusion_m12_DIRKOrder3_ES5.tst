<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> 2D unsteady DG implicit diffusion </description>
    <executable>ADRSolver</executable>
    <parameters> ImDiffusion_m12_DIRKOrder3_ES5.xml</parameters>
    <files>
        <file description="Session File"> ImDiffusion_m12_DIRKOrder3_ES5.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-10"> 1.60552e-08 </value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-10"> 3.73785e-07 </value>
        </metric>
    </metrics>
</test>
