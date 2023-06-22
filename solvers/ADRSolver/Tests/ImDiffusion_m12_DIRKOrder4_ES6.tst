<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> 2D unsteady DG implicit diffusion </description>
    <executable>ADRSolver</executable>
    <parameters> ImDiffusion_m12_DIRKOrder4_ES6.xml</parameters>
    <files>
        <file description="Session File"> ImDiffusion_m12_DIRKOrder4_ES6.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-08"> 1.21277e-06 </value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-08"> 2.07721e-05 </value>
        </metric>
    </metrics>
</test>
