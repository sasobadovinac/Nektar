<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>2D Unsteady Advection-Diffusion IMEXdirk (2,2,2)</description>
    <executable>ADRSolver</executable>
    <parameters>UnsteadyAdvectionDiffusion2D_IMEXdirk_2_2_2.xml</parameters>
    <files>
        <file description="Session File">UnsteadyAdvectionDiffusion2D_IMEXdirk_2_2_2.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-08">2.29902e-06</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-08">2.28871e-06</value>
        </metric>
    </metrics>
</test>
