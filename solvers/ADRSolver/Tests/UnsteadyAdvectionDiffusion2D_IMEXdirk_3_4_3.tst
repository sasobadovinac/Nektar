<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>2D Unsteady Advection-Diffusion IMEXdirk (3,4,3)</description>
    <executable>ADRSolver</executable>
    <parameters>UnsteadyAdvectionDiffusion2D_IMEXdirk_3_4_3.xml</parameters>
    <files>
        <file description="Session File">UnsteadyAdvectionDiffusion2D_IMEXdirk_3_4_3.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-08">4.10992e-08</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-08">2.56861e-07</value>
        </metric>
    </metrics>
</test>
