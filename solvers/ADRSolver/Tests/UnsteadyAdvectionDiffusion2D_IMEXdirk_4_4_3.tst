<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>2D Unsteady Advection-Diffusion IMEXdirk (4,4,3)</description>
    <executable>ADRSolver</executable>
    <parameters>UnsteadyAdvectionDiffusion2D_IMEXdirk_4_4_3.xml</parameters>
    <files>
        <file description="Session File">UnsteadyAdvectionDiffusion2D_IMEXdirk_4_4_3.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-08">3.41659e-08</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-08">1.46835e-07</value>
        </metric>
    </metrics>
</test>
