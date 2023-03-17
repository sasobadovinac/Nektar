<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>2D Unsteady Advection-Diffusion IMEXdirk (2,3,3)</description>
    <executable>ADRSolver</executable>
    <parameters>UnsteadyAdvectionDiffusion2D_IMEXdirk_2_3_3.xml</parameters>
    <files>
        <file description="Session File">UnsteadyAdvectionDiffusion2D_IMEXdirk_2_3_3.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-08">1.85112e-07</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-08">7.82281e-07</value>
        </metric>
    </metrics>
</test>
