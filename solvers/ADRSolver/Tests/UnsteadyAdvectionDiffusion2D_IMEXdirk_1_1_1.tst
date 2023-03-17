<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>2D Unsteady Advection-Diffusion IMEXdirk (1,1,1)</description>
    <executable>ADRSolver</executable>
    <parameters>UnsteadyAdvectionDiffusion2D_IMEXdirk_1_1_1.xml</parameters>
    <files>
        <file description="Session File">UnsteadyAdvectionDiffusion2D_IMEXdirk_1_1_1.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-08">0.00150483</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-08">0.00169437</value>
        </metric>
    </metrics>
</test>
