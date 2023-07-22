<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> 2D unsteady DG implicit diffusion, order 1, P=12 </description>
    <executable>ADRSolver</executable>
    <parameters>--npt 10 PararealDriverImDiffusion_m12_iter4.xml</parameters>
    <processes>10</processes>
    <files>
        <file description="Session File"> PararealDriverImDiffusion_m12_iter4.xml </file>
        <file description="Session File"> PararealDriverImDiffusion_m12_iter4_coarseSolver.xml </file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="5e-10">7.70078e-05</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="5e-10">7.85498e-05</value>
        </metric>
    </metrics>
</test>
