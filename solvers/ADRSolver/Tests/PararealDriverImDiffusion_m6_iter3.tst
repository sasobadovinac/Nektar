<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> 2D unsteady DG implicit diffusion, order 1, P=6 </description>
    <executable>ADRSolver</executable>
    <parameters>--npt 16 PararealDriverImDiffusion_m6_iter3.xml</parameters>
    <processes>16</processes>
    <files>
        <file description="Session File"> PararealDriverImDiffusion_m6_iter3.xml </file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-08"> 0.000226594 </value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-08"> 0.000582864 </value>
        </metric>
    </metrics>
</test>
