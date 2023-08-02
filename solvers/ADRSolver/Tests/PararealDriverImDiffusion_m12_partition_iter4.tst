<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> 2D unsteady DG implicit diffusion, order 1, P=12 </description>
    <executable>ADRSolver</executable>
    <parameters>--npt 5 PararealDriverImDiffusion_m12_iter4_xml</parameters>
    <processes>10</processes>
    <files>
        <file description="Session File"> PararealDriverImDiffusion_m12_iter4_xml </file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="5e-10">7.5238e-05</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="5e-10">7.53325e-05</value>
        </metric>
    </metrics>
</test>
