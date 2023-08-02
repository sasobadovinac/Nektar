<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> 1D unsteady DG explicit advection</description>
    <executable>ADRSolver</executable>
    <parameters>--npt 8 PFASSTDriverUnsteadyAdvection1D.xml</parameters>
    <processes>8</processes>
    <files>
        <file description="Session File"> PFASSTDriverUnsteadyAdvection1D.xml </file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="5e-11">2.06495e-06</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="5e-10">2.71381e-05</value>
        </metric>
    </metrics>
</test>
