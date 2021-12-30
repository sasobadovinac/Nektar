<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>NS, Couette flow, mixed bcs, WeakDG advection and IP diffusion, Implicit, Adiabatic wall</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>Couette_WeakDG_IP_MODIFIED_Adiabatic.xml</parameters>
    <files>
        <file description="Session File">Couette_WeakDG_IP_MODIFIED_Adiabatic.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-6">0.001587</value>
            <value variable="rhou" tolerance="1e-7">0.00982231</value>
            <value variable="rhov" tolerance="1e-7">0.0476072</value>
            <value variable="E" tolerance="1e-4">62.4086</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-7">0.00648099</value>
            <value variable="rhou" tolerance="1e-6">0.0121131</value>
            <value variable="rhov" tolerance="1e-7">0.0438114</value>
            <value variable="E" tolerance="1e-4">51.3568</value>
        </metric>
    </metrics>
</test>
