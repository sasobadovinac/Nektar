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
            <value variable="rhou" tolerance="1e-7">0.00988025</value>
            <value variable="rhov" tolerance="1e-7">0.0476664</value>
            <value variable="E" tolerance="1e-4">61.7316</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-7">0.00648549</value>
            <value variable="rhou" tolerance="1e-6">0.0120572</value>
            <value variable="rhov" tolerance="1e-7">0.0450055</value>
            <value variable="E" tolerance="1e-4">52.0678</value>
        </metric>
    </metrics>
</test>
