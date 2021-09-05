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
            <value variable="rho" tolerance="1e-7">0.281714</value>
            <value variable="rhou" tolerance="1e-6">26.2196</value>
            <value variable="rhov" tolerance="1e-8">1.00695</value>
            <value variable="E" tolerance="1e-12">1323.78</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">0.251024</value>
            <value variable="rhou" tolerance="1e-12">31.4528</value>
            <value variable="rhov" tolerance="1e-8">1.49992</value>
            <value variable="E" tolerance="1e-12">1920.9</value>
        </metric>
    </metrics>
</test>
