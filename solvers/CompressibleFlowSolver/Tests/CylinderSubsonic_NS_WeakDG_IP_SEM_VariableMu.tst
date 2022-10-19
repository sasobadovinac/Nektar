<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>NS, Subsonic Cylinder, mixed bcs, WeakDG advection and IP diffusion, variable viscosity</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>CylinderSubsonic_NS_WeakDG_IP_SEM_VariableMu.xml</parameters>
    <files>
        <file description="Session File">CylinderSubsonic_NS_WeakDG_IP_SEM_VariableMu.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-6">5.50084</value>
            <value variable="rhou" tolerance="1e-6">2016.91</value>
            <value variable="rhov" tolerance="1e-4">73.9407</value>
            <value variable="E" tolerance="1e-1">6.26943e+06</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-6">0.249876</value>
            <value variable="rhou" tolerance="1e-4">110.272</value>
            <value variable="rhov" tolerance="1e-4">59.0975</value>
            <value variable="E" tolerance="1e-1">261932</value>
        </metric>
    </metrics>
</test>


