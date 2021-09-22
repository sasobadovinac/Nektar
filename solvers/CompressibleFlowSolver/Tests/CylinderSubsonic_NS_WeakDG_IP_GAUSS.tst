<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>NS, Subsonic Cylinder, mixed bcs, WeakDG advection and LDG diffusion, SEM</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>CylinderSubsonic_NS_WeakDG_IP_GAUSS.xml</parameters>
    <files>
        <file description="Session File">CylinderSubsonic_NS_WeakDG_IP_GAUSS.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">5.50779</value>
            <value variable="rhou" tolerance="1e-12">2016.32</value>
            <value variable="rhov" tolerance="5e-4">24.7023</value>
            <value variable="E" tolerance="1e-12">6.27023e+06</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">0.32958</value>
            <value variable="rhou" tolerance="2e-3">88.9862</value>
            <value variable="rhov" tolerance="1e-8">18.4571</value>
            <value variable="E" tolerance="1e-12">287048</value>
        </metric>
    </metrics>
</test>


