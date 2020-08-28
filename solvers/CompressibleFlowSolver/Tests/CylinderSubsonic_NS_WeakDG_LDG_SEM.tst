<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>NS, Subsonic Cylinder, mixed bcs, WeakDG advection and LDG diffusion, SEM</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>CylinderSubsonic_NS_WeakDG_LDG_SEM.xml</parameters>
    <files>
        <file description="Session File">CylinderSubsonic_NS_WeakDG_LDG_SEM.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-6">5.50797</value>
            <value variable="rhou" tolerance="1e-3">2016.24</value>
            <value variable="rhov" tolerance="1e-5">25.6454</value>
            <value variable="E" tolerance="1e-1">6.27025e+06</value>
            <value variable="rhoC" tolerance="1e-5">17.4286</value>
            <value variable="rhoD" tolerance="1e-5">17.4286</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-7">0.515737</value>
            <value variable="rhou" tolerance="1e-5">89.7422</value>
            <value variable="rhov" tolerance="1e-4">52.3932</value>
            <value variable="E" tolerance="1e-1">342899</value>
            <value variable="rhoC" tolerance="1e-5">1</value>
            <value variable="rhoD" tolerance="1e-5">1</value>
        </metric>
    </metrics>
</test>


