<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Euler, Subsonic Cylinder, WeakDG advection, Implicit solver</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>CylinderSubsonic_WeakDG_Implicit.xml</parameters>
    <files>
        <file description="Session File">CylinderSubsonic_WeakDG_Implicit.xml</file>
    </files>
    <metrics>  
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-4">10.7744</value>
            <value variable="rhou" tolerance="1e-5">1.09375</value>
            <value variable="rhov" tolerance="2e-7">0.0768964</value>
            <value variable="E" tolerance="1e-12">2.228e+06</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">1.22501</value>
            <value variable="rhou" tolerance="2e-5">0.246677</value>
            <value variable="rhov" tolerance="1e-5">0.134325</value>
            <value variable="E" tolerance="1e-12">253314</value>
        </metric>
    </metrics>
</test>


