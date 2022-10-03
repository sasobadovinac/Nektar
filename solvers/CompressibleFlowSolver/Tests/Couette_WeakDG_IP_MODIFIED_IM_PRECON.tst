<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>NS, Couette flow, mixed bcs, WeakDG advection and IP diffusion, Implicit</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>Couette_WeakDG_IP_MODIFIED_IM_PRECON.xml</parameters>
    <files>
        <file description="Session File">Couette_WeakDG_IP_MODIFIED_IM_PRECON.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-8">2.25256e-03</value>
            <value variable="rhou" tolerance="1e-8">2.17252e-03</value>
            <value variable="rhov" tolerance="1e-8">1.61124e-03</value>
            <value variable="E" tolerance="1e-8">5.26196e-03</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-8">2.49404e-03</value>
            <value variable="rhou" tolerance="1e-8">3.66045e-03</value>
            <value variable="rhov" tolerance="1e-8">2.10457e-03</value>
            <value variable="E" tolerance="1e-8">4.80471e-03</value>
        </metric>
    </metrics>
</test>
