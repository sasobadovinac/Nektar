<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>NS, Couette flow, mixed bcs, WeakDG advection and IP diffusion, Implicit</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>Couette_WeakDG_IP_MODIFIED_IM_PRECON.xml</parameters>
    <processes>4</processes>
    <files>
        <file description="Session File">Couette_WeakDG_IP_MODIFIED_IM_PRECON.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="2e-9">3.52778e-04</value>
            <value variable="rhou" tolerance="2e-6">1.80406e-01</value>
            <value variable="rhov" tolerance="2e-7">3.48969e-02</value>
            <value variable="E" tolerance="5e-4">3.24753e+01</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="5e-9">2.71061e-04</value>
            <value variable="rhou" tolerance="2e-7">1.54715e-01</value>
            <value variable="rhov" tolerance="1e-9">4.43422e-02</value>
            <value variable="E" tolerance="5e-2">3.39670e+01</value>
        </metric>
    </metrics>
</test>
