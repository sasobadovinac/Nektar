<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>NS, Couette flow, mixed bcs, WeakDG advection and IP diffusion, Implicit parallel</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>Couette_WeakDG_IP_MODIFIED_IM.xml</parameters>
    <processes>2</processes>
    <files>
        <file description="Session File">Couette_WeakDG_IP_MODIFIED_IM.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-9">0.000352778</value>
            <value variable="rhou" tolerance="1e-6">0.180406</value>
            <value variable="rhov" tolerance="1e-7">0.0348969</value>
            <value variable="E" tolerance="1e-4">32.4753</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-9">0.000271061</value>
            <value variable="rhou" tolerance="1e-6">0.154715</value>
            <value variable="rhov" tolerance="1e-7">0.0443422</value>
            <value variable="E" tolerance="1e-3">33.967</value>
        </metric>
    </metrics>
</test>


