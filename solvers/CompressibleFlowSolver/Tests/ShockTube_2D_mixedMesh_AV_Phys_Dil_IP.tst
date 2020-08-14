<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>N-S 2D shocktube with mixed mesh, physical AV, dilatation sensor and interior penalty</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>ShockTube_2D_mixedMesh.xml ShockTube_2D_mixedMesh_AV_Phys_Dil_IP.xml</parameters>
    <files>
        <file description="Mesh File">ShockTube_2D_mixedMesh.xml</file>
        <file description="Session File">ShockTube_2D_mixedMesh_AV_Phys_Dil_IP.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-9">1.59707e-05</value>
            <value variable="rhou" tolerance="1e-5">1.29832</value>
            <value variable="rhov" tolerance="1e-6">0.146938</value>
            <value variable="E" tolerance="1e-2">657.021</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-7">0.0016932</value>
            <value variable="rhou" tolerance="1e-4">33.4157</value>
            <value variable="rhov" tolerance="1e-4">16.3785</value>
            <value variable="E" tolerance="1e-1">31917.8</value>
        </metric>
    </metrics>
</test>