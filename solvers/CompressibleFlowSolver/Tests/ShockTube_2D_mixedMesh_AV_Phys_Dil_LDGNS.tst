<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>N-S 2D shocktube with mixed mesh, physical AV, dilatation sensor and LDGNS</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>ShockTube_2D_mixedMesh.xml ShockTube_2D_mixedMesh_AV_Phys_Dil_LDGNS.xml</parameters>
    <files>
        <file description="Mesh File">ShockTube_2D_mixedMesh.xml</file>
        <file description="Session File">ShockTube_2D_mixedMesh_AV_Phys_Dil_LDGNS.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-8">0.00121163</value>
            <value variable="rhou" tolerance="1e-5">1.50394</value>
            <value variable="rhov" tolerance="1e-6">0.170155</value>
            <value variable="E" tolerance="1e-2">878.249</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-7">0.078705</value>
            <value variable="rhou" tolerance="1e-4">43.735</value>
            <value variable="rhov" tolerance="1e-4">23.2031</value>
            <value variable="E" tolerance="1e-1">71223.8</value>
        </metric>
    </metrics>
</test>