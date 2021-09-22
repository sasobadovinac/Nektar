<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Euler 2D shocktube with Quad mesh and Physical AV and Implicit solver</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>ShockTube_2D_QuadMesh.xml ShockTube_2D_QuadMesh_PhyIm.xml</parameters>
    <files>
        <file description="Mesh File">ShockTube_2D_QuadMesh.xml</file>
        <file description="Session File">ShockTube_2D_QuadMesh_PhyIm.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-9">0.000565113</value>
            <value variable="rhou" tolerance="1e-7">1.92406</value>
            <value variable="rhov" tolerance="1e-7">3.08407e-08</value>
            <value variable="E" tolerance="1e-2">1823.07</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-7">0.0111643</value>
            <value variable="rhou" tolerance="1e-6">32.8239</value>
            <value variable="rhov" tolerance="1e-6">4.19949e-06</value>
            <value variable="E" tolerance="1e-2">40356.8</value>
        </metric>
    </metrics>
</test>