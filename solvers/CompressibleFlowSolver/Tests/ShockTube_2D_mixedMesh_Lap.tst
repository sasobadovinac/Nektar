<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Euler 2D shocktube with mixed mesh and Laplacian AV</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>ShockTube_2D_mixedMesh.xml ShockTube_2D_mixedMesh_Lap.xml</parameters>
    <files>
        <file description="Mesh File">ShockTube_2D_mixedMesh.xml</file>
        <file description="Session File">ShockTube_2D_mixedMesh_Lap.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-9">0.00577162</value>
            <value variable="rhou" tolerance="1e-7">2.08174</value>
            <value variable="rhov" tolerance="1e-7">0.0251614</value>
            <value variable="E" tolerance="1e-2">15546</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-7">0.117437</value>
            <value variable="rhou" tolerance="1e-6">36.2498</value>
            <value variable="rhov" tolerance="1e-6">2.45249</value>
            <value variable="E" tolerance="1e-2">311418</value>
        </metric>
    </metrics>
</test>