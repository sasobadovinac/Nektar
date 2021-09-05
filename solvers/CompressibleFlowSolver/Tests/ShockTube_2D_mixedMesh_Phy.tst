<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Euler 2D shocktube with mixed mesh and Physical AV</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>ShockTube_2D_mixedMesh.xml ShockTube_2D_mixedMesh_Phy.xml</parameters>
    <files>
        <file description="Mesh File">ShockTube_2D_mixedMesh.xml</file>
        <file description="Session File">ShockTube_2D_mixedMesh_Phy.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-9">0.000164192</value>
            <value variable="rhou" tolerance="1e-7">6.2884</value>
            <value variable="rhov" tolerance="1e-7">0.90211</value>
            <value variable="E" tolerance="1e-2">3185.61</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-7">0.00829796</value>
            <value variable="rhou" tolerance="1e-6">165.971</value>
            <value variable="rhov" tolerance="1e-6">41.7404</value>
            <value variable="E" tolerance="1e-2">92663.6</value>
        </metric>
    </metrics>
</test>