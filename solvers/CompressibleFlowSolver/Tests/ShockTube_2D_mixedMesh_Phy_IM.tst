<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Euler 2D shocktube with mixed mesh and Physical AV using Implicit CFS solver</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>ShockTube_2D_mixedMesh.xml ShockTube_2D_mixedMesh_Phy_IM.xml</parameters>
    <files>
        <file description="Mesh File">ShockTube_2D_mixedMesh.xml</file>
        <file description="Session File">ShockTube_2D_mixedMesh_Phy_IM.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-9">0.0213636</value>
            <value variable="rhou" tolerance="1e-7">61.7548</value>
            <value variable="rhov" tolerance="1e-7">3.76726</value>
            <value variable="E" tolerance="1e-2">80698.8</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-7">0.524274</value>
            <value variable="rhou" tolerance="1e-6">931.907</value>
            <value variable="rhov" tolerance="1e-6">128.928</value>
            <value variable="E" tolerance="1e-2">1.76616e+06</value>
        </metric>
    </metrics>
</test>