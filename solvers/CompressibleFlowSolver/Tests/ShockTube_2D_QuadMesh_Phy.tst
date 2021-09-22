<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Euler 2D shocktube with Quad mesh and Physical AV</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>ShockTube_2D_QuadMesh.xml ShockTube_2D_QuadMesh_Phy.xml</parameters>
    <files>
        <file description="Mesh File">ShockTube_2D_QuadMesh.xml</file>
        <file description="Session File">ShockTube_2D_QuadMesh_Phy.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-9">4.78027e-05</value>
            <value variable="rhou" tolerance="1e-7">3.27121</value>
            <value variable="rhov" tolerance="1e-7">1.95781e-11</value>
            <value variable="E" tolerance="1e-2">1791.44</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-7">0.00195937</value>
            <value variable="rhou" tolerance="1e-6">54.6089</value>
            <value variable="rhov" tolerance="1e-6">7.89271e-10</value>
            <value variable="E" tolerance="1e-2">28228.8</value>
        </metric>
    </metrics>
</test>