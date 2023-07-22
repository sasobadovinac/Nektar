<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>N-S 2D shocktube with mixed mesh, HLLC Riemann solver, physical AV, dilatation sensor and interior penalty</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>ShockTube_2D_mixedMesh.xml ShockTube_2D_mixedMesh_AV_Phys_Dil_IP_HLLC.xml</parameters>
    <files>
        <file description="Mesh File">ShockTube_2D_mixedMesh.xml</file>
        <file description="Session File">ShockTube_2D_mixedMesh_AV_Phys_Dil_IP_HLLC.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-10">2.111e-05</value>
            <value variable="rhou" tolerance="1e-5">1.33915</value>
            <value variable="rhov" tolerance="1e-7">0.0983775</value>
            <value variable="E" tolerance="1e-3">662.788</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-8">0.00249791</value>
            <value variable="rhou" tolerance="1e-4">35.3069</value>
            <value variable="rhov" tolerance="1e-5">8.69836</value>
            <value variable="E" tolerance="1e-1">32353.9</value>
        </metric>
    </metrics>
</test>
