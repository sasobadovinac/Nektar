<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>N-S 2D shocktube with mixed mesh, AUSM2 Riemann solver, physical AV, dilatation sensor and interior penalty</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>ShockTube_2D_mixedMesh.xml ShockTube_2D_mixedMesh_AV_Phys_Dil_IP_AUSM2.xml</parameters>
    <files>
        <file description="Mesh File">ShockTube_2D_mixedMesh.xml</file>
        <file description="Session File">ShockTube_2D_mixedMesh_AV_Phys_Dil_IP_AUSM2.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-10">1.99286e-05</value>
            <value variable="rhou" tolerance="1e-5">1.33882</value>
            <value variable="rhov" tolerance="1e-7">0.097965</value>
            <value variable="E" tolerance="1e-3">661.684</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-8">0.00259851</value>
            <value variable="rhou" tolerance="1e-4">35.1446</value>
            <value variable="rhov" tolerance="1e-5">8.60605</value>
            <value variable="E" tolerance="1e-1">30638.7</value>
        </metric>
    </metrics>
</test>
