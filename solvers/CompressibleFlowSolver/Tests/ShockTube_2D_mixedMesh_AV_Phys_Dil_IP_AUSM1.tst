<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>N-S 2D shocktube with mixed mesh, AUSM1 Riemann solver, physical AV, dilatation sensor and interior penalty</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>ShockTube_2D_mixedMesh.xml ShockTube_2D_mixedMesh_AV_Phys_Dil_IP_AUSM1.xml</parameters>
    <files>
        <file description="Mesh File">ShockTube_2D_mixedMesh.xml</file>
        <file description="Session File">ShockTube_2D_mixedMesh_AV_Phys_Dil_IP_AUSM1.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-10">1.95652e-05</value>
            <value variable="rhou" tolerance="1e-5">1.33845</value>
            <value variable="rhov" tolerance="1e-7">0.0981891</value>
            <value variable="E" tolerance="1e-3">661.044</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-8">0.00273414</value>
            <value variable="rhou" tolerance="1e-4">34.9565</value>
            <value variable="rhov" tolerance="1e-5">8.66406</value>
            <value variable="E" tolerance="1e-1">28724.5</value>
        </metric>
    </metrics>
</test>
