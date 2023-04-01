<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>N-S 2D shocktube with mixed mesh, LaxFriedrichs Riemann solver, physical AV, dilatation sensor and interior penalty</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>ShockTube_2D_mixedMesh.xml ShockTube_2D_mixedMesh_AV_Phys_Dil_IP_LaxFriedrichs.xml</parameters>
    <files>
        <file description="Mesh File">ShockTube_2D_mixedMesh.xml</file>
        <file description="Session File">ShockTube_2D_mixedMesh_AV_Phys_Dil_IP_LaxFriedrichs.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-10">6.25963e-05</value>
            <value variable="rhou" tolerance="1e-5">1.32182</value>
            <value variable="rhov" tolerance="1e-7">0.0935274</value>
            <value variable="E" tolerance="1e-3">685.187</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-8">0.00952261</value>
            <value variable="rhou" tolerance="1e-4">32.1626</value>
            <value variable="rhov" tolerance="1e-5">8.69787</value>
            <value variable="E" tolerance="1e-1">36595.9</value>
        </metric>
    </metrics>
</test>
