<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>N-S 2D shocktube with mixed mesh, physical AV, modal sensor and LDGNS</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>ShockTube_2D_mixedMesh.xml ShockTube_2D_mixedMesh_AV_Phys_Mod_LDGNS.xml</parameters>
    <files>
        <file description="Mesh File">ShockTube_2D_mixedMesh.xml</file>
        <file description="Session File">ShockTube_2D_mixedMesh_AV_Phys_Mod_LDGNS.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-8">0.00121235</value>
            <value variable="rhou" tolerance="1e-5">2.68086</value>
            <value variable="rhov" tolerance="1e-6">0.170225</value>
            <value variable="E" tolerance="1e-2">1494.76</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-7">0.0816055</value>
            <value variable="rhou" tolerance="1e-4">74.3684</value>
            <value variable="rhov" tolerance="1e-4">22.4862</value>
            <value variable="E" tolerance="1e-1">59111.8</value>
        </metric>
    </metrics>
</test>