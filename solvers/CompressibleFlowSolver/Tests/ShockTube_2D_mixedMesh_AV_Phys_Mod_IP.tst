<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>N-S 2D shocktube with mixed mesh, physical AV, modal sensor and interior penalty</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>ShockTube_2D_mixedMesh.xml ShockTube_2D_mixedMesh_AV_Phys_Mod_IP.xml</parameters>
    <files>
        <file description="Mesh File">ShockTube_2D_mixedMesh.xml</file>
        <file description="Session File">ShockTube_2D_mixedMesh_AV_Phys_Mod_IP.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-9">1.95874e-05</value>
            <value variable="rhou" tolerance="1e-5">2.40445</value>
            <value variable="rhov" tolerance="1e-6">0.255749</value>
            <value variable="E" tolerance="1e-1">1264.6</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-7">0.00170242</value>
            <value variable="rhou" tolerance="1e-4">55.1967</value>
            <value variable="rhov" tolerance="1e-4">31.0243</value>
            <value variable="E" tolerance="1e-1">24283.2</value>
        </metric>
    </metrics>
</test>