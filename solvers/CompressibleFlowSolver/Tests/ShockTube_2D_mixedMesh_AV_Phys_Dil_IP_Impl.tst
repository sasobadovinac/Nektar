<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>
        N-S 2D shocktube with mixed mesh, physical AV, dilatation sensor and interior penalty

        This test case is identical to ShockTube_2D_mixedMesh_AV_Phys_Dil_IP, except that it
        uses the implicit Navier Stokes solver.
    </description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>ShockTube_2D_mixedMesh.xml ShockTube_2D_mixedMesh_AV_Phys_Dil_IP_Impl.xml</parameters>
    <files>
        <file description="Mesh File">ShockTube_2D_mixedMesh.xml</file>
        <file description="Session File">ShockTube_2D_mixedMesh_AV_Phys_Dil_IP_Impl.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-8">0.00121095</value>
            <value variable="rhou" tolerance="5e-5">1.32668</value>
            <value variable="rhov" tolerance="2e-6">0.0617728</value>
            <value variable="E" tolerance="2e-2">3398.85</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="4e-7">0.0733394</value>
            <value variable="rhou" tolerance="1e-3">28.0426</value>
            <value variable="rhov" tolerance="3e-4">7.09624</value>
            <value variable="E" tolerance="1e-0">212302</value>
        </metric>
    </metrics>
</test>
