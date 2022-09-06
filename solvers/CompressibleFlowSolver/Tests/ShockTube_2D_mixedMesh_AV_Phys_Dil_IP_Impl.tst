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
            <value variable="rho" tolerance="1e-8">0.00121093</value>
            <value variable="rhou" tolerance="1e-5">1.34602</value>
            <value variable="rhov" tolerance="1e-7">0.0671994</value>
            <value variable="E" tolerance="1e-2">3398.66</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-7">0.0733791</value>
            <value variable="rhou" tolerance="1e-4">27.7611</value>
            <value variable="rhov" tolerance="1e-5">9.42978</value>
            <value variable="E" tolerance="1e-0">212621</value>
        </metric>
    </metrics>
</test>