<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>
        Navier Stokes 2D shocktube with mixed mesh and Laplacian AV.

        This test is meant to verify that the NonSmooth shock capture implementation
        works correctly for the explicit Navier Stokes solver. Therefore, this test
        is identical to the ShockTube_2D_mixedMesh_AV_Lap test case, except that we 
        are using the explicit Navier Stokes solver instead of the explicit Euler solver. 
        To make sure both tests give the same results, the viscosity and the thermal 
        conductivity are set to 0.
    </description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>ShockTube_2D_mixedMesh.xml ShockTube_2D_mixedMesh_AV_Lap_NS.xml</parameters>
    <files>
        <file description="Mesh File">ShockTube_2D_mixedMesh.xml</file>
        <file description="Session File">ShockTube_2D_mixedMesh_AV_Lap_NS.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-9">0.000395188</value>
            <value variable="rhou" tolerance="1e-7">0.110602</value>
            <value variable="rhov" tolerance="1e-7">0.00309297</value>
            <value variable="E" tolerance="1e-2">1077.22</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-7">0.0196245</value>
            <value variable="rhou" tolerance="1e-6">2.70304</value>
            <value variable="rhov" tolerance="1e-6">0.406501</value>
            <value variable="E" tolerance="1e-2">55681.2</value>
        </metric>
    </metrics>
</test>