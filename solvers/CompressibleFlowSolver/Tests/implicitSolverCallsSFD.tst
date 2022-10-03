<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Test the selective frequency damping function in the implicit solver</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>implicitSolverCallsSFD_mesh.xml implicitSolverCallsSFD_session.xml</parameters>
    <files>
        <file description="Mesh File">implicitSolverCallsSFD_mesh.xml</file>
        <file description="Session File">implicitSolverCallsSFD_session.xml</file>
        <file description="Restart File">implicitSolverCallsSFD.fld</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho"  tolerance="1e-8">0.00707092</value>
            <value variable="rhou" tolerance="1e-8">0.00657595</value>
            <value variable="rhov" tolerance="1e-8">1.43331e-06</value>
            <value variable="E"    tolerance="1e-3">0.0535753</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho"  tolerance="1e-8">1.00018</value>
            <value variable="rhou" tolerance="1e-8">1.00178</value>
            <value variable="rhov" tolerance="1e-8">0.000987263</value>
            <value variable="E"    tolerance="1e-3">7.64495</value>
        </metric>
    </metrics>
</test>
