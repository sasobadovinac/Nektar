<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Embedded cylinder flow (2D)</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>EmbededCylinder2D_mesh.xml EmbededCylinder2D_session.xml</parameters>
    <files>
        <file description="Mesh File">EmbededCylinder2D_mesh.xml</file>
        <file description="Session File">EmbededCylinder2D_session.xml</file>
        <file description="Boundary Data File">EmbededCylinder2D_bnd_inflow.bc</file>
        <file description="Boundary Data File">EmbededCylinder2D_bnd_outflow.bc</file>
        <file description="Simulation Restart File">EmbededCylinder2D_restart_P2.fld</file>
        <file description="Filter Restart File and Body-fitted Coordinate File">EmbededCylinder2D_bfc.fld</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="x"              tolerance="1e-6">3.13063</value>
            <value variable="y"              tolerance="1e-6">2.39221</value>
            <value variable="rho"            tolerance="1e-6">2.63641</value>
            <value variable="rhou"           tolerance="1e-6">2.64688</value>
            <value variable="rhov"           tolerance="1e-6">0.585958</value>
            <value variable="E"              tolerance="1e-6">53.6529</value>
            <value variable="u"              tolerance="1e-6">2.68507</value>
            <value variable="v"              tolerance="1e-6">0.597729</value>
            <value variable="p"              tolerance="1e-6">20.9044</value>
            <value variable="T"              tolerance="1e-6">2.64096</value>
            <value variable="s"              tolerance="1e-10">0.000805546</value>
            <value variable="a"              tolerance="1e-6">8.80904</value>
            <value variable="Mach"           tolerance="1e-6">0.827863</value>
            <value variable="Sensor"         tolerance="1e-6">9.96263</value>
            <value variable="distanceToWall" tolerance="1e-6">2.75772</value>
            <value variable="u_bfc"          tolerance="1e-6">2.20496</value>
            <value variable="v_bfc"          tolerance="1e-6">1.64472</value>  
        </metric>
        <metric type="Linf" id="2">
            <value variable="x"              tolerance="1e-4">2</value>
            <value variable="y"              tolerance="1e-6">1.93221</value>
            <value variable="rho"            tolerance="1e-6">1.04662</value>
            <value variable="rhou"           tolerance="1e-6">1.81238</value>
            <value variable="rhov"           tolerance="1e-6">0.960331</value>
            <value variable="E"              tolerance="1e-4">21.1317</value>
            <value variable="u"              tolerance="1e-6">2.14061</value>
            <value variable="v"              tolerance="1e-6">1.00556</value>
            <value variable="p"              tolerance="1e-6">8.45278</value>
            <value variable="T"              tolerance="1e-6">1.01809</value>
            <value variable="s"              tolerance="1e-7">0.0026417</value>
            <value variable="a"              tolerance="1e-6">3.36336</value>
            <value variable="Mach"           tolerance="1e-6">0.663917</value>
            <value variable="Sensor"         tolerance="1e-6">-2.84266</value>
            <value variable="distanceToWall" tolerance="1e-5">2.3285</value>
            <value variable="u_bfc"          tolerance="1e-6">2.13868</value>
            <value variable="v_bfc"          tolerance="1e-6">0.933669</value>
        </metric>
        
    </metrics>
</test>
