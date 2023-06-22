<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Laminar Channel Flow 3D homogeneous 1D, p-Refinenement tag</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>ChanFlow_3DH1D_pRef.xml</parameters>
    <files>
        <file description="Session File">ChanFlow_3DH1D_pRef.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-8">0.00617619</value>
            <value variable="v" tolerance="1e-8">0.000200183</value>
            <value variable="w" tolerance="1e-12">0</value>
            <value variable="p" tolerance="1e-7">0.0101901</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-7">0.0165127</value>
            <value variable="v" tolerance="1e-8">0.000954683</value>
            <value variable="w" tolerance="1e-12">0</value>
            <value variable="p" tolerance="1e-7">0.0575085</value>
        </metric>
    </metrics>
</test>


