<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>3D flexible cylinder flow simulation using homogeneous strip modeling</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>CylFlow_HomoStrip.xml</parameters>
    <files>
        <file description="Session File">CylFlow_HomoStrip.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">9.81073</value>
            <value variable="v" tolerance="1e-12">0.660518</value>
            <value variable="w" tolerance="1e-12">0.0</value>
            <value variable="p" tolerance="1e-12">29.4734</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">1.44442</value>
            <value variable="v" tolerance="1e-12">0.660518</value>
            <value variable="w" tolerance="1e-12">0.0</value>
            <value variable="p" tolerance="1e-12">5.61629</value>
        </metric>
    </metrics>
</test>
