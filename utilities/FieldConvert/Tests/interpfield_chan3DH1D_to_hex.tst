<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Interpolate field from chan3DH1D to Hex_channel </description>
    <executable>FieldConvert</executable>
    <parameters> -f -e -m interpfield:fromxml=chan3DH1D.xml:fromfld=chan3DH1D.fld Hex_channel_C0helmsmoothing.xml chan.plt
    </parameters>
    <files>
        <file description="From Mesh File 3DH1D">chan3DH1D.xml</file>
        <file description="From Field File">chan3DH1D.fld</file>
        <file description="To Mesh File Hex">Hex_channel_C0helmsmoothing.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="x" tolerance="1e-6">0.57735</value>
            <value variable="y" tolerance="1e-6">0.57735</value>
            <value variable="z" tolerance="1e-6">0.57735</value>
            <value variable="u" tolerance="1e-6">0.183207</value>
            <value variable="v" tolerance="1e-12">2.47059e-16</value>
            <value variable="w" tolerance="1e-6">0</value>
            <value variable="p" tolerance="1e-6">1.1547</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="x" tolerance="1e-6">1</value>
            <value variable="y" tolerance="1e-6">1</value>
            <value variable="z" tolerance="1e-6">1</value>
            <value variable="u" tolerance="1e-6">0.25</value>
            <value variable="v" tolerance="1e-12">9.70661e-16</value>
            <value variable="w" tolerance="1e-12">8.87226e-18</value>
            <value variable="p" tolerance="1e-6">2</value>
        </metric>
    </metrics>
</test>