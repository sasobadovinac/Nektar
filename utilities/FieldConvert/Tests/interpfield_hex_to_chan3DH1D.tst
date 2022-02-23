<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Interpolate field from Hex_channel to chan3DH1D </description>
    <executable>FieldConvert</executable>
    <parameters> -f -e -m interpfield:fromxml=Hex_channel_C0helmsmoothing.xml:fromfld=Hex_channel_C0helmsmoothing.fld chan3DH1D.xml chan.plt
    </parameters>
    <files>
        <file description="From Mesh File Hex">Hex_channel_C0helmsmoothing.xml</file>
        <file description="From Field File Hex">Hex_channel_C0helmsmoothing.fld</file>
        <file description="To Mesh File 3DH1D">chan3DH1D.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="x" tolerance="1e-6">0.408248</value>
            <value variable="y" tolerance="1e-6">0.408248</value>
            <value variable="z" tolerance="1e-6">0.369755</value>
            <value variable="u" tolerance="1e-6">0.662361</value>
            <value variable="v" tolerance="1e-6">0.662361</value>
            <value variable="w" tolerance="1e-6">0.62619</value>
            <value variable="p" tolerance="1e-6">1.12748</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="x" tolerance="1e-6">1</value>
            <value variable="y" tolerance="1e-6">1</value>
            <value variable="z" tolerance="1e-6">0.875</value>
            <value variable="u" tolerance="1e-6">1.16466</value>
            <value variable="v" tolerance="1e-6">1.16466</value>
            <value variable="w" tolerance="1e-6">1.03862</value>
            <value variable="p" tolerance="1e-6">3.36795</value>
        </metric>
    </metrics>
</test>