<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Interpolate field from chan3DH1D to chan3DH1D </description>
    <executable>FieldConvert</executable>
    <parameters> -f -e -m interpfield:fromxml=chan3DH1D.xml:fromfld=chan3DH1D.fld chan3DH1D.xml chan.plt
    </parameters>
    <files>
        <file description="From and To Mesh File 3DH1D">chan3DH1D.xml</file>
        <file description="From Field File">chan3DH1D.fld</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="x" tolerance="1e-6">0.408248</value>
            <value variable="y" tolerance="1e-6">0.408248</value>
            <value variable="z" tolerance="1e-6">0.369755</value>
            <value variable="u" tolerance="1e-6">0.129547</value>
            <value variable="v" tolerance="1e-12">1.98701e-16</value>
            <value variable="w" tolerance="1e-6">0</value>
            <value variable="p" tolerance="1e-6">0.816497</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="x" tolerance="1e-6">1</value>
            <value variable="y" tolerance="1e-6">1</value>
            <value variable="z" tolerance="1e-6">0.875</value>
            <value variable="u" tolerance="1e-6">0.25</value>
            <value variable="v" tolerance="1e-12">9.70677e-16</value>
            <value variable="w" tolerance="1e-12">1.26243e-17</value>
            <value variable="p" tolerance="1e-6">2</value>
        </metric>
    </metrics>
</test>