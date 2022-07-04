<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Extract an isocontour</description>
    <executable>FieldConvert</executable>
    <parameters> -f -e -m CFL chan3D.xml chan3D.fld chan3D_cfl.fld</parameters>
    <files>
        <file description="Session File">chan3D.xml</file>
        <file description="Session File">chan3D.fld</file>
    </files>

     <metrics>
        <metric type="L2" id="1">
            <value variable="x" tolerance="1e-5">1.63299</value>
            <value variable="y" tolerance="1e-5">1.63299</value>
            <value variable="z" tolerance="1e-5">1.63299</value>
            <value variable="u" tolerance="1e-5">2.06559</value>
            <value variable="v" tolerance="1e-5">0</value>
            <value variable="w" tolerance="1e-5">0</value>
            <value variable="p" tolerance="1e-5">6.53197</value>
            <value variable="CFL" tolerance="1e-7">0.00619922</value>
        </metric>
    </metrics>
</test>

