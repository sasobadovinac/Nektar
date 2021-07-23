<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Process surfDistance </description>
    <executable>FieldConvert</executable>
    <parameters> -f -e -m surfdistance:bnd=4 surfDist.xml surfDist.fld </parameters>
    <files>
        <file description="Session File">surfDist.xml</file>
    </files>
     <metrics>
        <metric type="L2" id="1">
            <value variable="x"   tolerance="1e-4">0.827414</value>
            <value variable="y"   tolerance="1e-4">0.0932301</value>
            <value variable="dist"   tolerance="1e-4">0.0483366</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="x"   tolerance="1e-4">0.999987</value>
            <value variable="y" tolerance="1e-4">0.0660123</value>
            <value variable="dist" tolerance="1e-4">0.0401285</value>
        </metric>
    </metrics>
</test>

