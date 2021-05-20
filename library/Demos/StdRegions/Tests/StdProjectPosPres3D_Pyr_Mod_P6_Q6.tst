<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>StdProjectPositivityPres3D  Pyr  Modified basis P=6 Q=6</description>
    <executable>StdProjectPositivityPres3D</executable>
    <parameters>-s pyramid -b Modified_A Modified_A ModifiedPyr_C -o 6 6 6 -p 6 6 6  -z</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-12">5.37715e-16</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-12">1.16573e-15</value>
        </metric>
    </metrics>
</test>