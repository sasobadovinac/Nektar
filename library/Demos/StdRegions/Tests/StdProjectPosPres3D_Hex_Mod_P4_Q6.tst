<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>StdProjectPositivityPres3D  Hex  Modified basis P=4 Q=6</description>
    <executable>StdProjectPositivityPres3D</executable>
    <parameters>-s hexahedron -b Modified_A Modified_A Modified_A -o 4 4 4 -p 6 6 6 -z</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-12">5.37715e-16</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-12">1.16573e-15</value>
        </metric>
    </metrics>
</test>