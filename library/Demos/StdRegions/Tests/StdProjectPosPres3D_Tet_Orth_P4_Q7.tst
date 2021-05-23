<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>StdProjectPositivityPres3D  Hex  Orthonormal basis P=4 Q=7</description>
    <executable>StdProjectPositivityPres3D</executable>
    <parameters>-s tetrahedron -b Ortho_A Ortho_B Ortho_C -o 4 4 4 -p 7 7 7 -z</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-12">5.37715e-16</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-12">1.16573e-15</value>
        </metric>
    </metrics>
</test>
