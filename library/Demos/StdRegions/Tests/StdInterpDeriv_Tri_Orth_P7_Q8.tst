<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>StdInterpDeriv Tri Ortho basis P=7 Q=8</description>
    <executable>StdInterpDeriv</executable>
    <parameters>-s triangle -b Ortho_A Ortho_B -o 7 7 -p 8 8</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-12">5.480074e-14</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-12">1.793232e-12</value>
        </metric>
    </metrics>
</test>