<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>StdInterpDeriv Tri Modified basis P=7 Q=8</description>
    <executable>StdInterpDeriv</executable>
    <parameters>-s triangle -b Modified_A Modified_B -o 7 7 -p 8 8 -P GaussLobattoLegendre GaussRadauMAlpha1Beta0</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-7">0</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-7">0</value>
        </metric>
    </metrics>
</test>
