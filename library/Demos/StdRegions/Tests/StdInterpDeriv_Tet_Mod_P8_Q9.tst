<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>StdInterpDeriv Tet Mod basis P=8 Q=9</description>
    <executable>StdInterpDeriv</executable>
    <parameters> -s tetrahedron -b Modified_A Modified_B Modified_C -o 8 8 8 -p 9 9 9 -P GaussLobattoLegendre GaussRadauMAlpha1Beta0 GaussRadauMAlpha2Beta0</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-7">0</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-7">0</value>
        </metric>
    </metrics>
</test>