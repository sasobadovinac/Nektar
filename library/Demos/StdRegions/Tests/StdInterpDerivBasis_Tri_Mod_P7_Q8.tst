<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>StdInterpDerivBasis Tri Lagrange basis P=6 Q=7</description>
    <executable>StdInterpDerivBasis</executable>
    <parameters>-s triangle -b Modified_A Modified_B -o 7 7 -p 8 8</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-12">0</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-12">0</value>
        </metric>
    </metrics>
</test>
