<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>StdInterpDerivBasis Tri Ortho basis P=7 Q=8</description>
    <executable>StdInterpDerivBasis</executable>
    <parameters> -s triangle -b Ortho_A Ortho_B -o 7 7 -p 8 8</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-12">0</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-12">0</value>
        </metric>
    </metrics>
</test>
