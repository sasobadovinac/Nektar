<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>StdInterpDerivBasis Pyr Mod basis P=7 Q=8</description>
    <executable>StdInterpDerivBasis</executable>
    <parameters> -s pyramid -b Modified_A Modified_A ModifiedPyr_C -o 7 7 7 -p 8 8 8 -P GaussGaussLegendre GaussGaussLegendre GaussGaussLegendre</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-12">0</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-12">0</value>
        </metric>
    </metrics>
</test>