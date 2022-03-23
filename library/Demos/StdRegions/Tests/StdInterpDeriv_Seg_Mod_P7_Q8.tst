<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>StdInterpDeriv Seg Mod basis P=7 Q=8</description>
    <executable>StdInterpDeriv</executable>
    <parameters> -s segment -b Modified_A -o 7 -p 8 -P GaussGaussLegendre</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-12">0</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-12">0</value>
        </metric>
    </metrics>
</test>