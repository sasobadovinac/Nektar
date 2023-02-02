<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>StdInterpDeriv Seg Lagrange basis P=7 Q=8</description>
    <executable>StdInterpDeriv</executable>
    <parameters> -s segment -b GLL_Lagrange -o 7 -p 8 -P GaussLobattoLegendre</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-7">0</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-7">0</value>
        </metric>
    </metrics>
</test>