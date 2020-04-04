<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>BaryInterpBasis Quadrilateral Lagrange basis P=6 Q=7</description>
    <executable>BaryInterpBasis</executable>
    <parameters>-s quadrilateral -b GLL_Lagrange GLL_Lagrange -o 6 6 -p 7 7</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-12">0.00024527</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-12">0.000360131</value>
        </metric>
    </metrics>
</test>


