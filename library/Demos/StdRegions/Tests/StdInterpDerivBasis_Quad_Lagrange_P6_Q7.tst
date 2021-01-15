<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>StdInterpDerivBasis Quadrilateral Lagrange basis P=62D Q=7</description>
    <executable>StdInterpDerivBasis</executable>
    <parameters>-s quadrilateral -b GLL_Lagrange GLL_Lagrange -o 6 6 -p 7 7 </parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-12"> -0.000000e+00</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-12"> 0.000000e+00</value>
        </metric>
    </metrics>
</test>