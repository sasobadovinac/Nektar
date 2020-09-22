<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>StdInterpSecondDeriv Seg Ortho_A basis P=6 Q=7</description>
    <executable>StdInterpSecondDeriv</executable>
    <parameters>-s segment -b Ortho_A -o 6 -p 7 -P GaussGaussLegendre</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-12">1e-12</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-12">1e-12</value>
        </metric>
    </metrics>
</test>
