<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>StdInterpDeriv Pyr Mod basis P=7 Q=8</description>
    <executable>StdInterpDeriv</executable>
    <parameters> -s pyramid -b Modified_A Modified_A ModifiedPyr_C -o 7 7 7 -p
8 8 8 -P GaussGaussChebyshev GaussGaussChebyshev GaussGaussChebyshev</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-11">1e-15</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-11">1e-14</value>
        </metric>
    </metrics>
</test>