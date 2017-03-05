<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Option Pricing Solver Test</description>
    <executable>OptionPricingSolver</executable>
    <parameters>BlackScholes1DTest.xml</parameters>
    <files>
        <file description="Session File">BlackScholes1DTest.xml</file>
    </files>
    <metrics>
       <metric type="L2" id="1">
            <value variable="u" tolerance="1e-8">5.1197e-08</value>
            </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-7">1.59721e-07</value>
            </metric>
    </metrics>
</test>
