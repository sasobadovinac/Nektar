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
            <value variable="u" tolerance="1e-6">2.91653e-06</value>
            </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-6">3.58428e-06</value>
            </metric>
    </metrics>
</test>
