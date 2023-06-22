<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Test for time integration schemes</description>
    <executable>TimeIntegrationDemo</executable>
    <parameters>--dof 10000 --timesteps 10 --method ExtrapolationMethod --variant ImplicitMidpoint --order 2</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-12">0.00160331</value>
        </metric>
    </metrics>
</test>
