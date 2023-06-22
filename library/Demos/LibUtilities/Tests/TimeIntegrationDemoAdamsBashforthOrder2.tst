<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Test for time integration schemes</description>
    <executable>TimeIntegrationDemo</executable>
    <parameters>--dof 100 --timesteps 200 --method AdamsBashforth --order 2</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-12">0.00175531</value>
        </metric>
    </metrics>
</test>
