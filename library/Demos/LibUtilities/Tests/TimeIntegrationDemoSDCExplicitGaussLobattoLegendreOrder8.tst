<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Test for time integration schemes</description>
    <executable>TimeIntegrationDemo</executable>
    <parameters>--dof 800 --timesteps 400 --method ExplicitSDC --variant GaussLobattoLegendre --order 8 --param 1.0 5</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-10">6.48390e-05</value>
        </metric>
    </metrics>
</test>
