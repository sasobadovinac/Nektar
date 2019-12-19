<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>LSIAC 1D test on equal spaced mesh</description>
    <executable>test1DSIAC</executable>
    <parameters> P2_mesh_1D_80.xml 2</parameters>
    <files>
        <file description="Session File">P2_mesh_1D_80.xml</file>
    </files>

    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-9">3.98247e-07</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-9">2.55821e-06</value>
        </metric>
    </metrics>
</test>


