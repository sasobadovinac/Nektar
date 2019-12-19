<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>LSIAC 1D test on equal spaced mesh</description>
    <executable>test1DSIAC</executable>
    <parameters> P3_mesh_1D_80.xml 3</parameters>
    <files>
        <file description="Session File">P3_mesh_1D_80.xml</file>
    </files>

    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-11">1.61077e-09</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-11">2.19729e-08</value>
        </metric>
    </metrics>
</test>


