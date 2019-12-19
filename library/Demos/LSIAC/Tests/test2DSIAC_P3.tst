<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>LSIAC 2D test on equal spaced mesh</description>
    <executable>test2DSIAC</executable>
    <parameters> P3_quadSDExtra21.xml 3</parameters>
    <files>
        <file description="Session File">P3_quadSDExtra21.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-6">1.9969e-04</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-6">8.96892e-04</value>
        </metric>
    </metrics>
</test>


