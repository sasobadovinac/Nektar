<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>LSIAC 2D test on equal spaced mesh</description>
    <executable>test2DSIAC</executable>
    <parameters> P2_quadSDExtra21.xml 2</parameters>
    <files>
        <file description="Session File">P2_quadSDExtra21.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-5">8.17919e-03</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-5">8.64424e-03</value>
        </metric>
    </metrics>
</test>
