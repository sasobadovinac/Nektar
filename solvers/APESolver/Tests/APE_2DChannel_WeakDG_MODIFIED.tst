<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>2D channel with base flow, P=5</description>
    <executable>APESolver</executable>
    <parameters>APE_2DChannel_WeakDG_MODIFIED.xml</parameters>
    <files>
        <file description="Session File">APE_2DChannel_WeakDG_MODIFIED.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="p" tolerance="1e-12">12.8017</value>
            <value variable="u" tolerance="1e-12">0.0314268</value>
            <value variable="v" tolerance="1e-12">0.00291624</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="p" tolerance="1e-12">57.9942</value>
            <value variable="u" tolerance="1e-12">0.180743</value>
            <value variable="v" tolerance="1e-12">0.174942</value>
        </metric>
    </metrics>
</test>
