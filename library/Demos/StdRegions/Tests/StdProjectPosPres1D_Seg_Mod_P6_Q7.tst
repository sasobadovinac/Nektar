<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>StdProjectPositivityPres1D Segment Modified basis P=6 Q=7</description>
    <executable>StdProjectPositivityPres1D</executable>
    <parameters>-s Segment -b Modified_A -o 6 -p 7 -z</parameters>
    <metrics>
	<metric type="L2" id="1">
            <value tolerance="1e-12">5.37715e-16</value>
	</metric>
	<metric type="Linf" id="2">
            <value tolerance="1e-12">1.16573e-15</value>
	</metric>
    </metrics>
</test>
