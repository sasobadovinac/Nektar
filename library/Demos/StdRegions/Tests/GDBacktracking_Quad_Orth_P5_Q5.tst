<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>GDBackTrackingTest  Quad Ortho basis P=5 Q=5 </description>
    <executable>GDBackTrackingTest</executable>
    <parameters>-s quadrilateral -b Ortho_A Ortho_A -o 5 5 -p 5 5 -g 0.7 -c 0.4 </parameters>
    <metrics>
	<metric type="L2" id="1">
            <value tolerance="1e-2"> 1e-5 </value>
	</metric>
	<metric type="Linf" id="2">
            <value tolerance="1e-2"> 1e-5 </value>
	</metric>
    </metrics>
</test>

