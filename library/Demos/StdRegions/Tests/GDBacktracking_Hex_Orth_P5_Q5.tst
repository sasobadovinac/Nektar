<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>GDBackTrackingTest Hex  Ortho basis P=5 Q=5 </description>
    <executable>GDBackTrackingTest</executable>
    <parameters>-s hexahedron -b Ortho_A Ortho_A Ortho_A -o 5 5 5   -p 5 5 5   -g 0.7 -c 0.4 </parameters>
    <metrics>
	<metric type="L2" id="1">
            <value tolerance="1e-1"> 1e-4 </value>
	</metric>
	<metric type="Linf" id="2">
            <value tolerance="1e-1"> 1e-3 </value>
	</metric>
    </metrics>
</test>

