<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>GDBackTrackingTest Tet Ortho basis P=5 Q=5 f= -2x + (x+0.6)^3 + y^2 - 0.2</description>
    <executable>GDBackTrackingTest</executable>
    <parameters>-s tetrahedron -b Ortho_A Ortho_B Ortho_C -o 6 6 6  -p 6 6 6  -P GaussGaussLegendre GaussGaussLegendre GaussGaussLegendre -c 0.4 -g 0.7</parameters>
    <metrics>
	<metric type="L2" id="1">
            <value tolerance="1e-2">1e-5</value>
	</metric>
	<metric type="Linf" id="2">
            <value tolerance="1e-2">1e-5</value>
	</metric>
    </metrics>
</test>



