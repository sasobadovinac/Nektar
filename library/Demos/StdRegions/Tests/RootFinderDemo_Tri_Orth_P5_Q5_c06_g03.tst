<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>RootFinderTest Tri Ortho basis P=5 Q=5 f= -2x + (x+0.6)^3 + y^2 - 0.2</description>
    <executable>RootFinderTest</executable>
    <parameters>-s triangle -b Ortho_A Ortho_B -o 5 5 -p 5 5 -g 0.3 -c 0.6 -z</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-2">1e-5</value>
        </metric>
        <metric type="Linf" id="2">
 	    <value tolerance="1e-2">1e-5</value>
        </metric>  
    </metrics>
</test>