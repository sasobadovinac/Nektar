<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Moving Frame of Reference test </description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>MovingRefFrame_translation.xml</parameters>
    <files>
        <file description="Session File">MovingRefFrame_translation.xml</file>
	<file description="Restart File">MovingRefFrame_translation.rst</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="5e-5">5.04747</value>
            <value variable="v" tolerance="5e-5">1.2054</value>
	        <value variable="p" tolerance="5e-5">1.78602</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="5e-5">1.99607</value>
            <value variable="v" tolerance="5e-5">0.799828</value>
	        <value variable="p" tolerance="5e-5">0.960228</value>
        </metric>
    </metrics>
</test>

