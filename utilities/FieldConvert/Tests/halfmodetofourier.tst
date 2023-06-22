<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Process halfmodetofourier with w being transformed to imaginary mode </description>
    <executable>FieldConvert</executable>
    <parameters> -f -e -m halfmodetofourier:realmodetoimag=2  TriQuad.xml TriQuad.fld  TriQuad.plt</parameters>
    <files>
        <file description="Session File">TriQuad.xml</file>
	<file description="input fld File">TriQuad.fld</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="x"   tolerance="1e-4">0.816497</value>
            <value variable="y"   tolerance="1e-4">0.816497</value>
            <value variable="z"   tolerance="1e-4">0.661438</value>
            <value variable="u"   tolerance="1e-4">2.12687</value>
            <value variable="v"   tolerance="1e-4">0.455007</value>
            <value variable="w"   tolerance="1e-4">0.57735</value>
            <value variable="p"   tolerance="1e-4">1.35401</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="x"   tolerance="1e-4">1.</value>
            <value variable="y"   tolerance="1e-4">1.</value>
            <value variable="z"   tolerance="1e-4">0.75</value>
            <value variable="u"   tolerance="1e-4">6.28319</value>
            <value variable="v"   tolerance="1e-4">1.</value>
            <value variable="w"   tolerance="1e-4">1.</value>
            <value variable="p"   tolerance="1e-4">4.</value>
        </metric>
    </metrics>
</test>

