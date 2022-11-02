<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Unit test of the Python interface for the Nektar::NekMesh::Element class.</description>
    <executable python="true">test_nekmesh_element.py</executable>
    <parameters>-v</parameters>
    <metrics>
        <metric type="pyunittest" id="1">
            <function>testElementGetDim</function>
            <function>testElementGetId</function>
            <function>testElementGetShapeType</function>
            <function>testElementGetTag</function>
            <function>testElmtConfigConstructor</function>
            <function>testNoDefaultConstructor</function>
        </metric>
    </metrics>
</test>
