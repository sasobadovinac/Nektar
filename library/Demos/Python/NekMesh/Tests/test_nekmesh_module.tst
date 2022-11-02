<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Unit test of the Python interface for the Nektar::NekMesh::Module class.</description>
    <executable python="true">test_nekmesh_module.py</executable>
    <parameters>-v</parameters>
    <metrics>
        <metric type="pyunittest" id="1">
            <function>testCreateExceptionUnknownModule</function>
            <function>testCreateExceptionWrongArgs</function>
            <function>testCreateModuleUnknownParameter</function>
            <function>testExceptionNoMesh</function>
            <function>testInheritFromInputModuleTest</function>
            <function>testInheritFromOutputModuleTest</function>
            <function>testInheritFromProcessModuleTest</function>
            <function>testModuleProcessRuntimeError</function>
        </metric>
    </metrics>
</test>
