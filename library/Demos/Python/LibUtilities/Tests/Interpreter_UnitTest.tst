<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Unit test of the Python interface for the Nektar::LibUtilities::Interpreter class.</description>
    <executable python="true">Interpreter_UnitTest.py</executable>
    <parameters>-v</parameters>
    <metrics>
        <metric type="pyunittest" id="1">
            <function>testInterpreterEvaluate</function>
            <function>testInterpreterGetConstant</function>
            <function>testInterpreterGetConstants</function>
            <function>testInterpreterSetParameter</function>
            <function>testInterpreterSetParameters</function>
        </metric>
    </metrics>
</test>
