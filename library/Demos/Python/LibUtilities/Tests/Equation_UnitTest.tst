<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Unit test of the Python interface for the Nektar::LibUtilities::Equation class.</description>
    <executable python="true">Equation_UnitTest.py</executable>
    <parameters>-v</parameters>
    <metrics>
        <metric type="pyunittest" id="1">
            <function>testEquationEvaluate</function>
            <function>testEquationGetExpression</function>
            <function>testEquationGetVlist</function>
            <function>testEquationSetConstants</function>
            <function>testEquationSetParameter</function>
        </metric>
    </metrics>
</test>
