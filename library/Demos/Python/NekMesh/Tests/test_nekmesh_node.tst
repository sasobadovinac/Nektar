<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Unit test of the Python interface for the Nektar::NekMesh::Node class.</description>
    <executable python="true">test_nekmesh_node.py</executable>
    <parameters>-v</parameters>
    <metrics>
        <metric type="pyunittest" id="1">
            <function>testNodeAbs2</function>
            <function>testNodeConstructor</function>
            <function>testNodeDistance</function>
            <function>testNodeFieldAccess</function>
            <function>testNodeGetID</function>
            <function>testNodeGetLoc</function>
            <function>testNodeSetID</function>
            <function>testNodeSetAdd</function>
            <function>testNodeSetClear</function>
            <function>testNodeSet__contains__</function>
            <function>testNodeSet__iter__</function>
            <function>testNodeSet__len__</function>
        </metric>
    </metrics>
</test>
