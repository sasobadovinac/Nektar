<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Unit test of the Python interface for the Nektar::NekMesh::Mesh class.</description>
    <executable python="true">test_nekmesh_mesh.py</executable>
    <parameters>-v</parameters>
    <metrics>
        <metric type="pyunittest" id="1">
            <function>testMeshConstructor</function>
            <function>testMeshFieldAccess</function>
        </metric>
    </metrics>
</test>
