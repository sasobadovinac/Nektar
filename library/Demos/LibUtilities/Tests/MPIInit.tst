<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Test MPI Initialisation</description>
    <executable>MPIInit</executable>
    <parameters />
    <processes>3</processes>
    <metrics>
        <metric type="regex" id="1">
            <regex>^Obtained reduce = (\d+)</regex>
            <matches>
                <match>
                    <field id="0">3</field>
                </match>
            </matches>
        </metric>
    </metrics>
</test>
