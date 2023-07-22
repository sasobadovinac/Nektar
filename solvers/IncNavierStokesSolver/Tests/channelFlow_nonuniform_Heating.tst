<?xml version="1.0" encoding="utf-8" ?>
<test>
   <description> Constant flow-rate Channel Flow with temperature field P=3 </description>
   <executable>IncNavierStokesSolver</executable>
   <parameters> channelFlow_nonuniform_Heating.xml</parameters>
   <files>
        <file description ="Session File"> channelFlow_nonuniform_Heating.xml </file>
   </files>
   <metrics>
     <metric type="L2" id="1">
       <value variable ="x" tolerance="1e-5">2.3094</value>
       <value variable ="y" tolerance="1e-5">1.1547 </value>
       <value variable ="u" tolerance="1e-5">1.46059</value>
       <value variable ="v" tolerance="1e-10">0 </value>
       <value variable ="teta" tolerance="1e-5">1.1547 </value>
       <value variable ="p" tolerance="1e-4">205.079 </value>
       <value variable ="uu" tolerance="1e-10">0 </value>
       <value variable ="uv" tolerance="1e-10">0 </value>
       <value variable ="uteta" tolerance="1e-10">0</value>
       <value variable ="vv" tolerance="1e-10">0 </value>
       <value variable ="vteta" tolerance="1e-10">0 </value>
       <value variable ="tetateta" tolerance="1e-10">0 </value>
       <value variable ="u" tolerance="1e-10">0</value>
       <value variable ="v" tolerance="1e-10">0 </value>
       <value variable ="teta" tolerance="1e-5">0 </value>
       <value variable ="p" tolerance="1e-4">0 </value>
     </metric>
      <metric type="Linf" id="1">
       <value variable ="x" tolerance="1e-5">2</value>
       <value variable ="y" tolerance="1e-5">1 </value>
       <value variable ="u" tolerance="1e-5">1</value>
       <value variable ="v" tolerance="1e-10">0 </value>
       <value variable ="teta" tolerance="1e-5">1 </value>
       <value variable ="p" tolerance="1e-4">140.845 </value>
       <value variable ="uu" tolerance="1e-10">0 </value>
       <value variable ="uv" tolerance="1e-10">0 </value>
       <value variable ="uteta" tolerance="1e-10">0</value>
       <value variable ="vv" tolerance="1e-10">0 </value>
       <value variable ="vteta" tolerance="1e-10">0 </value>
       <value variable ="tetateta" tolerance="1e-10">0 </value>
       <value variable ="u" tolerance="1e-10">0 </value>
       <value variable ="v" tolerance="1e-10">0 </value>
       <value variable ="teta" tolerance="1e-10">0 </value>
       <value variable ="p" tolerance="1e-10">0</value>
     </metric>
   </metrics>
</test>