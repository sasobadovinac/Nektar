<?xml version="1.0" encoding="utf-8" ?>
<test>
   <description> Constant flow-rate Channel Flow with temperature field P=3 </description>
   <executable>IncNavierStokesSolver</executable>
   <parameters> ChannelFlow_Heating_FlowRate.xml</parameters>
   <files>
        <file description ="Session File"> ChannelFlow_Heating_FlowRate.xml </file>
   </files>
   <metrics>
     <metric type="L2" id="1">
       <value variable ="u" tolerance="1e-10">3.46569e-15</value>
       <value variable ="v" tolerance="1e-10">4.42773e-15 </value>
       <value variable ="teta" tolerance="1e-5">5.07433e-15 </value>
       <value variable ="p" tolerance="1e-4">4.40726e-13 </value>
     </metric>
      <metric type="Linf" id="1">
       <value variable ="u" tolerance="1e-10">6.77236e-15 </value>
       <value variable ="v" tolerance="1e-10">5.89886e-15 </value>
       <value variable ="teta" tolerance="1e-5">4.21885e-15 </value>
       <value variable ="p" tolerance="1e-4">4.64944e-13</value>
     </metric>
   </metrics>
</test>