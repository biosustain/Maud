<?xml version="1.0" encoding="UTF-8"?>
<!-- generated with COPASI 4.24 (Build 197) (http://www.copasi.org) at 2019-03-08 15:07:34 UTC -->
<?oxygen RNGSchema="http://www.copasi.org/static/schema/CopasiML.rng" type="xml"?>
<COPASI xmlns="http://www.copasi.org/static/schema" versionMajor="4" versionMinor="24" versionDevel="197" copasiSourcesModified="0">
  <ListOfFunctions>
    <Function key="Function_14" name="Mass action (reversible)" type="MassAction" reversible="true">
      <MiriamAnnotation>
<rdf:RDF xmlns:CopasiMT="http://www.copasi.org/RDF/MiriamTerms#" xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
   <rdf:Description rdf:about="#Function_14">
   <CopasiMT:is rdf:resource="urn:miriam:obo.sbo:SBO:0000042" />
   </rdf:Description>
   </rdf:RDF>
      </MiriamAnnotation>
      <Comment>
        <body xmlns="http://www.w3.org/1999/xhtml">
<b>Mass action rate law for reversible reactions</b>
<p>
Reaction scheme where the products are created from the reactants and the change of a product quantity is proportional to the product of reactant activities. The reaction scheme does include a reverse process that creates the reactants from the products.
</p>
</body>
      </Comment>
      <Expression>
        k1*PRODUCT&lt;substrate_i>-k2*PRODUCT&lt;product_j>
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_69" name="k1" order="0" role="constant"/>
        <ParameterDescription key="FunctionParameter_68" name="substrate" order="1" role="substrate"/>
        <ParameterDescription key="FunctionParameter_78" name="k2" order="2" role="constant"/>
        <ParameterDescription key="FunctionParameter_79" name="product" order="3" role="product"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_41" name="Function for PGK" type="UserDefined" reversible="true">
      <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Function_41">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-26T08:49:30Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

      </MiriamAnnotation>
      <Expression>
        Vmax*(MgADP*BPG-MgATP*PGA3/Keq)/(KmADPMg*KmBPG)/(1+MgADP/KmADPMg+BPG/KmBPG+MgADP/KmADPMg*BPG/KmBPG+MgATP/KmATPMg+PGA3/KmPGA3+MgATP/KmATPMg*PGA3/KmPGA3)
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_317" name="BPG" order="0" role="substrate"/>
        <ParameterDescription key="FunctionParameter_316" name="Keq" order="1" role="constant"/>
        <ParameterDescription key="FunctionParameter_315" name="KmADPMg" order="2" role="constant"/>
        <ParameterDescription key="FunctionParameter_314" name="KmATPMg" order="3" role="constant"/>
        <ParameterDescription key="FunctionParameter_313" name="KmBPG" order="4" role="constant"/>
        <ParameterDescription key="FunctionParameter_312" name="KmPGA3" order="5" role="constant"/>
        <ParameterDescription key="FunctionParameter_311" name="MgADP" order="6" role="modifier"/>
        <ParameterDescription key="FunctionParameter_310" name="MgATP" order="7" role="modifier"/>
        <ParameterDescription key="FunctionParameter_309" name="PGA3" order="8" role="product"/>
        <ParameterDescription key="FunctionParameter_308" name="Vmax" order="9" role="constant"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_42" name="Function for GDH" type="UserDefined" reversible="true">
      <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Function_42">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-25T19:43:08Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

      </MiriamAnnotation>
      <Expression>
        Vmax*(P*GAP*NAD-BPG*NADH/Keq)/(KmP*KmGAP*KmNAD)/((1+P/KmP)*(1+GAP/KmGAP)*(1+NAD/KmNAD)+(1+BPG/KmBPG)*(1+NADH/KmNADH)-1)
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_298" name="BPG" order="0" role="product"/>
        <ParameterDescription key="FunctionParameter_299" name="GAP" order="1" role="substrate"/>
        <ParameterDescription key="FunctionParameter_300" name="Keq" order="2" role="constant"/>
        <ParameterDescription key="FunctionParameter_301" name="KmBPG" order="3" role="constant"/>
        <ParameterDescription key="FunctionParameter_302" name="KmGAP" order="4" role="constant"/>
        <ParameterDescription key="FunctionParameter_303" name="KmNAD" order="5" role="constant"/>
        <ParameterDescription key="FunctionParameter_304" name="KmNADH" order="6" role="constant"/>
        <ParameterDescription key="FunctionParameter_305" name="KmP" order="7" role="constant"/>
        <ParameterDescription key="FunctionParameter_306" name="NAD" order="8" role="substrate"/>
        <ParameterDescription key="FunctionParameter_307" name="NADH" order="9" role="product"/>
        <ParameterDescription key="FunctionParameter_297" name="P" order="10" role="substrate"/>
        <ParameterDescription key="FunctionParameter_296" name="Vmax" order="11" role="constant"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_44" name="Function for ENO" type="UserDefined" reversible="true">
      <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Function_44">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2019-03-08T15:58:14Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

      </MiriamAnnotation>
      <Expression>
        Vmax*(PGA2-PEP/Keq)/KmPGA2/(1+PGA2/KmPGA2+PEP/KmPEP)
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_295" name="Keq" order="0" role="constant"/>
        <ParameterDescription key="FunctionParameter_294" name="KmPEP" order="1" role="constant"/>
        <ParameterDescription key="FunctionParameter_293" name="KmPGA2" order="2" role="constant"/>
        <ParameterDescription key="FunctionParameter_292" name="PEP" order="3" role="product"/>
        <ParameterDescription key="FunctionParameter_291" name="PGA2" order="4" role="substrate"/>
        <ParameterDescription key="FunctionParameter_290" name="Vmax" order="5" role="constant"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_46" name="Function for FBP" type="UserDefined" reversible="false">
      <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Function_46">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-25T19:42:50Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

      </MiriamAnnotation>
      <Expression>
        Vmax*n*MgFDP/KirFDPMg/(1+KmrFDP/KirFDP*(MG/KmrMg)+P/KirP+P/KirP*(MG/KirPMg)+F6P/KirF6P+F6P/KirF6P*(MG/KirF6PMg)+P/KirP*(F6P/KirPF6P)+P/KirP*(F6P/KirPF6P)*(MG/KirPF6PMg)+(FDP-MgFDP)/KirFDP+KdFDPMg/KmrMg*(MgFDP/KirFDP)+AMP/KirAMP+MgFDP/KirFDPMg+MgFDP/KirFDPMg*(MG/KirFDPMgMg)+AMP/KirAMP*((FDP-MgFDP)/KirAMPFDP))/(1+L0*((1+KmtFDP/KitFDP*(MG/KmtMg)+P/KitP+P/KitP*(MG/KitPMg)+F6P/KitF6P+F6P/KitF6P*(MG/KitF6PMg)+P/KitP*(F6P/KitPF6P)+P/KitP*(F6P/KitPF6P)*(MG/KitPF6PMg)+(FDP-MgFDP)/KitFDP+KdFDPMg/KmtMg*(MgFDP/KitFDP)+AMP/KitAMP+MgFDP/KitFDPMg+MgFDP/KitFDPMg*(MG/KitFDPMgMg)+AMP/KitAMP*((FDP-MgFDP)/KitAMPFDP))/(1+KmrFDP/KirFDP*(MG/KmrMg)+P/KirP+P/KirP*(MG/KirPMg)+F6P/KirF6P+F6P/KirF6P*(MG/KirF6PMg)+P/KirP*(F6P/KirPF6P)+P/KirP*(F6P/KirPF6P)*(MG/KirPF6PMg)+(FDP-MgFDP)/KirFDP+KdFDPMg/KmrMg*(MgFDP/KirFDP)+AMP/KirAMP+MgFDP/KirFDPMg+MgFDP/KirFDPMg*(MG/KirFDPMgMg)+AMP/KirAMP*((FDP-MgFDP)/KirAMPFDP)))^n)
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_334" name="AMP" order="0" role="modifier"/>
        <ParameterDescription key="FunctionParameter_333" name="F6P" order="1" role="product"/>
        <ParameterDescription key="FunctionParameter_332" name="FDP" order="2" role="substrate"/>
        <ParameterDescription key="FunctionParameter_331" name="KdFDPMg" order="3" role="constant"/>
        <ParameterDescription key="FunctionParameter_330" name="KirAMP" order="4" role="constant"/>
        <ParameterDescription key="FunctionParameter_329" name="KirAMPFDP" order="5" role="constant"/>
        <ParameterDescription key="FunctionParameter_335" name="KirF6P" order="6" role="constant"/>
        <ParameterDescription key="FunctionParameter_336" name="KirF6PMg" order="7" role="constant"/>
        <ParameterDescription key="FunctionParameter_337" name="KirFDP" order="8" role="constant"/>
        <ParameterDescription key="FunctionParameter_338" name="KirFDPMg" order="9" role="constant"/>
        <ParameterDescription key="FunctionParameter_339" name="KirFDPMgMg" order="10" role="constant"/>
        <ParameterDescription key="FunctionParameter_340" name="KirP" order="11" role="constant"/>
        <ParameterDescription key="FunctionParameter_341" name="KirPF6P" order="12" role="constant"/>
        <ParameterDescription key="FunctionParameter_342" name="KirPF6PMg" order="13" role="constant"/>
        <ParameterDescription key="FunctionParameter_343" name="KirPMg" order="14" role="constant"/>
        <ParameterDescription key="FunctionParameter_344" name="KitAMP" order="15" role="constant"/>
        <ParameterDescription key="FunctionParameter_345" name="KitAMPFDP" order="16" role="constant"/>
        <ParameterDescription key="FunctionParameter_346" name="KitF6P" order="17" role="constant"/>
        <ParameterDescription key="FunctionParameter_347" name="KitF6PMg" order="18" role="constant"/>
        <ParameterDescription key="FunctionParameter_348" name="KitFDP" order="19" role="constant"/>
        <ParameterDescription key="FunctionParameter_349" name="KitFDPMg" order="20" role="constant"/>
        <ParameterDescription key="FunctionParameter_350" name="KitFDPMgMg" order="21" role="constant"/>
        <ParameterDescription key="FunctionParameter_351" name="KitP" order="22" role="constant"/>
        <ParameterDescription key="FunctionParameter_352" name="KitPF6P" order="23" role="constant"/>
        <ParameterDescription key="FunctionParameter_353" name="KitPF6PMg" order="24" role="constant"/>
        <ParameterDescription key="FunctionParameter_354" name="KitPMg" order="25" role="constant"/>
        <ParameterDescription key="FunctionParameter_355" name="KmrFDP" order="26" role="constant"/>
        <ParameterDescription key="FunctionParameter_356" name="KmrMg" order="27" role="constant"/>
        <ParameterDescription key="FunctionParameter_357" name="KmtFDP" order="28" role="constant"/>
        <ParameterDescription key="FunctionParameter_358" name="KmtMg" order="29" role="constant"/>
        <ParameterDescription key="FunctionParameter_359" name="L0" order="30" role="constant"/>
        <ParameterDescription key="FunctionParameter_360" name="MG" order="31" role="modifier"/>
        <ParameterDescription key="FunctionParameter_361" name="MgFDP" order="32" role="modifier"/>
        <ParameterDescription key="FunctionParameter_362" name="P" order="33" role="product"/>
        <ParameterDescription key="FunctionParameter_363" name="Vmax" order="34" role="constant"/>
        <ParameterDescription key="FunctionParameter_364" name="n" order="35" role="constant"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_50" name="Function for PPS" type="UserDefined" reversible="true">
      <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Function_50">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-25T19:43:41Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

      </MiriamAnnotation>
      <Expression>
        Vmax*(MgATP*PYR-AMP*PEP*P*MG/Keq)/(KmATPMg*KmPYR)/(MgATP/KmATPMg+alpha*(P/KdP)*(MgATP/KmATPMg)+alpha*(AMP/KdAMP)*(MgATP/KmATPMg)+alpha*(P/KdP)*(AMP/KdAMP)*(MgATP/KmATPMg)+alpha*(MG/KdMg)*(P/KmP)*(AMP/KdAMP)*(MgATP/KdATPMgPPS)/(W*(1+MG/KdMg))+MgATP/KmATPMg*(AKG/KefAKG)+(1+MG/KdMg)*(AKG/KefAKG)*(PEP/KmPEP)/W+MgATP/KmATPMg*(OAA/KefOAA)+(1+MG/KdMg)*(OAA/KefOAA)*(PEP/KmPEP)/W+MG/KdMg*(P/KmP)*(AMP/KdAMP)/W+alpha*(P/KdP)*(AMP/KdAMP)*(PEP/KmPEP)/W+alpha*(MG/KdMg)*(P/KmP)*(AMP/KdAMP)*(PEP/KmPEP)/W+alpha*(1+MG/KdMg)*(KmAMP/KdAMP*(P/KmP)*(PEP/KmPEP)+AMP/KdAMP*(PEP/KmPEP))/W+(1+MG/KdMg)*(PYR/KmPYR)+MgATP/KmATPMg*(PYR/KmPYR)+KdADPMg/KdMg*(P/KmP)*(MgADP/KefADP)*(AMP/KdAMP)/(W*(1+MG/KdMg))+(ADP-MgADP)/KefADP*(PYR/KmPYR)+KdATPMg/KdMg*(P/KmP)*(AMP/KdAMP)*(MgATP/KefATP)/(W*(1+MG/KdMg))+(ATP-MgATP)/KefATP*(PYR/KmPYR)+(1+MG/KdMg)*(PEP/KmPEP)/W+alpha*(1+MG/KdMg)*(PEP/KdPEP)*(PYR/KmPYR)+(1+MG/KdMg)*(PYR/KdPYR)*(PEP/KmPEP)/W)
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_420" name="ADP" order="0" role="modifier"/>
        <ParameterDescription key="FunctionParameter_421" name="AKG" order="1" role="modifier"/>
        <ParameterDescription key="FunctionParameter_422" name="AMP" order="2" role="product"/>
        <ParameterDescription key="FunctionParameter_423" name="ATP" order="3" role="substrate"/>
        <ParameterDescription key="FunctionParameter_424" name="KdADPMg" order="4" role="constant"/>
        <ParameterDescription key="FunctionParameter_425" name="KdAMP" order="5" role="constant"/>
        <ParameterDescription key="FunctionParameter_426" name="KdATPMg" order="6" role="constant"/>
        <ParameterDescription key="FunctionParameter_427" name="KdATPMgPPS" order="7" role="constant"/>
        <ParameterDescription key="FunctionParameter_428" name="KdMg" order="8" role="constant"/>
        <ParameterDescription key="FunctionParameter_429" name="KdP" order="9" role="constant"/>
        <ParameterDescription key="FunctionParameter_430" name="KdPEP" order="10" role="constant"/>
        <ParameterDescription key="FunctionParameter_443" name="KdPYR" order="11" role="constant"/>
        <ParameterDescription key="FunctionParameter_419" name="KefADP" order="12" role="constant"/>
        <ParameterDescription key="FunctionParameter_418" name="KefAKG" order="13" role="constant"/>
        <ParameterDescription key="FunctionParameter_417" name="KefATP" order="14" role="constant"/>
        <ParameterDescription key="FunctionParameter_416" name="KefOAA" order="15" role="constant"/>
        <ParameterDescription key="FunctionParameter_415" name="Keq" order="16" role="constant"/>
        <ParameterDescription key="FunctionParameter_414" name="KmAMP" order="17" role="constant"/>
        <ParameterDescription key="FunctionParameter_413" name="KmATPMg" order="18" role="constant"/>
        <ParameterDescription key="FunctionParameter_412" name="KmP" order="19" role="constant"/>
        <ParameterDescription key="FunctionParameter_411" name="KmPEP" order="20" role="constant"/>
        <ParameterDescription key="FunctionParameter_457" name="KmPYR" order="21" role="constant"/>
        <ParameterDescription key="FunctionParameter_458" name="MG" order="22" role="modifier"/>
        <ParameterDescription key="FunctionParameter_459" name="MgADP" order="23" role="modifier"/>
        <ParameterDescription key="FunctionParameter_460" name="MgATP" order="24" role="modifier"/>
        <ParameterDescription key="FunctionParameter_461" name="OAA" order="25" role="modifier"/>
        <ParameterDescription key="FunctionParameter_462" name="P" order="26" role="product"/>
        <ParameterDescription key="FunctionParameter_463" name="PEP" order="27" role="product"/>
        <ParameterDescription key="FunctionParameter_464" name="PYR" order="28" role="substrate"/>
        <ParameterDescription key="FunctionParameter_465" name="Vmax" order="29" role="constant"/>
        <ParameterDescription key="FunctionParameter_466" name="W" order="30" role="constant"/>
        <ParameterDescription key="FunctionParameter_467" name="alpha" order="31" role="constant"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_52" name="Function for FBA" type="UserDefined" reversible="true">
      <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Function_52">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2014-12-29T22:00:27Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

      </MiriamAnnotation>
      <Expression>
        Vmax*(FDP-DAP*GAP/Keq)/KmFDP/(1+FDP/KmFDP+DAP/KmDAP+DAP/KmDAP*(GAP/KmGAP)+PEP/KmPEP)
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_513" name="DAP" order="0" role="product"/>
        <ParameterDescription key="FunctionParameter_512" name="FDP" order="1" role="substrate"/>
        <ParameterDescription key="FunctionParameter_511" name="GAP" order="2" role="product"/>
        <ParameterDescription key="FunctionParameter_510" name="Keq" order="3" role="constant"/>
        <ParameterDescription key="FunctionParameter_509" name="KmDAP" order="4" role="constant"/>
        <ParameterDescription key="FunctionParameter_508" name="KmFDP" order="5" role="constant"/>
        <ParameterDescription key="FunctionParameter_507" name="KmGAP" order="6" role="constant"/>
        <ParameterDescription key="FunctionParameter_506" name="KmPEP" order="7" role="constant"/>
        <ParameterDescription key="FunctionParameter_505" name="PEP" order="8" role="modifier"/>
        <ParameterDescription key="FunctionParameter_504" name="Vmax" order="9" role="constant"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_54" name="Function for PFK" type="UserDefined" reversible="true">
      <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Function_54">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2014-12-29T22:03:23Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

      </MiriamAnnotation>
      <Expression>
        Vmax*n*(MgATP*F6P-MgADP*FDP/Keq)/(KirF6P*KmrATPMg)/(1+KmrFDP/KirFDP*(MgADP/KmrADP)+KmrF6P/KirF6P*(MgATP/KmrATPMg)+KmrFDP/KirFDP*(MgADP/KmrADP)*(F6P/KirF6P)+MgATP/KmrATPMg*(F6P/KirF6P)+MgADP/KirADP*(MgATP/KmrATPMg)*(F6P/KirF6P)+(1+(ATP-MgATP)/KirATP)*(F6P/KirF6P)+FDP/KirFDP+MgADP/KmrADP*(FDP/KirFDP)+KmrF6P/KirF6P*(MgATP/KmrATPMg)*(FDP/KirFDP)+Wr*(KmrF6P/KirF6P)*(MgADP/KirADP)*(MgATP/KmrATPMg)*(FDP/KmrFDP))/(1+L0*((1+KmtFDP/KitFDP*(MgADP/KmtADP)+KmtF6P/KitF6P*(MgATP/KmtATPMg)+KmtFDP/KitFDP*(MgADP/KmtADP)*(F6P/KitF6P)+MgATP/KmtATPMg*(F6P/KitF6P)+MgADP/KitADP*(MgATP/KmtATPMg)*(F6P/KitF6P)+(1+(ATP-MgATP)/KitATP)*(F6P/KitF6P)+FDP/KitFDP+MgADP/KmtADP*(FDP/KitFDP)+KmtF6P/KitF6P*(MgATP/KmtATPMg)*(FDP/KitFDP)+Wt*(KmtF6P/KitF6P)*(MgADP/KitADP)*(MgATP/KmtATPMg)*(FDP/KmtFDP))*(1+MgADP/KeftADP+PEP/KeftPEP+MgADP/KeftADP*(PEP/KeftPEP))/((1+KmrFDP/KirFDP*(MgADP/KmrADP)+KmrF6P*MgATP/(KirF6P*KmrATPMg)+KmrFDP/KirFDP*(MgADP/KmrADP)*(F6P/KirF6P)+MgATP/KmrATPMg*(F6P/KirF6P)+MgADP/KirADP*(MgATP/KmrATPMg)*(F6P/KirF6P)+(1+(ATP-MgATP)/KirATP)*(F6P/KirF6P)+FDP/KirFDP+MgADP/KmrADP*(FDP/KirFDP)+KmrF6P/KirF6P*(MgATP/KmrATPMg)*(FDP/KirFDP)+Wr*(KmrF6P/KirF6P)*(MgADP/KirADP)*(MgATP/KmrATPMg)*(FDP/KmrFDP))*(1+MgADP/KefrADP+PEP/KefrPEP+MgADP/KefrADP*(PEP/KefrPEP))))^n)
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_475" name="ATP" order="0" role="substrate"/>
        <ParameterDescription key="FunctionParameter_474" name="F6P" order="1" role="substrate"/>
        <ParameterDescription key="FunctionParameter_503" name="FDP" order="2" role="product"/>
        <ParameterDescription key="FunctionParameter_502" name="KefrADP" order="3" role="constant"/>
        <ParameterDescription key="FunctionParameter_501" name="KefrPEP" order="4" role="constant"/>
        <ParameterDescription key="FunctionParameter_500" name="KeftADP" order="5" role="constant"/>
        <ParameterDescription key="FunctionParameter_476" name="KeftPEP" order="6" role="constant"/>
        <ParameterDescription key="FunctionParameter_514" name="Keq" order="7" role="constant"/>
        <ParameterDescription key="FunctionParameter_515" name="KirADP" order="8" role="constant"/>
        <ParameterDescription key="FunctionParameter_516" name="KirATP" order="9" role="constant"/>
        <ParameterDescription key="FunctionParameter_517" name="KirF6P" order="10" role="constant"/>
        <ParameterDescription key="FunctionParameter_518" name="KirFDP" order="11" role="constant"/>
        <ParameterDescription key="FunctionParameter_519" name="KitADP" order="12" role="constant"/>
        <ParameterDescription key="FunctionParameter_520" name="KitATP" order="13" role="constant"/>
        <ParameterDescription key="FunctionParameter_521" name="KitF6P" order="14" role="constant"/>
        <ParameterDescription key="FunctionParameter_522" name="KitFDP" order="15" role="constant"/>
        <ParameterDescription key="FunctionParameter_523" name="KmrADP" order="16" role="constant"/>
        <ParameterDescription key="FunctionParameter_524" name="KmrATPMg" order="17" role="constant"/>
        <ParameterDescription key="FunctionParameter_525" name="KmrF6P" order="18" role="constant"/>
        <ParameterDescription key="FunctionParameter_526" name="KmrFDP" order="19" role="constant"/>
        <ParameterDescription key="FunctionParameter_527" name="KmtADP" order="20" role="constant"/>
        <ParameterDescription key="FunctionParameter_528" name="KmtATPMg" order="21" role="constant"/>
        <ParameterDescription key="FunctionParameter_529" name="KmtF6P" order="22" role="constant"/>
        <ParameterDescription key="FunctionParameter_530" name="KmtFDP" order="23" role="constant"/>
        <ParameterDescription key="FunctionParameter_531" name="L0" order="24" role="constant"/>
        <ParameterDescription key="FunctionParameter_532" name="MgADP" order="25" role="modifier"/>
        <ParameterDescription key="FunctionParameter_533" name="MgATP" order="26" role="modifier"/>
        <ParameterDescription key="FunctionParameter_534" name="PEP" order="27" role="modifier"/>
        <ParameterDescription key="FunctionParameter_535" name="Vmax" order="28" role="constant"/>
        <ParameterDescription key="FunctionParameter_536" name="Wr" order="29" role="constant"/>
        <ParameterDescription key="FunctionParameter_537" name="Wt" order="30" role="constant"/>
        <ParameterDescription key="FunctionParameter_538" name="n" order="31" role="constant"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_55" name="Function for PYK" type="UserDefined" reversible="false">
      <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Function_55">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-26T08:51:41Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

      </MiriamAnnotation>
      <Expression>
        Vmax*n*PEP*MgADP/(KirPEP*KmrADPMg)/(1+KmrPEP/KirPEP*(MgADP/KmrADPMg)+MgATP/KirATP+MgADP/KmrADPMg*(PEP/KirPEP)+KmrADPMg/KmrADPMg*(1+(ADP-MgADP)/KirADP)*(PEP/KirPEP)+PYR/KirPYR+MgATP/KirPyrATP*(PYR/KirPYR))/(1+L0*((1+KmtPEP/KitPEP*(MgADP/KmtADPMg)+MgATP/KitATP+MgADP*PEP/(KitPEP*KmtADPMg)+(1+(ADP-MgADP)/KitADP)*(PEP/KitPEP)+PYR/KitPYR+MgATP/KitPyrATP*(PYR/KitPYR))*(1+SUCCOA/KeftSUCCOA+MgATP*SUCCOA/(KeftATP*KeftSUCCOA))/((1+KmrPEP/KirPEP*(MgADP/KmrADPMg)+MgATP/KirATP+MgADP/KmrADPMg*(PEP/KirPEP)+(1+(ADP-MgADP)/KirADP)*(PEP/KirPEP)+PYR/KirPYR+MgATP/KirPyrATP*(PYR/KirPYR))*(1+FDP/KefrFDP+G6P/KefrG6P+GL6P/KefrGL6P+R5P/KefrR5P+RU5P/KefrRU5P+S7P/KefrS7P+X5P/KefrX5P)))^n)
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_570" name="ADP" order="0" role="substrate"/>
        <ParameterDescription key="FunctionParameter_569" name="FDP" order="1" role="modifier"/>
        <ParameterDescription key="FunctionParameter_568" name="G6P" order="2" role="modifier"/>
        <ParameterDescription key="FunctionParameter_567" name="GL6P" order="3" role="modifier"/>
        <ParameterDescription key="FunctionParameter_566" name="KefrFDP" order="4" role="constant"/>
        <ParameterDescription key="FunctionParameter_565" name="KefrG6P" order="5" role="constant"/>
        <ParameterDescription key="FunctionParameter_564" name="KefrGL6P" order="6" role="constant"/>
        <ParameterDescription key="FunctionParameter_563" name="KefrR5P" order="7" role="constant"/>
        <ParameterDescription key="FunctionParameter_562" name="KefrRU5P" order="8" role="constant"/>
        <ParameterDescription key="FunctionParameter_561" name="KefrS7P" order="9" role="constant"/>
        <ParameterDescription key="FunctionParameter_560" name="KefrX5P" order="10" role="constant"/>
        <ParameterDescription key="FunctionParameter_559" name="KeftATP" order="11" role="constant"/>
        <ParameterDescription key="FunctionParameter_558" name="KeftSUCCOA" order="12" role="constant"/>
        <ParameterDescription key="FunctionParameter_557" name="KirADP" order="13" role="constant"/>
        <ParameterDescription key="FunctionParameter_556" name="KirATP" order="14" role="constant"/>
        <ParameterDescription key="FunctionParameter_555" name="KirPEP" order="15" role="constant"/>
        <ParameterDescription key="FunctionParameter_554" name="KirPYR" order="16" role="constant"/>
        <ParameterDescription key="FunctionParameter_553" name="KirPyrATP" order="17" role="constant"/>
        <ParameterDescription key="FunctionParameter_552" name="KitADP" order="18" role="constant"/>
        <ParameterDescription key="FunctionParameter_551" name="KitATP" order="19" role="constant"/>
        <ParameterDescription key="FunctionParameter_550" name="KitPEP" order="20" role="constant"/>
        <ParameterDescription key="FunctionParameter_549" name="KitPYR" order="21" role="constant"/>
        <ParameterDescription key="FunctionParameter_548" name="KitPyrATP" order="22" role="constant"/>
        <ParameterDescription key="FunctionParameter_547" name="KmrADPMg" order="23" role="constant"/>
        <ParameterDescription key="FunctionParameter_546" name="KmrPEP" order="24" role="constant"/>
        <ParameterDescription key="FunctionParameter_545" name="KmtADPMg" order="25" role="constant"/>
        <ParameterDescription key="FunctionParameter_544" name="KmtPEP" order="26" role="constant"/>
        <ParameterDescription key="FunctionParameter_543" name="L0" order="27" role="constant"/>
        <ParameterDescription key="FunctionParameter_542" name="MgADP" order="28" role="modifier"/>
        <ParameterDescription key="FunctionParameter_541" name="MgATP" order="29" role="modifier"/>
        <ParameterDescription key="FunctionParameter_540" name="PEP" order="30" role="substrate"/>
        <ParameterDescription key="FunctionParameter_539" name="PYR" order="31" role="product"/>
        <ParameterDescription key="FunctionParameter_571" name="R5P" order="32" role="modifier"/>
        <ParameterDescription key="FunctionParameter_572" name="RU5P" order="33" role="modifier"/>
        <ParameterDescription key="FunctionParameter_573" name="S7P" order="34" role="modifier"/>
        <ParameterDescription key="FunctionParameter_574" name="SUCCOA" order="35" role="modifier"/>
        <ParameterDescription key="FunctionParameter_575" name="Vmax" order="36" role="constant"/>
        <ParameterDescription key="FunctionParameter_576" name="X5P" order="37" role="modifier"/>
        <ParameterDescription key="FunctionParameter_577" name="n" order="38" role="constant"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_56" name="Function for TPI" type="UserDefined" reversible="true">
      <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Function_56">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2019-03-08T16:05:29Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

      </MiriamAnnotation>
      <Expression>
        Vmax*(DAP-GAP/Keq)/KmDAP/(1+DAP/KmDAP+GAP/KmGAP)
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_616" name="DAP" order="0" role="substrate"/>
        <ParameterDescription key="FunctionParameter_615" name="GAP" order="1" role="product"/>
        <ParameterDescription key="FunctionParameter_614" name="Keq" order="2" role="constant"/>
        <ParameterDescription key="FunctionParameter_613" name="KmDAP" order="3" role="constant"/>
        <ParameterDescription key="FunctionParameter_612" name="KmGAP" order="4" role="constant"/>
        <ParameterDescription key="FunctionParameter_611" name="Vmax" order="5" role="constant"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_61" name="Function for F6P_GAP_TAL" type="UserDefined" reversible="true">
      <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Function_61">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2014-12-29T22:00:25Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

      </MiriamAnnotation>
      <Expression>
        kcat*(GAP*talC3-F6P*tal/Keq)
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_617" name="F6P" order="0" role="product"/>
        <ParameterDescription key="FunctionParameter_578" name="GAP" order="1" role="substrate"/>
        <ParameterDescription key="FunctionParameter_579" name="Keq" order="2" role="constant"/>
        <ParameterDescription key="FunctionParameter_580" name="kcat" order="3" role="constant"/>
        <ParameterDescription key="FunctionParameter_581" name="tal" order="4" role="product"/>
        <ParameterDescription key="FunctionParameter_582" name="talC3" order="5" role="substrate"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_62" name="Function for GPM" type="UserDefined" reversible="true">
      <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Function_62">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-13T10:20:43Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

      </MiriamAnnotation>
      <Expression>
        Vmax*(PGA3-PGA2/Keq)/KmPGA3/(1+PGA3/KmPGA3+PGA2/KmPGA2)
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_623" name="Keq" order="0" role="constant"/>
        <ParameterDescription key="FunctionParameter_622" name="KmPGA2" order="1" role="constant"/>
        <ParameterDescription key="FunctionParameter_621" name="KmPGA3" order="2" role="constant"/>
        <ParameterDescription key="FunctionParameter_620" name="PGA2" order="3" role="product"/>
        <ParameterDescription key="FunctionParameter_619" name="PGA3" order="4" role="substrate"/>
        <ParameterDescription key="FunctionParameter_618" name="Vmax" order="5" role="constant"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_70" name="Function for PTS_0" type="UserDefined" reversible="true">
      <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Function_70">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2015-01-22T14:03:00Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

      </MiriamAnnotation>
      <Expression>
        kF*ei*PEP^2/(KmPEP^2+PEP^2)-kR*eiP*PYR^2/(KmPYR^2+PYR^2)
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_701" name="KmPEP" order="0" role="constant"/>
        <ParameterDescription key="FunctionParameter_700" name="KmPYR" order="1" role="constant"/>
        <ParameterDescription key="FunctionParameter_699" name="PEP" order="2" role="substrate"/>
        <ParameterDescription key="FunctionParameter_698" name="PYR" order="3" role="product"/>
        <ParameterDescription key="FunctionParameter_697" name="ei" order="4" role="substrate"/>
        <ParameterDescription key="FunctionParameter_696" name="eiP" order="5" role="product"/>
        <ParameterDescription key="FunctionParameter_695" name="kF" order="6" role="constant"/>
        <ParameterDescription key="FunctionParameter_694" name="kR" order="7" role="constant"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_73" name="Function for GLC_FEED" type="UserDefined" reversible="false">
      <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Function_73">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2014-12-29T21:58:51Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

      </MiriamAnnotation>
      <Expression>
        cell*FEED/extracellular
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_719" name="FEED" order="0" role="constant"/>
        <ParameterDescription key="FunctionParameter_718" name="cell" order="1" role="volume"/>
        <ParameterDescription key="FunctionParameter_708" name="extracellular" order="2" role="volume"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_85" name="Function for PGI_1" type="UserDefined" reversible="true">
      <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Function_85">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2019-03-08T16:02:36Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

      </MiriamAnnotation>
      <Expression>
        Vmax*(G6P-F6P/Keq)/KmG6P/(1+F6P/KmF6P+G6P/KmG6P+PEP/KmPEP+PGN/KmPGN)
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_838" name="F6P" order="0" role="product"/>
        <ParameterDescription key="FunctionParameter_837" name="G6P" order="1" role="substrate"/>
        <ParameterDescription key="FunctionParameter_836" name="Keq" order="2" role="constant"/>
        <ParameterDescription key="FunctionParameter_835" name="KmF6P" order="3" role="constant"/>
        <ParameterDescription key="FunctionParameter_834" name="KmG6P" order="4" role="constant"/>
        <ParameterDescription key="FunctionParameter_833" name="KmPEP" order="5" role="constant"/>
        <ParameterDescription key="FunctionParameter_832" name="KmPGN" order="6" role="constant"/>
        <ParameterDescription key="FunctionParameter_831" name="PEP" order="7" role="modifier"/>
        <ParameterDescription key="FunctionParameter_839" name="PGN" order="8" role="modifier"/>
        <ParameterDescription key="FunctionParameter_840" name="Vmax" order="9" role="constant"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_92" name="Function for PTS_4_1" type="UserDefined" reversible="true">
      <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Function_92">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2019-03-08T16:04:21Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

      </MiriamAnnotation>
      <Expression>
        cell*(kF*eiicbP*GLCp/(KmGLC+GLCp)-kR*eiicb*G6P/(KmG6P+G6P))
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_897" name="G6P" order="0" role="product"/>
        <ParameterDescription key="FunctionParameter_896" name="GLCp" order="1" role="substrate"/>
        <ParameterDescription key="FunctionParameter_895" name="KmG6P" order="2" role="constant"/>
        <ParameterDescription key="FunctionParameter_894" name="KmGLC" order="3" role="constant"/>
        <ParameterDescription key="FunctionParameter_893" name="cell" order="4" role="volume"/>
        <ParameterDescription key="FunctionParameter_892" name="eiicb" order="5" role="product"/>
        <ParameterDescription key="FunctionParameter_891" name="eiicbP" order="6" role="substrate"/>
        <ParameterDescription key="FunctionParameter_898" name="kF" order="7" role="constant"/>
        <ParameterDescription key="FunctionParameter_899" name="kR" order="8" role="constant"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_97" name="Function for ATP_MAINTENANCE_1" type="UserDefined" reversible="true">
      <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Function_97">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2019-03-08T15:57:33Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

      </MiriamAnnotation>
      <Expression>
        Vmax*(ATP-ADP*P/Keq)
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_976" name="ADP" order="0" role="product"/>
        <ParameterDescription key="FunctionParameter_975" name="ATP" order="1" role="substrate"/>
        <ParameterDescription key="FunctionParameter_974" name="Keq" order="2" role="constant"/>
        <ParameterDescription key="FunctionParameter_973" name="P" order="3" role="product"/>
        <ParameterDescription key="FunctionParameter_972" name="Vmax" order="4" role="constant"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_98" name="Function for XCH_RMM_1" type="UserDefined" reversible="true">
      <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Function_98">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2019-03-08T16:05:41Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

      </MiriamAnnotation>
      <Expression>
        Vmax*(GLCx/Km-GLCp/Km)/(1+GLCx/Km+GLCp/Km)
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_967" name="GLCp" order="0" role="product"/>
        <ParameterDescription key="FunctionParameter_968" name="GLCx" order="1" role="substrate"/>
        <ParameterDescription key="FunctionParameter_969" name="Km" order="2" role="constant"/>
        <ParameterDescription key="FunctionParameter_970" name="Vmax" order="3" role="constant"/>
      </ListOfParameterDescriptions>
    </Function>
    <Function key="Function_103" name="Function for GL6P_HYDROLYSIS_1" type="UserDefined" reversible="true">
      <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Function_103">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2019-03-08T15:58:38Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

      </MiriamAnnotation>
      <Expression>
        KGl6Phydrol*(GL6P-PGN/KeqGl6Phydrol)
      </Expression>
      <ListOfParameterDescriptions>
        <ParameterDescription key="FunctionParameter_980" name="GL6P" order="0" role="substrate"/>
        <ParameterDescription key="FunctionParameter_979" name="KGl6Phydrol" order="1" role="constant"/>
        <ParameterDescription key="FunctionParameter_978" name="KeqGl6Phydrol" order="2" role="constant"/>
        <ParameterDescription key="FunctionParameter_977" name="PGN" order="3" role="product"/>
      </ListOfParameterDescriptions>
    </Function>
  </ListOfFunctions>
  <Model key="Model_1" name="Millard2016 - E. coli central carbon and energy metabolism" simulationType="time" timeUnit="s" volumeUnit="l" areaUnit="m²" lengthUnit="m" quantityUnit="mmol" type="deterministic" avogadroConstant="6.0221408570000002e+23">
    <MiriamAnnotation>
<rdf:RDF
   xmlns:CopasiMT="http://www.copasi.org/RDF/MiriamTerms#"
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"
   xmlns:vCard="http://www.w3.org/2001/vcard-rdf/3.0#">
  <rdf:Description rdf:about="#Model_1">
    <dcterms:bibliographicCitation>
      <rdf:Bag>
        <rdf:li>
          <rdf:Description>
            <CopasiMT:isDescribedBy rdf:resource="http://identifiers.org/doi/10.1371/journal.pcbi.1005396"/>
          </rdf:Description>
        </rdf:li>
      </rdf:Bag>
    </dcterms:bibliographicCitation>
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2015-05-10T22:00:37Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
    <dcterms:creator>
      <rdf:Bag>
        <rdf:li>
          <rdf:Description>
            <vCard:EMAIL>millard@insa-toulouse.fr</vCard:EMAIL>
            <vCard:N>
              <rdf:Description>
                <vCard:Family>Millard</vCard:Family>
                <vCard:Given>Pierre</vCard:Given>
              </rdf:Description>
            </vCard:N>
            <vCard:ORG>
              <rdf:Description>
                <vCard:Orgname>INRA</vCard:Orgname>
              </rdf:Description>
            </vCard:ORG>
          </rdf:Description>
        </rdf:li>
        <rdf:li>
          <rdf:Description>
            <vCard:EMAIL>kieran.smallbone@manchester.ac.uk</vCard:EMAIL>
            <vCard:N>
              <rdf:Description>
                <vCard:Family>Smallbone</vCard:Family>
                <vCard:Given>Kieran</vCard:Given>
              </rdf:Description>
            </vCard:N>
            <vCard:ORG>
              <rdf:Description>
                <vCard:Orgname>University of Manchester</vCard:Orgname>
              </rdf:Description>
            </vCard:ORG>
          </rdf:Description>
        </rdf:li>
        <rdf:li>
          <rdf:Description>
            <vCard:EMAIL>pedro.mendes@manchester.ac.uk</vCard:EMAIL>
            <vCard:N>
              <rdf:Description>
                <vCard:Family>Mendes</vCard:Family>
                <vCard:Given>Pedro</vCard:Given>
              </rdf:Description>
            </vCard:N>
            <vCard:ORG>
              <rdf:Description>
                <vCard:Orgname>University of Manchester</vCard:Orgname>
              </rdf:Description>
            </vCard:ORG>
          </rdf:Description>
        </rdf:li>
      </rdf:Bag>
    </dcterms:creator>
    <dcterms:modified>
      <rdf:Description>
        <dcterms:W3CDTF>2015-05-10T22:00:37Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:modified>
    <dcterms:modified>
      <rdf:Description>
        <dcterms:W3CDTF>2017-02-22T10:26:11Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:modified>
    <CopasiMT:is>
      <rdf:Bag>
        <rdf:li rdf:resource="http://identifiers.org/biomodels.db/MODEL1505110000"/>
      </rdf:Bag>
    </CopasiMT:is>
  </rdf:Description>
</rdf:RDF>

    </MiriamAnnotation>
    <Comment>
      <body xmlns="http://www.w3.org/1999/xhtml">
    <div class="dc:title">Millard2016 - E. coli central carbon and
energy metabolism</div>
    <div class="dc:bibliographicCitation">
      <p>This model is described in the article:</p>
      <div class="bibo:title">
        <a href="http://identifiers.org/doi/10.1371/journal.pcbi.1005396" title="Access to this publication">Metabolic regulation is
    sufficient for global and robust coordination of glucose
    uptake, catabolism, energy production and growth in Escherichia
    coli</a>
      </div>
      <div class="bibo:authorList">Pierre Millard , Kieran Smallbone,
  Pedro Mendes</div>
      <div class="bibo:Journal">PLoS ONE</div>
      <p>Abstract:</p>
      <div class="bibo:abstract">
        <p>The metabolism of microorganisms is regulated through two
    main mechanisms: changes of enzyme capacities as a consequence
    of gene expression modulation (“hierarchical
    control”) and changes of enzyme activities through
    metabolite-enzyme interactions. An increasing body of evidence
    indicates that hierarchical control is insufficient to explain
    metabolic behaviors, but the system-wide impact of metabolic
    regulation remains largely uncharacterized. To clarify its
    role, we developed and validated a detailed kinetic model of
    Escherichia coli central metabolism that links growth to
    environment. Metabolic control analyses confirm that the
    control is widely distributed across the network and highlight
    strong interconnections between all the pathways. Exploration
    of the model solution space reveals that several robust
    properties emerge from metabolic regulation, from the molecular
    level (e.g. homeostasis of total metabolite pool) to the
    overall cellular physiology (e.g. coordination of carbon
    uptake, catabolism, energy and redox production, and growth),
    while allowing a large degree of flexibility at most individual
    metabolic steps. These properties have important physiological
    implications for E. coli and significantly expand the
    self-regulating capacities of its metabolism.</p>
      </div>
    </div>
    <div class="dc:publisher">
      <p>This model is hosted on 
  <a href="http://www.ebi.ac.uk/biomodels/">BioModels Database</a>
  and identified by: 
  <a href="http://identifiers.org/biomodels.db/MODEL1505110000">MODEL1505110000</a>.</p>
      <p>To cite BioModels Database, please use: 
  <a href="http://identifiers.org/pubmed/20587024" title="Latest BioModels Database publication">BioModels Database:
  An enhanced, curated and annotated resource for published
  quantitative kinetic models</a>.</p>
    </div>
    <div class="dc:license">
      <p>To the extent possible under law, all copyright and related or
  neighbouring rights to this encoded model have been dedicated to
  the public domain worldwide. Please refer to 
  <a href="http://creativecommons.org/publicdomain/zero/1.0/" title="Access to: CC0 1.0 Universal (CC0 1.0), Public Domain Dedication">CC0
  Public Domain Dedication</a> for more information.</p>
    </div>
  </body>
    </Comment>
    <ListOfCompartments>
      <Compartment key="Compartment_0" name="cell_cytoplasm" simulationType="fixed" dimensionality="3" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Compartment_0">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-26T16:22:10Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

        </MiriamAnnotation>
      </Compartment>
      <Compartment key="Compartment_1" name="extracellular" simulationType="fixed" dimensionality="3" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Compartment_1">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-26T16:22:10Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

        </MiriamAnnotation>
      </Compartment>
      <Compartment key="Compartment_2" name="cell_periplasm" simulationType="fixed" dimensionality="3" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Compartment_2">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-16T10:08:42Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

        </MiriamAnnotation>
      </Compartment>
    </ListOfCompartments>
    <ListOfMetabolites>
      <Metabolite key="Metabolite_3" name="AKG" simulationType="reactions" compartment="Compartment_0" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_3">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-17T08:08:54Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_4" name="BPG" simulationType="reactions" compartment="Compartment_0" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_4">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-01T14:57:13Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_6" name="DAP" simulationType="reactions" compartment="Compartment_0" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_6">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-15T10:03:10Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_8" name="F6P" simulationType="reactions" compartment="Compartment_0" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_8">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-26T16:25:24Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_9" name="FDP" simulationType="reactions" compartment="Compartment_0" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_9">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-30T23:49:14Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_11" name="G6P" simulationType="reactions" compartment="Compartment_0" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_11">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2019-03-06T17:04:27Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_12" name="GAP" simulationType="reactions" compartment="Compartment_0" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_12">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2015-01-03T00:12:37Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_13" name="GL6P" simulationType="reactions" compartment="Compartment_0" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_13">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-30T23:53:03Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_18" name="NAD" simulationType="reactions" compartment="Compartment_0" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_18">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-25T23:44:07Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_19" name="NADH" simulationType="reactions" compartment="Compartment_0" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_19">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-25T23:44:56Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_22" name="OAA" simulationType="reactions" compartment="Compartment_0" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_22">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-25T19:45:19Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_23" name="PEP" simulationType="reactions" compartment="Compartment_0" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_23">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-24T16:37:45Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_24" name="PGA2" simulationType="reactions" compartment="Compartment_0" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_24">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2019-03-08T15:55:12Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_25" name="PGA3" simulationType="reactions" compartment="Compartment_0" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_25">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2019-03-08T15:55:14Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_26" name="PGN" simulationType="reactions" compartment="Compartment_0" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_26">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2019-03-08T15:49:00Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_27" name="PYR" simulationType="reactions" compartment="Compartment_0" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_27">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-31T12:41:19Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_30" name="R5P" simulationType="reactions" compartment="Compartment_0" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_30">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2019-03-08T15:55:32Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_31" name="RU5P" simulationType="reactions" compartment="Compartment_0" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_31">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2019-03-08T15:55:35Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_32" name="S7P" simulationType="reactions" compartment="Compartment_0" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_32">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-24T10:25:24Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_34" name="SUCCOA" simulationType="reactions" compartment="Compartment_0" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_34">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-14T13:47:22Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_35" name="X5P" simulationType="reactions" compartment="Compartment_0" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_35">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2015-01-22T20:51:59Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_36" name="ei" simulationType="reactions" compartment="Compartment_0" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_36">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-15T14:56:48Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_37" name="eiP" simulationType="reactions" compartment="Compartment_0" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_37">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2019-03-08T15:53:27Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_38" name="eiia" simulationType="reactions" compartment="Compartment_0" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_38">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-16T09:44:24Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_39" name="eiiaP" simulationType="reactions" compartment="Compartment_0" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_39">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2019-03-06T17:04:23Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_40" name="eiicb" simulationType="reactions" compartment="Compartment_0" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_40">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2019-03-08T15:53:26Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_41" name="eiicbP" simulationType="reactions" compartment="Compartment_0" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_41">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-15T11:48:22Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_42" name="hpr" simulationType="reactions" compartment="Compartment_0" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_42">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-25T19:47:23Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_43" name="hprP" simulationType="reactions" compartment="Compartment_0" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_43">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-10T10:15:23Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_46" name="tal" simulationType="reactions" compartment="Compartment_0" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_46">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-01T10:37:13Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_47" name="talC3" simulationType="reactions" compartment="Compartment_0" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_47">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-01T10:37:13Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_50" name="ADP" simulationType="reactions" compartment="Compartment_0" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_50">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-25T19:47:02Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_51" name="AMP" simulationType="reactions" compartment="Compartment_0" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_51">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-25T19:47:05Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_52" name="ATP" simulationType="reactions" compartment="Compartment_0" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_52">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-25T19:47:08Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_56" name="P" simulationType="reactions" compartment="Compartment_0" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_56">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-25T19:45:13Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_57" name="MG" simulationType="fixed" compartment="Compartment_0" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_57">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-30T23:50:03Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_58" name="MgADP" simulationType="assignment" compartment="Compartment_0" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_58">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-01T00:28:29Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

        </MiriamAnnotation>
        <Expression>
          &lt;CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[MG],Reference=Concentration>*&lt;CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[ADP],Reference=Concentration>/(&lt;CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Values[KdADPMg],Reference=Value>+&lt;CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[MG],Reference=Concentration>)
        </Expression>
      </Metabolite>
      <Metabolite key="Metabolite_59" name="MgATP" simulationType="assignment" compartment="Compartment_0" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_59">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2015-01-23T11:39:48Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

        </MiriamAnnotation>
        <Expression>
          &lt;CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[MG],Reference=Concentration>*&lt;CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[ATP],Reference=Concentration>/(&lt;CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Values[KdATPMg],Reference=Value>+&lt;CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[MG],Reference=Concentration>)
        </Expression>
      </Metabolite>
      <Metabolite key="Metabolite_60" name="MgFDP" simulationType="assignment" compartment="Compartment_0" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_60">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-01T00:28:01Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

        </MiriamAnnotation>
        <Expression>
          &lt;CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[MG],Reference=Concentration>*&lt;CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[FDP],Reference=Concentration>/(&lt;CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Values[KdFDPMg],Reference=Value>+&lt;CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[MG],Reference=Concentration>)
        </Expression>
      </Metabolite>
      <Metabolite key="Metabolite_70" name="GLCx" simulationType="reactions" compartment="Compartment_1" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_70">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-30T23:49:09Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

        </MiriamAnnotation>
      </Metabolite>
      <Metabolite key="Metabolite_74" name="GLCp" simulationType="reactions" compartment="Compartment_2" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Metabolite_74">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-23T16:11:40Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

        </MiriamAnnotation>
      </Metabolite>
    </ListOfMetabolites>
    <ListOfModelValues>
      <ModelValue key="ModelValue_0" name="FEED" simulationType="fixed" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#ModelValue_0">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-01T10:22:21Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </ModelValue>
      <ModelValue key="ModelValue_1" name="KdADPMg" simulationType="fixed" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#ModelValue_1">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-30T19:54:50Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </ModelValue>
      <ModelValue key="ModelValue_2" name="KdATPMg" simulationType="fixed" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#ModelValue_2">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-20T18:47:20Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </ModelValue>
      <ModelValue key="ModelValue_3" name="KdFDPMg" simulationType="fixed" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#ModelValue_3">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-20T18:47:20Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </ModelValue>
      <ModelValue key="ModelValue_4" name="KmICIT_ACN" simulationType="fixed" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#ModelValue_4">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-13T02:18:45Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </ModelValue>
      <ModelValue key="ModelValue_5" name="KmCIT_ACN" simulationType="fixed" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#ModelValue_5">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-13T02:18:59Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </ModelValue>
      <ModelValue key="ModelValue_6" name="KmACO_ACN" simulationType="fixed" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#ModelValue_6">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-13T02:19:05Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </ModelValue>
      <ModelValue key="ModelValue_7" name="KeqNDH" simulationType="fixed" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF xmlns:dcterms="http://purl.org/dc/terms/" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#ModelValue_7">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-25T10:29:33Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>
        </MiriamAnnotation>
      </ModelValue>
    </ListOfModelValues>
    <ListOfReactions>
      <Reaction key="Reaction_0" name="PGI" reversible="true" fast="false" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_0">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-30T23:08:22Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_11" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_8" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfModifiers>
          <Modifier metabolite="Metabolite_23" stoichiometry="1"/>
          <Modifier metabolite="Metabolite_26" stoichiometry="1"/>
        </ListOfModifiers>
        <ListOfConstants>
          <Constant key="Parameter_4999" name="Keq" value="0.36"/>
          <Constant key="Parameter_4998" name="KmF6P" value="0.147"/>
          <Constant key="Parameter_4997" name="KmG6P" value="0.28"/>
          <Constant key="Parameter_4996" name="KmPEP" value="1.999"/>
          <Constant key="Parameter_4995" name="Vmax" value="2.32456"/>
          <Constant key="Parameter_4994" name="KmPGN" value="0.515958"/>
        </ListOfConstants>
        <KineticLaw function="Function_85" unitType="Default" scalingCompartment="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm]">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_838">
              <SourceParameter reference="Metabolite_8"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_837">
              <SourceParameter reference="Metabolite_11"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_836">
              <SourceParameter reference="Parameter_4999"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_835">
              <SourceParameter reference="Parameter_4998"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_834">
              <SourceParameter reference="Parameter_4997"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_833">
              <SourceParameter reference="Parameter_4996"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_832">
              <SourceParameter reference="Parameter_4994"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_831">
              <SourceParameter reference="Metabolite_23"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_839">
              <SourceParameter reference="Metabolite_26"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_840">
              <SourceParameter reference="Parameter_4995"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_1" name="PFK" reversible="true" fast="false" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_1">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-07T10:52:04Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_52" stoichiometry="1"/>
          <Substrate metabolite="Metabolite_8" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_50" stoichiometry="1"/>
          <Product metabolite="Metabolite_9" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfModifiers>
          <Modifier metabolite="Metabolite_23" stoichiometry="1"/>
          <Modifier metabolite="Metabolite_58" stoichiometry="1"/>
          <Modifier metabolite="Metabolite_59" stoichiometry="1"/>
        </ListOfModifiers>
        <ListOfConstants>
          <Constant key="Parameter_4993" name="KefrADP" value="0.0735264"/>
          <Constant key="Parameter_4992" name="KefrPEP" value="19.98"/>
          <Constant key="Parameter_4991" name="KeftADP" value="9.009"/>
          <Constant key="Parameter_4990" name="KeftPEP" value="0.26026"/>
          <Constant key="Parameter_4989" name="Keq" value="1998"/>
          <Constant key="Parameter_4988" name="KirADP" value="54.945"/>
          <Constant key="Parameter_4987" name="KirATP" value="2.4975e-05"/>
          <Constant key="Parameter_4986" name="KirF6P" value="1.84615"/>
          <Constant key="Parameter_4985" name="KirFDP" value="0.045954"/>
          <Constant key="Parameter_4984" name="KitADP" value="80.08"/>
          <Constant key="Parameter_4983" name="KitATP" value="0.014014"/>
          <Constant key="Parameter_4982" name="KitF6P" value="0.00856856"/>
          <Constant key="Parameter_4981" name="KitFDP" value="50.5505"/>
          <Constant key="Parameter_4980" name="KmrADP" value="0.690009"/>
          <Constant key="Parameter_4979" name="KmrATPMg" value="8.12187e-05"/>
          <Constant key="Parameter_4978" name="KmrF6P" value="2.05205e-05"/>
          <Constant key="Parameter_4977" name="KmrFDP" value="10.01"/>
          <Constant key="Parameter_4976" name="KmtADP" value="2.002"/>
          <Constant key="Parameter_4975" name="KmtATPMg" value="3.34334"/>
          <Constant key="Parameter_4974" name="KmtF6P" value="32.967"/>
          <Constant key="Parameter_4973" name="KmtFDP" value="9.99"/>
          <Constant key="Parameter_4972" name="L0" value="14.0851"/>
          <Constant key="Parameter_4971" name="Vmax" value="0.185253"/>
          <Constant key="Parameter_4970" name="Wr" value="0.0237041"/>
          <Constant key="Parameter_4969" name="Wt" value="0.146735"/>
          <Constant key="Parameter_4968" name="n" value="4"/>
        </ListOfConstants>
        <KineticLaw function="Function_54" unitType="Default" scalingCompartment="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm]">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_475">
              <SourceParameter reference="Metabolite_52"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_474">
              <SourceParameter reference="Metabolite_8"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_503">
              <SourceParameter reference="Metabolite_9"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_502">
              <SourceParameter reference="Parameter_4993"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_501">
              <SourceParameter reference="Parameter_4992"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_500">
              <SourceParameter reference="Parameter_4991"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_476">
              <SourceParameter reference="Parameter_4990"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_514">
              <SourceParameter reference="Parameter_4989"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_515">
              <SourceParameter reference="Parameter_4988"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_516">
              <SourceParameter reference="Parameter_4987"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_517">
              <SourceParameter reference="Parameter_4986"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_518">
              <SourceParameter reference="Parameter_4985"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_519">
              <SourceParameter reference="Parameter_4984"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_520">
              <SourceParameter reference="Parameter_4983"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_521">
              <SourceParameter reference="Parameter_4982"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_522">
              <SourceParameter reference="Parameter_4981"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_523">
              <SourceParameter reference="Parameter_4980"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_524">
              <SourceParameter reference="Parameter_4979"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_525">
              <SourceParameter reference="Parameter_4978"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_526">
              <SourceParameter reference="Parameter_4977"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_527">
              <SourceParameter reference="Parameter_4976"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_528">
              <SourceParameter reference="Parameter_4975"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_529">
              <SourceParameter reference="Parameter_4974"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_530">
              <SourceParameter reference="Parameter_4973"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_531">
              <SourceParameter reference="Parameter_4972"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_532">
              <SourceParameter reference="Metabolite_58"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_533">
              <SourceParameter reference="Metabolite_59"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_534">
              <SourceParameter reference="Metabolite_23"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_535">
              <SourceParameter reference="Parameter_4971"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_536">
              <SourceParameter reference="Parameter_4970"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_537">
              <SourceParameter reference="Parameter_4969"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_538">
              <SourceParameter reference="Parameter_4968"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_2" name="FBA" reversible="true" fast="false" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_2">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-07T10:49:12Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_9" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_6" stoichiometry="1"/>
          <Product metabolite="Metabolite_12" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfModifiers>
          <Modifier metabolite="Metabolite_23" stoichiometry="1"/>
        </ListOfModifiers>
        <ListOfConstants>
          <Constant key="Parameter_4967" name="Keq" value="0.18981"/>
          <Constant key="Parameter_4966" name="KmDAP" value="0.13001"/>
          <Constant key="Parameter_4965" name="KmFDP" value="0.12012"/>
          <Constant key="Parameter_4964" name="KmGAP" value="0.13001"/>
          <Constant key="Parameter_4963" name="KmPEP" value="0.5"/>
          <Constant key="Parameter_4962" name="Vmax" value="21.6978"/>
        </ListOfConstants>
        <KineticLaw function="Function_52" unitType="Default" scalingCompartment="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm]">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_513">
              <SourceParameter reference="Metabolite_6"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_512">
              <SourceParameter reference="Metabolite_9"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_511">
              <SourceParameter reference="Metabolite_12"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_510">
              <SourceParameter reference="Parameter_4967"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_509">
              <SourceParameter reference="Parameter_4966"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_508">
              <SourceParameter reference="Parameter_4965"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_507">
              <SourceParameter reference="Parameter_4964"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_506">
              <SourceParameter reference="Parameter_4963"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_505">
              <SourceParameter reference="Metabolite_23"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_504">
              <SourceParameter reference="Parameter_4962"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_3" name="TPI" reversible="true" fast="false" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_3">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-07T10:49:41Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_6" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_12" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfConstants>
          <Constant key="Parameter_4961" name="Keq" value="0.270203"/>
          <Constant key="Parameter_4960" name="KmDAP" value="0.01"/>
          <Constant key="Parameter_4959" name="KmGAP" value="1.89301"/>
          <Constant key="Parameter_4958" name="Vmax" value="24.1843"/>
        </ListOfConstants>
        <KineticLaw function="Function_56" unitType="Default" scalingCompartment="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm]">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_616">
              <SourceParameter reference="Metabolite_6"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_615">
              <SourceParameter reference="Metabolite_12"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_614">
              <SourceParameter reference="Parameter_4961"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_613">
              <SourceParameter reference="Parameter_4960"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_612">
              <SourceParameter reference="Parameter_4959"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_611">
              <SourceParameter reference="Parameter_4958"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_4" name="GDH" reversible="true" fast="false" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_4">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-01T17:03:28Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_12" stoichiometry="1"/>
          <Substrate metabolite="Metabolite_18" stoichiometry="1"/>
          <Substrate metabolite="Metabolite_56" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_4" stoichiometry="1"/>
          <Product metabolite="Metabolite_19" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfConstants>
          <Constant key="Parameter_4957" name="Keq" value="20"/>
          <Constant key="Parameter_4956" name="KmBPG" value="0.2"/>
          <Constant key="Parameter_4955" name="KmGAP" value="2.47265"/>
          <Constant key="Parameter_4954" name="KmNAD" value="0.0110454"/>
          <Constant key="Parameter_4953" name="KmNADH" value="3.69797"/>
          <Constant key="Parameter_4952" name="KmP" value="0.017"/>
          <Constant key="Parameter_4951" name="Vmax" value="8.66573"/>
        </ListOfConstants>
        <KineticLaw function="Function_42" unitType="Default" scalingCompartment="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm]">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_298">
              <SourceParameter reference="Metabolite_4"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_299">
              <SourceParameter reference="Metabolite_12"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_300">
              <SourceParameter reference="Parameter_4957"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_301">
              <SourceParameter reference="Parameter_4956"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_302">
              <SourceParameter reference="Parameter_4955"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_303">
              <SourceParameter reference="Parameter_4954"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_304">
              <SourceParameter reference="Parameter_4953"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_305">
              <SourceParameter reference="Parameter_4952"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_306">
              <SourceParameter reference="Metabolite_18"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_307">
              <SourceParameter reference="Metabolite_19"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_297">
              <SourceParameter reference="Metabolite_56"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_296">
              <SourceParameter reference="Parameter_4951"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_5" name="PGK" reversible="true" fast="false" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_5">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-14T17:46:47Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_50" stoichiometry="1"/>
          <Substrate metabolite="Metabolite_4" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_52" stoichiometry="1"/>
          <Product metabolite="Metabolite_25" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfModifiers>
          <Modifier metabolite="Metabolite_58" stoichiometry="1"/>
          <Modifier metabolite="Metabolite_59" stoichiometry="1"/>
        </ListOfModifiers>
        <ListOfConstants>
          <Constant key="Parameter_4950" name="Keq" value="99.9925"/>
          <Constant key="Parameter_4949" name="KmADPMg" value="0.085416"/>
          <Constant key="Parameter_4948" name="KmATPMg" value="3.47737"/>
          <Constant key="Parameter_4947" name="KmBPG" value="0.0113296"/>
          <Constant key="Parameter_4946" name="KmPGA3" value="2.45722"/>
          <Constant key="Parameter_4945" name="Vmax" value="16.1089"/>
        </ListOfConstants>
        <KineticLaw function="Function_41" unitType="Default" scalingCompartment="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm]">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_317">
              <SourceParameter reference="Metabolite_4"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_316">
              <SourceParameter reference="Parameter_4950"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_315">
              <SourceParameter reference="Parameter_4949"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_314">
              <SourceParameter reference="Parameter_4948"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_313">
              <SourceParameter reference="Parameter_4947"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_312">
              <SourceParameter reference="Parameter_4946"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_311">
              <SourceParameter reference="Metabolite_58"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_310">
              <SourceParameter reference="Metabolite_59"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_309">
              <SourceParameter reference="Metabolite_25"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_308">
              <SourceParameter reference="Parameter_4945"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_6" name="GPM" reversible="true" fast="false" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_6">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-24T12:56:58Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_25" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_24" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfConstants>
          <Constant key="Parameter_4944" name="Keq" value="0.565818"/>
          <Constant key="Parameter_4943" name="KmPGA2" value="1.9153"/>
          <Constant key="Parameter_4942" name="KmPGA3" value="0.115"/>
          <Constant key="Parameter_4941" name="Vmax" value="10.9934"/>
        </ListOfConstants>
        <KineticLaw function="Function_62" unitType="Default" scalingCompartment="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm]">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_623">
              <SourceParameter reference="Parameter_4944"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_622">
              <SourceParameter reference="Parameter_4943"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_621">
              <SourceParameter reference="Parameter_4942"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_620">
              <SourceParameter reference="Metabolite_24"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_619">
              <SourceParameter reference="Metabolite_25"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_618">
              <SourceParameter reference="Parameter_4941"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_7" name="ENO" reversible="true" fast="false" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_7">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-25T20:00:11Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_24" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_23" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfConstants>
          <Constant key="Parameter_4940" name="Keq" value="3"/>
          <Constant key="Parameter_4939" name="KmPEP" value="0.1"/>
          <Constant key="Parameter_4937" name="KmPGA2" value="0.1"/>
          <Constant key="Parameter_4938" name="Vmax" value="11.7189"/>
        </ListOfConstants>
        <KineticLaw function="Function_44" unitType="Default" scalingCompartment="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm]">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_295">
              <SourceParameter reference="Parameter_4940"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_294">
              <SourceParameter reference="Parameter_4939"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_293">
              <SourceParameter reference="Parameter_4937"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_292">
              <SourceParameter reference="Metabolite_23"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_291">
              <SourceParameter reference="Metabolite_24"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_290">
              <SourceParameter reference="Parameter_4938"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_8" name="PYK" reversible="false" fast="false" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_8">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-08T11:11:05Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_50" stoichiometry="1"/>
          <Substrate metabolite="Metabolite_23" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_52" stoichiometry="1"/>
          <Product metabolite="Metabolite_27" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfModifiers>
          <Modifier metabolite="Metabolite_9" stoichiometry="1"/>
          <Modifier metabolite="Metabolite_11" stoichiometry="1"/>
          <Modifier metabolite="Metabolite_13" stoichiometry="1"/>
          <Modifier metabolite="Metabolite_58" stoichiometry="1"/>
          <Modifier metabolite="Metabolite_59" stoichiometry="1"/>
          <Modifier metabolite="Metabolite_30" stoichiometry="1"/>
          <Modifier metabolite="Metabolite_31" stoichiometry="1"/>
          <Modifier metabolite="Metabolite_32" stoichiometry="1"/>
          <Modifier metabolite="Metabolite_34" stoichiometry="1"/>
          <Modifier metabolite="Metabolite_35" stoichiometry="1"/>
        </ListOfModifiers>
        <ListOfConstants>
          <Constant key="Parameter_4936" name="KefrFDP" value="0.449149"/>
          <Constant key="Parameter_4935" name="KefrG6P" value="0.158746"/>
          <Constant key="Parameter_4934" name="KefrGL6P" value="0.150482"/>
          <Constant key="Parameter_4933" name="KefrR5P" value="9.33254"/>
          <Constant key="Parameter_4932" name="KefrRU5P" value="1.53591"/>
          <Constant key="Parameter_4931" name="KefrS7P" value="0.0785955"/>
          <Constant key="Parameter_4930" name="KefrX5P" value="0.677374"/>
          <Constant key="Parameter_4929" name="KeftATP" value="3.69117"/>
          <Constant key="Parameter_4928" name="KeftSUCCOA" value="8.26406"/>
          <Constant key="Parameter_4927" name="KirADP" value="0.517585"/>
          <Constant key="Parameter_4926" name="KirATP" value="96.0333"/>
          <Constant key="Parameter_4925" name="KirPEP" value="0.181056"/>
          <Constant key="Parameter_4924" name="KirPYR" value="15.1403"/>
          <Constant key="Parameter_4923" name="KirPyrATP" value="230.781"/>
          <Constant key="Parameter_4922" name="KitADP" value="0.224911"/>
          <Constant key="Parameter_4921" name="KitATP" value="0.039564"/>
          <Constant key="Parameter_4918" name="KitPEP" value="0.465672"/>
          <Constant key="Parameter_4919" name="KitPYR" value="0.2499"/>
          <Constant key="Parameter_4920" name="KitPyrATP" value="11.3691"/>
          <Constant key="Parameter_4917" name="KmrADPMg" value="0.326144"/>
          <Constant key="Parameter_4916" name="KmrPEP" value="5.56368e-07"/>
          <Constant key="Parameter_4915" name="KmtADPMg" value="0.054678"/>
          <Constant key="Parameter_4914" name="KmtPEP" value="0.11475"/>
          <Constant key="Parameter_4913" name="L0" value="50.4818"/>
          <Constant key="Parameter_4912" name="Vmax" value="0.74716"/>
          <Constant key="Parameter_4911" name="n" value="4"/>
        </ListOfConstants>
        <KineticLaw function="Function_55" unitType="Default" scalingCompartment="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm]">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_570">
              <SourceParameter reference="Metabolite_50"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_569">
              <SourceParameter reference="Metabolite_9"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_568">
              <SourceParameter reference="Metabolite_11"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_567">
              <SourceParameter reference="Metabolite_13"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_566">
              <SourceParameter reference="Parameter_4936"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_565">
              <SourceParameter reference="Parameter_4935"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_564">
              <SourceParameter reference="Parameter_4934"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_563">
              <SourceParameter reference="Parameter_4933"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_562">
              <SourceParameter reference="Parameter_4932"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_561">
              <SourceParameter reference="Parameter_4931"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_560">
              <SourceParameter reference="Parameter_4930"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_559">
              <SourceParameter reference="Parameter_4929"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_558">
              <SourceParameter reference="Parameter_4928"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_557">
              <SourceParameter reference="Parameter_4927"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_556">
              <SourceParameter reference="Parameter_4926"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_555">
              <SourceParameter reference="Parameter_4925"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_554">
              <SourceParameter reference="Parameter_4924"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_553">
              <SourceParameter reference="Parameter_4923"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_552">
              <SourceParameter reference="Parameter_4922"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_551">
              <SourceParameter reference="Parameter_4921"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_550">
              <SourceParameter reference="Parameter_4918"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_549">
              <SourceParameter reference="Parameter_4919"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_548">
              <SourceParameter reference="Parameter_4920"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_547">
              <SourceParameter reference="Parameter_4917"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_546">
              <SourceParameter reference="Parameter_4916"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_545">
              <SourceParameter reference="Parameter_4915"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_544">
              <SourceParameter reference="Parameter_4914"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_543">
              <SourceParameter reference="Parameter_4913"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_542">
              <SourceParameter reference="Metabolite_58"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_541">
              <SourceParameter reference="Metabolite_59"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_540">
              <SourceParameter reference="Metabolite_23"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_539">
              <SourceParameter reference="Metabolite_27"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_571">
              <SourceParameter reference="Metabolite_30"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_572">
              <SourceParameter reference="Metabolite_31"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_573">
              <SourceParameter reference="Metabolite_32"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_574">
              <SourceParameter reference="Metabolite_34"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_575">
              <SourceParameter reference="Parameter_4912"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_576">
              <SourceParameter reference="Metabolite_35"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_577">
              <SourceParameter reference="Parameter_4911"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_17" name="F6P_GAP_TAL" reversible="true" fast="false" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_17">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-27T11:53:43Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_12" stoichiometry="1"/>
          <Substrate metabolite="Metabolite_47" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_8" stoichiometry="1"/>
          <Product metabolite="Metabolite_46" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfConstants>
          <Constant key="Parameter_4866" name="Keq" value="0.11011"/>
          <Constant key="Parameter_4865" name="kcat" value="119.992"/>
        </ListOfConstants>
        <KineticLaw function="Function_61" unitType="Default" scalingCompartment="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm]">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_617">
              <SourceParameter reference="Metabolite_8"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_578">
              <SourceParameter reference="Metabolite_12"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_579">
              <SourceParameter reference="Parameter_4866"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_580">
              <SourceParameter reference="Parameter_4865"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_581">
              <SourceParameter reference="Metabolite_46"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_582">
              <SourceParameter reference="Metabolite_47"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_19" name="FBP" reversible="false" fast="false" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_19">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-01T17:07:07Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_9" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_8" stoichiometry="1"/>
          <Product metabolite="Metabolite_56" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfModifiers>
          <Modifier metabolite="Metabolite_51" stoichiometry="1"/>
          <Modifier metabolite="Metabolite_57" stoichiometry="1"/>
          <Modifier metabolite="Metabolite_60" stoichiometry="1"/>
        </ListOfModifiers>
        <ListOfConstants>
          <Constant key="Parameter_4862" name="KirAMP" value="0.00122122"/>
          <Constant key="Parameter_4861" name="KirAMPFDP" value="0.256256"/>
          <Constant key="Parameter_4860" name="KirF6P" value="1.12112"/>
          <Constant key="Parameter_4859" name="KirF6PMg" value="0.384615"/>
          <Constant key="Parameter_4858" name="KirFDP" value="1.35327"/>
          <Constant key="Parameter_4857" name="KirFDPMg" value="0.75924"/>
          <Constant key="Parameter_4856" name="KirFDPMgMg" value="0.356356"/>
          <Constant key="Parameter_4855" name="KirP" value="3.16316"/>
          <Constant key="Parameter_4854" name="KirPF6P" value="6.60538"/>
          <Constant key="Parameter_4853" name="KirPF6PMg" value="48.4484"/>
          <Constant key="Parameter_4852" name="KirPMg" value="0.856"/>
          <Constant key="Parameter_4851" name="KitAMP" value="0.000255"/>
          <Constant key="Parameter_4850" name="KitAMPFDP" value="690"/>
          <Constant key="Parameter_4849" name="KitF6P" value="0.304"/>
          <Constant key="Parameter_4848" name="KitF6PMg" value="315"/>
          <Constant key="Parameter_4847" name="KitFDP" value="0.043101"/>
          <Constant key="Parameter_4846" name="KitFDPMg" value="0.00642"/>
          <Constant key="Parameter_4845" name="KitFDPMgMg" value="100"/>
          <Constant key="Parameter_4844" name="KitP" value="0.642"/>
          <Constant key="Parameter_4843" name="KitPF6P" value="0.00689"/>
          <Constant key="Parameter_4842" name="KitPF6PMg" value="16.5"/>
          <Constant key="Parameter_4841" name="KitPMg" value="539"/>
          <Constant key="Parameter_4840" name="KmrFDP" value="0.0636141"/>
          <Constant key="Parameter_4839" name="KmrMg" value="0.039039"/>
          <Constant key="Parameter_4838" name="KmtFDP" value="1e-05"/>
          <Constant key="Parameter_4837" name="KmtMg" value="55.055"/>
          <Constant key="Parameter_4836" name="L0" value="0.000815"/>
          <Constant key="Parameter_4835" name="Vmax" value="0.215583"/>
          <Constant key="Parameter_4834" name="n" value="4"/>
          <Constant key="Parameter_4833" name="KdFDPMg" value="5.81"/>
        </ListOfConstants>
        <KineticLaw function="Function_46" unitType="Default" scalingCompartment="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm]">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_334">
              <SourceParameter reference="Metabolite_51"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_333">
              <SourceParameter reference="Metabolite_8"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_332">
              <SourceParameter reference="Metabolite_9"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_331">
              <SourceParameter reference="ModelValue_3"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_330">
              <SourceParameter reference="Parameter_4862"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_329">
              <SourceParameter reference="Parameter_4861"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_335">
              <SourceParameter reference="Parameter_4860"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_336">
              <SourceParameter reference="Parameter_4859"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_337">
              <SourceParameter reference="Parameter_4858"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_338">
              <SourceParameter reference="Parameter_4857"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_339">
              <SourceParameter reference="Parameter_4856"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_340">
              <SourceParameter reference="Parameter_4855"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_341">
              <SourceParameter reference="Parameter_4854"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_342">
              <SourceParameter reference="Parameter_4853"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_343">
              <SourceParameter reference="Parameter_4852"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_344">
              <SourceParameter reference="Parameter_4851"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_345">
              <SourceParameter reference="Parameter_4850"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_346">
              <SourceParameter reference="Parameter_4849"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_347">
              <SourceParameter reference="Parameter_4848"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_348">
              <SourceParameter reference="Parameter_4847"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_349">
              <SourceParameter reference="Parameter_4846"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_350">
              <SourceParameter reference="Parameter_4845"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_351">
              <SourceParameter reference="Parameter_4844"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_352">
              <SourceParameter reference="Parameter_4843"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_353">
              <SourceParameter reference="Parameter_4842"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_354">
              <SourceParameter reference="Parameter_4841"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_355">
              <SourceParameter reference="Parameter_4840"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_356">
              <SourceParameter reference="Parameter_4839"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_357">
              <SourceParameter reference="Parameter_4838"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_358">
              <SourceParameter reference="Parameter_4837"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_359">
              <SourceParameter reference="Parameter_4836"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_360">
              <SourceParameter reference="Metabolite_57"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_361">
              <SourceParameter reference="Metabolite_60"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_362">
              <SourceParameter reference="Metabolite_56"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_363">
              <SourceParameter reference="Parameter_4835"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_364">
              <SourceParameter reference="Parameter_4834"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_22" name="PPS" reversible="true" fast="false" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_22">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-01T13:57:27Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_52" stoichiometry="1"/>
          <Substrate metabolite="Metabolite_27" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_51" stoichiometry="1"/>
          <Product metabolite="Metabolite_23" stoichiometry="1"/>
          <Product metabolite="Metabolite_56" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfModifiers>
          <Modifier metabolite="Metabolite_50" stoichiometry="1"/>
          <Modifier metabolite="Metabolite_3" stoichiometry="1"/>
          <Modifier metabolite="Metabolite_22" stoichiometry="1"/>
          <Modifier metabolite="Metabolite_57" stoichiometry="1"/>
          <Modifier metabolite="Metabolite_58" stoichiometry="1"/>
          <Modifier metabolite="Metabolite_59" stoichiometry="1"/>
        </ListOfModifiers>
        <ListOfConstants>
          <Constant key="Parameter_4791" name="KdAMP" value="1480"/>
          <Constant key="Parameter_4790" name="KdATPMgPPS" value="0.0549"/>
          <Constant key="Parameter_4789" name="KdMg" value="36.9"/>
          <Constant key="Parameter_4788" name="KdP" value="346"/>
          <Constant key="Parameter_4787" name="KdPEP" value="95.7"/>
          <Constant key="Parameter_4786" name="KdPYR" value="2740"/>
          <Constant key="Parameter_4785" name="KefADP" value="0.0283"/>
          <Constant key="Parameter_4784" name="KefAKG" value="0.274"/>
          <Constant key="Parameter_4783" name="KefATP" value="0.000628"/>
          <Constant key="Parameter_4782" name="KefOAA" value="0.796"/>
          <Constant key="Parameter_4781" name="Keq" value="200000"/>
          <Constant key="Parameter_4780" name="KmAMP" value="0.000384"/>
          <Constant key="Parameter_4779" name="KmATPMg" value="0.0549"/>
          <Constant key="Parameter_4778" name="KmP" value="84.4"/>
          <Constant key="Parameter_4777" name="KmPEP" value="20.7"/>
          <Constant key="Parameter_4776" name="KmPYR" value="0.229"/>
          <Constant key="Parameter_4775" name="Vmax" value="0.0163772"/>
          <Constant key="Parameter_4774" name="W" value="10"/>
          <Constant key="Parameter_4773" name="alpha" value="38900"/>
          <Constant key="Parameter_4772" name="KdADPMg" value="1.27771"/>
          <Constant key="Parameter_4771" name="KdATPMg" value="0.0847634"/>
        </ListOfConstants>
        <KineticLaw function="Function_50" unitType="Default" scalingCompartment="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm]">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_420">
              <SourceParameter reference="Metabolite_50"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_421">
              <SourceParameter reference="Metabolite_3"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_422">
              <SourceParameter reference="Metabolite_51"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_423">
              <SourceParameter reference="Metabolite_52"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_424">
              <SourceParameter reference="ModelValue_1"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_425">
              <SourceParameter reference="Parameter_4791"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_426">
              <SourceParameter reference="ModelValue_2"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_427">
              <SourceParameter reference="Parameter_4790"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_428">
              <SourceParameter reference="Parameter_4789"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_429">
              <SourceParameter reference="Parameter_4788"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_430">
              <SourceParameter reference="Parameter_4787"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_443">
              <SourceParameter reference="Parameter_4786"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_419">
              <SourceParameter reference="Parameter_4785"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_418">
              <SourceParameter reference="Parameter_4784"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_417">
              <SourceParameter reference="Parameter_4783"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_416">
              <SourceParameter reference="Parameter_4782"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_415">
              <SourceParameter reference="Parameter_4781"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_414">
              <SourceParameter reference="Parameter_4780"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_413">
              <SourceParameter reference="Parameter_4779"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_412">
              <SourceParameter reference="Parameter_4778"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_411">
              <SourceParameter reference="Parameter_4777"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_457">
              <SourceParameter reference="Parameter_4776"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_458">
              <SourceParameter reference="Metabolite_57"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_459">
              <SourceParameter reference="Metabolite_58"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_460">
              <SourceParameter reference="Metabolite_59"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_461">
              <SourceParameter reference="Metabolite_22"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_462">
              <SourceParameter reference="Metabolite_56"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_463">
              <SourceParameter reference="Metabolite_23"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_464">
              <SourceParameter reference="Metabolite_27"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_465">
              <SourceParameter reference="Parameter_4775"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_466">
              <SourceParameter reference="Parameter_4774"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_467">
              <SourceParameter reference="Parameter_4773"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_50" name="PTS_0" reversible="true" fast="false" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_50">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-30T12:14:31Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_36" stoichiometry="1"/>
          <Substrate metabolite="Metabolite_23" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_37" stoichiometry="1"/>
          <Product metabolite="Metabolite_27" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfConstants>
          <Constant key="Parameter_4608" name="KmPEP" value="0.6"/>
          <Constant key="Parameter_4607" name="KmPYR" value="1"/>
          <Constant key="Parameter_4606" name="kF" value="12000"/>
          <Constant key="Parameter_4605" name="kR" value="8000"/>
        </ListOfConstants>
        <KineticLaw function="Function_70" unitType="Default" scalingCompartment="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm]">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_701">
              <SourceParameter reference="Parameter_4608"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_700">
              <SourceParameter reference="Parameter_4607"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_699">
              <SourceParameter reference="Metabolite_23"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_698">
              <SourceParameter reference="Metabolite_27"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_697">
              <SourceParameter reference="Metabolite_36"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_696">
              <SourceParameter reference="Metabolite_37"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_695">
              <SourceParameter reference="Parameter_4606"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_694">
              <SourceParameter reference="Parameter_4605"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_51" name="PTS_1" reversible="true" fast="false" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_51">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2015-01-03T00:49:35Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_42" stoichiometry="1"/>
          <Substrate metabolite="Metabolite_37" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_43" stoichiometry="1"/>
          <Product metabolite="Metabolite_36" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfConstants>
          <Constant key="Parameter_4604" name="k1" value="200000"/>
          <Constant key="Parameter_4603" name="k2" value="8000"/>
        </ListOfConstants>
        <KineticLaw function="Function_14" unitType="Default" scalingCompartment="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm]">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_69">
              <SourceParameter reference="Parameter_4604"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_68">
              <SourceParameter reference="Metabolite_42"/>
              <SourceParameter reference="Metabolite_37"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_78">
              <SourceParameter reference="Parameter_4603"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_79">
              <SourceParameter reference="Metabolite_43"/>
              <SourceParameter reference="Metabolite_36"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_52" name="PTS_2" reversible="true" fast="false" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_52">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-15T10:24:08Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_38" stoichiometry="1"/>
          <Substrate metabolite="Metabolite_43" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_39" stoichiometry="1"/>
          <Product metabolite="Metabolite_42" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfConstants>
          <Constant key="Parameter_4602" name="k1" value="61000"/>
          <Constant key="Parameter_4601" name="k2" value="47000"/>
        </ListOfConstants>
        <KineticLaw function="Function_14" unitType="Default" scalingCompartment="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm]">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_69">
              <SourceParameter reference="Parameter_4602"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_68">
              <SourceParameter reference="Metabolite_38"/>
              <SourceParameter reference="Metabolite_43"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_78">
              <SourceParameter reference="Parameter_4601"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_79">
              <SourceParameter reference="Metabolite_39"/>
              <SourceParameter reference="Metabolite_42"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_53" name="PTS_3" reversible="true" fast="false" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_53">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-30T12:14:41Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_40" stoichiometry="1"/>
          <Substrate metabolite="Metabolite_39" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_41" stoichiometry="1"/>
          <Product metabolite="Metabolite_38" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfConstants>
          <Constant key="Parameter_4600" name="k1" value="11000"/>
          <Constant key="Parameter_4599" name="k2" value="4000"/>
        </ListOfConstants>
        <KineticLaw function="Function_14" unitType="Default" scalingCompartment="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm]">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_69">
              <SourceParameter reference="Parameter_4600"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_68">
              <SourceParameter reference="Metabolite_40"/>
              <SourceParameter reference="Metabolite_39"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_78">
              <SourceParameter reference="Parameter_4599"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_79">
              <SourceParameter reference="Metabolite_41"/>
              <SourceParameter reference="Metabolite_38"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_54" name="PTS_4" reversible="true" fast="false" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_54">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-30T12:14:34Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_74" stoichiometry="1"/>
          <Substrate metabolite="Metabolite_41" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_11" stoichiometry="1"/>
          <Product metabolite="Metabolite_40" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfConstants>
          <Constant key="Parameter_4598" name="KmG6P" value="2125.91"/>
          <Constant key="Parameter_4597" name="KmGLC" value="0.02"/>
          <Constant key="Parameter_4596" name="kF" value="4000"/>
          <Constant key="Parameter_4595" name="kR" value="1e-05"/>
        </ListOfConstants>
        <KineticLaw function="Function_92" unitType="Default">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_897">
              <SourceParameter reference="Metabolite_11"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_896">
              <SourceParameter reference="Metabolite_74"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_895">
              <SourceParameter reference="Parameter_4598"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_894">
              <SourceParameter reference="Parameter_4597"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_893">
              <SourceParameter reference="Compartment_0"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_892">
              <SourceParameter reference="Metabolite_40"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_891">
              <SourceParameter reference="Metabolite_41"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_898">
              <SourceParameter reference="Parameter_4596"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_899">
              <SourceParameter reference="Parameter_4595"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_55" name="_GLC_FEED" reversible="false" fast="false" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_55">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-09-30T23:49:02Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

        </MiriamAnnotation>
        <ListOfProducts>
          <Product metabolite="Metabolite_70" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfConstants>
          <Constant key="Parameter_4594" name="FEED" value="0.23"/>
        </ListOfConstants>
        <KineticLaw function="Function_73" unitType="Default" scalingCompartment="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[extracellular]">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_719">
              <SourceParameter reference="ModelValue_0"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_718">
              <SourceParameter reference="Compartment_0"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_708">
              <SourceParameter reference="Compartment_1"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_60" name="ATP_MAINTENANCE" reversible="true" fast="false" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_60">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-21T14:43:37Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_52" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_50" stoichiometry="1"/>
          <Product metabolite="Metabolite_56" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfConstants>
          <Constant key="Parameter_4572" name="Vmax" value="1.30166"/>
          <Constant key="Parameter_4571" name="Keq" value="3.63369"/>
        </ListOfConstants>
        <KineticLaw function="Function_97" unitType="Default" scalingCompartment="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm]">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_976">
              <SourceParameter reference="Metabolite_50"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_975">
              <SourceParameter reference="Metabolite_52"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_974">
              <SourceParameter reference="Parameter_4571"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_973">
              <SourceParameter reference="Metabolite_56"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_972">
              <SourceParameter reference="Parameter_4572"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_61" name="XCH_GLC" reversible="true" fast="false" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_61">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2013-10-23T16:07:48Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_70" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_74" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfConstants>
          <Constant key="Parameter_4570" name="Vmax" value="100"/>
          <Constant key="Parameter_4569" name="Km" value="10"/>
        </ListOfConstants>
        <KineticLaw function="Function_98" unitType="Default">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_967">
              <SourceParameter reference="Metabolite_74"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_968">
              <SourceParameter reference="Metabolite_70"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_969">
              <SourceParameter reference="Parameter_4569"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_970">
              <SourceParameter reference="Parameter_4570"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
      <Reaction key="Reaction_67" name="GL6P_HYDROLYSIS" reversible="true" fast="false" addNoise="false">
        <MiriamAnnotation>
<rdf:RDF
   xmlns:dcterms="http://purl.org/dc/terms/"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
  <rdf:Description rdf:about="#Reaction_67">
    <dcterms:created>
      <rdf:Description>
        <dcterms:W3CDTF>2015-01-22T11:03:40Z</dcterms:W3CDTF>
      </rdf:Description>
    </dcterms:created>
  </rdf:Description>
</rdf:RDF>

        </MiriamAnnotation>
        <ListOfSubstrates>
          <Substrate metabolite="Metabolite_13" stoichiometry="1"/>
        </ListOfSubstrates>
        <ListOfProducts>
          <Product metabolite="Metabolite_26" stoichiometry="1"/>
        </ListOfProducts>
        <ListOfConstants>
          <Constant key="Parameter_4557" name="KGl6Phydrol" value="0.000167"/>
          <Constant key="Parameter_4556" name="KeqGl6Phydrol" value="42.8"/>
        </ListOfConstants>
        <KineticLaw function="Function_103" unitType="Default" scalingCompartment="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm]">
          <ListOfCallParameters>
            <CallParameter functionParameter="FunctionParameter_980">
              <SourceParameter reference="Metabolite_13"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_979">
              <SourceParameter reference="Parameter_4557"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_978">
              <SourceParameter reference="Parameter_4556"/>
            </CallParameter>
            <CallParameter functionParameter="FunctionParameter_977">
              <SourceParameter reference="Metabolite_26"/>
            </CallParameter>
          </ListOfCallParameters>
        </KineticLaw>
      </Reaction>
    </ListOfReactions>
    <ListOfModelParameterSets activeSet="ModelParameterSet_1">
      <ModelParameterSet key="ModelParameterSet_1" name="Initial State">
        <ModelParameterGroup cn="String=Initial Time" type="Group">
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism" value="0" type="Model" simulationType="time"/>
        </ModelParameterGroup>
        <ModelParameterGroup cn="String=Initial Compartment Sizes" type="Group">
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm]" value="1" type="Compartment" simulationType="fixed"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[extracellular]" value="100" type="Compartment" simulationType="fixed"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_periplasm]" value="0.25" type="Compartment" simulationType="fixed"/>
        </ModelParameterGroup>
        <ModelParameterGroup cn="String=Initial Species Values" type="Group">
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[AKG]" value="3.6004593044112428e+20" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[BPG]" value="3.9391164070065283e+19" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[DAP]" value="2.6322420431957089e+20" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[F6P]" value="1.5763946066185716e+20" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[FDP]" value="1.6970884333545056e+20" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[G6P]" value="5.1858431598197925e+20" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[GAP]" value="7.0569319205180695e+19" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[GL6P]" value="1.9642127972845123e+18" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[NAD]" value="8.5005146998199353e+20" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[NADH]" value="9.5424644986276282e+19" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[OAA]" value="7.6986801648223224e+19" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[PEP]" value="6.0043053493910241e+20" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[PGA2]" value="2.2781584503402476e+20" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[PGA3]" value="4.193061098189998e+20" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[PGN]" value="7.9251241914975207e+19" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[PYR]" value="1.4265932662081769e+20" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[R5P]" value="6.4341569172230717e+19" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[RU5P]" value="2.0585278619139365e+20" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[S7P]" value="8.5505452688743842e+19" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[SUCCOA]" value="2.474363866819711e+19" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[X5P]" value="3.0473116183618054e+20" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[ei]" value="2.0114735423457258e+17" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[eiP]" value="3.8428260841912832e+18" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[eiia]" value="8.5525575608716513e+18" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[eiiaP]" value="2.8940514631167993e+20" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[eiicb]" value="28419065176204708" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[eiicbP]" value="2.0633142966785149e+17" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[hpr]" value="1.1515043650212757e+17" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[hprP]" value="3.1621128357967836e+18" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[tal]" value="1.6747235514063363e+18" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[talC3]" value="3.4458121569610027e+19" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[ADP]" value="3.6031354180883579e+20" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[AMP]" value="1.1216394702819669e+20" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[ATP]" value="1.5490159511281688e+21" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[P]" value="5.8799869312313674e+21" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[MG]" value="6.0221408570000002e+20" type="Species" simulationType="fixed"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[MgADP]" value="1.581911401402443e+20" type="Species" simulationType="assignment"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[MgATP]" value="1.42797586195125e+21" type="Species" simulationType="assignment"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_cytoplasm],Vector=Metabolites[MgFDP]" value="2.4920534997863522e+19" type="Species" simulationType="assignment"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[extracellular],Vector=Metabolites[GLCx]" value="3.8140444595268039e+20" type="Species" simulationType="reactions"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Compartments[cell_periplasm],Vector=Metabolites[GLCp]" value="6.0723801560420096e+17" type="Species" simulationType="reactions"/>
        </ModelParameterGroup>
        <ModelParameterGroup cn="String=Initial Global Quantities" type="Group">
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Values[FEED]" value="0.23000000000000001" type="ModelValue" simulationType="fixed"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Values[KdADPMg]" value="1.2777099999999999" type="ModelValue" simulationType="fixed"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Values[KdATPMg]" value="0.084763400000000003" type="ModelValue" simulationType="fixed"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Values[KdFDPMg]" value="5.8099999999999996" type="ModelValue" simulationType="fixed"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Values[KmICIT_ACN]" value="9.3135200000000005" type="ModelValue" simulationType="fixed"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Values[KmCIT_ACN]" value="0.062888200000000005" type="ModelValue" simulationType="fixed"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Values[KmACO_ACN]" value="0.02001" type="ModelValue" simulationType="fixed"/>
          <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Values[KeqNDH]" value="27.619299999999999" type="ModelValue" simulationType="fixed"/>
        </ModelParameterGroup>
        <ModelParameterGroup cn="String=Kinetic Parameters" type="Group">
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PGI]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PGI],ParameterGroup=Parameters,Parameter=Keq" value="0.35999999999999999" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PGI],ParameterGroup=Parameters,Parameter=KmF6P" value="0.14699999999999999" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PGI],ParameterGroup=Parameters,Parameter=KmG6P" value="0.28000000000000003" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PGI],ParameterGroup=Parameters,Parameter=KmPEP" value="1.9990000000000001" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PGI],ParameterGroup=Parameters,Parameter=Vmax" value="2.32456" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PGI],ParameterGroup=Parameters,Parameter=KmPGN" value="0.51595800000000003" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PFK]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PFK],ParameterGroup=Parameters,Parameter=KefrADP" value="0.073526400000000006" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PFK],ParameterGroup=Parameters,Parameter=KefrPEP" value="19.98" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PFK],ParameterGroup=Parameters,Parameter=KeftADP" value="9.0090000000000003" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PFK],ParameterGroup=Parameters,Parameter=KeftPEP" value="0.26025999999999999" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PFK],ParameterGroup=Parameters,Parameter=Keq" value="1998" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PFK],ParameterGroup=Parameters,Parameter=KirADP" value="54.945" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PFK],ParameterGroup=Parameters,Parameter=KirATP" value="2.4975000000000001e-05" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PFK],ParameterGroup=Parameters,Parameter=KirF6P" value="1.84615" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PFK],ParameterGroup=Parameters,Parameter=KirFDP" value="0.045954000000000002" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PFK],ParameterGroup=Parameters,Parameter=KitADP" value="80.079999999999998" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PFK],ParameterGroup=Parameters,Parameter=KitATP" value="0.014014" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PFK],ParameterGroup=Parameters,Parameter=KitF6P" value="0.0085685599999999994" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PFK],ParameterGroup=Parameters,Parameter=KitFDP" value="50.5505" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PFK],ParameterGroup=Parameters,Parameter=KmrADP" value="0.69000899999999998" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PFK],ParameterGroup=Parameters,Parameter=KmrATPMg" value="8.1218700000000002e-05" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PFK],ParameterGroup=Parameters,Parameter=KmrF6P" value="2.0520500000000001e-05" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PFK],ParameterGroup=Parameters,Parameter=KmrFDP" value="10.01" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PFK],ParameterGroup=Parameters,Parameter=KmtADP" value="2.0019999999999998" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PFK],ParameterGroup=Parameters,Parameter=KmtATPMg" value="3.34334" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PFK],ParameterGroup=Parameters,Parameter=KmtF6P" value="32.966999999999999" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PFK],ParameterGroup=Parameters,Parameter=KmtFDP" value="9.9900000000000002" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PFK],ParameterGroup=Parameters,Parameter=L0" value="14.085100000000001" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PFK],ParameterGroup=Parameters,Parameter=Vmax" value="0.185253" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PFK],ParameterGroup=Parameters,Parameter=Wr" value="0.023704099999999999" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PFK],ParameterGroup=Parameters,Parameter=Wt" value="0.146735" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PFK],ParameterGroup=Parameters,Parameter=n" value="4" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBA]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBA],ParameterGroup=Parameters,Parameter=Keq" value="0.18981000000000001" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBA],ParameterGroup=Parameters,Parameter=KmDAP" value="0.13000999999999999" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBA],ParameterGroup=Parameters,Parameter=KmFDP" value="0.12012" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBA],ParameterGroup=Parameters,Parameter=KmGAP" value="0.13000999999999999" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBA],ParameterGroup=Parameters,Parameter=KmPEP" value="0.5" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBA],ParameterGroup=Parameters,Parameter=Vmax" value="21.697800000000001" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[TPI]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[TPI],ParameterGroup=Parameters,Parameter=Keq" value="0.27020300000000003" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[TPI],ParameterGroup=Parameters,Parameter=KmDAP" value="0.01" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[TPI],ParameterGroup=Parameters,Parameter=KmGAP" value="1.8930100000000001" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[TPI],ParameterGroup=Parameters,Parameter=Vmax" value="24.1843" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GDH]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GDH],ParameterGroup=Parameters,Parameter=Keq" value="20" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GDH],ParameterGroup=Parameters,Parameter=KmBPG" value="0.20000000000000001" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GDH],ParameterGroup=Parameters,Parameter=KmGAP" value="2.4726499999999998" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GDH],ParameterGroup=Parameters,Parameter=KmNAD" value="0.0110454" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GDH],ParameterGroup=Parameters,Parameter=KmNADH" value="3.6979700000000002" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GDH],ParameterGroup=Parameters,Parameter=KmP" value="0.017000000000000001" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GDH],ParameterGroup=Parameters,Parameter=Vmax" value="8.6657299999999999" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PGK]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PGK],ParameterGroup=Parameters,Parameter=Keq" value="99.992500000000007" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PGK],ParameterGroup=Parameters,Parameter=KmADPMg" value="0.085416000000000006" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PGK],ParameterGroup=Parameters,Parameter=KmATPMg" value="3.4773700000000001" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PGK],ParameterGroup=Parameters,Parameter=KmBPG" value="0.0113296" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PGK],ParameterGroup=Parameters,Parameter=KmPGA3" value="2.45722" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PGK],ParameterGroup=Parameters,Parameter=Vmax" value="16.108899999999998" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GPM]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GPM],ParameterGroup=Parameters,Parameter=Keq" value="0.56581800000000004" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GPM],ParameterGroup=Parameters,Parameter=KmPGA2" value="1.9153" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GPM],ParameterGroup=Parameters,Parameter=KmPGA3" value="0.115" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GPM],ParameterGroup=Parameters,Parameter=Vmax" value="10.993399999999999" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ENO]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ENO],ParameterGroup=Parameters,Parameter=Keq" value="3" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ENO],ParameterGroup=Parameters,Parameter=KmPEP" value="0.10000000000000001" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ENO],ParameterGroup=Parameters,Parameter=KmPGA2" value="0.10000000000000001" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ENO],ParameterGroup=Parameters,Parameter=Vmax" value="11.7189" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PYK]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PYK],ParameterGroup=Parameters,Parameter=KefrFDP" value="0.44914900000000002" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PYK],ParameterGroup=Parameters,Parameter=KefrG6P" value="0.158746" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PYK],ParameterGroup=Parameters,Parameter=KefrGL6P" value="0.150482" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PYK],ParameterGroup=Parameters,Parameter=KefrR5P" value="9.3325399999999998" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PYK],ParameterGroup=Parameters,Parameter=KefrRU5P" value="1.5359100000000001" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PYK],ParameterGroup=Parameters,Parameter=KefrS7P" value="0.078595499999999999" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PYK],ParameterGroup=Parameters,Parameter=KefrX5P" value="0.67737400000000003" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PYK],ParameterGroup=Parameters,Parameter=KeftATP" value="3.6911700000000001" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PYK],ParameterGroup=Parameters,Parameter=KeftSUCCOA" value="8.2640600000000006" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PYK],ParameterGroup=Parameters,Parameter=KirADP" value="0.51758499999999996" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PYK],ParameterGroup=Parameters,Parameter=KirATP" value="96.033299999999997" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PYK],ParameterGroup=Parameters,Parameter=KirPEP" value="0.18105599999999999" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PYK],ParameterGroup=Parameters,Parameter=KirPYR" value="15.1403" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PYK],ParameterGroup=Parameters,Parameter=KirPyrATP" value="230.78100000000001" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PYK],ParameterGroup=Parameters,Parameter=KitADP" value="0.224911" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PYK],ParameterGroup=Parameters,Parameter=KitATP" value="0.039564000000000002" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PYK],ParameterGroup=Parameters,Parameter=KitPEP" value="0.46567199999999997" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PYK],ParameterGroup=Parameters,Parameter=KitPYR" value="0.24990000000000001" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PYK],ParameterGroup=Parameters,Parameter=KitPyrATP" value="11.3691" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PYK],ParameterGroup=Parameters,Parameter=KmrADPMg" value="0.32614399999999999" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PYK],ParameterGroup=Parameters,Parameter=KmrPEP" value="5.5636800000000005e-07" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PYK],ParameterGroup=Parameters,Parameter=KmtADPMg" value="0.054677999999999997" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PYK],ParameterGroup=Parameters,Parameter=KmtPEP" value="0.11475" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PYK],ParameterGroup=Parameters,Parameter=L0" value="50.4818" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PYK],ParameterGroup=Parameters,Parameter=Vmax" value="0.74716000000000005" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PYK],ParameterGroup=Parameters,Parameter=n" value="4" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[F6P_GAP_TAL]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[F6P_GAP_TAL],ParameterGroup=Parameters,Parameter=Keq" value="0.11011" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[F6P_GAP_TAL],ParameterGroup=Parameters,Parameter=kcat" value="119.992" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBP]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBP],ParameterGroup=Parameters,Parameter=KirAMP" value="0.0012212200000000001" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBP],ParameterGroup=Parameters,Parameter=KirAMPFDP" value="0.25625599999999998" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBP],ParameterGroup=Parameters,Parameter=KirF6P" value="1.1211199999999999" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBP],ParameterGroup=Parameters,Parameter=KirF6PMg" value="0.38461499999999998" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBP],ParameterGroup=Parameters,Parameter=KirFDP" value="1.35327" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBP],ParameterGroup=Parameters,Parameter=KirFDPMg" value="0.75924000000000003" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBP],ParameterGroup=Parameters,Parameter=KirFDPMgMg" value="0.35635600000000001" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBP],ParameterGroup=Parameters,Parameter=KirP" value="3.16316" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBP],ParameterGroup=Parameters,Parameter=KirPF6P" value="6.6053800000000003" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBP],ParameterGroup=Parameters,Parameter=KirPF6PMg" value="48.448399999999999" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBP],ParameterGroup=Parameters,Parameter=KirPMg" value="0.85599999999999998" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBP],ParameterGroup=Parameters,Parameter=KitAMP" value="0.00025500000000000002" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBP],ParameterGroup=Parameters,Parameter=KitAMPFDP" value="690" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBP],ParameterGroup=Parameters,Parameter=KitF6P" value="0.30399999999999999" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBP],ParameterGroup=Parameters,Parameter=KitF6PMg" value="315" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBP],ParameterGroup=Parameters,Parameter=KitFDP" value="0.043101" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBP],ParameterGroup=Parameters,Parameter=KitFDPMg" value="0.0064200000000000004" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBP],ParameterGroup=Parameters,Parameter=KitFDPMgMg" value="100" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBP],ParameterGroup=Parameters,Parameter=KitP" value="0.64200000000000002" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBP],ParameterGroup=Parameters,Parameter=KitPF6P" value="0.0068900000000000003" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBP],ParameterGroup=Parameters,Parameter=KitPF6PMg" value="16.5" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBP],ParameterGroup=Parameters,Parameter=KitPMg" value="539" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBP],ParameterGroup=Parameters,Parameter=KmrFDP" value="0.063614100000000007" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBP],ParameterGroup=Parameters,Parameter=KmrMg" value="0.039038999999999997" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBP],ParameterGroup=Parameters,Parameter=KmtFDP" value="1.0000000000000001e-05" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBP],ParameterGroup=Parameters,Parameter=KmtMg" value="55.055" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBP],ParameterGroup=Parameters,Parameter=L0" value="0.00081499999999999997" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBP],ParameterGroup=Parameters,Parameter=Vmax" value="0.215583" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBP],ParameterGroup=Parameters,Parameter=n" value="4" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[FBP],ParameterGroup=Parameters,Parameter=KdFDPMg" value="5.8099999999999996" type="ReactionParameter" simulationType="assignment">
              <InitialExpression>
                &lt;CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Values[KdFDPMg],Reference=InitialValue>
              </InitialExpression>
            </ModelParameter>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPS]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPS],ParameterGroup=Parameters,Parameter=KdAMP" value="1480" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPS],ParameterGroup=Parameters,Parameter=KdATPMgPPS" value="0.054899999999999997" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPS],ParameterGroup=Parameters,Parameter=KdMg" value="36.899999999999999" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPS],ParameterGroup=Parameters,Parameter=KdP" value="346" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPS],ParameterGroup=Parameters,Parameter=KdPEP" value="95.700000000000003" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPS],ParameterGroup=Parameters,Parameter=KdPYR" value="2740" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPS],ParameterGroup=Parameters,Parameter=KefADP" value="0.028299999999999999" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPS],ParameterGroup=Parameters,Parameter=KefAKG" value="0.27400000000000002" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPS],ParameterGroup=Parameters,Parameter=KefATP" value="0.00062799999999999998" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPS],ParameterGroup=Parameters,Parameter=KefOAA" value="0.79600000000000004" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPS],ParameterGroup=Parameters,Parameter=Keq" value="200000" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPS],ParameterGroup=Parameters,Parameter=KmAMP" value="0.00038400000000000001" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPS],ParameterGroup=Parameters,Parameter=KmATPMg" value="0.054899999999999997" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPS],ParameterGroup=Parameters,Parameter=KmP" value="84.400000000000006" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPS],ParameterGroup=Parameters,Parameter=KmPEP" value="20.699999999999999" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPS],ParameterGroup=Parameters,Parameter=KmPYR" value="0.22900000000000001" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPS],ParameterGroup=Parameters,Parameter=Vmax" value="0.016377200000000001" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPS],ParameterGroup=Parameters,Parameter=W" value="10" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPS],ParameterGroup=Parameters,Parameter=alpha" value="38900" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPS],ParameterGroup=Parameters,Parameter=KdADPMg" value="1.2777099999999999" type="ReactionParameter" simulationType="assignment">
              <InitialExpression>
                &lt;CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Values[KdADPMg],Reference=InitialValue>
              </InitialExpression>
            </ModelParameter>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PPS],ParameterGroup=Parameters,Parameter=KdATPMg" value="0.084763400000000003" type="ReactionParameter" simulationType="assignment">
              <InitialExpression>
                &lt;CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Values[KdATPMg],Reference=InitialValue>
              </InitialExpression>
            </ModelParameter>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PTS_0]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PTS_0],ParameterGroup=Parameters,Parameter=KmPEP" value="0.59999999999999998" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PTS_0],ParameterGroup=Parameters,Parameter=KmPYR" value="1" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PTS_0],ParameterGroup=Parameters,Parameter=kF" value="12000" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PTS_0],ParameterGroup=Parameters,Parameter=kR" value="8000" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PTS_1]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PTS_1],ParameterGroup=Parameters,Parameter=k1" value="200000" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PTS_1],ParameterGroup=Parameters,Parameter=k2" value="8000" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PTS_2]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PTS_2],ParameterGroup=Parameters,Parameter=k1" value="61000" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PTS_2],ParameterGroup=Parameters,Parameter=k2" value="47000" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PTS_3]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PTS_3],ParameterGroup=Parameters,Parameter=k1" value="11000" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PTS_3],ParameterGroup=Parameters,Parameter=k2" value="4000" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PTS_4]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PTS_4],ParameterGroup=Parameters,Parameter=KmG6P" value="2125.9099999999999" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PTS_4],ParameterGroup=Parameters,Parameter=KmGLC" value="0.02" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PTS_4],ParameterGroup=Parameters,Parameter=kF" value="4000" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[PTS_4],ParameterGroup=Parameters,Parameter=kR" value="1.0000000000000001e-05" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[_GLC_FEED]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[_GLC_FEED],ParameterGroup=Parameters,Parameter=FEED" value="0.23000000000000001" type="ReactionParameter" simulationType="assignment">
              <InitialExpression>
                &lt;CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Values[FEED],Reference=InitialValue>
              </InitialExpression>
            </ModelParameter>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ATP_MAINTENANCE]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ATP_MAINTENANCE],ParameterGroup=Parameters,Parameter=Vmax" value="1.30166" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[ATP_MAINTENANCE],ParameterGroup=Parameters,Parameter=Keq" value="3.6336900000000001" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[XCH_GLC]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[XCH_GLC],ParameterGroup=Parameters,Parameter=Vmax" value="100" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[XCH_GLC],ParameterGroup=Parameters,Parameter=Km" value="10" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
          <ModelParameterGroup cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GL6P_HYDROLYSIS]" type="Reaction">
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GL6P_HYDROLYSIS],ParameterGroup=Parameters,Parameter=KGl6Phydrol" value="0.00016699999999999999" type="ReactionParameter" simulationType="fixed"/>
            <ModelParameter cn="CN=Root,Model=Millard2016 - E. coli central carbon and energy metabolism,Vector=Reactions[GL6P_HYDROLYSIS],ParameterGroup=Parameters,Parameter=KeqGl6Phydrol" value="42.799999999999997" type="ReactionParameter" simulationType="fixed"/>
          </ModelParameterGroup>
        </ModelParameterGroup>
      </ModelParameterSet>
    </ListOfModelParameterSets>
    <StateTemplate>
      <StateTemplateVariable objectReference="Model_1"/>
      <StateTemplateVariable objectReference="Metabolite_52"/>
      <StateTemplateVariable objectReference="Metabolite_12"/>
      <StateTemplateVariable objectReference="Metabolite_8"/>
      <StateTemplateVariable objectReference="Metabolite_23"/>
      <StateTemplateVariable objectReference="Metabolite_56"/>
      <StateTemplateVariable objectReference="Metabolite_6"/>
      <StateTemplateVariable objectReference="Metabolite_38"/>
      <StateTemplateVariable objectReference="Metabolite_70"/>
      <StateTemplateVariable objectReference="Metabolite_24"/>
      <StateTemplateVariable objectReference="Metabolite_36"/>
      <StateTemplateVariable objectReference="Metabolite_11"/>
      <StateTemplateVariable objectReference="Metabolite_25"/>
      <StateTemplateVariable objectReference="Metabolite_26"/>
      <StateTemplateVariable objectReference="Metabolite_42"/>
      <StateTemplateVariable objectReference="Metabolite_74"/>
      <StateTemplateVariable objectReference="Metabolite_4"/>
      <StateTemplateVariable objectReference="Metabolite_40"/>
      <StateTemplateVariable objectReference="Metabolite_50"/>
      <StateTemplateVariable objectReference="Metabolite_46"/>
      <StateTemplateVariable objectReference="Metabolite_27"/>
      <StateTemplateVariable objectReference="Metabolite_9"/>
      <StateTemplateVariable objectReference="Metabolite_43"/>
      <StateTemplateVariable objectReference="Metabolite_13"/>
      <StateTemplateVariable objectReference="Metabolite_47"/>
      <StateTemplateVariable objectReference="Metabolite_39"/>
      <StateTemplateVariable objectReference="Metabolite_51"/>
      <StateTemplateVariable objectReference="Metabolite_37"/>
      <StateTemplateVariable objectReference="Metabolite_19"/>
      <StateTemplateVariable objectReference="Metabolite_18"/>
      <StateTemplateVariable objectReference="Metabolite_41"/>
      <StateTemplateVariable objectReference="Metabolite_58"/>
      <StateTemplateVariable objectReference="Metabolite_59"/>
      <StateTemplateVariable objectReference="Metabolite_60"/>
      <StateTemplateVariable objectReference="Metabolite_57"/>
      <StateTemplateVariable objectReference="Metabolite_3"/>
      <StateTemplateVariable objectReference="Metabolite_22"/>
      <StateTemplateVariable objectReference="Metabolite_30"/>
      <StateTemplateVariable objectReference="Metabolite_31"/>
      <StateTemplateVariable objectReference="Metabolite_32"/>
      <StateTemplateVariable objectReference="Metabolite_34"/>
      <StateTemplateVariable objectReference="Metabolite_35"/>
      <StateTemplateVariable objectReference="Compartment_0"/>
      <StateTemplateVariable objectReference="Compartment_1"/>
      <StateTemplateVariable objectReference="Compartment_2"/>
      <StateTemplateVariable objectReference="ModelValue_0"/>
      <StateTemplateVariable objectReference="ModelValue_1"/>
      <StateTemplateVariable objectReference="ModelValue_2"/>
      <StateTemplateVariable objectReference="ModelValue_3"/>
      <StateTemplateVariable objectReference="ModelValue_4"/>
      <StateTemplateVariable objectReference="ModelValue_5"/>
      <StateTemplateVariable objectReference="ModelValue_6"/>
      <StateTemplateVariable objectReference="ModelValue_7"/>
    </StateTemplate>
    <InitialState type="initialState">
      0 1.5490159511281688e+21 7.0569319205180695e+19 1.5763946066185716e+20 6.0043053493910241e+20 5.8799869312313674e+21 2.6322420431957089e+20 8.5525575608716513e+18 3.8140444595268039e+20 2.2781584503402476e+20 2.0114735423457258e+17 5.1858431598197925e+20 4.193061098189998e+20 7.9251241914975207e+19 1.1515043650212757e+17 6.0723801560420096e+17 3.9391164070065283e+19 28419065176204708 3.6031354180883579e+20 1.6747235514063363e+18 1.4265932662081769e+20 1.6970884333545056e+20 3.1621128357967836e+18 1.9642127972845123e+18 3.4458121569610027e+19 2.8940514631167993e+20 1.1216394702819669e+20 3.8428260841912832e+18 9.5424644986276282e+19 8.5005146998199353e+20 2.0633142966785149e+17 1.581911401402443e+20 1.42797586195125e+21 2.4920534997863522e+19 6.0221408570000002e+20 3.6004593044112428e+20 7.6986801648223224e+19 6.4341569172230717e+19 2.0585278619139365e+20 8.5505452688743842e+19 2.474363866819711e+19 3.0473116183618054e+20 1 100 0.25 0.23000000000000001 1.2777099999999999 0.084763400000000003 5.8099999999999996 9.3135200000000005 0.062888200000000005 0.02001 27.619299999999999 
    </InitialState>
  </Model>
  <ListOfTasks>
    <Task key="Task_14" name="Steady-State" type="steadyState" scheduled="false" updateModel="false">
      <Report reference="Report_9" target="" append="1" confirmOverwrite="1"/>
      <Problem>
        <Parameter name="JacobianRequested" type="bool" value="1"/>
        <Parameter name="StabilityAnalysisRequested" type="bool" value="1"/>
      </Problem>
      <Method name="Enhanced Newton" type="EnhancedNewton">
        <Parameter name="Resolution" type="unsignedFloat" value="1.0000000000000001e-09"/>
        <Parameter name="Derivation Factor" type="unsignedFloat" value="0.001"/>
        <Parameter name="Use Newton" type="bool" value="1"/>
        <Parameter name="Use Integration" type="bool" value="1"/>
        <Parameter name="Use Back Integration" type="bool" value="0"/>
        <Parameter name="Accept Negative Concentrations" type="bool" value="0"/>
        <Parameter name="Iteration Limit" type="unsignedInteger" value="50"/>
        <Parameter name="Maximum duration for forward integration" type="unsignedFloat" value="1000000000"/>
        <Parameter name="Maximum duration for backward integration" type="unsignedFloat" value="1000000"/>
      </Method>
    </Task>
    <Task key="Task_15" name="Time-Course" type="timeCourse" scheduled="false" updateModel="false">
      <Problem>
        <Parameter name="AutomaticStepSize" type="bool" value="0"/>
        <Parameter name="StepNumber" type="unsignedInteger" value="100"/>
        <Parameter name="StepSize" type="float" value="0.01"/>
        <Parameter name="Duration" type="float" value="1"/>
        <Parameter name="TimeSeriesRequested" type="bool" value="1"/>
        <Parameter name="OutputStartTime" type="float" value="0"/>
        <Parameter name="Output Event" type="bool" value="0"/>
        <Parameter name="Start in Steady State" type="bool" value="0"/>
      </Problem>
      <Method name="Deterministic (LSODA)" type="Deterministic(LSODA)">
        <Parameter name="Integrate Reduced Model" type="bool" value="0"/>
        <Parameter name="Relative Tolerance" type="unsignedFloat" value="9.9999999999999995e-07"/>
        <Parameter name="Absolute Tolerance" type="unsignedFloat" value="9.9999999999999998e-13"/>
        <Parameter name="Max Internal Steps" type="unsignedInteger" value="10000"/>
        <Parameter name="Max Internal Step Size" type="unsignedFloat" value="0"/>
      </Method>
    </Task>
    <Task key="Task_16" name="Scan" type="scan" scheduled="false" updateModel="false">
      <Problem>
        <Parameter name="Subtask" type="unsignedInteger" value="1"/>
        <ParameterGroup name="ScanItems">
        </ParameterGroup>
        <Parameter name="Output in subtask" type="bool" value="1"/>
        <Parameter name="Adjust initial conditions" type="bool" value="0"/>
      </Problem>
      <Method name="Scan Framework" type="ScanFramework">
      </Method>
    </Task>
    <Task key="Task_17" name="Elementary Flux Modes" type="fluxMode" scheduled="false" updateModel="false">
      <Report reference="Report_10" target="" append="1" confirmOverwrite="1"/>
      <Problem>
      </Problem>
      <Method name="EFM Algorithm" type="EFMAlgorithm">
      </Method>
    </Task>
    <Task key="Task_18" name="Optimization" type="optimization" scheduled="false" updateModel="false">
      <Report reference="Report_11" target="" append="1" confirmOverwrite="1"/>
      <Problem>
        <Parameter name="Subtask" type="cn" value="CN=Root,Vector=TaskList[Steady-State]"/>
        <ParameterText name="ObjectiveExpression" type="expression">
          
        </ParameterText>
        <Parameter name="Maximize" type="bool" value="0"/>
        <Parameter name="Randomize Start Values" type="bool" value="0"/>
        <Parameter name="Calculate Statistics" type="bool" value="1"/>
        <ParameterGroup name="OptimizationItemList">
        </ParameterGroup>
        <ParameterGroup name="OptimizationConstraintList">
        </ParameterGroup>
      </Problem>
      <Method name="Random Search" type="RandomSearch">
        <Parameter name="Log Verbosity" type="unsignedInteger" value="0"/>
        <Parameter name="Number of Iterations" type="unsignedInteger" value="100000"/>
        <Parameter name="Random Number Generator" type="unsignedInteger" value="1"/>
        <Parameter name="Seed" type="unsignedInteger" value="0"/>
      </Method>
    </Task>
    <Task key="Task_19" name="Parameter Estimation" type="parameterFitting" scheduled="false" updateModel="false">
      <Report reference="Report_12" target="" append="1" confirmOverwrite="1"/>
      <Problem>
        <Parameter name="Maximize" type="bool" value="0"/>
        <Parameter name="Randomize Start Values" type="bool" value="0"/>
        <Parameter name="Calculate Statistics" type="bool" value="1"/>
        <ParameterGroup name="OptimizationItemList">
        </ParameterGroup>
        <ParameterGroup name="OptimizationConstraintList">
        </ParameterGroup>
        <Parameter name="Steady-State" type="cn" value="CN=Root,Vector=TaskList[Steady-State]"/>
        <Parameter name="Time-Course" type="cn" value="CN=Root,Vector=TaskList[Time-Course]"/>
        <Parameter name="Create Parameter Sets" type="bool" value="0"/>
        <ParameterGroup name="Experiment Set">
        </ParameterGroup>
        <ParameterGroup name="Validation Set">
          <Parameter name="Weight" type="unsignedFloat" value="1"/>
          <Parameter name="Threshold" type="unsignedInteger" value="5"/>
        </ParameterGroup>
      </Problem>
      <Method name="Evolutionary Programming" type="EvolutionaryProgram">
        <Parameter name="Log Verbosity" type="unsignedInteger" value="0"/>
        <Parameter name="Number of Generations" type="unsignedInteger" value="200"/>
        <Parameter name="Population Size" type="unsignedInteger" value="20"/>
        <Parameter name="Random Number Generator" type="unsignedInteger" value="1"/>
        <Parameter name="Seed" type="unsignedInteger" value="0"/>
        <Parameter name="Stop after # Stalled Generations" type="unsignedInteger" value="0"/>
      </Method>
    </Task>
    <Task key="Task_20" name="Metabolic Control Analysis" type="metabolicControlAnalysis" scheduled="false" updateModel="false">
      <Report reference="Report_13" target="" append="1" confirmOverwrite="1"/>
      <Problem>
        <Parameter name="Steady-State" type="key" value="Task_14"/>
      </Problem>
      <Method name="MCA Method (Reder)" type="MCAMethod(Reder)">
        <Parameter name="Modulation Factor" type="unsignedFloat" value="1.0000000000000001e-09"/>
        <Parameter name="Use Reder" type="bool" value="1"/>
        <Parameter name="Use Smallbone" type="bool" value="1"/>
      </Method>
    </Task>
    <Task key="Task_21" name="Lyapunov Exponents" type="lyapunovExponents" scheduled="false" updateModel="false">
      <Report reference="Report_14" target="" append="1" confirmOverwrite="1"/>
      <Problem>
        <Parameter name="ExponentNumber" type="unsignedInteger" value="3"/>
        <Parameter name="DivergenceRequested" type="bool" value="1"/>
        <Parameter name="TransientTime" type="float" value="0"/>
      </Problem>
      <Method name="Wolf Method" type="WolfMethod">
        <Parameter name="Orthonormalization Interval" type="unsignedFloat" value="1"/>
        <Parameter name="Overall time" type="unsignedFloat" value="1000"/>
        <Parameter name="Relative Tolerance" type="unsignedFloat" value="9.9999999999999995e-07"/>
        <Parameter name="Absolute Tolerance" type="unsignedFloat" value="9.9999999999999998e-13"/>
        <Parameter name="Max Internal Steps" type="unsignedInteger" value="10000"/>
      </Method>
    </Task>
    <Task key="Task_22" name="Time Scale Separation Analysis" type="timeScaleSeparationAnalysis" scheduled="false" updateModel="false">
      <Report reference="Report_15" target="" append="1" confirmOverwrite="1"/>
      <Problem>
        <Parameter name="StepNumber" type="unsignedInteger" value="100"/>
        <Parameter name="StepSize" type="float" value="0.01"/>
        <Parameter name="Duration" type="float" value="1"/>
        <Parameter name="TimeSeriesRequested" type="bool" value="1"/>
        <Parameter name="OutputStartTime" type="float" value="0"/>
      </Problem>
      <Method name="ILDM (LSODA,Deuflhard)" type="TimeScaleSeparation(ILDM,Deuflhard)">
        <Parameter name="Deuflhard Tolerance" type="unsignedFloat" value="0.0001"/>
      </Method>
    </Task>
    <Task key="Task_23" name="Sensitivities" type="sensitivities" scheduled="false" updateModel="false">
      <Report reference="Report_16" target="" append="1" confirmOverwrite="1"/>
      <Problem>
        <Parameter name="SubtaskType" type="unsignedInteger" value="1"/>
        <ParameterGroup name="TargetFunctions">
          <Parameter name="SingleObject" type="cn" value=""/>
          <Parameter name="ObjectListType" type="unsignedInteger" value="7"/>
        </ParameterGroup>
        <ParameterGroup name="ListOfVariables">
          <ParameterGroup name="Variables">
            <Parameter name="SingleObject" type="cn" value=""/>
            <Parameter name="ObjectListType" type="unsignedInteger" value="41"/>
          </ParameterGroup>
          <ParameterGroup name="Variables">
            <Parameter name="SingleObject" type="cn" value=""/>
            <Parameter name="ObjectListType" type="unsignedInteger" value="0"/>
          </ParameterGroup>
        </ParameterGroup>
      </Problem>
      <Method name="Sensitivities Method" type="SensitivitiesMethod">
        <Parameter name="Delta factor" type="unsignedFloat" value="0.001"/>
        <Parameter name="Delta minimum" type="unsignedFloat" value="9.9999999999999998e-13"/>
      </Method>
    </Task>
    <Task key="Task_24" name="Moieties" type="moieties" scheduled="false" updateModel="false">
      <Problem>
      </Problem>
      <Method name="Householder Reduction" type="Householder">
      </Method>
    </Task>
    <Task key="Task_25" name="Cross Section" type="crosssection" scheduled="false" updateModel="false">
      <Problem>
        <Parameter name="AutomaticStepSize" type="bool" value="0"/>
        <Parameter name="StepNumber" type="unsignedInteger" value="100"/>
        <Parameter name="StepSize" type="float" value="0.01"/>
        <Parameter name="Duration" type="float" value="1"/>
        <Parameter name="TimeSeriesRequested" type="bool" value="1"/>
        <Parameter name="OutputStartTime" type="float" value="0"/>
        <Parameter name="Output Event" type="bool" value="0"/>
        <Parameter name="Start in Steady State" type="bool" value="0"/>
        <Parameter name="LimitCrossings" type="bool" value="0"/>
        <Parameter name="NumCrossingsLimit" type="unsignedInteger" value="0"/>
        <Parameter name="LimitOutTime" type="bool" value="0"/>
        <Parameter name="LimitOutCrossings" type="bool" value="0"/>
        <Parameter name="PositiveDirection" type="bool" value="1"/>
        <Parameter name="NumOutCrossingsLimit" type="unsignedInteger" value="0"/>
        <Parameter name="LimitUntilConvergence" type="bool" value="0"/>
        <Parameter name="ConvergenceTolerance" type="float" value="9.9999999999999995e-07"/>
        <Parameter name="Threshold" type="float" value="0"/>
        <Parameter name="DelayOutputUntilConvergence" type="bool" value="0"/>
        <Parameter name="OutputConvergenceTolerance" type="float" value="9.9999999999999995e-07"/>
        <ParameterText name="TriggerExpression" type="expression">
          
        </ParameterText>
        <Parameter name="SingleVariable" type="cn" value=""/>
      </Problem>
      <Method name="Deterministic (LSODA)" type="Deterministic(LSODA)">
        <Parameter name="Integrate Reduced Model" type="bool" value="0"/>
        <Parameter name="Relative Tolerance" type="unsignedFloat" value="9.9999999999999995e-07"/>
        <Parameter name="Absolute Tolerance" type="unsignedFloat" value="9.9999999999999998e-13"/>
        <Parameter name="Max Internal Steps" type="unsignedInteger" value="10000"/>
        <Parameter name="Max Internal Step Size" type="unsignedFloat" value="0"/>
      </Method>
    </Task>
    <Task key="Task_26" name="Linear Noise Approximation" type="linearNoiseApproximation" scheduled="false" updateModel="false">
      <Report reference="Report_17" target="" append="1" confirmOverwrite="1"/>
      <Problem>
        <Parameter name="Steady-State" type="key" value=""/>
      </Problem>
      <Method name="Linear Noise Approximation" type="LinearNoiseApproximation">
      </Method>
    </Task>
  </ListOfTasks>
  <ListOfReports>
    <Report key="Report_9" name="Steady-State" taskType="steadyState" separator="&#x09;" precision="6">
      <Comment>
        Automatically generated report.
      </Comment>
      <Footer>
        <Object cn="CN=Root,Vector=TaskList[Steady-State]"/>
      </Footer>
    </Report>
    <Report key="Report_10" name="Elementary Flux Modes" taskType="fluxMode" separator="&#x09;" precision="6">
      <Comment>
        Automatically generated report.
      </Comment>
      <Footer>
        <Object cn="CN=Root,Vector=TaskList[Elementary Flux Modes],Object=Result"/>
      </Footer>
    </Report>
    <Report key="Report_11" name="Optimization" taskType="optimization" separator="&#x09;" precision="6">
      <Comment>
        Automatically generated report.
      </Comment>
      <Header>
        <Object cn="CN=Root,Vector=TaskList[Optimization],Object=Description"/>
        <Object cn="String=\[Function Evaluations\]"/>
        <Object cn="Separator=&#x09;"/>
        <Object cn="String=\[Best Value\]"/>
        <Object cn="Separator=&#x09;"/>
        <Object cn="String=\[Best Parameters\]"/>
      </Header>
      <Body>
        <Object cn="CN=Root,Vector=TaskList[Optimization],Problem=Optimization,Reference=Function Evaluations"/>
        <Object cn="Separator=&#x09;"/>
        <Object cn="CN=Root,Vector=TaskList[Optimization],Problem=Optimization,Reference=Best Value"/>
        <Object cn="Separator=&#x09;"/>
        <Object cn="CN=Root,Vector=TaskList[Optimization],Problem=Optimization,Reference=Best Parameters"/>
      </Body>
      <Footer>
        <Object cn="String=&#x0a;"/>
        <Object cn="CN=Root,Vector=TaskList[Optimization],Object=Result"/>
      </Footer>
    </Report>
    <Report key="Report_12" name="Parameter Estimation" taskType="parameterFitting" separator="&#x09;" precision="6">
      <Comment>
        Automatically generated report.
      </Comment>
      <Header>
        <Object cn="CN=Root,Vector=TaskList[Parameter Estimation],Object=Description"/>
        <Object cn="String=\[Function Evaluations\]"/>
        <Object cn="Separator=&#x09;"/>
        <Object cn="String=\[Best Value\]"/>
        <Object cn="Separator=&#x09;"/>
        <Object cn="String=\[Best Parameters\]"/>
      </Header>
      <Body>
        <Object cn="CN=Root,Vector=TaskList[Parameter Estimation],Problem=Parameter Estimation,Reference=Function Evaluations"/>
        <Object cn="Separator=&#x09;"/>
        <Object cn="CN=Root,Vector=TaskList[Parameter Estimation],Problem=Parameter Estimation,Reference=Best Value"/>
        <Object cn="Separator=&#x09;"/>
        <Object cn="CN=Root,Vector=TaskList[Parameter Estimation],Problem=Parameter Estimation,Reference=Best Parameters"/>
      </Body>
      <Footer>
        <Object cn="String=&#x0a;"/>
        <Object cn="CN=Root,Vector=TaskList[Parameter Estimation],Object=Result"/>
      </Footer>
    </Report>
    <Report key="Report_13" name="Metabolic Control Analysis" taskType="metabolicControlAnalysis" separator="&#x09;" precision="6">
      <Comment>
        Automatically generated report.
      </Comment>
      <Header>
        <Object cn="CN=Root,Vector=TaskList[Metabolic Control Analysis],Object=Description"/>
      </Header>
      <Footer>
        <Object cn="String=&#x0a;"/>
        <Object cn="CN=Root,Vector=TaskList[Metabolic Control Analysis],Object=Result"/>
      </Footer>
    </Report>
    <Report key="Report_14" name="Lyapunov Exponents" taskType="lyapunovExponents" separator="&#x09;" precision="6">
      <Comment>
        Automatically generated report.
      </Comment>
      <Header>
        <Object cn="CN=Root,Vector=TaskList[Lyapunov Exponents],Object=Description"/>
      </Header>
      <Footer>
        <Object cn="String=&#x0a;"/>
        <Object cn="CN=Root,Vector=TaskList[Lyapunov Exponents],Object=Result"/>
      </Footer>
    </Report>
    <Report key="Report_15" name="Time Scale Separation Analysis" taskType="timeScaleSeparationAnalysis" separator="&#x09;" precision="6">
      <Comment>
        Automatically generated report.
      </Comment>
      <Header>
        <Object cn="CN=Root,Vector=TaskList[Time Scale Separation Analysis],Object=Description"/>
      </Header>
      <Footer>
        <Object cn="String=&#x0a;"/>
        <Object cn="CN=Root,Vector=TaskList[Time Scale Separation Analysis],Object=Result"/>
      </Footer>
    </Report>
    <Report key="Report_16" name="Sensitivities" taskType="sensitivities" separator="&#x09;" precision="6">
      <Comment>
        Automatically generated report.
      </Comment>
      <Header>
        <Object cn="CN=Root,Vector=TaskList[Sensitivities],Object=Description"/>
      </Header>
      <Footer>
        <Object cn="String=&#x0a;"/>
        <Object cn="CN=Root,Vector=TaskList[Sensitivities],Object=Result"/>
      </Footer>
    </Report>
    <Report key="Report_17" name="Linear Noise Approximation" taskType="linearNoiseApproximation" separator="&#x09;" precision="6">
      <Comment>
        Automatically generated report.
      </Comment>
      <Header>
        <Object cn="CN=Root,Vector=TaskList[Linear Noise Approximation],Object=Description"/>
      </Header>
      <Footer>
        <Object cn="String=&#x0a;"/>
        <Object cn="CN=Root,Vector=TaskList[Linear Noise Approximation],Object=Result"/>
      </Footer>
    </Report>
  </ListOfReports>
  <GUI>
  </GUI>
  <SBMLReference file="MODEL1505110000_url.xml">
    <SBMLMap SBMLid="ADP" COPASIkey="Metabolite_50"/>
    <SBMLMap SBMLid="AKG" COPASIkey="Metabolite_3"/>
    <SBMLMap SBMLid="AMP" COPASIkey="Metabolite_51"/>
    <SBMLMap SBMLid="ATP" COPASIkey="Metabolite_52"/>
    <SBMLMap SBMLid="ATP_MAINTENANCE" COPASIkey="Reaction_60"/>
    <SBMLMap SBMLid="BPG" COPASIkey="Metabolite_4"/>
    <SBMLMap SBMLid="DAP" COPASIkey="Metabolite_6"/>
    <SBMLMap SBMLid="ENO" COPASIkey="Reaction_7"/>
    <SBMLMap SBMLid="F6P" COPASIkey="Metabolite_8"/>
    <SBMLMap SBMLid="F6P_GAP_TAL" COPASIkey="Reaction_17"/>
    <SBMLMap SBMLid="FBA" COPASIkey="Reaction_2"/>
    <SBMLMap SBMLid="FBP" COPASIkey="Reaction_19"/>
    <SBMLMap SBMLid="FDP" COPASIkey="Metabolite_9"/>
    <SBMLMap SBMLid="FEED" COPASIkey="ModelValue_0"/>
    <SBMLMap SBMLid="Function_for_ENO" COPASIkey="Function_44"/>
    <SBMLMap SBMLid="Function_for_F6P_GAP_TAL" COPASIkey="Function_61"/>
    <SBMLMap SBMLid="Function_for_FBA" COPASIkey="Function_52"/>
    <SBMLMap SBMLid="Function_for_FBP" COPASIkey="Function_46"/>
    <SBMLMap SBMLid="Function_for_GDH" COPASIkey="Function_42"/>
    <SBMLMap SBMLid="Function_for_GLC_FEED" COPASIkey="Function_73"/>
    <SBMLMap SBMLid="Function_for_GPM" COPASIkey="Function_62"/>
    <SBMLMap SBMLid="Function_for_PFK" COPASIkey="Function_54"/>
    <SBMLMap SBMLid="Function_for_PGK" COPASIkey="Function_41"/>
    <SBMLMap SBMLid="Function_for_PPS" COPASIkey="Function_50"/>
    <SBMLMap SBMLid="Function_for_PTS_0" COPASIkey="Function_70"/>
    <SBMLMap SBMLid="Function_for_PYK" COPASIkey="Function_55"/>
    <SBMLMap SBMLid="Function_for_TPI" COPASIkey="Function_56"/>
    <SBMLMap SBMLid="G6P" COPASIkey="Metabolite_11"/>
    <SBMLMap SBMLid="GAP" COPASIkey="Metabolite_12"/>
    <SBMLMap SBMLid="GDH" COPASIkey="Reaction_4"/>
    <SBMLMap SBMLid="GL6P" COPASIkey="Metabolite_13"/>
    <SBMLMap SBMLid="GL6P_HYDROLYSIS" COPASIkey="Reaction_67"/>
    <SBMLMap SBMLid="GLC_feed" COPASIkey="Reaction_55"/>
    <SBMLMap SBMLid="GLCp" COPASIkey="Metabolite_74"/>
    <SBMLMap SBMLid="GLCx" COPASIkey="Metabolite_70"/>
    <SBMLMap SBMLid="GPM" COPASIkey="Reaction_6"/>
    <SBMLMap SBMLid="KdADPMg" COPASIkey="ModelValue_1"/>
    <SBMLMap SBMLid="KdATPMg" COPASIkey="ModelValue_2"/>
    <SBMLMap SBMLid="KdFDPMg" COPASIkey="ModelValue_3"/>
    <SBMLMap SBMLid="KeqNDH" COPASIkey="ModelValue_7"/>
    <SBMLMap SBMLid="KmACO_ACN" COPASIkey="ModelValue_6"/>
    <SBMLMap SBMLid="KmCIT_ACN" COPASIkey="ModelValue_5"/>
    <SBMLMap SBMLid="KmICIT_ACN" COPASIkey="ModelValue_4"/>
    <SBMLMap SBMLid="MG" COPASIkey="Metabolite_57"/>
    <SBMLMap SBMLid="MgADP" COPASIkey="Metabolite_58"/>
    <SBMLMap SBMLid="MgATP" COPASIkey="Metabolite_59"/>
    <SBMLMap SBMLid="MgFDP" COPASIkey="Metabolite_60"/>
    <SBMLMap SBMLid="NAD" COPASIkey="Metabolite_18"/>
    <SBMLMap SBMLid="NADH" COPASIkey="Metabolite_19"/>
    <SBMLMap SBMLid="OAA" COPASIkey="Metabolite_22"/>
    <SBMLMap SBMLid="P" COPASIkey="Metabolite_56"/>
    <SBMLMap SBMLid="PEP" COPASIkey="Metabolite_23"/>
    <SBMLMap SBMLid="PFK" COPASIkey="Reaction_1"/>
    <SBMLMap SBMLid="PGA2" COPASIkey="Metabolite_24"/>
    <SBMLMap SBMLid="PGA3" COPASIkey="Metabolite_25"/>
    <SBMLMap SBMLid="PGI" COPASIkey="Reaction_0"/>
    <SBMLMap SBMLid="PGK" COPASIkey="Reaction_5"/>
    <SBMLMap SBMLid="PGN" COPASIkey="Metabolite_26"/>
    <SBMLMap SBMLid="PPS" COPASIkey="Reaction_22"/>
    <SBMLMap SBMLid="PTS_0" COPASIkey="Reaction_50"/>
    <SBMLMap SBMLid="PTS_1" COPASIkey="Reaction_51"/>
    <SBMLMap SBMLid="PTS_2" COPASIkey="Reaction_52"/>
    <SBMLMap SBMLid="PTS_3" COPASIkey="Reaction_53"/>
    <SBMLMap SBMLid="PTS_4" COPASIkey="Reaction_54"/>
    <SBMLMap SBMLid="PYK" COPASIkey="Reaction_8"/>
    <SBMLMap SBMLid="PYR" COPASIkey="Metabolite_27"/>
    <SBMLMap SBMLid="R5P" COPASIkey="Metabolite_30"/>
    <SBMLMap SBMLid="RU5P" COPASIkey="Metabolite_31"/>
    <SBMLMap SBMLid="S7P" COPASIkey="Metabolite_32"/>
    <SBMLMap SBMLid="SUCCOA" COPASIkey="Metabolite_34"/>
    <SBMLMap SBMLid="TPI" COPASIkey="Reaction_3"/>
    <SBMLMap SBMLid="X5P" COPASIkey="Metabolite_35"/>
    <SBMLMap SBMLid="XCH_GLC" COPASIkey="Reaction_61"/>
    <SBMLMap SBMLid="cell" COPASIkey="Compartment_0"/>
    <SBMLMap SBMLid="cell_periplasm" COPASIkey="Compartment_2"/>
    <SBMLMap SBMLid="ei" COPASIkey="Metabolite_36"/>
    <SBMLMap SBMLid="eiP" COPASIkey="Metabolite_37"/>
    <SBMLMap SBMLid="eiia" COPASIkey="Metabolite_38"/>
    <SBMLMap SBMLid="eiiaP" COPASIkey="Metabolite_39"/>
    <SBMLMap SBMLid="eiicb" COPASIkey="Metabolite_40"/>
    <SBMLMap SBMLid="eiicbP" COPASIkey="Metabolite_41"/>
    <SBMLMap SBMLid="extracellular" COPASIkey="Compartment_1"/>
    <SBMLMap SBMLid="hpr" COPASIkey="Metabolite_42"/>
    <SBMLMap SBMLid="hprP" COPASIkey="Metabolite_43"/>
    <SBMLMap SBMLid="tal" COPASIkey="Metabolite_46"/>
    <SBMLMap SBMLid="talC3" COPASIkey="Metabolite_47"/>
  </SBMLReference>
  <ListOfUnitDefinitions>
    <UnitDefinition key="Unit_0" name="meter" symbol="m">
      <Expression>
        m
      </Expression>
    </UnitDefinition>
    <UnitDefinition key="Unit_2" name="second" symbol="s">
      <Expression>
        s
      </Expression>
    </UnitDefinition>
    <UnitDefinition key="Unit_6" name="Avogadro" symbol="Avogadro">
      <Expression>
        Avogadro
      </Expression>
    </UnitDefinition>
    <UnitDefinition key="Unit_8" name="item" symbol="#">
      <Expression>
        #
      </Expression>
    </UnitDefinition>
    <UnitDefinition key="Unit_17" name="liter" symbol="l">
      <Expression>
        0.001*m^3
      </Expression>
    </UnitDefinition>
    <UnitDefinition key="Unit_20" name="mole" symbol="mol">
      <Expression>
        Avogadro*#
      </Expression>
    </UnitDefinition>
  </ListOfUnitDefinitions>
</COPASI>
