classdef hNRUplinkGrantFormat
%hNRUplinkGrantFormat Represents the uplink grant information
% The uplink assignment information does not include all the fields of
% downlink control information (DCI) format 0_1 as per 3GPP standard. Grant
% packets containing only the information fields which feature as member of
% this class are assumed to be exchanged with UEs.

%   Copyright 2019-2020 The MathWorks, Inc.

%#codegen

   properties

      %SlotOffset Offset of the allocated slot from the current slot (k2)
      SlotOffset = []

      %RBGAllocationBitmap Resource block group(RBG) allocation represented as bit vector
      RBGAllocationBitmap = []

      %StartSymbol Location of first symbol
      StartSymbol = []

      %NumSymbols Number of symbols
      NumSymbols = []

      %MCS Modulation and coding scheme
      MCS = []

      %NDI New data indicator flag
      NDI = []

      %RV Redundancy version sequence number
      RV = []

      %HARQID HARQ process identifier
      HARQID = []

      %DMRSLength DM-RS length 
      DMRSLength = []

      %MappingType Mapping type 
      MappingType = []
      
      %NumLayers Number of layers (1...4)
      NumLayers = 1
   end
end