%c
%c                                                       ..'´`'..'´`..'´`..                                                   
%c       File: Data_Structure_ONDA
%c
%c     A comprehensive description of ONDA's data structure. 
%c
%c                                           by M.Segatto
%c                                           15/04/2015
%c                                           segatto@ele.ufes.br
%c 
%c     References:
%c       [1] http://en.wikipedia.org/wiki/Data_structure.
%c
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c                    G L O B A L  V A R I A B L E S                      %
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c   None by now.
%c
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c                   L I S T  O F  V A R I A B L E S                      %
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c   
%c   alpha             : Fibre attenuation                           [1/km]
%c   alphadB           : Fibre attenuation in dB                    [dB/km]
%c   Aeff              : Area efetiva da fibra                        [m^2]
%c   betha2            : 2nd order dispersion coefficient           [s^2/m]
%c   betha3            : 3nd order dispersion coefficient           [s^3/m]
%c   Bit_Rate          : Transmission rate                         [Gbit/s]
%c   gamma             : Fibre nonlinear coefficient [W/??] ## REVER
%c   n2                : Indice de refracao nao linear              [m^2/W]
%c   frequency         : Frequency                                     [Hz]
%c   frequencynorm     : Normalised frequency (see optilux)
%c   lambda            : Wavelength                                     [m]
%c   mod_type          : Modulation format (m-ASK, m-PSK, m-QAM,...)
%c   n2                : Fibre nonlinear refractive index           [m^2/W]
%c   NBITS             : Number of bits
%c   NSYMBOLS          : Number of symbos
%c   omega             : 2 *pi*frequency                            [rad/s]
%c   OSNR              : Opitical signal to noise ratio                [dB]
%c   Sample_Rate       : Number of point per bit
%c   SNR               : Electrical signal to noise ratio              [dB]
%c   Symbol_Rate       : Symbol rate                                [GBaud]
%c   t                 : Time                                           [s]




%c
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c                             S T R U C T S                              %
%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%c
%c   Amplifier.
%c                   .id
%c                   .id_Node
%c                   .Gain
%c                   .Noise_Figure
%c                   .Product_Code
%c                   .type          : Must be 'electrical' or 'optical'
%c   Band.
%c   Optical_Channel.
%c                   .id
%c                   .id_Route
%c                   .lambda_used
%c                   .lambda.viable
%c   DCM.
%c   Demand.
%c                   .id
%c                   .id_Node_Source
%c                   .id_Node_Destination
%c                   .id_Optical_Channel
%c                   .type           : IMDD, Coherent, EOFDM, OOFDM 
%c                                     or FlexGrid.
%c                                     IMDD -  Intensity-Modulation 
%c                                             Direct Detection
%c                                     EOFDM - Electrically generated 
%c                                             optical OFDM
%c                                     OOFDM - Optically generated 
%c                                             optical OFDM
%c                   .value          : [Gbit/s]
%c                   
%c
%c   Fibre_Segment. 
%c   Interface.
%c   Link.
%c   Network.   
%c   Node.
%c   OADM_Terminal.
%c   Route.