# Tunnel-Ionization
=================================

Program Description:

    1) This code is based on Dr. Tenio Popmintchev's PhD thesis
    2) Most of the work accompolished by Peter Frame, Siyang Wang ,
    Shiwen An from Tenio's Group
    3) Code Structure:
    main.m => ArgonPlotFunc
                 => HeliumPlotFunc
                 => NeonPlotFunc
    WIth specific function use under each function
    Will automatically draw the function with certain range based on
    [Wavelength, Intensity, Gas/Ion, Cycles, Sample rate]
    
    4) Problems happened and wrong doings:
       a)Sample Rate too low: The diagram will be distorted under this
       circumstance because the 8th~9th even higher derivative might
       require higher sample rate to get the 9th derivative numerically.
       b)Keldysh Parameter or Cutoff energy distorted on the diagram
       c) Similar or the same color among various plot [still fixing it]
    
Contact Info:
    http://popmintchev.ucsd.edu/team/
