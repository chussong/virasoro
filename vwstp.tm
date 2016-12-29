void RunFromComponents P((char*, char*, char*, char*, char*));

:Begin:
:Function:		RunFromComponents
:Pattern:		VPass[c_String, hl_String, hh_String, hp_String, maxOrder_String]
:Arguments:		{c, hl, hh, hp, maxOrder}
:ArgumentTypes:	{String, String, String, String, String}
:ReturnType:	Manual
:End:

:Evaluate: VPass::usage = "VPass[c, hl, hh, hp, maxOrder] sends the run to Virasoro for computation, returning the Mathematica object that would be given by VRead[]. Arguments must be strings!"
