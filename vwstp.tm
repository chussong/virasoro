void RunFromComponents P((const char*, const char*, const char*, const char*, const char*));

:Begin:
:Function:		RunFromComponents
:Pattern:		VPass[c_String, hl_String, hh_String, hp_String, maxOrder_String]
:Arguments:		{c, hl, hh, hp, maxOrder}
:ArgumentTypes:	{String, String, String, String, String}
:ReturnType:	Manual
:End:

:Evaluate: VPass::usage = "VPass[c, hl, hh, hp, maxOrder] or sends the run to Virasoro for computation, returning the Mathematica object that would be given by VRead[]. Arguments must be strings!"

void RunFromFile P((const char*));

:Begin:
:Function:		RunFromFile
:Pattern:		VPassFilename[filename_String]
:Arguments:		{filename}
:ArgumentTypes:	{String}
:ReturnType:	Manual
:End:

:Evaluate: VPassFilename::usage = "VPassFilename[filename] tells Virasoro to compute the given runfile, returning the Mathematica object that would be given by VRead[]. Argument must be a string!"
