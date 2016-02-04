{ //LIBRARY_STRATEGY[text()] }	{ //LIBRARY_LAYOUT/* }	{ //SPOT_LENGTH[text()] }	{ //XREF_LINK/DB[text()=pubmed]/../ID[text()] }	{
for $x in //SAMPLE_ATTRIBUTE/*
return { TAG[@text] }={ VALUE[@text] }; 
}
