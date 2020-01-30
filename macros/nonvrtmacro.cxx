void nonvrtmacro()
{	
	//Selecting which plot to make	
	int i;
	cout 	<< "Choose:" 					<< endl
			<< "10: (x,y,z) of primaries"	<< endl
			<< "11: (x,y,z) of hits"		<< endl
			<< "20: (x,y) of primaries"		<< endl
			<< "21: (x,y) of hits"			<< endl
			<< "30: r of primaries"	 		<< endl
			<< "31: r of hits"		 		<< endl
			<< "40: (r,z) of primaries"		<< endl
			<< "41: (r,z) of hits"			<< endl
			<< "42: (r2,z) of hits"			<< endl
			<< "50: energy spectrum"		<< endl;
	cin		>> i;
	cout	<< endl;

	//Selecting if cuts have to be made (no choices on which cut yet)
	TString input_cut;
	cout	<< "For FV-cut choose 'y', else choose 'n':" << endl;
	cin		>> input_cut;
	cout	<< endl;

	//Setting the cut
	TString cut;
	if (input_cut == 'y')
	{	TString c_edep   = "ed > 0";
		TString c_fidz = "abs(zp_pri) < 670";
        TString c_fidr = "sqrt(xp^2 + yp^2) < 570";
	
		cut  = "( " + c_edep + " && " + c_fidz + " && " + c_fidr + " )";
		cout << "Cut: " << cut << endl;
	}
    else
    {   cut = "ed > 0";
    }
    
	
	//Radial distance of particle and primary
	TString rp 		= "sqrt(xp^2 + yp^2)";
	TString rp2     = "xp^2 + yp^2";
	TString rp_pri 	= "sqrt(xp_pri^2 + yp_pri^2)";
	
	//Plotting	
	switch (i)
	{
		case 10: evt->Draw("xp_pri:yp_pri:zp_pri");				break;
		case 11: evt->Draw("xp:yp:zp");				 			break;
		case 20: evt->Draw("xp_pri:yp_pri",cut); 				break;
		case 21: evt->Draw("xp:yp",cut); 						break;
		case 30: evt->Draw(rp_pri,cut); 						break;
		case 31: evt->Draw(rp,cut);		 						break;
		case 40: evt->Draw("zp_pri:"+rp_pri,cut,"colz");	    break;
		case 41: evt->Draw("zp:"+rp,cut,"colz");			    break;
		case 42: evt->Draw("zp:"+rp2,cut,"colz");		        break;
		case 50: evt->Draw("ed",cut,"hist");					break;
	}
}


