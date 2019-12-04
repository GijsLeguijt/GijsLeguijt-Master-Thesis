void rootmacro()
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
			<< "41: (r,z) of hits"			<< endl;
	cin		>> i;
	cout	<< endl;

	//Selecting if cuts have to be made (no choices on which cut yet)
	TString input_cut;
	cout	<< "For cut choose 'y', else choose 'n':" << endl;
	cin		>> input_cut;
	cout	<< endl;

	//Setting the cut
	TString cut;
	if (input_cut == 'y')
	{	TString c_edep   = "etot > 0";
		TString c_nodisc = "abs(zp_pri) < 750";
	
		cut  = "w_pri * ( " + c_edep + " && " + c_nodisc + " )";
		cout << "Cut: " << cut << endl;
	}
	else
		cut  = "w_pri";
	
	//Radial distance of particle and primary
	TString rp 		= "sqrt(xp^2 + yp^2)";
	TString rp_pri 	= "sqrt(xp_pri^2 + yp_pri^2)";
	
	//Plotting	
	switch (i)
	{
		case 10: evt_vrt->Draw("xp_pri:yp_pri:zp_pri");						break;
		case 11: evt_vrt->Draw("xp:yp:zp");				 					break;
		case 20: evt_vrt->Draw("xp_pri:yp_pri",cut); 						break;
		case 21: evt_vrt->Draw("xp:yp",cut); 								break;
		case 30: evt_vrt->Draw(rp_pri,cut); 						        break;
		case 31: evt_vrt->Draw(rp,cut);		 						        break;
		case 40: evt_vrt->Draw("zp_pri:"+rp_pri,"w_pri * (ed > 0)","colz");	break;
		case 41: evt_vrt->Draw("zp:"+rp,"w_pri * (ed > 0)","colz");			break;
	}
}


