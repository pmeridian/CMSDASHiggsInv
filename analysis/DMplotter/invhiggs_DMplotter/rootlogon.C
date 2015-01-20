{
    // Set up FW Lite for automatic loading of CMS libraries
    // and data formats.   As you may have other user-defined setup
    // in your rootlogon.C, the CMS setup is executed only if the CMS
    // environment is set up.
    //

    
/*
      TString cmsswbase = getenv("CMSSW_BASE");
      if (cmsswbase.Length() > 0) {
        //
        // The CMSSW environment is defined (this is true even for FW Lite)
        // so set up the rest.
        //
        cout << "Loading FW Lite setup." << endl;
        gSystem->Load("libFWCoreFWLite.so");
        AutoLibraryLoader::enable();
        gSystem->Load("libDataFormatsFWLite.so");
        gSystem->Load("libDataFormatsPatCandidates.so");

 	gSystem->Load("libRooFit") ;
	using namespace RooFit ;
       }  

*/
    //gSystem->Load("libRooFit");
    //using namespace RooFit;

    // if root complains about headers not found, try "gSystem->SetIncludePath("-I$ROOFITSYS/include")"
    //gSystem->SetIncludePath("-I$ROOFITSYS/include");

    //rootlogon.C
    printf("\033[32m\n");
    printf("######################################### \n");
    printf("#  WELCOME to RENJIE WANG's rootlogo    # \n");
    printf("#          renjie.wang@cern.ch          # \n");
    printf("#             Apr. 10, 2013             # \n");
    printf("######################################### \n");
    printf("\033[0m\n");
    //printf("\n >>>>>> WELCOME to RENJIE WANG rootlogon <<<<<< \n\n");     // print a message


    gROOT->SetStyle("Plain");                       // plain histogram style
    gStyle->SetOptStat("nemruoi");                  // expanded stats box
    gStyle->SetPadTickX(1);                         // tic marks on all axes
    gStyle->SetPadTickY(1);                         //

//gStyle->SetOptFit(1111);                        // the results of the fits
//gStyle->SetOptFit();
    gStyle->SetFillStyle(0);
    gStyle->SetTitleBorderSize(0);
    gStyle->SetLegendBorderSize(0);
//gStyle->SetLineWidth(1.0);

//gStyle->SetPadGridX(kTRUE);                     // draw horizontal and vertical grids
//gStyle->SetPadGridY(kTRUE);
    gStyle->SetPalette(1,0);                        // pretty and useful palette
    gStyle->SetHistLineWidth(2);                    // a thicker histogram line
    gStyle->SetFrameFillColor(10);                  // a different frame colour
    gStyle->SetTitleFillColor(33);                 // title colour to highlight it
    gStyle->SetTitleW(.77);                         // Title Width
    gStyle->SetTitleH(.07);                        // Title height
    gStyle->SetHistMinimumZero(true);               // Suppresses Histogram Zero Supression

    //gROOT->Reset();
    gROOT->LoadMacro("$HOME/tdrstyle.C");
    setTDRStyle();

}
