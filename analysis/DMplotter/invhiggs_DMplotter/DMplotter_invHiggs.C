//
// Renjie Wang 
// (renjie.wang@cern.ch)
// Aug 2013
//
// Fig13 in arXiv:1404.1344
//
// usesage: root -l -b -q DMplotter_invHiggs.C 
// 	    plug rootlogon.C in ~/.rootrc  

#include <stdio>
#include <stdlib>
#include <iostream>
#include <fstream>
#include <string>
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TString.h"


using namespace std;
double const Mn = 0.93895;
double const Width_SM = 0.00407;
double const v = 174;

/*FIXME*/
double const CLval=90;
//published
//double const BRinv = 0.51;
//parked
double const BRinv = 0.40;

//double const CLval=95;
//double const BRinv = 0.58;

double const Mh = 125;
/*FIXME*/


//double const BRinv = 0.58;
//double const Mh = 105;
//double const BRinv = 0.61;
//double const Mh = 115;
//double const BRinv = 0.80;
//double const Mh = 135;
//double const BRinv = 0.84;
//double const Mh = 145;


double scalar(double *x, double fn)
{
    if(1 < 2*x[0]/Mh) return 0;
    double c1 = (Width_SM*BRinv)/(1.-BRinv);
    double beta = sqrt(1-4*pow(x[0]/Mh,2));
    double c2 = 4*c1*pow(Mn,4)*fn*fn;
    c2 /= v*v*beta*pow(Mh,3)*(x[0]+Mn)*(x[0]+Mn);
    return 0.3894*c2*1.0e+9;
}


double vector(double *x, double fn)
{
    if(1 < 2*x[0]/Mh) return 0;
    double c1 = (Width_SM*BRinv)/(1.-BRinv);
    double beta = sqrt(1-4*pow(x[0]/Mh,2));
    double c2 = c1*16*pow(x[0]*Mn,4)*fn*fn;
    c2 /= v*v*beta*pow(Mh,7)*(x[0]+Mn)*(x[0]+Mn);
    c2 /= 1-4*pow(x[0]/Mh,2)+12*pow(x[0]/Mh,4);
    return 0.3894*c2*1.0e+9;
}


double fermion(double *x, double fn)
{
    if(1 < 2*x[0]/Mh) return 0;
    double c1 = (Width_SM*BRinv)/(1.-BRinv);
    double beta = sqrt(1-4*pow(x[0]/Mh,2));
    double c2 = c1*8*pow(x[0]*Mn*fn,2)*Mn*Mn;
    c2 /= v*v*pow(beta,3)*pow(Mh,5)*(x[0]+Mn)*(x[0]+Mn);
    return 0.3894*c2*1.0e+9;
}


enum ModelType {SCALAR,VECTOR,FERMION};
TGraph *MakeGraph(int Type=SCALAR, double fn)
{

    double const k = 2;
    double const numb = Mh*5*k;
    Double_t x[numb], y[numb];
    Int_t n = numb;
    for (Int_t i=0; i<n; i++) {
        x[i] = i*0.1*1/k;
        if(Type == SCALAR) 	y[i] = scalar(&x[i], fn);
        else if (Type == VECTOR) 	y[i] = vector(&x[i], fn);
        else if (Type == FERMION)  y[i] = fermion(&x[i], fn);
        else 			y[i] = 0;
    }
    TGraph *gr = new TGraph(n,x,y);
    gr->GetXaxis()->SetLimits(4,1000);
    return gr;
}

void DMplotter_invHiggs()
{

    TCanvas *canv = new TCanvas("canv", "limits canvas", 800., 600.);
    TPad* t1d = new TPad();
    t1d = new TPad("t1d","t1d", 0.0, 0.0, 1.0, 1.0);
    t1d->Draw();
    t1d->SetRightMargin(0.03);
    t1d->cd();

    //t1d->SetGridx(1);
    //t1d->SetGridy(1);
    t1d->SetLogy();
    t1d->SetLogx();

    //double const fn = 0.629;//0.326;
    TGraph *h_scalar_min = MakeGraph(SCALAR,0.260);
    TGraph *h_scalar_lat = MakeGraph(SCALAR,0.326);
    TGraph *h_scalar_max = MakeGraph(SCALAR,0.629);
    TGraph *h_vector_min = MakeGraph(VECTOR,0.260);
    TGraph *h_vector_lat = MakeGraph(VECTOR,0.326);
    TGraph *h_vector_max = MakeGraph(VECTOR,0.629);
    TGraph *h_fermion_min = MakeGraph(FERMION,0.260);
    TGraph *h_fermion_lat = MakeGraph(FERMION,0.326);
    TGraph *h_fermion_max = MakeGraph(FERMION,0.629);

//////////////////////////////////////////////////
    TGraph *g1 = new TGraph("CRESST_II_2sigma_pt2.dat","%lg %lg");
    TGraph *g2 = new TGraph("CRESST_II_1sigma.dat","%lg %lg");
    TGraph *g3 = new TGraph("CRESST_II_2sigma_pt1.dat","%lg %lg");
    
    TGraph *h1 = new TGraph("XENON100.dat","%lg %lg");
    TGraph *h2 = new TGraph("XENON10_2011.dat","%lg %lg");

    TGraph *f1 = new TGraph("DAMA_08_Savage_3sig_pt1.dat","%lg %lg");
    TGraph *f2 = new TGraph("DAMA_08_Savage_3sig_pt2.dat","%lg %lg");

    TGraph *k1 = new TGraph("CoGeNT_90CL.dat","%lg %lg");
    TGraph *k2 = new TGraph("CoGeNT_99CL.dat","%lg %lg");

    TGraph *j1 = new TGraph("CDMS_2010.dat","%lg %lg");
    TGraph *j2 = new TGraph("CDMS_2011.dat","%lg %lg");

    TGraph *cdms1 = new TGraph("95CLCDMS_2013.dat","%lg %lg");

    TGraph *p1 = new TGraph("COUPP_exp.dat","%lg %lg");
    TGraph *p2 = new TGraph("COUPP_flat.dat","%lg %lg");
    TGraphErrors *p3 = new TGraphErrors("COUPP.dat","%lg %lg %lg");


    TGraph *z1 = new TGraph("LUX_90CL.dat","%lg %lg");
    //TGraph *z2 = new TGraph("LUX_95CL.dat","%lg %lg");


//    TLegend *leg2 = new TLegend(0.65, 0.15, 0.93, 0.42);
    TLegend *leg2 = new TLegend(0.70, 0.15, 0.97, 0.42);
    leg2->AddEntry(g2,"CRESST 1#sigma","F");
    leg2->AddEntry(g1,"CRESST 2#sigma","F");
    leg2->AddEntry(h1,"XENON100(2012)","L");
    leg2->AddEntry(h2,"XENON10(2011)","L");
    leg2->AddEntry(f1,"DAMA/LIBRA","F");
    leg2->AddEntry(k1,"CoGeNT(2013)/90\%CL","F");
    leg2->AddEntry(k2,"CoGeNT(2013)/99\%CL","F");
    //leg2->AddEntry(j1,"CDMS(2010/11)","L");
    leg2->AddEntry(cdms1,"CDMS(2013)/95\%CL","F");
    leg2->AddEntry(p3,"COUPP(2012)","F");
    leg2->AddEntry(z1,"LUX(90\%CL)","L");
    //leg2->AddEntry(z2,"LUX(95\%CL)","L");
//////////////////////////////////////////////////

    h_scalar_min->SetMinimum(1.0e-13);
    h_scalar_min->SetMaximum(1.0e-1);
    h_scalar_min->SetLineColor(3);
    h_scalar_min->SetLineStyle(7);
    h_scalar_min->SetLineWidth(3);
    h_scalar_min->GetXaxis()->SetTitleOffset(1.03);
    h_scalar_min->GetXaxis()->SetTitle("DM Mass #it{M_{#chi}} [GeV]");
    h_scalar_min->GetYaxis()->SetTitle("DM-nucleon cross section #it{#sigma_{#chi-N}^{SI}} [pb]");
    h_scalar_lat->SetLineColor(3);
    h_scalar_lat->SetLineStyle(1);
    h_scalar_lat->SetLineWidth(3);
    h_scalar_max->SetLineColor(3);
    h_scalar_max->SetLineStyle(5);
    h_scalar_max->SetLineWidth(3);


    h_vector_min->SetLineColor(kBlue);
    h_vector_min->SetLineStyle(7);
    h_vector_min->SetLineWidth(3);
    h_vector_lat->SetLineColor(kBlue);
    h_vector_lat->SetLineStyle(1);
    h_vector_lat->SetLineWidth(3);
    h_vector_max->SetLineColor(kBlue);
    h_vector_max->SetLineStyle(5);
    h_vector_max->SetLineWidth(3);


    h_fermion_min->SetLineColor(kRed);
    h_fermion_min->SetLineStyle(7);
    h_fermion_min->SetLineWidth(3);
    h_fermion_lat->SetLineColor(kRed);
    h_fermion_lat->SetLineStyle(1);
    h_fermion_lat->SetLineWidth(3);
    h_fermion_max->SetLineColor(kRed);
    h_fermion_max->SetLineStyle(5);
    h_fermion_max->SetLineWidth(3);


    h_scalar_min->Draw("AC");
    h_scalar_lat->Draw("C");
    h_scalar_max->Draw("C");
    h_vector_min->Draw("C");
    h_vector_lat->Draw("C");
    h_vector_max->Draw("C");
    //h_fermion_min->Draw("C");
    //h_fermion_lat->Draw("C");
    //h_fermion_max->Draw("C");

///////////////////////////////////
    cdms1->SetFillColor(kViolet-9);
    cdms1->Draw("F");

    g1->SetFillColor(kCyan);
    g1->Draw("F");

    g3->SetFillColor(kCyan);
    g3->Draw("F");

    g2->SetFillColor(kPink+1);
    g2->Draw("F");


    h1->SetLineWidth(2);
    h1->SetLineColor(kViolet);
    h1->SetLineStyle(9);
    h1->Draw("L");


    z1->SetLineWidth(2);
    z1->SetLineColor(kRed+4);
    z1->SetLineStyle(9);
    z1->Draw("L");

    //z2->SetLineWidth(2);
    //z2->SetLineColor(kRed);
    //z2->SetLineStyle(9);
    //z2->Draw("L");

    h2->SetLineWidth(2);
    h2->SetLineColor(kOrange+8);
    h2->SetLineStyle(3);
    h2->Draw("L");
   

    f1->SetFillColor(kOrange);
    f1->Draw("F");
    f2->SetFillColor(kOrange);
    f2->Draw("F");


    k2->SetFillColor(kRed-9);
    k2->Draw("F");

    k1->SetFillColor(kRed);
    k1->Draw("F");

    j1->SetLineWidth(2);
    j1->SetLineColor(kOrange+4);
    j1->SetLineStyle(6);
    //j1->Draw("L");

    j2->SetLineWidth(2);
    j2->SetLineColor(kOrange+4);
    j2->SetLineStyle(6);
    //j2->Draw("L");


    p1->SetLineWidth(2);
    p2->SetLineWidth(2);
    p1->Draw("L"); p2->Draw("L");
    p3->SetFillColor(kBlack);
    p3->SetFillStyle(3005);
    p3->DrawClone("E3");


    h_fermion_min->Draw("C");
    h_fermion_lat->Draw("C");
    h_fermion_max->Draw("C");

    leg2->Draw();
///////////////////////////////////

    TPaveText* T = new TPaveText(0.1,0.995,0.9,0.95, "NDC");
    T->SetFillColor(0);
    T->SetFillStyle(0);
    T->SetLineColor(0);
    T->SetTextAlign(22);
    T->SetTextFont(42);
    char Buffer[1024];
    //double iEcm  = systpostfix.Contains('8') ? 8.0 : 7.0;
    //double iLumi = systpostfix.Contains('8') ? 19577 : 5051;
    double iEcm_8  = 8.0;
    double iEcm_7  = 7.0;
    double iLumi_8 = 19712;
    double iLumi_7 = 4900;//5051;


    //T = new TPaveText(0.80,0.92,0.95,0.82, "NDC");
    //sprintf(Buffer, "#splitline{#bf{CMS}}{Preliminary}");
    T = new TPaveText(0.80,0.92,0.95,0.87, "NDC");
    sprintf(Buffer, "#bf{CMS}");
    T->AddText(Buffer);
    T->SetTextFont(42);
    T->Draw("same");
    T->SetBorderSize(0);


    T = new TPaveText(0.35,0.95,0.6,0.85, "NDC");
    sprintf(Buffer, "#splitline{Combination of VBF and}{ZH, H #rightarrow invisible}"); 
    T->AddText(Buffer);
    T->SetTextFont(42);
    T->Draw("same");
    T->SetBorderSize(0);

    T = new TPaveText(0.70,0.82,0.95,0.72, "NDC");
    sprintf(Buffer, "#splitline{#it{B(H#rightarrow inv) < %.2f @%3d%% CL}}{#it{m_{H} = %d GeV}}",BRinv,CLval,Mh);
    T->AddText(Buffer);
    T->SetTextFont(42);
    T->Draw("same");
    T->SetBorderSize(0);

    //T = new TPaveText(0.345,0.89,0.66,0.76, "NDC");
    T = new TPaveText(0.345-0.002,0.89+0.01,0.68+0.005-0.002+0.06,0.73+0.01-0.01, "NDC");
    T = new TPaveText(0.345-0.008,0.86,0.75,0.78, "NDC");
    T->AddText("#sqrt{s} = 8.0 TeV, L = 18.9-19.7 fb^{-1} (VBF+ZH)");
    T->AddText("#sqrt{s} = 7.0 TeV, L = 4.9 fb^{-1} (ZH)");

    //sprintf(Buffer, "#splitline{#sqrt{s} = %.1f TeV, L = 18.9-%.1f fb^{-1} (VBF+ZH)}{#sqrt{s} = %.1f TeV, L = %.1f fb^{-1} (ZH)}", iEcm_8, iLumi_8/1000,iEcm_7, iLumi_7/1000);
    //T->AddText(Buffer);
    T->SetTextAlign(12);
    T->SetTextFont(42);
    T->Draw("same");
    T->SetBorderSize(0);


    T = new TPaveText(0.35,0.265+0.04,0.56,0.215+0.04, "NDC");
    sprintf(Buffer, "#it{vector}");
    T->AddText(Buffer);
    T->SetTextAlign(12);
    T->SetTextColor(kBlue);
    T->SetTextFont(42);
    T->Draw("same");
    T->SetBorderSize(0);

    T = new TPaveText(0.18-0.01,0.58-0.04,0.39-0.01,0.53-0.04, "NDC");
    sprintf(Buffer, "#it{scalar}");
    T->AddText(Buffer);
    T->SetTextAlign(12);
    T->SetTextColor(kGreen);
    T->SetTextFont(42);
    T->Draw("same");
    T->SetBorderSize(0);

    T = new TPaveText(0.2,0.42-0.03,0.41,0.37-0.03, "NDC");
    sprintf(Buffer, "#it{fermion}");
    T->AddText(Buffer);
    T->SetTextAlign(12);
    T->SetTextColor(kRed);
    T->SetTextFont(42);
    T->Draw("same");
    T->SetBorderSize(0);


    TGraph *h_vector_min_leg = new TGraph();
    TGraph *h_vector_lat_leg = new TGraph();
    TGraph *h_vector_max_leg = new TGraph();
    h_vector_min_leg->SetLineStyle(7);
    h_vector_lat_leg->SetLineStyle(1);
    h_vector_max_leg->SetLineStyle(5);
    h_vector_min_leg->SetLineWidth(3);
    h_vector_lat_leg->SetLineWidth(3);
    h_vector_max_leg->SetLineWidth(3);

    TLegend *leg = new TLegend(0.55, 0.2-0.02, 0.83, 0.32-0.02);
    leg->AddEntry(h_vector_min_leg,"Min","L");
    leg->AddEntry(h_vector_lat_leg,"Lattice","L");
    leg->AddEntry(h_vector_max_leg,"Max","L");
    leg->Draw();

    if(Mh == 105) canv->SaveAs("DM_nucleon_XS_ZH105.pdf");
    if(Mh == 115) canv->SaveAs("DM_nucleon_XS_ZH115.pdf");
    if(Mh !=0) { canv->SaveAs("DM_nucleon_XS_ZH125.eps");
		    canv->SaveAs("DM_nucleon_XS_ZH125.pdf");
		    canv->SaveAs("DM_nucleon_XS_ZH125.png");
		    /*canv->SaveAs("DM_nucleon_XS_ZH125.jpg");*/}

    if(Mh == 135) canv->SaveAs("DM_nucleon_XS_ZH135.pdf");
    if(Mh == 145) canv->SaveAs("DM_nucleon_XS_ZH145.pdf");

    //if(Mh == 105) canv->SaveAs("DM_nucleon_XS_ZH105.png");
    //if(Mh == 115) canv->SaveAs("DM_nucleon_XS_ZH115.png");
    //if(Mh == 125) canv->SaveAs("DM_nucleon_XS_ZH125.png");
    //if(Mh == 135) canv->SaveAs("DM_nucleon_XS_ZH135.png");
    //if(Mh == 145) canv->SaveAs("DM_nucleon_XS_ZH145.png");



}
