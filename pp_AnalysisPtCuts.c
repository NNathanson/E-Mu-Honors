
#include <iostream>
#include <cstdio>
#include <tuple>
#include <vector>
#include <stdlib.h>
#include "cmath"
#include <string>

void pp_Analysis2()
{
    gROOT->Macro("loadPythia.C");

////////////Reading from root file
    TFile *f = new TFile("ntuples10Mill_MoreInfo.root");

    TNtuple *pairsClean = (TNtuple*)f->Get("ntAll_clean");
    TNtuple *pairsDirty = (TNtuple*)f->Get("ntAll_dirty");
    TNtuple *pairsDetector = (TNtuple*)f->Get("ntAll_detector");
    TNtuple *elecsNoBkg = (TNtuple*)f->Get("noisyElecs");
    TNtuple *elecsBkg = (TNtuple*)f->Get("cleanElecs");
    TNtuple *ntMuons = (TNtuple*)f->Get("ntMuons");

    //Defining variables as described in the paper
    Double_t Ntrig = elecsNoBkg->GetEntries("ptE>3");
    Double_t Nmu = pairsClean->GetEntries("ptE>3 && ptMu>3");
    Double_t sigma = 3.04e-3;



////////////Drawing graphs

// 2: Fit total dPhi graphs
    TF1 *fSpectrum = new TF1("fSpectrum","gaus+pol0(3)",0.4,5.6);
    Int_t bins = 40;
   /* Double_t width = 5.2/(Double_t)bins;
    

    //Fit all potential pairs ptE>1 for visual reference
    TCanvas* fitsAll = new TCanvas("fitsAll","dPhi Fits");

    fitsAll->cd(1);
    TH1F* pcPhiAll = new TH1F("pcPhiAll", "clean dPhi, ptE>3", bins, 0, 2*M_PI);
    pcPhiAll->SetLineColor(kBlue);
    pairsClean->Draw("dPhi>>pcPhiAll","ptE>3 && ptMu>3"); 

    pcPhiAll->Scale(1/(width));
    pcPhiAll->Scale(1/Ntrig);

    fSpectrum->SetParameters(0.03,3,0.6,0.006);
    pcPhiAll->Fit("fSpectrum","","",0.4,5.6);

    pcPhiAll->SetXTitle("#Delta#phi");
    pcPhiAll->SetTitleSize(0.04,"X");
    pcPhiAll->SetYTitle("(1/N_{e})dN_{#mu}/d#Delta #phi");
    pcPhiAll->SetTitleSize(0.04,"Y");
    fitsAll->Modified(); 
    fitsAll->Update();

    Double_t fwhmC = fSpectrum->GetParameter(2)*2.35482;
    Double_t fwhmCErr = fSpectrum->GetParError(2)*2.35482;

    cout<<"Width (confirmed quark parent): "<<fwhmC<<" +/- "<<fwhmCErr<<endl;*/

//
    TCanvas* cuts = new TCanvas("cuts","Checking different cuts");

    TGraphErrors* cutsA = new TGraphErrors(3);
    TGraphErrors* cutsB = new TGraphErrors(3);
    TGraphErrors* cutsC = new TGraphErrors(3);
    TGraphErrors* cutsD = new TGraphErrors(3);
    TMultiGraph *mg = new TMultiGraph();

  //cut 1: ptE>1
    //Double_t    = 1/(2*M_PI*1)
    Double_t scale1 = elecsNoBkg->GetEntries("ptE>1"); //check this

    TH1F* cut1a = new TH1F("cut1a", "ptE>1 ptMu>0.3", bins, 0.4,5.6);
    pairsClean->Draw("dPhi>>cut1a","ptE>1 && ptMu>0.3");

    fSpectrum->SetParameters(0.03,3,0.6,0.006);
    cut1a->Scale(1/scale1,"width");
    TFitResultPtr r = cut1a->Fit("fSpectrum","S","",0.4,5.6);

    //cutsA->SetPoint(0, 1, fSpectrum->GetParameter(2)*2.35482);
    ////cutsA->SetPointError(0, 0, fSpectrum->GetParError(2)*2.35482);
    Double_t bkgInteg = fSpectrum->GetParameter(3)*5.2;
    cutsA->SetPoint(0, 1, fSpectrum->Integral(0.4,5.6)-bkgInteg);
    cutsA->SetPointError(0, 0, fSpectrum->IntegralError(0.4,5.6,r->GetParams(), r->GetCovarianceMatrix().GetMatrixArray()));



    TH1F* cut1b = new TH1F("cut1b", "ptE>1 ptMu>1", bins, 0.4,5.6);
    pairsClean->Draw("dPhi>>cut1b","ptE>1 && ptMu>1");

    fSpectrum->SetParameters(0.03,3,0.6,0.006);
    cut1b->Scale(1/scale1,"width");
    r = cut1b->Fit("fSpectrum","S","",0.4,5.6);

    //cutsB->SetPoint(0, 1, fSpectrum->GetParameter(2)*2.35482);
    ////cutsB->SetPointError(0, 0, fSpectrum->GetParError(2)*2.35482);
    bkgInteg = fSpectrum->GetParameter(3)*5.2;
    cutsB->SetPoint(0, 1, fSpectrum->Integral(0.4,5.6)-bkgInteg);
    cutsB->SetPointError(0, 0, fSpectrum->IntegralError(0.4,5.6,r->GetParams(), r->GetCovarianceMatrix().GetMatrixArray()));


    TH1F* cut1c = new TH1F("cut1c", "ptE>1 ptMu>3", bins, 0.4,5.6);
    pairsClean->Draw("dPhi>>cut1c","ptE>1 && ptMu>3");

    fSpectrum->SetParameters(0.03,3,0.6,0.006);
    cut1c->Scale(1/scale1,"width");
    r = cut1c->Fit("fSpectrum","S","",0.4,5.6);

    //cutsC->SetPoint(0, 1, fSpectrum->GetParameter(2)*2.35482);
    ////cutsC->SetPointError(0, 0, fSpectrum->GetParError(2)*2.35482);
    bkgInteg = fSpectrum->GetParameter(3)*5.2;
    cutsC->SetPoint(0, 1, fSpectrum->Integral(0.4,5.6)-bkgInteg);
    cutsC->SetPointError(0, 0, fSpectrum->IntegralError(0.4,5.6,r->GetParams(), r->GetCovarianceMatrix().GetMatrixArray()));


    TH1F* cut1d = new TH1F("cut1d", "ptE>1 ptMu>5", bins, 0.4,5.6);
    pairsClean->Draw("dPhi>>cut1d","ptE>1 && ptMu>5");

    fSpectrum->SetParameters(0.03,3,0.6,0.006);
    cut1d->Scale(1/scale1,"width");
    r = cut1d->Fit("fSpectrum","S","",0.4,5.6);

    //cutsD->SetPoint(0, 1, fSpectrum->GetParameter(2)*2.35482);
    ////cutsD->SetPointError(0, 0, fSpectrum->GetParError(2)*2.35482);
    bkgInteg = fSpectrum->GetParameter(3)*5.2;
    cutsD->SetPoint(0, 1, fSpectrum->Integral(0.4,5.6)-bkgInteg);
    cutsD->SetPointError(0, 0, fSpectrum->IntegralError(0.4,5.6,r->GetParams(), r->GetCovarianceMatrix().GetMatrixArray()));

  //cut 2: ptE>3
    Double_t scale2 = elecsNoBkg->GetEntries("ptE>3"); //check this

    TH1F* cut2a = new TH1F("cut2a", "ptE>3 ptMu>0.3", bins, 0.4,5.6);
    pairsClean->Draw("dPhi>>cut2a","ptE>3 && ptMu>0.3");

    fSpectrum->SetParameters(0.03,3,0.6,0.006);
    cut2a->Scale(1/scale2,"width");
    r = cut2a->Fit("fSpectrum","S","",0.4,5.6);

    //cutsA->SetPoint(1, 3, fSpectrum->GetParameter(2)*2.35482);
    ////cutsA->SetPointError(1, 0, fSpectrum->GetParError(2)*2.35482);
    bkgInteg = fSpectrum->GetParameter(3)*5.2;
    cutsA->SetPoint(1, 3, fSpectrum->Integral(0.4,5.6)-bkgInteg);
    cutsA->SetPointError(1, 0, fSpectrum->IntegralError(0.4,5.6,r->GetParams(), r->GetCovarianceMatrix().GetMatrixArray()));


    TH1F* cut2b = new TH1F("cut2b", "ptE>3 ptMu>1", bins, 0.4,5.6);
    pairsClean->Draw("dPhi>>cut2b","ptE>3 && ptMu>1");

    fSpectrum->SetParameters(0.03,3,0.6,0.006);
    cut2b->Scale(1/scale2,"width");
    r = cut2b->Fit("fSpectrum","S","",0.4,5.6);

    //cutsB->SetPoint(1, 3, fSpectrum->GetParameter(2)*2.35482);
    ////cutsB->SetPointError(1, 0, fSpectrum->GetParError(2)*2.35482);
    bkgInteg = fSpectrum->GetParameter(3)*5.2;
    cutsB->SetPoint(1, 3, fSpectrum->Integral(0.4,5.6)-bkgInteg);
    cutsB->SetPointError(1, 0,fSpectrum->IntegralError(0.4,5.6,r->GetParams(), r->GetCovarianceMatrix().GetMatrixArray()));



    TH1F* cut2c = new TH1F("cut2c", "ptE>3 ptMu>3", bins, 0.4,5.6);
    pairsClean->Draw("dPhi>>cut2c","ptE>3 && ptMu>3");

    fSpectrum->SetParameters(0.03,3,0.6,0.006);
    cut2c->Scale(1/scale2,"width");
    r = cut2c->Fit("fSpectrum","S","",0.4,5.6);

    //cutsC->SetPoint(1, 3, fSpectrum->GetParameter(2)*2.35482);
    ////cutsC->SetPointError(1, 0, fSpectrum->GetParError(2)*2.35482);
    bkgInteg = fSpectrum->GetParameter(3)*5.2;
    cutsC->SetPoint(1, 3, fSpectrum->Integral(0.4,5.6)-bkgInteg);
    cutsC->SetPointError(1, 0, fSpectrum->IntegralError(0.4,5.6,r->GetParams(), r->GetCovarianceMatrix().GetMatrixArray()));


    TH1F* cut2d = new TH1F("cut2d", "ptE>3 ptMu>5", bins, 0.4,5.6);
    pairsClean->Draw("dPhi>>cut2d","ptE>3 && ptMu>5");

    fSpectrum->SetParameters(0.03,3,0.6,0.006);
    cut2d->Scale(1/scale2,"width");
    r = cut2d->Fit("fSpectrum","S","",0.4,5.6);

    //cutsD->SetPoint(1, 3, fSpectrum->GetParameter(2)*2.35482);
    ////cutsD->SetPointError(1, 0, fSpectrum->GetParError(2)*2.35482);
    bkgInteg = fSpectrum->GetParameter(3)*5.2;
    cutsD->SetPoint(1, 3, fSpectrum->Integral(0.4,5.6)-bkgInteg);
    cutsD->SetPointError(1, 0, fSpectrum->IntegralError(0.4,5.6,r->GetParams(), r->GetCovarianceMatrix().GetMatrixArray()));

  //cut 3: ptE>5
    Double_t scale3 = elecsNoBkg->GetEntries("ptE>5"); //check this

    TH1F* cut3a = new TH1F("cut3a", "ptE>5 ptMu>0.3", bins, 0.4,5.6);
    pairsClean->Draw("dPhi>>cut3a","ptE>5 && ptMu>0.3");

    fSpectrum->SetParameters(0.03,3,0.6,0.006);
    cut3a->Scale(1/scale3,"width");
    r = cut3a->Fit("fSpectrum","S","",0.4,5.6);

    //cutsA->SetPoint(2, 5, fSpectrum->GetParameter(2)*2.35482);
    ////cutsA->SetPointError(2, 0, fSpectrum->GetParError(2)*2.35482);
    bkgInteg = fSpectrum->GetParameter(3)*5.2;
    cutsA->SetPoint(2, 5, fSpectrum->Integral(0.4,5.6)-bkgInteg);
    cutsA->SetPointError(2, 0, fSpectrum->IntegralError(0.4,5.6,r->GetParams(), r->GetCovarianceMatrix().GetMatrixArray()));



    TH1F* cut3b = new TH1F("cut3b", "ptE>5 ptMu>1", bins, 0.4,5.6);
    pairsClean->Draw("dPhi>>cut3b","ptE>5 && ptMu>1");

    fSpectrum->SetParameters(0.03,3,0.6,0.006);
    cut3b->Scale(1/scale3,"width");
    r = cut3b->Fit("fSpectrum","S","",0.4,5.6);

    //cutsB->SetPoint(2, 5, fSpectrum->GetParameter(2)*2.35482);
    ////cutsB->SetPointError(2, 0, fSpectrum->GetParError(2)*2.35482);
    bkgInteg = fSpectrum->GetParameter(3)*5.2;
    cutsB->SetPoint(2, 5, fSpectrum->Integral(0.4,5.6)-bkgInteg);
    cutsB->SetPointError(2, 0, fSpectrum->IntegralError(0.4,5.6,r->GetParams(), r->GetCovarianceMatrix().GetMatrixArray()));



    TH1F* cut3c = new TH1F("cut3c", "ptE>5 ptMu>3", bins, 0.4,5.6);
    pairsClean->Draw("dPhi>>cut3c","ptE>5 && ptMu>3");

    fSpectrum->SetParameters(0.03,3,0.6,0.006);
    cut3c->Scale(1/scale3,"width");
    r = cut3c->Fit("fSpectrum","S","",0.4,5.6);

    //cutsC->SetPoint(2, 5, fSpectrum->GetParameter(2)*2.35482);
    ////cutsC->SetPointError(2, 0, fSpectrum->GetParError(2)*2.35482);
    bkgInteg = fSpectrum->GetParameter(3)*5.2;
    cutsC->SetPoint(2, 5, fSpectrum->Integral(0.4,5.6)-bkgInteg);
    cutsC->SetPointError(2, 0, fSpectrum->IntegralError(0.4,5.6,r->GetParams(), r->GetCovarianceMatrix().GetMatrixArray()));


    TH1F* cut3d = new TH1F("cut3d", "ptE>5 ptMu>5", bins, 0.4,5.6);
    pairsClean->Draw("dPhi>>cut3d","ptE>5 && ptMu>5");

    fSpectrum->SetParameters(0.03,3,0.6,0.006);
    cut3d->Scale(1/scale3,"width");
    r = cut3d->Fit("fSpectrum","S","",0.4,5.6);

    //cutsD->SetPoint(2, 5, fSpectrum->GetParameter(2)*2.35482);
    ////cutsD->SetPointError(2, 0, fSpectrum->GetParError(2)*2.35482);
    bkgInteg = fSpectrum->GetParameter(3)*5.2;
    cutsD->SetPoint(2, 5, fSpectrum->Integral(0.4,5.6)-bkgInteg);
    cutsD->SetPointError(2, 0, fSpectrum->IntegralError(0.4,5.6,r->GetParams(), r->GetCovarianceMatrix().GetMatrixArray()));
    cout<<"!!!!!"<< fSpectrum->IntegralError(0.4,5.6,r->GetParams(), r->GetCovarianceMatrix().GetMatrixArray())<<endl;


  //

    TCanvas* cutsMG = new TCanvas("cutsMG","Multigraph of ptE vs ptMu cut");
    cutsMG->DrawFrame(0,0,6,0.06);
    gPad->SetLogy();
    cutsA->SetMarkerColor(kRed);
    cutsA->SetLineColor(kRed);
    cutsA->SetMarkerStyle(20);

    cutsB->SetMarkerColor(kViolet);
    cutsB->SetLineColor(kViolet);
    cutsB->SetMarkerStyle(21);

    cutsC->SetMarkerColor(kBlue);
    cutsC->SetLineColor(kBlue);
    cutsC->SetMarkerStyle(22);

    cutsD->SetMarkerColor(8);
    cutsD->SetLineColor(8);
    cutsD->SetMarkerStyle(23);

    mg->Add(cutsA);
    mg->Add(cutsB);
    mg->Add(cutsC);
    mg->Add(cutsD);

    
    mg->Draw("APC");
    mg->GetXaxis()->SetRangeUser(0,6);
    mg->GetXaxis()->SetLimits(0,6);
    mg->GetXaxis()->SetTitle("p_{t}^{e} (GeV)");
    mg->GetYaxis()->SetTitle("(1/N_{trig})dN/d#Delta #phi (Integral)");
    mg->GetXaxis()->SetTitleSize(0.03);
    mg->GetYaxis()->SetTitleSize(0.03);
    cutsMG->Modified(); 
    cutsMG->Update();

    TLegend* key = new TLegend(1,1,1,1);
    key->AddEntry(cutsA, "p_{t}(#mu)>0.3", "P");
    key->AddEntry(cutsB, "p_{t}(#mu)>1", "P");
    key->AddEntry(cutsC, "p_{t}(#mu)>3", "P");
    key->AddEntry(cutsD, "p_{t}(#mu)>5", "P");
    key->Draw();


}