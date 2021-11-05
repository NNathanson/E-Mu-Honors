
#include <iostream>
#include <cstdio>
#include <tuple>
#include <vector>
#include <stdlib.h>
#include "cmath"
#include <string>

void pp_Analysis_CS()
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
    Double_t sigma = 3.04e-3;
    Double_t xArr[18] = {1.25,1.75,2.25,2.75,3.25,3.75,4.25,4.75,5.25,5.75,6.25,6.75,7.25,7.75,8.25,8.75,9.25,9.75};
    Double_t xErr[18] = {0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5};


////// a: Charm Quark Parent
    TCanvas* rAll = new TCanvas("rAll","dPhi Fit Ranges (charm quark parent)");
    cout<<endl<<endl<<"CHARM QUARK PARENT"<<endl<<endl;
    Double_t bins = 85;
    TF1 *fSpectrum = new TF1("fSpectrum","gaus+pol0(3)",0.4,5.6);
/*
////Muons
    //define and fill histogram objects to use in the fits
    TH1F* pcPhi1 = new TH1F("pcPhi1", "clean dPhi, 1<ptMu<1.5", bins, 0.4,5.6);
    TH1F* pcPhi2 = new TH1F("pcPhi2", "clean dPhi, 1.5<ptMu<2", bins, 0.4,5.6);
    TH1F* pcPhi3 = new TH1F("pcPhi3", "clean dPhi, 2<ptMu<2.5", bins, 0.4,5.6);
    TH1F* pcPhi4 = new TH1F("pcPhi4", "clean dPhi, 2.5<ptMu<3", bins, 0.4,5.6);
    TH1F* pcPhi5 = new TH1F("pcPhi5", "clean dPhi, 3<ptMu<3.5", bins, 0.4,5.6);
    TH1F* pcPhi6 = new TH1F("pcPhi6", "clean dPhi, 3.5<ptMu<4", bins, 0.4,5.6);
    TH1F* pcPhi7 = new TH1F("pcPhi7", "clean dPhi, 4<ptMu<4.5", bins, 0.4,5.6);
    TH1F* pcPhi8 = new TH1F("pcPhi8", "clean dPhi, 4.5<ptMu<5", bins, 0.4,5.6);
    TH1F* pcPhi9 = new TH1F("pcPhi9", "clean dPhi, 5<ptMu<5.5", bins, 0.4,5.6);
    TH1F* pcPhi10 = new TH1F("pcPhi10", "clean dPhi, 5.5<ptMu<6", bins, 0.4,5.6);
    TH1F* pcPhi11 = new TH1F("pcPhi11", "clean dPhi, 6<ptMu<6.5", bins, 0.4,5.6);
    TH1F* pcPhi12 = new TH1F("pcPhi12", "clean dPhi, 6.5<ptMu<7", bins, 0.4,5.6);
    TH1F* pcPhi13 = new TH1F("pcPhi13", "clean dPhi, 7<ptMu<7.5", bins, 0.4,5.6);
    TH1F* pcPhi14 = new TH1F("pcPhi14", "clean dPhi, 7.5<ptMu<8", bins, 0.4,5.6);
    TH1F* pcPhi15 = new TH1F("pcPhi15", "clean dPhi, 8<ptMu<8.5", bins, 0.4,5.6);
    TH1F* pcPhi16 = new TH1F("pcPhi16", "clean dPhi, 8.5<ptMu<9", bins, 0.4,5.6);
    TH1F* pcPhi17 = new TH1F("pcPhi17", "clean dPhi, 9<ptMu<9.5", bins, 0.4,5.6);
    TH1F* pcPhi18 = new TH1F("pcPhi18", "clean dPhi, 9.5<ptMu<10", bins, 0.4,5.6);

    pairsClean->Draw("dPhi>>pcPhi1","ptE>3 && ptMu>1 && ptMu<1.5 && 0.4<dPhi && 5.6>dPhi");
    pairsClean->Draw("dPhi>>pcPhi2","ptE>3 && ptMu>1.5 && ptMu<2 && 0.4<dPhi && 5.6>dPhi");
    pairsClean->Draw("dPhi>>pcPhi3","ptE>3 && ptMu>2 && ptMu<2.5 && 0.4<dPhi && 5.6>dPhi");
    pairsClean->Draw("dPhi>>pcPhi4","ptE>3 && ptMu>2.5 && ptMu<3 && 0.4<dPhi && 5.6>dPhi");
    pairsClean->Draw("dPhi>>pcPhi5","ptE>3 && ptMu>3 && ptMu<3.5 && 0.4<dPhi && 5.6>dPhi");
    pairsClean->Draw("dPhi>>pcPhi6","ptE>3 && ptMu>3.5 && ptMu<4 && 0.4<dPhi && 5.6>dPhi");
    pairsClean->Draw("dPhi>>pcPhi7","ptE>3 && ptMu>4 && ptMu<4.5 && 0.4<dPhi && 5.6>dPhi");
    pairsClean->Draw("dPhi>>pcPhi8","ptE>3 && ptMu>4.5 && ptMu<5 && 0.4<dPhi && 5.6>dPhi");
    pairsClean->Draw("dPhi>>pcPhi9","ptE>3 && ptMu>5 && ptMu<5.5 && 0.4<dPhi && 5.6>dPhi");
    pairsClean->Draw("dPhi>>pcPhi10","ptE>3 && ptMu>5.5 && ptMu<6 && 0.4<dPhi && 5.6>dPhi");
    pairsClean->Draw("dPhi>>pcPhi11","ptE>3 && ptMu>6 && ptMu<6.5 && 0.4<dPhi && 5.6>dPhi");
    pairsClean->Draw("dPhi>>pcPhi12","ptE>3 && ptMu>6.5 && ptMu<7 && 0.4<dPhi && 5.6>dPhi");
    pairsClean->Draw("dPhi>>pcPhi13","ptE>3 && ptMu>7 && ptMu<7.5 && 0.4<dPhi && 5.6>dPhi");
    pairsClean->Draw("dPhi>>pcPhi14","ptE>3 && ptMu>7.5 && ptMu<8 && 0.4<dPhi && 5.6>dPhi");
    pairsClean->Draw("dPhi>>pcPhi15","ptE>3 && ptMu>8 && ptMu<8.5 && 0.4<dPhi && 5.6>dPhi");
    pairsClean->Draw("dPhi>>pcPhi16","ptE>3 && ptMu>8.5 && ptMu<9 && 0.4<dPhi && 5.6>dPhi");
    pairsClean->Draw("dPhi>>pcPhi17","ptE>3 && ptMu>9 && ptMu<9.5 && 0.4<dPhi && 5.6>dPhi");
    pairsClean->Draw("dPhi>>pcPhi18","ptE>3 && ptMu>9.5 && ptMu<10 && 0.4<dPhi && 5.6>dPhi");
*/
////Electrons
    TH1F* pcPhi1 = new TH1F("pcPhi1", "clean dPhi, 1<pTE<1.5", bins, 0.4,5.6);
    TH1F* pcPhi2 = new TH1F("pcPhi2", "clean dPhi, 1.5<pTE<2", bins, 0.4,5.6);
    TH1F* pcPhi3 = new TH1F("pcPhi3", "clean dPhi, 2<pTE<2.5", bins, 0.4,5.6);
    TH1F* pcPhi4 = new TH1F("pcPhi4", "clean dPhi, 2.5<pTE<3", bins, 0.4,5.6);
    TH1F* pcPhi5 = new TH1F("pcPhi5", "clean dPhi, 3<pTE<3.5", bins, 0.4,5.6);
    TH1F* pcPhi6 = new TH1F("pcPhi6", "clean dPhi, 3.5<pTE<4", bins, 0.4,5.6);
    TH1F* pcPhi7 = new TH1F("pcPhi7", "clean dPhi, 4<pTE<4.5", bins, 0.4,5.6);
    TH1F* pcPhi8 = new TH1F("pcPhi8", "clean dPhi, 4.5<pTE<5", bins, 0.4,5.6);
    TH1F* pcPhi9 = new TH1F("pcPhi9", "clean dPhi, 5<pTE<5.5", bins, 0.4,5.6);
    TH1F* pcPhi10 = new TH1F("pcPhi10", "clean dPhi, 5.5<pTE<6", bins, 0.4,5.6);
    TH1F* pcPhi11 = new TH1F("pcPhi11", "clean dPhi, 6<pTE<6.5", bins, 0.4,5.6);
    TH1F* pcPhi12 = new TH1F("pcPhi12", "clean dPhi, 6.5<pTE<7", bins, 0.4,5.6);
    TH1F* pcPhi13 = new TH1F("pcPhi13", "clean dPhi, 7<pTE<7.5", bins, 0.4,5.6);
    TH1F* pcPhi14 = new TH1F("pcPhi14", "clean dPhi, 7.5<pTE<8", bins, 0.4,5.6);
    TH1F* pcPhi15 = new TH1F("pcPhi15", "clean dPhi, 8<pTE<8.5", bins, 0.4,5.6);
    TH1F* pcPhi16 = new TH1F("pcPhi16", "clean dPhi, 8.5<pTE<9", bins, 0.4,5.6);
    TH1F* pcPhi17 = new TH1F("pcPhi17", "clean dPhi, 9<pTE<9.5", bins, 0.4,5.6);
    TH1F* pcPhi18 = new TH1F("pcPhi18", "clean dPhi, 9.5<pTE<10", bins, 0.4,5.6);

    pairsClean->Draw("dPhi>>pcPhi1","ptMu>3 && ptE>1 && ptE<1.5 && 0.4<dPhi && 5.6>dPhi");
    pairsClean->Draw("dPhi>>pcPhi2","ptMu>3 && ptE>1.5 && ptE<2 && 0.4<dPhi && 5.6>dPhi");
    pairsClean->Draw("dPhi>>pcPhi3","ptMu>3 && ptE>2 && ptE<2.5 && 0.4<dPhi && 5.6>dPhi");
    pairsClean->Draw("dPhi>>pcPhi4","ptMu>3 && ptE>2.5 && ptE<3 && 0.4<dPhi && 5.6>dPhi");
    pairsClean->Draw("dPhi>>pcPhi5","ptMu>3 && ptE>3 && ptE<3.5 && 0.4<dPhi && 5.6>dPhi");
    pairsClean->Draw("dPhi>>pcPhi6","ptMu>3 && ptE>3.5 && ptE<4 && 0.4<dPhi && 5.6>dPhi");
    pairsClean->Draw("dPhi>>pcPhi7","ptMu>3 && ptE>4 && ptE<4.5 && 0.4<dPhi && 5.6>dPhi");
    pairsClean->Draw("dPhi>>pcPhi8","ptMu>3 && ptE>4.5 && ptE<5 && 0.4<dPhi && 5.6>dPhi");
    pairsClean->Draw("dPhi>>pcPhi9","ptMu>3 && ptE>5 && ptE<5.5 && 0.4<dPhi && 5.6>dPhi");
    pairsClean->Draw("dPhi>>pcPhi10","ptMu>3 && ptE>5.5 && ptE<6 && 0.4<dPhi && 5.6>dPhi");
    pairsClean->Draw("dPhi>>pcPhi11","ptMu>3 && ptE>6 && ptE<6.5 && 0.4<dPhi && 5.6>dPhi");
    pairsClean->Draw("dPhi>>pcPhi12","ptMu>3 && ptE>6.5 && ptE<7 && 0.4<dPhi && 5.6>dPhi");
    pairsClean->Draw("dPhi>>pcPhi13","ptMu>3 && ptE>7 && ptE<7.5 && 0.4<dPhi && 5.6>dPhi");
    pairsClean->Draw("dPhi>>pcPhi14","ptMu>3 && ptE>7.5 && ptE<8 && 0.4<dPhi && 5.6>dPhi");
    pairsClean->Draw("dPhi>>pcPhi15","ptMu>3 && ptE>8 && ptE<8.5 && 0.4<dPhi && 5.6>dPhi");
    pairsClean->Draw("dPhi>>pcPhi16","ptMu>3 && ptE>8.5 && ptE<9 && 0.4<dPhi && 5.6>dPhi");
    pairsClean->Draw("dPhi>>pcPhi17","ptMu>3 && ptE>9 && ptE<9.5 && 0.4<dPhi && 5.6>dPhi");
    pairsClean->Draw("dPhi>>pcPhi18","ptMu>3 && ptE>9.5 && ptE<10 && 0.4<dPhi && 5.6>dPhi");

    std::vector<Double_t> widthsC;
    std::vector<Double_t> widthsCErr;
    std::vector<Double_t> integC;
    std::vector<Double_t> integCErr;
    Double_t bkgInteg = 0;
    int count = 0;

    TH1F* allCleanHists[18] = {pcPhi1,pcPhi2,pcPhi3,pcPhi4,pcPhi5,pcPhi6,pcPhi7,pcPhi8,pcPhi9,pcPhi10,pcPhi11,pcPhi12,pcPhi13,pcPhi14,pcPhi15,pcPhi16,pcPhi17,pcPhi18};

    for (TH1F* hist: allCleanHists) 
    {
        cout<<endl<<count+1<<endl;
    
        hist->Draw();
        fSpectrum->SetParameters(10,3,0.6,11);
        TFitResultPtr r = hist->Fit("fSpectrum","S","same",0.4,5.6);

        widthsC.push_back(std::abs(fSpectrum->GetParameter(2)*2.35482));
        widthsCErr.push_back(fSpectrum->GetParError(2)*2.35482);

        fSpectrum->SetParameter(2,std::abs(fSpectrum->GetParameter(2)));
        bkgInteg = fSpectrum->GetParameter(3)*5.2;

        Double_t integ = fSpectrum->Integral(0.4,5.6)-bkgInteg;
        Double_t fac = (1/(2*M_PI*xArr[count])) * (sigma/0.5);

        integC.push_back(fac*integ);
        integCErr.push_back(fac*fSpectrum->IntegralError(0.4,5.6,r->GetParams(), r->GetCovarianceMatrix().GetMatrixArray()));
        count++;
    }

    Double_t wcArr[widthsC.size()];
    std::copy(widthsC.begin(), widthsC.end(), wcArr);
    Double_t wcErrArr[widthsCErr.size()];
    std::copy(widthsCErr.begin(), widthsCErr.end(), wcErrArr);

    TGraphErrors* widthCompClean = new TGraphErrors(18,xArr,wcArr,0,wcErrArr);

    Double_t icArr[integC.size()];
    std::copy(integC.begin(), integC.end(), icArr);
    Double_t icErrArr[integCErr.size()];
    std::copy(integCErr.begin(), integCErr.end(), icErrArr);

    TGraphErrors* integCompClean = new TGraphErrors(18,xArr,icArr,0,icErrArr);
    cout<<endl<<"Quark Parent"<<endl;
    for(int i=0; i<18; i++){cout<<widthsC[i]<<" +/- "<<widthsCErr[i]<<"  \t\t"<<integC[i]<<" +/- "<<integCErr[i]<<endl;}
    
////// b: No Background
    //TCanvas* rD = new TCanvas("rD","dPhi Fit Ranges (no background)");
    cout<<endl<<endl<<"NO BACKGROUND"<<endl<<endl;
/*
////Muons
    //define and fill histogram objects to use in the fits
    TH1F* pdPhi1 = new TH1F("pdPhi1", "no bkg dPhi, 1<ptMu<1.5", bins, 0.4,5.6);
    TH1F* pdPhi2 = new TH1F("pdPhi2", "no bkg dPhi, 1.5<ptMu<2", bins, 0.4,5.6);
    TH1F* pdPhi3 = new TH1F("pdPhi3", "no bkg dPhi, 2<ptMu<2.5", bins, 0.4,5.6);
    TH1F* pdPhi4 = new TH1F("pdPhi4", "no bkg dPhi, 2.5<ptMu<3", bins, 0.4,5.6);
    TH1F* pdPhi5 = new TH1F("pdPhi5", "no bkg dPhi, 3<ptMu<3.5", bins, 0.4,5.6);
    TH1F* pdPhi6 = new TH1F("pdPhi6", "no bkg dPhi, 3.5<ptMu<4", bins, 0.4,5.6);
    TH1F* pdPhi7 = new TH1F("pdPhi7", "no bkg dPhi, 4<ptMu<4.5", bins, 0.4,5.6);
    TH1F* pdPhi8 = new TH1F("pdPhi8", "no bkg dPhi, 4.5<ptMu<5", bins, 0.4,5.6);
    TH1F* pdPhi9 = new TH1F("pdPhi9", "no bkg dPhi, 5<ptMu<5.5", bins, 0.4,5.6);
    TH1F* pdPhi10 = new TH1F("pdPhi10", "no bkg dPhi, 5.5<ptMu<6", bins, 0.4,5.6);
    TH1F* pdPhi11 = new TH1F("pdPhi11", "no bkg dPhi, 6<ptMu<6.5", bins, 0.4,5.6);
    TH1F* pdPhi12 = new TH1F("pdPhi12", "no bkg dPhi, 6.5<ptMu<7", bins, 0.4,5.6);
    TH1F* pdPhi13 = new TH1F("pdPhi13", "no bkg dPhi, 7<ptMu<7.5", bins, 0.4,5.6);
    TH1F* pdPhi14 = new TH1F("pdPhi14", "no bkg dPhi, 7.5<ptMu<8", bins, 0.4,5.6);
    TH1F* pdPhi15 = new TH1F("pdPhi15", "no bkg dPhi, 8<ptMu<8.5", bins, 0.4,5.6);
    TH1F* pdPhi16 = new TH1F("pdPhi16", "no bkg dPhi, 8.5<ptMu<9", bins, 0.4,5.6);
    TH1F* pdPhi17 = new TH1F("pdPhi17", "no bkg dPhi, 9<ptMu<9.5", bins, 0.4,5.6);
    TH1F* pdPhi18 = new TH1F("pdPhi18", "no bkg dPhi, 9.5<ptMu<10", bins, 0.4,5.6);

    pairsDirty->Draw("dPhi>>pdPhi1","ptE>3 && ptMu>1 && ptMu<1.5 && 0.4<dPhi && 5.6>dPhi");
    pairsDirty->Draw("dPhi>>pdPhi2","ptE>3 && ptMu>1.5 && ptMu<2 && 0.4<dPhi && 5.6>dPhi");
    pairsDirty->Draw("dPhi>>pdPhi3","ptE>3 && ptMu>2 && ptMu<2.5 && 0.4<dPhi && 5.6>dPhi");
    pairsDirty->Draw("dPhi>>pdPhi4","ptE>3 && ptMu>2.5 && ptMu<3 && 0.4<dPhi && 5.6>dPhi");
    pairsDirty->Draw("dPhi>>pdPhi5","ptE>3 && ptMu>3 && ptMu<3.5 && 0.4<dPhi && 5.6>dPhi");
    pairsDirty->Draw("dPhi>>pdPhi6","ptE>3 && ptMu>3.5 && ptMu<4 && 0.4<dPhi && 5.6>dPhi");
    pairsDirty->Draw("dPhi>>pdPhi7","ptE>3 && ptMu>4 && ptMu<4.5 && 0.4<dPhi && 5.6>dPhi");
    pairsDirty->Draw("dPhi>>pdPhi8","ptE>3 && ptMu>4.5 && ptMu<5 && 0.4<dPhi && 5.6>dPhi");
    pairsDirty->Draw("dPhi>>pdPhi9","ptE>3 && ptMu>5 && ptMu<5.5 && 0.4<dPhi && 5.6>dPhi");
    pairsDirty->Draw("dPhi>>pdPhi10","ptE>3 && ptMu>5.5 && ptMu<6 && 0.4<dPhi && 5.6>dPhi");
    pairsDirty->Draw("dPhi>>pdPhi11","ptE>3 && ptMu>6 && ptMu<6.5 && 0.4<dPhi && 5.6>dPhi");
    pairsDirty->Draw("dPhi>>pdPhi12","ptE>3 && ptMu>6.5 && ptMu<7 && 0.4<dPhi && 5.6>dPhi");
    pairsDirty->Draw("dPhi>>pdPhi13","ptE>3 && ptMu>7 && ptMu<7.5 && 0.4<dPhi && 5.6>dPhi");
    pairsDirty->Draw("dPhi>>pdPhi14","ptE>3 && ptMu>7.5 && ptMu<8 && 0.4<dPhi && 5.6>dPhi");
    pairsDirty->Draw("dPhi>>pdPhi15","ptE>3 && ptMu>8 && ptMu<8.5 && 0.4<dPhi && 5.6>dPhi");
    pairsDirty->Draw("dPhi>>pdPhi16","ptE>3 && ptMu>8.5 && ptMu<9 && 0.4<dPhi && 5.6>dPhi");
    pairsDirty->Draw("dPhi>>pdPhi17","ptE>3 && ptMu>9 && ptMu<9.5 && 0.4<dPhi && 5.6>dPhi");
    pairsDirty->Draw("dPhi>>pdPhi18","ptE>3 && ptMu>9.5 && ptMu<10 && 0.4<dPhi && 5.6>dPhi");
*/
////Electrons
    TH1F* pdPhi1 = new TH1F("pdPhi1", "no bkg dPhi, 1<pTE<1.5", bins, 0.4,5.6);
    TH1F* pdPhi2 = new TH1F("pdPhi2", "no bkg dPhi, 1.5<pTE<2", bins, 0.4,5.6);
    TH1F* pdPhi3 = new TH1F("pdPhi3", "no bkg dPhi, 2<pTE<2.5", bins, 0.4,5.6);
    TH1F* pdPhi4 = new TH1F("pdPhi4", "no bkg dPhi, 2.5<pTE<3", bins, 0.4,5.6);
    TH1F* pdPhi5 = new TH1F("pdPhi5", "no bkg dPhi, 3<pTE<3.5", bins, 0.4,5.6);
    TH1F* pdPhi6 = new TH1F("pdPhi6", "no bkg dPhi, 3.5<pTE<4", bins, 0.4,5.6);
    TH1F* pdPhi7 = new TH1F("pdPhi7", "no bkg dPhi, 4<pTE<4.5", bins, 0.4,5.6);
    TH1F* pdPhi8 = new TH1F("pdPhi8", "no bkg dPhi, 4.5<pTE<5", bins, 0.4,5.6);
    TH1F* pdPhi9 = new TH1F("pdPhi9", "no bkg dPhi, 5<pTE<5.5", bins, 0.4,5.6);
    TH1F* pdPhi10 = new TH1F("pdPhi10", "no bkg dPhi, 5.5<pTE<6", bins, 0.4,5.6);
    TH1F* pdPhi11 = new TH1F("pdPhi11", "no bkg dPhi, 6<pTE<6.5", bins, 0.4,5.6);
    TH1F* pdPhi12 = new TH1F("pdPhi12", "no bkg dPhi, 6.5<pTE<7", bins, 0.4,5.6);
    TH1F* pdPhi13 = new TH1F("pdPhi13", "no bkg dPhi, 7<pTE<7.5", bins, 0.4,5.6);
    TH1F* pdPhi14 = new TH1F("pdPhi14", "no bkg dPhi, 7.5<pTE<8", bins, 0.4,5.6);
    TH1F* pdPhi15 = new TH1F("pdPhi15", "no bkg dPhi, 8<pTE<8.5", bins, 0.4,5.6);
    TH1F* pdPhi16 = new TH1F("pdPhi16", "no bkg dPhi, 8.5<pTE<9", bins, 0.4,5.6);
    TH1F* pdPhi17 = new TH1F("pdPhi17", "no bkg dPhi, 9<pTE<9.5", bins, 0.4,5.6);
    TH1F* pdPhi18 = new TH1F("pdPhi18", "no bkg dPhi, 9.5<pTE<10", bins, 0.4,5.6);

    pairsDirty->Draw("dPhi>>pdPhi1","ptMu>3 && ptE>1 && ptE<1.5 && 0.4<dPhi && 5.6>dPhi");
    pairsDirty->Draw("dPhi>>pdPhi2","ptMu>3 && ptE>1.5 && ptE<2 && 0.4<dPhi && 5.6>dPhi");
    pairsDirty->Draw("dPhi>>pdPhi3","ptMu>3 && ptE>2 && ptE<2.5 && 0.4<dPhi && 5.6>dPhi");
    pairsDirty->Draw("dPhi>>pdPhi4","ptMu>3 && ptE>2.5 && ptE<3 && 0.4<dPhi && 5.6>dPhi");
    pairsDirty->Draw("dPhi>>pdPhi5","ptMu>3 && ptE>3 && ptE<3.5 && 0.4<dPhi && 5.6>dPhi");
    pairsDirty->Draw("dPhi>>pdPhi6","ptMu>3 && ptE>3.5 && ptE<4 && 0.4<dPhi && 5.6>dPhi");
    pairsDirty->Draw("dPhi>>pdPhi7","ptMu>3 && ptE>4 && ptE<4.5 && 0.4<dPhi && 5.6>dPhi");
    pairsDirty->Draw("dPhi>>pdPhi8","ptMu>3 && ptE>4.5 && ptE<5 && 0.4<dPhi && 5.6>dPhi");
    pairsDirty->Draw("dPhi>>pdPhi9","ptMu>3 && ptE>5 && ptE<5.5 && 0.4<dPhi && 5.6>dPhi");
    pairsDirty->Draw("dPhi>>pdPhi10","ptMu>3 && ptE>5.5 && ptE<6 && 0.4<dPhi && 5.6>dPhi");
    pairsDirty->Draw("dPhi>>pdPhi11","ptMu>3 && ptE>6 && ptE<6.5 && 0.4<dPhi && 5.6>dPhi");
    pairsDirty->Draw("dPhi>>pdPhi12","ptMu>3 && ptE>6.5 && ptE<7 && 0.4<dPhi && 5.6>dPhi");
    pairsDirty->Draw("dPhi>>pdPhi13","ptMu>3 && ptE>7 && ptE<7.5 && 0.4<dPhi && 5.6>dPhi");
    pairsDirty->Draw("dPhi>>pdPhi14","ptMu>3 && ptE>7.5 && ptE<8 && 0.4<dPhi && 5.6>dPhi");
    pairsDirty->Draw("dPhi>>pdPhi15","ptMu>3 && ptE>8 && ptE<8.5 && 0.4<dPhi && 5.6>dPhi");
    pairsDirty->Draw("dPhi>>pdPhi16","ptMu>3 && ptE>8.5 && ptE<9 && 0.4<dPhi && 5.6>dPhi");
    pairsDirty->Draw("dPhi>>pdPhi17","ptMu>3 && ptE>9 && ptE<9.5 && 0.4<dPhi && 5.6>dPhi");
    pairsDirty->Draw("dPhi>>pdPhi18","ptMu>3 && ptE>9.5 && ptE<10 && 0.4<dPhi && 5.6>dPhi");


    std::vector<Double_t> widthsD;
    std::vector<Double_t> widthsDErr;
    std::vector<Double_t> integD;
    std::vector<Double_t> integDErr;
    count = 0;

    TH1F* allDirtyHists[18] = {pdPhi1,pdPhi2,pdPhi3,pdPhi4,pdPhi5,pdPhi6,pdPhi7,pdPhi8,pdPhi9,pdPhi10,pdPhi11,pdPhi12,pdPhi13,pdPhi14,pdPhi15,pdPhi16,pdPhi17,pdPhi18};

    for (TH1F* hist: allDirtyHists) 
    {   
        cout<<endl<<count+1<<endl;


        hist->Draw();
        fSpectrum->SetParameters(100,3,0.6,11);
        TFitResultPtr r = hist->Fit("fSpectrum","S","same",0.4,5.6);

        widthsD.push_back(std::abs(fSpectrum->GetParameter(2)*2.35482));
        widthsDErr.push_back(fSpectrum->GetParError(2)*2.35482);

        fSpectrum->SetParameter(2,std::abs(fSpectrum->GetParameter(2)));
        bkgInteg = fSpectrum->GetParameter(3)*5.2;
        
        Double_t integ = (fSpectrum->Integral(0.4,5.6))-bkgInteg;
        Double_t fac = (1/(2*M_PI*xArr[count])) * (sigma/0.5);

        integD.push_back(fac*integ);
        integDErr.push_back(fac*fSpectrum->IntegralError(0.4,5.6,r->GetParams(), r->GetCovarianceMatrix().GetMatrixArray()));
        count++;
    }

    Double_t wdArr[widthsD.size()];
    std::copy(widthsD.begin(), widthsD.end(), wdArr);
    Double_t wdErrArr[widthsDErr.size()];
    std::copy(widthsDErr.begin(), widthsDErr.end(), wdErrArr);

    TGraphErrors* widthCompDirty = new TGraphErrors(18,xArr,wdArr,0,wdErrArr);

    Double_t idArr[integD.size()];
    std::copy(integD.begin(), integD.end(), idArr);
    Double_t idErrArr[integDErr.size()];
    std::copy(integDErr.begin(), integDErr.end(), idErrArr);

    TGraphErrors* integCompDirty = new TGraphErrors(18,xArr,idArr,0,idErrArr);

////// c: With Background
    //TCanvas* rDet = new TCanvas("rDet","dPhi Fit Ranges (with background)");
    cout<<endl<<endl<<"WITH BACKGROUND"<<endl<<endl;
/*
////Muons
    //define and fill histogram objects to use in the fits
    TH1F* pdetPhi1 = new TH1F("pdetPhi1", "det dPhi, 1<ptMu<1.5", bins, 0.4,5.6);
    TH1F* pdetPhi2 = new TH1F("pdetPhi2", "det dPhi, 1.5<ptMu<2", bins, 0.4,5.6);
    TH1F* pdetPhi3 = new TH1F("pdetPhi3", "det dPhi, 2<ptMu<2.5", bins, 0.4,5.6);
    TH1F* pdetPhi4 = new TH1F("pdetPhi4", "det dPhi, 2.5<ptMu<3", bins, 0.4,5.6);
    TH1F* pdetPhi5 = new TH1F("pdetPhi5", "det dPhi, 3<ptMu<3.5", bins, 0.4,5.6);
    TH1F* pdetPhi6 = new TH1F("pdetPhi6", "det dPhi, 3.5<ptMu<4", bins, 0.4,5.6);
    TH1F* pdetPhi7 = new TH1F("pdetPhi7", "det dPhi, 4<ptMu<4.5", bins, 0.4,5.6);
    TH1F* pdetPhi8 = new TH1F("pdetPhi8", "det dPhi, 4.5<ptMu<5", bins, 0.4,5.6);
    TH1F* pdetPhi9 = new TH1F("pdetPhi9", "det dPhi, 5<ptMu<5.5", bins, 0.4,5.6);
    TH1F* pdetPhi10 = new TH1F("pdetPhi10", "det dPhi, 5.5<ptMu<6", bins, 0.4,5.6);
    TH1F* pdetPhi11 = new TH1F("pdetPhi11", "det dPhi, 6<ptMu<6.5", bins, 0.4,5.6);
    TH1F* pdetPhi12 = new TH1F("pdetPhi12", "det dPhi, 6.5<ptMu<7", bins, 0.4,5.6);
    TH1F* pdetPhi13 = new TH1F("pdetPhi13", "det dPhi, 7<ptMu<7.5", bins, 0.4,5.6);
    TH1F* pdetPhi14 = new TH1F("pdetPhi14", "det dPhi, 7.5<ptMu<8", bins, 0.4,5.6);
    TH1F* pdetPhi15 = new TH1F("pdetPhi15", "det dPhi, 8<ptMu<8.5", bins, 0.4,5.6);
    TH1F* pdetPhi16 = new TH1F("pdetPhi16", "det dPhi, 8.5<ptMu<9", bins, 0.4,5.6);
    TH1F* pdetPhi17 = new TH1F("pdetPhi17", "det dPhi, 9<ptMu<9.5", bins, 0.4,5.6);
    TH1F* pdetPhi18 = new TH1F("pdetPhi18", "det dPhi, 9.5<ptMu<10", bins, 0.4,5.6);

    pairsDetector->Draw("dPhi>>pdetPhi1","ptE>3 && ptMu>1 && ptMu<1.5 && 0.4<dPhi && 5.6>dPhi");
    pairsDetector->Draw("dPhi>>pdetPhi2","ptE>3 && ptMu>1.5 && ptMu<2 && 0.4<dPhi && 5.6>dPhi");
    pairsDetector->Draw("dPhi>>pdetPhi3","ptE>3 && ptMu>2 && ptMu<2.5 && 0.4<dPhi && 5.6>dPhi");
    pairsDetector->Draw("dPhi>>pdetPhi4","ptE>3 && ptMu>2.5 && ptMu<3 && 0.4<dPhi && 5.6>dPhi");
    pairsDetector->Draw("dPhi>>pdetPhi5","ptE>3 && ptMu>3 && ptMu<3.5 && 0.4<dPhi && 5.6>dPhi");
    pairsDetector->Draw("dPhi>>pdetPhi6","ptE>3 && ptMu>3.5 && ptMu<4 && 0.4<dPhi && 5.6>dPhi");
    pairsDetector->Draw("dPhi>>pdetPhi7","ptE>3 && ptMu>4 && ptMu<4.5 && 0.4<dPhi && 5.6>dPhi");
    pairsDetector->Draw("dPhi>>pdetPhi8","ptE>3 && ptMu>4.5 && ptMu<5 && 0.4<dPhi && 5.6>dPhi");
    pairsDetector->Draw("dPhi>>pdetPhi9","ptE>3 && ptMu>5 && ptMu<5.5 && 0.4<dPhi && 5.6>dPhi");
    pairsDetector->Draw("dPhi>>pdetPhi10","ptE>3 && ptMu>5.5 && ptMu<6 && 0.4<dPhi && 5.6>dPhi");
    pairsDetector->Draw("dPhi>>pdetPhi11","ptE>3 && ptMu>6 && ptMu<6.5 && 0.4<dPhi && 5.6>dPhi");
    pairsDetector->Draw("dPhi>>pdetPhi12","ptE>3 && ptMu>6.5 && ptMu<7 && 0.4<dPhi && 5.6>dPhi");
    pairsDetector->Draw("dPhi>>pdetPhi13","ptE>3 && ptMu>7 && ptMu<7.5 && 0.4<dPhi && 5.6>dPhi");
    pairsDetector->Draw("dPhi>>pdetPhi14","ptE>3 && ptMu>7.5 && ptMu<8 && 0.4<dPhi && 5.6>dPhi");
    pairsDetector->Draw("dPhi>>pdetPhi15","ptE>3 && ptMu>8 && ptMu<8.5 && 0.4<dPhi && 5.6>dPhi");
    pairsDetector->Draw("dPhi>>pdetPhi16","ptE>3 && ptMu>8.5 && ptMu<9 && 0.4<dPhi && 5.6>dPhi");
    pairsDetector->Draw("dPhi>>pdetPhi17","ptE>3 && ptMu>9 && ptMu<9.5 && 0.4<dPhi && 5.6>dPhi");
    pairsDetector->Draw("dPhi>>pdetPhi18","ptE>3 && ptMu>9.5 && ptMu<10 && 0.4<dPhi && 5.6>dPhi");
*/
////Electrons
    TH1F* pdetPhi1 = new TH1F("pdetPhi1", "det dPhi, 1<pTE<1.5", bins, 0.4,5.6);
    TH1F* pdetPhi2 = new TH1F("pdetPhi2", "det dPhi, 1.5<pTE<2", bins, 0.4,5.6);
    TH1F* pdetPhi3 = new TH1F("pdetPhi3", "det dPhi, 2<pTE<2.5", bins, 0.4,5.6);
    TH1F* pdetPhi4 = new TH1F("pdetPhi4", "det dPhi, 2.5<pTE<3", bins, 0.4,5.6);
    TH1F* pdetPhi5 = new TH1F("pdetPhi5", "det dPhi, 3<pTE<3.5", bins, 0.4,5.6);
    TH1F* pdetPhi6 = new TH1F("pdetPhi6", "det dPhi, 3.5<pTE<4", bins, 0.4,5.6);
    TH1F* pdetPhi7 = new TH1F("pdetPhi7", "det dPhi, 4<pTE<4.5", bins, 0.4,5.6);
    TH1F* pdetPhi8 = new TH1F("pdetPhi8", "det dPhi, 4.5<pTE<5", bins, 0.4,5.6);
    TH1F* pdetPhi9 = new TH1F("pdetPhi9", "det dPhi, 5<pTE<5.5", bins, 0.4,5.6);
    TH1F* pdetPhi10 = new TH1F("pdetPhi10", "det dPhi, 5.5<pTE<6", bins, 0.4,5.6);
    TH1F* pdetPhi11 = new TH1F("pdetPhi11", "det dPhi, 6<pTE<6.5", bins, 0.4,5.6);
    TH1F* pdetPhi12 = new TH1F("pdetPhi12", "det dPhi, 6.5<pTE<7", bins, 0.4,5.6);
    TH1F* pdetPhi13 = new TH1F("pdetPhi13", "det dPhi, 7<pTE<7.5", bins, 0.4,5.6);
    TH1F* pdetPhi14 = new TH1F("pdetPhi14", "det dPhi, 7.5<pTE<8", bins, 0.4,5.6);
    TH1F* pdetPhi15 = new TH1F("pdetPhi15", "det dPhi, 8<pTE<8.5", bins, 0.4,5.6);
    TH1F* pdetPhi16 = new TH1F("pdetPhi16", "det dPhi, 8.5<pTE<9", bins, 0.4,5.6);
    TH1F* pdetPhi17 = new TH1F("pdetPhi17", "det dPhi, 9<pTE<9.5", bins, 0.4,5.6);
    TH1F* pdetPhi18 = new TH1F("pdetPhi18", "det dPhi, 9.5<pTE<10", bins, 0.4,5.6);

    pairsDetector->Draw("dPhi>>pdetPhi1","ptMu>3 && ptE>1 && ptE<1.5 && 0.4<dPhi && 5.6>dPhi");
    pairsDetector->Draw("dPhi>>pdetPhi2","ptMu>3 && ptE>1.5 && ptE<2 && 0.4<dPhi && 5.6>dPhi");
    pairsDetector->Draw("dPhi>>pdetPhi3","ptMu>3 && ptE>2 && ptE<2.5 && 0.4<dPhi && 5.6>dPhi");
    pairsDetector->Draw("dPhi>>pdetPhi4","ptMu>3 && ptE>2.5 && ptE<3 && 0.4<dPhi && 5.6>dPhi");
    pairsDetector->Draw("dPhi>>pdetPhi5","ptMu>3 && ptE>3 && ptE<3.5 && 0.4<dPhi && 5.6>dPhi");
    pairsDetector->Draw("dPhi>>pdetPhi6","ptMu>3 && ptE>3.5 && ptE<4 && 0.4<dPhi && 5.6>dPhi");
    pairsDetector->Draw("dPhi>>pdetPhi7","ptMu>3 && ptE>4 && ptE<4.5 && 0.4<dPhi && 5.6>dPhi");
    pairsDetector->Draw("dPhi>>pdetPhi8","ptMu>3 && ptE>4.5 && ptE<5 && 0.4<dPhi && 5.6>dPhi");
    pairsDetector->Draw("dPhi>>pdetPhi9","ptMu>3 && ptE>5 && ptE<5.5 && 0.4<dPhi && 5.6>dPhi");
    pairsDetector->Draw("dPhi>>pdetPhi10","ptMu>3 && ptE>5.5 && ptE<6 && 0.4<dPhi && 5.6>dPhi");
    pairsDetector->Draw("dPhi>>pdetPhi11","ptMu>3 && ptE>6 && ptE<6.5 && 0.4<dPhi && 5.6>dPhi");
    pairsDetector->Draw("dPhi>>pdetPhi12","ptMu>3 && ptE>6.5 && ptE<7 && 0.4<dPhi && 5.6>dPhi");
    pairsDetector->Draw("dPhi>>pdetPhi13","ptMu>3 && ptE>7 && ptE<7.5 && 0.4<dPhi && 5.6>dPhi");
    pairsDetector->Draw("dPhi>>pdetPhi14","ptMu>3 && ptE>7.5 && ptE<8 && 0.4<dPhi && 5.6>dPhi");
    pairsDetector->Draw("dPhi>>pdetPhi15","ptMu>3 && ptE>8 && ptE<8.5 && 0.4<dPhi && 5.6>dPhi");
    pairsDetector->Draw("dPhi>>pdetPhi16","ptMu>3 && ptE>8.5 && ptE<9 && 0.4<dPhi && 5.6>dPhi");
    pairsDetector->Draw("dPhi>>pdetPhi17","ptMu>3 && ptE>9 && ptE<9.5 && 0.4<dPhi && 5.6>dPhi");
    pairsDetector->Draw("dPhi>>pdetPhi18","ptMu>3 && ptE>9.5 && ptE<10 && 0.4<dPhi && 5.6>dPhi");
    
    std::vector<Double_t> widthsDet;
    std::vector<Double_t> widthsDetErr;
    std::vector<Double_t> integDet;
    std::vector<Double_t> integDetErr;
    count=0;

    TH1F* allDetHists[18] = {pdetPhi1,pdetPhi2,pdetPhi3,pdetPhi4,pdetPhi5,pdetPhi6,pdetPhi7,pdetPhi8,pdetPhi9,pdetPhi10,pdetPhi11,pdetPhi12,pdetPhi13,pdetPhi14,pdetPhi15,pdetPhi16,pdetPhi17,pdetPhi18};

    for (TH1F* hist: allDetHists) 
    {   
        cout<<endl<<count+1<<endl;


        hist->Draw();
        fSpectrum->SetParameters(100,3,0.6,11);
        TFitResultPtr r = hist->Fit("fSpectrum","S","same",0.4,5.6);

        widthsDet.push_back(std::abs(fSpectrum->GetParameter(2)*2.35482));
        widthsDetErr.push_back(fSpectrum->GetParError(2)*2.35482);

        fSpectrum->SetParameter(2,std::abs(fSpectrum->GetParameter(2)));
        bkgInteg = fSpectrum->GetParameter(3)*5.2;

        Double_t integ = fSpectrum->Integral(0.4,5.6)-bkgInteg;
        Double_t fac = (1/(2*M_PI*xArr[count])) * (sigma/0.5);

        integDet.push_back(fac*integ);
        integDetErr.push_back(fac*fSpectrum->IntegralError(0.4,5.6,r->GetParams(), r->GetCovarianceMatrix().GetMatrixArray()));
        count++;
    }
    
    Double_t wdetArr[widthsDet.size()];
    std::copy(widthsDet.begin(), widthsDet.end(), wdetArr);
    Double_t wdetErrArr[widthsDetErr.size()];
    std::copy(widthsDetErr.begin(), widthsDetErr.end(), wdetErrArr);

    TGraphErrors* widthCompDet = new TGraphErrors(18,xArr,wdetArr,0,wdetErrArr);

    Double_t idetArr[integDet.size()];
    std::copy(integDet.begin(), integDet.end(), idetArr);
    Double_t idetErrArr[integDetErr.size()];
    std::copy(integDetErr.begin(), integDetErr.end(), idetErrArr);

    TGraphErrors* integCompDet = new TGraphErrors(18,xArr,idetArr,0,idetErrArr);

    cout<<endl<<"Quark Parent"<<endl;
    for(int i=0; i<18; i++){cout<<widthsC[i]<<" +/- "<<widthsCErr[i]<<"  \t\t"<<integC[i]<<" +/- "<<integCErr[i]<<endl;}
    cout<<endl<<"No Bkg"<<endl;
    for(int i=0; i<18; i++){cout<<widthsD[i]<<" +/- "<<widthsDErr[i]<<"  \t\t"<<integD[i]<<" +/- "<<integDErr[i]<<endl;}
    cout<<endl<<"Pion Bkg"<<endl;
    for(int i=0; i<18; i++){cout<<widthsDet[i]<<" +/- "<<widthsDetErr[i]<<"  \t\t"<<integDet[i]<<" +/- "<<integDetErr[i]<<endl;}

// 4: Draw the graph of histogram width vs momentum range of electron
    TCanvas* ranges = new TCanvas("ranges","Comparison: electron momentum range vs histogram");
    
    TMultiGraph *mg2 = new TMultiGraph();
    gPad->SetLogy();
    
    integCompDet->SetMarkerColor(kRed);
    integCompDet->SetLineColor(kRed);
    integFit->SetLineColor(kRed);
    integCompDet->SetMarkerStyle(21);

    integCompDirty->SetMarkerColor(kViolet);
    integCompDirty->SetLineColor(kViolet);
    integFit->SetLineColor(kViolet);
    integCompDirty->SetMarkerStyle(21);

    integCompClean->SetMarkerColor(kBlue);
    integCompClean->SetLineColor(kBlue);
    integFit->SetLineColor(kBlue);
    integCompClean->SetMarkerStyle(21);


    mg2->Add(integCompDet);
    mg2->Add(integCompDirty);
    mg2->Add(integCompClean);

    mg2->Draw("AP");
    //mg2->GetXaxis()->SetTitle("p_{t}^{#mu} (GeV)");
    mg2->GetXaxis()->SetTitle("p_{t}^{e} (GeV)");
    mg2->GetYaxis()->SetTitle("(1/2#pi p_{t}) d#sigma/dp_{t} (Integral)");

    ranges->Modified(); 
    ranges->Update();

    TLegend* keyI = new TLegend(1, 1, 1, 1);
    keyI->AddEntry(integCompDet, "Tier 3 Data", "P");
    keyI->AddEntry(integCompDirty, "Tier 2 Data", "P");
    keyI->AddEntry(integCompClean, "Tier 1 Data", "P");
    keyI->Draw();
 

}