
#include <iostream>
#include <cstdio>
#include <tuple>
#include <vector>
#include <stdlib.h>
#include "cmath"
#include <string>

void pp_Analysis()
{
    gROOT->Macro("loadPythia.C");

////////////Reading from root file
    TFile *f = new TFile("ntuples10Mill.root");

    TNtuple *pairsClean = (TNtuple*)f->Get("ntAll_clean");
    TNtuple *pairsDirty = (TNtuple*)f->Get("ntAll_dirty");
    TNtuple *pairsDetector = (TNtuple*)f->Get("ntAll_detector");
    TNtuple *elecsNoBkg = (TNtuple*)f->Get("noisyElecs");
    TNtuple *elecsBkg = (TNtuple*)f->Get("cleanElecs");
    TNtuple *ntMuons = (TNtuple*)f->Get("ntMuons");

////////////Drawing graphs

// 1: Momentum comparison of all electrons/muons
//ISSUE: 
    TCanvas* pTs = new TCanvas("pTs","Momentum Comparison (all electrons and muons)");
    pTs->Divide(1, 2);

    pTs->cd(1);
    gPad->SetLogy();
    elecsBkg->SetLineColor(kRed);
    //elecsBkg->Draw("ptE"); 

    TH1F* ptEBkg = new TH1F("ptE_Bkg", "Momenta of electrons with background", 85, 0, 250);
    elecsBkg->Draw("ptE>>ptEBkg"); 

    //ptEBkg->GetXaxis()->SetRangeUser(0,150);
    //ptEBkg->SetXTitle("p_{t}^{e}");
    //ptEBkg->SetTitleSize(0.06,"X");
    //pTs->Modified(); 
   // pTs->Update();

    elecsNoBkg->SetLineColor(kViolet);
    //elecsNoBkg->Draw("ptE","","same");
    TH1F* ptENoBkg = new TH1F("ptE_NoBkg", "Momenta of electrons with no background", 85, 0, 150);
    elecsNoBkg->Draw("ptE>>ptENoBkg","","same");

    pairsClean->SetLineColor(kBlue);
    //pairsClean->Draw("ptE","","same");
    TH1F* ptEPairs = new TH1F("ptE_Pairs", "Momenta of electrons paired to muons", 85, 0, 150);
    pairsClean->Draw("ptE>>ptEPairs","","same");

    //pTs->GetXaxis()->SetTitle("GeV");
    //elecsBkg->SetTitleSize(0.06,"X");
    //pTs->GetYaxis()->SetTitle("N");
    //elecsBkg->SetTitleSize(0.06,"Y");

    TLegend* key1 = new TLegend(1, 1, 1, 1);
    key1->AddEntry(elecsBkg, "All Electrons (with pion background)", "l");
    key1->AddEntry(elecsNoBkg, "All Electrons (no pion background)", "l");
    key1->AddEntry(pairsClean, "Paired Electrons", "l");
    key1->Draw();

    pTs->cd(2);
    gPad->SetLogy();

    ntMuons->SetLineColor(kRed);
    ntMuons->Draw("ptMu");

    pairsClean->SetLineColor(kBlue);
    pairsClean->Draw("ptMu","","same");

    TLegend* key2 = new TLegend(1, 1, 1, 1);
    key2->AddEntry(ntMuons, "All Muons", "l");
    key2->AddEntry(pairsClean, "Paired Muons", "l");
    key2->Draw();

// 2: Fit total dPhi graphs
    TF1 *fSpectrum = new TF1("fSpectrum","gaus+pol0(3)",0.4,5.6);
    Int_t bins = 70;

    //Fit all potential pairs ptE>1 for visual reference
    TCanvas* fitsAll = new TCanvas("fitsAll","dPhi Fits");
    fitsAll->Divide(1, 3);

    fitsAll->cd(1);
    TH1F* pcPhiAll = new TH1F("pcPhiAll", "Tier 1 Pairs", bins, 0, 2*M_PI);
    //p_{t}^{e}>1 GeV and p_{t}^{#mu}>1 GeV
    pcPhiAll->SetLineColor(kBlue);
    pairsClean->Draw("dPhi>>pcPhiAll","ptE>5 && ptMu>5"); 


    Int_t elecsC = elecsNoBkg->GetEntries("ptE>5"); 
    pcPhiAll->Scale(1./elecsC,"width");

    //fSpectrum->SetParameters(2,3,0.6,1);
    //pcPhiAll->Fit("fSpectrum","S","",0.4,5.6);
    fSpectrum->SetParameters(0.03,3,0.6,0.006);
    fSpectrum->SetLineColor(kBlue);
    //pcPhiAll->Fit("fSpectrum","","",0.4,5.6);
    //Double_t fwhmC = fSpectrum->GetParameter(2)*2.35482;
    //Double_t fwhmCErr = fSpectrum->GetParError(2)*2.35482;


    pcPhiAll->SetXTitle("#Delta#phi (rad)");
    pcPhiAll->SetTitleSize(0.06,"X");
    pcPhiAll->SetYTitle("(1/N_{trig})dN/d#Delta #phi");
    pcPhiAll->SetTitleSize(0.06,"Y");
    gStyle->SetOptStat("emr");
    fitsAll->Modified(); 
    fitsAll->Update();

    


    fitsAll->cd(2);
    TH1F* pdPhiAll = new TH1F("pdPhiAll", "Tier 2 Pairs", bins, 0, 2*M_PI);
    pdPhiAll->SetLineColor(kViolet);
    pairsDirty->Draw("dPhi>>pdPhiAll","ptE>5 && ptMu>5");  

    
    pdPhiAll->Scale(1./elecsC,"width"); 

    //fSpectrum->SetParameters(2,3,0.6,1);
    fSpectrum->SetParameters(0.03,3,0.6,0.006);
    fSpectrum->SetLineColor(kViolet);
    //pdPhiAll->Fit("fSpectrum","","",0.4,5.6);

    pdPhiAll->SetXTitle("#Delta#phi (rad)");
    pdPhiAll->SetTitleSize(0.06,"X");
    pdPhiAll->SetYTitle("(1/N_{trig})dN/d#Delta #phi");
    pdPhiAll->SetTitleSize(0.06,"Y");
    gStyle->SetOptStat("emr");
    fitsAll->Modified(); 
    fitsAll->Update();

    //Double_t fwhmD = fSpectrum->GetParameter(2)*2.35482;
    //Double_t fwhmDErr = fSpectrum->GetParError(2)*2.35482;


    fitsAll->cd(3);
    TH1F* pdetPhiAll = new TH1F("pdetPhiAll", "Tier 3 Pairs", bins, 0, 2*M_PI);
    pdetPhiAll->SetLineColor(kRed);
    pairsDetector->Draw("dPhi>>pdetPhiAll","ptE>5 && ptMu>5");   

    Int_t elecsDet = elecsBkg->GetEntries("ptE>5"); //change this to total number of electrons recorded once current run finishes
    pdetPhiAll->Scale(1./elecsDet,"width");

    //fSpectrum->SetParameters(2,3,0.6,1);
    //pdetPhiAll->Fit("fSpectrum","","",0.4,5.6);
    fSpectrum->SetParameters(0.0003,3,0.6,0.0001);
    fSpectrum->SetLineColor(kRed);
   // pdetPhiAll->Fit("fSpectrum","","",0.4,5.6);

    pdetPhiAll->SetXTitle("#Delta#phi (rad)");
    pdetPhiAll->SetTitleSize(0.06,"X");
    pdetPhiAll->SetYTitle("(1/N_{trig})dN/d#Delta #phi");
    pdetPhiAll->SetTitleSize(0.06,"Y");
    gStyle->SetOptStat("emr");
    fitsAll->Modified(); 
    fitsAll->Update();

    //Double_t fwhmDet = fSpectrum->GetParameter(2)*2.35482;
    //Double_t fwhmDetErr = fSpectrum->GetParError(2)*2.35482;

    //cout<<"Width (confirmed quark parent): "<<fwhmC<<" +/- "<<fwhmCErr<<endl;
    //cout<<"Width (no background): "<<fwhmD<<" +/- "<<fwhmDErr<<endl;
    //cout<<"Width (with background): "<<fwhmDet<<" +/- "<<fwhmDetErr<<endl;
//

    //TCanvas* c5 = new TCanvas("c5","dPhi Comparison");

    //pdetPhiAll->Draw();
    //pdPhiAll->Draw("same");
    //pcPhiAll->Draw("same");
    //pairsClean->Fit("gaus","dPhi");


// 3: Fit ranges dPhi graphs
    //TCanvas* rAll = new TCanvas("rAll","dPhi Fit Ranges (charm quark parent)");
    Double_t xArr[18] = {1.25,1.75,2.25,2.75,3.25,3.75,4.25,4.75,5.25,5.75,6.25,6.75,7.25,7.75,8.25,8.75,9.25,9.75};
    TGraphErrors* numElecs = new TGraphErrors(18);
    TGraphErrors* numBkgElecs = new TGraphErrors(18);

    Int_t numElecsArr[18];
    Int_t numBkgElecsArr[18];

    for (int i=0; i<18; i++) 
    {
        Double_t minCut = (Double_t)xArr[i]-0.25;
        Double_t maxCut = (Double_t)xArr[i]+0.25;
        string min_str = to_string(minCut);
        Int_t n = min_str.length();
        char min_ch[n + 1];
        strcpy(min_ch, min_str.c_str());
        char req[5+n];
        strcpy(req, "ptE>");
        strcat(req,min_ch);
        strcat(req,"&& ptE<");
        string max_str = to_string(maxCut);
        Int_t n2 = max_str.length();
        char max_ch[n2 + 1];
        strcpy(max_ch, max_str.c_str());
        strcat(req,max_ch);

        Int_t num = elecsNoBkg->GetEntries(req);
        Int_t numBkg = elecsBkg->GetEntries(req);
        numElecs->SetPoint(i,xArr[i],num);
        numBkgElecs->SetPoint(i,xArr[i],numBkg);
        numElecsArr[i] = num;
        numBkgElecsArr[i] = numBkg;
    }

////// a: Charm Quark Parent
    TCanvas* rAll = new TCanvas("rAll","dPhi Fit Ranges (charm quark parent)");
    cout<<endl<<endl<<"CHARM QUARK PARENT"<<endl<<endl;

    //define and fill histogram objects to use in the fits
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
        hist->Scale(1./numElecsArr[count],"width"); 
        fSpectrum->SetParameters(0.02,3,0.6,0.0001);
        //fSpectrum->SetParameters(10,3,0.6,11);
        //hist->Fit("fSpectrum","S","same",0.4,5.6);
        TFitResultPtr r = hist->Fit("fSpectrum","S","same",0.4,5.6);

        widthsC.push_back(std::abs(fSpectrum->GetParameter(2)*2.35482));
        widthsCErr.push_back(fSpectrum->GetParError(2)*2.35482);

        fSpectrum->SetParameter(2,std::abs(fSpectrum->GetParameter(2)));
        bkgInteg = fSpectrum->GetParameter(3)*5.2;
        integC.push_back(fSpectrum->Integral(0.4,5.6)-bkgInteg);
        //integCErr.push_back(fSpectrum->IntegralError(0.4,5.6,fSpectrum->GetParameters(),fSpectrum->GetParErrors()));
        integCErr.push_back(fSpectrum->IntegralError(0.4,5.6,r->GetParams(), r->GetCovarianceMatrix().GetMatrixArray()));
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

    //define and fill histogram objects to use in the fits
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
        hist->Scale(1./numElecsArr[count],"width"); 
        fSpectrum->SetParameters(0.02,3,0.6,0.0001);
        //fSpectrum->SetParameters(10,3,0.6,11);
        TFitResultPtr r = hist->Fit("fSpectrum","S","same",0.4,5.6);

        widthsD.push_back(std::abs(fSpectrum->GetParameter(2)*2.35482));
        widthsDErr.push_back(fSpectrum->GetParError(2)*2.35482);

        fSpectrum->SetParameter(2,std::abs(fSpectrum->GetParameter(2)));
        bkgInteg = fSpectrum->GetParameter(3)*5.2;
        integD.push_back(fSpectrum->Integral(0.4,5.6)-bkgInteg);
        //integCErr.push_back(fSpectrum->IntegralError(0.4,5.6,fSpectrum->GetParameters(),fSpectrum->GetParErrors()));
        integDErr.push_back(fSpectrum->IntegralError(0.4,5.6,r->GetParams(), r->GetCovarianceMatrix().GetMatrixArray()));
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

    //define and fill histogram objects to use in the fits
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
        hist->Scale(1./numBkgElecsArr[count],"width"); 
        fSpectrum->SetParameters(0.02,3,0.6,0.0001);
        //fSpectrum->SetParameters(10,3,0.6,11);
        TFitResultPtr r = hist->Fit("fSpectrum","S","same",0.4,5.6);

        widthsDet.push_back(std::abs(fSpectrum->GetParameter(2)*2.35482));
        widthsDetErr.push_back(fSpectrum->GetParError(2)*2.35482);

        fSpectrum->SetParameter(2,std::abs(fSpectrum->GetParameter(2)));
        bkgInteg = fSpectrum->GetParameter(3)*5.2;
        integDet.push_back(fSpectrum->Integral(0.4,5.6)-bkgInteg);
        //integCErr.push_back(fSpectrum->IntegralError(0.4,5.6,fSpectrum->GetParameters(),fSpectrum->GetParErrors()));
        integDetErr.push_back(fSpectrum->IntegralError(0.4,5.6,r->GetParams(), r->GetCovarianceMatrix().GetMatrixArray()));
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
    TF1 *integFit = new TF1("integFit","expo+pol0(2)",1,10);
    //TF1 *integFit = new TF1("integFit","pol1",1,10);
    TCanvas* ranges = new TCanvas("ranges","Comparison: electron momentum range vs histogram");
    
    ranges->Divide(1, 2);
    
////Widths
    ranges->cd(1);
    TMultiGraph *mg = new TMultiGraph();

    widthCompDet->SetMarkerColor(kRed);
    widthCompDet->SetLineColor(kRed);
    widthCompDet->SetMarkerStyle(21);

    widthCompDirty->SetMarkerColor(kViolet);
    widthCompDirty->SetLineColor(kViolet);
    widthCompDirty->SetMarkerStyle(21);

    widthCompClean->SetMarkerColor(kBlue);
    widthCompClean->SetLineColor(kBlue);
    widthCompClean->SetMarkerStyle(21);

    mg->Add(widthCompDet);
    mg->Add(widthCompDirty);
    mg->Add(widthCompClean);

    mg->Draw("AP");
    mg->GetXaxis()->SetTitle("p_{t}^{e}(GeV)");
    mg->GetYaxis()->SetTitle("#Delta #phi FWHM");
    mg->GetXaxis()->SetTitleSize(0.05);
    mg->GetYaxis()->SetTitleSize(0.05);
    ranges->Modified(); 
    ranges->Update();

    TLegend* key = new TLegend(1,1,1,1);
    key->AddEntry(widthCompDet, "Tier 3 Data", "P");
    key->AddEntry(widthCompDirty, "Tier 2 Data", "P");
    key->AddEntry(widthCompClean, "Tier 1 Data", "P");
    key->Draw();

////Integrals
    ranges->cd(2);
    //TCanvas* ranges2 = new TCanvas("ranges2","Comparison: electron momentum range vs histogram 2");
    TMultiGraph *mg2 = new TMultiGraph();
    //ranges->SetLogy();
    gPad->SetLogy();
    
    integCompDet->SetMarkerColor(kRed);
    integCompDet->SetLineColor(kRed);
    integFit->SetLineColor(kRed);
    integCompDet->SetMarkerStyle(21);
    integFit->SetParameters(0.000002,-1,0.000100);
    //integCompDet->Fit("integFit","","same",1,10);

    integCompDirty->SetMarkerColor(kViolet);
    integCompDirty->SetLineColor(kViolet);
    integFit->SetLineColor(kViolet);
    integCompDirty->SetMarkerStyle(21);
    //integCompDirty->Fit("integFit","","same",1,10);

    integCompClean->SetMarkerColor(kBlue);
    integCompClean->SetLineColor(kBlue);
    integFit->SetLineColor(kBlue);
    integCompClean->SetMarkerStyle(21);
    //integCompClean->Fit("integFit","","same",1,10);

    numElecs->SetMarkerColor(kGreen);
    numElecs->SetMarkerStyle(22);

    mg2->Add(integCompDet);
    mg2->Add(integCompDirty);
    mg2->Add(integCompClean);
    //mg2->Add(numElecs);

    mg2->Draw("AP");
    mg2->GetXaxis()->SetTitle("p_{t}^{e} (GeV)");
    mg2->GetYaxis()->SetTitle("(1/N_{trig})dN/d#Delta #phi (Integral)");
    mg2->GetXaxis()->SetTitleSize(0.05);
    mg2->GetYaxis()->SetTitleSize(0.05);
    ranges->Modified(); 
    ranges->Update();

    TLegend* keyI = new TLegend(1, 1, 1, 1);
    keyI->AddEntry(widthCompDet, "Tier 3 Data", "P");
    keyI->AddEntry(widthCompDirty, "Tier 2 Data", "P");
    keyI->AddEntry(widthCompClean, "Tier 1 Data", "P");
    //keyI->AddEntry(numElecs, "Total (unpaired) Electrons in p_{t} Range", "P");
    keyI->Draw();
 

// 6: Effect of ptE cut on number of pairs
/*
    TH2F* cutsC = new TH2F("cutsC","Effect of Electron and Muon Cuts on number of pairs",9,1,10,9,1,10);
    
    for(int i=1; i<10; i++)
    {
        for(int j=1; j<10; j++)
        {
            string cut_str = to_string(i);
            Int_t n = cut_str.length();
            char cut_ch[n + 1];
            strcpy(cut_ch, cut_str.c_str());
            char req[5+n];
            strcpy(req, "ptE>");
            strcat(req,cut_ch);
            strcat(req,"&& ptMu>");
            string cut_str2 = to_string(j);
            Int_t n2 = cut_str2.length();
            char cut_ch2[n2 + 1];
            strcpy(cut_ch2, cut_str2.c_str());
            strcat(req,cut_ch2);
            //cout<<req<<endl;

            cutsC->Fill(i,j,pairsClean->GetEntries(req));
            //cout<<countC[i-1]<<endl;
        }

    }
    TCanvas* cuts = new TCanvas("cuts","Effect of cuts");
    cutsC->Draw("colz text");
    cutsC->SetXTitle("Electron Cut (GeV)");
    cutsC->SetYTitle("Muon Cut (GeV)");*/

}