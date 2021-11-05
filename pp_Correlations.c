//
//  pp_Correlations.c
//  
// In Root, run ".x pp_Correlations.c"

#include <iostream>
#include <cstdio>
#include <tuple>
#include <vector>
#include <stdlib.h>
#include "cmath"
#include <string>

auto f = TFile::Open("ntuples.root","RECREATE");

//General Tracking Variables (over all events)
//(unideal setup, will fix)
double nElectron = 0;
double nMuon = 0;
double nChElectron = 0;
double nChMuon = 0;
double corrCount = 0;


TNtuple* noisyElecs = new TNtuple("noisyElecs","All electron momentums and some pion background","pdg:ptE");
TNtuple* cleanElecs = new TNtuple("cleanElecs","All electron momentums","pdg:ptE");
TNtuple* ntMuons = new TNtuple("ntMuons","All muon momentums","pdg:ptMu");

//Searches if specific particle has a charm quark parent, returns a boolean and the details of the mother quark
//Previous issue: only checks one of the two listed parent particles, this leads to it missing any cbar parents - correlations can't be found until this is fixed
//Current fix: loop through the first mother as always but store and check the second mother and if it is a quark/antiquark, return that
//This is based on my manual analysis of event listings - as far as I can tell if the second mother becomes relevant it is because that mother is immediately an antiquark

std::tuple<bool, TParticle*, Int_t> QuarkParent(Int_t quarkPDG, TParticle* part, TClonesArray* particles) 
{
    
    Int_t motherNo = part->GetMother(0);
    Int_t mother2No = part->GetMother(1);

    TParticle* mother = (TParticle*) particles->At(motherNo);
    TParticle* mother2 = (TParticle*) particles->At(motherNo); //this line is just to initialize the variable mother2 with something, I define it properly later down - I chose to initialize with motherNo for now to prevent the throwing of an error in the print statements even though that error is solved later on

    Int_t mPDG = mother->GetPdgCode();
    
    while(mPDG!=2212)
    {
        if(mother2No!=-1) //this check condition solves the aforementioned error by not performing operations on the second mother particle if it doesn't exist
        {
            mother2 = (TParticle*) particles->At(mother2No); //correct initialization of mother2

            if(mother2->GetPdgCode()==quarkPDG  or mother2->GetPdgCode()== -quarkPDG)
            {
                return std::tuple<bool, TParticle*, Int_t>{TRUE, mother2, mother2No};
            }
        }


        if((mPDG)==quarkPDG  or (mPDG)== -quarkPDG)
        {
            return std::tuple<bool, TParticle*, Int_t>{TRUE, mother, motherNo};
        }

        motherNo = mother->GetMother(0);
        mother2No = mother->GetMother(1);

        mother = (TParticle*) particles->At(motherNo);
        mother2 = (TParticle*) particles->At(motherNo); //wrong initialization for same reasons as above

        mPDG = mother->GetPdgCode();

        
    }
    return std::tuple<bool, TParticle*, Int_t>{FALSE, mother, motherNo};
}

vector<vector<TParticle*>> FindRelevantParticles(Int_t np, TClonesArray* particles, Int_t pdgQuark, Int_t evNo) 
{
    //Counters for my sanity, delete some later
    Int_t nElecPerEvent = 0;
    Int_t nElecPionPerEvent = 0;
    Int_t nChElecPerEvent = 0;
    Int_t nMuonPerEvent = 0;
    Int_t nChMuonPerEvent = 0;
    
    //PDG Codes
    Int_t pdgElectron= 11;
    Int_t pdgMuon= 13;

    Int_t pdgCharm = 4;
    Int_t pdgBottom = 5;

    Int_t pdgPion = 211;
    Int_t pdgPi0 = 111;
    
    /*//Array Initialization for printed out displays that can be used as checks
    int elecMothers[10] = {0,0,0,0,0,0,0,0,0,0};
    int muonMothers[10] = {0,0,0,0,0,0,0,0,0,0};
    
    int chElecs[10] = {0,0,0,0,0,0,0,0,0,0};
    int chMuons[10] = {0,0,0,0,0,0,0,0,0,0};*/

    //Useful objects - storage for electrons and muons in a vector
    //"rel" = has charm quark parent

    std::vector<TParticle*> dirtyElecs; //includes pion background

    std::vector<TParticle*> allElecs;
    std::vector<TParticle*> relElecs;
    std::vector<Int_t> relMothersE;

    std::vector<TParticle*> allMuons;
    std::vector<TParticle*> relMuons;
    std::vector<Int_t> relMothersM;
     
    bool foundM=FALSE;
    bool foundE=FALSE;

    //Correlation vectors - different levels of "accuracy" in the different vectors 
        // correlatedClean: only electrons and muons that have been confirmed as correlated by looking through pythia event record
        // correlatedDirty: electrons and muons are considered regardless of if they have a confirmed quark parent
        // correlatedDetector: rough reconstruction of what we'll get from the detector (1% pion background included)
    
    std::vector<TParticle*> correlatedClean;
    std::vector<TParticle*> correlatedDirty;
    std::vector<TParticle*> correlatedDetector;

    //Loop over all particles in event
    for (Int_t ip = 0; ip < np; ip++)
    {
        //"part": the specific particle we are examining from the event
        TParticle* part = (TParticle*) particles->At(ip);
        
        //checking if the particle is in its final state
        if(part->GetStatusCode()>0)
        {
    
    //Introducing pions for intentional background
            if (part->GetPdgCode() == pdgPion or part->GetPdgCode() == -pdgPion or part->GetPdgCode() == pdgPi0){
                //Generate a random number between 1 and 100
                int bgChance = rand()%100;

                //Add 1% of the pions into the data
                if(bgChance<=1)
                {
                    nElecPionPerEvent++;

                    //for momentum plot
                    noisyElecs->Fill(part->GetPdgCode(),part->Pt());

                    //for use when populating correlation vectors
                    dirtyElecs.push_back(part);
                }
                

            }
    
    //Electrons
            if (part->GetPdgCode() == pdgElectron or part->GetPdgCode() == -pdgElectron){
                nElecPerEvent++;
                nElecPionPerEvent++;
                nElectron++;

                //for momentum plot
                cleanElecs->Fill(part->GetPdgCode(),part->Pt());
                noisyElecs->Fill(part->GetPdgCode(),part->Pt());
                
                //for use when populating correlation vectors
                allElecs.push_back(part);
                dirtyElecs.push_back(part);

                //Search for Quark Parent
                std::tuple <bool, TParticle*, Int_t> found = QuarkParent(pdgQuark, part, particles);
                bool parentFound = std::get<0>(found);
                
                if (parentFound==TRUE){
                    TParticle* mother = (TParticle*) std::get<1>(found);
                    //chElecs[nChElecPerEvent] = part->GetPdgCode();
                    //elecMothers[nChElecPerEvent] = mother->GetPdgCode();
                    
                    relElecs.push_back(part);
                    relMothersE.push_back(std::get<2>(found));
                    
                    nChElectron++;
                    nChElecPerEvent++;
                    foundE=TRUE;
                }
            }     
            
    //Muons
            if (part->GetPdgCode() == pdgMuon or part->GetPdgCode() == -pdgMuon){
                nMuonPerEvent++; 
                nMuon++;

                //for momentum plot
                ntMuons->Fill(part->GetPdgCode(),part->Pt());

                //for use when populating correlation vectors
                allMuons.push_back(part);

                //Search for Quark Parent
                std::tuple <bool, TParticle*, Int_t> found = QuarkParent(pdgQuark, part, particles);
                bool parentFound = std::get<0>(found);

                //Check to see if particle in question has the desired quark parent
                if (parentFound==TRUE){
                    TParticle* mother = (TParticle*) std::get<1>(found);
                    //chMuons[nChMuonPerEvent] = part->GetPdgCode();
                    //muonMothers[nChMuonPerEvent] = mother->GetPdgCode();
            
                    relMuons.push_back(part);
                    relMothersM.push_back(std::get<2>(found));
                    
                    nChMuon++;
                    nChMuonPerEvent++;
                    foundM=TRUE;
                    
                }
            }
        }
    }

    //Print out the arrays for manual error checking:
    /*if(foundM==TRUE and foundE==TRUE)
    {
    
        cout<<endl;
        cout<<"Relevant particles (potentially) found in event "<<evNo<<endl;
        cout<<endl;
        //cout << nChElecPerEvent << " electrons with quark parents" << endl;
        //cout << nChMuonPerEvent << " muons with quark parents" << endl;
        
        cout<<"Potentially Relevant Electrons: "<<endl;
        for (int i=0;i<nChElecPerEvent;i++) {
              cout << relElecs[i]->GetPdgCode() << "  ";
           }
        cout<<endl;

        cout<<"Electron Mothers: "<<endl;
        for (int element: elecMothers) {
              cout << element << "  ";
           }

        cout<<endl<<"Indices:"<<endl;
        for (int i=0;i<nChElecPerEvent;i++) {
              cout << relMothersE[i] << "  ";
           }

        cout<<endl<<endl;
        
        cout<<"Potentially Relevant Muons: "<<endl;
        for (int i=0;i<nChMuonPerEvent;i++) {
              cout << relMuons[i]->GetPdgCode() << "  ";
           }
        cout<<endl;
        
        cout<<"Muon Mothers: "<<endl;
        for (int element: muonMothers) {
              cout << element << "  ";
           }

        cout<<endl<<endl;
        cout<<endl;
    }*/


    //Populate correlatedClean 
    if(foundM==TRUE and foundE==TRUE) //checks that electrons and muons with charm quark parents have both been found, hence a possible correlation exists
    {
        int count = 0;
        for (int i=0; i<nChMuonPerEvent; i++)
        {
            Int_t mSign = relMuons[i]->GetPdgCode();

            for (int j=0; j<nChElecPerEvent; j++)
            {
                Int_t eSign = relElecs[j]->GetPdgCode();

                if(mSign * eSign < 0) //correlation occurs when opposite sign pair is found
                {
                    //cout<<"potential correlation found"<<endl<<endl;
                    corrCount++;
                    count++;

                    correlatedClean.push_back(relElecs[j]);
                    correlatedClean.push_back(relMuons[i]);

                    //break; //1 electron can only be correlated to 1 muon, break the loop when correlation is found
                }
            }
        }
    }

    //Populate correlatedDirty
    if(nElecPerEvent>0 and nMuonPerEvent>0) //checks that both electrons and muons have been found this event, even if they don't come from relevant quark parents
    {
        int count = 0;
        for (int i=0; i<nMuonPerEvent; i++)
        {
            Int_t mSign = allMuons[i]->GetPdgCode();

            for (int j=0; j<nElecPerEvent; j++)
            {
                Int_t eSign = allElecs[j]->GetPdgCode();

                if(mSign * eSign < 0) //correlation occurs when opposite sign pair is found
                {
                    correlatedDirty.push_back(allElecs[j]);
                    correlatedDirty.push_back(allMuons[i]);

                    //break; //1 electron can only be correlated to 1 muon, break the loop when correlation is found
                }
            }
        }
    }

    //Populate correlatedDetector
    if(nElecPionPerEvent>0 and nMuonPerEvent>0) //checks that both "electrons" and muons have been found this event, with the chance a pion has been misread as an electron
    {
        int count = 0;
        for (int i=0; i<nMuonPerEvent; i++)
        {
            Int_t mSign = allMuons[i]->GetPdgCode();

            for (int j=0; j<nElecPionPerEvent; j++)
            {
                Int_t eSign = dirtyElecs[j]->GetPdgCode();

                if(mSign * eSign < 0) //correlation occurs when opposite sign pair is found
                {
                    correlatedDetector.push_back(dirtyElecs[j]);
                    correlatedDetector.push_back(allMuons[i]);

                    //break; //1 electron can only be correlated to 1 muon, break the loop when correlation is found
                }
            }
        }
    }

    std::vector<std::vector<TParticle*>> returnPairs;
    returnPairs.push_back(correlatedClean);
    returnPairs.push_back(correlatedDirty);
    returnPairs.push_back(correlatedDetector);

    return returnPairs;
}


void pp_Correlations(Float_t energy = 13500, Int_t nev  = 10000000, Int_t ndeb = 1)
{
    gROOT->Macro("loadPythia.C");
    
    // Array of particles
    TClonesArray* particles = new TClonesArray("TParticle", 5000);
    
    // Create pythia8 object
    TPythia8* pythia8 = new TPythia8();
    
    // Configure
    pythia8->ReadString("HardQCD:hardccbar = on");
    //pythia8->ReadString("HardQCD:hardbbbar = on");
    pythia8->ReadString("PhaseSpace:pTHatMin = 20"); //used to be 2
    pythia8->ReadString("MultipartonInteractions:processLevel = 1");
    //pythia8->ReadString("MultipartonInteractions:Kfactor = 3.5"); //potentially relevant pythia setting, also potentially redundant in pythia8, not using for now
    
    // Initialize collision
    pythia8->Initialize(2212 /* p */, 2212 /* p */, energy /* GeV */);

    //PDG Codes
    Int_t pdgCharm = 4;
    Int_t pdgBottom = 5;

//Ntuple initialization
    TNtuple* pairsClean = new TNtuple("ntAll_clean","Statistics of correlated pairs","etaE:etaMu:ptE:ptMu:phiE:phiMu:dPhi:ratio");
    TNtuple* pairsElecNClean = new TNtuple("ntEN_clean","Statistics of correlated pairs","etaE:etaMu:ptE:ptMu:phiE:phiMu:dPhi");
    TNtuple* pairsElecPClean = new TNtuple("ntEP_clean","Statistics of correlated pairs","etaE:etaMu:ptE:ptMu:phiE:phiMu:dPhi");

    TNtuple* pairsDirty = new TNtuple("ntAll_dirty","Statistics of correlated pairs","etaE:etaMu:ptE:ptMu:phiE:phiMu:dPhi:ratio");
    TNtuple* pairsElecNDirty = new TNtuple("ntEN_dirty","Statistics of correlated pairs","etaE:etaMu:ptE:ptMu:phiE:phiMu:dPhi");
    TNtuple* pairsElecPDirty = new TNtuple("ntEP_dirty","Statistics of correlated pairs","etaE:etaMu:ptE:ptMu:phiE:phiMu:dPhi");

    TNtuple* pairsDetector = new TNtuple("ntAll_detector","Statistics of correlated pairs","etaE:etaMu:ptE:ptMu:phiE:phiMu:dPhi:ratio");
    TNtuple* pairsElecNDetector = new TNtuple("ntEN_detector","Statistics of correlated pairs","etaE:etaMu:ptE:ptMu:phiE:phiMu:dPhi");
    TNtuple* pairsElecPDetector = new TNtuple("ntEP_detector","Statistics of correlated pairs","etaE:etaMu:ptE:ptMu:phiE:phiMu:dPhi");

    //TNtuple* noisyElecs = new TNtuple("noisyElecs","All electron momentums and some pion background","pdg:ptE");
    //TNtuple* cleanElecs = new TNtuple("cleanElecs","All electron momentums","pdg:ptE");
    //TNtuple* ntMuons = new TNtuple("ntMuons","All muon momentums","pdg:ptMu");
    // Event loop

    for (Int_t iev = 0; iev < nev; iev++)
    {
        //cout<<"----------------------------EV NO: " <<iev<< endl;
        //Event Tracking Variables
        Int_t nElecPerEvent = 0;
        Int_t nElecPionPerEvent = 0;
        Int_t countCorrEv = 0;
        
        
        pythia8->GenerateEvent();
        //pythia8->EventListing();
        
        pythia8->ImportParticles(particles,"All");
        Int_t np = particles->GetEntriesFast();
        
        
        std::vector<std::vector<TParticle*>> allCorrTypes = FindRelevantParticles(np,particles,pdgCharm,iev);
        //std::vector<std::vector<TParticle*>> allCorrTypes = FindRelevantParticles(np,particles,pdgBottom,iev);
        
        //Confirmed correlations, no background
        std::vector<TParticle*> corrClean = allCorrTypes[0];
        int num = corrClean.size(); 
        if(num>0)
        {
            for(int i=0; i<num; i=i+2)
            {
                Double_t ptE = corrClean[i]->Pt();
                Double_t ptMu = corrClean[i+1]->Pt();
                Double_t phiE = corrClean[i]->Phi();
                Double_t phiMu = corrClean[i+1]->Phi();
                Double_t dPhi = (phiE - phiMu);
                Double_t etaE = corrClean[i]->Eta();
                Double_t etaMu = corrClean[i+1]->Eta();
                Double_t ratio = (ptE-ptMu)/(ptE+ptMu);

                //if(dPhi<-0.5*M_PI) {dPhi = dPhi + 2.*M_PI;}

                //else if(dPhi>1.5*M_PI) {dPhi = dPhi - 2.*M_PI;}

                if(dPhi<0) {dPhi = dPhi + 2.*M_PI;}
                else if(dPhi>2*M_PI) {dPhi = dPhi - 2.*M_PI;}

                pairsClean->Fill(etaE,etaMu,ptE,ptMu,phiE,phiMu,dPhi,ratio);

                if(corrClean[i]->GetPdgCode() == 11) {pairsElecNClean->Fill(etaE,etaMu,ptE,ptMu,phiE,phiMu,dPhi); }

                else {pairsElecPClean->Fill(etaE,etaMu,ptE,ptMu,phiE,phiMu,dPhi);}
            }
        }

        //Unconfirmed correlations, no background
        std::vector<TParticle*> corrDirty = allCorrTypes[1];
        num = corrDirty.size();
        if(num>0)
        {
            for(int i=0; i<num; i=i+2)
            {
                Double_t ptE = corrDirty[i]->Pt();
                Double_t ptMu = corrDirty[i+1]->Pt();
                Double_t phiE = corrDirty[i]->Phi();
                Double_t phiMu = corrDirty[i+1]->Phi();
                Double_t dPhi = (phiE - phiMu);
                Double_t etaE = corrDirty[i]->Eta();
                Double_t etaMu = corrDirty[i+1]->Eta();
                Double_t ratio = (ptE-ptMu)/(ptE+ptMu);

                //if(dPhi<-0.5*M_PI) {dPhi = dPhi + 2.*M_PI;}

                //else if(dPhi>1.5*M_PI) {dPhi = dPhi - 2.*M_PI;}

                if(dPhi<0) {dPhi = dPhi + 2.*M_PI;}
                else if(dPhi>2*M_PI) {dPhi = dPhi - 2.*M_PI;}

                pairsDirty->Fill(etaE,etaMu,ptE,ptMu,phiE,phiMu,dPhi,ratio);

                if(corrDirty[i]->GetPdgCode() == 11) {pairsElecNDirty->Fill(etaE,etaMu,ptE,ptMu,phiE,phiMu,dPhi);}

                else {pairsElecPDirty->Fill(etaE,etaMu,ptE,ptMu,phiE,phiMu,dPhi);}
            }
        }

        //Unconfirmed correlations, pion background
        std::vector<TParticle*> corrDetector = allCorrTypes[2];
        num = corrDetector.size();
        if(num>0)
        {
            for(int i=0; i<num; i=i+2)
            {
                Double_t ptE = corrDetector[i]->Pt();
                Double_t ptMu = corrDetector[i+1]->Pt();
                Double_t phiE = corrDetector[i]->Phi();
                Double_t phiMu = corrDetector[i+1]->Phi();
                Double_t dPhi = (phiE - phiMu);
                Double_t etaE = corrDetector[i]->Eta();
                Double_t etaMu = corrDetector[i+1]->Eta();
                Double_t ratio = (ptE-ptMu)/(ptE+ptMu);

                //if(dPhi<-0.5*M_PI) {dPhi = dPhi + 2.*M_PI;}

                //else if(dPhi>1.5*M_PI) {dPhi = dPhi - 2.*M_PI;}

                if(dPhi<0) {dPhi = dPhi + 2.*M_PI;}
                else if(dPhi>2*M_PI) {dPhi = dPhi - 2.*M_PI;}

                pairsDetector->Fill(etaE,etaMu,ptE,ptMu,phiE,phiMu,dPhi,ratio);

                if(corrDetector[i]->GetPdgCode() == 11) {pairsElecNDetector->Fill(etaE,etaMu,ptE,ptMu,phiE,phiMu,dPhi);}

                else {pairsElecPDirty->Fill(etaE,etaMu,ptE,ptMu,phiE,phiMu,dPhi);}
            }

        }
    }
//
    pythia8->PrintStatistics();
    //Double_t crossErr = pythia8->ReadString("Info::sigmaErr()");
    //cout<<"CROSS SECTION: "<<cross<<" +/- "<<crossErr<<endl;

    pairsClean->Write();
    pairsDirty->Write();
    pairsDetector->Write();
    noisyElecs->Write();
    cleanElecs->Write();
    ntMuons->Write();

//Momentum comparison of all electrons/muons
    TCanvas* c1 = new TCanvas("c1","Momentum Comparison (all electrons and muons)");
    c1->Divide(1, 2);

    c1->cd(1);
    gPad->SetLogy();
    noisyElecs->SetLineColor(kRed);
    noisyElecs->Draw("ptE");

    cleanElecs->SetLineColor(kViolet);
    cleanElecs->Draw("ptE","","same");

    pairsClean->SetLineColor(kBlue);
    pairsClean->Draw("ptE","","same");

    TLegend* key1 = new TLegend(1, 1, 1, 1);
    key1->AddEntry(noisyElecs, "pT of all electrons and background", "l");
    key1->AddEntry(cleanElecs, "pT of all electrons", "l");
    key1->AddEntry(pairsClean, "pT of correlated electrons", "l");
    key1->Draw();

    c1->cd(2);
    gPad->SetLogy();

    ntMuons->SetLineColor(kRed);
    ntMuons->Draw("ptMu");

    pairsClean->SetLineColor(kBlue);
    pairsClean->Draw("ptMu","","same");

    TLegend* key2 = new TLegend(1, 1, 1, 1);
    key2->AddEntry(ntMuons, "pT of all muons", "l");
    key2->AddEntry(pairsClean, "pT of correlated muons", "l");
    key2->Draw();

//dPhi comparison
    TCanvas* c5 = new TCanvas("c5","dPhi Comparison");

    pairsDetector->SetLineColor(kRed);
    pairsDetector->Draw("dPhi","ptE>1");
    

    pairsDirty->SetLineColor(kViolet);
    pairsDirty->Draw("dPhi","ptE>1","same");

    pairsClean->SetLineColor(kBlue);
    pairsClean->Draw("dPhi","ptE>1","same");
    //pairsClean->Fit("gaus","dPhi");

    key1 = new TLegend(1, 1, 1, 1);
    key1->AddEntry(pairsDetector, "Unconfirmed correlation, pion background included", "l");
    key1->AddEntry(pairsDirty, "Unconfirmed correlation, no background included", "l");
    key1->AddEntry(pairsClean, "Confirmed correlation through pythia", "l");
    key1->Draw();

//dPhi fits 

    TF1 *fSpectrum = new TF1("fSpectrum","gaus+pol0(3)",0.4,5.6);
    Int_t bins = 85;

    //Fit all potential pairs ptE>1 for visual reference
    TCanvas* fitsAll = new TCanvas("fitsAll","dPhi Fits");
    fitsAll->Divide(1, 3);

    fitsAll->cd(1);
    TH1F* pcPhiAll = new TH1F("pcPhiAll", "clean dPhi, 1<pTE", bins, 0, 2*M_PI);
    pcPhiAll->SetLineColor(kBlue);
    pairsClean->Draw("dPhi>>pcPhiAll","ptE>1");    
    fSpectrum->SetParameters(1,1,1,1);
    pcPhiAll->Fit("fSpectrum","","",0.4,5.6);
    Double_t fwhm = fSpectrum->GetParameter(2)*2.35482;

    fitsAll->cd(2);
    TH1F* pdPhiAll = new TH1F("pdPhiAll", "no bkg dPhi, 1<pTE", bins, 0, 2*M_PI);
    pdPhiAll->SetLineColor(kViolet);
    pairsDirty->Draw("dPhi>>pdPhiAll","ptE>1");    
    fSpectrum->SetParameters(1,1,1,1);
    pdPhiAll->Fit("fSpectrum","","",0.4,5.6);
    fwhm = fSpectrum->GetParameter(2)*2.35482;

    fitsAll->cd(3);
    TH1F* pdetPhiAll = new TH1F("pdetPhiAll", "bkg dPhi, 1<pTE", bins, 0, 2*M_PI);
    pdetPhiAll->SetLineColor(kRed);
    pairsDetector->Draw("dPhi>>pdetPhiAll","ptE>1");    
    fSpectrum->SetParameters(1,1,1,1);
    pdetPhiAll->Fit("fSpectrum","","",0.4,5.6);
    fwhm = fSpectrum->GetParameter(2)*2.35482;

    cout<<"Width (confirmed quark parent): "<<fwhm<<endl;
    cout<<"Width (no background): "<<fwhm<<endl;
    cout<<"Width (with background): "<<fwhm<<endl;


//Calculated stats
    ULong64_t entriesAll = pairsClean->GetEntries();
    ULong64_t entriesP = pairsElecPClean->GetEntries();
    ULong64_t entriesN = pairsElecNClean->GetEntries();
    //cout << "Entries: " << entriesP << " + "<<entriesN<<" = "<<entriesAll << endl;

    cout<<endl;
    cout << "found " << nElectron << " electrons total" << endl;
    cout << "found " << nChElectron << " electrons with charm quark parents" <<endl;
    double percentE = (nChElectron/nElectron)*100;
    cout << percentE << "%" <<endl<<endl;

    cout << "found " << nMuon << " muons total" << endl;
    cout << "found " << nChMuon << " muons with charm quark parents" <<endl;
    double percentM = (nChMuon/nMuon)*100;
    cout << percentM << "%" <<endl<<endl;

    cout << "found " << corrCount << " correlated pairs from charm quark parents" << endl<<endl;
    double percentCorrE = (corrCount/nElectron)*100;
    double percentCorrM = (corrCount/nMuon)*100;
    double percentCorrChE = (corrCount/nChElectron)*100;
    double percentCorrChM = (corrCount/nChMuon)*100;

    cout<<"of all electrons found, "<<percentCorrE<<"% were correlated to a muon"<<endl;
    cout<< "this means "<<percentCorrChE<<"% of electrons with a charm quark were correlated to a muon"<<endl<<endl;

    cout<<"of all muons found, "<<percentCorrM<<"% were correlated to an electron"<<endl;
    cout<<"this means "<<percentCorrChM<<"% of muons with a charm quark parent were correlated to an electron"<<endl;

    /*cout<<endl;
    cout << "found " << nElectron << " electrons total" << endl;
    cout << "found " << nChElectron << " electrons with bottom quark parents" <<endl;
    double percentE = (nChElectron/nElectron)*100;
    cout << percentE << "%" <<endl<<endl;

    cout << "found " << nMuon << " muons total" << endl;
    cout << "found " << nChMuon << " muons with bottom quark parents" <<endl;
    double percentM = (nChMuon/nMuon)*100;
    cout << percentM << "%" <<endl<<endl;

    cout << "found " << corrCount << " correlated pairs from bottom quark parents" << endl<<endl;
    double percentCorrE = (corrCount/nElectron)*100;
    double percentCorrM = (corrCount/nMuon)*100;
    double percentCorrChE = (corrCount/nChElectron)*100;
    double percentCorrChM = (corrCount/nChMuon)*100;

    cout<<"of all electrons found, "<<percentCorrE<<"% were correlated to a muon"<<endl;
    cout<< "this means "<<percentCorrChE<<"% of electrons with a bottom quark were correlated to a muon"<<endl<<endl;

    cout<<"of all muons found, "<<percentCorrM<<"% were correlated to an electron"<<endl;
    cout<<"this means "<<percentCorrChM<<"% of muons with a bottom quark parent were correlated to an electron"<<endl;*/


}