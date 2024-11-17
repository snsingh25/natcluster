// Event generation for PHOTPRODUCTION process
//
// Authors: Ilkka Helenius <ilkka.m.helenius@jyu.fi>
// Siddharth Singh
///home/siddharth/HEP_Projects/Project1_PhotoProduction/PYTHIA_EventsJets/Events_PYHTIA
// Keywords: photon beam; photoproduction; photon-photon;
// In case of photon-photon interactions four different contributions are
// present:              ProcessType:
//  - resolved-resolved  1
//  - resolved-direct    2
//  - direct-resolved    3
//  - direct-direct      4
// 
// PYTHIA examples: main70.cc and main69.cc
//
// cd /home/siddharth/HEP_Projects/Project1_PhotoProduction/PYTHIA_EventsJets/Events_PYHTIA
//***************************************************************************************

#include "TH1.h"
#include "TH3.h"
#include "TFile.h"
#include "TStyle.h"
#include "TTree.h"
#include "TGraph.h"
#include "Pythia8/Pythia.h"
#include <iostream>
using namespace Pythia8;

int main() {
		
	// Number of events per run.
	int nEvent = 10000000;	//10e7 events for the final run

	// Generator.
	Pythia pythia;
	// Shorthand for some public members of pythia (also static ones)
	Settings& settings  = pythia.settings;
	const Info& info = pythia.info;

	// Decrease the output.
	pythia.readString("Init:showMultipartonInteractions = on");
	pythia.readString("Init:showChangedSettings = on");
	pythia.readString("Init:showChangedParticleData = on");
	pythia.readString("Next:numberCount = 0");
	pythia.readString("Next:numberShowInfo = 1");
	pythia.readString("Next:numberShowProcess = 0");
	pythia.readString("Next:numberShowEvent = 0");

	// Beam parameters : 63.2(10x100) 114.9(10x275) 141(18x275)(20x250)
	pythia.readString("Beams:frameType = 2");		
	pythia.readString("Beams:idA = 2212");			
	pythia.readString("Beams:idB = 11");				
pythia.readString("Beams:eA = 820");
pythia.readString("Beams:eB = 27.5");
	pythia.readString("PDF:beamB2gamma = on");	

	// Set Event Settings : 
	settings.mode("Photon:ProcessType", 0);													// Switch on : Automatic Mix
	pythia.readString("Photon:Q2max = 1.0");												// Maximal Q2
	pythia.readString("MultipartonInteractions:pT0Ref = 3");				// Use tuned pT0ref for photon-hadron (for LEP2, EIC = 3)
	pythia.readString("PhaseSpace:pTHatMin = 5");									// Limit partonic pThat.
	// pythia.readString("MultipartonInteractions:pTmin = 1.9");
  // pythia.settings.forceParm("PhaseSpace:pTHatMinDiverge",1.0);

	// pythia.readString("SpaceShower:pTmaxMatch = 1");								// Allow emissions up to the kinematical limit
	// pythia.readString("SpaceShower:dipoleRecoil = on");							// Set dipole recoil on. Necessary for DIS + shower.
	// pythia.settings.forceParm("PhaseSpace:pTHatMinDiverge",1.0);		// Extra pT cut to avoid the divergences in the limit pT→0
	// pythia.settings.parm("PhaseSpace:mHatMin",1.0);									// The minimum invariant mass.
	// pythia.readString("SpaceShower:pTmaxMatch = 2");
  // pythia.readString("SigmaElastic:Coulomb = on");  
	// pythia.readString("Photon:Wmin  = 134.0");			// invariant mass_min of gamma-hadron pair
	// pythia.readString("Photon:Wmax  = 277.0");			// invariant mass_max of gamma-hadron pair

	// Switch relevant processes on :
	pythia.readString("HardQCD:all = on");									//For All Resolved Process
	// pythia.readString("HardQCD:gg2gg = on");							//gg2gg - Gluon Induced Events (code 111)
	// pythia.readString("HardQCD:gg2qqbar = on");					//gg2qqbar - (code 112)
	// pythia.readString("HardQCD:qg2qg = on");							//qg2qg - (code 113)
	// pythia.readString("HardQCD:qq2qq = on");							//qq2qq - Quark Induced Events (code 114)
	// pythia.readString("HardQCD:qqbar2gg = on");					//qqbar2gg - (code 115)
	//pythia.readString("HardQCD:qqbar2qqbarNew = on");			//qqbar2qqbarNew (code 116)

	// For photon-parton interaction :
	pythia.readString("PhotonParton:all = on");							// For All Direct Process
	// pythia.readString("PhotonParton:ggm2qqbar = on");		// Scattering g gamma → q qbar, (q = u,d,s) Code 271 (281).
	// pythia.readString("PhotonParton:ggm2ccbar = off");		// Scattering g gamma → c cbar. Code 272 (282).
	// pythia.readString("PhotonParton:ggm2bbbar = off");		// Scattering g gamma → b bbar. Code 273 (283).
	// pythia.readString("PhotonParton:qgm2qg = on");				// Scattering q gamma → q g. Code 274 (284).
	// pythia.readString("PhotonParton:qgm2qgm = on");			// Scattering q gamma → q gamma. Code 275 (285).

	//ROOT for Storing Generated Events :
	TFile *dataevents = new TFile("dataevents.root","recreate");
	std::vector<Float_t> pxC,pyC,pzC,pTC,eC,eTC;												// COMBINED
	std::vector<Float_t> pxR,pyR,pzR,pTR,eR,eTR,outpartonpTR;						// RESOLVED
	std::vector<Float_t> pxD,pyD,pzD,pTD,eD,eTD,outpartonpTD;						// DIRECT
	std::vector<Float_t> bpx,bpy,bpz,be,beT, beta, bpT,inelasticity;		// BEAM ELECTRONS
	std::vector<Int_t> eventidC, eventidR, eventidD, scattering; 				// Event ID for C/R/D
	// TC stores the 4 momentum of all the final state particles belonging to res+dir events:
	TTree *TC = new TTree("eventdata", "Event Information");
	TC->Branch("eventidC", "std::vector<Int_t>",  &eventidC);
	TC->Branch("pxC", "std::vector<Float_t>",  &pxC);
	TC->Branch("pyC", "std::vector<Float_t>",  &pyC);
	TC->Branch("pzC", "std::vector<Float_t>",  &pzC);
	TC->Branch("pTC", "std::vector<Float_t>",  &pTC);
	TC->Branch("eC", "std::vector<Float_t>",  &eC);
	TC->Branch("eTC", "std::vector<Float_t>",  &eTC);
	// TR stores the 4 momentum of all the final state particles belonging to res events:
	// TTree *TR = new TTree("resolvedevents", "4-Momentum of all Charged Final State Particles (Resolved)");
	TC->Branch("eventidR", "std::vector<Int_t>", &eventidR);
	TC->Branch("pxR", "std::vector<Float_t>",  &pxR);
	TC->Branch("pyR", "std::vector<Float_t>",  &pyR);
	TC->Branch("pzR", "std::vector<Float_t>",  &pzR);
	TC->Branch("pTR", "std::vector<Float_t>",  &pTR);
	TC->Branch("eR", "std::vector<Float_t>",  &eR);
	TC->Branch("eTR", "std::vector<Float_t>",  &eTR);
	TC->Branch("outpartonpTR", "std::vector<Float_t>",  &outpartonpTR);
	// TD stores the 4 momentum of all the final state particles belonging to dir events:
	// TTree *TD = new TTree("directevents", "4-Momentum of all Charged Final State Particles (Direct)");
	TC->Branch("eventidD", "std::vector<Int_t>", &eventidD);
	TC->Branch("pxD", "std::vector<Float_t>",  &pxD);
	TC->Branch("pyD", "std::vector<Float_t>",  &pyD);
	TC->Branch("pzD", "std::vector<Float_t>",  &pzD);
	TC->Branch("pTD", "std::vector<Float_t>",  &pTD);
	TC->Branch("eD", "std::vector<Float_t>",  &eD);
	TC->Branch("eTD", "std::vector<Float_t>",  &eTD);
	TC->Branch("outpartonpTD", "std::vector<Float_t>",  &outpartonpTD);
	// TB stores the 4-Momentum of all scattered Beam electrons:
	// TTree *TB = new TTree("beamelectrons", "4-Momentum of all Charged Scattered Electrons");
	TC->Branch("bpx", "std::vector<Float_t>",  &bpx);
	TC->Branch("bpy", "std::vector<Float_t>",  &bpy);
	TC->Branch("bpz", "std::vector<Float_t>",  &bpz);
	TC->Branch("bpT", "std::vector<Float_t>",  &bpT);
	TC->Branch("be", "std::vector<Float_t>",  &be);
	TC->Branch("beta", "std::vector<Float_t>",  &beta);
	TC->Branch("beT", "std::vector<Float_t>",  &beT);
	TC->Branch("inelasticity", "std::vector<Float_t>",  &inelasticity);
	// type of subprocess/scatterring qq,qg or gg
	TC->Branch("scattering", "std::vector<Int_t>",  &scattering);

	double Xt=0;
	bool qq = false, qg = false, gg = false;
	// Parameters for histograms :
	double pTmin = 0.0;
	double pTmax = 5.0;
	int nBinsPT  = 100;
	//ROOT Histogram pT of Direct and Resolved Processes :
	TH1F *pTres = new TH1F("pTresresR","Resolved-resolved contribution for pT distribution; pT  Res-Res ;Events ", nBinsPT, pTmin, pTmax);
	TH1F *pTdir = new TH1F("pTdirresR","Direct-resolved contribution for pT distribution; pT  Dir-Res ;Events", nBinsPT, pTmin, pTmax);

	//Counters :
	int ev=0, rescount = 0, dircount = 0, beame = 0, c=0, chk=0, qcc=0, gcc=0, qgcc=0;
	int dsub=0,rsub=0,csub=0;
	bool dircheck=false, rescheck=false;
	int qq_events=0, gg_events=0, gq_events=0;
	float total_cs=0, qq_cs=0, gg_cs=0, gq_cs=0;
	//y = inelasticity = Beam energy (ben) / Photon Energy (gen)
	double ben, gen, y;

	// Initialize the generator.
	pythia.init();

	// Event Loop
	for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
		
		// Generate next event.
		if (!pythia.next()) continue;
		ev++;

		if (iEvent == 0) {
			pythia.process.list();
			pythia.event.list();
		}

		double weight = info.weight();				//used to Normalize Histograms

		// check the type of subprocess C/D/R
		if (pythia.info.code()!=0){		
			// Combined Process
			eventidC.push_back(1);

			// Resolved Process
			if (info.photonMode() == 1){
				// eventidR.push_back(1);	
				rescheck = true;
				rescount++;
			}
			// else
			// 	eventidR.push_back(0);
			
			// Direct Process
			if (info.photonMode() == 2){
				// eventidD.push_back(1);
				dircheck=true;
				dircount++;
			}
			// else
			// 	eventidD.push_back(0);	
		}

		if(dircheck) eventidD.push_back(1);
			else eventidD.push_back(0);
		if(rescheck) eventidR.push_back(1);
			else eventidR.push_back(0);

		// Cross Section:
		float cs = info.sigmaGen();
		total_cs+=cs;
		int codechk = pythia.info.code();

		// //////////////////// OLD SCATTERING SELECTIONS ///////////////////
		// // qq scattering
		// if (codechk == 112 || codechk == 114 || codechk == 116 ||
		// 		codechk == 281 || codechk == 282 || codechk == 283 || codechk == 285){
		// 	qq_cs+=cs;
		// 	qq_events++;
		// 	scattering.push_back(1);
		// }
		// // gg scattering
		// if (codechk == 111 || codechk == 115){
		// 	gg_cs+=cs;
		// 	gg_events++;
		// 	scattering.push_back(2);
		// }
		// // gq scattering
		// if (codechk == 113 || codechk == 284){
		// 	gq_cs+=cs;
		// 	gq_events++;
		// 	scattering.push_back(3);
		// }

		// qq scattering
		if (codechk == 112 || codechk == 114 || codechk == 116 || codechk == 121 || codechk == 122 || codechk == 123 || codechk == 124 || codechk == 271 || codechk == 272 || codechk == 273 || codechk == 281 || codechk == 282 || codechk == 283 || codechk == 284){
			qq_cs+=cs;
			qq_events++;
			scattering.push_back(1);
		}
		// gg scattering
		if (codechk == 111 || codechk == 115){
			gg_cs+=cs;
			gg_events++;
			scattering.push_back(2);
		}
		// gq scattering
		if (codechk == 113 || codechk == 274 || codechk == 284){
			gq_cs+=cs;
			gq_events++;
			scattering.push_back(3);
		}

		// Particle Loop
		for (int i = 0; i < pythia.event.size(); ++i){

			//Store incoming beam-photon energy (to calculate inelasticity):
			if (pythia.event[i].id()==22 && pythia.event[i].status()==-13) {
				gen = pythia.event[i].e();
			}

			//Skip for scattered Beam Electron   
			if(pythia.event[i].id() == 11 && pythia.event[i].mother1() == 0 && pythia.event[i].mother2() == 0){
				bpx.push_back(pythia.event[i].px());
				bpy.push_back(pythia.event[i].py());
				bpz.push_back(pythia.event[i].pz());
				beT.push_back(pythia.event[i].eT());
				bpT.push_back(pythia.event[i].pT());
				be.push_back(pythia.event[i].e());
				beta.push_back(pythia.event[i].eta()); 

				ben = pythia.event[i].e();
				y=gen/ben;
				inelasticity.push_back(y);
				beame++; continue;
			}
			
			if(dircheck && pythia.event[i].isParton() && pythia.event[i].status()==-71){
				outpartonpTD.push_back(pythia.event[i].pT());
			}
			if(rescheck && pythia.event[i].isParton() && pythia.event[i].status()==-71){
				outpartonpTR.push_back(pythia.event[i].pT());
			}

			//Store 4-momentum of the final state charged particles (Combined Events)
			if (pythia.event[i].isFinal()){ //  && pythia.event[i].isHadron() && pythia.event[i].isCharged()
				//Save the 4mom of all combined final state particles.
				pxC.push_back(pythia.event[i].px());
				pyC.push_back(pythia.event[i].py());
				pzC.push_back(pythia.event[i].pz());
				pTC.push_back(pythia.event[i].pT());
				eC.push_back(pythia.event[i].e());
				eTC.push_back(pythia.event[i].eT());
				csub++;
				// Fill Histograms :
				double pTch = pythia.event[i].pT();
				if (info.photonMode() == 1) pTres->Fill(pTch, weight);
				if (info.photonMode() == 3) pTdir->Fill(pTch, weight);
				// // Counters to check for unecessary events :
				// if (info.photonMode() == 2 || info.photonMode() == 4)

				//Store 4-momentum of the final state charged particles (Resolved Events)
				if (rescheck){
					pxR.push_back(pythia.event[i].px());
					pyR.push_back(pythia.event[i].py());
					pzR.push_back(pythia.event[i].pz());
					pTR.push_back(pythia.event[i].pT());
					eR.push_back(pythia.event[i].e());
					eTR.push_back(pythia.event[i].eT());
					rsub++;			
				}

				//Store 4-momentum of the final state charged particles (Direct Events)
				if (dircheck){
					pxD.push_back(pythia.event[i].px());
					pyD.push_back(pythia.event[i].py());
					pzD.push_back(pythia.event[i].pz());
					pTD.push_back(pythia.event[i].pT());
					eD.push_back(pythia.event[i].e());
					eTD.push_back(pythia.event[i].eT());
					dsub++;			
				}
			}
		}

		TC->Fill();	// TR->Fill();	TD->Fill();	TB->Fill();	

		// Clear Variables for next Run
		eventidC.clear(); eventidR.clear(); eventidD.clear(); 
		pxC.clear(); pyC.clear(); pzC.clear(); pTC.clear(); eC.clear(); eTC.clear();	// Combined Process
		pxR.clear(); pyR.clear(); pzR.clear(); pTR.clear(); eR.clear(); eTR.clear(); outpartonpTR.clear();	// Resolved Process
		pxD.clear(); pyD.clear(); pzD.clear(); pTD.clear(); eD.clear(); eTD.clear(); outpartonpTD.clear();	// Direct Process
		bpx.clear(); bpy.clear(); bpz.clear(); bpT.clear(); be.clear(); beT.clear(); beta.clear();	// Beam electrons
		inelasticity.clear(); scattering.clear();
		dircheck = false, rescheck = false; //flags
	}

	// Show statistics after each run.
	pythia.stat();

	// cout << "Fraction of Resolved and Direct Subprocesses (based on their Subprocess Codes) for "<< nEvent << " Events :" << endl;
	cout << "\nTotal Selected Events : " << ev << endl;
	cout << "Total Beam Remnant electrons (should be same as Events) : " << beame << endl;
	cout << "\nTotal Combined (Res+Dir) Events : " << ev << endl;
	// cout << "	Combined Constituents : " << csub << endl;
	cout << "Total Direct Events : " << dircount
			<< " | Fraction: " << fixed << (dircount/(double)ev)*100.0 << endl;
	// cout << "	Direct Constituents : " << dsub << endl;
	cout << "Total Resolved Events : " << rescount
			<< " | Fraction: " << fixed << (rescount/(double)ev)*100.0 << endl;
	// cout << "	Total Resolved Constituents : " << rsub << endl;

	//cout << fixed;	// To remove Scientific Notation.
	cout << "\n\nSelected Events: " << ev << endl;
	cout << fixed << "QQ Scattering: " << qq_events 
			<< "	|	QQ-Fraction of Total Events: " << fixed << (qq_events/(double)ev)*100.0 << endl;
	cout << fixed << "GG Scattering: " << gg_events
			<< "	|	GG-Fraction of Total Events: " << fixed << (gg_events/(double)ev)*100.0 << endl;
	cout << fixed << "GQ Scattering: " << gq_events
			<< "	|	GQ-Fraction of Total Events: " << fixed << (gq_events/(double)ev)*100.0 << endl;

	cout << "\nTotal Cross Section:" << scientific << (total_cs/ev) << endl;
	cout << scientific << "QQ Cross Section: " << (qq_cs/ev) << "	|	QQ-Fraction of Total Events: " 
			<< fixed << (qq_cs/total_cs)*100.0 << endl;
	cout << scientific << "GG Cross Section: " << (gg_cs/ev) << "	|	GG-Fraction of Total Events: "
			<< fixed << (gg_cs/total_cs)*100.0 << endl;
	cout << scientific << "GQ Cross Section: " << (gq_cs/ev) << "	|	GQ-Fraction of Total Events: "
			<< fixed << (gq_cs/total_cs)*100.0 << endl;

	dataevents->Write();
	dataevents->Close();
	delete dataevents;
}
