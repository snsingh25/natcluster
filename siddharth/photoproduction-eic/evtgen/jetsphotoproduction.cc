//*******************************************
// Photoproduction - Jets and SubJets
//
// Siddharth Singh
// Dr. Manjit Kaur, Dr. Ritu Aggarwal
//
// Creating Jets/Subjets at HERA energies
//
// cd /home/siddharth/HEP_Projects/Project1_PhotoProduction/PYTHIA_EventsJets/Jets_PYHTIA
//*******************************************

#include "importheaders.h"

//------------------Reading input data.root file------------------
TFile *myFile = new TFile("/Users/siddharthsingh/Analysis/photoproduction/EvtGen/finalevents/dataeventsEIC3.root", "READ");
// Combined Events :
TTreeReader myReader("eventdata", myFile);
TTreeReaderArray<Int_t> eventidC(myReader, "eventidC");
TTreeReaderArray<Float_t> pxC(myReader, "pxC");
TTreeReaderArray<Float_t> pyC(myReader, "pyC");
TTreeReaderArray<Float_t> pzC(myReader, "pzC");
TTreeReaderArray<Float_t> pTC(myReader, "pTC");
TTreeReaderArray<Float_t> eC(myReader, "eC");
TTreeReaderArray<Float_t> eTC(myReader, "eTC");
// Resolved Events
TTreeReaderArray<Int_t> eventidR(myReader, "eventidR");
TTreeReaderArray<Float_t> pxR(myReader, "pxR");
TTreeReaderArray<Float_t> pyR(myReader, "pyR");
TTreeReaderArray<Float_t> pzR(myReader, "pzR");
TTreeReaderArray<Float_t> pTR(myReader, "pTR");
TTreeReaderArray<Float_t> eR(myReader, "eR");
TTreeReaderArray<Float_t> eTR(myReader, "eTR");
// Direct Events
TTreeReaderArray<Int_t> eventidD(myReader, "eventidD");
TTreeReaderArray<Float_t> pxD(myReader, "pxD");
TTreeReaderArray<Float_t> pyD(myReader, "pyD");
TTreeReaderArray<Float_t> pzD(myReader, "pzD");
TTreeReaderArray<Float_t> pTD(myReader, "pTD");
TTreeReaderArray<Float_t> eD(myReader, "eD");
TTreeReaderArray<Float_t> eTD(myReader, "eTD");
// Reading Beam Electrons, inelasticity and scattering
TTreeReaderArray<Float_t> bpx(myReader, "bpx");
TTreeReaderArray<Float_t> bpy(myReader, "bpy");
TTreeReaderArray<Float_t> bpz(myReader, "bpz");
TTreeReaderArray<Float_t> be(myReader, "be");
TTreeReaderArray<Float_t> beT(myReader, "beT");
// TTreeReaderArray<Float_t> beta(myReader, "beta");
TTreeReaderArray<Float_t> bpT(myReader, "bpT");
TTreeReaderArray<Float_t> inelasticity(myReader, "inelasticity");
TTreeReaderArray<Int_t> scattering(myReader, "scattering");
//---------------------------------------------------------------

//Declaring the Functions :
void combined_jets (const vector<fastjet::PseudoJet> &);									// Jet for Combined (Res+Dir) Processes
void direct_jets (const vector<fastjet::PseudoJet> &);										// Jet for Direct Processes
void resolved_jets (const vector<fastjet::PseudoJet> &);									// Jet for Resolved Processes
void print_jets (const vector<fastjet::PseudoJet> &, int flag);						// Printing Jets and Subjets
void dijets (const vector<fastjet::PseudoJet> &);													// To study Dijets
double jet_shape_diff (const vector<fastjet::PseudoJet> &);								// Differential Jet Shape
double jet_shape_int (const vector<fastjet::PseudoJet> &, double inrad);	// Integrated Jet Shape, for quark and gluon jets
double average(const vector<double> &);																		// used to count the average vaue of an array
double x_obs(const vector<fastjet::PseudoJet> &);													// xgamma obs to differentiate between D&R jets
TVectorD transverse_sphericity(const vector<fastjet::PseudoJet> &);				// transverse sphericity, returns the eigenvalues
 
// Set the respective parameters :toggleVim
// double cc = 0.0005;
double EtMin_C=10.0, R_C=1.0, crapmin=-4.0, crapmax=4.0, ycut_C = 0.0005, dcut_C;				
double EtMin_D=10.0, R_D=1.0, drapmin=-4.0, drapmax=4.0, ycut_D = 0.0005, dcut_D;				
double EtMin_R=10.0, R_R=1.0, rrapmin=-4.0, rrapmax=4.0, ycut_R = 0.0005, dcut_R;				

// differential jet shape
double rr_diff = 0.6;
vector<double> rho_C, rho_Q, rho_G;

// Integrated jet shape
double rr_int = 0.6;
vector<double> psi_C, psi_D, psi_R, psi_Q, psi_G;

// Transverse Sphericity global-variables
double lambda1=0., lambda2=0., ST=0.;

// Define variables to store the mean subjet multiplicity vs ycut
vector<double> n_C, n_R, n_D, n_Q, n_G;
vector<double> etaq, etag;
vector<fastjet::PseudoJet> temp;

// Output jets.root file to store jet data
TFile *datajets = new TFile("datajets.root", "RECREATE");

// COMBINED EVENTS :
std::vector<Float_t> pxcj, pycj, pzcj, ecj, mcj, etcj, ptcj, phicj, etacj, ncjets;
std::vector<Double_t> mean_js_incl, mean_js_diff;
TTree *jetcombined = new TTree("jetcombined", "");
TH1F *ncJets = new TH1F("ncJets", "#Combined Jets; Number of Combined Jets; Events", 10, 0.5, 10);
TH1F *ncsubJets = new TH1F("ncsubJets", "#Combined <mean subJet multiplicity>; Combined nSubjets; Jets", 10, 0.5, 10);
TH2F *etaphi_planeC = new TH2F("etaphi_planeC", "eta-phi plane (Combined Events); eta; phi", 50, -5, 5, 50, 0, 6.282);
// Transverse Sphericity
TH1F *STplot = new TH1F("STplot", "Transverse Sphericity for integrated-jet Events", 50, 0, 1);
TH1F *STplot0 = new TH1F("STplot0", "Transverse Sphericity for 0-jet Events", 50, 0, 1);
TH1F *STplot1 = new TH1F("STplot1", "Transverse Sphericity for 1-jet Events", 50, 0, 1);
TH1F *STplot2 = new TH1F("STplot2", "Transverse Sphericity for 2-jet Events", 50, 0, 1);
TH1F *STplot3 = new TH1F("STplot3", "Transverse Sphericity for 3-jet Events", 50, 0, 1);
TH1F *STplot4 = new TH1F("STplot4", "Transverse Sphericity for 4-jet Events", 50, 0, 1);
TH1F *STlambda1 = new TH1F("STlambda1", "Eigenvalue Lambda 1", 20, 0, 1);
TH1F *STlambda2 = new TH1F("STlambda2", "Eigenvalue Lambda 2", 20, 0, 1);
// Single Jet Events & Jet Shapes
TH1F *singlejetshapeC = new TH1F("singlejetshapeC", "Jet Shape", 100, 0, 10);

// RESOLVED EVENTS :
std::vector<Float_t> pxrj, pyrj, pzrj, erj, mrj, etrj, phirj, etarj, nrjets;
TTree *jetresolved = new TTree("jetresolved", "");
TH1F *nrJets = new TH1F("nrJets", "#Resolved Jets; Number of Resolved Jets; Events", 10, 0.5, 10);
TH1F *nrsubJets = new TH1F("nrsubJets", "#Resolved <mean subJet multiplicity>; Resolved nSubjets; Jets", 10, 0.5, 10);
TH1F *rxobs = new TH1F("rxobs", "", 10, 0, 100);

// DIRECT EVENTS :
std::vector<Float_t> pxdj, pydj, pzdj, edj, mdj, etdj, phidj, etadj, ndjets;
TTree *jetdirect = new TTree("jetdirect", "");
TH1F *ndJets = new TH1F("ndJets", "#Direct Jets; Number of Direct Jets; Events", 10, 0.5, 10);
TH1F *ndsubJets = new TH1F("ndsubJets", "#Direct <mean subJet multiplicity>; Direct nSubjets; Jets", 10, 0.5, 10);
TH1F *dxobs = new TH1F("dxobs", "", 10, 0, 100);

// THICK JETS (Gluon Initiated):
std::vector<Float_t> gpxj, gpyj, gpzj, gej, gmj, getj, gphij, getaj, gnjets;
TTree *jetthick = new TTree("jetthick", "");
TH1F *ngthicksubJets = new TH1F("ngthicksubJets", "#Thick subJets; Number of Thick subJets; Events", 10, 0.5, 10); 
TH2F *etaphi_planethk = new TH2F("etaphi_planethk", "eta-phi plane (Thick Jet Events); eta; phi", 50, -5, 3, 50, 0, 6.282);

// THIN JETS (Quark Initiated):
std::vector<Float_t> qpxj, qpyj, qpzj, qej, qmj, qetj, qphij, qetaj, qnjets; 
TTree *jetthin = new TTree("jetthin", "");
TH1F *nqthinsubJets = new TH1F("nqthinsubJets", "#Thin subJets; Number of Thin subJets; Events", 10, 0.5, 10); 
TH2F *etaphi_planethn = new TH2F("etaphi_planethn", "eta-phi plane (Thin Jet Events); eta; phi", 50, -5, 3, 50, 0, 6.282);

int ncombjets=0, ndirjets=0, nresjets=0;					// Number of Jets
int ncombsubjets=0, ndirsubjets=0, nressubjets=0;	// Number of subJets
int ev=0;		// Number of Events
double in_y, beam_eT, beam_e, beam_pT; float beam_eta;

// variables to know about the type of scattering process
int qq=0, gg=0, gq=0;
bool chk_qq=false, chk_gg=false, chk_gq=false;
// variables for thiick and thin jets scattering process
int thkqq=0, thkgg=0, thkgq=0;
int thnqq=0, thngg=0, thngq=0;
int zz=0, yy=0, count1jet=0, count2jet=0, count3jet=0;
vector<double> dijet1,dijet2,jet1;

int main(){

	vector<fastjet::PseudoJet> input_particles_combined, input_particles_direct, input_particles_resolved;

	//int ev=0;	
	int c=0, d=0, r=0, cs=0, ds=0, rs=0, g=0, q=0, qs=0, gs=0;
	
	//------------------------------------DECLARING JET PARAMETERS------------------------------------
	//COMBINED ENTRIES:
	jetcombined->Branch("combined_px", "std::vector<Float_t>",  &pxcj);
	jetcombined->Branch("combined_py", "std::vector<Float_t>",  &pycj);
	jetcombined->Branch("combined_pz", "std::vector<Float_t>",  &pzcj);
	jetcombined->Branch("combined_E", "std::vector<Float_t>",  &ecj);
	jetcombined->Branch("combined_mass", "std::vector<Float_t>",  &mcj);
	jetcombined->Branch("combined_Et", "std::vector<Float_t>",  &etcj);
	jetcombined->Branch("combined_pt", "std::vector<Float_t>",  &ptcj);
	jetcombined->Branch("combined_phi", "std::vector<Float_t>",  &phicj);
	jetcombined->Branch("combined_eta", "std::vector<Float_t>",  &etacj);
	jetcombined->Branch("combined_njets", "std::vector<Float_t>",  &ncjets);
	//DIRECT ENTRIES:
	jetdirect->Branch("direct_px", "std::vector<Float_t>",  &pxdj);
	jetdirect->Branch("direct_py", "std::vector<Float_t>",  &pydj);
	jetdirect->Branch("direct_pz", "std::vector<Float_t>",  &pzdj);
	jetdirect->Branch("direct_E", "std::vector<Float_t>",  &edj);
	jetdirect->Branch("direct_mass", "std::vector<Float_t>",  &mdj);
	jetdirect->Branch("direct_Et", "std::vector<Float_t>",  &etdj);
	jetdirect->Branch("direct_phi", "std::vector<Float_t>",  &phidj);
	jetdirect->Branch("direct_eta", "std::vector<Float_t>",  &etadj);
	jetdirect->Branch("direct_njets", "std::vector<Float_t>",  &ndjets);
	//RESOLVED ENTRIES:
	jetresolved->Branch("resolved_px", "std::vector<Float_t>",  &pxrj);
	jetresolved->Branch("resolved_py", "std::vector<Float_t>",  &pyrj);
	jetresolved->Branch("resolved_pz", "std::vector<Float_t>",  &pzrj);
	jetresolved->Branch("resolved_E", "std::vector<Float_t>",  &erj);
	jetresolved->Branch("resolved_mass", "std::vector<Float_t>",  &mrj);
	jetresolved->Branch("resolved_Et", "std::vector<Float_t>",  &etrj);
	jetresolved->Branch("resolved_phi", "std::vector<Float_t>",  &phirj);
	jetresolved->Branch("resolved_eta", "std::vector<Float_t>",  &etarj);
	jetresolved->Branch("resolved_njets", "std::vector<Float_t>",  &nrjets);
	//THIN JETS:
	jetthin->Branch("jetthin_px", "std::vector<Float_t>",  &qpxj);
	jetthin->Branch("jetthin_py", "std::vector<Float_t>",  &qpyj);
	jetthin->Branch("jetthin_pz", "std::vector<Float_t>",  &qpzj);
	jetthin->Branch("jetthin_E", "std::vector<Float_t>",  &qej);
	jetthin->Branch("jetthin_mass", "std::vector<Float_t>",  &qmj);
	jetthin->Branch("jetthin_Et", "std::vector<Float_t>",  &qetj);
	jetthin->Branch("jetthin_phi", "std::vector<Float_t>",  &qphij);
	jetthin->Branch("jetthin_eta", "std::vector<Float_t>",  &qetaj);
	jetthin->Branch("jetthin_njets", "std::vector<Float_t>",  &qnjets);
	//THICK JETS:
	jetthick->Branch("jetthick_px", "std::vector<Float_t>",  &gpxj);
	jetthick->Branch("jetthick_py", "std::vector<Float_t>",  &gpyj);
	jetthick->Branch("jetthick_pz", "std::vector<Float_t>",  &gpzj);
	jetthick->Branch("jetthick_E", "std::vector<Float_t>",  &gej);
	jetthick->Branch("jetthick_mass", "std::vector<Float_t>",  &gmj);
	jetthick->Branch("jetthick_Et", "std::vector<Float_t>",  &getj);
	jetthick->Branch("jetthick_phi", "std::vector<Float_t>",  &gphij);
	jetthick->Branch("jetthick_eta", "std::vector<Float_t>",  &getaj);
	jetthick->Branch("jetthick_njets", "std::vector<Float_t>",  &gnjets);

	while (myReader.Next()){
		for (int i = 0; i < eventidC.GetSize(); i++){
			ev++;
			
			for (int x = 0; x < scattering.GetSize(); x++){
				if (scattering[x] == 1) { qq++; chk_qq=true; chk_gg=false; chk_gq=false; }
				if (scattering[x] == 2) { gg++; chk_qq=false; chk_gg=true; chk_gq=false; }
				if (scattering[x] == 3) { gq++; chk_qq=false; chk_gg=false; chk_gq=true; }
			}

			// Combined (Res+Dir) Events
			for (int x = 0; x < pxC.GetSize(); x++){
				input_particles_combined.push_back(fastjet::PseudoJet(pxC[x],pyC[x],pzC[x],eC[x]));
				cs++;
			}

			// Resolved Events
			if (eventidR[i]>0){
				for (int x = 0; x < pxR.GetSize(); x++){
					input_particles_resolved.push_back(fastjet::PseudoJet(pxR[x],pyR[x],pzR[x],eR[x]));
					rs++;
				}
				r++;
			}
				
			// Direct Events
			if (eventidD[i]>0){
				for (int x = 0; x < pxD.GetSize(); x++){
					input_particles_direct.push_back(fastjet::PseudoJet(pxD[x],pyD[x],pzD[x],eD[x]));
					ds++;
				}
				d++;
			}

			// Beam Electrons
			in_y = inelasticity[i];
			beam_eT = beT[i];
			beam_e = be[i];
			beam_pT = bpT[i];
			// beam_eta = beta[i];
		}
			//Calling Functions for Jets/Subjets :
			combined_jets(input_particles_combined);
			resolved_jets(input_particles_resolved);
			direct_jets(input_particles_direct);
			
			//Clearing Particles for the next run :
			input_particles_combined.clear();
			input_particles_resolved.clear();
			input_particles_direct.clear();
	}
		cout << "\n=======================================================================================\n" << endl;
		cout << "Number of Combined Events Read : " << ev << " (Constituents : " << cs << ")" << endl;
		cout << "Number of Resolved Events : " << r << " (Constituents : " << rs << ")" << endl;
		cout << "Number of Direct Events : " << d << " (Constituents : " << ds << ")" << endl;
		cout << "\n=======================================================================================\n" << endl;
		cout << "#Jets formed (with Combined Jet Constituents : " << cs << ") : " << 	ncombjets << endl;
		cout << "		Minimum Etmin Combined: " << EtMin_C << " GeV" << endl;
		cout << "		R-Combined: " << R_C << endl;
		cout << "		PseudoRapidity Range : " << crapmin << " to " << crapmax << endl;
		cout << "		Mean SubJet Multiplicity for Combined Jets : " << ncombsubjets << endl;
		
		cout << "\n#Jets formed (with Resolved Jet Events, Constituents : " << rs << ") : " << 	nresjets << endl;
		cout << "		Minimum Etmin Resolved: " << EtMin_R << " GeV" << endl;
		cout << "		R-Resolved: " << R_R << endl;
		cout << "		PseudoRapidity Range : " << rrapmin << " to " << rrapmax << endl;
		cout << "		Mean SubJet Multiplicity for Resolved Jets : " << nressubjets << endl;
		
		cout << "\n#Jets formed (with Direct Jet Events) Constituents : " << ds << ") : " << 	ndirjets << endl;
		cout << "		Minimum Etmin Direct: " << EtMin_D << " GeV" << endl;
		cout << "		R-Direct: " << R_D << endl;
		cout << "		PseudoRapidity Range : " << drapmin << " to " << drapmax << endl;
		cout << "		Mean SubJet Multiplicity for Direct Jets : " << ndirsubjets << endl;
		cout << "\n=======================================================================================\n" << endl;
		cout << fixed << "For Thick Jets:" << endl;
		double t1 = thkqq+thkgg+thkgq;
		cout << fixed << "	QQ Scattering: " << thkqq << "	|	Fraction: " << fixed << (thkqq/t1)*100.0 << endl;
		cout << fixed << "	GG Scattering: " << thkgg << "	|	Fraction: " << fixed << (thkgg/t1)*100.0 << endl;
		cout << fixed << "	GQ Scattering: " << thkgq << "	|	Fraction: " << fixed << (thkgq/t1)*100.0 << endl;
		cout << fixed << "For Thin Jets:" << endl;
		double t2 = thnqq+thngg+thngq;
		cout << fixed << "	QQ Scattering: " << thnqq << "	|	Fraction: " << fixed << (thnqq/(double)t2)*100.0 << endl;
		cout << fixed << "	GG Scattering: " << thngg << "	|	Fraction: " << fixed << (thngg/(double)t2)*100.0 << endl;
		cout << fixed << "	GQ Scattering: " << thngq << "	|	Fraction: " << fixed << (thngq/(double)t2)*100.0 << endl;

		double rho_C_mean = average(rho_C);
		double rho_Q_mean = average(rho_Q);
		double rho_G_mean = average(rho_G);

		double psi_C_mean = average(psi_C); double psi_R_mean = average(psi_R); double psi_D_mean = average(psi_D);
		double psi_Q_mean = average(psi_Q); double psi_G_mean = average(psi_G);

		double etaq_mean = average(etaq);
		double etag_mean = average(etag);
		double dijet1_mean = average(dijet1);
		double dijet2_mean = average(dijet2);
		double jet1_mean = average(jet1);

		cout << "\n=======================================================================================\n" << endl;
		cout << "Differential Shape for r(" << rr_diff << ") EtJetMin(" <<  EtMin_C << ")" << endl;
		cout << fixed << "Combined Jet : " << rho_C_mean << endl;
		cout << fixed << "Quark Jet : " << rho_Q_mean << endl;
		cout << fixed << "Gluon Jet : " << rho_G_mean << endl;
		cout << "\n=======================================================================================\n" << endl;
		cout << "Integrated Shape for r(" << rr_int << ") EtJetMin(" <<  EtMin_C << ")" << endl;
		cout << fixed << "Combined Jet : " << psi_C_mean << endl;
		cout << fixed << "Resolved Jet : " << psi_R_mean << endl;
		cout << fixed << "Direct Jet : " << psi_D_mean << endl;
		cout << "-------------------------------------------------------------" << endl;
		cout << fixed << "Quark Jet : " << psi_Q_mean << endl;
		cout << fixed << "Gluon Jet : " << psi_G_mean << endl;
		cout << "\n=======================================================================================\n" << endl;
		cout << "1 Jet Events : " << count1jet << endl;
		cout << "	with mean ET1 : " << jet1_mean << endl;
		cout << "2 Jet Events : " << count2jet << endl;
		cout << "	with mean ET1 and ET2 : " << dijet1_mean << "	" << dijet2_mean << endl;
		cout << "3+Jet Events : " << count3jet << endl;
		cout << "\n=======================================================================================\n" << endl;
		cout << zz << endl;
		cout << yy << endl;
		cout << etaq_mean <<"	"<<etag_mean;
		
		datajets->Write();
		datajets->Close();
		delete datajets;
}

//---------------------------------------COMBINED JETS-------------------------------------------
void combined_jets (const vector<fastjet::PseudoJet> & input_particles_combined){
  
	fastjet::JetDefinition jet_def(kt_algorithm,1.0);
	fastjet::ClusterSequence clust_seq(input_particles_combined, jet_def);
	// vector<PseudoJet> inclusive_jets_combined = sorted_by_pt(clust_seq.inclusive_jets());
	Selector jet_selector = SelectorEtMin(EtMin_C) && SelectorRapRange(crapmin, crapmax);
	vector<PseudoJet> inclusive_jets_combined = sorted_by_pt(jet_selector(clust_seq.inclusive_jets()));
	temp=inclusive_jets_combined;

	if(inclusive_jets_combined.size()!=0){
		// Number of Jets
		ncjets.push_back(inclusive_jets_combined.size());
		ncJets->Fill(inclusive_jets_combined.size());
		ncombjets+=inclusive_jets_combined.size();

		// Jet shape
		rho_C.push_back(jet_shape_diff(temp));					// Diff
		psi_C.push_back(jet_shape_int(temp, rr_int));		// Int

		// if (jet_shape_int(temp, 0.3) > 0.8){
		// 	rho_Q.push_back(jet_shape_diff(inclusive_jets_combined));
		// 	psi_Q.push_back(jet_shape_int(temp, rr_int));
		// }
		// if (jet_shape_int(temp, 0.3) < 0.6){
		// 	rho_G.push_back(jet_shape_diff(inclusive_jets_combined));
		// 	psi_G.push_back(jet_shape_int(temp, rr_int));
		// }

		TVectorD lambda = transverse_sphericity(inclusive_jets_combined);
		lambda1 = lambda[0];	lambda2 = lambda[1];	ST = abs((2*lambda2)/(lambda1+lambda2));
		STplot->Fill(ST);
		double STCheck = floorf(ST * 10e7) / 10e7;

		if(STCheck==0){
			singlejetshapeC->Fill(jet_shape_int(inclusive_jets_combined, 1.0));
			zz++;

			if (jet_shape_int(temp, 0.3) > 0.8){
				rho_Q.push_back(jet_shape_diff(inclusive_jets_combined));
				psi_Q.push_back(jet_shape_int(temp, rr_int));
			}
			if (jet_shape_int(temp, 0.3) < 0.6){
				rho_G.push_back(jet_shape_diff(inclusive_jets_combined));
				psi_G.push_back(jet_shape_int(temp, rr_int));
			}
		}
		
		if (inclusive_jets_combined.size() == 1){
			// cout << std::fixed << std::setprecision(4);
			// cout << "\nEigenvalue 1 : " << lambda1 << " && Eigenvalue 2 : " << lambda2 << endl;
			// cout << "Lambda1+Lambda2 : " << lambda1+lambda2 << endl;
			// zz++;
			// cout << "Transverse Sphericity : " << lround(abs(ST)) << "	" << zz << endl;
			// cout << "\nBeam Et, E, pT, eta : " << beam_eT << " " << beam_e << " " << beam_pT << " " << beam_eta << endl;
			// cout << "1-Jet Et, E, pT, eta : " << temp[0].Et() << " " << temp[0].e() << " " << temp[0].pt() << " " << temp[0].eta() << endl;
			count1jet++;
			if ((inclusive_jets_combined[0].Et()>0))
			jet1.push_back(inclusive_jets_combined[0].Et());
			STplot1->Fill(ST);
			yy++;
		}
		
		if (inclusive_jets_combined.size() == 2){
			cout << std::fixed << std::setprecision(4);
			// cout << "\nEigenvalue 1 : " << lambda1 << " && Eigenvalue 2 : " << lambda2 << endl;
			// cout << "Lambda1+Lambda2 : " << lambda1+lambda2 << endl;
			// cout << "Transverse Sphericity : " << ST << endl;
			// cout << inclusive_jets_combined[0].eta() << "	" << inclusive_jets_combined[1].eta() << endl;
			dijets(inclusive_jets_combined);
			count2jet++;
			if (inclusive_jets_combined[1].Et()>0){
				dijet1.push_back(inclusive_jets_combined[0].Et());
				dijet2.push_back(inclusive_jets_combined[1].Et());
			}
			STlambda1->Fill(lambda1);
			STlambda2->Fill(lambda2);
			STplot2->Fill(ST);
			// x_obs(inclusive_jets_combined);
		}
		if (inclusive_jets_combined.size() > 2){
			count3jet++;
			STplot3->Fill(ST);
		}
		// if (inclusive_jets_combined.size() == 4){
		// 	STplot4->Fill(ST);
		// }
	}

	for (int i = 0; i < inclusive_jets_combined.size(); i++){
		
		// if(inclusive_jets_combined[i].Et()<14) 
		// 	continue;

		pxcj.push_back(inclusive_jets_combined[i].px());
		pycj.push_back(inclusive_jets_combined[i].py());
		pzcj.push_back(inclusive_jets_combined[i].pz());
		mcj.push_back(inclusive_jets_combined[i].m());
		ecj.push_back(inclusive_jets_combined[i].e());
		etcj.push_back(inclusive_jets_combined[i].Et());
		ptcj.push_back(inclusive_jets_combined[i].perp());
		phicj.push_back(inclusive_jets_combined[i].phi());
		etacj.push_back(inclusive_jets_combined[i].eta());
		etaphi_planeC->Fill(inclusive_jets_combined[i].eta(), inclusive_jets_combined[i].phi());
		if (jet_shape_int(temp, 0.3) > 0.8)	etaq.push_back(inclusive_jets_combined[i].eta());
		if (jet_shape_int(temp, 0.3) < 0.6)	etag.push_back(inclusive_jets_combined[i].eta());

		// int n_constituents = inclusive_jets_combined[i].constituents().size();
		// if(n_constituents > 20)
		// print_jets(inclusive_jets_combined,1);

		//--------------------------------------SUBJETS for COMBINED EVENTS--------------------------------------
		dcut_C =  ycut_C*pow(inclusive_jets_combined[i].Et(),2);
		vector<PseudoJet> subJets_comb = inclusive_jets_combined[i].exclusive_subjets(dcut_C);

		if (subJets_comb.size() != 0){
			ncsubJets->Fill(subJets_comb.size());
			ncombsubjets+=subJets_comb.size();
			n_C.push_back(subJets_comb.size());

			// Thin Jets gives Quark initiated jets
			if (jet_shape_int(temp, 0.3) > 0.8){
				nqthinsubJets->Fill(subJets_comb.size());
				// etaphi_planethn->Fill(inclusive_jets_combined[i].eta(), inclusive_jets_combined[i].phi());
				// n_Q.push_back(subJets_comb.size());
				// rho_Q.push_back(jet_shape_diff(inclusive_jets_combined));
				// psi_Q.push_back(jet_shape_int(inclusive_jets_combined, 1.0));
				if(chk_qq==true) thnqq++;
				if(chk_gg==true) thngg++;
				if(chk_gq==true) thngq++;
				for (int i = 0; i < subJets_comb.size(); i++){
					qpxj.push_back(subJets_comb[i].px());
					qpyj.push_back(subJets_comb[i].py());
					qpzj.push_back(subJets_comb[i].pz());
					qmj.push_back(subJets_comb[i].m());
					qej.push_back(subJets_comb[i].e());
					qetj.push_back(subJets_comb[i].Et());
					qphij.push_back(subJets_comb[i].phi());
					qetaj.push_back(subJets_comb[i].eta());
				}
			}

			// Thick Jets gives Gluon initiated jets
			if (jet_shape_int(temp, 0.3) < 0.6){
				ngthicksubJets->Fill(subJets_comb.size());
				// etaphi_planethk->Fill(inclusive_jets_combined[i].eta(), inclusive_jets_combined[i].phi());
				// n_G.push_back(subJets_comb.size());
				// rho_G.push_back(jet_shape_diff(inclusive_jets_combined));
				// psi_G.push_back(jet_shape_int(inclusive_jets_combined, 1.0));
				if(chk_qq==true) thkqq++;
				if(chk_gg==true) thkgg++;
				if(chk_gq==true) thkgq++;
				for (int i = 0; i < subJets_comb.size(); i++){
					gpxj.push_back(subJets_comb[i].px());
					gpyj.push_back(subJets_comb[i].py());
					gpzj.push_back(subJets_comb[i].pz());
					gmj.push_back(subJets_comb[i].m());
					gej.push_back(subJets_comb[i].e());
					getj.push_back(subJets_comb[i].Et());
					gphij.push_back(subJets_comb[i].phi());
					getaj.push_back(subJets_comb[i].eta());
				}
			}
		}
	}
	
	jetcombined->Fill(); 
	pxcj.clear(); pycj.clear(); pzcj.clear(); mcj.clear(); 
	ecj.clear(); etcj.clear(); ptcj.clear(); phicj.clear(); etacj.clear(); ncjets.clear();
	jetthin->Fill();
  qpxj.clear(); qpyj.clear(); qpzj.clear(); qmj.clear(); 
	qej.clear(); qetj.clear(); qphij.clear(); qetaj.clear(); qnjets.clear();
	mean_js_diff.clear();
	jetthick->Fill();
  gpxj.clear(); gpyj.clear(); gpzj.clear(); gmj.clear(); 
	gej.clear(); getj.clear(); gphij.clear(); getaj.clear(); gnjets.clear();
}

//---------------------------------------RESOLVED JETS-------------------------------------------
void resolved_jets (const vector<fastjet::PseudoJet> & input_particles_resolved){
  
  fastjet::JetDefinition jet_def(kt_algorithm,R_R);
	fastjet::ClusterSequence clust_seq(input_particles_resolved, jet_def);
	Selector jet_selector = SelectorEtMin(EtMin_R) && SelectorRapRange(rrapmin, rrapmax); 
	vector<PseudoJet> inclusive_jets_resolved = sorted_by_pt(jet_selector(clust_seq.inclusive_jets()));
	
	if(inclusive_jets_resolved.size()!=0){
		nrjets.push_back(inclusive_jets_resolved.size());
		nrJets->Fill(inclusive_jets_resolved.size());
		nresjets+=inclusive_jets_resolved.size();

		// Jet Shape
		psi_R.push_back(jet_shape_int(inclusive_jets_resolved, rr_int));
	}

	if(inclusive_jets_resolved.size()==2){
		rxobs->Fill(x_obs(inclusive_jets_resolved));
	}

	for (int i = 0; i < inclusive_jets_resolved.size(); i++){
		pxrj.push_back(inclusive_jets_resolved[i].px());
		pyrj.push_back(inclusive_jets_resolved[i].py());
		pzrj.push_back(inclusive_jets_resolved[i].pz());
		erj.push_back(inclusive_jets_resolved[i].e());
		etrj.push_back(inclusive_jets_resolved[i].Et());
		phirj.push_back(inclusive_jets_resolved[i].phi());
		etarj.push_back(inclusive_jets_resolved[i].eta());

		//--------------------------------------SUBJETS for RESOLVED EVENTS--------------------------------------
		dcut_R =  ycut_R*pow(inclusive_jets_resolved[i].Et(),2);
		vector<PseudoJet> subJets_res = inclusive_jets_resolved[i].exclusive_subjets(dcut_R);
		
		if (subJets_res.size() != 0){
			//nsubcjets.push_back(subJets_quark.size());
			nrsubJets->Fill(subJets_res.size());
			nressubjets+=subJets_res.size();
			n_R.push_back(subJets_res.size());
		}
	}
	
	jetresolved->Fill();
	pxrj.clear(); pyrj.clear(); pzrj.clear(); 
	erj.clear(); etrj.clear(); phirj.clear(); etarj.clear(); nrjets.clear();
}

//---------------------------------------DIRECT JETS---------------------------------------------
void direct_jets (const vector<fastjet::PseudoJet> & input_particles_direct){
  
  fastjet::JetDefinition jet_def(kt_algorithm,R_D);
	fastjet::ClusterSequence clust_seq(input_particles_direct, jet_def);
	Selector jet_selector = SelectorEtMin(EtMin_D) && SelectorRapRange(drapmin, drapmax); 
	vector<PseudoJet> inclusive_jets_direct = sorted_by_pt(jet_selector(clust_seq.inclusive_jets()));
	
	if(inclusive_jets_direct.size()!=0){
		ndjets.push_back(inclusive_jets_direct.size());
		ndJets->Fill(inclusive_jets_direct.size());
		ndirjets+=inclusive_jets_direct.size();

		// Jet Shape
		psi_D.push_back(jet_shape_int(inclusive_jets_direct, rr_int));
	}
	if(inclusive_jets_direct.size()==2){
		dxobs->Fill(x_obs(inclusive_jets_direct));
	}
		
	for (int i = 0; i < inclusive_jets_direct.size(); i++){
		pxdj.push_back(inclusive_jets_direct[i].px());
		pydj.push_back(inclusive_jets_direct[i].py());
		pzdj.push_back(inclusive_jets_direct[i].pz());
		edj.push_back(inclusive_jets_direct[i].e());
		etdj.push_back(inclusive_jets_direct[i].Et());
		phidj.push_back(inclusive_jets_direct[i].phi());
		etadj.push_back(inclusive_jets_direct[i].eta());

		//--------------------------------------SUBJETS for DIRECT EVENTS--------------------------------------
		dcut_D =  ycut_D*pow(inclusive_jets_direct[i].Et(),2);
		vector<PseudoJet> subJets_dir = inclusive_jets_direct[i].exclusive_subjets(dcut_D);
		
		if (subJets_dir.size() != 0){
			//nsubcjets.push_back(subJets_quark.size());
			ndsubJets->Fill(subJets_dir.size());
			ndirsubjets+=subJets_dir.size();
			n_D.push_back(subJets_dir.size());
		}
	}
		
	jetdirect->Fill(); gStyle->SetLineWidth(2);
	pxdj.clear(); pydj.clear(); pzdj.clear(); 
	edj.clear(); etdj.clear(); phidj.clear(); etadj.clear(); ndjets.clear();
}

//------------------------------------------Dijets-----------------------------------------------
void dijets (const vector<fastjet::PseudoJet> & inclusive_dijets){
	// if (inclusive_dijets[1].Et()>12)
	// cout << 1;
}

//------------------------------------------X-OBS------------------------------------------------
double x_obs(const vector<fastjet::PseudoJet> & inclusive_jets_combined){
	double xobs=0., pt1=0., eta1=0., pt2=0., eta2=0.;
	eta1 = inclusive_jets_combined[0].eta();
	pt1 = inclusive_jets_combined[0].pt();
	eta2 = inclusive_jets_combined[1].eta();
	pt2 = inclusive_jets_combined[1].pt();
	double num = (pt1*exp(-eta1))+(pt2*exp(-eta2));
	double denom = 2*beam_e;
	xobs = num/denom;
	// cout << xobs << endl;
	return xobs;
}

//---------------------------------------TRANSVERSE SPHERICITY-----------------------------------
TVectorD transverse_sphericity(const vector<fastjet::PseudoJet> & inclusive_jets_combined){
	TMatrixD Si(2,2), Siall(2,2);
	Si.Zero(); Siall.Zero();
	double pxi=0,pyi=0,pTi=0,pTall=0;
	for (int i = 0; i < inclusive_jets_combined.size(); i++){
		pxi=inclusive_jets_combined[i].px();
		pyi=inclusive_jets_combined[i].py();
		pTi=inclusive_jets_combined[i].pt();
		pTall+=pTi;
		Si[0][0]+= pow(pxi,2)/pTi;
		Si[0][1]+= (pxi*pyi)/pTi;
		Si[1][0]+= (pxi*pyi)/pTi;
		Si[1][1]+= pow(pyi,2)/pTi;
	}
	Siall[0][0]= Si[0][0]/pTall;
	Siall[0][1]= Si[0][1]/pTall;
	Siall[1][0]= Si[1][0]/pTall;
	Siall[1][1]= Si[1][1]/pTall;
	TVectorD eigenVal;
	const TMatrixD eigenVec = Siall.EigenVectors(eigenVal);
	return eigenVal;
}

//------------------------------------Jet shape (differential)-----------------------------------
double jet_shape_diff (const vector<fastjet::PseudoJet> & inclusive_jets){
	
	double r=rr_diff;								// select r for R_jet = 1
	double delta_r = 0.10;			 		// R2-R1
	double R1 = r-(delta_r/2);			// inner radius
	double R2 = r+(delta_r/2);			// outer radius
	double select_const_r;					// stores the contituent radius
	double select_const_Et = 0;			// stores the contituent Et
	double select_sum = 0;					// summation
	double rho = 0;									// differential 

	// Loop over all inclusive jets
	for (int i = 0; i < inclusive_jets.size(); i++)
	{
		// store the i-jet constituents
		vector <fastjet::PseudoJet> constituents = inclusive_jets[i].constituents();

		// loop over the constituents and select the ones which lie between R1 and R2
		for (unsigned int j = 0; j < constituents.size(); j++){

			select_const_r= pow((inclusive_jets[i].phi()-constituents[j].phi()),2) + pow((inclusive_jets[i].eta()-constituents[j].eta()),2);
			select_const_r= pow(select_const_r,0.5);

			// add Et of constituents in the doughnut of R2-R1
			if( select_const_r > R1 && select_const_r < R2){
				select_const_Et = select_const_Et + constituents[j].Et();
			}
			// divide the total selected jet contituents' Et with the total jet Et 
			select_sum = select_const_Et/inclusive_jets[i].Et();
		}
		// calculate rho, ie the jet shape for the inclusive jets
		rho = rho + (select_sum/delta_r);
		mean_js_diff.push_back(rho);

		// reset variable for another r run 
		select_const_Et=0;
	}
	return rho;
}

//------------------------------------Jet shape (integrated)-------------------------------------
double jet_shape_int (const vector<fastjet::PseudoJet> & inclusive_jets, double rr_int){

	// stores the constituent r, Et
	double select_const_r=0.0, select_const_Et=0.0;
	double select_sum=0.0;
	double psi = 0.0, t1=0.0, t2=0.0, r=rr_int;
	int Njets = inclusive_jets.size();

	for (int i = 0; i < Njets; i++){
		// store the i-jet constituents
		vector <fastjet::PseudoJet> constituents = inclusive_jets[i].constituents();

		// loop over the constituents and select the ones which lie < rr_inclusive 
		for (unsigned int j = 0; j < constituents.size(); j++){

			t1 = pow((inclusive_jets[i].phi()-constituents[j].phi()),2);
			t2 = pow((inclusive_jets[i].rapidity()-constituents[j].rapidity()),2);
			select_const_r = pow((t1 + t2),0.50);

			// add Et of constituents < r_inclusive
			if(select_const_r < r)	select_const_Et = select_const_Et + constituents[j].Et();
		}

		// divide the total selected jet contituents' Et with the total jet Et 
		select_sum = select_const_Et/inclusive_jets[i].Et();

		// calculate psi, ie the jet shape for the inclusive jet sample 
		psi = psi + (select_sum);
		// mean_js_incl.push_back(psi);

		// reset variable for another r run 
		select_const_Et=0.;
		select_sum=0.;
	}
	//cout << " " << psi << endl; 
	return (psi/Njets);
}

//------------------------------------PRINT JETS and SUBJETS-------------------------------------
void print_jets (const vector<fastjet::PseudoJet> & jets, int flag) {
  // sort jets into increasing pt
  vector<fastjet::PseudoJet> sorted_jets = sorted_by_pt(jets);  
	cout << "Event :: " << ev << endl;
  // label the columns
  if (flag == 1)
  printf("%5s %15s %15s %15s %15s %15s\n","Jet", "rapidity", "phi", "pt", "Et", "constituents");
	if (flag == 2)
	printf("%5s %15s %15s %15s %15s %15s\n","sJet", "rapidity", "phi", "pt", "Et", "constituents");
  
  // print out the details for each jet
  for (unsigned int i = 0; i < sorted_jets.size(); i++) {
    int n_constituents = sorted_jets[i].constituents().size();
		vector<PseudoJet> constituents = sorted_jets[i].constituents();

		// print Jets (rapidity, phi_std (−π to π), pT and eT)
    printf("%5u %15.8f %15.8f %15.8f %15.8f %8u\n",
	   i, sorted_jets[i].rap(), sorted_jets[i].phi(),
	   sorted_jets[i].pt(), sorted_jets[i].Et(), n_constituents);
		cout << "		Constituents :" << endl;
		// print the constituents (rapidity, phi_std (−π to π), pT and eT)
		for (unsigned int j = 0; j < n_constituents; j++){

				cout<< "						" << constituents[j].rap() << "				"
						<< constituents[j].phi() << "				"
						<< constituents[j].pt() << "					"
						<< constituents[j].Et() << endl;
		}
		cout << "\n";
  }
}

// Auxiliary: To calculate the average of a std::vector<double>
double average(const vector<double> &v){
    if(v.empty())	return 0;
    double avg = accumulate( v.begin(), v.end(), 0.0)/v.size();
    return avg;
}
