/*
    File: chromoCell_main.cpp
    Function: simulate chromosome dynamics in the context of nuclear envelope and nucleolus
    Model: chromoCell
    Created: 7 December, 2017 (XF)
    ---------------------------------------------
    format: executible argv[1]
    argv[1] is the file name of data
    ---------------------------------------------
    Notes:
    Week of 2018.01.08
        (1) DONE: output a distance matrix from simulation
    Week of 2018.01.15
        (1) DONE: construct a random chain initial configuration (pymol doesn't recognize properly PDB file with >9999 atoms)
            --> summary: potential problem is initial self-entanglement!
        (2) DONE: introduce core functions for Brownian dynamics
            --> summary: tension, attraction, repulsion, entropic
        (3) DONE: introduce condensin tension
            --> summary: good to see condensin binding sites cluster into one without dissociation or maximum number restraint;
                         now testing these restraints
        (4) NOT DOING: Hilbert curve filling

    Week of 2018.01.22 & Week of 2018.01.29
        (1) DONE: change to (overdamped) Langevin equation
            --> thought: proper parameters for random noise
        (2) NOT DONE: apply cylindrical restraint for random walk and simulation
            --> thought: may or may not change VoxMap dimension
        (3) DONE: implement loop extrusion
            --> [Y] thought: color the loop region using a different color (different protein type)
            --> [Y] thought: diffusion capture of binders instead of condensin binding sites?
        (4) NOT DONE: allow dynamic pairing in diffusion capture

    Week of 2018.02.05 & Week of 2018.02.12
        (1) DONE: combine loop extrusion and diffusion capture
            --> [Y] thought: stop binder walking when two binders meet!
        (2) DONE: openmp implementation
            --> summary: repulsion forces
        (3) DONE: introduce static road blocks -- random distribution
        (4) DONE: introduec static road blocks -- based on chipseq cohesin
        (5) DONE: modify initial configuration to bined random walk
            --> thought: fix centromere coordinate location & grow chain step by step

    Week of 2018.02.19 & Week of 2018.02.26
        (1) DONE: initial configuration of random chain binned tether
            --> thought: fix both centromere and telomere coordinates
            --> thought: number of beads per bin scaled by volume?
        (2) DONE: put fluorophores on the chromosome
            --> thought: how many possible configurations based on Haering paper?
        (3) DONE: experiment: diffusion capture only, loop extrusion only, both
        (4) DONE: generate multiple relaxed configurations (interphase) and plot frequency over genomic separation

    Week of 2018.03.05 & Week of 2018.03.12
        (1) DONE: initial configuration
            --> plan: binned chain in a cylinder; scaled by volume
        (2) TO DO: bubble plot & HiC-map like figure
        (3) DONE: make MSD plot
            --> summary: slope = 0.5 in log-log plot makes sense, but the absolute MSD is not good; need smaller viscosity?

    Week of 2018.03.19 & Week of 2018.03.26
        (1) TO DO: correct the mistake in cylindrical initial configuration
        (2) TO DO: tune viscosity for MSD

        (?) TO DO: initialize chain configuration "segment" by "segment" according to fluorophores
        (?) DONE: save locations of binders over time
        (?) DONE: save parameter list
        (?) DONE: add PID in output file name (in preparation for CAMP simulations)
        (?) DONE: plot binder lifetime (residence time), loop size distribution, etc
        (?) DONE: plot frequency over genomic separation
        (?) DONE: plot end-to-end distance over time
        (?) TO DO: plot distribution of binder locations (validation against chipseq condensin)
*/

#include "initConfig.hpp"
#include "dynamics.hpp"


void writeParameterInfo(void);  // write model parameters into file
int main(int argc, char** argv)
{
    // write all parameter into file
    writeParameterInfo();

    // Set up random number generators
    int seed = chrono::system_clock::now().time_since_epoch().count();
    mt19937 generator (seed);
    uniform_real_distribution<double> uniform01(0.0, 1.0);
    //normal_distribution<double> normal0(0.0, sqrt(F_ENTROPIC_VARIANCE));
    normal_distribution<double> normal0(0.0, sqrt(DT));

    // initialize variables
    unordered_map<string, double> energy;
    energy["total"] = 0.;
    energy["kinetic"] = 0.;
    energy["potential"] = 0.;

    // initialize voxMap -- if using VoxMap method
    //const int nSubdiv = 3;
    //const int dimX = pow(2, nSubdiv);
    const int dimX = DIM_VOXMAP;
    const iVec dim = {dimX, dimX, dimX};
    unordered_map<int, vector<int>> mapBeadIds;
    unordered_map<int, vector<int>> mapBeadIdsHalo;
    VoxMap voxMap = VoxMap(dim, mapBeadIds, mapBeadIdsHalo);

    // initialize octree -- if using Octree method
    /*
    int depth = 0;
    vector<Leaf> vecLeaf;
    Octree octree = Octree(depth, vecLeaf);
    */

    // initialize vectors
    vector<Chromosome> vecChromosome;           // collect all Chromosome objects
    vector<Binder>     vecBinder;               // collect all Binder objects
    vector<Node>       vecNode;                 // collect all Node objects (Beads, Vertics)
    vector<int>        vecFluoroSite;           // collect all fluorophore sites
    vector<int>        vecCondensinSite;        // collect all condensin sites              // change to Domain object some time?
    vector<int>        vecBinderBlockSite;      // collect all sites that block binders     // change to Domain object some time?
    vector<int>        vecBinderBlockSiteMitosis;      // collect all sites that block binders     // change to Domain object some time?
   
    //double **matChromoContFreq;                 // contact frequency matrix for intra- and inter-chromosome interactions
    vector<iPair> vecBeadPairInteract;          //
    vector<double> vecBeadPairInteractFreq;
    //
    vector<int> vecLoadedBinders;
    vector<int> vecCondensinSiteEmpty;
    vector<int> vecBinderIds2;
    vector<int> vecUnloadedBinders;
    vector<int> vecUnLoadedBinders;
    vector<int> vecCondensinSiteEmpty2;
    vector<int> vecDomains_DC;
    vector<int> vecDomains_DCMitosis;
    //


/////////////////////////////////////////////////////////////////////////
	vector<iPair> capturer;       // monitor DC sites


// ----------------------------------- 03-08-2018 -------------------------------------
	vector<int> hit_valence;
	
/////////////////////////////////////////////////////////////////////////

    // initialize chromosomes (beads)
    if (SIM_TYPE == "RECONSTRUCT")
        initChromosome(vecChromosome, vecNode, vecBeadPairInteract, vecBeadPairInteractFreq, argc, argv);  // pass the reference to vectors
    else if (SIM_TYPE == "SIMULATION")
        initChromosomeBrownianSimulation(vecChromosome, vecNode);

    // initialize fluorophores & condensin sites (binders starting sites) & block sites (static road blocks)
    // ... need to group into functions
    if (SIM_TYPE == "SIMULATION")
    {   

        cout << ">> initializing fluorophore sites (Petrova2013 paper figure 2B-C)" << endl;		
        vector<int> vecDistFluoroCentromere;    // distance measured in number of beads
        vecDistFluoroCentromere.push_back(50);  //push_back(100000/(BP_SEP_NUCLEOSOME*COARSE_FACTOR));
        vecDistFluoroCentromere.push_back(350); 	//push_back(700000/(BP_SEP_NUCLEOSOME*COARSE_FACTOR));
        vecDistFluoroCentromere.push_back(600); 	//push_back(1200000/(BP_SEP_NUCLEOSOME*COARSE_FACTOR));
        vecDistFluoroCentromere.push_back(850);	//push_back(1700000/(BP_SEP_NUCLEOSOME*COARSE_FACTOR));
        vecDistFluoroCentromere.push_back(1100);	//push_back(2200000/(BP_SEP_NUCLEOSOME*COARSE_FACTOR));
        vecDistFluoroCentromere.push_back(900);	//push_back(1803687/(BP_SEP_NUCLEOSOME*COARSE_FACTOR)); // !!!!!!!!  Yasu's fluorophore chrI 1.95Mb, CEN1=3 753 687 -> 3 789 421, distance_from_CEN1 = 3753687 - 1950000 = 1803687
        
        cout << vecDistFluoroCentromere[1] << endl;
        cout << vecDistFluoroCentromere[2] << endl;
        cout << vecDistFluoroCentromere[3] << endl;
        cout << vecDistFluoroCentromere[4] << endl;
        cout << vecDistFluoroCentromere[5] << endl;
        cout << vecDistFluoroCentromere[6] << endl;        
        
        for (int i = 0; i < vecNode.size(); i ++)
        {
            int dist2cen = vecNode.size() - i;
            if (i == 0 || i == vecNode.size()-1 \
                || find(vecDistFluoroCentromere.begin(), vecDistFluoroCentromere.end(), dist2cen) != vecDistFluoroCentromere.end())
            {
                vecFluoroSite.push_back(i);
                cout << "... beadId = " << i << " is labeled by fluorophore ..." << endl;
            }
        }
        cout << "    done!" << endl;
        
        if (BINDER_POSITION == "uniform")
        {

			cout << ">> initializing **uniform** condensin binding sites (start sites for condensin translocation) ..." << endl;
			for (int i = 0; i < vecNode.size(); i ++)
			{
				if (i != 0 && i % BINDER_SPACING == 0)
    
                //if (i != 0 && i % ((int) ((double) NUM_BEAD/NUM_CONDENSIN_PEAKS)) == 0 && vecCondensinSite.size() < NUM_CONDENSIN_PEAKS)
	//			if (i != 0 && i % ((int) ((double) BP_SEP_CONDENSIN_BINDING/(BP_SEP_NUCLEOSOME*10))) == 0)
                                //if (i != 0 && i % ((int) ((double) NUM_BEAD/NUM_CONDENSIN_PEAKS)) == 0 && vecCondensinSite.size() < NUM_CONDENSIN_PEAKS)
	//			if (i != 0 && i % ((int) ((double) BP_SEP_CONDENSIN_BINDING/(BP_SEP_NUCLEOSOME*10))) == 0 && vecCondensinSite.size() < NUM_CONDENSIN_PEAKS)

				{
					vecCondensinSite.push_back(i);
					cout << "... beadId = " << i << " is a condensin BINDing site ..." << endl;
				}
			}
			cout << "    done!" << endl;
	    }
        
        
        //----------------------16-07-2019- CHIP-corrected vector -------------------
        if (BINDER_POSITION == "CHIP-corrected")
        {

			cout << ">> initializing **uniform** condensin binding sites (start sites for condensin translocation) ..." << endl;
			vector<int> vecCHIP_seq_condensin_spombe{46,96,102,103,123,126,146,170,183,223,228,229,259,264,266,267,268,270,279,318,352,359,363,365,366,370,371,449,479,480,484,494,499,501,548,575,592,598,601,609,618,623,643,658,669,690,700,711,725,740,743,756,764,785,786,798,823,834,852,856,867,873,879,886,911,915,927,949,955,959,962,981,987,1054,1067,1107,1108,1116,1121,1124,1125,1146,1150,1165,1166,1183,1211,1217,1220,1221,1255,1274,1276,1281,1289,1313,1334,1339,1340,1342,1349,1353,1356,1367,1373,1378,1396,1400,1401,1403,1404,1405,1422,1424,1425,1427,1440,1463,1488,1504,1510,1524,1535,1543,1582,1583,1620,1621,1645,1649,1651,1655,1657,1693,1698,1722,1730,1774,1780,1792,1796,1809,1817,1828,1834,1848,1850,1851,1855,1862,1867,1868,1873,1874,1876,1877,1878,1880}; 
			
			for (int i = 0; i < vecCHIP_seq_condensin_spombe.size(); i ++)
			{
				{
					vecCondensinSite.push_back(vecCHIP_seq_condensin_spombe[i]-1);
					cout << "... beadId = " << vecCHIP_seq_condensin_spombe[i]-1 << " is a condensin BINDing site ..." << endl;
				}
			}
			cout << "    done!" << endl;
	    }        
        
        
                    // ---------------- 08-08-2018 ChIP-seq loading -------------------
            
        if (BINDER_POSITION == "CHIP-seq")
        {
            cout << ">> initializing **" << BINDER_POSITION << "** ChIP-seq binder loading sites ..." << endl;
            string iFilename_chipseq_condesin;
            vector<int> vecCondensinBP_CG;
            //vector<int> vecCondensinSiteOccupied;
            if (argc == 1)
            {
                cout << "\nPlease input cohesin chipseq filename (without suffix) as argument! \n";
                exit(0);
            }
            if (argc == 2)
            {
                iFilename_chipseq_condesin = argv[1]+SUFFIX_chipseq;
                cout << ">> reading file " << iFilename_chipseq_condesin << " ... "<< endl;
                ifstream file_chipseq_condesin (iFilename_chipseq_condesin);
                string line_chipseq_condesin;
                
                if (file_chipseq_condesin.is_open())
                {
                    int cnt = 0;
                    getline (file_chipseq_condesin, line_chipseq_condesin); // skip the first line
                    while (getline (file_chipseq_condesin, line_chipseq_condesin))
                    {
                        stringstream ss (line_chipseq_condesin);
                        vector<string> items;
                        string buf;
                        while (ss >> buf)
                            items.push_back(buf);
                        //int chromoId, loc_bp, loc_bead;
                        int loc_bp_start, loc_bp_end, loc_bp, loc_bead, loc_diff;
                        //chromoId = stoi(items[0]);
                        loc_bp_start = stoi(items[0]);
                        loc_bp_end = stoi(items[1]);
                        loc_diff   = (loc_bp_end - loc_bp_start)/2;
                        loc_bp = loc_bp_start + loc_diff;  // unit is bp
                        //loc_bead = loc_bp/(N_BP_PER_NUCLEO);
                        loc_bead = loc_bp/(BP_SEP_NUCLEOSOME)*0.1;
                        //vecCondensinBP_CG.push_back(loc_bead);
                        cout << " start " << loc_bp_start << endl;
                        cout << " loc_bp = " << loc_bp << endl;


                            // push_back locations on left arm to vecBinderBlockSite
                        //if (chromoId == 1 && loc_bead < NUM_BEAD)
                        //int NUM_BEAD_BP = NUM_BEAD * COARSE_FACTOR * (N_BP_PER_NUCLEO + BP_SEP_NUCLEOSOME); 
                        
                        if (loc_bead < NUM_BEAD && find(vecCondensinSite.begin(), vecCondensinSite.end(), loc_bead) == vecCondensinSite.end())
                        {
                            vecCondensinSite.push_back(loc_bead);
                           // vecCondensinSiteOccupied.push_back(loc_bead);
                            cnt ++;
                            cout << "... beadId = " << loc_bead << " is a condensin LOADING site ... " << endl;
                        }
                        
                    }
                    cout << "total number : " << cnt << endl;
                }

				for (int i=0; i < vecCondensinSite.size(); i++)
				{	
					cout << vecCondensinSite[i] << " " << endl;
				}
				cout << vecCondensinSite.size() << " CondensinSites" << endl;
				
            }
            cout << "    done!" << endl;
        }    
        
        
        // --------------------------------------08-08-2018 RANDOM distribution pre-defined ------------------------
        
          
        if (BINDER_POSITION == "random1")
        {
            cout << ">> initializing **" << BINDER_POSITION << "** pre-defined random binder loading sites ..." << endl;
            string iFilename_random_distrib;
            vector<int> vecCondensinBP_CG;
            //vector<int> vecCondensinSiteOccupied;
            if (argc == 1)
            {
                cout << "\nPlease input cohesin chipseq filename (without suffix) as argument! \n";
                exit(0);
            }
            if (argc == 2)
            {
                iFilename_random_distrib = argv[1]+SUFFIX_random;
                cout << ">> reading file " << iFilename_random_distrib << " ... "<< endl;
                ifstream file_random_distrib (iFilename_random_distrib);
                string line_random_distrib;
                
                if (file_random_distrib.is_open())
                {
                    int cnt = 0;
                    //getline (file_random_distrib, line_random_distrib); // skip the first line
                    while (getline (file_random_distrib, line_random_distrib))
                    {
                        stringstream ss (line_random_distrib);
                        vector<string> items;
                        string buf;
                        while (ss >> buf)
                            items.push_back(buf);
                        //int chromoId, loc_bp, loc_bead;
                        int loc_bead;
                        //chromoId = stoi(items[0]);
                        loc_bead = stoi(items[0]);
                        cout << " checkpoint1" << endl;
                        //vecCondensinBP_CG.push_back(loc_bead);
                        cout << " loc_bead = " << loc_bead << endl;


                            // push_back locations on left arm to vecBinderBlockSite
                        //if (chromoId == 1 && loc_bead < NUM_BEAD)
                        //int NUM_BEAD_BP = NUM_BEAD * COARSE_FACTOR * (N_BP_PER_NUCLEO + BP_SEP_NUCLEOSOME); 
                        
                        if (loc_bead < NUM_BEAD)
                        {
                            vecCondensinSite.push_back(loc_bead);
                           // vecCondensinSiteOccupied.push_back(loc_bead);
                            cnt ++;
                            cout << "... beadId = " << loc_bead << " is a condensin LOADING site ... " << endl;
                        }
                        
                    }
                    cout << "total number : " << cnt << endl;
                }
/*
				for (int i=0; i < vecCondensinSite.size(); i++)
				{	
					cout << vecCondensinSite[i] << " " << endl;
				}
				cout << vecCondensinSite.size() << " CondensinSites" << endl;
*/				
            }
            cout << "    done!" << endl;
        }         
               
            
        // ----------------------------------------------------------------

        // block sites between condensin sites -- assume 0/1/2 random block sites between condensin sites
        if (STATIC_BINDER_BLOCK_ON == true)
        {
            if (STATIC_BINDER_BLOCK_TYPE == "random")
            {
                cout << ">> initializing **" << STATIC_BINDER_BLOCK_TYPE << "** static binder block sites ..." << endl;
                int cnt = 0;
                for (int j = 0; j < vecCondensinSite.size()-1; j ++)
                {
                    double a  =  rand()/(double) RAND_MAX;
                    int interval = vecCondensinSite[j+1] - vecCondensinSite[j];
                    // default values to give about 111 BLOCKing sites
                    if ( a <= 0.7 )    // no block sites; default 0.7
                        continue;
                    else if (a <= 0.9) // 1 block site;   default 0.9
                    {
                        int b = rand()%interval + vecCondensinSite[j]+1;
                        vecBinderBlockSite.push_back(b);
                        cnt ++;
                        cout << "... beadId = " << b << " is a condensin BLOCKing site ... between ("
                             << vecCondensinSite[j] << ", " << vecCondensinSite[j+1] << ")" << endl;
                    }
                    else                // 2 block sites
                    {
                        int b, c, small, big;
                        b = rand()%interval + vecCondensinSite[j]+1; c = b;
                        while (c == b)
                            c = rand()%interval + vecCondensinSite[j]+1;
                        if (c < b)
                        {
                            small = c; big = b;
                        }
                        else
                        {
                            small = b; big = c;
                        }
                        vecBinderBlockSite.push_back(small);
                        cnt ++;
                        cout << "... beadId = " << small << " is a condensin BLOCKing site ... between ("
                             << vecCondensinSite[j] << ", " << vecCondensinSite[j+1] << ")" << endl;
                        vecBinderBlockSite.push_back(big);
                        cnt ++;
                        cout << "... beadId = " << big << " is a condensin BLOCKing site ... between ("
                             << vecCondensinSite[j] << ", " << vecCondensinSite[j+1] << ")" << endl;
                    }
                }
                cout << "total number : " << cnt << endl;
                cout << "    done!" << endl;
            }
 
 
 // --------------------------- 15-10-2018  TAD-Interphase boundaries ~ 267kb -------------------------
 
            if (STATIC_BINDER_BLOCK_TYPE == "TAD-Interphase")
            {
                cout << ">> initializing **" << STATIC_BINDER_BLOCK_TYPE << "** static binder block sites ..." << endl;
                int cnt = 0;
                
				for (int i = 0; i < vecNode.size(); i ++)
				{
					if (i != 0 && i % INTERPHASE_TAD_SIZE == 0)
					{
						vecBinderBlockSite.push_back(i);
						vecDomains_DC.push_back(i); 
						cout << "... beadId = " << i << " is a condensin BLOCK based on Interphase TAD-size ..." << endl;
					}
				}
				cout << "    done!" << endl;
            }
 
 // --------------------------- 15-10-2018  TAD-Mitosis boundaries ~ 480kb -------------------------

            if (STATIC_BINDER_BLOCK_TYPE == "TAD-Mitosis")
            {
                cout << ">> initializing **" << STATIC_BINDER_BLOCK_TYPE << "** static binder block sites ..." << endl;
                int cnt = 0;
                
				for (int i = 0; i < vecNode.size(); i ++)
				{
					if (i != 0 && i % MITOSIS_TAD_SIZE == 0)
					{
						vecBinderBlockSite.push_back(i);
						vecDomains_DCMitosis.push_back(i);
						cout << "... beadId = " << i << " is a condensin BLOCK based on Mitosis TAD-size ..." << endl;
					}
				}
				cout << "    done!" << endl;
            }
            
// -------------------------------------19-10-2018 Interphase-Mitosis transition for TAD boundaries -----------------------------------------------------------------------------------------------

            if (STATIC_BINDER_BLOCK_TYPE == "Interphase-TAD-Mitosis")
            {
                cout << ">> initializing **" << STATIC_BINDER_BLOCK_TYPE << "** static binder block sites in interphase..." << endl;
                int cnt = 0;
                
				for (int i = 0; i < vecNode.size(); i ++)
				{
					if (i != 0 && i % INTERPHASE_TAD_SIZE == 0)
					{
						vecBinderBlockSite.push_back(i);
						vecDomains_DC.push_back(i); 
						cout << "... beadId = " << i << " is a condensin BLOCK based on Interphase TAD-size ..." << endl;
					}
				}
				cout << "    done!" << endl;
                cout << ">> initializing **" << STATIC_BINDER_BLOCK_TYPE << "** static binder block sites in mitosis ..." << endl;
                cnt = 0;
                
				for (int i = 0; i < vecNode.size(); i ++)
				{
					if (i != 0 && i % MITOSIS_TAD_SIZE == 0)
					{
						vecBinderBlockSiteMitosis.push_back(i);
						vecDomains_DCMitosis.push_back(i);
						cout << "... beadId = " << i << " is a condensin BLOCK based on Mitosis TAD-size ..." << endl;
					}
				}
				cout << "    done!" << endl;
            }


//------------------------------------ 17-10-2018 Domain boundaries by Yasu -------------------------------------------------------------

            if (STATIC_BINDER_BLOCK_TYPE == "domains")
            {
                cout << ">> initializing **" << STATIC_BINDER_BLOCK_TYPE << "** static binder block sites ..." << endl;
                string iFilename_domains;
                if (argc == 1)
                {
                    cout << "\nPlease input condesin chipseq filename (without suffix) as argument! \n";
                    exit(0);
                }
                if (argc == 2)
                {
                    iFilename_domains = argv[1]+SUFFIX_domains;
                    //cout << ">> reading file " << iFilename_domains << " ... "<< endl;
                    ifstream file_domains (iFilename_domains);
                    string line_domains;
                    if (file_domains.is_open())
                    {
                        int cnt = 0;
                        //getline (file_chipseq_cohesin, line_chipseq_cohesin); // skip the first line
                        while (getline (file_domains, line_domains))
                        {
                            stringstream ss (line_domains);
                            vector<string> items;
                            string buf;
                            while (ss >> buf)
                                items.push_back(buf);
                            //int chromoId, loc_bp, loc_bead;
                            
							int loc_bp_start, loc_bp_end, loc_bp, loc_bead, loc_diff;
							loc_bp_start = stoi(items[0]);
							loc_bp_end = stoi(items[1]);
							loc_diff   = (loc_bp_end - loc_bp_start)/2;
							loc_bp = loc_bp_start + loc_diff;  // unit is bp
							loc_bead = loc_bp/(BP_SEP_NUCLEOSOME)*0.1;
							cout << " start " << loc_bp_start << endl;
							cout << " loc_bp = " << loc_bp << endl;

                        
							if (loc_bead < NUM_BEAD && find(vecBinderBlockSite.begin(), vecBinderBlockSite.end(), loc_bead) == vecBinderBlockSite.end())
							{
								vecBinderBlockSite.push_back(loc_bead);
								vecDomains_DC.push_back(loc_bead);                         
                                cnt ++;
                                cout << "... beadId = " << loc_bead << " is a condensin BLOCKing site ... " << endl;
                            }
                        }
                        cout << "total number of BLOCKs : " << cnt << endl;
                        //cout << " total number of BLOCKs2: " << vecBinderBlocksSite.size() << endl;

                    }
                    //cout << " total number of BLOCKs2: " << vecBinderBlocksSite.size() << endl;
                    
                }
                cout << "    done!" << endl;
            }


//------------------------------------ 18-10-2018 Domain boundaries by Yasu -------------------------------------------------------------

            if (STATIC_BINDER_BLOCK_TYPE == "Interphase-domains-Mitosis")
            {
                cout << ">> initializing **" << STATIC_BINDER_BLOCK_TYPE << "** static binder block sites ..." << endl;
                string iFilename_domains;
                if (argc == 1)
                {
                    cout << "\nPlease input condesin chipseq filename (without suffix) as argument! \n";
                    exit(0);
                }
                if (argc == 3)
                {
                    iFilename_domains = argv[1]+SUFFIX_domains1;
                    //cout << ">> reading file " << iFilename_domains << " ... "<< endl;
                    ifstream file_domains (iFilename_domains);
                    string line_domains;
                    if (file_domains.is_open())
                    {
                        int cnt = 0;
                        //getline (file_chipseq_cohesin, line_chipseq_cohesin); // skip the first line
                        while (getline (file_domains, line_domains))
                        {
                            stringstream ss (line_domains);
                            vector<string> items;
                            string buf;
                            while (ss >> buf)
                                items.push_back(buf);
                            //int chromoId, loc_bp, loc_bead;
                            
							int loc_bp_start, loc_bp_end, loc_bp, loc_bead, loc_diff;
							loc_bp_start = stoi(items[0]);
							loc_bp_end = stoi(items[1]);
							loc_diff   = (loc_bp_end - loc_bp_start)/2;
							loc_bp = loc_bp_start + loc_diff;  // unit is bp
							loc_bead = loc_bp/(BP_SEP_NUCLEOSOME)*0.1;
							cout << " start " << loc_bp_start << endl;
							cout << " loc_bp = " << loc_bp << endl;

                        
							if (loc_bead < NUM_BEAD && find(vecBinderBlockSite.begin(), vecBinderBlockSite.end(), loc_bead) == vecBinderBlockSite.end())
							{
								vecBinderBlockSite.push_back(loc_bead);
								vecDomains_DC.push_back(loc_bead);                         
                                cnt ++;
                                cout << "... beadId = " << loc_bead << " is a condensin BLOCKing site ... " << endl;
                            }
                        }
                        cout << "total number of BLOCKs : " << cnt << endl;
                        //cout << " total number of BLOCKs2: " << vecBinderBlocksSite.size() << endl;

                    }
                    //cout << " total number of BLOCKs2: " << vecBinderBlocksSite.size() << endl;
                    
                }
                cout << "    done!" << endl;
                
                string iFilename_domains2;

                if (argc == 3)
                {
                    iFilename_domains2 = argv[2]+SUFFIX_domains2;
                    //cout << ">> reading file " << iFilename_domains << " ... "<< endl;
                    ifstream file_domains2 (iFilename_domains2);
                    string line_domains2;
                    if (file_domains2.is_open())
                    {
                        int cnt = 0;
                        //getline (file_chipseq_cohesin, line_chipseq_cohesin); // skip the first line
                        while (getline (file_domains2, line_domains2))
                        {
                            stringstream ss (line_domains2);
                            vector<string> items;
                            string buf;
                            while (ss >> buf)
                                items.push_back(buf);
                            //int chromoId, loc_bp, loc_bead;
                            
							int loc_bp_start, loc_bp_end, loc_bp, loc_bead, loc_diff;
							loc_bp_start = stoi(items[0]);
							loc_bp_end = stoi(items[1]);
							loc_diff   = (loc_bp_end - loc_bp_start)/2;
							loc_bp = loc_bp_start + loc_diff;  // unit is bp
							loc_bead = loc_bp/(BP_SEP_NUCLEOSOME)*0.1;
							cout << " start " << loc_bp_start << endl;
							cout << " loc_bp = " << loc_bp << endl;

                        
							if (loc_bead < NUM_BEAD && find(vecBinderBlockSiteMitosis.begin(), vecBinderBlockSiteMitosis.end(), loc_bead) == vecBinderBlockSiteMitosis.end())
							{
								vecBinderBlockSiteMitosis.push_back(loc_bead);
								vecDomains_DCMitosis.push_back(loc_bead);                         
                                cnt ++;
                                cout << "... beadId = " << loc_bead << " is a condensin BLOCKing site ... " << endl;
                            }
                        }
                        cout << "total number of BLOCKs : " << cnt << endl;
                        //cout << " total number of BLOCKs2: " << vecBinderBlocksSite.size() << endl;

                    }
                    //cout << " total number of BLOCKs2: " << vecBinderBlocksSite.size() << endl;
                    
                }
                cout << "    done!" << endl;
            }


//---------------------------------------------------------------------------------------------------------------------------------------
            
            if (STATIC_BINDER_BLOCK_TYPE == "chipseq")
            {
                cout << ">> initializing **" << STATIC_BINDER_BLOCK_TYPE << "** static binder block sites ..." << endl;
                string iFilename_chipseq_cohesin;
                if (argc == 1)
                {
                    cout << "\nPlease input condesin chipseq filename (without suffix) as argument! \n";
                    exit(0);
                }
                if (argc == 2)
                {
                    iFilename_chipseq_cohesin = argv[1]+SUFFIX_chipseq;
                    //cout << ">> reading file " << iFilename_chipseq_cohesin << " ... "<< endl;
                    ifstream file_chipseq_cohesin (iFilename_chipseq_cohesin);
                    string line_chipseq_cohesin;
                    if (file_chipseq_cohesin.is_open())
                    {
                        int cnt = 0;
                        //getline (file_chipseq_cohesin, line_chipseq_cohesin); // skip the first line
                        while (getline (file_chipseq_cohesin, line_chipseq_cohesin))
                        {
                            stringstream ss (line_chipseq_cohesin);
                            vector<string> items;
                            string buf;
                            while (ss >> buf)
                                items.push_back(buf);
                            //int chromoId, loc_bp, loc_bead;
                            
							int loc_bp_start, loc_bp_end, loc_bp, loc_bead, loc_diff;
							loc_bp_start = stoi(items[0]);
							loc_bp_end = stoi(items[1]);
							loc_diff   = (loc_bp_end - loc_bp_start)/2;
							loc_bp = loc_bp_start + loc_diff;  // unit is bp
							loc_bead = loc_bp/(BP_SEP_NUCLEOSOME)*0.1;
							cout << " start " << loc_bp_start << endl;
							cout << " loc_bp = " << loc_bp << endl;

                        
							if (loc_bead < NUM_BEAD && find(vecBinderBlockSite.begin(), vecBinderBlockSite.end(), loc_bead) == vecBinderBlockSite.end())
							{
								vecBinderBlockSite.push_back(loc_bead);                         
                                                                cnt ++;
                                                                cout << "... beadId = " << loc_bead << " is a condensin BLOCKing site ... " << endl;
                                                        }
                        }
                        cout << "total number of BLOCKs : " << cnt << endl;
                       // cout << " total number of BLOCKs2: " << vecBinderBlocksSite.size() << endl;

                    }
                }
                cout << "    done!" << endl;
            }
        }
        
        else
            cout << ">> NO static binder block sites ..." << endl;
    }

    // --------------------------- 18-10-2018 Initialize DOMAINS predicted by Yasu log2 direcitonality ----------------------------
    
    //unordered_map<int, vector<int>> Domain_CSB;
    
    //if( vecDomains_DC.size() > 0)
    //{
		map<int, vector<int>> Domain_CSB;
		vector<int>vecCSB_domain;
		cout << vecDomains_DC.size() << " number of domains in interphase " << endl;
		
		for (int j = 1; j <= vecDomains_DC.size(); j++)
		{
			for (int i = 0; i < vecCondensinSite.size(); i++)
			{
				int bIdCSB = vecCondensinSite[i];
				if (bIdCSB >= vecDomains_DC[j-1] && bIdCSB < vecDomains_DC[j])
				{
					vecCSB_domain.push_back(vecCondensinSite[i]);
					Domain_CSB[j].push_back(vecCondensinSite[i]);
				}
			}
			//Domain_CSB[j].push_back(vecCSB_domain);
		}
		for ( auto map_iter = Domain_CSB.begin(); map_iter != Domain_CSB.end(); ++map_iter)
		{
			cout << "key: " << map_iter->first << "    value: [";
			for (auto vec_iter = map_iter->second.begin() ; vec_iter != map_iter->second.end(); ++vec_iter)
			{
				cout << *vec_iter << " ";
			}
			cout << "]\n" ;
		}
	//}
	
	//if( vecDomains_DCMitosis.size() > 0)
	//{
			
		map<int, vector<int>> Domain_CSBMitosis;
		vector<int>vecCSB_domainMitosis;
		cout << vecDomains_DCMitosis.size() << " number of domains in mitosis " << endl;
		
		for (int j = 1; j <= vecDomains_DCMitosis.size(); j++)
		{
			for (int i = 0; i < vecCondensinSite.size(); i++)
			{
				int bIdCSB = vecCondensinSite[i];
				if (bIdCSB >= vecDomains_DC[j-1] && bIdCSB < vecDomains_DCMitosis[j])
				{
					vecCSB_domainMitosis.push_back(vecCondensinSite[i]);
					Domain_CSBMitosis[j].push_back(vecCondensinSite[i]);
				}
			}
			//Domain_CSB[j].push_back(vecCSB_domain);
		}
	
		for ( auto map_iter = Domain_CSBMitosis.begin(); map_iter != Domain_CSBMitosis.end(); ++map_iter)
		{
			cout << "key: " << map_iter->first << "    value: [";
			for (auto vec_iter = map_iter->second.begin() ; vec_iter != map_iter->second.end(); ++vec_iter)
			{
				cout << *vec_iter << " ";
			}
			cout << "]\n" ;
		}
	//}
	
	
    // initialize genes

    // initialize bead coordinate (only Chromosome beads)
    vector<double> randNum01;
    for (int i = 0; i < vecNode.size()*9; i ++)
        randNum01.push_back(uniform01(generator));  // used for coord and veloc initiation
    if (SIM_TYPE == "RECONSTRUCT")
    {
        string initType = "random";
        initDynamics(vecChromosome, vecNode, argc, argv, initType, randNum01);
    }
    else if (SIM_TYPE == "SIMULATION")
    {
        //string initType = "randomChain";        // this is random walk chain
        //string initType = "randomChainBinned";    // this is cylindrical extending from centromere to center of mass
        //[NOT USED] string initType = "randomChainUniformTether";
        //string initType = "randomChainBinaryInsertion";  // this is a chain formed by recursively inserting a bead in between
        string initType = "randomChainBinnedCylinder";    // this is cylindrical extending from centromere to center of mass
        //string initType = "loadPDB";    // this is for checkpoint or decondensation - starting from a PDB-defined configuration

        initDynamicsBrownianSimulation(vecChromosome, vecNode, vecFluoroSite, initType, randNum01, argc, argv);
    }

    // initialize loop-extrusion binders (both vector of binders and coordinates)
    if (SIM_TYPE == "SIMULATION")
    {
        initBindersBrownianSimulation(vecBinder, vecNode, vecCondensinSite, vecChromosome);
        //cout << vecBinder.size() << endl;
    }


    /* ---------- simulation starts here ----------- */

    double t_now = -T_PREP-T_INIT;
    bool flagIsPrep = true;
    bool flagIsInitConfig = true;
    long nIter = 0;

    int cntWriteMSD = 0;
    bool flagWriteMSD = false;
    double tWriteMSD = 0;
    
    // --------------------------
    
    bool flagBindersLoading = false;
    bool flagBindersUnloading = false;
    
    
    //------------------------------------------------
    bool flagDC_Interphase = true;

    //-------------------- 19-10-2018 -------------------
    bool flagLE_ON = false;
    bool flagIsInterphase = true;
    
    // -------------------------- 10-08-2018 ---------------
    
    bool flagBindersDynamicsCoupled = false;
    /*
    int cntWriteDynamicBinderCoupled = 0;
    bool flagWriteDynamicBinderCoupled = false;
    double tWriteDynamicBinderCoupled = 0
    */
    
/////////////////////////////////////////////////////////////////////////
    int cntWriteDC = 0;
    bool flagWriteDC = true;
    double tWriteDC = 0;
/////////////////////////////////////////////////////////////////////////

// ------------------------------------01-08-2018 ---------------------------

    int cntWriteDynamicBinder = 0;
    bool flagWriteDynamicBinder = false;
    double tWriteDynamicBinder = 0;
    
// ---------------------------------------------------------------------------

    while (t_now <= T)
    {
        // single iteration of simulation
        energy["total"] = 0; energy["kinetic"] = 0; energy["potential"] = 0;
        bool FLAG_INTERACT_ON = false;
        if (nIter >= CNT_START_INTERACT)
            FLAG_INTERACT_ON = true;

        if (SIM_TYPE == "RECONSTRUCT")
        {
			vector<double> randNormNum0;
            for (int i = 0; i < vecNode.size()*6; i ++)
            {
                randNormNum0.push_back(normal0(generator));
                //cout << randNormNum0[randNormNum0.size()-1] << endl;
            }
            oneIter(vecChromosome, vecNode, vecBeadPairInteract, vecBeadPairInteractFreq, vecFluoroSite, energy, voxMap, FLAG_INTERACT_ON, randNormNum0, flagIsPrep); // if using voxMap
            //oneIter(vecChromosome, vecNode, vecBeadPairInteract, energy, octree, FLAG_INTERACT_ON); // if using octree
		}
		
//----------------------------------


		else if (SIM_TYPE == "SIMULATION")
				{
					vector<double> randNormNum0;
					for (int i = 0; i < vecNode.size()*6; i ++)
					{
						randNormNum0.push_back(normal0(generator));
						//cout << randNormNum0[randNormNum0.size()-1] << endl;
					}
					//oneIterBrownianSimulation(vecChromosome, vecNode, vecCondensinSite, energy, voxMap, randNormNum0);
					
				   // ---------- 21.6.2018 --------------------------- 
					if (nIter != 1 && flagBindersLoading == true && flagIsPrep == false && nIter % T_RESIDENCE == 0  && (nIter*DT) == (SWITCH_INTERPHASE_MITOSIS + T_PREP))
					{
						initBindersLoadingBrownianSimulation(vecBinder, vecNode, vecCondensinSite, vecChromosome, vecLoadedBinders, vecCondensinSiteEmpty,vecBinderIds2);
						//cout << "Dynamic loading of Binders....." << endl;
						////cout << "Empty sites: " << vecCondensinSiteEmpty.size() << endl;
						//cout <<  vecLoadedBinders.size()  << " Binders loaded. " << endl;
						////cout << "Empty sites: " << vecCondensinSiteEmpty.size() << endl;
						////cout << "vecCondensinSite: " << endl;
						/*for (int i = 0; i < vecCondensinSite.size(); i ++)
						{
							cout << vecCondensinSite[i] << endl;
						}
						*/
						//
						/*for (int j = 0; j < vecCondensinSite.size(); j ++)
						{
							cout << vecBinderIds2[j] << endl;
						}
						*/
					
					}
					if (nIter != 1 && flagBindersUnloading == true && flagIsPrep == false && nIter % T_RESIDENCE == 0 && (nIter*DT) == (SWITCH_INTERPHASE_MITOSIS + T_PREP))
					{
						initBindersUnloadingBrownianSimulation(vecBinder, vecNode, vecCondensinSite, vecChromosome, vecLoadedBinders);
						//cout << "Unloading of Binders....." << endl;
						//cout << "Binders to be deleted: " << vecUnloadedBinders.size() << endl;						
						//cout << "Empty sites after dynamic loading and unloading: " << vecCondensinSiteEmpty2.size() << endl;
						////cout << "Number of Binders before unloading: " << vecBinderIds2.size() / 2 << endl;
						//cout << "Number of Binders before after dynamic loading and unloading: " << vecBinder.size() << endl;
									

					}
					
					if (flagIsPrep == false && (nIter*DT) == (SWITCH_INTERPHASE_MITOSIS + T_PREP) && SWITCH_INTERPHASE_MITOSIS > 0)
					{
						flagDC_Interphase = false;
						cout << "BLAAAAA ******************************************************** SWITCH*********" << endl;
						vecBinderBlockSite.clear();
						vecBinderBlockSite = vecBinderBlockSiteMitosis;
						
											
						//map<int, vector<int>> Domain_CSB;
						//vector<int>vecCSB_domain;
						Domain_CSB.clear();
						
						Domain_CSB = Domain_CSBMitosis;
						cout << vecDomains_DCMitosis.size() << " number of domains in mitosis " << endl;
						
						

						for ( auto map_iter = Domain_CSB.begin(); map_iter != Domain_CSB.end(); ++map_iter)
						{
							cout << "key: " << map_iter->first << "    value: [";
							for (auto vec_iter = map_iter->second.begin() ; vec_iter != map_iter->second.end(); ++vec_iter)
							{
								cout << *vec_iter << " ";
							}
							cout << "]\n" ;
						}
						flagLE_ON = true;
						flagIsInterphase = false;
						
						/*for (int i = 0; i <= vecBinderBlockSite.size()+1; i=i+2)
						{
								vecBinderBlockSite.erase(vecBinderBlockSite.begin() + (i-1));
		
						} */
					}

					oneIterBrownianSimulationOverdamp(vecChromosome, vecNode, vecCondensinSite, vecBinderBlockSite, vecFluoroSite, vecBinder, energy, voxMap, randNormNum0, flagIsPrep, flagIsInitConfig, capturer, hit_valence, Domain_CSB, flagLE_ON, flagIsInterphase);
				}
				// print information


		
		
//------------------------------------
		//---- static implementation ------------
      /*  else if (SIM_TYPE == "SIMULATION")
        {
            vector<double> randNormNum0;
            for (int i = 0; i < vecNode.size()*6; i ++)
            {
                randNormNum0.push_back(normal0(generator));
                //cout << randNormNum0[randNormNum0.size()-1] << endl;
            }
            //oneIterBrownianSimulation(vecChromosome, vecNode, vecCondensinSite, energy, voxMap, randNormNum0);
            
            oneIterBrownianSimulationOverdamp(vecChromosome, vecNode, vecCondensinSite, vecBinderBlockSite, vecFluoroSite, vecBinder, energy, voxMap, randNormNum0, flagIsPrep, flagIsInitConfig, capturer);
        }
        */
        
        
        // print information
        if (nIter % FREQ_PRINT == 0)
        {
            cout << "\n --------- t_now = " << t_now << " ---------" << endl;
            /*
            cout << " --------- Energy (x1e-18 J/kg ) = " << energy["total"] << endl;
            cout << " ---------        (potential   ) = " << energy["potential"] << endl;
            cout << " ---------        (kinetic     ) = " << energy["kinetic"] << endl;
            */
            printChromoInfo(vecChromosome, vecNode, vecBinder);
        }

        // write coordinate to file
        if (nIter % FREQ_WRITE == 0)
        {
            writeChromoInfo(vecChromosome, vecNode, vecFluoroSite, vecCondensinSite, vecBinderBlockSite, vecBinder, voxMap, t_now);
        }
        
        if (nIter % MSD_FREQ_WRITE == 0)
        {
            flagWriteMSD = true;
            tWriteMSD = t_now;
        }
        if (flagWriteMSD)
        {
            if ( cntWriteMSD % MSD_WRITE_INTERVAL == 0 )
                writeMSD(vecNode, t_now);
            cntWriteMSD ++;

            if ( cntWriteMSD > MSD_WRITE_DURATION )
            {
                flagWriteMSD = false;
                cntWriteMSD = 0;
            }
        }
 ////////////////////////////////////////////////////////////////////////  
       /*
        if (nIter % DC_FREQ_WRITE == 0)
        {
			flagWriteDC = true;
			tWriteDC = t_now;
		}
		
		if (flagWriteDC)
		{
			if ( cntWriteDC % DC_WRITE_INTERVAL == 0 )
				//writeDC2(capturer, t_now);
				writeDC2_stats(capturer, hit_valence, t_now);
			cntWriteDC ++;
			
			
			if ( cntWriteDC > DC_WRITE_DURATION )
			{
				flagWriteDC = false;
				cntWriteDC = 0;
			}
			
			
		}
		*/
		
// ----------------------------------- 03-08-2018 -----------------------
		if (nIter % DC_FREQ_WRITE == 0)
        {
			flagWriteDC = true;
			tWriteDC = t_now;
		}
		
		if (flagWriteDC == true)
		{
			//writeDC2(capturer, t_now);
			writeDC2_stats(capturer, hit_valence, t_now);
		}

// ----------------------------------- 01-08-2018 -----------------------				
		//if (flagWriteDynamicBinder == true && nIter % DB_FREQ_WRITE == 0)
		if (flagWriteDynamicBinder == true && nIter % (FREQ_WRITE) == 0)  // FREQ_WRITE/10000

		{

			writeDynamicsBinders(vecBinder, vecLoadedBinders, \
				                 vecUnloadedBinders, vecCondensinSiteEmpty, t_now, vecUnLoadedBinders);
			/*cntWriteDynamicBinder ++;
			
			if ( cntWriteDynamicBinder > DC_WRITE_DURATION )
			{
				flagWriteDC = false;
				cntWriteDC = 0;
			}
			*/ 	
		}
/////////////////////////////////////////////////////////////////////////
		//writeDC2(capturer, t_now);
        /*if (nIter % FREQ_WRITE == 0)
        {
            //writeDC2(capturer, t_now);
        }
        */
/////////////////////////////////////////////////////////////////////////

        t_now += DT;
        nIter ++;

        if (t_now >= 0)
            flagIsPrep = false;
        if (t_now >= -T_PREP)
            flagIsInitConfig = false;
    }

    return 0;
}



/* ------------------ writing model parameters -------------------*/
void writeParameterInfo(void )  // all parameters should be defined as global parameters
{
    cout << "-------------------------------------------------------------------" << endl;
    cout << "Writing parameter information ... " << endl;
    cout << "-------------------------------------------------------------------" << endl;

    ofstream pFile;
    string pFileName;
    stringstream stream;

    stringstream PROC_ID_SS; PROC_ID_SS << PROC_ID;

    pFileName = "_PID_"+ PROC_ID_SS.str() +"_all_parameters_o.csv";

    pFile.open(pFileName, ios::out | ios::trunc); pFile.close();    // delete existing content

    pFile.open(pFileName, ios::app | ios::binary);

    if (pFile.is_open())
    {
        // header file
        pFile << "FLAG"         << "\t"
              << "Parameter"    << "\t"
              << "Value"        << "\t"
              << "Unit"         << "\n";
        // parameters
        pFile << "0" << "\t" << "SIM_TYPE"     << "\t" << SIM_TYPE     << "\t" << "N/A"        << "\n";
        pFile << "0" << "\t" << "DIFFUSION_CAPTURE_ON"                 << "\t" << DIFFUSION_CAPTURE_ON     << "\t" << "N/A"        << "\n";
        pFile << "0" << "\t" << "LOOP_EXTRUSION_ON"                    << "\t" << LOOP_EXTRUSION_ON        << "\t" << "N/A"        << "\n";

        pFile << "0" << "\t" << "T"            << "\t" << T            << "\t" << "seconds"    << "\n";
        pFile << "0" << "\t" << "DT"           << "\t" << DT           << "\t" << "seconds"    << "\n";
        pFile << "0" << "\t" << "T_PREP"       << "\t" << T_PREP       << "\t" << "seconds"    << "\n";

        pFile << "0" << "\t" << "DIM_VOXMAP"   << "\t" << DIM_VOXMAP   << "\t" << "N/A"        << "\n";
        pFile << "0" << "\t" << "OMP_ON"       << "\t" << OMP_ON       << "\t" << "N/A"        << "\n";
        pFile << "0" << "\t" << "N_THREAD"     << "\t" << NUM_THREADS  << "\t" << "1"          << "\n";

        // constants
        pFile << "0" << "\t" << "NUM_BEAD"     << "\t" << NUM_BEAD     << "\t" << "1"          << "\n";
        pFile << "0" << "\t" << "N_BP_PER_NUCLEO"      << "\t" << N_BP_PER_NUCLEO              << "\t" << "bp"          << "\n";
        pFile << "0" << "\t" << "BP_SEP_NUCLEOSOME"    << "\t" << BP_SEP_NUCLEOSOME            << "\t" << "bp"          << "\n";
        pFile << "0" << "\t" << "BP_SEP_CONDENSIN_BINDING"     << "\t" << BP_SEP_CONDENSIN_BINDING     << "\t" << "bp"  << "\n";
        pFile << "0" << "\t" << "COARSE_FACTOR"         << "\t" << COARSE_FACTOR         << "\t" << "1"          << "\n";


        // entropic
        pFile << "1" << "\t" << "GAMMA_VISCOUS"<< "\t" << GAMMA_VISCOUS<< "\t" << "kg/s"       << "\n";
        pFile << "1" << "\t" << "F_ENTROPIC_VARIANCE"  << "\t" << F_ENTROPIC_VARIANCE          << "\t" << "pN*pN"       << "\n";

        // tension
        pFile << "1" << "\t" << "K_TENSION"    << "\t" << K_TENSION    << "\t" << "pN/nm"      << "\n";
        pFile << "1" << "\t" << "L0_TENSION"   << "\t" << L0_TENSION   << "\t" << "nm"         << "\n";

        // repulsion
        pFile << "1" << "\t" << "F_REPULSION"  << "\t" << F_REPULSION  << "\t" << "pN"         << "\n";
        pFile << "1" << "\t" << "LTHR1_REPULSION"      << "\t" << LTHR1_REPULSION              << "\t" << "nm"         << "\n";

        // attraction
        pFile << "1" << "\t" << "K_ATTRACTION" << "\t" << K_ATTRACTION << "\t" << "pN/nm"      << "\n";

        // diffusion capture dissociation
        pFile << "1" << "\t" << "P_DISSOCIATE_CONDENSIN" << "\t" << P_DISSOCIATE_CONDENSIN << "\t" << "1"      << "\n";
        
        // loop extrusion rate ~ walking probability
        pFile << "1" << "\t" << "P_BINDER_SLIDING"         << "\t" << P_BINDER_SLIDING         << "\t" << "1"          << "\n";


        // binder tension
        pFile << "1" << "\t" << "K_TENSION_BINDER"    << "\t" << K_TENSION_BINDER    << "\t" << "pN/nm"      << "\n";
        pFile << "1" << "\t" << "L0_TENSION_BINDER"   << "\t" << L0_TENSION_BINDER   << "\t" << "nm"         << "\n";
        pFile << "1" << "\t" << "K_TENSION_BINDER_BEAD"    << "\t" << K_TENSION_BINDER_BEAD    << "\t" << "pN/nm"      << "\n";
        pFile << "1" << "\t" << "L0_TENSION_BINDER_BEAD"   << "\t" << L0_TENSION_BINDER_BEAD   << "\t" << "nm"         << "\n";
        
        pFile << "0" << "\t" << "STATIC_BINDER_BLOCK_ON"   << "\t" << STATIC_BINDER_BLOCK_ON   << "\t" << "N/A"        << "\n";
        pFile << "0" << "\t" << "STATIC_BINDER_BLOCK_TYPE" << "\t" << STATIC_BINDER_BLOCK_TYPE << "\t" << "N/A"        << "\n";
        
        // mechanism 
        pFile << "0" << "\t" << "SEPARATE_DYNAMICS_BINDER" << "\t" << SEPARATE_DYNAMICS_BINDER << "\t" << "N/A"        << "\n";
        pFile << "0" << "\t" << "COUPLED_DYNAMICS_BINDER" << "\t" << COUPLED_DYNAMICS_BINDER << "\t" << "N/A"        << "\n";
        pFile << "0" << "\t" << "P_ACTIVATION_BINDER" << "\t" << P_ACTIVATION_BINDER << "\t" << "N/A"        << "\n";
        pFile << "0" << "\t" << "P_DISSOCIATE_BINDER" << "\t" << P_DISSOCIATE_BINDER << "\t" << "N/A"        << "\n";
        
        pFile << "0" << "\t" << "DELETE_BINDER" << "\t" << DELETE_BINDER << "\t" << "*100 = %"        << "\n";
        pFile << "0" << "\t" << "BINDER_POSITION" << "\t" << BINDER_POSITION << "\t" << "N/A"        << "\n";
        pFile << "0" << "\t" << "NUM_CONDENSIN_PEAKS" << "\t" << NUM_CONDENSIN_PEAKS << "\t" << "N/A"        << "\n";
        
        pFile << "0" << "\t" << "DIFFUSION_CAPTURE" << "\t" << DIFFUSION_CAPTURE << "\t" << "N/A"        << "\n";
        pFile << "0" << "\t" << "P_DC_PROBAB" << "\t" << P_DC_PROBAB << "\t" << "N/A"        << "\n";

        pFile << "0" << "\t" << "LOOP_EXTRUSION" << "\t" << LOOP_EXTRUSION << "\t" << "N/A"        << "\n";
        pFile << "0" << "\t" << "P_LE_PROBAB" << "\t" << P_DC_PROBAB << "\t" << "N/A"        << "\n";

    }

    pFile.close();

}
/* ---------------------------------------------------------------*/


/* -------------- Notes ---------------
    Week of 2017.12.18
    * Achievements
      (1) first stage modelling dynamics according to Hi-C matrix
    * Expectations
      (1) output simulated connectivity/contact matrix
      (2) output coordinates of chromosomes in separate files (e.g., PDB format, for rod/stick representation in visualisation)

   ------------------------------------ */
