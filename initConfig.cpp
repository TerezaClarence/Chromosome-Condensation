/*
    File: initConfig.cpp
    Function: Initialize chromosomes and nuclear envelope
    Model: chromoCell
    Created: 7 December, 2017 (XF)
*/

#include "initConfig.hpp"

// ===================== Classes =========================

/* -------------- Node ---------------- */
Node::Node() {}
Node::Node(int id, int hostChromoId)
{
    this->id = id;
    this->hostChromoId = hostChromoId;
}
Node::Node(int id, int hostChromoId, dVec coord, dVec veloc, dVec accel, dVec accelPrev)
{
    this->id = id;
    this->hostChromoId = hostChromoId;
    this->coord = coord;
    this->veloc = veloc;
    this->accel = accel;
    this->accelPrev = accelPrev;
}
Node::~Node() {}

/* -------------- Chromosome ---------------- */
Chromosome::Chromosome() {}
Chromosome::Chromosome(int id, vector<int> vecBeadIds)
{
    this->id = id;
    this->vecBeadIds = vecBeadIds;
}
Chromosome::Chromosome(int id, vector<int> vecBeadIds, vector<int> vecGeneIds, vector<int> vecBinderIds)
{
    this->id = id;
    this->vecBeadIds = vecBeadIds;
    this->vecGeneIds = vecGeneIds;
    this->vecBinderIds = vecBinderIds;
}
Chromosome::~Chromosome() {}

/* -------------- Binder ---------------- */
Binder::Binder() {}
Binder::Binder(int id, int attachChromoId, vector<bool> canSlide, vector<int> vecBeadIds, vector<int> vecAttachBeadIds)
{
    this->id = id;
    this->attachChromoId = attachChromoId;
    this->canSlide = canSlide;
    this->vecBeadIds = vecBeadIds;
    this->vecAttachBeadIds = vecAttachBeadIds;
}
Binder::~Binder() {}

/* -------------- Domain ---------------- */
Domain::Domain() {}
Domain::Domain(string name, int hostChromoId, iPair location)
{
    this->name = name;
    this->hostChromoId = hostChromoId;
    this->location = location;
}
Domain::~Domain() {}

/* -------------- VoxMap ---------------- */
VoxMap::VoxMap() {}
VoxMap::VoxMap(iVec dim, unordered_map<int, vector<int>> mapBeadIds, unordered_map<int, vector<int>> mapBeadIdsHalo)
{
    this->dim = dim;
    this->mapBeadIds = mapBeadIds;
    this->mapBeadIdsHalo = mapBeadIdsHalo;
}
VoxMap::~VoxMap() {}

/* -------------- Octree ---------------- */
Leaf::Leaf() {}
Leaf::Leaf(vector<bVec> position, vector<int> vecBeadIds, vector<Leaf> vecLeaf, bool flagBranched)
{
    this->position = position;
    this->vecBeadIds = vecBeadIds;
    this->vecLeaf = vecLeaf;
    this->flagBranched = flagBranched;
}
Leaf::~Leaf() {}
Octree::Octree() {}
Octree::Octree(int depth, vector<Leaf> vecLeaf)
{
    this->depth = depth;
    this->vecLeaf = vecLeaf;
}
Octree::~Octree() {}

// ===================== Functions =========================
// --- following functions are for RECONSTRUCTION objective ---
void initChromosome(vector<Chromosome> & vecChromosome, vector<Node> & vecNode, \
                    vector<iPair> & vecBeadPairInteract, vector<double> & vecBeadPairInteractFreq, int argc, char ** argv)
{
    double **matChromoContFreq = NULL;
    string iFilename_normMatrix, iFilename_numberBead;
    if (argc == 1)
    {
        cout << "\nPlease input filename (without suffix) as argument! \n";
        exit(0);
    }
    if (argc == 2)
    {
        iFilename_normMatrix = argv[1]+SUFFIX_normMatrix;
        iFilename_numberBead = argv[1]+SUFFIX_numberBead;

        cout << ">> reading file " << iFilename_numberBead << " ... "<< endl;
        ifstream file_numberBead (iFilename_numberBead);
        string line_numberBead;
        if (file_numberBead.is_open())
        {
            getline (file_numberBead, line_numberBead); // skip the first line
            while (getline (file_numberBead, line_numberBead))
            {
                stringstream ss (line_numberBead);
                vector<string> items;
                string buf;
                while (ss >> buf)
                    items.push_back(buf);
                //cout << line_numberBead << endl;
                //cout << items[0] << ", " << items[1] << endl;   // items[1] is chromosome id, items[2] is number of beads

                /* ---------- create Node and Chromosome ----------- */
                int chromoId, nBead;
                chromoId = stoi(items[0]);
                nBead    = stoi(items[1]);
                //cout << chromoId << ", " << nBead << endl;

                int currBeadId = vecNode.size();
                vector<int> vecBeadIds;
                for (int i = 0; i < nBead; i ++)
                {
                    Node bead0 = Node(currBeadId, chromoId);
                    vecNode.push_back(bead0);
                    vecBeadIds.push_back(currBeadId);
                    currBeadId ++;
                }
                Chromosome chromo0 = Chromosome(chromoId, vecBeadIds);
                vecChromosome.push_back(chromo0);
            }
            file_numberBead.close();
            cout << "    done!" << endl;
        }
        else
        {
            cout << "Unable to open file " << iFilename_numberBead << endl;
            exit(1);
        }

        // allocate memory to contact frequency map
        const int nBeadTotal = vecNode.size();
        cout << ">> allocating memory to temporary matChromoContFreq with dimension " << nBeadTotal << " by " << nBeadTotal << endl;
        matChromoContFreq = new double* [nBeadTotal];
        for (int i = 0; i<nBeadTotal; i++)
            matChromoContFreq[i] = new double[nBeadTotal];
        for (int j = 0; j<nBeadTotal; j++)
        {
            for (int k = 0; k<nBeadTotal; k++)
                matChromoContFreq[j][k] = 0;
        }

        cout << "    done!" << endl;

        cout << ">> reading file " << iFilename_normMatrix << " ... "<< endl;
        ifstream file_normMatrix (iFilename_normMatrix);
        string line_normMatrix;
        if (file_normMatrix.is_open())
        {
            getline (file_normMatrix, line_normMatrix); // skip the first line
            int currRow = 0;
            double freqMax = 0;
            int cntOverCutoff = 0;
            while(getline (file_normMatrix, line_normMatrix))
            {
                stringstream ss (line_normMatrix);
                vector<string> items;
                string buf;
                while (ss >> buf)
                    items.push_back(buf);
                /* --------- initiate contact frequency matrix --------- */
                for (vector<string>::iterator it = items.begin(); it != items.end(); it++)
                {
                    double freq = stod (*it);
                    int currCol = it-items.begin();
                    matChromoContFreq[currRow][currCol] = freq;
                    if (freq > freqMax)
                        freqMax = freq;
                    if (freq != 0 && freq >= CUTOFF_NORMMAT && currCol > currRow + GENOMIC_LEAST_SEP_INTRA)
                    {
                        cntOverCutoff ++;
                        cout << "Bead id pairs with interaction: " << currRow << " " << currCol << endl;
                        vecBeadPairInteract.push_back({currRow, currCol});
                        vecBeadPairInteractFreq.push_back(freq);
                    }
                }
                currRow ++;
            }
            cout << "freqMax = " << freqMax << "; " << cntOverCutoff << " over cutoff (" << CUTOFF_NORMMAT << ")" << endl;
            //exit(88);
            file_normMatrix.close();

            cout << "    done!" << endl;
        }
        else
        {
            cout << "Unable to open file " << iFilename_normMatrix << endl;
            exit(1);
        }
    }
}
void initDynamics(vector<Chromosome> & vecChromosome, vector<Node> & vecNode, int argc, char ** argv,
                    string initType, vector<double> & randNum01)
{
    int k = 0;  // position in randNum01
    cout << ">> initializing bead coordinates with *" << initType << "* method ... " << endl;
    if (initType == "reconstruct")
    {

    }
    else if (initType == "random")
    {
        double rInit = 0.75*RAD_NUCLEUS;
        double r, theta, phi;
        for (vector<Node>::iterator it = vecNode.begin(); it != vecNode.end(); it ++)
        {
            // ref1: https://math.stackexchange.com/questions/87230/picking-random-points-in-the-volume-of-sphere-with-uniform-probability
            // ref2: http://corysimon.github.io/articles/uniformdistn-on-sphere/
            r     = rInit * cbrt(randNum01[k++]);
            theta = acos(2.* randNum01[k++] - 1);
            phi   = 2.* PI * randNum01[k++];
            dVec coord = {r*sin(theta)*cos(phi), r*sin(theta)*sin(phi), r*cos(theta)};
            it->set_coord(coord);
        }
    }
    cout << "    done!" << endl;

    cout << ">> initializing bead velocities & accelerations ... " << endl;
    for (vector<Node>::iterator it = vecNode.begin(); it != vecNode.end(); it ++)
    {
        //double s = 10.; // nm/s
        //dVec veloc = {s*randNum01[k++], s*randNum01[k++], s*randNum01[k++]};
        dVec veloc = {0., 0., 0.};
        dVec accel = {0., 0., 0.}, accelPrev = {0., 0., 0.};
        it->set_veloc(veloc);
        it->set_accel(accel);
        it->set_accelPrev(accelPrev);
    }
    cout << "    done!" << endl;
}
// ------------------------------------------------------------

// --- following functions are for Brownian SIMULATION objective ---
void initChromosomeBrownianSimulation(vector<Chromosome> & vecChromosome, vector<Node> & vecNode)
{
    /* ---------- create Node and Chromosome ----------- */
    int chromoId = 1, nBead = NUM_BEAD; // NUM_BEAD is derived from Hi-C data
    int currBeadId = vecNode.size();
    vector<int> vecBeadIds;
    for (int i = 0; i < nBead; i ++)
    {
        Node bead0 = Node(currBeadId, chromoId);
        vecNode.push_back(bead0);
        vecBeadIds.push_back(currBeadId);
        currBeadId ++;
    }
    Chromosome chromo0 = Chromosome(chromoId, vecBeadIds);
    vecChromosome.push_back(chromo0);
}
void initDynamicsBrownianSimulation(vector<Chromosome> & vecChromosome, vector<Node> & vecNode, vector<int> & vecFluoroSite, string initType, vector<double> & randNum01, int argc, char ** argv)
{
    int k = 0;
    cout << ">> initializing bead coordinates with *" << initType << "* method ... " << endl;
    if (initType == "hilbertChain")
    {

    }
    else if (initType == "randomChain") // need to rule out Binder nodes in the loop
    {
        double rInit = 0.95*RAD_NUCLEUS;
        double l = L0_TENSION, theta, phi;
        double x_prev = 0, y_prev = 0, z_prev = 0;
        double theta_prev = acos(2.* randNum01[k++] - 1);
        double phi_prev = 2.* PI * randNum01[k++];
        dVec coord;
        /* start from the center, then extend the chain by l with randomly chosen theta and phi */
        for (vector<Node>::iterator it = vecNode.begin(); it != vecNode.end(); it ++)
        {
            // random chain
            theta    = acos(2.* randNum01[k++] - 1);
            phi      = 2.* PI * randNum01[k++];

            // "persistent" random chain
            //double d_theta = PI*0.1 * randNum01[k++];
            //double d_phi   = PI*0.1 * randNum01[k++];
            //theta    = theta_prev + d_theta;
            //phi      = phi_prev   + d_phi;

            if (it-vecNode.begin() == 0)
            {
                double r = rInit * cbrt(randNum01[k++]);
                coord = {r*sin(theta)*cos(phi), r*sin(theta)*sin(phi), r*cos(theta)};
            }
            else
            {
                coord = {x_prev+l*sin(theta)*cos(phi), y_prev+l*sin(theta)*sin(phi), z_prev+l*cos(theta)};
                if (coord.x*coord.x+coord.y*coord.y+coord.z*coord.z >= rInit*rInit) // if the chain extrude from allowed scope
                    coord = {x_prev-l*sin(theta)*cos(phi), y_prev-l*sin(theta)*sin(phi), z_prev-l*cos(theta)};
            }
            it->set_coord(coord);
            x_prev = coord.x;
            y_prev = coord.y;
            z_prev = coord.z;
            theta_prev = theta;
            phi_prev = phi;
        }

    }
    else if (initType == "randomChainBinned")
    {
        /*
            (1) put centromere at (0, 0, RAD_NUCLEUS)
            (2) calculate distance per bead according to total centromere-telomere distance
            (3) every growth step, pick randomly a phi
        */
        double l = L0_TENSION, rMax = l*0.5;
        double theta_cent = 0, theta_telo = 2.*PI/3;    // assume phi_telo = 0
        double x_cent = RAD_NUCLEUS*sin(theta_cent), y_cent = 0, z_cent = RAD_NUCLEUS*cos(theta_cent);
        double x_telo = RAD_NUCLEUS*sin(theta_telo), y_telo = 0, z_telo = RAD_NUCLEUS*cos(theta_telo);

        int nBead = vecNode.size();
        vecNode[0].set_coord({x_telo, y_telo, z_telo});       // telomere
        vecNode[nBead-1].set_coord({x_cent, y_cent, z_cent}); // centromere

        double x_prev = x_cent, y_prev = y_cent, z_prev = z_cent;

        double Ltot = z_cent-z_telo, step_z = Ltot/nBead;
        double r = sqrt(l*l - step_z*step_z), r2;

        string typeShift = "linear";
        typeShift = "squared";
        typeShift = "squared_sine";

        dVec shiftXYZ;
        if (typeShift == "linear")
            shiftXYZ = {(x_telo-x_cent)/nBead, (y_telo-y_cent)/nBead, 0};
        if (typeShift == "squared")
            shiftXYZ = {(x_telo-x_cent)/nBead/nBead, (y_telo-y_cent)/nBead/nBead, 0};
        if (typeShift == "squared_sine")
            shiftXYZ = {(x_telo-x_cent)/nBead/nBead, 0, 0};

        double x, y, z, phi;
        double noise_x, noise_y, noise_size = 10, buffer = 25;

        for (vector<Node>::iterator it = vecNode.end(); it != vecNode.begin(); it --) // assume centromere is the last bead
        {
            if (it-vecNode.end() == 0)
            {
                continue;
            }
            else
            {
                phi      = 2.* PI * randNum01[k++];
                r2       = rMax * randNum01[k++];

                noise_x  = (2 * randNum01[k++] - 1) * noise_size;
                noise_y  = (2 * randNum01[k++] - 1) * noise_size;

                if (typeShift == "linear")
                {
                    x = r2*cos(phi) + shiftXYZ.x*(vecNode.end() - it);
                    y = r2*sin(phi) + shiftXYZ.y*(vecNode.end() - it);
                }
                if (typeShift == "squared")
                {
                    x = r2*cos(phi) + shiftXYZ.x*(vecNode.end() - it)*(vecNode.end() - it);
                    y = r2*sin(phi) + shiftXYZ.y*(vecNode.end() - it)*(vecNode.end() - it);
                }
                if (typeShift == "squared_sine")
                {
                    x = r2*cos(phi) + shiftXYZ.x*(vecNode.end() - it)*(vecNode.end() - it);
                    int N_PERIOD = 3;
                    double AMP = 0.9 * sqrt(RAD_NUCLEUS*RAD_NUCLEUS - z_prev*z_prev );
                    y = r2*sin(phi) + AMP*sin( (float)(vecNode.end() - it)*N_PERIOD/nBead *2*PI );
                }

                if (abs(x) > RAD_NUCLEUS-buffer)
                    noise_x = 0;
                if (abs(y) > RAD_NUCLEUS-buffer)
                    noise_y = 0;
                if (abs(z) > RAD_NUCLEUS-buffer)
                {
                    noise_x = 0;
                    noise_y = 0;
                }

                x += noise_x;
                y += noise_y;
                z = z_prev - step_z;

                dVec coord = {x, y, z};
                it->set_coord(coord);
                x_prev = coord.x;
                y_prev = coord.y;
                z_prev = coord.z;
            }

        }
    }
    else if (initType == "randomChainBinnedCylinder")
    {
        /*
            (1) put centromere at (0, 0, RAD_NUCLEUS)
            (2) calculate radius of cylinder by accounting for volume fraction
            (3) calculate distance per bead according to total centromere-telomere distance
            (4) every growth step, pick randomly a phi
        */

        /*
        spherical cone: v = 1/3 * PI * R^3 * sin(theta)^2 * cos(theta)
        spherical cap : v = 1/3 * PI * R^3 * (2 - 3*cos(theta) + cos(theta)^3 )
        ----------- volume of nucleus excluding nucleolus -----------
            Type I : concave interface
            volume = 1 cone vol + 2 small cone cap vol + 1 big cone cap   !!!!!!!! THIS IS WRONG !!!!!!
        */

        double cos_theta1 = RAD_NUCLEOLUS / RAD_NUCLEUS * 0.5;
        double theta1     = acos(cos_theta1);
        double theta2	  = 2*theta1;
        double cos_theta2 = cos(theta2);
        //double sin_theta1 = sqrt(1 - cos_theta1 * cos_theta1);
        //double cos_theta2 = sin_theta1;
        
        double A_factor = 3.0; //scaling factor
        double vol_sphere = 4*PI/3. * RAD_NUCLEUS*RAD_NUCLEUS*RAD_NUCLEUS;
        double vol_cap    = PI/3. * RAD_NUCLEOLUS*RAD_NUCLEOLUS*RAD_NUCLEOLUS * (2 - 3*cos_theta1 + cos_theta1*cos_theta1*cos_theta1);
        double vol_cap_s  = PI/3. * RAD_NUCLEUS*RAD_NUCLEUS*RAD_NUCLEUS * (2 - 3*cos_theta2 + cos_theta2*cos_theta2*cos_theta2);
        double vol_total = vol_sphere - (vol_cap_s - vol_cap);
        double vol_occupy = vol_total * VOL_FRAC_OF_GENOME * A_factor;
        //double vol_occupy = vol_total * VOL_FRAC_OF_GENOME;

        //double vol_cone   = PI/3. * RAD_NUCLEOLUS*RAD_NUCLEOLUS*RAD_NUCLEOLUS * sin_theta1*sin_theta1 * cos_theta1;
        //double vol_total  = vol_cone + vol_cap + 2*vol_cap_s;
        //double vol_occupy = vol_total * VOL_FRAC_OF_GENOME;
        
        double Ltot       = RAD_NUCLEOLUS;
        double rMax       = sqrt(vol_occupy / PI / Ltot); // this is a cylinder with vol fraction equal to genome size fraction
        //double rMax       = 0.95*RAD_NUCLEUS;
        //double rMax       = 110; //Yasu's data [nm]

        cout << "... whole available space excluding nucleolus      (um^3): " << vol_total  *1E-9  << endl;
        cout << "... initializing chromosome in a cylindical volume (um^3): " << vol_occupy *1E-9 << endl;
        cout << "      Length = " << Ltot << " ; Radius = " << rMax << endl;

        double l = L0_TENSION;
        double x_cent = 0, y_cent = 0, z_cent = RAD_NUCLEUS;
        double x_telo = 0, y_telo = 0, z_telo = RAD_NUCLEUS-RAD_NUCLEOLUS;

        int nBead = vecNode.size();

        double x_prev = x_cent, y_prev = y_cent, z_prev = z_cent;
        double phi_prev = 0;

        double step_z = Ltot/nBead, step_xy = sqrt(l*l - step_z*step_z);

        double x, y, z, phi;
        dVec coord;

        for (vector<Node>::iterator it = vecNode.end(); it != vecNode.begin(); it --) // assume centromere is the last bead
        {
            // random chain
            //theta    = acos(2.* randNum01[k++] - 1);
            //phi      = 2.* PI * randNum01[k++];

            // "persistent" random chain
            double d_phi   = PI*0.1 * (2*randNum01[k++]-1);
            phi      = phi_prev   + d_phi;

            if (it-vecNode.end() == 0)
            {
                coord = {x_cent, y_cent, z_cent};
            }
            else
            {
                coord = {x_prev+step_xy*cos(phi), y_prev+step_xy*sin(phi), z_prev-step_z};
                if ( (coord.x*coord.x+coord.y*coord.y+coord.z*coord.z >= RAD_NUCLEUS*RAD_NUCLEUS)
                    || (coord.x*coord.x+coord.y*coord.y >= rMax*rMax) )
                {
                    //coord = {x_prev-step_xy*cos(phi), y_prev-step_xy*sin(phi), z_prev-step_z};
                    phi += PI;
                    coord = {x_prev+step_xy*cos(phi), y_prev+step_xy*sin(phi), z_prev-step_z};
                }

            }
            it->set_coord(coord);
            x_prev = coord.x;
            y_prev = coord.y;
            z_prev = coord.z;
            phi_prev = phi;
            //cout << z_prev << endl;
        }

        //cout << z_telo << endl;

        vecNode[0].set_coord({x_prev, y_prev, z_telo});
    }
    else if (initType == "randomChainBinaryInsertion")
    {
        /*
            (1) put centromere at (0, 0, RAD_NUCLEUS)
            (2) put telomere at (RAD_NUCLEUS*sin(theta)*cos(phi), RAD_NUCLEUS*sin(theta)*sin(phi), RAD_NUCLEUS*cos(theta) )  theta = PI*2/3,
            (3) introduce a new bead in the middle of the chain recursively
        */
        bool flagFluoroConstraint = true;

        unordered_map<int, double> mapInterFluoroDist;
        if (flagFluoroConstraint == true)   // Petrova2013 Figure 1C
        {
            cout << "flagFluoroConstraint = TRUE" << endl;
            mapInterFluoroDist[vecFluoroSite[1]] = 2.0;
            mapInterFluoroDist[vecFluoroSite[2]] = 1.84;
            mapInterFluoroDist[vecFluoroSite[3]] = 1.67;
            mapInterFluoroDist[vecFluoroSite[4]] = 1.5;
            mapInterFluoroDist[vecFluoroSite[5]] = 0.5;
            mapInterFluoroDist[vecFluoroSite[6]] = 0.996;  // Yasu's paper chromosome I 1.95Mb LacO, chromosome II 3.66Mb
        }
        else
            cout << "flagFluoroConstraint = FALSE" << endl;

        int nBead = NUM_BEAD, nBeadCreated = 0;
        double theta_cent = 0, theta_telo = 2.*PI/3;    // assume phi_telo = 0
        double x_cent = RAD_NUCLEUS*sin(theta_cent), y_cent = 0, z_cent = RAD_NUCLEUS*cos(theta_cent);
        double x_telo = RAD_NUCLEUS*sin(theta_telo), y_telo = 0, z_telo = RAD_NUCLEUS*cos(theta_telo);

        vector<dVec> vecCoord, vecCoord_tmp;
        vecCoord.push_back({x_telo, y_telo, z_telo});
        vecCoord.push_back({x_cent, y_cent, z_cent});
        nBeadCreated = 2;
        double sep = RAD_NUCLEUS, fracShorten = 0.625;    // 0.8 too spread/diffuse; 0.6 too compact
        sep = RAD_NUCLEUS;
        fracShorten = 0.64;
        cout << "Initial spacing for insertion = " << sep << endl;
        cout << "fracShorten = " << fracShorten << endl;

        while (nBeadCreated < nBead)
        {
            vecCoord_tmp.clear();
            vecCoord_tmp.push_back(vecCoord[0]);
            bool flagDoneCreation = false;
            for (int k = 0; k < vecCoord.size()-1; k ++)
            {
                if (flagDoneCreation == false)
                {
                    double xm = (vecCoord[k].x + vecCoord[k+1].x) * 0.5;
                    double ym = (vecCoord[k].y + vecCoord[k+1].y) * 0.5;
                    double zm = (vecCoord[k].z + vecCoord[k+1].z) * 0.5;

                    bool flagOutsideNucleus = true;
                    double dx, dy, dz;
                    while (flagOutsideNucleus)
                    {
                        dx = (2*random()/(float)RAND_MAX-1) * sep;
                        dy = (2*random()/(float)RAND_MAX-1) * sep;
                        dz = (2*random()/(float)RAND_MAX-1) * sep;
                        if (  (xm+dx)*(xm+dx) + (ym+dy)*(ym+dy) + (zm+dz)*(zm+dz) >= RAD_NUCLEUS*RAD_NUCLEUS )
                            flagOutsideNucleus = true;
                        else
                            flagOutsideNucleus = false;
                    }

                    vecCoord_tmp.push_back({xm+dx, ym+dy, zm+dz});
                }

                vecCoord_tmp.push_back(vecCoord[k+1]);
                nBeadCreated ++;
                if (nBeadCreated >= nBead)
                    flagDoneCreation = true;
            }

            sep *= fracShorten;
            vecCoord.clear();
            vecCoord = vecCoord_tmp;
        }
        cout << vecCoord.size() << " " << nBead << endl;
        for (vector<dVec>::iterator it = vecCoord.begin(); it != vecCoord.end(); it ++)
        {
            dVec coord = *it;
            int beadId = it-vecCoord.begin();

            if (flagFluoroConstraint == true)
            {
                if ( beadId != 0 && beadId != nBead-1 && find(vecFluoroSite.begin(), vecFluoroSite.end(), beadId) != vecFluoroSite.end() ) // re-assign the position
                {
                    // find the closest point on the spherical surface
                    cout << beadId << endl;
                    dVec coord_tmp;
                    double dist2cent = sqrt(  (coord.x-x_cent)*(coord.x-x_cent) +
                                              (coord.y-y_cent)*(coord.y-y_cent) +
                                              (coord.z-z_cent)*(coord.z-z_cent)   ) / 1000.;    // um
                    double dist2cent_exp = mapInterFluoroDist[beadId];

                    coord_tmp = { x_cent + (coord.x-x_cent)*dist2cent_exp/dist2cent,
                                  y_cent + (coord.y-y_cent)*dist2cent_exp/dist2cent,
                                  z_cent + (coord.z-z_cent)*dist2cent_exp/dist2cent};

                    // check if the coord is outside nucleus
                    bool flagOutsideNucleus = true;
                    double dist2cent_exp_flex = dist2cent_exp;                          // um
                    while (flagOutsideNucleus)
                    {
                        if ( coord_tmp.x*coord_tmp.x + coord_tmp.y*coord_tmp.y + coord_tmp.z*coord_tmp.z < RAD_NUCLEUS*RAD_NUCLEUS )
                            flagOutsideNucleus = false;
                        else
                        {
                            dist2cent_exp_flex -= RAD_BEAD / 1000.;
                            coord_tmp = { x_cent + (coord.x-x_cent)*dist2cent_exp_flex/dist2cent,
                                          y_cent + (coord.y-y_cent)*dist2cent_exp_flex/dist2cent,
                                          z_cent + (coord.z-z_cent)*dist2cent_exp_flex/dist2cent};
                        }
                    }

                    cout << dist2cent << " (exp: " << dist2cent_exp << "; exp_flex: " << dist2cent_exp_flex << ")" << endl;
                    coord = coord_tmp;
                }
            }

            vecNode[beadId].set_coord(coord);
        }
    }
    else if (initType == "randomChainUniformTether")
    {
        /*
            (1) put centromere at (0, 0, RAD_NUCLEUS)
            (2) put telomere at (RAD_NUCLEUS*sin(theta)*cos(phi), RAD_NUCLEUS*sin(theta)*sin(phi), RAD_NUCLEUS*cos(theta) )  theta = PI*2/3,
            (3) uniformly initialize beads
            (4) rank all coordinates by z and assign to beads
        */
        int nBead = NUM_BEAD;
        double theta_cent = 0, theta_telo = 2.*PI/3;    // assume phi_telo = 0
        double x_cent = RAD_NUCLEUS*sin(theta_cent), y_cent = 0, z_cent = RAD_NUCLEUS*cos(theta_cent);
        double x_telo = RAD_NUCLEUS*sin(theta_telo), y_telo = 0, z_telo = RAD_NUCLEUS*cos(theta_telo);

        vecNode[0].set_coord({x_telo, y_telo, z_telo});       // telomere
        vecNode[nBead-1].set_coord({x_cent, y_cent, z_cent}); // centromere

        vector<dVec> vecCoord, vecCoord_tmp;  // from SMALL to LARGE z

        double rInit = 0.99*RAD_NUCLEUS;
        double r, theta, phi;
        int cnt = 0;
        while (vecCoord.size() < nBead - 2)
        {
            double xx, yy, zz, zz_eff;

            r     = rInit * cbrt(randNum01[k++]);
            theta = acos(2.* randNum01[k++] - 1);
            phi   = 2.* PI * randNum01[k++];

            xx = r*sin(theta)*cos(phi);
            yy = r*sin(theta)*sin(phi);
            zz = r*cos(theta);
            zz_eff = r*cos(theta+PI/2-theta_telo/2);

            if (k >= randNum01.size()-1)
            {
                cout << "random numbers depleted! " << endl;
                exit(4);
            }

            // in the future, account for nucleolus!
            //if (zz > z_telo)
            if (true)
            {
                if (vecCoord.size() == 1)
                {
                    dVec coord0 = vecCoord[0];
                    double r0 = sqrt( coord0.x*coord0.x + coord0.y*coord0.y + coord0.z*coord0.z );
                    double theta0 = acos(coord0.z / r0);
                    double z0_eff = r0*cos(theta0+PI/2-theta_telo/2);

                    //if (zz >= value.z)            // large zz in the end
                    if (zz_eff >= z0_eff)
                        vecCoord.push_back({xx,yy,zz});
                    //if (zz < value.z)                            // small zz in the beginning
                    if (zz_eff < z0_eff)
                        vecCoord.insert(vecCoord.begin(), {xx,yy,zz});
                }
                if (vecCoord.size() > 1)    // SORT by z' = r*cos(theta_telo/2 + theta)
                {
                    vecCoord_tmp.clear();
                    bool flagIncludedNow = false;
                    for (int p = 0; p < vecCoord.size(); p ++)
                    {
                        dVec coord0 = vecCoord[p];
                        double r0 = sqrt( coord0.x*coord0.x + coord0.y*coord0.y + coord0.z*coord0.z );
                        double theta0 = acos(coord0.z / r0);
                        double z0_eff = r0*cos(theta0+PI/2-theta_telo/2);

                        //if (coord0.z > zz && flagIncludedNow == false)
                        if (z0_eff > zz_eff && flagIncludedNow == false)
                        {
                            vecCoord_tmp.push_back({xx, yy, zz});
                            flagIncludedNow = true;
                        }
                        vecCoord_tmp.push_back(coord0);
                        //cout << coord0.z << "\t";
                    }
                    //cout << endl;
                    vecCoord.clear();
                    vecCoord = vecCoord_tmp;

                }
                if (vecCoord.size() == 0)
                    vecCoord.push_back({xx,yy,zz});
            }
        }
        cout << vecCoord.size() << " " << nBead << endl;

        for (vector<dVec>::iterator it = vecCoord.begin(); it != vecCoord.end(); it ++)
        {
            dVec coord = *it;
            int beadId = 1+(it-vecCoord.begin());
            vecNode[beadId].set_coord(coord);
        }
    }
    /*cout << "    done!" << endl;

    cout << ">> initializing bead accelerations ... " << endl;
    for (vector<Node>::iterator it = vecNode.begin(); it != vecNode.end(); it ++)
    {
        dVec accel = {0., 0., 0.}, accelPrev = {0., 0., 0.};
        it->set_accel(accel);
        it->set_accelPrev(accelPrev);
    }
    cout << "    done!" << endl;
    */
//------------------------------------------------------------------------------------------------------------------------------   
    else if (initType == "loadPDB")
    {
        int nBead = NUM_BEAD;
        //vecNode[0].set_coord({x_telo, y_telo, z_telo});       // telomere
        //vecNode[nBead-1].set_coord({x_cent, y_cent, z_cent}); // centromere

        vector<dVec> vecCoord; 
		string iFilename_PDB;
		if (argc == 1)
		{
			cout << "\nPlease input PDB filename (without suffix) as argument! \n";
			exit(0);
		}
        if (argc == 2)
        {
            iFilename_PDB = argv[1]+SUFFIX_PDB;
            cout << ">> reading PDB file " << iFilename_PDB << " ... "<< endl;
            ifstream file_PDB (iFilename_PDB);
            cout << " ifstream " << endl; // 1.checkpoint
            
            string line_PDB;
            
            if (file_PDB.is_open())
            {
                int cnt = 0;
                //getline (file_PDB, line_PDB); // skip the first line
                while (getline (file_PDB, line_PDB))
                {
                    stringstream ss (line_PDB);
                    cout << line_PDB << endl;  // 2.checkpoint
                                        
                    vector<string> items;
                    string buf;
                    while (ss >> buf)
                    {
                        items.push_back(buf);
                        
						//cout << "3" << endl;
						/*
                        for (vector<string>::const_iterator i = items.begin(); i != items.end(); ++i) // 4.checkpoint
                        {
							cout << *i << ' ';
                        }
                        */
                        //cout << "3" << endl;
						cout << items.size() << endl;
						cout << buf << endl;
					}
					//cout << items[0] << " " << items[1] << " " << items[2] << " 4.checkpoint" << endl;
						

                    double X_coord, Y_coord, Z_coord;
                    int NUM_LINE;
                    X_coord = stod(items[0]); //7th column in a PDB file
					//cout << "5.checkpoint " << X_coord << endl;
                        
                    Y_coord   = stoi(items[1]);
                    Z_coord = stoi(items[2]);

                    vecCoord.push_back({X_coord,Y_coord,Z_coord});
                        // cout << "... beadId = bead/file line = " << NUM_LINE << " ..and coordinates " << vecCoord[cnt]  << endl;
                    cnt ++;     
                        
                    //exit(124);
                        //cout << "total number : " << cnt << endl;
                }
            }
                //cout << "    done!" << endl; 
				//cout << vecCoord.size() << " = vectorCoord size" << nBead << " = bead number " << endl;
                 
        }
        //cout << vecCoord.size() << " = vectorCoord size" << nBead << " = bead number " << endl;

        for (vector<dVec>::iterator it = vecCoord.begin(); it != vecCoord.end(); it ++)
        {
            dVec coord = *it;
            int beadId = 1+(it-vecCoord.begin());
            vecNode[beadId].set_coord(coord);
        }
	//------------------------------------------------------------------------------------------------------------------------------   

    //cout << "    done!" << endl;

    //cout << ">> initializing bead accelerations ... " << endl;
		for (vector<Node>::iterator it = vecNode.begin(); it != vecNode.end(); it ++)
		{
			dVec accel = {0., 0., 0.}, accelPrev = {0., 0., 0.};
			it->set_accel(accel);
			it->set_accelPrev(accelPrev);
		}
		cout << "    done!" << endl;
    }			
	
}


void initBindersBrownianSimulation(vector<Binder> & vecBinder, vector<Node> & vecNode, vector<int> & vecCondensinSite, \
                                   vector<Chromosome> & vecChromosome)
{
    /* ---------- create Node and Binder ----------- */
    int currBeadId   = vecNode.size();
    int currBinderId = vecBinder.size();
    int nBinder = vecCondensinSite.size();  // assume that each condensin binding site is associated with a binder
    for (int i = 0; i < nBinder; i ++)
    {
        // get coordinate of the base condensin binding site
        int beadIdCBS = vecCondensinSite[i];
        int hostChromoIdCBS = vecNode[beadIdCBS].get_hostChromoId(); // get chromosome id
        dVec coordCBS = vecNode[beadIdCBS].get_coord();
        dVec coordCBS_next = vecNode[beadIdCBS+1].get_coord();
        dVec coordCBS_prev = vecNode[beadIdCBS-1].get_coord();

        vector<int> vecBeadIds, vecAttachBeadIds;
        vecAttachBeadIds.push_back(beadIdCBS);  // assume two feet both attach CBS to start with
        vecAttachBeadIds.push_back(beadIdCBS);

        Node bead_front = Node(currBeadId, -1); // assume the binder doesn't belong to any Chromosome
        // simple coord method
        dVec coord_front = {coordCBS.x*2/5.+coordCBS_next.x*2/5.+coordCBS_prev.x*1/5.,
                            coordCBS.y*2/5.+coordCBS_next.y*2/5.+coordCBS_prev.y*1/5.,
                            coordCBS.z*2/5.+coordCBS_next.z*2/5.+coordCBS_prev.z*1/5.};
        /*
        if (abs(coord_front.x-coordCBS.x) > 20 ||
            abs(coord_front.y-coordCBS.y) > 20 ||
            abs(coord_front.z-coordCBS.z) > 20)
        {
            cout << "coordCBS:    " << coordCBS.x << ", " << coordCBS.y << ", " << coordCBS.z << endl;
            cout << "coord_front: " << coord_front.x << ", " << coord_front.y << ", " << coord_front.z << endl;
        }
        */
        //dVec coord_front = {coordCBS.x+1, coordCBS.y, coordCBS.z};
        bead_front.set_coord(coord_front);
        bead_front.set_veloc({0,0,0});
        vecNode.push_back(bead_front);
        vecBeadIds.push_back(currBeadId);
        currBeadId++;

        Node bead_rear  = Node(currBeadId, -1); // assume the binder doesn't belong to any Chromosome
        // simple coord method
        dVec coord_rear  = {coordCBS.x*2/5.+coordCBS_next.x*1/5.+coordCBS_prev.x*2/5.,
                            coordCBS.y*2/5.+coordCBS_next.y*1/5.+coordCBS_prev.y*2/5.,
                            coordCBS.z*2/5.+coordCBS_next.z*1/5.+coordCBS_prev.z*2/5.};
        //dVec coord_rear = {coordCBS.x-1, coordCBS.y, coordCBS.z};
        bead_rear.set_coord(coord_rear);
        bead_rear.set_veloc({0,0,0});
        vecNode.push_back(bead_rear);
        vecBeadIds.push_back(currBeadId);
        currBeadId++;

        vector<bool> canSlide;
        canSlide.push_back(true);   // front feet
        canSlide.push_back(true);   // rear  feet
        Binder binder0 = Binder(currBinderId, hostChromoIdCBS, canSlide, vecBeadIds, vecAttachBeadIds);
        vecBinder.push_back(binder0);

        // push_back currBinderId to Chromosome
        vector<int> vecBinderIdsInChromo = vecChromosome[hostChromoIdCBS-1].get_vecBinderIds();
        vecBinderIdsInChromo.push_back(currBinderId);
        vecChromosome[hostChromoIdCBS-1].set_vecBinderIds(vecBinderIdsInChromo);

        currBinderId++;
    }
    //  ----- 03-05-2018 ----  'simple' koff rates, changing the amount of condensin molecules by random removal from vecBinder //

    if (DELETION_BINDER_TYPE == "random")
    {
		if (DELETE_BINDER < 1.0)
		{
			int CONC_trivial = DELETE_BINDER*vecBinder.size();        // 0.1 * vecBinder.size();
			//cout << CONC_trivial << "deleted BINDERS!!!!!!!!!!!!!!!!!!!" << endl; //		 number of binders to be deleted
			for (int i = 1; i < CONC_trivial; i++)  
			{
				vecBinder.erase(vecBinder.begin() + rand() % vecBinder.size());    // deletes a random element
			}
		}
		if (DELETE_BINDER == 1.0)
		{
			vecBinder.clear();
		}			
    }
    //vecBinder.clear();
  /*  
    if (DELETION_BINDER_TYPE == "uniform")
    {
		//int CONC_trivial = DELETE_BINDER*vecBinder.size();        // 0.1 * vecBinder.size();
		//cout << CONC_trivial << "deleted BINDERS!!!!!!!!!!!!!!!!!!!" << endl; //		 number of binders to be deleted
		for (int i = 1; i < vecBinder.size(); i++)  
		{
			cout << vecBinderIds(i) << endl;
			if ( i != 0 && i % Nth_REMOVE == 0)
			{
				vecBinder.erase(vecBinder.begin() + (i-1));    // deletes a random element
			}
		}
    }
    */
    int NUM_Binders_before_dynamics = vecBinder.size();
    //cout << NUM_Binders_before_dynamics << endl;
    //

}
// ----------------------------------------------------------------- DYNAMIC LOADING AND UNLOADING ------------------------------------------------------------------------------------------------------------

// ---------------- uploaded from old integration: 31/07/2018

void initBindersLoadingBrownianSimulation(vector<Binder> & vecBinder, vector<Node> & vecNode, vector<int> & vecCondensinSite, \
                                   vector<Chromosome> & vecChromosome, vector<int> & vecLoadedBinders, vector<int> & vecCondensinSiteEmpty, vector<int> & vecBinderIds2)
{
	//vector<int> vecCondensinSiteEmpty;
	vecCondensinSiteEmpty.clear();
	vecLoadedBinders.clear();
	vecBinderIds2.clear();
	
	
	//vector<int> vecBinderIds2;
	for (vector<Binder>::iterator it = vecBinder.begin(); it != vecBinder.end(); it ++)
	{
		//int ID = it->get_id();
		//int ID = it->get_attachChromoId();
		//vecBinderIds2.push_back(ID);
		
		vector<int> vecAttachBeadIds = it->get_vecAttachBeadIds();	
		vecBinderIds2.push_back(vecAttachBeadIds[0]);	
		vecBinderIds2.push_back(vecAttachBeadIds[1]);	
		//cout << "Front" << vecAttachBeadIds[0] << endl;
		//cout << "Rear" << vecAttachBeadIds[1]<< endl;
		/*
		for (int j = 0; j < vecBinderIds2.size(); j ++)
		{
			cout << vecBinderIds2[j] << endl;
		}
		*/
		//vector<int> vecBeadIds = it->get_vecBeadIds();
		//int bIdF = vecBeadIds[0];
		//vecBinderIds2.push_back(vecBeadIds[0]);
		//vecBinderIds2.push_back(vecBeadIds[1]);

	
	}
	
	//exit(0); // to check what part of the code works
	 // if ( false && bla  // to check what part of the code works
	
	//for (vector<CondensinSite>::iterator it = vecCondensinSite.begin(); it != vecCondensinSite.end(); it ++)
	for (int j = 0; j < vecCondensinSite.size(); j ++)
	{
		int CSId = vecCondensinSite[j];
		
        if (find (vecBinderIds2.begin(), vecBinderIds2.end(), (CSId)) == vecBinderIds2.end()) 
		{
			vecCondensinSiteEmpty.push_back(CSId);
			double bId_random = rand() / double(RAND_MAX);
			
			if ( bId_random < P_ACTIVATION_BINDER )
			{
				int BINDER = 1;
				vecLoadedBinders.push_back(BINDER);

				
				int currBeadId   = vecNode.size();
				int currBinderId = vecBinder.size();
				int nBinder = 1;  // 50 % on CondensinSites occupied with a Binder 						// assume that each condensin binding site is associated with a binder
				for (int i = 0; i < nBinder; i ++)
				{
					// get coordinate of the base condensin binding site
					int beadIdCBS = vecCondensinSite[j];
					int hostChromoIdCBS = vecNode[beadIdCBS].get_hostChromoId(); // get chromosome id
					dVec coordCBS = vecNode[beadIdCBS].get_coord();
					dVec coordCBS_next = vecNode[beadIdCBS+1].get_coord();
					dVec coordCBS_prev = vecNode[beadIdCBS-1].get_coord();

					vector<int> vecBeadIds, vecAttachBeadIds;
					vecAttachBeadIds.push_back(beadIdCBS);  // assume two feet both attach CBS to start with
					vecAttachBeadIds.push_back(beadIdCBS);

					Node bead_front = Node(currBeadId, -1); // assume the binder doesn't belong to any Chromosome
					// simple coord method
					dVec coord_front = {coordCBS.x*2/5.+coordCBS_next.x*2/5.+coordCBS_prev.x*1/5.,
										coordCBS.y*2/5.+coordCBS_next.y*2/5.+coordCBS_prev.y*1/5.,
										coordCBS.z*2/5.+coordCBS_next.z*2/5.+coordCBS_prev.z*1/5.};
					/*
					if (abs(coord_front.x-coordCBS.x) > 20 ||
						abs(coord_front.y-coordCBS.y) > 20 ||
						abs(coord_front.z-coordCBS.z) > 20)
					{
						cout << "coordCBS:    " << coordCBS.x << ", " << coordCBS.y << ", " << coordCBS.z << endl;
						cout << "coord_front: " << coord_front.x << ", " << coord_front.y << ", " << coord_front.z << endl;
					}
					*/
					//dVec coord_front = {coordCBS.x+1, coordCBS.y, coordCBS.z};
					bead_front.set_coord(coord_front);
					bead_front.set_veloc({0,0,0});
					vecNode.push_back(bead_front);
					vecBeadIds.push_back(currBeadId);
					currBeadId++;

					Node bead_rear  = Node(currBeadId, -1); // assume the binder doesn't belong to any Chromosome
					// simple coord method
					dVec coord_rear  = {coordCBS.x*2/5.+coordCBS_next.x*1/5.+coordCBS_prev.x*2/5.,
										coordCBS.y*2/5.+coordCBS_next.y*1/5.+coordCBS_prev.y*2/5.,
										coordCBS.z*2/5.+coordCBS_next.z*1/5.+coordCBS_prev.z*2/5.};
					//dVec coord_rear = {coordCBS.x-1, coordCBS.y, coordCBS.z};
					bead_rear.set_coord(coord_rear);
					bead_rear.set_veloc({0,0,0});
					vecNode.push_back(bead_rear);
					vecBeadIds.push_back(currBeadId);
					currBeadId++;

					vector<bool> canSlide;
					canSlide.push_back(true);   // front feet
					canSlide.push_back(true);   // rear  feet
					Binder binder0 = Binder(currBinderId, hostChromoIdCBS, canSlide, vecBeadIds, vecAttachBeadIds);
					vecBinder.push_back(binder0);

					// push_back currBinderId to Chromosome
					vector<int> vecBinderIdsInChromo = vecChromosome[hostChromoIdCBS-1].get_vecBinderIds();
					vecBinderIdsInChromo.push_back(currBinderId);
					vecChromosome[hostChromoIdCBS-1].set_vecBinderIds(vecBinderIdsInChromo);

					currBinderId++;
				}
							
			
		    }
			
		}
		
	}	
	//cout << " Binders LOADED = " << vecLoadedBinders.size() << endl;	
	//cout << "Number of Binders after loading = " << vecBinder.size() << endl;
	//cout << " ----------------" << endl;
}	
		
		
		
		
		//int CondensinSiteId = it->get_id();
		
		/*for (vector<Binder>::iterator it = vecBinder.begin(); it != vecBinder.end(); it ++)
		{
			int BinderId = it->get_id();
			if ( BinderId != CSId )
			{
				vecCondensinSiteEmpty.push_back(CSId);	
			}
		}
		*/
		//created vector of empty CondensinSites
	
	
	/*for (int i = 0; i < vecCondensinSiteEmpty.size(); i ++)
	{
		//bIdEmpty = vecCondensinSiteEmpty[i];
		double bId_random = rand() / double(RAND_MAX);
		if ( bId_random < P_ACTIVATION_BINDER )
		{
			// record the loaded Binder
			int currBeadId   = vecNode.size();
			int currBinderId = vecBinder.size();
			int nBinder = 1;  // 
			vecLoadedBinders.push_back(nBinder);
			
		// get coordinate of the base condensin binding site
			int beadIdCBS = vecCondensinSiteEmpty[i];
			int hostChromoIdCBS = vecNode[beadIdCBS].get_hostChromoId(); // get chromosome id
			dVec coordCBS = vecNode[beadIdCBS].get_coord();
			dVec coordCBS_next = vecNode[beadIdCBS+1].get_coord();
			dVec coordCBS_prev = vecNode[beadIdCBS-1].get_coord();

			vector<int> vecBeadIds, vecAttachBeadIds;
			vecAttachBeadIds.push_back(beadIdCBS);  // assume two feet both attach CBS to start with
			vecAttachBeadIds.push_back(beadIdCBS);

			Node bead_front = Node(currBeadId, -1); // assume the binder doesn't belong to any Chromosome
        // simple coord method
			dVec coord_front = {coordCBS.x*2/5.+coordCBS_next.x*2/5.+coordCBS_prev.x*1/5.,
								coordCBS.y*2/5.+coordCBS_next.y*2/5.+coordCBS_prev.y*1/5.,
								coordCBS.z*2/5.+coordCBS_next.z*2/5.+coordCBS_prev.z*1/5.};

        //dVec coord_front = {coordCBS.x+1, coordCBS.y, coordCBS.z};
			bead_front.set_coord(coord_front);
			bead_front.set_veloc({0,0,0});
			vecNode.push_back(bead_front);
			vecBeadIds.push_back(currBeadId);
			currBeadId++;

			Node bead_rear  = Node(currBeadId, -1); // assume the binder doesn't belong to any Chromosome
			// simple coord method
			dVec coord_rear  = {coordCBS.x*2/5.+coordCBS_next.x*1/5.+coordCBS_prev.x*2/5.,
								coordCBS.y*2/5.+coordCBS_next.y*1/5.+coordCBS_prev.y*2/5.,
								coordCBS.z*2/5.+coordCBS_next.z*1/5.+coordCBS_prev.z*2/5.};
        //dVec coord_rear = {coordCBS.x-1, coordCBS.y, coordCBS.z};
			bead_rear.set_coord(coord_rear);
			bead_rear.set_veloc({0,0,0});
			vecNode.push_back(bead_rear);
			vecBeadIds.push_back(currBeadId);
			currBeadId++;

			vector<bool> canSlide;
			canSlide.push_back(true);   // front feet
			canSlide.push_back(true);   // rear  feet
			Binder binder0 = Binder(currBinderId, hostChromoIdCBS, canSlide, vecBeadIds, vecAttachBeadIds);
			vecBinder.push_back(binder0);

			// push_back currBinderId to Chromosome
			vector<int> vecBinderIdsInChromo = vecChromosome[hostChromoIdCBS-1].get_vecBinderIds();
			vecBinderIdsInChromo.push_back(currBinderId);
			vecChromosome[hostChromoIdCBS-1].set_vecBinderIds(vecBinderIdsInChromo);

			//currBinderId++;
		}
	}
	* /
		// iterated across empty CondensinSites and check their probability for being loaded with a new Binder, if so -> add a new Binder

}
*/
// -----------------------------------------------------------------

// ------ 21.6.2018 ---------------- unloading of Binders
void initBindersUnloadingBrownianSimulation(vector<Binder> & vecBinder, vector<Node> & vecNode, vector<int> & vecCondensinSite, vector<Chromosome> & vecChromosome, vector<int> & vecLoadedBinders)
{
	/* check vecCondensinSite if the position is empty 
	 * if busy, then check probability for Poff
	 * if Poff, then remove the Binder
	 */
	 

	
	int NUM_UNLOADED_BINDERS = vecLoadedBinders.size();                  // example ; vecBinder.size() - 4;              //  vecLoadedBinders.size(); not working!
	//cout << vecLoadedBinders.size() << "NUMBER LOADED" << endl;
	//cout << "Number of Binders = " << vecBinder.size() << " after loading." << endl;
	//cout << " ................UNLOADING..............." << endl;
	//cout << NUM_UNLOADED_BINDERS << " Binders to be unloaded! " << endl;
	
															// 0.1 * vecBinder.size();		 number of binders to be deleted
    for (int i = 1; i <= NUM_UNLOADED_BINDERS; i++)  
    {
    		vecBinder.erase(vecBinder.begin() + rand() % vecBinder.size());    // deletes a random element
    }
    
    //cout << "Number of Binders = " << vecBinder.size() << " after unloading." << endl;
	//cout << " //////////////////////" << endl;
	//cout << " //////////////////////" << endl;
	//cout << " //////////////////////" << endl;


	/*
	//vector<int> vecBinderIds2;
	for (vector<Binder>::iterator it = vecBinder.begin(); it != vecBinder.end(); it ++)
	{
		//int ID = it->get_id();
		//int ID = it->get_attachChromoId();
		//vecBinderIds2.push_back(ID);
		
		vector<int> vecAttachBeadIds = it->get_vecAttachBeadIds();	
		vecBinderIds2.push_back(vecAttachBeadIds[0]);	
		vecBinderIds2.push_back(vecAttachBeadIds[1]);	
		
		//vector<int> vecBeadIds = it->get_vecBeadIds();
		//int bIdF = vecBeadIds[0];
		//vecBinderIds2.push_back(vecBeadIds[0]);
		//vecBinderIds2.push_back(vecBeadIds[1]);
	}
	
	for (int j = 0; j < vecBinder.size(); j ++)
	{	
		double bId_random = rand() / double(RAND_MAX);
			
		if ( bId_random < P_DISSOCIATE_BINDER)
		{
			int BINDER_del = 1;
			vecUnloadedBinders.push_back(BINDER_del);
			vecBinder.erase(vecBinder.begin() + j);
		}
	
	}
	
	*/
	//exit(0); // to check what part of the code works
	 // if ( false && bla  // to check what part of the code works
			
}


// ------------------------------- END 31/07/2018 ------------------------------------


// -------------------------------- COUPLED DYNAMIC LOADING AND UNLOADING OF BINDERS (Mirny) 10/08/2018 ---------------------------------------



// -----------------------------------------------------------------
void printChromoInfo(vector<Chromosome> & vecChromosome, vector<Node> & vecNode, vector<Binder> & vecBinder)
{
    cout << "Chromosome_id\tBead_id\tSpeed_avg(nm/s)\tSpeed_max(nm/s)" << endl;
    for (vector<Chromosome>::iterator it_c = vecChromosome.begin(); it_c != vecChromosome.end(); it_c++)
    {
        vector<int> vecBeadIds = it_c->get_vecBeadIds();
        int chromoId = it_c->get_id();

        int beadIdM;
        dVec coordM, velocM, accelM;
        double speedSqM = -1.;
        double speedSqAvg = 0;
        double speedAvg = 0;
        for (vector<int>::iterator it_n = vecBeadIds.begin(); it_n != vecBeadIds.end(); it_n++)
        {
            int beadId = *it_n;
            dVec veloc = vecNode[beadId].get_veloc();
            double speedSq = veloc.x*veloc.x+veloc.y*veloc.y+veloc.z*veloc.z;
            //speedSqAvg += speedSq;
            speedAvg   += sqrt(speedSq);
            if ( speedSq > speedSqM )
            {
                speedSqM = speedSq;
                beadIdM = beadId;
                coordM = vecNode[beadId].get_coord();
                velocM = veloc;
                accelM = vecNode[beadId].get_accel();
            }
        }

        //cout << chromoId << "\t" << beadIdM << "\t" << sqrt(speedSqAvg/vecBeadIds.size()) << "\t" << sqrt(speedSqM) << endl;
        cout << chromoId << "\t" << beadIdM << "\t" << speedAvg/vecBeadIds.size() << "\t" << sqrt(speedSqM) << endl;
        /*
        cout << "Bead id = " << beadIdM << ": c = {" << coordM.x << ", " << coordM.y << ", " << coordM.z << "}; "
                                        <<  " v = {" << velocM.x << ", " << velocM.y << ", " << velocM.z << "}; "
                                        <<  " a = {" << accelM.x << ", " << accelM.y << ", " << accelM.z << "}; "
                                        << endl;
        */
    }
/*
    cout << "Binder_id\tFront(canSlide?)\tRear_Foot(canSlide?)" << endl;
    for (vector<Binder>::iterator it_b = vecBinder.begin(); it_b != vecBinder.end(); it_b++)
    {
        int binderId = it_b->get_id();
        vector<int> vecAttachBeadIds = it_b->get_vecAttachBeadIds();
        vector<bool> canSlide        = it_b->get_canSlide();
        string printCanSlideF = "(Y)", printCanSlideR = "(Y)";
        if (canSlide[0] == false)
            printCanSlideF = "(N)";
        if (canSlide[1] == false)
            printCanSlideR = "(N)";
        cout << binderId << "\t" << vecAttachBeadIds[0] << printCanSlideF << "\t" << vecAttachBeadIds[1] << printCanSlideR << endl;
    }
*/
}

void writeChromoInfo(vector<Chromosome> & vecChromosome, vector<Node> & vecNode, vector<int> & vecFluoroSite, \
                     vector<int> & vecCondensinSite, vector<int> & vecBinderBlockSite, \
                     vector<Binder> & vecBinder, VoxMap & voxMap, double t_now)
{
    /*
        WARNING: several places need modification for multi-chromosome simulations!!!
    */

    cout << ">> writing bead information at t_now = " << t_now << " ... " << endl;

    // write in x, y, z format (maybe with color information)
    cout << ">> >> coordinate & color ..." << endl;
    ofstream beadFile;
    string beadFileName;
    stringstream ss;
    ss << fixed << setprecision(1) << t_now;
    //beadFileName = "_bead_coord_color_" + ss.str() + ".txt";
    //if (t_now == -T_PREP)
    //    beadFile.open(beadFileName, ios::out | ios::trunc); beadFile.close();

    ofstream binderFile;
    string binderFileName;
    binderFileName = "_binder_attach_site.txt";
    if (t_now == T)
    {
        binderFile.open(binderFileName, ios::out | ios::trunc); binderFile.close();
        binderFile.open(binderFileName, ios::app | ios::binary);
        binderFile << "time(s)\tbinderId\trear\tfront" << endl;
        binderFile.close();
    }

    // write in PDB format
    // ....... format ........ see https://www.ichemlabs.com/166
    // ATOM   {atom_id as global bead id} CA  ASP {chain_id as A, B, C, ...}  {residue_id as bead_id in Chromosome} {x} {y} {z} 1.00 0.00
    // TER  {atom_id_prev + 1}  {empty} ASP {chain_id}  {residue_id}
    bool flagWritePDB = true;
    if (flagWritePDB)   // has an error about large residue number on chain
    {
        /* find out which segments are loops (assuming just one chromosome) */
        vector<int> vecBeadIdInLoop;
        binderFile.open(binderFileName, ios::app | ios::binary);
        
        /*              -------------- old version -----------------------------
		for (vector<Binder>::iterator it_b = vecBinder.begin(); it_b != vecBinder.end(); it_b++)
        {
            vector<int> vecAttachBeadIds = it_b->get_vecAttachBeadIds();
            binderFile <<  t_now << "\t" <<  (it_b-vecBinder.begin()) << "\t"
                       <<  vecAttachBeadIds[1] << "\t"
                       <<  vecAttachBeadIds[0] << endl;
            for (int beadId = vecAttachBeadIds[1]; beadId <= vecAttachBeadIds[0]; beadId++)
            {
                if (find(vecCondensinSite.begin(), vecCondensinSite.end(), beadId) == vecCondensinSite.end()) // not a condensin binding site
                    vecBeadIdInLoop.push_back(beadId);
            }
        }
        binderFile.close();
        
        
        */
        //----------------- 16-08-2018 -----------------------------------------------
        binderFile << t_now << "\t";
        
        for (vector<Binder>::iterator it_b = vecBinder.begin(); it_b != vecBinder.end(); it_b++)
        {
            vector<int> vecAttachBeadIds = it_b->get_vecAttachBeadIds();
            //binderFile <<  t_now << "\t" <<  (it_b-vecBinder.begin()) << "\t" 
            
            
            //binderFile <<  (it_b-vecBinder.begin()) << "\t"                       
            binderFile <<  vecAttachBeadIds[1] << "\t";
            binderFile <<  vecAttachBeadIds[0] << "\t";
            
            for (int beadId = vecAttachBeadIds[1]; beadId <= vecAttachBeadIds[0]; beadId++)
            {
                if (find(vecCondensinSite.begin(), vecCondensinSite.end(), beadId) == vecCondensinSite.end()) // not a condensin binding site
                    vecBeadIdInLoop.push_back(beadId);
            }
            
        }
        binderFile << endl;
        binderFile.close();
        /* ---------------------------------------------------------------- */

        cout << ">> >> PDB file ..." << endl;
        ofstream beadFilePDB;
        string beadFileNamePDB;
        stringstream ssPDB;
        ssPDB << fixed << setprecision(6) << t_now;
        beadFileNamePDB = "_chromoPDB_" + ss.str() + ".pdb";
        //beadFileNamePDB = "PID_"+ PROC_ID_SS.str() + "_chromoPDB_" + ss.str() + ".pdb";
        if (t_now == -T_PREP)
            beadFilePDB.open(beadFileNamePDB, ios::out | ios::trunc); beadFilePDB.close();

        //beadFile.open(beadFileName, ios::app | ios::binary);
        beadFilePDB.open(beadFileNamePDB, ios::app | ios::binary);
        int nChromo = vecChromosome.size();
        unordered_map<int, string> int2letter;
        int2letter[1] = "A"; int2letter[2] = "B"; int2letter[3] = "C"; int2letter[4] = "D"; int2letter[5] = "E"; int2letter[6] = "F";
        int2letter[7] = "G"; int2letter[8] = "H"; int2letter[9] = "I"; int2letter[10] = "J"; int2letter[11] = "K"; int2letter[12] = "L";
        int2letter[13] = "M"; int2letter[14] = "N"; int2letter[15] = "O"; int2letter[16] = "P"; int2letter[17] = "Q"; int2letter[18] = "R";
        int cntAll = 0;
        for (vector<Chromosome>::iterator it_c = vecChromosome.begin(); it_c != vecChromosome.end(); it_c++)
        {
            vector<int> vecBeadIds = it_c->get_vecBeadIds();
            int chromoId = it_c->get_id();
            double cx = (float)chromoId/nChromo, cz = 0, cy = 1-cx;

            int chromoSegLength = 9999; // this is used for visualization purpose because the length is too long for nucleosome level simulation
            int segId = 0;
            for (vector<int>::iterator it_n = vecBeadIds.begin(); it_n != vecBeadIds.end(); it_n++)
            {
                int beadId = *it_n;
                bool flagIsCondensinSite = false, flagIsBinderBlockSite = false;
                bool flagIsInLoop = false;
                bool flagIsFluoroSite = false;
                if (SIM_TYPE == "SIMULATION")
                {
                    if (find(vecFluoroSite.begin(), vecFluoroSite.end(), beadId) != vecFluoroSite.end())
                        flagIsFluoroSite = true;
                    if (find(vecCondensinSite.begin(), vecCondensinSite.end(), beadId) != vecCondensinSite.end())
                        flagIsCondensinSite = true;
                    if (find(vecBinderBlockSite.begin(), vecBinderBlockSite.end(), beadId) != vecBinderBlockSite.end())
                        flagIsBinderBlockSite = true;
                    if (find(vecBeadIdInLoop.begin(), vecBeadIdInLoop.end(), beadId) != vecBeadIdInLoop.end())
                        flagIsInLoop = true;
                }

                dVec coord = vecNode[beadId].get_coord();
                double x = coord.x, y = coord.y, z = coord.z;

                //beadFile <<  x << "\t" <<  y << "\t" <<  z << "\t"
                //         << cx << "\t" << cy << "\t" << cz << endl;

                stringstream ss0, ss1;
                stringstream ssx, ssy, ssz;
                ss0 << setw(5) << beadId+chromoId;
                ss1 << setw(4) << 1+(it_n-vecBeadIds.begin());
                //ssx << fixed << setw(8) << setprecision(3) << x;
                //ssy << fixed << setw(8) << setprecision(3) << y;
                //ssz << fixed << setw(8) << setprecision(3) << z;
                if (x<-1000)
                    ssx << fixed << setw(8) << setprecision(2) << x;
                else
                    ssx << fixed << setw(8) << setprecision(3) << x;
                if (y<-1000)
                    ssy << fixed << setw(8) << setprecision(2) << y;
                else
                    ssy << fixed << setw(8) << setprecision(3) << y;
                if (z<-1000)
                    ssz << fixed << setw(8) << setprecision(2) << z;
                else
                    ssz << fixed << setw(8) << setprecision(3) << z;

                if (SIM_TYPE == "RECONSTRUCT")
                    beadFilePDB << "ATOM  " << ss0.str() << "  CA" << "  " << "ASP" << " " << int2letter[chromoId] << ss1.str() << "    "
                            << ssx.str() << ssy.str() << ssz.str() << "  1.00" << "  0.00" << endl;
                if (SIM_TYPE == "SIMULATION") // Warning: may need to modify for more than one chromosomes!!!
                {
                    int beadRank = (it_n-vecBeadIds.begin());
                    segId = beadRank / chromoSegLength;
                    stringstream ss0_new, ss1_new;
                    ss0_new << setw(5) << beadId+chromoId+segId;
                    ss1_new << setw(4) << 1+beadRank%chromoSegLength;   // starting with 1

                    if ( beadRank >= chromoSegLength && beadRank % chromoSegLength == 0 )  // break the chromosome into "segments"
                    {
                        stringstream ss2;
                        ss2 << setw(5) << beadRank+chromoId;
                        beadFilePDB << "TER   " << ss2.str() << "    " << "  " << "ASP" << " " << int2letter[chromoId] << endl;
                    }

                    if (flagIsFluoroSite == false && flagIsCondensinSite == false && flagIsBinderBlockSite == false && flagIsInLoop == false)
                        beadFilePDB << "ATOM  " << ss0_new.str() << "  CA" << "  " << "ASP" << " " << int2letter[chromoId+segId] << ss1_new.str() << "    "
                                << ssx.str() << ssy.str() << ssz.str() << "  1.00" << "  0.00" << endl;
                    if (flagIsFluoroSite == false && flagIsInLoop == true)
                        beadFilePDB << "ATOM  " << ss0_new.str() << "  CA" << "  " << "HIS" << " " << int2letter[chromoId+segId] << ss1_new.str() << "    "
                                << ssx.str() << ssy.str() << ssz.str() << "  1.00" << "  0.00" << endl;
                    if (flagIsCondensinSite == true)
                        beadFilePDB << "ATOM  " << ss0_new.str() << "  CA" << "  " << "CYS" << " " << int2letter[chromoId+segId] << ss1_new.str() << "    "
                                << ssx.str() << ssy.str() << ssz.str() << "  1.00" << "  0.00" << endl;
                    if (flagIsBinderBlockSite == true)
                        beadFilePDB << "ATOM  " << ss0_new.str() << "  CA" << "  " << "LYS" << " " << int2letter[chromoId+segId] << ss1_new.str() << "    "
                                << ssx.str() << ssy.str() << ssz.str() << "  1.00" << "  0.00" << endl;
                    if (flagIsFluoroSite == true)
                        beadFilePDB << "ATOM  " << ss0_new.str() << "  CA" << "  " << "PHE" << " " << int2letter[chromoId+segId] << ss1_new.str() << "    "
                                << ssx.str() << ssy.str() << ssz.str() << "  1.00" << "  0.00" << endl;
                }
            }
            stringstream ss2;
            cntAll = vecBeadIds[vecBeadIds.size()-1]+1+chromoId;    // used to derive atomId for Binders below
            ss2 << setw(5) << vecBeadIds[vecBeadIds.size()-1]+1+chromoId;
            beadFilePDB << "TER   " << ss2.str() << "    " << "  " << "ASP" << " " << int2letter[chromoId] << endl;
        }
        // .. write all Binders as a separate "chain"
        cntAll ++;
        int cntBinderBead = 1;
        for (vector<Binder>::iterator it_b = vecBinder.begin(); it_b != vecBinder.end(); it_b++)
        {
            vector<int> vecBeadIds = it_b->get_vecBeadIds();
            for (vector<int>::iterator it_n = vecBeadIds.begin(); it_n != vecBeadIds.end(); it_n++)
            {
                int beadId = *it_n;

                dVec coord = vecNode[beadId].get_coord();
                double x = coord.x, y = coord.y, z = coord.z;

                stringstream ss0, ss1;
                stringstream ssx, ssy, ssz;
                ss0 << setw(5) << cntAll++;
                ss1 << setw(4) << cntBinderBead++;
                //ssx << fixed << setw(8) << setprecision(3) << x;
                //ssy << fixed << setw(8) << setprecision(3) << y;
                //ssz << fixed << setw(8) << setprecision(3) << z;

                if (x<-1000)
                    ssx << fixed << setw(8) << setprecision(2) << x;
                else
                    ssx << fixed << setw(8) << setprecision(3) << x;
                if (y<-1000)
                    ssy << fixed << setw(8) << setprecision(2) << y;
                else
                    ssy << fixed << setw(8) << setprecision(3) << y;
                if (z<-1000)
                    ssz << fixed << setw(8) << setprecision(2) << z;
                else
                    ssz << fixed << setw(8) << setprecision(3) << z;

                if (SIM_TYPE == "SIMULATION") // Warning: may need to modify for more than one chromosomes!!!
                {
                    beadFilePDB << "ATOM  " << ss0.str() << "  CA" << "  " << "SER" << " " << int2letter[18] << ss1.str() << "    "
                            << ssx.str() << ssy.str() << ssz.str() << "  1.00" << "  0.00" << endl;
                }
            }
        }
        stringstream ss3;
        ss3 << setw(5) << cntAll++;
        beadFilePDB << "TER   " << ss3.str() << "    " << "  " << "SER" << " " << int2letter[18] << endl;

        //beadFile.close();
        beadFilePDB.close();
    }

    // write VoxMap occupancy
    bool flagWriteVoxMapOcc = true;
    if (flagWriteVoxMapOcc)
    {
        ofstream voxMapOccFile;
        string voxMapOccFileName;
        voxMapOccFileName = "_voxMap_occupancy.txt";

        ofstream voxMapOccFileNEW;
        string voxMapOccFileNameNEW;
        voxMapOccFileNameNEW = "_voxMap_occupancy2new.txt";

        if (t_now == -T_PREP)
        {
            voxMapOccFile.open(voxMapOccFileName, ios::out | ios::trunc); voxMapOccFile.close();
            voxMapOccFile.open(voxMapOccFileName, ios::app | ios::binary);
            //voxMapOccFile << "time(s)\tnBead_in_allVoxBoxIds" << endl;
            voxMapOccFile << "time(s)\tnOccupiedVoxBoxes" << endl;
            voxMapOccFile.close();
            
            voxMapOccFileNEW.open(voxMapOccFileNameNEW, ios::out | ios::trunc); voxMapOccFileNEW.close();
            voxMapOccFileNEW.open(voxMapOccFileNameNEW, ios::app | ios::binary);
            voxMapOccFileNEW << "time(s)\tnBead_in_allVoxBoxIds" << endl;
            //voxMapOccFile << "time(s)\tnOccupiedVoxBoxes" << endl;
            voxMapOccFileNEW.close();
        }
        voxMapOccFile.open(voxMapOccFileName, ios::app | ios::binary);
        voxMapOccFile << t_now << "\t";

        voxMapOccFileNEW.open(voxMapOccFileNameNEW, ios::app | ios::binary);
        voxMapOccFileNEW << t_now << "\t";

        iVec dim = voxMap.get_dim();
        int nMax = dim.i*dim.j*dim.k;
        vector<int> OccupiedVoxBox;
        unordered_map<int, vector<int>> mapBeadIds = voxMap.get_mapBeadIds();
        for (int boxId = 0; boxId < nMax; boxId ++)
        {
            int nBead = 0;
            unordered_map<int, vector<int>>::const_iterator got = mapBeadIds.find(boxId);
            if ( got != mapBeadIds.end() )  // this box has Bead inside
                nBead = mapBeadIds[boxId].size();
                // ------------------ 10-08-2018 ------------------------------------------------
                // only print the voxel with non-zero number of beads inside (nBead != 0)
            if ( nBead > 0 )
            {
				voxMapOccFileNEW << nBead << "\t";
                                cout << nBead << endl;
				OccupiedVoxBox.push_back(nBead);
	    }
        }
        voxMapOccFile << OccupiedVoxBox.size() << endl;
        /*for (int i=0; i < OccupiedVoxBox.size(); i++)
            cout << OccupiedVoxBox[i] << "iterate" << endl;
            cout << "end" << endl;
*/
        voxMapOccFileNEW << endl;
        voxMapOccFile.close();
        voxMapOccFileNEW.close();
        
    }

    // [NOT USED] write pairwise distance
    bool flagWriteDistMat = false;
    if (flagWriteDistMat)
    {
        cout << ">> >> bead pairwise distance ..." << endl;
        ofstream beadDistFile;
        string beadDistFileName;
        stringstream ss3;
        ss3 << fixed << setprecision(1) << t_now;
        beadDistFileName = "_bead_dist_matrix_" + ss3.str() + ".txt";
        if (t_now == 0)
            beadDistFile.open(beadDistFileName, ios::out | ios::trunc); beadDistFile.close();
        beadDistFile.open(beadDistFileName, ios::app | ios::binary);

        int nBeadTotal = vecNode.size();
        double **matBeadDist = NULL;
        matBeadDist = new double* [nBeadTotal];
        for (int i = 0; i<nBeadTotal; i++)
            matBeadDist[i] = new double[nBeadTotal];
        for (int j = 0; j<nBeadTotal; j++)
        {
            for (int k = 0; k<nBeadTotal; k++)
                matBeadDist[j][k] = 0;
        }
        cout << "    ____calculating____" << endl;
        int cnt = 0;
        for (int i = 0; i < nBeadTotal; i ++)
        {
            for (int j = i+1; j < nBeadTotal; j ++)
            {
                dVec coord1 = vecNode[i].get_coord();
                dVec coord2 = vecNode[j].get_coord();
                matBeadDist[i][j] = sqrt( (coord1.x-coord2.x)*(coord1.x-coord2.x) + (coord1.y-coord2.y)*(coord1.y-coord2.y) + (coord1.z-coord2.z)*(coord1.z-coord2.z) );
                matBeadDist[j][i] = matBeadDist[i][j];
                cnt ++;
                if (cnt % ((int) nBeadTotal*nBeadTotal/40) == 0)
                    cout << "     ... " << (int) (cnt*200./nBeadTotal/nBeadTotal) << "\% finished ..." << endl;
            }
        }
        cout << "    ____writing____" << endl;
        int cnt2 = 0;
        for (int i = 0; i < nBeadTotal; i ++)
        {
            for (int j = 0; j < nBeadTotal; j ++)
            {
                beadDistFile <<  matBeadDist[i][j] << "\t";
                cnt2 ++;
                if (cnt2 % ((int) nBeadTotal*nBeadTotal/40) == 0)
                    cout << "     ... " << (int) (cnt2*200./nBeadTotal/nBeadTotal) << "\% finished ..." << endl;
            }
            beadDistFile << endl;
        }
        beadDistFile.close();
    }

    cout << "    done!" << endl;
}

void writeMSD(vector<Node> & vecNode, double t_now)
{
    cout << ">> writing bead MSD information at t_now = " << t_now << " ... " << endl;

    // write in x, y, z format (maybe with color information)
    ofstream beadFile;
    string beadFileName;
    //stringstream ss;
    //ss << fixed << setprecision(3) << t_now;
    beadFileName = "_msd_bead_coord.txt";
    if (t_now == -T_PREP)
    {
        beadFile.open(beadFileName, ios::out | ios::trunc); beadFile.close();
        beadFile.open(beadFileName, ios::app | ios::binary);
        beadFile << "time(s)\tbeadId\tx(nm)\ty(nm)\tz(nm)\tType(bead)" << endl;
        beadFile.close();
    }

    beadFile.open(beadFileName, ios::app | ios::binary);
    for (vector<Node>::iterator it = vecNode.begin(); it != vecNode.end(); it++)
    {
        if (it - vecNode.begin() < NUM_BEAD)    // limit calculations to chromatin beads
        {
            dVec coord = it->get_coord();
            beadFile <<  t_now << "\t" <<  (it-vecNode.begin()) << "\t"
                     <<  coord.x << "\t" <<  coord.y << "\t" <<  coord.z << "\t" << "N" << endl;      // N = Nucleosome
        }
        
        else if (it - vecNode.begin() > NUM_BEAD)    // limit calculations to chromatin beads
        {
            dVec coord = it->get_coord();
            beadFile <<  t_now << "\t" <<  (it-vecNode.begin()) << "\t"
                     <<  coord.x << "\t" <<  coord.y << "\t" <<  coord.z << "\t" << "B" << endl;       // B = Binder
        }
        
    }
    
    beadFile.close();
}

// ---------------------------------------------------------------------
void writeDC2(vector<iPair> & capturer, double t_now)
{
	//cout << ">> writing difussion capturers at t_now = " << t_now << " ... " << endl;
	// write out
	ofstream capturersFile;
	string capturersFileName;
	capturersFileName = "_DC_capturers_index2.txt";
	
	
	/*if (t_now == T)  // is this correct??
	{
		capturersFile.open(capturersFileName, ios::out | ios::trunc); capturersFile.close();
		capturersFile.open(capturersFileName, ios::out | ios::binary);
		capturersFile << "time(s)\tnCapturerId" << endl;
		capturersFile.close();
	}
	*/
	capturersFile.open(capturersFileName, ios::app | ios::binary);
	capturersFile << t_now << "\t";
	
	//for (vector<iPair>::iterator it = capturer.begin(); it != capturer.end(); it++)  // another option how to iterate?

	for (int i=0; i < capturer.size(); i++)
    {
		iPair cap0 = capturer[i];
		int id1 = cap0.i, id2 = cap0.j;
        //capturersFile << it << endl;
        //capturersFile << "writing capturers DC	"; //<< endl;
		capturersFile << id1 << " ";
		capturersFile << id2 << " ";
    }
    capturersFile << endl;
	capturersFile.close();	
}
// ---------------------------------------------------------------------


// -------------------------------------------------------- 03-08-2018 ------------------------------------------

void writeDC2_stats(vector<iPair> & capturer, vector<int> & hit_valence, double t_now)
{
	//cout << ">> writing difussion capturers at t_now = " << t_now << " ... " << endl;
	// write out
	ofstream capturers2File;
	string capturers2FileName;
	capturers2FileName = "_DC_stats_index2.txt";
	
	capturers2File.open(capturers2FileName, ios::app | ios::binary);
	capturers2File << t_now << "\t"; // current time
	capturers2File << capturer.size() << "\t"; // number of DC pairs
	/*

	for (int i=0; i < hit_valence.size(); i++)
    {
		capturers2File << hit_valence[i] << " "; // valence at t_now
    }
    */
    capturers2File << endl;
	capturers2File.close();	
}

//	---------------------------------------------------31-07-2018 and finished 01-08-2018


void writeDynamicsBinders(vector<Binder> & vecBinder, vector<int> & vecLoadedBinders, vector<int> & vecUnloadedBinders, \
                          vector<int> & vecCondensinSiteEmpty2, double t_now, vector<int> & vecUnLoadedBinders)
{
	if ( SEPARATE_DYNAMICS_BINDER == true )
	{
	
		cout << ">> writing information about dynamics Binders at t_now = " << t_now << " ... " << endl;
		// write out
		ofstream DynamicBinderFile;
		string DynamicBinderFileName;
		DynamicBinderFileName = "_dynamic_separate_Binder.txt";
		
		DynamicBinderFile.open(DynamicBinderFileName, ios::app | ios::binary);
		DynamicBinderFile << t_now << "\t";  // 1st column t_now values
		
		//DynamicBinderFile << NUM_Binders_before_dynamics << " ";
		//DynamicBinderFile << vecUnloadedBinders.size() << " "; // Binders to be deleted (UNLOADED)
		DynamicBinderFile << vecLoadedBinders.size()  << " ";  // Number of Binders LOADED
		
		//DynamicBinderFile << (vecBinder.size() - vecLoadedBinders.size() + vecUnloadedBinders.size()) << " "; // Number of Binders after loading and BEFORE unloading
		DynamicBinderFile << vecBinder.size() << " ";  // Number of Binders after all
		//DynamicBinderFile << vecCondensinSiteEmpty2.size() << " "; // Empty sites after all
		
		DynamicBinderFile << endl;
		DynamicBinderFile.close();	
	}
		
	if ( COUPLED_DYNAMICS_BINDER == true )
	{
		cout << ">> writing information about dynamics Binders at t_now = " << t_now << " ... " << endl;
		// write out
		ofstream DynamicBinderFile;
		string DynamicBinderFileName;
		DynamicBinderFileName = "_dynamic_coupled_Binder.txt";
		
		DynamicBinderFile.open(DynamicBinderFileName, ios::app | ios::binary);
		DynamicBinderFile << t_now << "\t";  // 1st column t_now values
		
		//DynamicBinderFile << NUM_Binders_before_dynamics << " ";
		//DynamicBinderFile << vecUnLoadedBinders.size() << " "; // Binders to be deleted (UNLOADED)
		DynamicBinderFile << vecLoadedBinders.size() << " "; //Binders to be loaded ( and then reloaded)
		DynamicBinderFile << vecBinder.size() << " ";  // Number of Binders after all
	
		
		DynamicBinderFile << endl;
		DynamicBinderFile.close();	
	}
	 
}


/*


*/

////////////////////////////////////////////////////////////////////////// 31/07/2018
/*
void writeDC(vector<iPair> & capturer, double t_now)
{
	cout << ">> writing difussion capturers at t_now = " << t_now << " ... " << endl;
	// write out
	ofstream capturersFile;
	string capturersFileName;
	capturersFileName = "_DC_capturers_index.txt";
	
	if (t_now != -T_PREP)  // is this correct??
	{
		capturersFile.open(capturersFileName, ios::out | ios::trunc); capturersFile.close();
		capturersFile.open(capturersFileName, ios::out | ios::binary);
		capturersFile << "time(s)\tnCapturerId" << endl;
		capturersFile.close();
	}
	capturersFile.open(capturersFileName, ios::app | ios::binary);
	capturersFile << t_now << "\t";
	
	//for (vector<iPair>::iterator it = capturer.begin(); it != capturer.end(); it++)  // another option how to iterate?

	for (int i=0; i < capturer.size(); i++)
    {
		//iPair cap0 = capturer[i];
		//int id1 = cap0.i, id2 = cap2.j;
        //capturersFile << it << endl;
		capturersFile << capturer[i].i << " ";
		capturersFile << capturer[i].j << " ";
    }
    cout << endl;
	capturersFile.close();	
}
*/
/////////////////////////////////////////////////////////////////////////
