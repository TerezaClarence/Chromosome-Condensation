/*
    File: dynamics.cpp
    Function: Simulate dynamics of chromosomes and nuclear envelope
    Model: chromoCell
    Created: 11 December, 2017 (XF)
*/

#include "dynamics.hpp"

// ===================== Functions =========================
void oneIter(vector<Chromosome> & vecChromosome, vector<Node> & vecNode, vector<iPair> vecBeadPairInteract,
             vector<double> & vecBeadPairInteractFreq, vector<int> & vecFluoroSite,
             unordered_map<string, double> & energy, VoxMap & voxMap, bool FLAG_INTERACT_ON, vector<double> & randNormNum0, bool flagIsPrep)
{
    /* --- Apply VoxMap algorithm to update the voxMap for detecting collision
    STEP 1 : assign each bead into VoxMap by doing division and modulo
    --- */
    updateVoxMap(voxMap, vecNode);
    //exit(123);

    /* --- Apply Octree algorithm to update the octree for detecting collision
    STEP 1 : assign each bead into Octree Leaf by doing iterative division and modulo
    STEP 2 : check if a certain Leaf is overloaded
    STEP 3 : branch overloaded Leaf and reassign beads in this Leaf to child Leaf
    --- */
    //updateOctree(octree, vecNode, 3);

    /* --- Apply verlet integration (https://en.wikipedia.org/wiki/Verlet_integration#Velocity_Verlet)
    STEP 1 : x(t+dt) = x(t) + v(t)*dt + 0.5*a(t)*dt*dt
    STEP 2 : a(t+dt) updated according to x(t+dt)
    STEP 3 : v(t+dt) = v(t) + 0.5* ( a(t)+a(t+dt) ) * dt
    --- */

    /* -- STEP 1 -- */
    bool flagFreezeFluoro = false;
    bool flagTetherEnds = false;
    updateCoord(vecNode, randNormNum0, flagTetherEnds, flagFreezeFluoro, vecFluoroSite, flagIsPrep);
    //updateCoord_OMP(vecNode);

    /* -- STEP 2 -- */
    updateAccel(vecChromosome, vecNode, vecBeadPairInteract, vecBeadPairInteractFreq, voxMap, energy, FLAG_INTERACT_ON);

    /* -- STEP 3 -- */
    //updateVeloc(vecNode, energy);
    updateVelocLinearDamp(vecNode, energy);
}
void oneIterBrownianSimulation(vector<Chromosome> & vecChromosome, vector<Node> & vecNode,
                               vector<int> & vecCondensinSite, vector<int> & vecFluoroSite,
                               unordered_map<string, double> & energy, VoxMap & voxMap,
                               vector<double> & randNormNum0, bool flagIsPrep)
{
    /* --- Apply VoxMap algorithm to update the voxMap for detecting collision
    STEP 1 : assign each bead into VoxMap by doing division and modulo
    --- */
    updateVoxMap(voxMap, vecNode);

    /* --- Apply verlet integration (https://en.wikipedia.org/wiki/Verlet_integration#Velocity_Verlet)
    STEP 1 : x(t+dt) = x(t) + v(t)*dt + 0.5*a(t)*dt*dt
    STEP 2 : a(t+dt) updated according to x(t+dt)
    STEP 3 : v(t+dt) = v(t) + 0.5* ( a(t)+a(t+dt) ) * dt
    --- */

    /* -- STEP 1 -- */
    bool flagFreezeFluoro = false;
    bool flagTetherEnds = false;
    updateCoord(vecNode, randNormNum0, flagTetherEnds, flagFreezeFluoro, vecFluoroSite, flagIsPrep);

    /* -- STEP 2 -- */
    updateAccelBrownianSimulation(vecChromosome, vecNode, vecCondensinSite, voxMap, energy, randNormNum0);

    /* -- STEP 3 -- */
    //updateVeloc(vecNode, energy);
    updateVelocLinearDamp(vecNode, energy);

}
void oneIterBrownianSimulationOverdamp(vector<Chromosome> & vecChromosome, vector<Node> & vecNode,
                                       vector<int> & vecCondensinSite, vector<int> & vecBinderBlockSite, vector<int> & vecFluoroSite,
                                       vector<Binder> & vecBinder,
                                       unordered_map<string, double> & energy, VoxMap & voxMap,
                                       vector<double> & randNormNum0, bool flagIsPrep, bool flagIsInitConfig, vector<iPair> & capturer, vector<int> & hit_valence, map<int, vector<int>> & Domain_CSB, bool flagLE_ON, bool flagIsInterphase)
{
    /* --- Overdamped approximation
       0 = -e*v(t) + f(t)
    */

    /* --- Apply VoxMap algorithm to update the voxMap for detecting collision
    STEP 1 : assign each bead into VoxMap by doing division and modulo
    --- */
    updateVoxMap(voxMap, vecNode);

    /* --- Apply Euler integration
    STEP 2 : x(t) = dt * v(t)
    STEP 1 : v(t+dt) updated according to x(t)
    */

    /* -- STEP 1 -- */
    bool flagFreezeFluoro = false;
    if (flagIsInitConfig)
        flagFreezeFluoro = true;
    //bool flagTetherEnds = false;
    //if (flagIsPrep)
      //  flagTetherEnds = true;
    bool flagTetherEnds = true;
    updateCoord(vecNode, randNormNum0, flagTetherEnds, flagFreezeFluoro, vecFluoroSite, flagIsPrep);

    /* -- STEP 2 -- */
    updateVelocBrownianSimulationOverdamp(vecChromosome, vecNode, vecCondensinSite, vecBinderBlockSite, vecBinder, voxMap, energy, randNormNum0, flagIsPrep, capturer, hit_valence, Domain_CSB, flagLE_ON, flagIsInterphase);
}

// voxmap
void updateVoxMap(VoxMap & voxMap, vector<Node> & vecNode)
{
    /* indexing i = nz*(dim.i*dim.j) + ny*(dim.k) + nx */
    iVec dim = voxMap.get_dim();
    int nMax = dim.i*dim.j*dim.k;
    double rMax = RAD_NUCLEUS;
    double aVoxX = rMax*2./dim.i;
    double aVoxY = rMax*2./dim.j;
    double aVoxZ = rMax*2./dim.k;

    unordered_map<int, vector<int>> mapBeadIds;
    unordered_map<int, vector<int>> mapBeadIdsHalo;
    for (vector<Node>::iterator it = vecNode.begin(); it != vecNode.end(); it ++)
    {
        int bId = it->get_id();
        dVec coord = it->get_coord();

        double xx, yy, zz;
        int nx, ny, nz;
        xx = (coord.x+rMax) / aVoxX;
        yy = (coord.y+rMax) / aVoxY;
        zz = (coord.z+rMax) / aVoxZ;

        if (xx < 0 || yy < 0 || zz < 0) // warn negative values
        {
            cout << "Please offset coordinates so that x, y, z values are positive !" << endl;
            exit(234);
        }

        nx = (int) floor( xx );
        ny = (int) floor( yy );
        nz = (int) floor( zz );
        int n = nz*(dim.i*dim.j) + ny*dim.i + nx;
        mapBeadIds[n].push_back(bId);

        // OFFSET in space to examine if this Node is within Halo space of neighboring 26 boxes
        int nx2, ny2, nz2;
        double aOffsetX = 2.*RAD_BEAD/aVoxX, aOffsetY = 2.*RAD_BEAD/aVoxY, aOffsetZ = 2.*RAD_BEAD/aVoxZ;
        if (SIM_TYPE == "SIMULATION")
        {
            aOffsetX = LTHR2_REPULSION;
            aOffsetY = LTHR2_REPULSION;
            aOffsetZ = LTHR2_REPULSION;
        }

        // TYPE 1 : 6 faces -------------------------- HALO repaired ----------------------------------
        if (true)
        {
            nx2 = (int) floor( xx+aOffsetX );         // plus x
            if (nx2 == nx+1 && n+1 < nMax)                                    // !!!!!!!!!!!!!!!!!!!!!!!!!! FIXED by Xiao 1_06_2018
                mapBeadIdsHalo[n+1].push_back(bId);
            nx2 = (int) floor( xx-aOffsetX );         // minus x
            if (nx2 == nx-1 && n-1 >= 0)
                mapBeadIdsHalo[n-1].push_back(bId);
            // ---------------------------------------------------
            ny2 = (int) floor( yy+aOffsetY );         // plus y
            if (ny2 == ny+1 && n+dim.i < nMax)
                mapBeadIdsHalo[n+dim.i].push_back(bId);
            ny2 = (int) floor( yy-aOffsetY );         // minus y
            if (ny2 == ny-1 && n-dim.i >= 0)
                mapBeadIdsHalo[n-dim.i].push_back(bId);
            // ---------------------------------------------------
            nz2 = (int) floor( zz+aOffsetZ );         // plus y
            if (nz2 == nz+1 && n+dim.i*dim.j < nMax)
                mapBeadIdsHalo[n+dim.i*dim.j].push_back(bId);
            nz2 = (int) floor( zz-aOffsetZ );         // minus y
            if (nz2 == nz-1 && n-dim.i*dim.j >= 0)
                mapBeadIdsHalo[n-dim.i*dim.j].push_back(bId);
        }        
        // TYPE 2 : 12 edges
        if (true)
        {
            nx2 = (int) floor( xx+aOffsetX );         // plus x, plus y
            ny2 = (int) floor( yy+aOffsetY );
            if (nx2 == nx+1 && ny2 == ny+1 && n+1+dim.i < nMax)
                mapBeadIdsHalo[n+1+dim.i].push_back(bId);
            ny2 = (int) floor( yy-aOffsetY );         // plus x, minus y
            if (nx2 == nx+1 && ny2 == ny-1 && n+1-dim.i >= 0)
                mapBeadIdsHalo[n+1-dim.i].push_back(bId);
            nx2 = (int) floor( xx-aOffsetX );         // minus x, plus y
            ny2 = (int) floor( yy+aOffsetY );
            if (nx2 == nx-1 && ny2 == ny+1 && n-1+dim.i < nMax)
                mapBeadIdsHalo[n-1+dim.i].push_back(bId);
            ny2 = (int) floor( yy-aOffsetY );         // minus x, minus y
            if (nx2 == nx-1 && ny2 == ny-1 && n-1-dim.i >= 0)
                mapBeadIdsHalo[n-1-dim.i].push_back(bId);
            // ----------------------------------------------------------
            nx2 = (int) floor( xx+aOffsetX );         // plus x, plus z
            nz2 = (int) floor( zz+aOffsetZ );
            if (nx2 == nx+1 && nz2 == nz+1 && n+1+dim.i*dim.j < nMax)
                mapBeadIdsHalo[n+1+dim.i*dim.j].push_back(bId);
            nz2 = (int) floor( zz-aOffsetZ );         // plus x, minus z
            if (nx2 == nx+1 && nz2 == nz-1 && n+1-dim.i*dim.j >= 0)
                mapBeadIdsHalo[n+1-dim.i*dim.j].push_back(bId);
            nx2 = (int) floor( xx-aOffsetX );         // minus x, plus z
            nz2 = (int) floor( zz+aOffsetZ );
            if (nx2 == nx-1 && nz2 == nz+1 && n-1+dim.i*dim.j < nMax)
                mapBeadIdsHalo[n-1+dim.i*dim.j].push_back(bId);
            nz2 = (int) floor( zz-aOffsetZ );         // minus x, minus z
            if (nx2 == nx-1 && nz2 == nz-1 && n-1-dim.i*dim.j >= 0)
                mapBeadIdsHalo[n-1-dim.i*dim.j].push_back(bId);
            // ----------------------------------------------------------
            ny2 = (int) floor( yy+aOffsetY );         // plus y, plus z
            nz2 = (int) floor( zz+aOffsetZ );
            if (ny2 == ny+1 && nz2 == nz+1 && n+dim.i+dim.i*dim.j < nMax)
                mapBeadIdsHalo[n+dim.i+dim.i*dim.j].push_back(bId);
            nz2 = (int) floor( zz-aOffsetZ );         // plus y, minus z
            if (ny2 == ny+1 && nz2 == nz-1 && n+dim.i-dim.i*dim.j >= 0)
                mapBeadIdsHalo[n+dim.i-dim.i*dim.j].push_back(bId);
            ny2 = (int) floor( yy-aOffsetY );         // minus y, plus z
            nz2 = (int) floor( zz+aOffsetZ );
            if (ny2 == ny-1 && nz2 == nz+1 && n-dim.i+dim.i*dim.j < nMax)
                mapBeadIdsHalo[n-dim.i+dim.i*dim.j].push_back(bId);
            nz2 = (int) floor( zz-aOffsetZ );         // minus y, minus z
            if (ny2 == ny-1 && nz2 == nz-1 && n-dim.i-dim.i*dim.j >= 0)
                mapBeadIdsHalo[n-dim.i-dim.i*dim.j].push_back(bId);
        }
        // TYPE 3 : 8 corners
        if (true)
        {
            nx2 = (int) floor( xx+aOffsetX );         // plus x, plus y, plus z
            ny2 = (int) floor( yy+aOffsetY );
            nz2 = (int) floor( zz+aOffsetZ );
            if (nx2 == nx+1 && ny2 == ny+1 && nz2 == nz+1 && n+1+dim.i+dim.i*dim.j < nMax)
                mapBeadIdsHalo[n+1+dim.i+dim.i*dim.j].push_back(bId);
            nz2 = (int) floor( zz-aOffsetZ );         // plus x, plus y, minus z
            if (nx2 == nx+1 && ny2 == ny+1 && nz2 == nz-1 && n+1+dim.i-dim.i*dim.j >= 0)
                mapBeadIdsHalo[n+1+dim.i-dim.i*dim.j].push_back(bId);
            ny2 = (int) floor( yy-aOffsetY );         // plus x, minus y, plus z
            nz2 = (int) floor( zz+aOffsetZ );
            if (nx2 == nx+1 && ny2 == ny-1 && nz2 == nz+1 && n+1-dim.i+dim.i*dim.j < nMax)
                mapBeadIdsHalo[n+1-dim.i+dim.i*dim.j].push_back(bId);
            nz2 = (int) floor( zz-aOffsetZ );         // plus x, minus y, minus z
            if (nx2 == nx+1 && ny2 == ny-1 && nz2 == nz-1 && n+1-dim.i-dim.i*dim.j >= 0)
                mapBeadIdsHalo[n+1-dim.i-dim.i*dim.j].push_back(bId);
            // ----------------------------------------------------------
            nx2 = (int) floor( xx-aOffsetX );         // minus x, plus y, plus z
            ny2 = (int) floor( yy+aOffsetY );
            nz2 = (int) floor( zz+aOffsetZ );
            if (nx2 == nx-1 && ny2 == ny+1 && nz2 == nz+1 && n-1+dim.i+dim.i*dim.j < nMax)
                mapBeadIdsHalo[n-1+dim.i+dim.i*dim.j].push_back(bId);
            nz2 = (int) floor( zz-aOffsetZ );         // minus x, plus y, minus z
            if (nx2 == nx-1 && ny2 == ny+1 && nz2 == nz-1 && n-1+dim.i-dim.i*dim.j >= 0)
                mapBeadIdsHalo[n-1+dim.i-dim.i*dim.j].push_back(bId);
            ny2 = (int) floor( yy-aOffsetY );         // minus x, minus y, plus z
            nz2 = (int) floor( zz+aOffsetZ );
            if (nx2 == nx-1 && ny2 == ny-1 && nz2 == nz+1 && n-1-dim.i+dim.i*dim.j < nMax)
                mapBeadIdsHalo[n-1-dim.i+dim.i*dim.j].push_back(bId);
            nz2 = (int) floor( zz-aOffsetZ );         // minus x, minus y, minus z
            if (nx2 == nx-1 && ny2 == ny-1 && nz2 == nz-1 && n-1-dim.i-dim.i*dim.j >= 0)
                mapBeadIdsHalo[n-1-dim.i-dim.i*dim.j].push_back(bId);
        }
    }
    voxMap.set_mapBeadIds(mapBeadIds);
    voxMap.set_mapBeadIdsHalo(mapBeadIdsHalo);

    // printing
    bool flagPrint = false;
    if (flagPrint)
    {
        for (unordered_map<int, vector<int>>::iterator it = mapBeadIds.begin(); it != mapBeadIds.end(); it ++)
        {
            int boxId = it->first;
            vector<int> vecBeadIds = it->second;
            cout << "boxId = " << boxId << " contains beadIds: ";
            for (vector<int>::iterator it2 = vecBeadIds.begin(); it2 != vecBeadIds.end(); it2 ++)
                cout << *it2 << "\t";

            unordered_map<int, vector<int>>::const_iterator got = mapBeadIdsHalo.find(boxId);
            if ( got != mapBeadIdsHalo.end() )
            {
                cout << "with Halo beadIds: ";
                for (vector<int>::iterator it3 = mapBeadIdsHalo[boxId].begin(); it3 != mapBeadIdsHalo[boxId].end(); it3 ++)
                    cout << *it3 << "\t";
            }

            cout << "\n" << endl;
        }
    }
}
// [NOT USED] octree
void updateOctree(Octree & octree, vector<Node> & vecNode, int depth_init)
{
    int depth = octree.get_depth();
    if (depth == 0) // initialize Leaf at first generation
    {
        vector<Leaf> vecLeaf = octree.get_vecLeaf();
        vector<bVec> position_child0;
        position_child0.push_back({false, false, false});   // this indicates child Leaf has depth 1

        branchLeaf(vecLeaf, position_child0);
        octree.set_vecLeaf(vecLeaf);
        octree.set_depth(depth+1);
    }
    if (depth > 0)
    {
        vector<Leaf> vecLeaf = octree.get_vecLeaf();

    }
}
void branchLeaf(vector<Leaf> & vecLeaf, vector<bVec> position_child0)
{
    vector<bVec> positions = OCTREE_POSTIIONS_CHILD;
    vector<bVec> position_child = position_child0;
    vector<int> vecBeadIds;
    bool flagBranched = false;
    // pop last item
    position_child.pop_back();
    for (vector<bVec>::iterator it = positions.begin(); it != positions.end(); it++)
    {
        position_child.push_back(*it);
        Leaf leaf = Leaf(position_child, vecBeadIds, vecLeaf, flagBranched);
        vecLeaf.push_back(leaf);
    }
}
void fillBeads(Octree & octree, vector<Node> & vecNode)
{
    for (vector<Node>::iterator it = vecNode.begin(); it != vecNode.end(); it ++)
    {
        dVec coord = it->get_coord();
        bVec currLeaf;

        double sizeBox = RAD_NUCLEUS;     // this is the size of first Leaf
        bool flagBranched = true;         // this is the first braching from root
        double fracX = coord.x/sizeBox, fracY = coord.y/sizeBox, fracZ = coord.z/sizeBox;


        while (true)
        {
            bool posX = coord2position(fracX);
            bool posY = coord2position(fracY);
            bool posZ = coord2position(fracZ);
        }
        // function ... convert the coordinate to bVec

        // put the bead into Leaf of Octree (iterative division and modulo)
    }
}
bool coord2position(double fracX)
{
    bool posX = false;
    if ( fracX >= 1 )
        posX = true;
    return posX;
}

// dynamics
void updateCoord(vector<Node> & vecNode, vector<double> & randNormNum0, bool flagTetherEnds, bool flagFreezeFluoro, vector<int> & vecFluoroSite, bool flagIsPrep)
{
    int k = 0;
    for (vector<Node>::iterator it = vecNode.begin(); it != vecNode.end(); it++)
    {
        int beadId = it-vecNode.begin();
        dVec coord = it->get_coord(), veloc = it->get_veloc(), accel = it->get_accel();
        dVec coord_new;
        
        double dx_brownian=0, dy_brownian=0, dz_brownian=0;
        if (SIM_TYPE == "SIMULATION")
        {
			double SCALE_FACTOR = 1.0;
			dx_brownian = randNormNum0[k++]*D_BROWNIAN_SQRT*SCALE_FACTOR;  // nm
            dy_brownian = randNormNum0[k++]*D_BROWNIAN_SQRT*SCALE_FACTOR;
            dz_brownian = randNormNum0[k++]*D_BROWNIAN_SQRT*SCALE_FACTOR;
		}
        
        coord_new = { coord.x + veloc.x*DT + 0.5*accel.x*DT*DT + dx_brownian,
                      coord.y + veloc.y*DT + 0.5*accel.y*DT*DT + dy_brownian,
                      coord.z + veloc.z*DT + 0.5*accel.z*DT*DT + dz_brownian};
  /*      
        if (beadId < 10)
        {
			
			cout << ">> beadId = " << beadId << endl;
			cout << "x: " << coord.x << " y: " << coord.y << "z: " << coord.z << endl;
			cout << "new coord: " << coord_new.x << " " << coord_new.y << " " << coord_new.z << endl;

			cout << "VELOCITY x: " << veloc.x << " y: " << veloc.y << "z: " << veloc.z << endl;
			cout << "BRWN coord: " << dx_brownian << " "<< dy_brownian << " "<< dz_brownian <<  " " << endl;
			
		}  
*/
        if (SIM_TYPE == "RECONSTRUCT")
            it->set_coord(coord_new);

        if (SIM_TYPE == "SIMULATION")
        {
			//cout << "RAD_NUCL_sqrt: " << RAD_NUCLEUS*RAD_NUCLEUS << endl;
            //it->set_coord(coord_new);

            // check confinement by rigid nucleus && fix the coord of centromere           
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////  here calculate distance to axis of the cylinder  - but it's pretty rough 
            double dist2centerSq = coord_new.x*coord_new.x + coord_new.y*coord_new.y + coord_new.z*coord_new.z;
            if (dist2centerSq >= RAD_NUCLEUS*RAD_NUCLEUS)
            {
			//	cout << "dist2centerSq = " << dist2centerSq << endl;
				continue;
			}
                
            //else if (it == vecNode.begin() && flagTetherEnds) // telomere
              //  continue;
            
            //else if ((it - vecNode.begin()) == NUM_BEAD-1 && flagTetherEnds) // centromere
				//continue;
          
            else if (flagFreezeFluoro &&  (find(vecFluoroSite.begin(), vecFluoroSite.end(), beadId)) != vecFluoroSite.end() )
                continue;
            
            // 30_07_2018 - rigid boundary of cylinder; if outside, do not update the coordinates
            // else if 
            /*
            else if (flagIsPrep == true)
            {
				//int bId2 = ; // id of each single bead in a cylinder
				
				//----parameters for the cylindrical confinement
				double cos_theta1 = RAD_NUCLEOLUS / RAD_NUCLEUS * 0.5;
				double theta1     = acos(cos_theta1);
				double theta2	  = 2*theta1;
				double cos_theta2 = cos(theta2);
				double A_factor = 3.0; //scaling factor
				double vol_sphere = 4*PI/3. * RAD_NUCLEUS*RAD_NUCLEUS*RAD_NUCLEUS;
				double vol_cap    = PI/3. * RAD_NUCLEOLUS*RAD_NUCLEOLUS*RAD_NUCLEOLUS * (2 - 3*cos_theta1 + cos_theta1*cos_theta1*cos_theta1);
				double vol_cap_s  = PI/3. * RAD_NUCLEUS*RAD_NUCLEUS*RAD_NUCLEUS * (2 - 3*cos_theta2 + cos_theta2*cos_theta2*cos_theta2);
				double vol_total = vol_sphere - (vol_cap_s - vol_cap);
				double vol_occupy = vol_total * VOL_FRAC_OF_GENOME * A_factor;
				//double vol_occupy = vol_total * VOL_FRAC_OF_GENOME;
				double Ltot       = RAD_NUCLEOLUS;
				//double rMax       = sqrt(vol_occupy / PI / Ltot); // this is a cylinder with vol fraction equal to genome size fraction
				double rMax       = 110; // nm based on Yasu

				//----parameters for pulling force
				double RAD_CYLINDER = RAD_NUCLEUS/3; //110; //Yasu's data d=0.22um  //rMax default value   //10.0; // ?? !! need an estimation 
				double RAD_CYLINDER_sq = RAD_CYLINDER*RAD_CYLINDER;
				double k = 10.0*K_TENSION;       // 1.0. *K_TENSION previous default?? !! need an estimation
				//dVec coord1 = vecNode[bId1].get_coord();
				dVec coord0 = {0, 0, coord_new.z};  //  !!!!!!!!! just for simple case, when using initial configuration with binned cylinder !!!!!!!!!!!!!!!!! if tilted, then modify!
				//dVec ORIGIN = coord0;
				double L = dVecDist(coord0, coord_new);
				//double L_sq = coord2.x*coord2.x + coord2.y*coord3.y + coord2.z*coord2.z;
				double L_inv;
				//dVec veloc1 = vecNode[bId1].get_veloc(), veloc1_new;
				dVec unitVec12 = {0,0,0};

				if (L > 0)
				{
					L_inv = 1./L;  // 1/nm
					unitVec12 = {(coord_new.x-coord0.x)*L_inv, (coord_new.y-coord0.y)*L_inv, (coord_new.z-coord0.z)*L_inv};
				}
				//if ( L_sq > RAD_CYLINDER_sq)
				if ( L > RAD_CYLINDER)  // or if ( L_sq > RAD_CYLINDER_sq)
				{
					
					//double v_pull_size = k*(RAD_CYLINDER-L)*gamma_inv*1E-3;
					//dVec veloc1_brw_temp = veloc1_new;
					//dVec v_pull = {v_pull_size*unitVec12.x, v_pull_size*unitVec12.y, v_pull_size*unitVec12.z};
					//veloc1_new = {  veloc1_brw_temp.x + v_pull.x, veloc1_brw_temp.y + v_pull.y, veloc1_brw_temp.z + v_pull.z};	
					//vecNode[bId1].set_veloc(veloc1_new);
					
					continue;
				}
				else
					it->set_coord(coord_new);
			}
            */
            else
                it->set_coord(coord_new);

        }
    }
}
void updateCoord_OMP(vector<Node> & vecNode)    // not as fast as serial code
{
    /*
    int thread_rank, nNodePerThread = vecNode.size() / NUM_THREADS;
    #pragma omp parallel num_threads(NUM_THREADS) \
                         shared(vecNode, nNodePerThread) \
                         private(thread_rank)
    {
        thread_rank  = omp_get_thread_num();
        int thread_count = omp_get_num_threads();
        int nodeIdOffset = thread_rank*nNodePerThread;

        vector<Node>::iterator first = vecNode.begin() + thread_rank*nNodePerThread;
        vector<Node>::iterator last  = vecNode.begin() + thread_rank*nNodePerThread + nNodePerThread;
        if (thread_rank == thread_count-1)
            last = vecNode.end();

        //cout << "thread_rank = " << thread_rank << " handling " << (last-first) << " Beads " << endl;

        for (vector<Node>::iterator it = first; it != last; it++)  // iterate over thread-specific vector
        //for (vector<Node>::iterator it = vecNode.begin(); it != vecNode.end(); it++)
        {
            dVec coord = it->get_coord(), veloc = it->get_veloc(), accel = it->get_accel();
            dVec coord_new;
            coord_new = { coord.x + veloc.x*DT + 0.5*accel.x*DT*DT,
                          coord.y + veloc.y*DT + 0.5*accel.y*DT*DT,
                          coord.z + veloc.z*DT + 0.5*accel.z*DT*DT };
            it->set_coord(coord_new);
        }
    }
    */

    #pragma omp parallel for num_threads(NUM_THREADS) shared(vecNode) schedule(dynamic)
    for (int i = 0; i < vecNode.size(); i++)
    {
        dVec coord = vecNode[i].get_coord(), veloc = vecNode[i].get_veloc(), accel = vecNode[i].get_accel();
        dVec coord_new;
        coord_new = { coord.x + veloc.x*DT + 0.5*accel.x*DT*DT,
                      coord.y + veloc.y*DT + 0.5*accel.y*DT*DT,
                      coord.z + veloc.z*DT + 0.5*accel.z*DT*DT };
        vecNode[i].set_coord(coord_new);
    }
}

// --- following functions are for RECONSTRUCTION objective ---
void updateAccel(vector<Chromosome> & vecChromosome, vector<Node> & vecNode,
                 vector<iPair> vecBeadPairInteract, vector<double> & vecBeadPairInteractFreq,
                 VoxMap & voxMap, unordered_map<string, double> & energy, bool FLAG_INTERACT_ON)
{
    // (0) clear accel values in Nodes
    for (vector<Node>::iterator it = vecNode.begin(); it != vecNode.end(); it++)
    {
        it->set_accelPrev(it->get_accel()); // save previous accel for verlet integration
        it->set_accel({0, 0, 0});
    }

    // (1) beads connected at nearby genomic site
    vector<int> vecBeadIdLocked;
    if (true)   // interaction between beads close in genomic position; note that the Hi-C map includes the nearest interactions!
    {
        for (vector<Chromosome>::iterator it = vecChromosome.begin(); it != vecChromosome.end(); it++)
        {
            vector<int> vecBeadIds = it->get_vecBeadIds();  // temporary
            // bi-node interaction (harmonic spring)
            for (int i = 0; i < vecBeadIds.size()-1; i ++)
            {
                int bId1 = vecBeadIds[i], bId2 = vecBeadIds[i+1];   // consecutive ids of Bead
                if (i == 0)
                    vecBeadIdLocked.push_back(bId1);
                dVec coord1 = vecNode[bId1].get_coord(), coord2 = vecNode[bId2].get_coord();
                dVec accel1 = vecNode[bId1].get_accel(), accel2 = vecNode[bId2].get_accel(), accel1_new, accel2_new;
                double L0 = 2.*RAD_BEAD, L = dVecDist(coord1, coord2), L_inv;   // nm
                double k = K_NEAR;  // nN/nm     (N = kg*m/s/s;    nN = kg*nm/s/s)
                double m_inv = 1./MASS_BEAD;    // 1/kg, ~ 1E+23
                dVec unitVec12;
                if (L > 0)
                    L_inv = 1./L;   // 1/nm
                else if (L == 0)
                    L_inv = 10.;
                unitVec12 = {(coord2.x-coord1.x)*L_inv, (coord2.y-coord1.y)*L_inv, (coord2.z-coord1.z)*L_inv};
                double a = k*(L-L0)*m_inv;  //  nm /s/s
                //cout << "a = " << a << endl;
                double ax = a * unitVec12.x, ay = a * unitVec12.y, az = a * unitVec12.z;

                // !!! Assume that 1st bead of each chromosome is locked in position !!!
                //if (i != 0)
                if (true)
                {
                    accel1_new = {  accel1.x + ax, accel1.y + ay, accel1.z + az };
                    vecNode[bId1].set_accel(accel1_new);
                }
                accel2_new = {  accel2.x - ax, accel2.y - ay, accel2.z - az };
                vecNode[bId2].set_accel(accel2_new);

                double ePot12 = 0.5*k*(L-L0)*(L-L0);
                energy["potential"] += ePot12;
                energy["total"]     += ePot12;
            }
            // tri-node interaction (confine angle)
            if (false)
            {
                for (int i = 0; i < vecBeadIds.size()-2; i ++)
                {
                    int bId1 = vecBeadIds[i], bId2 = vecBeadIds[i+2];   // consecutive ids of Bead
                    if (i == 0)
                        vecBeadIdLocked.push_back(bId1);
                    dVec coord1 = vecNode[bId1].get_coord(), coord2 = vecNode[bId2].get_coord();
                    dVec accel1 = vecNode[bId1].get_accel(), accel2 = vecNode[bId2].get_accel(), accel1_new, accel2_new;
                    double angle = PI*2./3.;
                    double L0 = 2.*RAD_BEAD*sin(angle), L = dVecDist(coord1, coord2), L_inv;   // nm
                    double k = K_NEAR;  // nN/nm     (N = kg*m/s/s;    nN = kg*nm/s/s)
                    double m_inv = 1./MASS_BEAD;    // 1/kg, ~ 1E+23
                    dVec unitVec12;
                    if (L > 0)
                        L_inv = 1./L;   // 1/nm
                    else if (L == 0)
                        L_inv = 10.;
                    unitVec12 = {(coord2.x-coord1.x)*L_inv, (coord2.y-coord1.y)*L_inv, (coord2.z-coord1.z)*L_inv};
                    double a = k*(L-L0)*m_inv;  //  nm /s/s
                    //cout << "a = " << a << endl;
                    double ax = a * unitVec12.x, ay = a * unitVec12.y, az = a * unitVec12.z;

                    // !!! Assume that 1st bead of each chromosome is locked in position !!!
                    //if (i != 0)
                    if (true)
                    {
                        accel1_new = {  accel1.x + ax, accel1.y + ay, accel1.z + az };
                        vecNode[bId1].set_accel(accel1_new);
                    }
                    accel2_new = {  accel2.x - ax, accel2.y - ay, accel2.z - az };
                    vecNode[bId2].set_accel(accel2_new);

                    double ePot12 = 0.5*k*(L-L0)*(L-L0);
                    energy["potential"] += ePot12;
                    energy["total"]     += ePot12;
                }
            }
        }

    }

    // (2) beads interacting according to Hi-C map
    if (FLAG_INTERACT_ON)   // interaction between beads close in physical position but not in genomic position
    {
        double freqMaxHiC_inv = 1./ *max_element(vecBeadPairInteractFreq.begin(), vecBeadPairInteractFreq.end());
        for (vector<iPair>::iterator it = vecBeadPairInteract.begin(); it != vecBeadPairInteract.end(); it++)
        {
            int pos = it-vecBeadPairInteract.begin();
            int bId1 = (*it).i, bId2 = (*it).j;
            int hostChromoId1 = vecNode[bId1].get_hostChromoId();
            int hostChromoId2 = vecNode[bId2].get_hostChromoId();
            if (hostChromoId1 == hostChromoId2) // this limit the following calculations to intra-chromosomal interactions
            {
                //cout << "beadId: " << bId1 << ", " << bId2 << endl;
                //cout << "chromoID: " << hostChromoId1 << ", " << hostChromoId2;
                //cout << " ... " << endl;
                dVec coord1 = vecNode[bId1].get_coord(), coord2 = vecNode[bId2].get_coord();
                dVec accel1 = vecNode[bId1].get_accel(), accel2 = vecNode[bId2].get_accel(), accel1_new, accel2_new;
                double L0 = 2.*RAD_BEAD, L = dVecDist(coord1, coord2), L_inv;   // nm

                double k = K_INTRA;  // nN/nm     (N = kg*m/s/s;    nN = kg*nm/s/s)
                //if (hostChromoId1 != hostChromoId2)
                //    k = K_INTER;
                k = k*freqMaxHiC_inv*vecBeadPairInteractFreq[pos];  // linearly relate k to the frequency in HiC map

                double m_inv = 1./MASS_BEAD;    // 1/kg, ~ 1E+23
                dVec unitVec12;
                if (L > 0)
                    L_inv = 1./L;   // 1/nm
                else if (L == 0)
                    L_inv = 10.;
                unitVec12 = {(coord2.x-coord1.x)*L_inv, (coord2.y-coord1.y)*L_inv, (coord2.z-coord1.z)*L_inv};
                double a = k*(L-L0)*m_inv;  //  nm /s/s
                //cout << "a = " << a << endl;
                double ax = a * unitVec12.x, ay = a * unitVec12.y, az = a * unitVec12.z;

                //auto findLocked = find(vecBeadIdLocked.begin(), vecBeadIdLocked.end(), bId1);   // temporary
                //if (findLocked == vecBeadIdLocked.end())
                if (true)
                {
                    accel1_new = {  accel1.x + ax, accel1.y + ay, accel1.z + az };
                    vecNode[bId1].set_accel(accel1_new);
                }
                //auto findLocked2 = find(vecBeadIdLocked.begin(), vecBeadIdLocked.end(), bId1);  // temporary
                //if (findLocked2 == vecBeadIdLocked.end())
                if (true)
                {
                    accel2_new = {  accel2.x - ax, accel2.y - ay, accel2.z - az };
                    vecNode[bId2].set_accel(accel2_new);
                }

                double ePot12 = 0.5*k*(L-L0)*(L-L0);
                energy["potential"] += ePot12;
                energy["total"]     += ePot12;
            }
        }
        //exit(123);
    }

    // (3) beads volume exclusion / collision detection
    if (true)   // interaction between beads that collide (physically too close)
    {
        iVec dim = voxMap.get_dim();
        int nMax = dim.i*dim.j*dim.k;
        unordered_map<int, vector<int>> mapBeadIds = voxMap.get_mapBeadIds(); // use reference to save memory
        unordered_map<int, vector<int>> mapBeadIdsHalo = voxMap.get_mapBeadIdsHalo();
        for (int boxId = 0; boxId < nMax; boxId ++)
        {
            vector<int> vecBeadIds0, vecBeadIds1;
            unordered_map<int, vector<int>>::const_iterator got = mapBeadIds.find(boxId);
            if ( got != mapBeadIds.end() )  // this box has Bead inside
            {
                for (vector<int>::iterator it = mapBeadIds[boxId].begin(); it != mapBeadIds[boxId].end(); it ++)
                    vecBeadIds0.push_back(*it);
            }
            unordered_map<int, vector<int>>::const_iterator got2 = mapBeadIdsHalo.find(boxId);
            if ( got2 != mapBeadIdsHalo.end() )  // this box has Bead inside
            {
                for (vector<int>::iterator it = mapBeadIdsHalo[boxId].begin(); it != mapBeadIdsHalo[boxId].end(); it ++)
                    vecBeadIds1.push_back(*it);
            }

            // STEP 1 : repulsion between Beads within the Box
            int nBeadInBox = vecBeadIds0.size();
            if (nBeadInBox > 1)
            {
                for (int l = 0; l < nBeadInBox-1; l++)
                {
                    for (int m = l+1; m < nBeadInBox; m++)
                    {
                        int bId1 = vecBeadIds0[l], bId2 = vecBeadIds0[m];
                        dVec coord1 = vecNode[bId1].get_coord(), coord2 = vecNode[bId2].get_coord();
                        double L0 = 2.*RAD_BEAD, L = dVecDist(coord1, coord2);   // nm
                        if (L < L0) // repulsion only when overlapping of Beads occurs
                        {
                            double L_inv;
                            double m_inv = 1./MASS_BEAD;    // 1/kg, ~ 1E+23
                            dVec unitVec12;
                            if (L > 0)
                                L_inv = 1./L;   // 1/nm
                            else if (L == 0)
                                L_inv = 10.;
                            double k = K_INTRA;
                            dVec accel1 = vecNode[bId1].get_accel(), accel2 = vecNode[bId2].get_accel(), accel1_new, accel2_new;
                            unitVec12 = {(coord2.x-coord1.x)*L_inv, (coord2.y-coord1.y)*L_inv, (coord2.z-coord1.z)*L_inv};
                            double a = k*(L-L0)*m_inv;  //  nm /s/s
                            double ax = a * unitVec12.x, ay = a * unitVec12.y, az = a * unitVec12.z;

                            accel1_new = {  accel1.x + ax, accel1.y + ay, accel1.z + az };
                            vecNode[bId1].set_accel(accel1_new);

                            accel2_new = {  accel2.x - ax, accel2.y - ay, accel2.z - az };
                            vecNode[bId2].set_accel(accel2_new);
                        }
                    }
                }
            }
            // STEP 2 : repulsion between one Bead within the Box and another Bead within Halo
            int nBeadInHalo = vecBeadIds1.size();
            if (nBeadInBox > 0 && nBeadInHalo > 0)
            {
                for (int l = 0; l < nBeadInBox-1; l++)
                {
                    for (int m = 0; m < nBeadInHalo; m++)
                    {
                        int bId1 = vecBeadIds0[l], bId2 = vecBeadIds1[m];
                        dVec coord1 = vecNode[bId1].get_coord(), coord2 = vecNode[bId2].get_coord();
                        double L0 = 2.*RAD_BEAD, L = dVecDist(coord1, coord2);   // nm
                        if (L < L0) // repulsion only when overlapping of Beads occurs
                        {
                            double L_inv;
                            double m_inv = 1./MASS_BEAD;    // 1/kg, ~ 1E+23
                            dVec unitVec12;
                            if (L > 0)
                                L_inv = 1./L;   // 1/nm
                            else if (L == 0)
                                L_inv = 10.;
                            double k = K_INTRA;
                            dVec accel1 = vecNode[bId1].get_accel(), accel2 = vecNode[bId2].get_accel(), accel1_new, accel2_new;
                            unitVec12 = {(coord2.x-coord1.x)*L_inv, (coord2.y-coord1.y)*L_inv, (coord2.z-coord1.z)*L_inv};
                            double a = k*(L-L0)*m_inv;  //  nm /s/s
                            double ax = a * unitVec12.x, ay = a * unitVec12.y, az = a * unitVec12.z;

                            accel1_new = {  accel1.x + ax, accel1.y + ay, accel1.z + az };
                            vecNode[bId1].set_accel(accel1_new);
                        }
                    }
                }
            }
        }
    }

    // (4) beads interacting with nuclear membrane

}

// --- following functions are for Brownian SIMULATION objective ---
void updateAccelBrownianSimulation(vector<Chromosome> & vecChromosome, vector<Node> & vecNode, vector<int> & vecCondensinSite,
                                   VoxMap & voxMap, unordered_map<string, double> & energy, vector<double> & randNormNum0)
{
    // (0) clear accel values in Nodes
    for (vector<Node>::iterator it = vecNode.begin(); it != vecNode.end(); it++)
    {
        it->set_accelPrev(it->get_accel()); // save previous accel for verlet integration
        it->set_accel({0, 0, 0});
    }

    // (1) entropic force & tension force (nearest bi-bead) & attraction force (nearest tri-bead)
    int k = 0;  // for using randNormNum0
    double m_inv = 1./MASS_BEAD;
    for (vector<Chromosome>::iterator it = vecChromosome.begin(); it != vecChromosome.end(); it++)
    {
        vector<int> vecBeadIds = it->get_vecBeadIds();
        for (int i = 0; i < vecBeadIds.size(); i ++)
        {
            //double m_inv = 1./MASS_BEAD;    // (1/kg) make sure MASS_BEAD is calculated correctly

            int bId1 = vecBeadIds[i], bId2, bId3;   // consecutive ids of Bead
            dVec coord1 = vecNode[bId1].get_coord(), coord2, coord3;
            dVec accel1 = vecNode[bId1].get_accel(), accel2, accel3;
            dVec accel1_new = accel1, accel2_new, accel3_new;

            /* --------- (1.1) entropic force 1st bead ---------- */
            if (true)
            {
                double theta    = acos(2.* ((double) rand() / (RAND_MAX)) - 1);
                double phi      = 2.* PI * ((double) rand() / (RAND_MAX));
                //double a_ent_size = F_ENTROPIC*m_inv;
                /* this is constant size of entropic force
                double a_ent_size = F_ENTROPIC_UNIT_MASS;                                  // pN/kg
                dVec a_ent = {a_ent_size*sin(theta)*cos(phi), a_ent_size*sin(theta)*sin(phi), a_ent_size*cos(theta)};
                accel1_new = {accel1.x+a_ent.x, accel1.y+a_ent.y, accel1.z+a_ent.z};
                */
                double a_ent_x = randNormNum0[k++]*m_inv*1E-3;
                double a_ent_y = randNormNum0[k++]*m_inv*1E-3;
                double a_ent_z = randNormNum0[k++]*m_inv*1E-3;
                accel1_new = {accel1.x+a_ent_x, accel1.y+a_ent_y, accel1.z+a_ent_z};
            }



////////////////////////////////////////////////////////////////////////
                    /* ------------------- (1.0) repulsion for beads outside a cylider ------------------------------------- */
                    
////////////////////////////////////////////////////////////////////////



            /* --------- (1.2) tension force between 1st & 2nd bead
               &         (1.3) attraction force between 1st & 3rd bead ----------- */
            if (true)
            {
                if (i < vecBeadIds.size()-1)
                {
                    bId2 = vecBeadIds[i+1];
                    coord2 = vecNode[bId2].get_coord();
                    accel2 = vecNode[bId2].get_accel();
                    double L0 = L0_TENSION, L = dVecDist(coord1, coord2);   // nm
                    //double k = K_TENSION;                                               // pN/nm     (N = kg*m/s/s; nN = kg*nm/s/s; pN = kg*pm/s/s)
                    double k = K_TENSION_UNIT_MASS;
                    dVec unitVec12 = {0,0,0};
                    if (L > 0)
                    {
                        double L_inv = 1./L;  // 1/nm
                        unitVec12 = {(coord2.x-coord1.x)*L_inv, (coord2.y-coord1.y)*L_inv, (coord2.z-coord1.z)*L_inv};
                    }

                    if ( abs(L - L0) > L_DIFF_RESOLUTION) // non-zero tension
                    {
                        //double a_ten_size = k*(L-L0)*m_inv;
                        double a_ten_size = k*(L-L0);   // nN/nm/kg
                        //cout << "L = " << L << " , L0 = " << L0 << " , abs(L-L0) = " << abs(L-L0) << endl;
                        //cout << "a_ten_size = " << a_ten_size << endl;
                        dVec a_ten = {a_ten_size*unitVec12.x, a_ten_size*unitVec12.y, a_ten_size*unitVec12.z};
                        dVec accel1_tmp1 = accel1_new;
                        accel1_new = {  accel1_tmp1.x + a_ten.x, accel1_tmp1.y + a_ten.y, accel1_tmp1.z + a_ten.z };
                        accel2_new = {  accel2.x      - a_ten.x, accel2.y      - a_ten.y, accel2.z      - a_ten.z };
                        // update acceleration
                        vecNode[bId2].set_accel(accel2_new);
                        // update potential & total energy
                        double ePot12 = 0.5*k*(L-L0)*(L-L0);
                        energy["potential"] += ePot12;
                        //cout << "potential energy: add = " << ePot12 << " , sum = " << energy["potential"] << endl;
                        energy["total"]     += ePot12;
                    }

                    if (i < vecBeadIds.size()-2)
                    {
                        bId3 = vecBeadIds[i+2];
                        coord3 = vecNode[bId3].get_coord();
                        accel3 = vecNode[bId3].get_accel();
                        double L = dVecDist(coord1, coord3), L_inv;
                        //double k = K_ATTRACTION;
                        double k = K_ATTRACTION_UNIT_MASS;
                        dVec unitVec13 = {0,0,0};
                        if (L > 0)
                        {
                            L_inv = 1./L;
                            unitVec13 = {(coord3.x-coord1.x)*L_inv, (coord3.y-coord1.y)*L_inv, (coord3.z-coord1.z)*L_inv};
                        }
                        //double a_att_size = k*L*m_inv;
                        double a_att_size = k*L;
                        dVec a_att = {a_att_size*unitVec13.x, a_att_size*unitVec13.y, a_att_size*unitVec13.z};
                        dVec accel1_tmp2 = accel1_new;
                        accel1_new = {  accel1_tmp2.x + a_att.x, accel1_tmp2.y + a_att.y, accel1_tmp2.z + a_att.z };
                        accel3_new = {  accel3.x - a_att.x, accel3.y - a_att.y, accel3.z - a_att.z };
                        // update acceleration
                        vecNode[bId3].set_accel(accel3_new);
                        // update potential & total energy
                        double ePot13 = 0.5*k*L*L;
                        energy["potential"] += ePot13;
                        energy["total"]     += ePot13;
                    }
                }
                // update acceleration
                vecNode[bId1].set_accel(accel1_new);
            }
        }
    }

    // (2) repulsion force
    if (true)   // interaction between beads that collide (physically too close)
    {
        iVec dim = voxMap.get_dim();
        int nMax = dim.i*dim.j*dim.k;
        unordered_map<int, vector<int>> mapBeadIds = voxMap.get_mapBeadIds(); // use reference to save memory
        unordered_map<int, vector<int>> mapBeadIdsHalo = voxMap.get_mapBeadIdsHalo();
        for (int boxId = 0; boxId < nMax; boxId ++)
        {
            vector<int> vecBeadIds0, vecBeadIds1;
            unordered_map<int, vector<int>>::const_iterator got = mapBeadIds.find(boxId);
            if ( got != mapBeadIds.end() )  // this box has Bead inside
            {
                for (vector<int>::iterator it = mapBeadIds[boxId].begin(); it != mapBeadIds[boxId].end(); it ++)
                    vecBeadIds0.push_back(*it);
            }
            unordered_map<int, vector<int>>::const_iterator got2 = mapBeadIdsHalo.find(boxId);
            if ( got2 != mapBeadIdsHalo.end() )  // this box has Bead inside
            {
                for (vector<int>::iterator it = mapBeadIdsHalo[boxId].begin(); it != mapBeadIdsHalo[boxId].end(); it ++)
                    vecBeadIds1.push_back(*it);
            }

            // STEP 1 : repulsion between Beads within the Box
            int nBeadInBox = vecBeadIds0.size();
            if (nBeadInBox > 1)
            {
                for (int l = 0; l < nBeadInBox-1; l++)
                {
                    for (int m = l+1; m < nBeadInBox; m++)
                    {
                        int bId1 = vecBeadIds0[l], bId2 = vecBeadIds0[m];
                        dVec coord1 = vecNode[bId1].get_coord(), coord2 = vecNode[bId2].get_coord();
                        double LTHR1 = LTHR1_REPULSION, LTHR2 = LTHR2_REPULSION, L = dVecDist(coord1, coord2);   // nm
                        if (abs(bId1-bId2) >= 2 &&  L < LTHR2) // repulsion only when overlapping of Beads occurs && when two Beads are not connected
                        {
                            double L_inv;
                            //double m_inv = 1./MASS_BEAD;    // 1/kg
                            dVec unitVec12 = {0,0,0};
                            if (L > 0)
                            {
                                L_inv = 1./L;   // 1/nm
                                unitVec12 = {(coord2.x-coord1.x)*L_inv, (coord2.y-coord1.y)*L_inv, (coord2.z-coord1.z)*L_inv};
                            }
                            //double f = F_REPULSION;
                            double f = F_REPULSION_UNIT_MASS;   // nN/kg
                            dVec accel1 = vecNode[bId1].get_accel(), accel2 = vecNode[bId2].get_accel(), accel1_new, accel2_new;

                            //double a_rep_size = f*m_inv*pow(L_inv, 12);
                            double a_rep_size = f;
/*
                            if (L < LTHR1)
                                a_rep_size = f;
                            if (L >= LTHR1)
                                a_rep_size = f*pow(L_inv, 12);
*/
                            dVec a_rep = {a_rep_size*unitVec12.x, a_rep_size*unitVec12.y, a_rep_size*unitVec12.z};

                            accel1_new = {  accel1.x - a_rep.x, accel1.y - a_rep.y, accel1.z - a_rep.z };
                            vecNode[bId1].set_accel(accel1_new);

                            accel2_new = {  accel2.x + a_rep.x, accel2.y + a_rep.y, accel2.z + a_rep.z };
                            vecNode[bId2].set_accel(accel2_new);

                            double ePot12 = f*(LTHR2-L);
/*
                            if (L < LTHR1)
                                ePot12 = -f*L;
                            if (L >= LTHR1)
                                ePot12 = -f*pow(L_inv, 11)/11.;
*/
                            energy["potential"] += ePot12;
                            energy["total"]     += ePot12;
                        }
                    }
                }
            }
            // STEP 2 : repulsion between one Bead within the Box and another Bead within Halo
            int nBeadInHalo = vecBeadIds1.size();
            if (nBeadInBox > 0 && nBeadInHalo > 0)
            {
                for (int l = 0; l < nBeadInBox-1; l++)
                {
                    for (int m = 0; m < nBeadInHalo; m++)
                    {
                        int bId1 = vecBeadIds0[l], bId2 = vecBeadIds1[m];
                        if (bId1 < bId2 && abs(bId1-bId2) >= 2) // bId1 < bId2 : avoid calculating twice! && when two Beads are not connected
                        {
                            dVec coord1 = vecNode[bId1].get_coord(), coord2 = vecNode[bId2].get_coord();
                            double LTHR1 = LTHR1_REPULSION, LTHR2 = LTHR2_REPULSION, L = dVecDist(coord1, coord2);   // nm
                            if (L < LTHR2)
                            {
                                double L_inv;
                                //double m_inv = 1./MASS_BEAD;    // 1/kg
                                dVec unitVec12 = {0,0,0};
                                if (L > 0)
                                {
                                    L_inv = 1./L;   // 1/nm
                                    unitVec12 = {(coord2.x-coord1.x)*L_inv, (coord2.y-coord1.y)*L_inv, (coord2.z-coord1.z)*L_inv};
                                }
                                //double f = F_REPULSION;
                                double f = F_REPULSION_UNIT_MASS;
                                dVec accel1 = vecNode[bId1].get_accel(), accel2 = vecNode[bId2].get_accel(), accel1_new, accel2_new;
                                //double a_rep_size = f*m_inv*pow(L_inv, 12);  //  nm /s/s
                                double a_rep_size = f;
/*
                                if (L < LTHR1)
                                    a_rep_size = f;
                                    //a_rep_size = f*m_inv;
                                if (L >= LTHR1)
                                    a_rep_size = f*pow(L_inv, 12);
*/
                                dVec a_rep = {a_rep_size*unitVec12.x, a_rep_size*unitVec12.y, a_rep_size*unitVec12.z};

                                accel1_new = {  accel1.x - a_rep.x, accel1.y - a_rep.y, accel1.z - a_rep.z };
                                vecNode[bId1].set_accel(accel1_new);

                                accel2_new = {  accel2.x + a_rep.x, accel2.y + a_rep.y, accel2.z + a_rep.z };
                                vecNode[bId2].set_accel(accel2_new);

                                double ePot12 = f*(LTHR2-L);
/*
                                if (L < LTHR1)
                                    ePot12 = -f*L;
                                if (L >= LTHR1)
                                    ePot12 = -f*pow(L_inv, 11)/11.;
*/
                                energy["potential"] += ePot12;
                                energy["total"]     += ePot12;
                            }
                        }
                    }
                }
            }
        }
    }

    // (3) condensin attraction
    if (true)
    {
        vector<int> vecCondensinSiteBusy;   // contain bead ids already interacting with another
        for (int i = 0; i < vecCondensinSite.size(); i ++)
        {
            for (int j = i+1; j < vecCondensinSite.size(); j ++)
            {
                if ( (rand() / double(RAND_MAX)) > P_DISSOCIATE_CONDENSIN )
                {
                    int bId1 = vecCondensinSite[i], bId2 = vecCondensinSite[j];
                    dVec coord1 = vecNode[bId1].get_coord(), coord2 = vecNode[bId2].get_coord();
                    double L = dVecDist(coord1, coord2);
                    if (   L < LTHR_TENSION_CONDENSIN
                        && find(vecCondensinSiteBusy.begin(), vecCondensinSiteBusy.end(), bId1) == vecCondensinSiteBusy.end()
                        && find(vecCondensinSiteBusy.begin(), vecCondensinSiteBusy.end(), bId2) == vecCondensinSiteBusy.end())  // dissociation & avoid more than one contact
                    {
                        vecCondensinSiteBusy.push_back(bId1);
                        vecCondensinSiteBusy.push_back(bId2);
                        //cout << "forming condensin association!" << endl;
                        dVec accel1 = vecNode[bId1].get_accel(), accel2 = vecNode[bId2].get_accel(), accel1_new, accel2_new;
                        double L0 = L0_TENSION_CONDENSIN, L_inv;
                        double k = K_TENSION_CONDENSIN_UNIT_MASS;
                        dVec unitVec12 = {0,0,0};
                        if (L > 0)
                        {
                            L_inv = 1./L;  // 1/nm
                            unitVec12 = {(coord2.x-coord1.x)*L_inv, (coord2.y-coord1.y)*L_inv, (coord2.z-coord1.z)*L_inv};
                        }

                        if ( abs(L - L0) > L_DIFF_RESOLUTION) // non-zero tension
                        {
                            double a_ten_con_size = k*(L-L0);
                            dVec a_ten_con = {a_ten_con_size*unitVec12.x, a_ten_con_size*unitVec12.y, a_ten_con_size*unitVec12.z};
                            accel1_new = {  accel1.x + a_ten_con.x, accel1.y + a_ten_con.y, accel1.z + a_ten_con.z };
                            accel2_new = {  accel2.x - a_ten_con.x, accel2.y - a_ten_con.y, accel2.z - a_ten_con.z };
                            // update acceleration
                            vecNode[bId1].set_accel(accel1_new);
                            vecNode[bId2].set_accel(accel2_new);
                            // update potential & total energy
                            double ePot12 = 0.5*k*(L-L0)*(L-L0);
                            energy["potential"] += ePot12;
                            energy["total"]     += ePot12;
                        }
                    }
                }

            }
        }
    }

    // (4) loop extrusion

}
void updateVelocBrownianSimulationOverdamp(vector<Chromosome> & vecChromosome, vector<Node> & vecNode, vector<int> & vecCondensinSite,
                                           vector<int> & vecBinderBlockSite, vector<Binder> & vecBinder,
                                           VoxMap & voxMap, unordered_map<string, double> & energy, vector<double> & randNormNum0, bool flagIsPrep, vector<iPair> & capturer, vector<int> & hit_valence, map<int, vector<int>> & Domain_CSB, bool flagLE_ON, bool flagIsInterphase)
{
    // (0) clear veloc values in Nodes
    for (vector<Node>::iterator it = vecNode.begin(); it != vecNode.end(); it++)
    {
        it->set_veloc({0, 0, 0});
    }

    // (1) entropic force & tension force (nearest bi-bead) & attraction force (nearest tri-bead)
    int k = 0;  // for using randNormNum0
    double m_inv = 1./MASS_BEAD;            // 1/kg
    double gamma_inv = 1./GAMMA_VISCOUS;    // s/kg
    for (vector<Chromosome>::iterator it = vecChromosome.begin(); it != vecChromosome.end(); it++)
    {
        vector<int> vecBeadIds = it->get_vecBeadIds();
        for (int i = 0; i < vecBeadIds.size(); i ++)
        {
            //double m_inv = 1./MASS_BEAD;    // (1/kg) make sure MASS_BEAD is calculated correctly

            int bId1 = vecBeadIds[i], bId2, bId3;   // consecutive ids of Bead
            dVec coord1 = vecNode[bId1].get_coord(), coord2, coord3;
            dVec veloc1 = vecNode[bId1].get_veloc(), veloc2, veloc3;
            dVec veloc1_new = veloc1, veloc2_new, veloc3_new;

            /* --------- [NO LONGER USED] (1.1) entropic force 1st bead ---------- */
            if (false)
            {
                //double theta    = acos(2.* ((double) rand() / (RAND_MAX)) - 1);
                //double phi      = 2.* PI * ((double) rand() / (RAND_MAX));

                double v_ent_x = randNormNum0[k++]*gamma_inv*1E-3;  // nm/s
                double v_ent_y = randNormNum0[k++]*gamma_inv*1E-3;
                double v_ent_z = randNormNum0[k++]*gamma_inv*1E-3;
                veloc1_new = {veloc1.x+v_ent_x, veloc1.y+v_ent_y, veloc1.z+v_ent_z};
            }
            
            
            
            
            
            
            /* --------------- ( ** 1.1/2 ) checking confiniment  ------------------------------------------------ */
            // exchange the bId2 etc!!!!!!!!!!!!!!!!!!!!!
/*
            if (true)
            {
				//int bId2 = ; // id of each single bead in a cylinder
				
				//----parameters for the cylindrical confinement
				double cos_theta1 = RAD_NUCLEOLUS / RAD_NUCLEUS * 0.5;
				double theta1     = acos(cos_theta1);
				double theta2	  = 2*theta1;
				double cos_theta2 = cos(theta2);
				double A_factor = 1.0; //scaling factor
				double vol_sphere = 4*PI/3. * RAD_NUCLEUS*RAD_NUCLEUS*RAD_NUCLEUS;
				double vol_cap    = PI/3. * RAD_NUCLEOLUS*RAD_NUCLEOLUS*RAD_NUCLEOLUS * (2 - 3*cos_theta1 + cos_theta1*cos_theta1*cos_theta1);
				double vol_cap_s  = PI/3. * RAD_NUCLEUS*RAD_NUCLEUS*RAD_NUCLEUS * (2 - 3*cos_theta2 + cos_theta2*cos_theta2*cos_theta2);
				double vol_total = vol_sphere - (vol_cap_s - vol_cap);
				double vol_occupy = vol_total * VOL_FRAC_OF_GENOME * A_factor;
				//double vol_occupy = vol_total * VOL_FRAC_OF_GENOME;
				double Ltot       = RAD_NUCLEOLUS;
				double rMax       = sqrt(vol_occupy / PI / Ltot); // this is a cylinder with vol fraction equal to genome size fraction
				
				//----parameters for pulling force
				double RAD_CYLINDER = RAD_NUCLEUS/3;//rMax;     //10.0; // ?? !! need an estimation 
				double RAD_CYLINDER_sq = RAD_CYLINDER*RAD_CYLINDER;
				double k = 1.0*K_TENSION;       // 1.0. *K_TENSION previous default?? !! need an estimation
				//dVec coord1 = vecNode[bId1].get_coord();
				dVec coord0 = {0, 0, coord1.z};  //  !!!!!!!!! just for simple case, when using initial configuration with binned cylinder !!!!!!!!!!!!!!!!! if tilted, then modify!
				//dVec ORIGIN = coord0;
				double L = dVecDist(coord0, coord1);
				//double L_sq = coord2.x*coord2.x + coord2.y*coord3.y + coord2.z*coord2.z;
				double L_inv;
				dVec veloc1 = vecNode[bId1].get_veloc(), veloc1_new;
				dVec unitVec12 = {0,0,0};

				if (L > 0)
				{
					L_inv = 1./L;  // 1/nm
					unitVec12 = {(coord1.x-coord0.x)*L_inv, (coord1.y-coord0.y)*L_inv, (coord1.z-coord0.z)*L_inv};
				}
				//if ( L_sq > RAD_CYLINDER_sq)
				if ( L > RAD_CYLINDER)  // or if ( L_sq > RAD_CYLINDER_sq)
				{
					
					double v_pull_size = k*(RAD_CYLINDER-L)*gamma_inv*1E-3;
					dVec veloc1_brw_temp = veloc1_new;
					dVec v_pull = {v_pull_size*unitVec12.x, v_pull_size*unitVec12.y, v_pull_size*unitVec12.z};
					veloc1_new = {  veloc1_brw_temp.x + v_pull.x, veloc1_brw_temp.y + v_pull.y, veloc1_brw_temp.z + v_pull.z};	
					vecNode[bId1].set_veloc(veloc1_new);
					
					//continue;
				}
				else
					
					dVec veloc1_new = veloc1; 
			}
			
	*/		
			
			
			

            /* --------- (1.2) tension force between 1st & 2nd bead
               &         (1.3) attraction force between 1st & 3rd bead ----------- */
            if (true)
            {
                if (i < vecBeadIds.size()-1)
                {
                    bId2 = vecBeadIds[i+1];
                    coord2 = vecNode[bId2].get_coord();
                    veloc2 = vecNode[bId2].get_veloc();
                    double L0 = L0_TENSION, L = dVecDist(coord1, coord2);   // nm
                    double k = K_TENSION;                                               // pN/nm     (N = kg*m/s/s; nN = kg*nm/s/s; pN = kg*pm/s/s)
                    dVec unitVec12 = {0,0,0};
                    if (L > 0)
                    {
                        double L_inv = 1./L;  // 1/nm
                        unitVec12 = {(coord2.x-coord1.x)*L_inv, (coord2.y-coord1.y)*L_inv, (coord2.z-coord1.z)*L_inv};
                    }

                    if ( abs(L - L0) > L_DIFF_RESOLUTION) // non-zero tension
                    {
                        double v_ten_size = k*(L-L0)*gamma_inv*1E-3;   // nm/s
                        dVec v_ten = {v_ten_size*unitVec12.x, v_ten_size*unitVec12.y, v_ten_size*unitVec12.z};
                        dVec veloc1_tmp1 = veloc1_new;
                        veloc1_new = {  veloc1_tmp1.x + v_ten.x, veloc1_tmp1.y + v_ten.y, veloc1_tmp1.z + v_ten.z };
                        veloc2_new = {  veloc2.x      - v_ten.x, veloc2.y      - v_ten.y, veloc2.z      - v_ten.z };
                        // update velocity
                        vecNode[bId2].set_veloc(veloc2_new);
                        // update potential & total energy
                        //double ePot12 = 0.5*k*(L-L0)*(L-L0);
                        //energy["potential"] += ePot12;
                        //energy["total"]     += ePot12;
                    }

                    if (i < vecBeadIds.size()-2)
                    {
                        bId3 = vecBeadIds[i+2];
                        coord3 = vecNode[bId3].get_coord();
                        veloc3 = vecNode[bId3].get_veloc();
                        double L = dVecDist(coord1, coord3), L_inv;
                        double k = K_ATTRACTION;
                        dVec unitVec13 = {0,0,0};
                        if (L > 0)
                        {
                            L_inv = 1./L;
                            unitVec13 = {(coord3.x-coord1.x)*L_inv, (coord3.y-coord1.y)*L_inv, (coord3.z-coord1.z)*L_inv};
                        }
                        double v_att_size = k*L*gamma_inv*1E-3;
                        dVec v_att = {v_att_size*unitVec13.x, v_att_size*unitVec13.y, v_att_size*unitVec13.z};
                        dVec veloc1_tmp2 = veloc1_new;
                        veloc1_new = {  veloc1_tmp2.x + v_att.x, veloc1_tmp2.y + v_att.y, veloc1_tmp2.z + v_att.z };
                        veloc3_new = {  veloc3.x - v_att.x, veloc3.y - v_att.y, veloc3.z - v_att.z };
                        // update acceleration
                        vecNode[bId3].set_veloc(veloc3_new);
                        // update potential & total energy
                        //double ePot13 = 0.5*k*L*L;
                        //energy["potential"] += ePot13;
                        //energy["total"]     += ePot13;
                    }
                }
                
                // update acceleration
                vecNode[bId1].set_veloc(veloc1_new);
            }
        }
    }

    // (2) repulsion force
    if (true)   // interaction between beads that collide (physically too close)
    {
        iVec dim = voxMap.get_dim();
        int nMax = dim.i*dim.j*dim.k;
        unordered_map<int, vector<int>> mapBeadIds = voxMap.get_mapBeadIds();
        unordered_map<int, vector<int>> mapBeadIdsHalo = voxMap.get_mapBeadIdsHalo();

        # pragma omp parallel for num_threads(NUM_THREADS) schedule(dynamic) \
        shared(mapBeadIds, mapBeadIdsHalo, vecNode)
        for (int boxId = 0; boxId < nMax; boxId ++)
        {
            vector<int> vecBeadIds0, vecBeadIds1;
            unordered_map<int, vector<int>>::const_iterator got = mapBeadIds.find(boxId);
            if ( got != mapBeadIds.end() )  // this box has Bead inside
            {
                for (vector<int>::iterator it = mapBeadIds[boxId].begin(); it != mapBeadIds[boxId].end(); it ++)
                    vecBeadIds0.push_back(*it);
            }
            unordered_map<int, vector<int>>::const_iterator got2 = mapBeadIdsHalo.find(boxId);
            if ( got2 != mapBeadIdsHalo.end() )  // this box has Bead inside
            {
                for (vector<int>::iterator it = mapBeadIdsHalo[boxId].begin(); it != mapBeadIdsHalo[boxId].end(); it ++)
                    vecBeadIds1.push_back(*it);
            }

            // STEP 1 : repulsion between Beads within the Box
            int nBeadInBox = vecBeadIds0.size();
            if (nBeadInBox > 1)
            {
                for (int l = 0; l < nBeadInBox-1; l++)
                {
                    for (int m = l+1; m < nBeadInBox; m++)
                    {
                        int bId1 = vecBeadIds0[l], bId2 = vecBeadIds0[m];
                        dVec coord1 = vecNode[bId1].get_coord(), coord2 = vecNode[bId2].get_coord();
                        double LTHR1 = LTHR1_REPULSION, LTHR2 = LTHR2_REPULSION, L = dVecDist(coord1, coord2);   // nm
                        if (abs(bId1-bId2) >= 2 &&  L < LTHR2) // repulsion only when overlapping of Beads occurs && when two Beads are not connected
                        {
                            double L_inv;
                            dVec unitVec12 = {0,0,0};
                            if (L > 0)
                            {
                                L_inv = 1./L;   // 1/nm
                                unitVec12 = {(coord2.x-coord1.x)*L_inv, (coord2.y-coord1.y)*L_inv, (coord2.z-coord1.z)*L_inv};
                            }
                            double f = F_REPULSION;
                            dVec veloc1 = vecNode[bId1].get_veloc(), veloc2 = vecNode[bId2].get_veloc(), veloc1_new, veloc2_new;

                            //double a_rep_size = f*m_inv*pow(L_inv, 12);
                            double v_rep_size = f*gamma_inv*1E-3;
/*
                            if (L < LTHR1)
                                a_rep_size = f;
                            if (L >= LTHR1)
                                a_rep_size = f*pow(L_inv, 12);
*/
                            dVec v_rep = {v_rep_size*unitVec12.x, v_rep_size*unitVec12.y, v_rep_size*unitVec12.z};

                            veloc1_new = {  veloc1.x - v_rep.x, veloc1.y - v_rep.y, veloc1.z - v_rep.z };
                            vecNode[bId1].set_veloc(veloc1_new);

                            veloc2_new = {  veloc2.x + v_rep.x, veloc2.y + v_rep.y, veloc2.z + v_rep.z };
                            vecNode[bId2].set_veloc(veloc2_new);

                            //double ePot12 = f*(LTHR2-L);
/*
                            if (L < LTHR1)
                                ePot12 = -f*L;
                            if (L >= LTHR1)
                                ePot12 = -f*pow(L_inv, 11)/11.;
*/
                            //energy["potential"] += ePot12;
                            //energy["total"]     += ePot12;
                        }
                    }
                }
            }
            // STEP 2 : repulsion between one Bead within the Box and another Bead within Halo
            int nBeadInHalo = vecBeadIds1.size();
            if (nBeadInBox > 0 && nBeadInHalo > 0)
            {
                for (int l = 0; l < nBeadInBox-1; l++)
                {
                    for (int m = 0; m < nBeadInHalo; m++)
                    {
                        int bId1 = vecBeadIds0[l], bId2 = vecBeadIds1[m];
                        if (bId1 < bId2 && abs(bId1-bId2) >= 2) // bId1 < bId2 : avoid calculating twice! && when two Beads are not connected
                        {
                            dVec coord1 = vecNode[bId1].get_coord(), coord2 = vecNode[bId2].get_coord();
                            double LTHR1 = LTHR1_REPULSION, LTHR2 = LTHR2_REPULSION, L = dVecDist(coord1, coord2);   // nm
                            if (L < LTHR2)
                            {
                                double L_inv;
                                dVec unitVec12 = {0,0,0};
                                if (L > 0)
                                {
                                    L_inv = 1./L;   // 1/nm
                                    unitVec12 = {(coord2.x-coord1.x)*L_inv, (coord2.y-coord1.y)*L_inv, (coord2.z-coord1.z)*L_inv};
                                }
                                double f = F_REPULSION;
                                dVec veloc1 = vecNode[bId1].get_veloc(), veloc2 = vecNode[bId2].get_veloc(), veloc1_new, veloc2_new;
                                //double a_rep_size = f*m_inv*pow(L_inv, 12);  //  nm /s/s
                                double v_rep_size = f*gamma_inv*1E-3;
/*
                                if (L < LTHR1)
                                    a_rep_size = f;
                                    //a_rep_size = f*m_inv;
                                if (L >= LTHR1)
                                    a_rep_size = f*pow(L_inv, 12);
*/
                                dVec v_rep = {v_rep_size*unitVec12.x, v_rep_size*unitVec12.y, v_rep_size*unitVec12.z};

                                veloc1_new = {  veloc1.x - v_rep.x, veloc1.y - v_rep.y, veloc1.z - v_rep.z };
                                vecNode[bId1].set_veloc(veloc1_new);

                                veloc2_new = {  veloc2.x + v_rep.x, veloc2.y + v_rep.y, veloc2.z + v_rep.z };
                                vecNode[bId2].set_veloc(veloc2_new);

                                //double ePot12 = f*(LTHR2-L);
/*
                                if (L < LTHR1)
                                    ePot12 = -f*L;
                                if (L >= LTHR1)
                                    ePot12 = -f*pow(L_inv, 11)/11.;
*/
                                //energy["potential"] += ePot12;
                                //energy["total"]     += ePot12;
                            }
                        }
                    }
                }
            }
        }
    }

//////
	// (3) condensin attraction (diffusion capture of condensin binding sites; see below for alternative method: diffusion capture of binders)
	//if (DIFFUSION_CAPTURE_ON == true && flagIsPrep == false && flagDC_Interphase == true)	
	if (DIFFUSION_CAPTURE_ON == true && flagIsPrep == false && flagIsInterphase == true && SWITCH_INTERPHASE_MITOSIS > 0)
	{
		//cout << " **                    NEW iteration                 **" << endl;
		/*
		vector<int> vecTADInterphaseBoundary;
		vector<int> vecTADMitosisBoundary;
		vecTADInterphaseBoundary.push_back(0);
		vecTADMitosisBoundary.push_back(0);
		
		
		int cnt = 0;
		for (int i = 1; i <= NUM_TAD_INTERPHASE; i++)
		{
			cnt = cnt + INTERPHASE_TAD_SIZE;
			vecTADInterphaseBoundary.push_back(cnt);
		}
		
		for (int i = 1; i <= NUM_TAD_MITOSIS; i++)
		{
			cnt = cnt + MITOSIS_TAD_SIZE;
			vecTADMitosisBoundary.push_back(cnt);
		}
		
		//cout << vecTADMitosisBoundary.size() << endl;
		//cout << vecTADInterphaseBoundary.size() << endl;
		
		*/
		
		
		vector<int> vecCondensinSiteBusy; // contain bead ids already interacting with another
		capturer.clear();
		hit_valence.clear();
		vector<int> vecDC_InDomain;
		
		for ( auto map_iter = Domain_CSB.begin(); map_iter != Domain_CSB.end(); ++map_iter)
		{
			vecDC_InDomain.clear();
			//cout << " Diffusion capture for domain key: " << map_iter->first << "    and CSB values: [";
			//vector<int> vecDC_InDomain;
			for (auto vec_iter = map_iter->second.begin() ; vec_iter != map_iter->second.end(); ++vec_iter)
			{
				//cout << *vec_iter << " ";
				int vecbId = *vec_iter;
				//cout << " vecbId "<< vecbId << "-- ";
				vecDC_InDomain.push_back(vecbId);
				
			}
			//cout << "]\n" ;
			
			//cout << "CSB values for actual domain SIZE " << vecDC_InDomain.size() << " // " << endl; 
			//vecDC_InDomain.clear();
				/*for (int k = 0; vecDC_InDomain.size(); k++)
				{
					cout << vecDC_InDomain[k] << " ";
				}
				cout << " ]" << endl ;
				*/
		
			for (int i = 0; i < vecDC_InDomain.size(); i++)
			{
				for (int j = i+1; j < vecDC_InDomain.size(); j++) 
				{
					if ((rand() / double(RAND_MAX)) > P_DISSOCIATE_CONDENSIN) 
					{
						int bId1 = vecDC_InDomain[i], bId2 = vecDC_InDomain[j];
						dVec coord1 = vecNode[bId1].get_coord(), coord2 = vecNode[bId2].get_coord();
						double L = dVecDist(coord1, coord2);

						if (DIFFUSION_CAPTURE == "probabilistic")
						{
							if ((rand() / double(RAND_MAX)) > P_DC_PROBAB)
							{ 
					// multivalent diffusion capture
								if (L < LTHR_TENSION_CONDENSIN)
								{
									vecCondensinSiteBusy.push_back(bId1);
									vecCondensinSiteBusy.push_back(bId2);

									//cout << "forming condensin association!" << endl;
									dVec veloc1 = vecNode[bId1].get_veloc(), veloc2 = vecNode[bId2].get_veloc(), veloc1_new, veloc2_new;
									double L0 = L0_TENSION_CONDENSIN, L_inv;
									double k = K_TENSION_CONDENSIN;
									dVec unitVec12 = { 0, 0, 0 };
									if (L > 0)
									{
										L_inv = 1. / L;  // 1/nm
										unitVec12 = {(coord2.x-coord1.x)*L_inv, (coord2.y-coord1.y)*L_inv, (coord2.z-coord1.z)*L_inv};
									}

									if (abs(L - L0) > L_DIFF_RESOLUTION) // non-zero tension
									{
										double v_ten_con_size = k * (L - L0)* gamma_inv * 1E-3;
										dVec v_ten_con = { v_ten_con_size* unitVec12.x, v_ten_con_size* unitVec12.y, v_ten_con_size* unitVec12.z };
										veloc1_new = {veloc1.x + v_ten_con.x, veloc1.y + v_ten_con.y, veloc1.z + v_ten_con.z};
										veloc2_new = {veloc2.x - v_ten_con.x, veloc2.y - v_ten_con.y, veloc2.z - v_ten_con.z};
										// update acceleration
										vecNode[bId1].set_veloc(veloc1_new);
										vecNode[bId2].set_veloc(veloc2_new);
										// update potential & total energy
										//double ePot12 = 0.5*k*(L-L0)*(L-L0);
										//energy["potential"] += ePot12;
										//energy["total"]     += ePot12;
									}
									capturer.push_back( {bId1, bId2});
								}
								
							}
														
					// bivalent diffusion capture
							else 
							{
								if (L < LTHR_TENSION_CONDENSIN && find(vecCondensinSiteBusy.begin(), vecCondensinSiteBusy.end(), bId1) == vecCondensinSiteBusy.end() && find(vecCondensinSiteBusy.begin(), vecCondensinSiteBusy.end(), bId2) == vecCondensinSiteBusy.end())
								//&& find(vecCondensinSiteBusy.begin(),vecCondensinSiteBusy.end(), bId1) == vecCondensinSiteBusy.end()
								//&& find(vecCondensinSiteBusy.begin(),vecCondensinSiteBusy.end(), bId2) == vecCondensinSiteBusy.end()) // dissociation & avoid more than one contact
								{
									vecCondensinSiteBusy.push_back(bId1);
									vecCondensinSiteBusy.push_back(bId2);

									//cout << "forming condensin association!" << endl;
									dVec veloc1 = vecNode[bId1].get_veloc(), veloc2 = vecNode[bId2].get_veloc(),veloc1_new, veloc2_new;
									double L0 = L0_TENSION_CONDENSIN, L_inv;
									double k = K_TENSION_CONDENSIN;
									dVec unitVec12 = { 0, 0, 0 };
									if (L > 0)
									{
										L_inv = 1. / L;  // 1/nm
										unitVec12 = {(coord2.x-coord1.x)*L_inv, (coord2.y-coord1.y)*L_inv, (coord2.z-coord1.z)*L_inv};
									}

									if ( abs(L - L0) > L_DIFF_RESOLUTION) // non-zero tension
									{
										double v_ten_con_size = k*(L-L0)*gamma_inv*1E-3;
										dVec v_ten_con = {v_ten_con_size*unitVec12.x, v_ten_con_size*unitVec12.y, v_ten_con_size*unitVec12.z};
										veloc1_new = {veloc1.x + v_ten_con.x, veloc1.y + v_ten_con.y, veloc1.z + v_ten_con.z};
										veloc2_new = {veloc2.x - v_ten_con.x, veloc2.y - v_ten_con.y, veloc2.z - v_ten_con.z};
										// update acceleration
										vecNode[bId1].set_veloc(veloc1_new);
										vecNode[bId2].set_veloc(veloc2_new);
										// update potential & total energy
										//double ePot12 = 0.5*k*(L-L0)*(L-L0);
										//energy["potential"] += ePot12;
										//energy["total"]     += ePot12;
									}
								}	capturer.push_back( {bId1, bId2});
								
							}
							
						}

	//					if (DIFFUSION_CAPTURE_PROBAB == true)


						if (DIFFUSION_CAPTURE == "bivalent")
						{
							
							if (   L < LTHR_TENSION_CONDENSIN
								&& find(vecCondensinSiteBusy.begin(), vecCondensinSiteBusy.end(), bId1) == vecCondensinSiteBusy.end()
								&& find(vecCondensinSiteBusy.begin(), vecCondensinSiteBusy.end(), bId2) == vecCondensinSiteBusy.end())  // dissociation & avoid more than one contact
								//if ( L < LTHR_TENSION_CONDENSIN )
							{
									vecCondensinSiteBusy.push_back(bId1);
									vecCondensinSiteBusy.push_back(bId2);
								
															//cout << "forming condensin association!" << endl;
									dVec veloc1 = vecNode[bId1].get_veloc(), veloc2 = vecNode[bId2].get_veloc(), veloc1_new, veloc2_new;
									double L0 = L0_TENSION_CONDENSIN, L_inv;
									double k = K_TENSION_CONDENSIN;
									dVec unitVec12 = {0,0,0};
									if (L > 0)
									{
										L_inv = 1./L;  // 1/nm
										unitVec12 = {(coord2.x-coord1.x)*L_inv, (coord2.y-coord1.y)*L_inv, (coord2.z-coord1.z)*L_inv};
									}

									if ( abs(L - L0) > L_DIFF_RESOLUTION) // non-zero tension
									{
										double v_ten_con_size = k*(L-L0)*gamma_inv*1E-3;
										dVec v_ten_con = {v_ten_con_size*unitVec12.x, v_ten_con_size*unitVec12.y, v_ten_con_size*unitVec12.z};
										veloc1_new = {veloc1.x + v_ten_con.x, veloc1.y + v_ten_con.y, veloc1.z + v_ten_con.z};
										veloc2_new = {veloc2.x - v_ten_con.x, veloc2.y - v_ten_con.y, veloc2.z - v_ten_con.z};
																// update acceleration
										vecNode[bId1].set_veloc(veloc1_new);
										vecNode[bId2].set_veloc(veloc2_new);
																// update potential & total energy
																//double ePot12 = 0.5*k*(L-L0)*(L-L0);
																//energy["potential"] += ePot12;
																//energy["total"]     += ePot12;
									}
									capturer.push_back( {bId1, bId2});
								//}
								
							}
						}
							
						
						if (DIFFUSION_CAPTURE == "multivalent")
						{	
							//if (   L < LTHR_TENSION_CONDENSIN
								//&& find(vecCondensinSiteBusy.begin(), vecCondensinSiteBusy.end(), bId1) == vecCondensinSiteBusy.end()
								//&& find(vecCondensinSiteBusy.begin(), vecCondensinSiteBusy.end(), bId2) == vecCondensinSiteBusy.end())  // dissociation & avoid more than one contact
							if ( L < LTHR_TENSION_CONDENSIN )							
							{
								/*	
								int bId1_TADid;
								int bId2atta_TADid;
								
								for (int i = 0; i < vecTADInterphaseBoundary.size(); i++)
								{
									if ( bId1 >= vecTADInterphaseBoundary[i] && bId1 < vecTADInterphaseBoundary[i+1])
									{
										bId1_TADid = i;
									}
									if ( bId2 >= vecTADInterphaseBoundary[i] && bId2 < vecTADInterphaseBoundary[i+1])
									{
										bId2atta_TADid = i;
									}								
								}
								//cout << bId1_TADid << "  beadId1_TADid" << endl;
								//cout << bId2atta_TADid << "  beadId2_TADid" << endl;
								//cout << "blabla" << endl;
								
								if ( bId1_TADid == bId2atta_TADid)							
								{
									*/
									//cout << "working" << endl;
									vecCondensinSiteBusy.push_back(bId1);
									vecCondensinSiteBusy.push_back(bId2);
								
															//cout << "forming condensin association!" << endl;
									dVec veloc1 = vecNode[bId1].get_veloc(), veloc2 = vecNode[bId2].get_veloc(), veloc1_new, veloc2_new;
									double L0 = L0_TENSION_CONDENSIN, L_inv;
									double k = K_TENSION_CONDENSIN;
									dVec unitVec12 = {0,0,0};
									if (L > 0)
									{
										L_inv = 1./L;  // 1/nm
										unitVec12 = {(coord2.x-coord1.x)*L_inv, (coord2.y-coord1.y)*L_inv, (coord2.z-coord1.z)*L_inv};
									}

									if ( abs(L - L0) > L_DIFF_RESOLUTION) // non-zero tension
									{
										double v_ten_con_size = k*(L-L0)*gamma_inv*1E-3;
										dVec v_ten_con = {v_ten_con_size*unitVec12.x, v_ten_con_size*unitVec12.y, v_ten_con_size*unitVec12.z};
										veloc1_new = {veloc1.x + v_ten_con.x, veloc1.y + v_ten_con.y, veloc1.z + v_ten_con.z};
										veloc2_new = {veloc2.x - v_ten_con.x, veloc2.y - v_ten_con.y, veloc2.z - v_ten_con.z};
																// update acceleration
										vecNode[bId1].set_veloc(veloc1_new);
										vecNode[bId2].set_veloc(veloc2_new);
																// update potential & total energy
																//double ePot12 = 0.5*k*(L-L0)*(L-L0);
																//energy["potential"] += ePot12;
																//energy["total"]     += ePot12;
									}
									capturer.push_back( {bId1, bId2});
								//}
								
							}
						}						
		
					}
				}
			}
		
		}
	}

	if (DIFFUSION_CAPTURE_ON == true && flagIsPrep == false && flagIsInterphase == false && SWITCH_INTERPHASE_MITOSIS > 0)
	{
		//cout << " **                    NEW iteration                 **" << endl;
		/*
		vector<int> vecTADInterphaseBoundary;
		vector<int> vecTADMitosisBoundary;
		vecTADInterphaseBoundary.push_back(0);
		vecTADMitosisBoundary.push_back(0);
		
		
		int cnt = 0;
		for (int i = 1; i <= NUM_TAD_INTERPHASE; i++)
		{
			cnt = cnt + INTERPHASE_TAD_SIZE;
			vecTADInterphaseBoundary.push_back(cnt);
		}
		
		for (int i = 1; i <= NUM_TAD_MITOSIS; i++)
		{
			cnt = cnt + MITOSIS_TAD_SIZE;
			vecTADMitosisBoundary.push_back(cnt);
		}
		
		//cout << vecTADMitosisBoundary.size() << endl;
		//cout << vecTADInterphaseBoundary.size() << endl;
		
		*/
		
		
		vector<int> vecCondensinSiteBusy; // contain bead ids already interacting with another
		capturer.clear();
		hit_valence.clear();
		vector<int> vecDC_InDomain;
		
		for ( auto map_iter = Domain_CSB.begin(); map_iter != Domain_CSB.end(); ++map_iter)
		{
			vecDC_InDomain.clear();
			//cout << " Diffusion capture for domain key: " << map_iter->first << "    and CSB values: [";
			//vector<int> vecDC_InDomain;
			for (auto vec_iter = map_iter->second.begin() ; vec_iter != map_iter->second.end(); ++vec_iter)
			{
				//cout << *vec_iter << " ";
				int vecbId = *vec_iter;
				//cout << " vecbId "<< vecbId << "-- ";
				vecDC_InDomain.push_back(vecbId);
				
			}
			//cout << "]\n" ;
			
			//cout << "CSB values for actual domain SIZE " << vecDC_InDomain.size() << " // " << endl; 
			//vecDC_InDomain.clear();
				/*for (int k = 0; vecDC_InDomain.size(); k++)
				{
					cout << vecDC_InDomain[k] << " ";
				}
				cout << " ]" << endl ;
				*/
		
			for (int i = 0; i < vecDC_InDomain.size(); i++)
			{
				for (int j = i+1; j < vecDC_InDomain.size(); j++) 
				{
					if ((rand() / double(RAND_MAX)) > P_DISSOCIATE_MITOSIS) 
					{
						int bId1 = vecDC_InDomain[i], bId2 = vecDC_InDomain[j];
						dVec coord1 = vecNode[bId1].get_coord(), coord2 = vecNode[bId2].get_coord();
						double L = dVecDist(coord1, coord2);

						if (DIFFUSION_CAPTURE_M == "probabilistic")
						{
							if ((rand() / double(RAND_MAX)) > P_DC_PROBAB)
							{ 
					// multivalent diffusion capture
								if (L < LTHR_TENSION_CONDENSIN)
								{
									vecCondensinSiteBusy.push_back(bId1);
									vecCondensinSiteBusy.push_back(bId2);

									//cout << "forming condensin association!" << endl;
									dVec veloc1 = vecNode[bId1].get_veloc(), veloc2 = vecNode[bId2].get_veloc(), veloc1_new, veloc2_new;
									double L0 = L0_TENSION_CONDENSIN, L_inv;
									double k = K_TENSION_CONDENSIN;
									dVec unitVec12 = { 0, 0, 0 };
									if (L > 0)
									{
										L_inv = 1. / L;  // 1/nm
										unitVec12 = {(coord2.x-coord1.x)*L_inv, (coord2.y-coord1.y)*L_inv, (coord2.z-coord1.z)*L_inv};
									}

									if (abs(L - L0) > L_DIFF_RESOLUTION) // non-zero tension
									{
										double v_ten_con_size = k * (L - L0)* gamma_inv * 1E-3;
										dVec v_ten_con = { v_ten_con_size* unitVec12.x, v_ten_con_size* unitVec12.y, v_ten_con_size* unitVec12.z };
										veloc1_new = {veloc1.x + v_ten_con.x, veloc1.y + v_ten_con.y, veloc1.z + v_ten_con.z};
										veloc2_new = {veloc2.x - v_ten_con.x, veloc2.y - v_ten_con.y, veloc2.z - v_ten_con.z};
										// update acceleration
										vecNode[bId1].set_veloc(veloc1_new);
										vecNode[bId2].set_veloc(veloc2_new);
										// update potential & total energy
										//double ePot12 = 0.5*k*(L-L0)*(L-L0);
										//energy["potential"] += ePot12;
										//energy["total"]     += ePot12;
									}
									capturer.push_back( {bId1, bId2});
								}
								
							}
														
					// bivalent diffusion capture
							else 
							{
								if (L < LTHR_TENSION_CONDENSIN && find(vecCondensinSiteBusy.begin(), vecCondensinSiteBusy.end(), bId1) == vecCondensinSiteBusy.end() && find(vecCondensinSiteBusy.begin(), vecCondensinSiteBusy.end(), bId2) == vecCondensinSiteBusy.end())
								//&& find(vecCondensinSiteBusy.begin(),vecCondensinSiteBusy.end(), bId1) == vecCondensinSiteBusy.end()
								//&& find(vecCondensinSiteBusy.begin(),vecCondensinSiteBusy.end(), bId2) == vecCondensinSiteBusy.end()) // dissociation & avoid more than one contact
								{
									vecCondensinSiteBusy.push_back(bId1);
									vecCondensinSiteBusy.push_back(bId2);

									//cout << "forming condensin association!" << endl;
									dVec veloc1 = vecNode[bId1].get_veloc(), veloc2 = vecNode[bId2].get_veloc(),veloc1_new, veloc2_new;
									double L0 = L0_TENSION_CONDENSIN, L_inv;
									double k = K_TENSION_CONDENSIN;
									dVec unitVec12 = { 0, 0, 0 };
									if (L > 0)
									{
										L_inv = 1. / L;  // 1/nm
										unitVec12 = {(coord2.x-coord1.x)*L_inv, (coord2.y-coord1.y)*L_inv, (coord2.z-coord1.z)*L_inv};
									}

									if ( abs(L - L0) > L_DIFF_RESOLUTION) // non-zero tension
									{
										double v_ten_con_size = k*(L-L0)*gamma_inv*1E-3;
										dVec v_ten_con = {v_ten_con_size*unitVec12.x, v_ten_con_size*unitVec12.y, v_ten_con_size*unitVec12.z};
										veloc1_new = {veloc1.x + v_ten_con.x, veloc1.y + v_ten_con.y, veloc1.z + v_ten_con.z};
										veloc2_new = {veloc2.x - v_ten_con.x, veloc2.y - v_ten_con.y, veloc2.z - v_ten_con.z};
										// update acceleration
										vecNode[bId1].set_veloc(veloc1_new);
										vecNode[bId2].set_veloc(veloc2_new);
										// update potential & total energy
										//double ePot12 = 0.5*k*(L-L0)*(L-L0);
										//energy["potential"] += ePot12;
										//energy["total"]     += ePot12;
									}
								}	capturer.push_back( {bId1, bId2});
								
							}
							
						}

	//					if (DIFFUSION_CAPTURE_PROBAB == true)


						if (DIFFUSION_CAPTURE_M == "bivalent")
						{
							
							if (   L < LTHR_TENSION_CONDENSIN
								&& find(vecCondensinSiteBusy.begin(), vecCondensinSiteBusy.end(), bId1) == vecCondensinSiteBusy.end()
								&& find(vecCondensinSiteBusy.begin(), vecCondensinSiteBusy.end(), bId2) == vecCondensinSiteBusy.end())  // dissociation & avoid more than one contact
								//if ( L < LTHR_TENSION_CONDENSIN )
							{
									vecCondensinSiteBusy.push_back(bId1);
									vecCondensinSiteBusy.push_back(bId2);
								
															//cout << "forming condensin association!" << endl;
									dVec veloc1 = vecNode[bId1].get_veloc(), veloc2 = vecNode[bId2].get_veloc(), veloc1_new, veloc2_new;
									double L0 = L0_TENSION_CONDENSIN, L_inv;
									double k = K_TENSION_CONDENSIN;
									dVec unitVec12 = {0,0,0};
									if (L > 0)
									{
										L_inv = 1./L;  // 1/nm
										unitVec12 = {(coord2.x-coord1.x)*L_inv, (coord2.y-coord1.y)*L_inv, (coord2.z-coord1.z)*L_inv};
									}

									if ( abs(L - L0) > L_DIFF_RESOLUTION) // non-zero tension
									{
										double v_ten_con_size = k*(L-L0)*gamma_inv*1E-3;
										dVec v_ten_con = {v_ten_con_size*unitVec12.x, v_ten_con_size*unitVec12.y, v_ten_con_size*unitVec12.z};
										veloc1_new = {veloc1.x + v_ten_con.x, veloc1.y + v_ten_con.y, veloc1.z + v_ten_con.z};
										veloc2_new = {veloc2.x - v_ten_con.x, veloc2.y - v_ten_con.y, veloc2.z - v_ten_con.z};
																// update acceleration
										vecNode[bId1].set_veloc(veloc1_new);
										vecNode[bId2].set_veloc(veloc2_new);
																// update potential & total energy
																//double ePot12 = 0.5*k*(L-L0)*(L-L0);
																//energy["potential"] += ePot12;
																//energy["total"]     += ePot12;
									}
									capturer.push_back( {bId1, bId2});
								//}
								
							}
						}
							
						
						if (DIFFUSION_CAPTURE_M == "multivalent")
						{	
							//if (   L < LTHR_TENSION_CONDENSIN
								//&& find(vecCondensinSiteBusy.begin(), vecCondensinSiteBusy.end(), bId1) == vecCondensinSiteBusy.end()
								//&& find(vecCondensinSiteBusy.begin(), vecCondensinSiteBusy.end(), bId2) == vecCondensinSiteBusy.end())  // dissociation & avoid more than one contact
							if ( L < LTHR_TENSION_CONDENSIN )							
							{
								/*	
								int bId1_TADid;
								int bId2atta_TADid;
								
								for (int i = 0; i < vecTADInterphaseBoundary.size(); i++)
								{
									if ( bId1 >= vecTADInterphaseBoundary[i] && bId1 < vecTADInterphaseBoundary[i+1])
									{
										bId1_TADid = i;
									}
									if ( bId2 >= vecTADInterphaseBoundary[i] && bId2 < vecTADInterphaseBoundary[i+1])
									{
										bId2atta_TADid = i;
									}								
								}
								//cout << bId1_TADid << "  beadId1_TADid" << endl;
								//cout << bId2atta_TADid << "  beadId2_TADid" << endl;
								//cout << "blabla" << endl;
								
								if ( bId1_TADid == bId2atta_TADid)							
								{
									*/
									//cout << "working" << endl;
									vecCondensinSiteBusy.push_back(bId1);
									vecCondensinSiteBusy.push_back(bId2);
								
															//cout << "forming condensin association!" << endl;
									dVec veloc1 = vecNode[bId1].get_veloc(), veloc2 = vecNode[bId2].get_veloc(), veloc1_new, veloc2_new;
									double L0 = L0_TENSION_CONDENSIN, L_inv;
									double k = K_TENSION_CONDENSIN;
									dVec unitVec12 = {0,0,0};
									if (L > 0)
									{
										L_inv = 1./L;  // 1/nm
										unitVec12 = {(coord2.x-coord1.x)*L_inv, (coord2.y-coord1.y)*L_inv, (coord2.z-coord1.z)*L_inv};
									}

									if ( abs(L - L0) > L_DIFF_RESOLUTION) // non-zero tension
									{
										double v_ten_con_size = k*(L-L0)*gamma_inv*1E-3;
										dVec v_ten_con = {v_ten_con_size*unitVec12.x, v_ten_con_size*unitVec12.y, v_ten_con_size*unitVec12.z};
										veloc1_new = {veloc1.x + v_ten_con.x, veloc1.y + v_ten_con.y, veloc1.z + v_ten_con.z};
										veloc2_new = {veloc2.x - v_ten_con.x, veloc2.y - v_ten_con.y, veloc2.z - v_ten_con.z};
																// update acceleration
										vecNode[bId1].set_veloc(veloc1_new);
										vecNode[bId2].set_veloc(veloc2_new);
																// update potential & total energy
																//double ePot12 = 0.5*k*(L-L0)*(L-L0);
																//energy["potential"] += ePot12;
																//energy["total"]     += ePot12;
									}
									capturer.push_back( {bId1, bId2});
								//}
								
							}
						}						
		
					}
				}
			}
		
		}
	}
		
//
	if (DIFFUSION_CAPTURE_ON == true && flagIsPrep == false && SWITCH_INTERPHASE_MITOSIS == 0 && DC_MECHANISM == "BF-BF")
	{
		//cout << "BF-BF"<< endl;
		vector<int> vecCondensinSiteBusy; // contain bead ids already interacting with another
		capturer.clear();
		hit_valence.clear();


		vector<int> vecAllFBinderBeadIds;
		for (vector<Binder>::iterator it = vecBinder.begin(); it != vecBinder.end(); it++) 
		{
			vector<int> vecBeadIds = it->get_vecBeadIds();
			int bIdF = vecBeadIds[0];
			vecAllFBinderBeadIds.push_back(vecBeadIds[0]);
		}

		for (int i = 0; i < vecAllFBinderBeadIds.size(); i++)
		{		
			for (int j = i+1; j < vecAllFBinderBeadIds.size(); j++) 
			{
				if ((rand() / double(RAND_MAX)) > P_DISSOCIATE_CONDENSIN) 
				{
					int bId1 = vecAllFBinderBeadIds[i], bId2 = vecAllFBinderBeadIds[j];
					dVec coord1 = vecNode[bId1].get_coord(), coord2 = vecNode[bId2].get_coord();
					double L = dVecDist(coord1, coord2);

					if (DIFFUSION_CAPTURE == "probabilistic")
					{
						if ((rand() / double(RAND_MAX)) > P_DC_PROBAB)
						{ 
				// multivalent diffusion capture
							if (L < LTHR_TENSION_CONDENSIN)
							{
								vecCondensinSiteBusy.push_back(bId1);
								vecCondensinSiteBusy.push_back(bId2);

								//cout << "forming condensin association!" << endl;
								dVec veloc1 = vecNode[bId1].get_veloc(), veloc2 = vecNode[bId2].get_veloc(), veloc1_new, veloc2_new;
								double L0 = L0_TENSION_CONDENSIN, L_inv;
								double k = K_TENSION_CONDENSIN;
								dVec unitVec12 = { 0, 0, 0 };
								if (L > 0)
								{
									L_inv = 1. / L;  // 1/nm
									unitVec12 = {(coord2.x-coord1.x)*L_inv, (coord2.y-coord1.y)*L_inv, (coord2.z-coord1.z)*L_inv};
								}

								if (abs(L - L0) > L_DIFF_RESOLUTION) // non-zero tension
								{
									double v_ten_con_size = k * (L - L0)* gamma_inv * 1E-3;
									dVec v_ten_con = { v_ten_con_size* unitVec12.x, v_ten_con_size* unitVec12.y, v_ten_con_size* unitVec12.z };
									veloc1_new = {veloc1.x + v_ten_con.x, veloc1.y + v_ten_con.y, veloc1.z + v_ten_con.z};
									veloc2_new = {veloc2.x - v_ten_con.x, veloc2.y - v_ten_con.y, veloc2.z - v_ten_con.z};
									// update acceleration
									vecNode[bId1].set_veloc(veloc1_new);
									vecNode[bId2].set_veloc(veloc2_new);
									// update potential & total energy
									//double ePot12 = 0.5*k*(L-L0)*(L-L0);
									//energy["potential"] += ePot12;
									//energy["total"]     += ePot12;
								}
								capturer.push_back( {bId1, bId2});
							}
							
						}
  													
			    // bivalent diffusion capture
						else 
						{
							if (L < LTHR_TENSION_CONDENSIN && find(vecCondensinSiteBusy.begin(), vecCondensinSiteBusy.end(), bId1) == vecCondensinSiteBusy.end() && find(vecCondensinSiteBusy.begin(), vecCondensinSiteBusy.end(), bId2) == vecCondensinSiteBusy.end())
							//&& find(vecCondensinSiteBusy.begin(),vecCondensinSiteBusy.end(), bId1) == vecCondensinSiteBusy.end()
							//&& find(vecCondensinSiteBusy.begin(),vecCondensinSiteBusy.end(), bId2) == vecCondensinSiteBusy.end()) // dissociation & avoid more than one contact
							{
								vecCondensinSiteBusy.push_back(bId1);
								vecCondensinSiteBusy.push_back(bId2);

								//cout << "forming condensin association!" << endl;
								dVec veloc1 = vecNode[bId1].get_veloc(), veloc2 = vecNode[bId2].get_veloc(),veloc1_new, veloc2_new;
								double L0 = L0_TENSION_CONDENSIN, L_inv;
								double k = K_TENSION_CONDENSIN;
								dVec unitVec12 = { 0, 0, 0 };
								if (L > 0)
								{
									L_inv = 1. / L;  // 1/nm
									unitVec12 = {(coord2.x-coord1.x)*L_inv, (coord2.y-coord1.y)*L_inv, (coord2.z-coord1.z)*L_inv};
								}

								if ( abs(L - L0) > L_DIFF_RESOLUTION) // non-zero tension
								{
									double v_ten_con_size = k*(L-L0)*gamma_inv*1E-3;
									dVec v_ten_con = {v_ten_con_size*unitVec12.x, v_ten_con_size*unitVec12.y, v_ten_con_size*unitVec12.z};
									veloc1_new = {veloc1.x + v_ten_con.x, veloc1.y + v_ten_con.y, veloc1.z + v_ten_con.z};
									veloc2_new = {veloc2.x - v_ten_con.x, veloc2.y - v_ten_con.y, veloc2.z - v_ten_con.z};
									// update acceleration
									vecNode[bId1].set_veloc(veloc1_new);
									vecNode[bId2].set_veloc(veloc2_new);
									// update potential & total energy
									//double ePot12 = 0.5*k*(L-L0)*(L-L0);
									//energy["potential"] += ePot12;
									//energy["total"]     += ePot12;
								}
							}	capturer.push_back( {bId1, bId2});
							
						}
						
					}

//					if (DIFFUSION_CAPTURE_PROBAB == true)

					if (DIFFUSION_CAPTURE == "bivalent")
					{
						
						if (   L < LTHR_TENSION_CONDENSIN
							&& find(vecCondensinSiteBusy.begin(), vecCondensinSiteBusy.end(), bId1) == vecCondensinSiteBusy.end()
							&& find(vecCondensinSiteBusy.begin(), vecCondensinSiteBusy.end(), bId2) == vecCondensinSiteBusy.end())  // dissociation & avoid more than one contact
							//if ( L < LTHR_TENSION_CONDENSIN )
						{
							vecCondensinSiteBusy.push_back(bId1);
							vecCondensinSiteBusy.push_back(bId2);

														//cout << "forming condensin association!" << endl;
							dVec veloc1 = vecNode[bId1].get_veloc(), veloc2 = vecNode[bId2].get_veloc(), veloc1_new, veloc2_new;
							double L0 = L0_TENSION_CONDENSIN, L_inv;
							double k = K_TENSION_CONDENSIN;
							dVec unitVec12 = {0,0,0};
							if (L > 0)
							{
								L_inv = 1./L;  // 1/nm
								unitVec12 = {(coord2.x-coord1.x)*L_inv, (coord2.y-coord1.y)*L_inv, (coord2.z-coord1.z)*L_inv};
							}

							if ( abs(L - L0) > L_DIFF_RESOLUTION) // non-zero tension
							{
								double v_ten_con_size = k*(L-L0)*gamma_inv*1E-3;
								dVec v_ten_con = {v_ten_con_size*unitVec12.x, v_ten_con_size*unitVec12.y, v_ten_con_size*unitVec12.z};
								veloc1_new = {veloc1.x + v_ten_con.x, veloc1.y + v_ten_con.y, veloc1.z + v_ten_con.z};
								veloc2_new = {veloc2.x - v_ten_con.x, veloc2.y - v_ten_con.y, veloc2.z - v_ten_con.z};
														// update acceleration
								vecNode[bId1].set_veloc(veloc1_new);
								vecNode[bId2].set_veloc(veloc2_new);
														// update potential & total energy
														//double ePot12 = 0.5*k*(L-L0)*(L-L0);
														//energy["potential"] += ePot12;
														//energy["total"]     += ePot12;
							}
							capturer.push_back( {bId1, bId2});
							
						}
					}
					//
					if (DIFFUSION_CAPTURE == "n-valent")
					{
							
							if (   L < LTHR_TENSION_CONDENSIN
								&& count(vecCondensinSiteBusy.begin(), vecCondensinSiteBusy.end(), bId1) < DC_VALENCE
								&& count(vecCondensinSiteBusy.begin(), vecCondensinSiteBusy.end(), bId2) < DC_VALENCE)  // dissociation & avoid more than 4 contactS
								//if ( L < LTHR_TENSION_CONDENSIN )
							{
								vecCondensinSiteBusy.push_back(bId1);
								vecCondensinSiteBusy.push_back(bId2);

															//cout << "forming condensin association!" << endl;
								dVec veloc1 = vecNode[bId1].get_veloc(), veloc2 = vecNode[bId2].get_veloc(), veloc1_new, veloc2_new;
								double L0 = L0_TENSION_CONDENSIN, L_inv;
								double k = K_TENSION_CONDENSIN;
								dVec unitVec12 = {0,0,0};
								if (L > 0)
								{
									L_inv = 1./L;  // 1/nm
									unitVec12 = {(coord2.x-coord1.x)*L_inv, (coord2.y-coord1.y)*L_inv, (coord2.z-coord1.z)*L_inv};
								}

								if ( abs(L - L0) > L_DIFF_RESOLUTION) // non-zero tension
								{
									double v_ten_con_size = k*(L-L0)*gamma_inv*1E-3;
									dVec v_ten_con = {v_ten_con_size*unitVec12.x, v_ten_con_size*unitVec12.y, v_ten_con_size*unitVec12.z};
									veloc1_new = {veloc1.x + v_ten_con.x, veloc1.y + v_ten_con.y, veloc1.z + v_ten_con.z};
									veloc2_new = {veloc2.x - v_ten_con.x, veloc2.y - v_ten_con.y, veloc2.z - v_ten_con.z};
															// update acceleration
									vecNode[bId1].set_veloc(veloc1_new);
									vecNode[bId2].set_veloc(veloc2_new);
															// update potential & total energy
															//double ePot12 = 0.5*k*(L-L0)*(L-L0);
															//energy["potential"] += ePot12;
															//energy["total"]     += ePot12;
								}
								capturer.push_back( {bId1, bId2});
								
							}
					}					
										
					
					
					//
					if (DIFFUSION_CAPTURE == "multivalent")
					{
						
						//if (   L < LTHR_TENSION_CONDENSIN
							//&& find(vecCondensinSiteBusy.begin(), vecCondensinSiteBusy.end(), bId1) == vecCondensinSiteBusy.end()
							//&& find(vecCondensinSiteBusy.begin(), vecCondensinSiteBusy.end(), bId2) == vecCondensinSiteBusy.end())  // dissociation & avoid more than one contact
						if ( L < LTHR_TENSION_CONDENSIN )
						{
							vecCondensinSiteBusy.push_back(bId1);
							vecCondensinSiteBusy.push_back(bId2);

														//cout << "forming condensin association!" << endl;
							dVec veloc1 = vecNode[bId1].get_veloc(), veloc2 = vecNode[bId2].get_veloc(), veloc1_new, veloc2_new;
							double L0 = L0_TENSION_CONDENSIN, L_inv;
							double k = K_TENSION_CONDENSIN;
							dVec unitVec12 = {0,0,0};
							if (L > 0)
							{
								L_inv = 1./L;  // 1/nm
								unitVec12 = {(coord2.x-coord1.x)*L_inv, (coord2.y-coord1.y)*L_inv, (coord2.z-coord1.z)*L_inv};
							}

							if ( abs(L - L0) > L_DIFF_RESOLUTION) // non-zero tension
							{
								double v_ten_con_size = k*(L-L0)*gamma_inv*1E-3;
								dVec v_ten_con = {v_ten_con_size*unitVec12.x, v_ten_con_size*unitVec12.y, v_ten_con_size*unitVec12.z};
								veloc1_new = {veloc1.x + v_ten_con.x, veloc1.y + v_ten_con.y, veloc1.z + v_ten_con.z};
								veloc2_new = {veloc2.x - v_ten_con.x, veloc2.y - v_ten_con.y, veloc2.z - v_ten_con.z};
														// update acceleration
								vecNode[bId1].set_veloc(veloc1_new);
								vecNode[bId2].set_veloc(veloc2_new);
														// update potential & total energy
														//double ePot12 = 0.5*k*(L-L0)*(L-L0);
														//energy["potential"] += ePot12;
														//energy["total"]     += ePot12;
							}
							capturer.push_back( {bId1, bId2});
							
						}
					}					
					
				}
			}
		}
	}	
	
	
	if (DIFFUSION_CAPTURE_ON == true && flagIsPrep == false && SWITCH_INTERPHASE_MITOSIS == 0 && DC_MECHANISM == "CSB-BF")
	{
		
		cout << "CSB-BF"<< endl;
		vector<int> vecCondensinSiteBusy; // contain bead ids already interacting with another
		capturer.clear();
		hit_valence.clear();

		for (int i = 0; i < vecCondensinSite.size(); i++)
		{
			vector<int> vecAllFBinderBeadIds;
			for (vector<Binder>::iterator it = vecBinder.begin(); it != vecBinder.end(); it++) 
			{
				vector<int> vecBeadIds = it->get_vecBeadIds();
				int bIdF = vecBeadIds[0];
				vecAllFBinderBeadIds.push_back(vecBeadIds[0]);
			}

			for (int j = 0; j < vecAllFBinderBeadIds.size(); j++) 
			{
				if ((rand() / double(RAND_MAX)) > P_DISSOCIATE_CONDENSIN) 
				{
					int bId1 = vecCondensinSite[i], bId2 = vecAllFBinderBeadIds[j];
					dVec coord1 = vecNode[bId1].get_coord(), coord2 = vecNode[bId2].get_coord();
					double L = dVecDist(coord1, coord2);

					if (DIFFUSION_CAPTURE == "probabilistic")
					{
						if ((rand() / double(RAND_MAX)) > P_DC_PROBAB)
						{ 
				// multivalent diffusion capture
							if (L < LTHR_TENSION_CONDENSIN)
							{
								vecCondensinSiteBusy.push_back(bId1);
								vecCondensinSiteBusy.push_back(bId2);

								//cout << "forming condensin association!" << endl;
								dVec veloc1 = vecNode[bId1].get_veloc(), veloc2 = vecNode[bId2].get_veloc(), veloc1_new, veloc2_new;
								double L0 = L0_TENSION_CONDENSIN, L_inv;
								double k = K_TENSION_CONDENSIN;
								dVec unitVec12 = { 0, 0, 0 };
								if (L > 0)
								{
									L_inv = 1. / L;  // 1/nm
									unitVec12 = {(coord2.x-coord1.x)*L_inv, (coord2.y-coord1.y)*L_inv, (coord2.z-coord1.z)*L_inv};
								}

								if (abs(L - L0) > L_DIFF_RESOLUTION) // non-zero tension
								{
									double v_ten_con_size = k * (L - L0)* gamma_inv * 1E-3;
									dVec v_ten_con = { v_ten_con_size* unitVec12.x, v_ten_con_size* unitVec12.y, v_ten_con_size* unitVec12.z };
									veloc1_new = {veloc1.x + v_ten_con.x, veloc1.y + v_ten_con.y, veloc1.z + v_ten_con.z};
									veloc2_new = {veloc2.x - v_ten_con.x, veloc2.y - v_ten_con.y, veloc2.z - v_ten_con.z};
									// update acceleration
									vecNode[bId1].set_veloc(veloc1_new);
									vecNode[bId2].set_veloc(veloc2_new);
									// update potential & total energy
									//double ePot12 = 0.5*k*(L-L0)*(L-L0);
									//energy["potential"] += ePot12;
									//energy["total"]     += ePot12;
								}
								capturer.push_back( {bId1, bId2});
							}
							
						}
  													
			    // bivalent diffusion capture
						else 
						{
							if (L < LTHR_TENSION_CONDENSIN && find(vecCondensinSiteBusy.begin(), vecCondensinSiteBusy.end(), bId1) == vecCondensinSiteBusy.end() && find(vecCondensinSiteBusy.begin(), vecCondensinSiteBusy.end(), bId2) == vecCondensinSiteBusy.end())
							//&& find(vecCondensinSiteBusy.begin(),vecCondensinSiteBusy.end(), bId1) == vecCondensinSiteBusy.end()
							//&& find(vecCondensinSiteBusy.begin(),vecCondensinSiteBusy.end(), bId2) == vecCondensinSiteBusy.end()) // dissociation & avoid more than one contact
							{
								vecCondensinSiteBusy.push_back(bId1);
								vecCondensinSiteBusy.push_back(bId2);

								//cout << "forming condensin association!" << endl;
								dVec veloc1 = vecNode[bId1].get_veloc(), veloc2 = vecNode[bId2].get_veloc(),veloc1_new, veloc2_new;
								double L0 = L0_TENSION_CONDENSIN, L_inv;
								double k = K_TENSION_CONDENSIN;
								dVec unitVec12 = { 0, 0, 0 };
								if (L > 0)
								{
									L_inv = 1. / L;  // 1/nm
									unitVec12 = {(coord2.x-coord1.x)*L_inv, (coord2.y-coord1.y)*L_inv, (coord2.z-coord1.z)*L_inv};
								}

								if ( abs(L - L0) > L_DIFF_RESOLUTION) // non-zero tension
								{
									double v_ten_con_size = k*(L-L0)*gamma_inv*1E-3;
									dVec v_ten_con = {v_ten_con_size*unitVec12.x, v_ten_con_size*unitVec12.y, v_ten_con_size*unitVec12.z};
									veloc1_new = {veloc1.x + v_ten_con.x, veloc1.y + v_ten_con.y, veloc1.z + v_ten_con.z};
									veloc2_new = {veloc2.x - v_ten_con.x, veloc2.y - v_ten_con.y, veloc2.z - v_ten_con.z};
									// update acceleration
									vecNode[bId1].set_veloc(veloc1_new);
									vecNode[bId2].set_veloc(veloc2_new);
									// update potential & total energy
									//double ePot12 = 0.5*k*(L-L0)*(L-L0);
									//energy["potential"] += ePot12;
									//energy["total"]     += ePot12;
								}
							}	capturer.push_back( {bId1, bId2});
							
						}
						
					}

//					if (DIFFUSION_CAPTURE_PROBAB == true)

					if (DIFFUSION_CAPTURE == "bivalent")
					{
						
						if (   L < LTHR_TENSION_CONDENSIN
							&& find(vecCondensinSiteBusy.begin(), vecCondensinSiteBusy.end(), bId1) == vecCondensinSiteBusy.end()
							&& find(vecCondensinSiteBusy.begin(), vecCondensinSiteBusy.end(), bId2) == vecCondensinSiteBusy.end())  // dissociation & avoid more than one contact
							//if ( L < LTHR_TENSION_CONDENSIN )
						{
							vecCondensinSiteBusy.push_back(bId1);
							vecCondensinSiteBusy.push_back(bId2);

														//cout << "forming condensin association!" << endl;
							dVec veloc1 = vecNode[bId1].get_veloc(), veloc2 = vecNode[bId2].get_veloc(), veloc1_new, veloc2_new;
							double L0 = L0_TENSION_CONDENSIN, L_inv;
							double k = K_TENSION_CONDENSIN;
							dVec unitVec12 = {0,0,0};
							if (L > 0)
							{
								L_inv = 1./L;  // 1/nm
								unitVec12 = {(coord2.x-coord1.x)*L_inv, (coord2.y-coord1.y)*L_inv, (coord2.z-coord1.z)*L_inv};
							}

							if ( abs(L - L0) > L_DIFF_RESOLUTION) // non-zero tension
							{
								double v_ten_con_size = k*(L-L0)*gamma_inv*1E-3;
								dVec v_ten_con = {v_ten_con_size*unitVec12.x, v_ten_con_size*unitVec12.y, v_ten_con_size*unitVec12.z};
								veloc1_new = {veloc1.x + v_ten_con.x, veloc1.y + v_ten_con.y, veloc1.z + v_ten_con.z};
								veloc2_new = {veloc2.x - v_ten_con.x, veloc2.y - v_ten_con.y, veloc2.z - v_ten_con.z};
														// update acceleration
								vecNode[bId1].set_veloc(veloc1_new);
								vecNode[bId2].set_veloc(veloc2_new);
														// update potential & total energy
														//double ePot12 = 0.5*k*(L-L0)*(L-L0);
														//energy["potential"] += ePot12;
														//energy["total"]     += ePot12;
							}
							capturer.push_back( {bId1, bId2});
							
						}
					}
					//
					if (DIFFUSION_CAPTURE == "n-valent")
					{
							
							if (   L < LTHR_TENSION_CONDENSIN
								&& count(vecCondensinSiteBusy.begin(), vecCondensinSiteBusy.end(), bId1) < DC_VALENCE
								&& count(vecCondensinSiteBusy.begin(), vecCondensinSiteBusy.end(), bId2) < DC_VALENCE)  // dissociation & avoid more than 4 contactS
								//if ( L < LTHR_TENSION_CONDENSIN )
							{
								vecCondensinSiteBusy.push_back(bId1);
								vecCondensinSiteBusy.push_back(bId2);

															//cout << "forming condensin association!" << endl;
								dVec veloc1 = vecNode[bId1].get_veloc(), veloc2 = vecNode[bId2].get_veloc(), veloc1_new, veloc2_new;
								double L0 = L0_TENSION_CONDENSIN, L_inv;
								double k = K_TENSION_CONDENSIN;
								dVec unitVec12 = {0,0,0};
								if (L > 0)
								{
									L_inv = 1./L;  // 1/nm
									unitVec12 = {(coord2.x-coord1.x)*L_inv, (coord2.y-coord1.y)*L_inv, (coord2.z-coord1.z)*L_inv};
								}

								if ( abs(L - L0) > L_DIFF_RESOLUTION) // non-zero tension
								{
									double v_ten_con_size = k*(L-L0)*gamma_inv*1E-3;
									dVec v_ten_con = {v_ten_con_size*unitVec12.x, v_ten_con_size*unitVec12.y, v_ten_con_size*unitVec12.z};
									veloc1_new = {veloc1.x + v_ten_con.x, veloc1.y + v_ten_con.y, veloc1.z + v_ten_con.z};
									veloc2_new = {veloc2.x - v_ten_con.x, veloc2.y - v_ten_con.y, veloc2.z - v_ten_con.z};
															// update acceleration
									vecNode[bId1].set_veloc(veloc1_new);
									vecNode[bId2].set_veloc(veloc2_new);
															// update potential & total energy
															//double ePot12 = 0.5*k*(L-L0)*(L-L0);
															//energy["potential"] += ePot12;
															//energy["total"]     += ePot12;
								}
								capturer.push_back( {bId1, bId2});
								
							}
					}					
										
					
					
					//
					if (DIFFUSION_CAPTURE == "multivalent")
					{
						
						//if (   L < LTHR_TENSION_CONDENSIN
							//&& find(vecCondensinSiteBusy.begin(), vecCondensinSiteBusy.end(), bId1) == vecCondensinSiteBusy.end()
							//&& find(vecCondensinSiteBusy.begin(), vecCondensinSiteBusy.end(), bId2) == vecCondensinSiteBusy.end())  // dissociation & avoid more than one contact
						if ( L < LTHR_TENSION_CONDENSIN )
						{
							vecCondensinSiteBusy.push_back(bId1);
							vecCondensinSiteBusy.push_back(bId2);

														//cout << "forming condensin association!" << endl;
							dVec veloc1 = vecNode[bId1].get_veloc(), veloc2 = vecNode[bId2].get_veloc(), veloc1_new, veloc2_new;
							double L0 = L0_TENSION_CONDENSIN, L_inv;
							double k = K_TENSION_CONDENSIN;
							dVec unitVec12 = {0,0,0};
							if (L > 0)
							{
								L_inv = 1./L;  // 1/nm
								unitVec12 = {(coord2.x-coord1.x)*L_inv, (coord2.y-coord1.y)*L_inv, (coord2.z-coord1.z)*L_inv};
							}

							if ( abs(L - L0) > L_DIFF_RESOLUTION) // non-zero tension
							{
								double v_ten_con_size = k*(L-L0)*gamma_inv*1E-3;
								dVec v_ten_con = {v_ten_con_size*unitVec12.x, v_ten_con_size*unitVec12.y, v_ten_con_size*unitVec12.z};
								veloc1_new = {veloc1.x + v_ten_con.x, veloc1.y + v_ten_con.y, veloc1.z + v_ten_con.z};
								veloc2_new = {veloc2.x - v_ten_con.x, veloc2.y - v_ten_con.y, veloc2.z - v_ten_con.z};
														// update acceleration
								vecNode[bId1].set_veloc(veloc1_new);
								vecNode[bId2].set_veloc(veloc2_new);
														// update potential & total energy
														//double ePot12 = 0.5*k*(L-L0)*(L-L0);
														//energy["potential"] += ePot12;
														//energy["total"]     += ePot12;
							}
							capturer.push_back( {bId1, bId2});
							
						}
					}					
					
				}
			}
		}
	}

	if (DIFFUSION_CAPTURE_ON == true && flagIsPrep == false && SWITCH_INTERPHASE_MITOSIS == 0 && DC_MECHANISM == "CSB-CSB")
	{
		
		cout << "CSB-CSB" << endl;
		vector<int> vecCondensinSiteBusy; // contain bead ids already interacting with another
		capturer.clear();
		hit_valence.clear();

		for (int i = 0; i < vecCondensinSite.size(); i++)
		{
			for (int j = i+1; j < vecCondensinSite.size(); j++) 
			{
				if ((rand() / double(RAND_MAX)) > P_DISSOCIATE_CONDENSIN) 
				{
					int bId1 = vecCondensinSite[i], bId2 = vecCondensinSite[j];
					dVec coord1 = vecNode[bId1].get_coord(), coord2 = vecNode[bId2].get_coord();
					double L = dVecDist(coord1, coord2);

					if (DIFFUSION_CAPTURE == "probabilistic")
					{
						if ((rand() / double(RAND_MAX)) > P_DC_PROBAB)
						{ 
				// multivalent diffusion capture
							if (L < LTHR_TENSION_CONDENSIN)
							{
								vecCondensinSiteBusy.push_back(bId1);
								vecCondensinSiteBusy.push_back(bId2);

								//cout << "forming condensin association!" << endl;
								dVec veloc1 = vecNode[bId1].get_veloc(), veloc2 = vecNode[bId2].get_veloc(), veloc1_new, veloc2_new;
								double L0 = L0_TENSION_CONDENSIN, L_inv;
								double k = K_TENSION_CONDENSIN;
								dVec unitVec12 = { 0, 0, 0 };
								if (L > 0)
								{
									L_inv = 1. / L;  // 1/nm
									unitVec12 = {(coord2.x-coord1.x)*L_inv, (coord2.y-coord1.y)*L_inv, (coord2.z-coord1.z)*L_inv};
								}

								if (abs(L - L0) > L_DIFF_RESOLUTION) // non-zero tension
								{
									double v_ten_con_size = k * (L - L0)* gamma_inv * 1E-3;
									dVec v_ten_con = { v_ten_con_size* unitVec12.x, v_ten_con_size* unitVec12.y, v_ten_con_size* unitVec12.z };
									veloc1_new = {veloc1.x + v_ten_con.x, veloc1.y + v_ten_con.y, veloc1.z + v_ten_con.z};
									veloc2_new = {veloc2.x - v_ten_con.x, veloc2.y - v_ten_con.y, veloc2.z - v_ten_con.z};
									// update acceleration
									vecNode[bId1].set_veloc(veloc1_new);
									vecNode[bId2].set_veloc(veloc2_new);
									// update potential & total energy
									//double ePot12 = 0.5*k*(L-L0)*(L-L0);
									//energy["potential"] += ePot12;
									//energy["total"]     += ePot12;
								}
								capturer.push_back( {bId1, bId2});
							}
						}
  													
			    // bivalent diffusion capture
						else 
						{
							if (L < LTHR_TENSION_CONDENSIN && find(vecCondensinSiteBusy.begin(), vecCondensinSiteBusy.end(), bId1) == vecCondensinSiteBusy.end() && find(vecCondensinSiteBusy.begin(), vecCondensinSiteBusy.end(), bId2) == vecCondensinSiteBusy.end())
							//&& find(vecCondensinSiteBusy.begin(),vecCondensinSiteBusy.end(), bId1) == vecCondensinSiteBusy.end()
							//&& find(vecCondensinSiteBusy.begin(),vecCondensinSiteBusy.end(), bId2) == vecCondensinSiteBusy.end()) // dissociation & avoid more than one contact
							{
								vecCondensinSiteBusy.push_back(bId1);
								vecCondensinSiteBusy.push_back(bId2);

								//cout << "forming condensin association!" << endl;
								dVec veloc1 = vecNode[bId1].get_veloc(), veloc2 = vecNode[bId2].get_veloc(),veloc1_new, veloc2_new;
								double L0 = L0_TENSION_CONDENSIN, L_inv;
								double k = K_TENSION_CONDENSIN;
								dVec unitVec12 = { 0, 0, 0 };
								if (L > 0)
								{
									L_inv = 1. / L;  // 1/nm
									unitVec12 = {(coord2.x-coord1.x)*L_inv, (coord2.y-coord1.y)*L_inv, (coord2.z-coord1.z)*L_inv};
								}

								if ( abs(L - L0) > L_DIFF_RESOLUTION) // non-zero tension
								{
									double v_ten_con_size = k*(L-L0)*gamma_inv*1E-3;
									dVec v_ten_con = {v_ten_con_size*unitVec12.x, v_ten_con_size*unitVec12.y, v_ten_con_size*unitVec12.z};
									veloc1_new = {veloc1.x + v_ten_con.x, veloc1.y + v_ten_con.y, veloc1.z + v_ten_con.z};
									veloc2_new = {veloc2.x - v_ten_con.x, veloc2.y - v_ten_con.y, veloc2.z - v_ten_con.z};
									// update acceleration
									vecNode[bId1].set_veloc(veloc1_new);
									vecNode[bId2].set_veloc(veloc2_new);
									// update potential & total energy
									//double ePot12 = 0.5*k*(L-L0)*(L-L0);
									//energy["potential"] += ePot12;
									//energy["total"]     += ePot12;
								}
							}	capturer.push_back( {bId1, bId2});
							
						}
						
					}

//					if (DIFFUSION_CAPTURE_PROBAB == true)

					if (DIFFUSION_CAPTURE == "bivalent")
					{
						
						if (   L < LTHR_TENSION_CONDENSIN
							&& find(vecCondensinSiteBusy.begin(), vecCondensinSiteBusy.end(), bId1) == vecCondensinSiteBusy.end()
							&& find(vecCondensinSiteBusy.begin(), vecCondensinSiteBusy.end(), bId2) == vecCondensinSiteBusy.end())  // dissociation & avoid more than one contact
							//if ( L < LTHR_TENSION_CONDENSIN )
						{
							vecCondensinSiteBusy.push_back(bId1);
							vecCondensinSiteBusy.push_back(bId2);

														//cout << "forming condensin association!" << endl;
							dVec veloc1 = vecNode[bId1].get_veloc(), veloc2 = vecNode[bId2].get_veloc(), veloc1_new, veloc2_new;
							double L0 = L0_TENSION_CONDENSIN, L_inv;
							double k = K_TENSION_CONDENSIN;
							dVec unitVec12 = {0,0,0};
							if (L > 0)
							{
								L_inv = 1./L;  // 1/nm
								unitVec12 = {(coord2.x-coord1.x)*L_inv, (coord2.y-coord1.y)*L_inv, (coord2.z-coord1.z)*L_inv};
							}

							if ( abs(L - L0) > L_DIFF_RESOLUTION) // non-zero tension
							{
								double v_ten_con_size = k*(L-L0)*gamma_inv*1E-3;
								dVec v_ten_con = {v_ten_con_size*unitVec12.x, v_ten_con_size*unitVec12.y, v_ten_con_size*unitVec12.z};
								veloc1_new = {veloc1.x + v_ten_con.x, veloc1.y + v_ten_con.y, veloc1.z + v_ten_con.z};
								veloc2_new = {veloc2.x - v_ten_con.x, veloc2.y - v_ten_con.y, veloc2.z - v_ten_con.z};
														// update acceleration
								vecNode[bId1].set_veloc(veloc1_new);
								vecNode[bId2].set_veloc(veloc2_new);
														// update potential & total energy
														//double ePot12 = 0.5*k*(L-L0)*(L-L0);
														//energy["potential"] += ePot12;
														//energy["total"]     += ePot12;
							}
							capturer.push_back( {bId1, bId2});
							
						}
					}
					//
					if (DIFFUSION_CAPTURE == "n-valent")
					{
							
							if (   L < LTHR_TENSION_CONDENSIN
								&& count(vecCondensinSiteBusy.begin(), vecCondensinSiteBusy.end(), bId1) < DC_VALENCE
								&& count(vecCondensinSiteBusy.begin(), vecCondensinSiteBusy.end(), bId2) < DC_VALENCE)  // dissociation & avoid more than 4 contactS
								//if ( L < LTHR_TENSION_CONDENSIN )
							{
								vecCondensinSiteBusy.push_back(bId1);
								vecCondensinSiteBusy.push_back(bId2);

															//cout << "forming condensin association!" << endl;
								dVec veloc1 = vecNode[bId1].get_veloc(), veloc2 = vecNode[bId2].get_veloc(), veloc1_new, veloc2_new;
								double L0 = L0_TENSION_CONDENSIN, L_inv;
								double k = K_TENSION_CONDENSIN;
								dVec unitVec12 = {0,0,0};
								if (L > 0)
								{
									L_inv = 1./L;  // 1/nm
									unitVec12 = {(coord2.x-coord1.x)*L_inv, (coord2.y-coord1.y)*L_inv, (coord2.z-coord1.z)*L_inv};
								}

								if ( abs(L - L0) > L_DIFF_RESOLUTION) // non-zero tension
								{
									double v_ten_con_size = k*(L-L0)*gamma_inv*1E-3;
									dVec v_ten_con = {v_ten_con_size*unitVec12.x, v_ten_con_size*unitVec12.y, v_ten_con_size*unitVec12.z};
									veloc1_new = {veloc1.x + v_ten_con.x, veloc1.y + v_ten_con.y, veloc1.z + v_ten_con.z};
									veloc2_new = {veloc2.x - v_ten_con.x, veloc2.y - v_ten_con.y, veloc2.z - v_ten_con.z};
															// update acceleration
									vecNode[bId1].set_veloc(veloc1_new);
									vecNode[bId2].set_veloc(veloc2_new);
															// update potential & total energy
															//double ePot12 = 0.5*k*(L-L0)*(L-L0);
															//energy["potential"] += ePot12;
															//energy["total"]     += ePot12;
								}
								capturer.push_back( {bId1, bId2});
								
							}
					}					
					
					//
					if (DIFFUSION_CAPTURE == "multivalent")
					{
						
						//if (   L < LTHR_TENSION_CONDENSIN
							//&& find(vecCondensinSiteBusy.begin(), vecCondensinSiteBusy.end(), bId1) == vecCondensinSiteBusy.end()
							//&& find(vecCondensinSiteBusy.begin(), vecCondensinSiteBusy.end(), bId2) == vecCondensinSiteBusy.end())  // dissociation & avoid more than one contact
						if ( L < LTHR_TENSION_CONDENSIN )
						{
							vecCondensinSiteBusy.push_back(bId1);
							vecCondensinSiteBusy.push_back(bId2);

														//cout << "forming condensin association!" << endl;
							dVec veloc1 = vecNode[bId1].get_veloc(), veloc2 = vecNode[bId2].get_veloc(), veloc1_new, veloc2_new;
							double L0 = L0_TENSION_CONDENSIN, L_inv;
							double k = K_TENSION_CONDENSIN;
							dVec unitVec12 = {0,0,0};
							if (L > 0)
							{
								L_inv = 1./L;  // 1/nm
								unitVec12 = {(coord2.x-coord1.x)*L_inv, (coord2.y-coord1.y)*L_inv, (coord2.z-coord1.z)*L_inv};
							}

							if ( abs(L - L0) > L_DIFF_RESOLUTION) // non-zero tension
							{
								double v_ten_con_size = k*(L-L0)*gamma_inv*1E-3;
								dVec v_ten_con = {v_ten_con_size*unitVec12.x, v_ten_con_size*unitVec12.y, v_ten_con_size*unitVec12.z};
								veloc1_new = {veloc1.x + v_ten_con.x, veloc1.y + v_ten_con.y, veloc1.z + v_ten_con.z};
								veloc2_new = {veloc2.x - v_ten_con.x, veloc2.y - v_ten_con.y, veloc2.z - v_ten_con.z};
														// update acceleration
								vecNode[bId1].set_veloc(veloc1_new);
								vecNode[bId2].set_veloc(veloc2_new);
														// update potential & total energy
														//double ePot12 = 0.5*k*(L-L0)*(L-L0);
														//energy["potential"] += ePot12;
														//energy["total"]     += ePot12;
							}
							capturer.push_back( {bId1, bId2});
							
						}
					}					
					
				}
			}
		}
	}
	
	/////
	if (DIFFUSION_CAPTURE_ON == true && flagIsPrep == false && SWITCH_INTERPHASE_MITOSIS == 0 && DC_MECHANISM == "BAll-BAll")	
	{
			cout << "BAll-BAll" << endl;
			vector<int> vecCondensinSiteBusy; // contain bead ids already interacting with another
			capturer.clear();
			hit_valence.clear();

			vector<int> vecAllFBinderBeadIds;
			vector<int> vecAllRBinderBeadIds;
			vector<int> vecAllBinderBeadIds;
			for (vector<Binder>::iterator it = vecBinder.begin(); it != vecBinder.end(); it++) 
			{
				vector<int> vecBeadIds = it->get_vecBeadIds();
				int bIdF = vecBeadIds[0];
				vecAllFBinderBeadIds.push_back(vecBeadIds[0]);
				int bIdR = vecBeadIds[1];
				vecAllRBinderBeadIds.push_back(vecBeadIds[1]);
				vecAllBinderBeadIds.push_back(vecBeadIds[0]);
				vecAllBinderBeadIds.push_back(vecBeadIds[1]);					
			}
			
			for (int i = 0; i < vecAllBinderBeadIds.size(); i++)
			{
				for (int j = i+1; j < vecAllBinderBeadIds.size(); j++) 
				{
					if ((rand() / double(RAND_MAX)) > P_DISSOCIATE_CONDENSIN) 
					{
						int bId1 = vecAllBinderBeadIds[i], bId2 = vecAllBinderBeadIds[j];
						dVec coord1 = vecNode[bId1].get_coord(), coord2 = vecNode[bId2].get_coord();
						double L = dVecDist(coord1, coord2);

						if (DIFFUSION_CAPTURE == "probabilistic")
						{
							if ((rand() / double(RAND_MAX)) > P_DC_PROBAB)
							{ 
					// multivalent diffusion capture
								if (L < LTHR_TENSION_CONDENSIN)
								{
									vecCondensinSiteBusy.push_back(bId1);
									vecCondensinSiteBusy.push_back(bId2);

									//cout << "forming condensin association!" << endl;
									dVec veloc1 = vecNode[bId1].get_veloc(), veloc2 = vecNode[bId2].get_veloc(), veloc1_new, veloc2_new;
									double L0 = L0_TENSION_CONDENSIN, L_inv;
									double k = K_TENSION_CONDENSIN;
									dVec unitVec12 = { 0, 0, 0 };
									if (L > 0)
									{
										L_inv = 1. / L;  // 1/nm
										unitVec12 = {(coord2.x-coord1.x)*L_inv, (coord2.y-coord1.y)*L_inv, (coord2.z-coord1.z)*L_inv};
									}

									if (abs(L - L0) > L_DIFF_RESOLUTION) // non-zero tension
									{
										double v_ten_con_size = k * (L - L0)* gamma_inv * 1E-3;
										dVec v_ten_con = { v_ten_con_size* unitVec12.x, v_ten_con_size* unitVec12.y, v_ten_con_size* unitVec12.z };
										veloc1_new = {veloc1.x + v_ten_con.x, veloc1.y + v_ten_con.y, veloc1.z + v_ten_con.z};
										veloc2_new = {veloc2.x - v_ten_con.x, veloc2.y - v_ten_con.y, veloc2.z - v_ten_con.z};
										// update acceleration
										vecNode[bId1].set_veloc(veloc1_new);
										vecNode[bId2].set_veloc(veloc2_new);
										// update potential & total energy
										//double ePot12 = 0.5*k*(L-L0)*(L-L0);
										//energy["potential"] += ePot12;
										//energy["total"]     += ePot12;
									}
									capturer.push_back( {bId1, bId2});
								}
								
							}
														
					// bivalent diffusion capture
							else 
							{
								if (L < LTHR_TENSION_CONDENSIN && find(vecCondensinSiteBusy.begin(), vecCondensinSiteBusy.end(), bId1) == vecCondensinSiteBusy.end() && find(vecCondensinSiteBusy.begin(), vecCondensinSiteBusy.end(), bId2) == vecCondensinSiteBusy.end())
								//&& find(vecCondensinSiteBusy.begin(),vecCondensinSiteBusy.end(), bId1) == vecCondensinSiteBusy.end()
								//&& find(vecCondensinSiteBusy.begin(),vecCondensinSiteBusy.end(), bId2) == vecCondensinSiteBusy.end()) // dissociation & avoid more than one contact
								{
									vecCondensinSiteBusy.push_back(bId1);
									vecCondensinSiteBusy.push_back(bId2);

									//cout << "forming condensin association!" << endl;
									dVec veloc1 = vecNode[bId1].get_veloc(), veloc2 = vecNode[bId2].get_veloc(),veloc1_new, veloc2_new;
									double L0 = L0_TENSION_CONDENSIN, L_inv;
									double k = K_TENSION_CONDENSIN;
									dVec unitVec12 = { 0, 0, 0 };
									if (L > 0)
									{
										L_inv = 1. / L;  // 1/nm
										unitVec12 = {(coord2.x-coord1.x)*L_inv, (coord2.y-coord1.y)*L_inv, (coord2.z-coord1.z)*L_inv};
									}

									if ( abs(L - L0) > L_DIFF_RESOLUTION) // non-zero tension
									{
										double v_ten_con_size = k*(L-L0)*gamma_inv*1E-3;
										dVec v_ten_con = {v_ten_con_size*unitVec12.x, v_ten_con_size*unitVec12.y, v_ten_con_size*unitVec12.z};
										veloc1_new = {veloc1.x + v_ten_con.x, veloc1.y + v_ten_con.y, veloc1.z + v_ten_con.z};
										veloc2_new = {veloc2.x - v_ten_con.x, veloc2.y - v_ten_con.y, veloc2.z - v_ten_con.z};
										// update acceleration
										vecNode[bId1].set_veloc(veloc1_new);
										vecNode[bId2].set_veloc(veloc2_new);
										// update potential & total energy
										//double ePot12 = 0.5*k*(L-L0)*(L-L0);
										//energy["potential"] += ePot12;
										//energy["total"]     += ePot12;
									}
								}	capturer.push_back( {bId1, bId2});
								
							}
							
						}

	//					if (DIFFUSION_CAPTURE_PROBAB == true)

						if (DIFFUSION_CAPTURE == "bivalent")
						{
							
							if (   L < LTHR_TENSION_CONDENSIN
								&& find(vecCondensinSiteBusy.begin(), vecCondensinSiteBusy.end(), bId1) == vecCondensinSiteBusy.end()
								&& find(vecCondensinSiteBusy.begin(), vecCondensinSiteBusy.end(), bId2) == vecCondensinSiteBusy.end())  // dissociation & avoid more than one contact
								//if ( L < LTHR_TENSION_CONDENSIN )
							{
								vecCondensinSiteBusy.push_back(bId1);
								vecCondensinSiteBusy.push_back(bId2);

															//cout << "forming condensin association!" << endl;
								dVec veloc1 = vecNode[bId1].get_veloc(), veloc2 = vecNode[bId2].get_veloc(), veloc1_new, veloc2_new;
								double L0 = L0_TENSION_CONDENSIN, L_inv;
								double k = K_TENSION_CONDENSIN;
								dVec unitVec12 = {0,0,0};
								if (L > 0)
								{
									L_inv = 1./L;  // 1/nm
									unitVec12 = {(coord2.x-coord1.x)*L_inv, (coord2.y-coord1.y)*L_inv, (coord2.z-coord1.z)*L_inv};
								}

								if ( abs(L - L0) > L_DIFF_RESOLUTION) // non-zero tension
								{
									double v_ten_con_size = k*(L-L0)*gamma_inv*1E-3;
									dVec v_ten_con = {v_ten_con_size*unitVec12.x, v_ten_con_size*unitVec12.y, v_ten_con_size*unitVec12.z};
									veloc1_new = {veloc1.x + v_ten_con.x, veloc1.y + v_ten_con.y, veloc1.z + v_ten_con.z};
									veloc2_new = {veloc2.x - v_ten_con.x, veloc2.y - v_ten_con.y, veloc2.z - v_ten_con.z};
															// update acceleration
									vecNode[bId1].set_veloc(veloc1_new);
									vecNode[bId2].set_veloc(veloc2_new);
															// update potential & total energy
															//double ePot12 = 0.5*k*(L-L0)*(L-L0);
															//energy["potential"] += ePot12;
															//energy["total"]     += ePot12;
								}
								capturer.push_back( {bId1, bId2});
								
							}
						}
	/////////////////////////////////////////////// 30-10-2018 tetravalent //////////////////////					
						
						if (DIFFUSION_CAPTURE == "n-valent")
						{
							
							if (   L < LTHR_TENSION_CONDENSIN
								&& count(vecCondensinSiteBusy.begin(), vecCondensinSiteBusy.end(), bId1) < DC_VALENCE
								&& count(vecCondensinSiteBusy.begin(), vecCondensinSiteBusy.end(), bId2) < DC_VALENCE)  // dissociation & avoid more than 4 contactS
								//if ( L < LTHR_TENSION_CONDENSIN )
							{
								vecCondensinSiteBusy.push_back(bId1);
								vecCondensinSiteBusy.push_back(bId2);

															//cout << "forming condensin association!" << endl;
								dVec veloc1 = vecNode[bId1].get_veloc(), veloc2 = vecNode[bId2].get_veloc(), veloc1_new, veloc2_new;
								double L0 = L0_TENSION_CONDENSIN, L_inv;
								double k = K_TENSION_CONDENSIN;
								dVec unitVec12 = {0,0,0};
								if (L > 0)
								{
									L_inv = 1./L;  // 1/nm
									unitVec12 = {(coord2.x-coord1.x)*L_inv, (coord2.y-coord1.y)*L_inv, (coord2.z-coord1.z)*L_inv};
								}

								if ( abs(L - L0) > L_DIFF_RESOLUTION) // non-zero tension
								{
									double v_ten_con_size = k*(L-L0)*gamma_inv*1E-3;
									dVec v_ten_con = {v_ten_con_size*unitVec12.x, v_ten_con_size*unitVec12.y, v_ten_con_size*unitVec12.z};
									veloc1_new = {veloc1.x + v_ten_con.x, veloc1.y + v_ten_con.y, veloc1.z + v_ten_con.z};
									veloc2_new = {veloc2.x - v_ten_con.x, veloc2.y - v_ten_con.y, veloc2.z - v_ten_con.z};
															// update acceleration
									vecNode[bId1].set_veloc(veloc1_new);
									vecNode[bId2].set_veloc(veloc2_new);
															// update potential & total energy
															//double ePot12 = 0.5*k*(L-L0)*(L-L0);
															//energy["potential"] += ePot12;
															//energy["total"]     += ePot12;
								}
								capturer.push_back( {bId1, bId2});
								
							}
						}
	//////////////////////////////////////////////////
						
						
						if (DIFFUSION_CAPTURE == "multivalent")
						{
							
							//if (   L < LTHR_TENSION_CONDENSIN
								//&& find(vecCondensinSiteBusy.begin(), vecCondensinSiteBusy.end(), bId1) == vecCondensinSiteBusy.end()
								//&& find(vecCondensinSiteBusy.begin(), vecCondensinSiteBusy.end(), bId2) == vecCondensinSiteBusy.end())  // dissociation & avoid more than one contact
							if ( L < LTHR_TENSION_CONDENSIN )
							{
								vecCondensinSiteBusy.push_back(bId1);
								vecCondensinSiteBusy.push_back(bId2);

															//cout << "forming condensin association!" << endl;
								dVec veloc1 = vecNode[bId1].get_veloc(), veloc2 = vecNode[bId2].get_veloc(), veloc1_new, veloc2_new;
								double L0 = L0_TENSION_CONDENSIN, L_inv;
								double k = K_TENSION_CONDENSIN;
								dVec unitVec12 = {0,0,0};
								if (L > 0)
								{
									L_inv = 1./L;  // 1/nm
									unitVec12 = {(coord2.x-coord1.x)*L_inv, (coord2.y-coord1.y)*L_inv, (coord2.z-coord1.z)*L_inv};
								}

								if ( abs(L - L0) > L_DIFF_RESOLUTION) // non-zero tension
								{
									double v_ten_con_size = k*(L-L0)*gamma_inv*1E-3;
									dVec v_ten_con = {v_ten_con_size*unitVec12.x, v_ten_con_size*unitVec12.y, v_ten_con_size*unitVec12.z};
									veloc1_new = {veloc1.x + v_ten_con.x, veloc1.y + v_ten_con.y, veloc1.z + v_ten_con.z};
									veloc2_new = {veloc2.x - v_ten_con.x, veloc2.y - v_ten_con.y, veloc2.z - v_ten_con.z};
															// update acceleration
									vecNode[bId1].set_veloc(veloc1_new);
									vecNode[bId2].set_veloc(veloc2_new);
															// update potential & total energy
															//double ePot12 = 0.5*k*(L-L0)*(L-L0);
															//energy["potential"] += ePot12;
															//energy["total"]     += ePot12;
								}
								capturer.push_back( {bId1, bId2});
								//cout << "multivalent difussion capture" << endl;
								
							}
						}					
						
					}
				}
			}
	}		
	////
	
	if (DIFFUSION_CAPTURE_ON == true && flagIsPrep == false && SWITCH_INTERPHASE_MITOSIS == 0 && DC_MECHANISM == "CSB-BAll")	
	{
			cout << "CSB-BAll" << endl;
			vector<int> vecCondensinSiteBusy; // contain bead ids already interacting with another
			capturer.clear();
			hit_valence.clear();

			vector<int> vecAllFBinderBeadIds;
			vector<int> vecAllRBinderBeadIds;
			vector<int> vecAllBinderBeadIds;
			for (vector<Binder>::iterator it = vecBinder.begin(); it != vecBinder.end(); it++) 
			{
				vector<int> vecBeadIds = it->get_vecBeadIds();
				int bIdF = vecBeadIds[0];
				vecAllFBinderBeadIds.push_back(vecBeadIds[0]);
				int bIdR = vecBeadIds[1];
				vecAllRBinderBeadIds.push_back(vecBeadIds[1]);
				vecAllBinderBeadIds.push_back(vecBeadIds[0]);
				vecAllBinderBeadIds.push_back(vecBeadIds[1]);					
			}
			
			for (int i = 0; i < vecCondensinSite.size(); i++)
			{
				for (int j = 0; j < vecAllBinderBeadIds.size(); j++) 
				{
					if ((rand() / double(RAND_MAX)) > P_DISSOCIATE_CONDENSIN) 
					{
						int bId1 = vecCondensinSite[i], bId2 = vecAllBinderBeadIds[j];
						dVec coord1 = vecNode[bId1].get_coord(), coord2 = vecNode[bId2].get_coord();
						double L = dVecDist(coord1, coord2);

						if (DIFFUSION_CAPTURE == "probabilistic")
						{
							if ((rand() / double(RAND_MAX)) > P_DC_PROBAB)
							{ 
					// multivalent diffusion capture
								if (L < LTHR_TENSION_CONDENSIN)
								{
									vecCondensinSiteBusy.push_back(bId1);
									vecCondensinSiteBusy.push_back(bId2);

									//cout << "forming condensin association!" << endl;
									dVec veloc1 = vecNode[bId1].get_veloc(), veloc2 = vecNode[bId2].get_veloc(), veloc1_new, veloc2_new;
									double L0 = L0_TENSION_CONDENSIN, L_inv;
									double k = K_TENSION_CONDENSIN;
									dVec unitVec12 = { 0, 0, 0 };
									if (L > 0)
									{
										L_inv = 1. / L;  // 1/nm
										unitVec12 = {(coord2.x-coord1.x)*L_inv, (coord2.y-coord1.y)*L_inv, (coord2.z-coord1.z)*L_inv};
									}

									if (abs(L - L0) > L_DIFF_RESOLUTION) // non-zero tension
									{
										double v_ten_con_size = k * (L - L0)* gamma_inv * 1E-3;
										dVec v_ten_con = { v_ten_con_size* unitVec12.x, v_ten_con_size* unitVec12.y, v_ten_con_size* unitVec12.z };
										veloc1_new = {veloc1.x + v_ten_con.x, veloc1.y + v_ten_con.y, veloc1.z + v_ten_con.z};
										veloc2_new = {veloc2.x - v_ten_con.x, veloc2.y - v_ten_con.y, veloc2.z - v_ten_con.z};
										// update acceleration
										vecNode[bId1].set_veloc(veloc1_new);
										vecNode[bId2].set_veloc(veloc2_new);
										// update potential & total energy
										//double ePot12 = 0.5*k*(L-L0)*(L-L0);
										//energy["potential"] += ePot12;
										//energy["total"]     += ePot12;
									}
									capturer.push_back( {bId1, bId2});
								}
								
							}
														
					// bivalent diffusion capture
							else 
							{
								if (L < LTHR_TENSION_CONDENSIN && find(vecCondensinSiteBusy.begin(), vecCondensinSiteBusy.end(), bId1) == vecCondensinSiteBusy.end() && find(vecCondensinSiteBusy.begin(), vecCondensinSiteBusy.end(), bId2) == vecCondensinSiteBusy.end())
								//&& find(vecCondensinSiteBusy.begin(),vecCondensinSiteBusy.end(), bId1) == vecCondensinSiteBusy.end()
								//&& find(vecCondensinSiteBusy.begin(),vecCondensinSiteBusy.end(), bId2) == vecCondensinSiteBusy.end()) // dissociation & avoid more than one contact
								{
									vecCondensinSiteBusy.push_back(bId1);
									vecCondensinSiteBusy.push_back(bId2);

									//cout << "forming condensin association!" << endl;
									dVec veloc1 = vecNode[bId1].get_veloc(), veloc2 = vecNode[bId2].get_veloc(),veloc1_new, veloc2_new;
									double L0 = L0_TENSION_CONDENSIN, L_inv;
									double k = K_TENSION_CONDENSIN;
									dVec unitVec12 = { 0, 0, 0 };
									if (L > 0)
									{
										L_inv = 1. / L;  // 1/nm
										unitVec12 = {(coord2.x-coord1.x)*L_inv, (coord2.y-coord1.y)*L_inv, (coord2.z-coord1.z)*L_inv};
									}

									if ( abs(L - L0) > L_DIFF_RESOLUTION) // non-zero tension
									{
										double v_ten_con_size = k*(L-L0)*gamma_inv*1E-3;
										dVec v_ten_con = {v_ten_con_size*unitVec12.x, v_ten_con_size*unitVec12.y, v_ten_con_size*unitVec12.z};
										veloc1_new = {veloc1.x + v_ten_con.x, veloc1.y + v_ten_con.y, veloc1.z + v_ten_con.z};
										veloc2_new = {veloc2.x - v_ten_con.x, veloc2.y - v_ten_con.y, veloc2.z - v_ten_con.z};
										// update acceleration
										vecNode[bId1].set_veloc(veloc1_new);
										vecNode[bId2].set_veloc(veloc2_new);
										// update potential & total energy
										//double ePot12 = 0.5*k*(L-L0)*(L-L0);
										//energy["potential"] += ePot12;
										//energy["total"]     += ePot12;
									}
								}	capturer.push_back( {bId1, bId2});
								
							}
							
						}

	//					if (DIFFUSION_CAPTURE_PROBAB == true)

						if (DIFFUSION_CAPTURE == "bivalent")
						{
							
							if (   L < LTHR_TENSION_CONDENSIN
								&& find(vecCondensinSiteBusy.begin(), vecCondensinSiteBusy.end(), bId1) == vecCondensinSiteBusy.end()
								&& find(vecCondensinSiteBusy.begin(), vecCondensinSiteBusy.end(), bId2) == vecCondensinSiteBusy.end())  // dissociation & avoid more than one contact
								//if ( L < LTHR_TENSION_CONDENSIN )
							{
								vecCondensinSiteBusy.push_back(bId1);
								vecCondensinSiteBusy.push_back(bId2);

															//cout << "forming condensin association!" << endl;
								dVec veloc1 = vecNode[bId1].get_veloc(), veloc2 = vecNode[bId2].get_veloc(), veloc1_new, veloc2_new;
								double L0 = L0_TENSION_CONDENSIN, L_inv;
								double k = K_TENSION_CONDENSIN;
								dVec unitVec12 = {0,0,0};
								if (L > 0)
								{
									L_inv = 1./L;  // 1/nm
									unitVec12 = {(coord2.x-coord1.x)*L_inv, (coord2.y-coord1.y)*L_inv, (coord2.z-coord1.z)*L_inv};
								}

								if ( abs(L - L0) > L_DIFF_RESOLUTION) // non-zero tension
								{
									double v_ten_con_size = k*(L-L0)*gamma_inv*1E-3;
									dVec v_ten_con = {v_ten_con_size*unitVec12.x, v_ten_con_size*unitVec12.y, v_ten_con_size*unitVec12.z};
									veloc1_new = {veloc1.x + v_ten_con.x, veloc1.y + v_ten_con.y, veloc1.z + v_ten_con.z};
									veloc2_new = {veloc2.x - v_ten_con.x, veloc2.y - v_ten_con.y, veloc2.z - v_ten_con.z};
															// update acceleration
									vecNode[bId1].set_veloc(veloc1_new);
									vecNode[bId2].set_veloc(veloc2_new);
															// update potential & total energy
															//double ePot12 = 0.5*k*(L-L0)*(L-L0);
															//energy["potential"] += ePot12;
															//energy["total"]     += ePot12;
								}
								capturer.push_back( {bId1, bId2});
								
							}
						}
	/////////////////////////////////////////////// 30-10-2018 tetravalent //////////////////////					
						
						if (DIFFUSION_CAPTURE == "n-valent")
						{
							
							if (   L < LTHR_TENSION_CONDENSIN
								&& count(vecCondensinSiteBusy.begin(), vecCondensinSiteBusy.end(), bId1) < DC_VALENCE
								&& count(vecCondensinSiteBusy.begin(), vecCondensinSiteBusy.end(), bId2) < DC_VALENCE)  // dissociation & avoid more than 4 contactS
								//if ( L < LTHR_TENSION_CONDENSIN )
							{
								vecCondensinSiteBusy.push_back(bId1);
								vecCondensinSiteBusy.push_back(bId2);

															//cout << "forming condensin association!" << endl;
								dVec veloc1 = vecNode[bId1].get_veloc(), veloc2 = vecNode[bId2].get_veloc(), veloc1_new, veloc2_new;
								double L0 = L0_TENSION_CONDENSIN, L_inv;
								double k = K_TENSION_CONDENSIN;
								dVec unitVec12 = {0,0,0};
								if (L > 0)
								{
									L_inv = 1./L;  // 1/nm
									unitVec12 = {(coord2.x-coord1.x)*L_inv, (coord2.y-coord1.y)*L_inv, (coord2.z-coord1.z)*L_inv};
								}

								if ( abs(L - L0) > L_DIFF_RESOLUTION) // non-zero tension
								{
									double v_ten_con_size = k*(L-L0)*gamma_inv*1E-3;
									dVec v_ten_con = {v_ten_con_size*unitVec12.x, v_ten_con_size*unitVec12.y, v_ten_con_size*unitVec12.z};
									veloc1_new = {veloc1.x + v_ten_con.x, veloc1.y + v_ten_con.y, veloc1.z + v_ten_con.z};
									veloc2_new = {veloc2.x - v_ten_con.x, veloc2.y - v_ten_con.y, veloc2.z - v_ten_con.z};
															// update acceleration
									vecNode[bId1].set_veloc(veloc1_new);
									vecNode[bId2].set_veloc(veloc2_new);
															// update potential & total energy
															//double ePot12 = 0.5*k*(L-L0)*(L-L0);
															//energy["potential"] += ePot12;
															//energy["total"]     += ePot12;
								}
								capturer.push_back( {bId1, bId2});
								
							}
						}
	//////////////////////////////////////////////////
						
						
						if (DIFFUSION_CAPTURE == "multivalent")
						{
							
							//if (   L < LTHR_TENSION_CONDENSIN
								//&& find(vecCondensinSiteBusy.begin(), vecCondensinSiteBusy.end(), bId1) == vecCondensinSiteBusy.end()
								//&& find(vecCondensinSiteBusy.begin(), vecCondensinSiteBusy.end(), bId2) == vecCondensinSiteBusy.end())  // dissociation & avoid more than one contact
							if ( L < LTHR_TENSION_CONDENSIN )
							{
								vecCondensinSiteBusy.push_back(bId1);
								vecCondensinSiteBusy.push_back(bId2);

															//cout << "forming condensin association!" << endl;
								dVec veloc1 = vecNode[bId1].get_veloc(), veloc2 = vecNode[bId2].get_veloc(), veloc1_new, veloc2_new;
								double L0 = L0_TENSION_CONDENSIN, L_inv;
								double k = K_TENSION_CONDENSIN;
								dVec unitVec12 = {0,0,0};
								if (L > 0)
								{
									L_inv = 1./L;  // 1/nm
									unitVec12 = {(coord2.x-coord1.x)*L_inv, (coord2.y-coord1.y)*L_inv, (coord2.z-coord1.z)*L_inv};
								}

								if ( abs(L - L0) > L_DIFF_RESOLUTION) // non-zero tension
								{
									double v_ten_con_size = k*(L-L0)*gamma_inv*1E-3;
									dVec v_ten_con = {v_ten_con_size*unitVec12.x, v_ten_con_size*unitVec12.y, v_ten_con_size*unitVec12.z};
									veloc1_new = {veloc1.x + v_ten_con.x, veloc1.y + v_ten_con.y, veloc1.z + v_ten_con.z};
									veloc2_new = {veloc2.x - v_ten_con.x, veloc2.y - v_ten_con.y, veloc2.z - v_ten_con.z};
															// update acceleration
									vecNode[bId1].set_veloc(veloc1_new);
									vecNode[bId2].set_veloc(veloc2_new);
															// update potential & total energy
															//double ePot12 = 0.5*k*(L-L0)*(L-L0);
															//energy["potential"] += ePot12;
															//energy["total"]     += ePot12;
								}
								capturer.push_back( {bId1, bId2});
								//cout << "multivalent difussion capture" << endl;
								
							}
						}					
						
					}
				}
			}
	}		
	
	
	
	
	
	
	
	
	
	
	
		/*
		 for (int i=0; i < capturer.size(); i++)
		 {
		 iPair CAP0 = capturer[i];
		 int id1 = CAP0.i, id2 = CAP0.j
		 int count_id1 = count( capturer.begin(), capturer.end(), id1);
		 int count_id2 = count( capturer.begin(), capturer.end(), id2);
		 hit_valence.push_back(count_id1);
		 hit_valence.push_back(count_id2);

		 }
		 */


    // (4) loop extrusion   -------   sliding binders
    if (true)
    {
        vector<int> vecAllAttachBeadIds;    // include all AttachBeadIds for condition for stopping binders
        for (vector<Binder>::iterator it = vecBinder.begin(); it != vecBinder.end(); it ++)
        {
            vector<int> vecBeadIds = it->get_vecBeadIds();
            vector<int> vecAttachBeadIds = it->get_vecAttachBeadIds();
            vecAllAttachBeadIds.push_back(vecAttachBeadIds[1]); // rear
            vecAllAttachBeadIds.push_back(vecAttachBeadIds[0]); // front
            int bIdF = vecBeadIds[0], bIdR = vecBeadIds[1], bIdFA = vecAttachBeadIds[0], bIdRA = vecAttachBeadIds[1];
            dVec coordF  = vecNode[bIdF].get_coord(),  coordR  = vecNode[bIdR].get_coord();
            dVec coordFA = vecNode[bIdFA].get_coord(), coordRA = vecNode[bIdRA].get_coord();
            dVec velocF  = vecNode[bIdF].get_veloc(),  velocR  = vecNode[bIdR].get_veloc();
            dVec velocFA = vecNode[bIdFA].get_veloc(), velocRA = vecNode[bIdRA].get_veloc();
            dVec velocF_new = velocF, velocR_new = velocR, velocFA_new = velocFA, velocRA_new = velocRA;

            if (true)
            {
                // (4.1) between front and rear feet beads
                double L0 = L0_TENSION_BINDER, L = dVecDist(coordF, coordR);   // nm
                double k = K_TENSION_BINDER;                                   // pN/nm     (N = kg*m/s/s; nN = kg*nm/s/s; pN = kg*pm/s/s)
                dVec unitVecRF = {0,0,0};
                if (L > 0)
                {
                    double L_inv = 1./L;  // 1/nm
                    unitVecRF = {(coordF.x-coordR.x)*L_inv, (coordF.y-coordR.y)*L_inv, (coordF.z-coordR.z)*L_inv};
                }
                double v_ten_size = k*(L-L0)*gamma_inv*1E-3;   // nm/s
                dVec v_ten = {v_ten_size*unitVecRF.x, v_ten_size*unitVecRF.y, v_ten_size*unitVecRF.z};
                dVec velocR_tmp1 = velocR_new, velocF_tmp1 = velocF_new;
                velocR_new = {  velocR_tmp1.x + v_ten.x, velocR_tmp1.y + v_ten.y, velocR_tmp1.z + v_ten.z };
                velocF_new = {  velocF_tmp1.x - v_ten.x, velocF_tmp1.y - v_ten.y, velocF_tmp1.z - v_ten.z };

                // (4.2) between rear feet bead and its attachment bead
                L0 = L0_TENSION_BINDER_BEAD;
                L = dVecDist(coordRA, coordR);   // nm
                k = K_TENSION_BINDER_BEAD;
                dVec unitVecRRA = {0,0,0};
                if (L > 0)
                {
                    double L_inv = 1./L;  // 1/nm
                    unitVecRRA = {(coordRA.x-coordR.x)*L_inv, (coordRA.y-coordR.y)*L_inv, (coordRA.z-coordR.z)*L_inv};
                }
                v_ten_size = k*(L-L0)*gamma_inv*1E-3;   // nm/s
                v_ten = {v_ten_size*unitVecRRA.x, v_ten_size*unitVecRRA.y, v_ten_size*unitVecRRA.z};
                dVec velocR_tmp2 = velocR_new, velocRA_tmp2 = velocRA_new;
                velocR_new  = {  velocR_tmp2.x  + v_ten.x, velocR_tmp2.y  + v_ten.y, velocR_tmp2.z  + v_ten.z };
                velocRA_new = {  velocRA_tmp2.x - v_ten.x, velocRA_tmp2.y - v_ten.y, velocRA_tmp2.z - v_ten.z };
                vecNode[bIdR].set_veloc(velocR_new);
                vecNode[bIdRA].set_veloc(velocRA_new);

                // (4.3) between front feer bead and its attachment bead
                L = dVecDist(coordFA, coordF);   // nm
                dVec unitVecFFA = {0,0,0};
                if (L > 0)
                {
                    double L_inv = 1./L;  // 1/nm
                    unitVecFFA = {(coordFA.x-coordF.x)*L_inv, (coordFA.y-coordF.y)*L_inv, (coordFA.z-coordF.z)*L_inv};
                }
                v_ten_size = k*(L-L0)*gamma_inv*1E-3;   // nm/s
                v_ten = {v_ten_size*unitVecFFA.x, v_ten_size*unitVecFFA.y, v_ten_size*unitVecFFA.z};
                dVec velocF_tmp2 = velocF_new, velocFA_tmp2 = velocFA_new;
                velocF_new  = {  velocF_tmp2.x  + v_ten.x, velocF_tmp2.y  + v_ten.y, velocF_tmp2.z  + v_ten.z };
                velocFA_new = {  velocFA_tmp2.x - v_ten.x, velocFA_tmp2.y - v_ten.y, velocFA_tmp2.z - v_ten.z };
                vecNode[bIdF].set_veloc(velocF_new);
                vecNode[bIdFA].set_veloc(velocFA_new);
            }

        }

        // also consider static road blocks in vecBinderBlockSite
        for (vector<int>::iterator it = vecBinderBlockSite.begin(); it != vecBinderBlockSite.end(); it ++)
            vecAllAttachBeadIds.push_back(*it);

        // print vecAttachBeadIds
        if (false)
        {
            cout << "--- Current vecAttachBeadIds ---\n";
            for (vector<int>::iterator it = vecAllAttachBeadIds.begin(); it != vecAllAttachBeadIds.end(); it ++)
            {
                cout << *it << "\t";
                if ( (it-vecAllAttachBeadIds.begin()) % 2 == 0  )
                 cout << "----\t";
            }
            cout << endl;
        }

        // (4.4) probabilistic walking -- simple implementation
        if (LOOP_EXTRUSION_ON == true && flagIsPrep == false && flagLE_ON == true)
        {
            for (vector<Binder>::iterator it = vecBinder.begin(); it != vecBinder.end(); it ++)
            {

                vector<bool> canSlide = it->get_canSlide(), canSlide_new = canSlide;
                vector<int> vecAttachBeadIds = it->get_vecAttachBeadIds();
                int binderId = it->get_id();
                int bIdFA = vecAttachBeadIds[0], bIdRA = vecAttachBeadIds[1];
                int bIdFA_new = bIdFA, bIdRA_new = bIdRA;
                
                
                // -------------------- sym ------------------------------
                if ( LOOP_EXTRUSION == "probabilistic") 
                {
					if (rand() / (double) RAND_MAX <= P_LE_PROBAB)
					{
						// front
						if (canSlide_new[0] == true)
						{
							if (find(vecAllAttachBeadIds.begin(), vecAllAttachBeadIds.end(), (bIdFA+1)) != vecAllAttachBeadIds.end()) // stop
							{
								cout << "Binder id = " << binderId << " FRONT stops translocation at Bead id = " << bIdFA << endl;
								canSlide_new[0] = false;
							}
							else if (rand() / (double) RAND_MAX <= P_BINDER_SLIDING && bIdFA < NUM_BEAD-1)                            // walk
								bIdFA_new += 1;
						}
						
						
						// rear
						if (canSlide_new[1] == true)
						{
							if ((find(vecAllAttachBeadIds.begin(), vecAllAttachBeadIds.end(), (bIdRA-1)) != vecAllAttachBeadIds.end())
							||  (find(vecAllAttachBeadIds.begin(), vecAllAttachBeadIds.end(), (bIdRA-2)) != vecAllAttachBeadIds.end())) // stop
							{
								cout << "Binder id = " << binderId << "  REAR stops translocation at Bead id = " << bIdRA << endl;
								canSlide_new[1] = false;
							}
							else if (rand() / (double) RAND_MAX <= P_BINDER_SLIDING && bIdRA > 0)                                     // walk
								bIdRA_new -= 1;
						}

						it->set_vecAttachBeadIds({bIdFA_new, bIdRA_new});
						it->set_canSlide(canSlide_new);
					}
					
					else 
					{
						// rear
						if (canSlide_new[1] == true)
						{
							if ((find(vecAllAttachBeadIds.begin(), vecAllAttachBeadIds.end(), (bIdRA-1)) != vecAllAttachBeadIds.end())
							||  (find(vecAllAttachBeadIds.begin(), vecAllAttachBeadIds.end(), (bIdRA-2)) != vecAllAttachBeadIds.end())) // stop
							{
								cout << "Binder id = " << binderId << "  REAR stops translocation at Bead id = " << bIdRA << endl;
								canSlide_new[1] = false;
							}
							else if (rand() / (double) RAND_MAX <= P_BINDER_SLIDING && bIdRA > 0)                                     // walk
								bIdRA_new -= 1;
						}


						it->set_vecAttachBeadIds({bIdFA_new, bIdRA_new});
						it->set_canSlide(canSlide_new);
					}	
				}
				
				
				
				
				// -------------------- sym ------------------------------

				if ( LOOP_EXTRUSION == "symmetric")
				{
					
				// front
					if (canSlide_new[0] == true)
					{
						if (find(vecAllAttachBeadIds.begin(), vecAllAttachBeadIds.end(), (bIdFA+1)) != vecAllAttachBeadIds.end()) // stop
						{
							cout << "Binder id = " << binderId << " FRONT stops translocation at Bead id = " << bIdFA << endl;
							canSlide_new[0] = false;
						}
						else if (rand() / (double) RAND_MAX <= P_BINDER_SLIDING && bIdFA < NUM_BEAD-1)                            // walk
							bIdFA_new += 1;
					}
				// rear
					if (canSlide_new[1] == true)
					{
						if ((find(vecAllAttachBeadIds.begin(), vecAllAttachBeadIds.end(), (bIdRA-1)) != vecAllAttachBeadIds.end())
						||  (find(vecAllAttachBeadIds.begin(), vecAllAttachBeadIds.end(), (bIdRA-2)) != vecAllAttachBeadIds.end())) // stop
						{
							cout << "Binder id = " << binderId << "  REAR stops translocation at Bead id = " << bIdRA << endl;
							canSlide_new[1] = false;
						}
						else if (rand() / (double) RAND_MAX <= P_BINDER_SLIDING && bIdRA > 0)                                     // walk
							bIdRA_new -= 1;
					}


					it->set_vecAttachBeadIds({bIdFA_new, bIdRA_new});
					it->set_canSlide(canSlide_new);
						
				}
				
				
				
				// -------------------- R asym ------------------------------
				if ( LOOP_EXTRUSION == "asymmetricR")
				{
				// rear
					if (canSlide_new[1] == true)
					{
						if ((find(vecAllAttachBeadIds.begin(), vecAllAttachBeadIds.end(), (bIdRA-1)) != vecAllAttachBeadIds.end())
						||  (find(vecAllAttachBeadIds.begin(), vecAllAttachBeadIds.end(), (bIdRA-2)) != vecAllAttachBeadIds.end())) // stop
						{
							cout << "Binder id = " << binderId << "  REAR stops translocation at Bead id = " << bIdRA << endl;
							canSlide_new[1] = false;
						}
						else if (rand() / (double) RAND_MAX <= P_BINDER_SLIDING && bIdRA > 0)                                     // walk
							bIdRA_new -= 1;
					}


					it->set_vecAttachBeadIds({bIdFA_new, bIdRA_new});
					it->set_canSlide(canSlide_new);
				}
				
				
				
				
				// -------------------- F sym ------------------------------
				if ( LOOP_EXTRUSION == "asymmetricF")
				{
				// front
					if (canSlide_new[0] == true)
					{
						if (find(vecAllAttachBeadIds.begin(), vecAllAttachBeadIds.end(), (bIdFA+1)) != vecAllAttachBeadIds.end()) // stop
						{
							cout << "Binder id = " << binderId << " FRONT stops translocation at Bead id = " << bIdFA << endl;
							canSlide_new[0] = false;
						}
						else if (rand() / (double) RAND_MAX <= P_BINDER_SLIDING && bIdFA < NUM_BEAD-1)                            // walk
							bIdFA_new += 1;
					}

					it->set_vecAttachBeadIds({bIdFA_new, bIdRA_new});
					it->set_canSlide(canSlide_new);
				}
				
				
				
				
				
				// -------------------- probability for F vs R asym ------------------------------
				if (LOOP_EXTRUSION == "asymprobabilistic")
				{
					if (rand() / (double) RAND_MAX <= P_LE_PROBAB)
					{
						// front

						if (canSlide_new[0] == true)
						{
							if (find(vecAllAttachBeadIds.begin(), vecAllAttachBeadIds.end(), (bIdFA+1)) != vecAllAttachBeadIds.end()) // stop
							{
								cout << "Binder id = " << binderId << " FRONT stops translocation at Bead id = " << bIdFA << endl;
								canSlide_new[0] = false;
							}
							else if (rand() / (double) RAND_MAX <= P_BINDER_SLIDING && bIdFA < NUM_BEAD-1)                            // walk
								bIdFA_new += 1;
						}

						it->set_vecAttachBeadIds({bIdFA_new, bIdRA_new});
						it->set_canSlide(canSlide_new);
					}
					
					else 
					{
						// rear
						if (canSlide_new[1] == true)
						{
							if ((find(vecAllAttachBeadIds.begin(), vecAllAttachBeadIds.end(), (bIdRA-1)) != vecAllAttachBeadIds.end())
							||  (find(vecAllAttachBeadIds.begin(), vecAllAttachBeadIds.end(), (bIdRA-2)) != vecAllAttachBeadIds.end())) // stop
							{
								cout << "Binder id = " << binderId << "  REAR stops translocation at Bead id = " << bIdRA << endl;
								canSlide_new[1] = false;
							}
							else if (rand() / (double) RAND_MAX <= P_BINDER_SLIDING && bIdRA > 0)                                     // walk
								bIdRA_new -= 1;
						}


						it->set_vecAttachBeadIds({bIdFA_new, bIdRA_new});
						it->set_canSlide(canSlide_new);
					}
				}                   

                // (3*) condensin attraction (diffusion capture of binders)
                if (false)
                {
                    it->set_canSlide({true, true}); // reset to true
                    vector<int> vecBeadIds = it->get_vecBeadIds();
                    int bIdF = vecBeadIds[0], bIdR = vecBeadIds[1];
                    dVec coordF = vecNode[bIdF].get_coord(), coordR = vecNode[bIdR].get_coord();
                    dVec coord = {(coordF.x+coordR.x)*0.5, (coordF.y+coordR.y)*0.5, (coordF.z+coordR.z)*0.5};
                    for (vector<Binder>::iterator it2 = it; it2 != vecBinder.end(); it2 ++)
                    {
                        vector<bool> canSlide2 = it2->get_canSlide(), canSlide2_new = canSlide2;
                        if ( (rand() / double(RAND_MAX)) > P_DISSOCIATE_CONDENSIN )
                        {
                            vector<int> vecBeadIds2 = it2->get_vecBeadIds();
                            int bIdF2 = vecBeadIds2[0], bIdR2 = vecBeadIds2[1];
                            dVec coordF2 = vecNode[bIdF2].get_coord(), coordR2 = vecNode[bIdR2].get_coord();
                            dVec coord2 = {(coordF2.x+coordR2.x)*0.5, (coordF2.y+coordR2.y)*0.5, (coordF2.z+coordR2.z)*0.5};
                            double L = dVecDist(coord, coord2);
                            if (   L < LTHR_TENSION_CONDENSIN)  //
                            {
                                dVec velocF  = vecNode[bIdF ].get_veloc(), velocR  = vecNode[bIdR ].get_veloc(), velocF_new,  velocR_new;
                                dVec velocF2 = vecNode[bIdF2].get_veloc(), velocR2 = vecNode[bIdR2].get_veloc(), velocF2_new, velocR2_new;
                                double L0 = L0_TENSION_CONDENSIN, L_inv;
                                double k = K_TENSION_CONDENSIN;
                                dVec unitVec12 = {0,0,0};
                                if (L > 0)
                                {
                                    L_inv = 1./L;  // 1/nm
                                    unitVec12 = {(coord2.x-coord.x)*L_inv, (coord2.y-coord.y)*L_inv, (coord2.z-coord.z)*L_inv};
                                }

                                if ( abs(L - L0) > L_DIFF_RESOLUTION) // non-zero tension
                                {
                                    double v_ten_con_size = k*(L-L0)*gamma_inv*1E-3;
                                    dVec v_ten_con = {v_ten_con_size*unitVec12.x, v_ten_con_size*unitVec12.y, v_ten_con_size*unitVec12.z};
                                    velocF_new  = {  velocF.x + v_ten_con.x,  velocF.y + v_ten_con.y,  velocF.z + v_ten_con.z };
                                    velocR_new  = {  velocR.x + v_ten_con.x,  velocR.y + v_ten_con.y,  velocR.z + v_ten_con.z };
                                    velocF2_new = {  velocF2.x - v_ten_con.x, velocF2.y - v_ten_con.y, velocF2.z - v_ten_con.z };
                                    velocR2_new = {  velocR2.x - v_ten_con.x, velocR2.y - v_ten_con.y, velocR2.z - v_ten_con.z };
                                    // update velocity
                                    vecNode[bIdF].set_veloc(velocF_new);
                                    vecNode[bIdR].set_veloc(velocR_new);
                                    vecNode[bIdF2].set_veloc(velocF2_new);
                                    vecNode[bIdR2].set_veloc(velocR2_new);
                                }

                                // update canSlide2 -- Note: if there is maximum number of assocition, check whether canSlide is already FALSE!!!
                                if (rand()/(double) RAND_MAX <= P_BINDER_STOP_UPON_CAPTURE)  // probabilistic stopping
                                {
                                    canSlide_new  = {false, false};
                                    canSlide2_new = {false, false};
                                }
                                else
                                {
                                    canSlide_new  = {true, true};
                                    canSlide2_new = {true, true};
                                }
                                it2->set_canSlide(canSlide2_new);
                            }
                        }
                    }

                    // update canSlide
                    it->set_canSlide(canSlide_new);
                }


            }
        }

    }

}
// -----------------------------------------------------------------

void updateVeloc(vector<Node> & vecNode, unordered_map<string, double> & energy)
{
    for (vector<Node>::iterator it = vecNode.begin(); it != vecNode.end(); it++)
    {
        dVec veloc = it->get_veloc(), accel = it->get_accel(), accelPrev = it->get_accelPrev();
        dVec veloc_new;
        veloc_new = { veloc.x + 0.5*(accel.x+accelPrev.x)*DT,
                      veloc.y + 0.5*(accel.y+accelPrev.y)*DT,
                      veloc.z + 0.5*(accel.z+accelPrev.z)*DT };
        it->set_veloc(veloc_new);

        // summ up kinetic energy
        double speedSq = veloc.x*veloc.x + veloc.y*veloc.y + veloc.z*veloc.z;
        //double speedSq = veloc_new.x*veloc_new.x + veloc_new.y*veloc_new.y + veloc_new.z*veloc_new.z;
        //double eKin = 0.5*MASS_BEAD*speedSq;
        double eKin = 0.5*speedSq;  // kinetic energy per unit mass
        energy["kinetic"] += eKin;
        //cout << "kinetic energy: add = " << eKin << " , sum = " << energy["kinetic"] << endl;
        energy["total"]   += eKin;
    }
}
void updateVelocLinearDamp(vector<Node> & vecNode, unordered_map<string, double> & energy)
{
    // refer to modified verlet method ( equation (34) in cmotion.pdf )
    // STEP 3 : v(t+dt) = ( v(t)*( 1 - 0.5*dt*gamma/m )   0.5* ( a(t)+a(t+dt) ) * dt ) / (1 + 0.5*dt*gamma/m)
    for (vector<Node>::iterator it = vecNode.begin(); it != vecNode.end(); it++)
    {
        dVec veloc = it->get_veloc(), accel = it->get_accel(), accelPrev = it->get_accelPrev();
        dVec veloc_new;
        double GAMMA = GAMMA_VISCOUS;
        double factor = 0.5*DT*GAMMA/MASS_BEAD;
        veloc_new = { ((1-factor)*veloc.x + 0.5*(accel.x+accelPrev.x)*DT)/(1+factor),
                      ((1-factor)*veloc.y + 0.5*(accel.y+accelPrev.y)*DT)/(1+factor),
                      ((1-factor)*veloc.z + 0.5*(accel.z+accelPrev.z)*DT)/(1+factor) };
        it->set_veloc(veloc_new);

        double speedSq = veloc.x*veloc.x + veloc.y*veloc.y + veloc.z*veloc.z;
        //double eKin = 0.5*MASS_BEAD*speedSq;
        double eKin = 0.5*speedSq;
        energy["kinetic"] += eKin;
        energy["total"]   += eKin;
    }
}

double dVecDist(dVec & p1, dVec & p2)
{
    return sqrt( (p1.x-p2.x)*(p1.x-p2.x) + (p1.y-p2.y)*(p1.y-p2.y) + (p1.z-p2.z)*(p1.z-p2.z) );
}
