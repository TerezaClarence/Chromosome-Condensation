/*
    File: dynamics.hpp
    Function: Simulate dynamics of chromosomes and nuclear envelope
    Model: chromoCell
    Created: 11 December, 2017 (XF)
*/

#ifndef DYNAMICS_HPP_INCLUDED
#define DYNAMICS_HPP_INCLUDED

#include "initConfig.hpp"

// ===================== Functions =========================
void oneIter(vector<Chromosome> & vecChromosome, vector<Node> & vecNode, vector<iPair> vecBeadPairInteract,
             vector<double> & vecBeadPairInteractFreq, vector<int> & vecFluoroSite,
             unordered_map<string, double> & energy, VoxMap & voxMap, bool FLAG_INTERACT_ON, vector<double> & randNormNum0, bool flagIsPrep);   // if using VoxMap
//void oneIter(vector<Chromosome> & vecChromosome, vector<Node> & vecNode, vector<iPair> vecBeadPairInteract,
//             unordered_map<string, double> & energy, Octree & octree, bool FLAG_INTERACT_ON);   // if using Octree
void oneIterBrownianSimulation(vector<Chromosome> & vecChromosome, vector<Node> & vecNode,
                               vector<int> & vecCondensinSite, vector<int> & vecFluoroSite,
                               unordered_map<string, double> & energy, VoxMap & voxMap,
                               vector<double> & randNormNum0, bool flagIsPrep);
void oneIterBrownianSimulationOverdamp(vector<Chromosome> & vecChromosome, vector<Node> & vecNode,
                                       vector<int> & vecCondensinSite, vector<int> & vecBinderBlockSite, vector<int> & vecFluoroSite,
                                       vector<Binder> & vecBinder,
                                       unordered_map<string, double> & energy, VoxMap & voxMap,
                                       vector<double> & randNormNum0, bool flagIsPrep, bool flagIsInitConfig, vector<iPair> & capturer, vector<int> & hit_valence, map<int, vector<int>> & Domain_CSB, bool flagLE_ON, bool flagIsInterphase);

// STEP 1
void updateCoord(vector<Node> & vecNode, vector<double> & randNormNum0, bool flagTetherEnds, bool flagFreezeFluoro, vector<int> & vecFluoroSite, bool flagIsPrep);
void updateCoord_OMP(vector<Node> & vecNode);   // not optimized
// STEP 2
// --- following functions are for RECONSTRUCTION objective ---
void updateAccel(vector<Chromosome> & vecChromosome, vector<Node> & vecNode,
                 vector<iPair> vecBeadPairInteract, vector<double> & vecBeadPairInteractFreq,
                 VoxMap & voxMap, unordered_map<string, double> & energy, bool FLAG_INTERACT_ON);

// --- following functions are for Brownian SIMULATION objective ---
void updateAccelBrownianSimulation(vector<Chromosome> & vecChromosome, vector<Node> & vecNode, vector<int> & vecCondensinSite,
                 VoxMap & voxMap, unordered_map<string, double> & energy, vector<double> & randNormNum0);
void updateVelocBrownianSimulationOverdamp(vector<Chromosome> & vecChromosome, vector<Node> & vecNode, vector<int> & vecCondensinSite,
                                           vector<int> & vecBinderBlockSite, vector<Binder> & vecBinder,
                                           VoxMap & voxMap, unordered_map<string, double> & energy, vector<double> & randNormNum0, bool flagIsPrep, vector<iPair> & capturer, vector<int> & hit_valence, map<int, vector<int>> & Domain_CSB, bool flagLE_ON, bool flagIsInterphase);

// STEP 3
void updateVeloc(vector<Node> & vecNode, unordered_map<string, double> & energy);
void updateVelocLinearDamp(vector<Node> & vecNode, unordered_map<string, double> & energy);

double dVecDist(dVec & p1, dVec & p2);

// VoxMap algorithm
void updateVoxMap(VoxMap & voxMap, vector<Node> & vecNode);
// Octree algorithm
void updateOctree(Octree & octree, vector<Node> & vecNode, int depth_init = 2);
void branchLeaf(vector<Leaf> & vecLeaf, vector<bVec> position_child0);
bool coord2position(double fracX);

#endif // DYNAMICS_HPP_INCLUDED
