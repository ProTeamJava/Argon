/*
 * argon.h
 *
 *  Created on: Oct 5, 2016
 *      Author: mongos
 */

#pragma once

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

struct coor
{
    double x;
    double y;
    double z;
};

struct atom
{
    coor r;
    coor p;
    coor E;
    coor F;

};

struct Parameters
{
    double n;
    double m;
    double e;
    double R;
    double f;
    double L;
    double a;
    double T_0;
    double tau;
    double S_0;
    double S_d;
    int S_out;
    int S_xyz;
};

struct State
{
    double V;
    double P;
    double H;
    double T;
    vector<atom> atoms;
};

Parameters getParameters(string inputFileName);

void setInitialState(Parameters &Parameters, State & state);
void setInitialLocation(Parameters &Parameters, vector<atom> &atoms);
void setInitialEnergies(Parameters &Parameters, vector<atom> &atoms);
void setInitialMomentum(Parameters &Parameters, vector<atom> &atoms);

void setPotentialForcesAndPressure(Parameters &Parameters, State & state);

double calcPotentialS(Parameters &Parameters, double ri);
double calcPotentialP(Parameters &Parameters, double rij);
coor calcForcesP(Parameters &Parameters, coor ri, coor rj);
coor calcForcesS(Parameters &Parameters, coor ri);

void simulate(Parameters &Parameters, State &state,
              ofstream &outputFileXYZ, ofstream &outputFileChar);

void updateState(Parameters &Parameters, State & state);
void setEnergyAndTemperature(Parameters &Parameters, State & state);

double getUniRandom();
int getPlusOrMinus();

double calcVectorModulus(coor vector);
coor subtractVectors(coor v1, coor v2);
coor addVectors(coor v1, coor v2);
coor calcOppositeVector(coor v);

void outputXYZ(vector<atom> &atoms, ofstream &outputFile);
void outputChar(State &state, ofstream &outputFile, double &t);
void outputMomentum(vector<atom> &atoms);

