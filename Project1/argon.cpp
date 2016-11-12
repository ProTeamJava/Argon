/*
 * argon.cpp
 *
 *  Created on: Oct 5, 2016
 *      Author: mongos
 */

#include "argon.h"

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstring>

#ifndef _WIN32
int fopen_s(FILE **f, const char *name, const char *mode) {
    int ret = 0;
    assert(f);
    *f = fopen(name, mode);
    if (!*f)
        ret = errno;
    return ret;
}
#define fscanf_s fscanf
#endif

/* physical constants */
namespace {

const double kB = 8.31e-3;

}

Parameters getParameters(std::string inputFileName)
{
    Parameters param = {};
    FILE *inputFile;
    fopen_s(&inputFile, inputFileName.c_str(), "r");

    if(inputFile)
    {
        fscanf_s(inputFile, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d",
                 &param.n,
                 &param.m,
                 &param.e,
                 &param.R,
                 &param.f,
                 &param.L,
                 &param.a,
                 &param.T_0,
                 &param.tau,
                 &param.S_0,
                 &param.S_d,
                 &param.S_out,
                 &param.S_xyz);
    }
    else
    {
        std::cerr << "\nError, unable to open products.txt.\n";
    }


    fclose(inputFile);

    return param;
}

void setInitialState(Parameters &parameters, State &state)
{
    setInitialLocation(parameters, state.atoms);
    setInitialEnergies(parameters, state.atoms);
    setInitialMomentum(parameters, state.atoms);
    setPotentialForcesAndPressure(parameters, state);
}

void setInitialEnergies(Parameters &parameters, vector<atom> & atoms)
{
    for(int i = 0; i < (int) atoms.size(); i++)
    {
        atoms[i].E ={- 0.5 * kB * parameters.T_0 * log(getUniRandom()) ,
                     - 0.5 * kB * parameters.T_0 * log(getUniRandom()),
                     - 0.5 * kB * parameters.T_0 * log(getUniRandom())};
    }
}

void setInitialLocation(Parameters &parameters, vector<atom> &atoms)
{
    coor b[3] = {{parameters.a, 0, 0},
                 {parameters.a / 2.0, parameters.a * sqrt(3.0) / 2.0, 0},
                 {parameters.a / 2.0, parameters.a * sqrt(3.0) / 6.0, parameters.a * sqrt(2.0 / 3.0)}};

    double N = pow(parameters.n, 3);

    atoms.resize((int)N);
    int i = 0;

    for(int i0 = 0; i0 < parameters.n; i0++)
        for(int i1 = 0; i1 < parameters.n; i1++)
            for(int i2 = 0; i2 < parameters.n; i2++, i++)
            {
                atoms[i].r = {(i0 - (parameters.n-1)/2.) * b[0].x
                              + (i1 - (parameters.n-1)/2.) * b[1].x
                              + (i2 - (parameters.n-1)/2.) * b[2].x,
                              (i0 - (parameters.n-1)/2.) * b[0].y
                              + (i1 - (parameters.n-1)/2.) * b[1].y
                              + (i2 - (parameters.n-1)/2.) * b[2].y,
                              (i0 - (parameters.n-1)/2.) * b[0].z
                              + (i1 - (parameters.n-1)/2.) * b[1].z
                              + (i2 - (parameters.n-1)/2.) * b[2].z};
            }

}

void setInitialMomentum(Parameters &parameters, vector<atom> &atoms)
{
    coor Psum = {0, 0, 0};

    for(int i = 0; i < (int) atoms.size(); i++)
    {
        atoms[i].p = {
            getPlusOrMinus() * sqrt(2.*parameters.m*atoms[i].E.x),
            getPlusOrMinus() * sqrt(2.*parameters.m*atoms[i].E.y),
            getPlusOrMinus() * sqrt(2.*parameters.m*atoms[i].E.z)};

        Psum.x += atoms[i].p.x;
        Psum.y += atoms[i].p.y;
        Psum.z += atoms[i].p.z;
    }

    for(int i = 0; i < (int) atoms.size(); i++)
    {
        atoms[i].p = {
            atoms[i].p.x - Psum.x / (double)atoms.size(),
            atoms[i].p.y - Psum.y / (double)atoms.size(),
            atoms[i].p.z - Psum.z / (double)atoms.size()};
    }
}


void setPotentialForcesAndPressure(Parameters &parameters, State &state)
{
    double N = pow(parameters.n, 3);

    state.V = 0;
    state.T = 0;
    state.P = 0;

    for(int i = 0; i < N; i++)
    {
        state.atoms[i].F = {0,0,0};
    }

    for(int i = 0; i < N; i++)
    {
        state.V += calcPotentialS(parameters, calcVectorModulus(state.atoms[i].r)); //potencjal od scianek (10)

        coor Fis = calcForcesS(parameters, state.atoms[i].r); //sily odpychania od scianek (14)
        state.atoms[i].F = addVectors(state.atoms[i].F, Fis ); //akumulacja do F
        state.P += calcVectorModulus(Fis) / 4. / M_PI / parameters.L / parameters.L; //akumulacja cisnienia chwilowego (15)

        if (i > 0)
        {
            for(int j = 0; j < i; j++)
            {
                double rij = calcVectorModulus( subtractVectors(state.atoms[i].r, state.atoms[j].r));
                state.V += calcPotentialP(parameters, rij); //obliczanie potencjalu par (9) i akumulacja do V
                coor Fip = calcForcesP(parameters, state.atoms[i].r, state.atoms[j].r); //obliczanie sil miedzyatomowych
                state.atoms[i].F = addVectors(state.atoms[i].F, Fip); //  akumulacja do Fi
                state.atoms[j].F = addVectors(state.atoms[j].F, calcOppositeVector(Fip)); // akuulacja do Fj
            }
        }

    }
}

double calcPotentialS(Parameters &parameters, double ri)
{
    if(ri <= parameters.L)
        return 0;

    return (0.5 * parameters.f * (ri - parameters.L) * (ri - parameters.L));
}

double calcPotentialP(Parameters &parameters, double rij)
{
    double rRatio = parameters.R / rij;

    return parameters.e * (pow( rRatio , 12 ) - 2 * pow(rRatio, 6) );
}

coor calcForcesP(Parameters &parameters, coor ri, coor rj)
{
    coor riMinusrj = subtractVectors(ri, rj);
    double rij = calcVectorModulus( riMinusrj );
    double rRatio = parameters.R/rij;

    double A = 12 * parameters.e * (pow( rRatio , 12 ) - pow(rRatio, 6) ) / (rij * rij);

    coor Fij = {A * riMinusrj.x,
                A * riMinusrj.y,
                A * riMinusrj.z};

    return Fij;
}

coor calcForcesS(Parameters &parameters, coor ri)
{
    double mri = calcVectorModulus(ri);
    coor fS = {0, 0, 0};

    if(mri > parameters.L)
    {
        double A  = parameters.f * (parameters.L - mri)  / mri;
        fS = {ri.x * A, ri.y * A, ri.z * A};
    }

    return fS;
}

void simulate(Parameters &parameters, State &state, ofstream & outputFileXYZ, ofstream & outputFileChar)
{
    double t = 0;

    for(int s = 1; s <= (parameters.S_0 + parameters.S_d); s++)
    {
        t = s*parameters.tau;
        updateState(parameters, state);

        if(!(s % parameters.S_out))
            outputChar(state, outputFileChar, t);

        if( !(s % parameters.S_xyz))
            outputXYZ(state.atoms, outputFileXYZ);
    }

}

void updateState(Parameters &parameters, State &state)
{

    for(int i = 0; i < (int) state.atoms.size(); i++)
    {
        state.atoms[i].p =
        {
            state.atoms[i].p.x + 0.5 * state.atoms[i].F.x * parameters.tau,
            state.atoms[i].p.y + 0.5 * state.atoms[i].F.y * parameters.tau,
            state.atoms[i].p.z + 0.5 * state.atoms[i].F.z * parameters.tau
        };

        state.atoms[i].r =
        {
            state.atoms[i].r.x + state.atoms[i].p.x * parameters.tau / parameters.m,
            state.atoms[i].r.y + state.atoms[i].p.y * parameters.tau / parameters.m,
            state.atoms[i].r.z + state.atoms[i].p.z * parameters.tau / parameters.m
        };
    }

    setPotentialForcesAndPressure(parameters, state);

    for(int i = 0; i < (int) state.atoms.size(); i++)
    {
        state.atoms[i].p =
        {
            state.atoms[i].p.x  + 0.5 * state.atoms[i].F.x * parameters.tau,
            state.atoms[i].p.y  + 0.5 * state.atoms[i].F.y * parameters.tau,
            state.atoms[i].p.z  + 0.5 * state.atoms[i].F.z * parameters.tau
        };
    }

    setEnergyAndTemperature(parameters, state);
}

void setEnergyAndTemperature(Parameters &parameters, State &state)
{
    double N = pow(parameters.n, 3);

    state.H = 0;

    for(int i = 0; i < N; i++)
    {
        state.H += calcVectorModulus(state.atoms[i].p) *
                calcVectorModulus(state.atoms[i].p) / 2. / parameters.m;
    }

    state.T = state.H * 2. / 3. / N / kB;
    state.H += state.V;
}

double getUniRandom()
{
    return (double)rand() / (RAND_MAX);
}

int getPlusOrMinus()
{
    return (rand() % 2) ? 1 : -1;
}

double calcVectorModulus(coor vector)
{
    return sqrt(vector.x * vector.x +
                vector.y * vector.y +
                vector.z * vector.z);
}

coor subtractVectors(coor v1, coor v2)
{
    return { v1.x - v2.x, v1.y - v2.y, v1.z - v2.z};
}

coor addVectors(coor v1, coor v2)
{
    return { v1.x + v2.x, v1.y + v2.y, v1.z + v2.z};
}

coor calcOppositeVector(coor v)
{
    return {-v.x, -v.y, -v.z};
}

void outputXYZ(vector<atom> & atoms, ofstream &outputFile)
{
    outputFile << atoms.size() << endl << "comment"<< endl;

    for(int i = 0; i < (int) atoms.size(); i++)
    {
        outputFile << "Ar " << atoms[i].r.x << " " << atoms[i].r.y << " "
                   << atoms[i].r.z << endl;
    }
}

void outputMomentum(vector<atom> &atoms)
{
    ofstream outputFile;
    outputFile.open("outputMomentum.txt");

    for(int i = 0; i < (int) atoms.size(); i++)
    {
        outputFile << atoms[i].p.x << "\t" << atoms[i].p.y << "\t"
                   << atoms[i].p.z << endl;
    }

    outputFile.close();
}

void outputChar(State & state, ofstream &outputFile, double &t)
{
    outputFile << t << "\t" << state.H << "\t" << state.V << "\t"
               << state.T << "\t" << state.P << endl;
}
