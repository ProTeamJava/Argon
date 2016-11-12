/*
 * main.cpp
 *
 *  Created on: Oct 5, 2016
 *      Author: mongos
 */

#include "argon.h"

int main(int /* argc */, char* /* argv */[])
{
    State state;
    Parameters parameters = getParameters("data.txt");

    setInitialState(parameters, state);

    ofstream outputFileXYZ;
    outputFileXYZ.open("output.xyz");

    ofstream outputFileChar;
    outputFileChar.open("output.txt");

    outputFileChar << "t" << "\t" << "H" << "\t" << "V" << "\t" << "T" << "\t " << "P" << endl;

    outputXYZ(state.atoms, outputFileXYZ);

    outputMomentum(state.atoms);

    simulate(parameters, state, outputFileXYZ, outputFileChar);

    /* Nie potrzebne wykonywane w destruktorze
     * outputFileXYZ.close();
     * outputFileChar.close();
    */

    return 0;
}


