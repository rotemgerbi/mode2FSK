#include <iostream>
#include <complex>
#include <vector>
#include <fstream>
void MyNewFile();


//This is function FSK modulation

///input /*only*/vector date

using namespace std;

void FskMod(int LengthFille) {


    ////???

    double M = 2;
    double Freq_sep = 5e3;
    double nSamp = 25;
    double Fs = 200e3;

    ifstream x("myfile.bin", ios::binary);
    int rows=LengthFille;
    /*int rows = x.tellg();*/
    vector<double> buf(rows);
    x.read(reinterpret_cast<char *>(buf.data()), buf.size());
    x.close();


    float Samptime = 1 / Fs;
    int kp = 0;
    //obtain the total number of channels

    vector<complex<float>> y(rows * nSamp);

    vector<vector<double>> OneSamples((nSamp, 1));
    for (int counter = 0; counter < nSamp; counter++) OneSamples[counter][1] == 1;

    //Initialize the phase increments and the oscillator phase for modulator with
    //discontinuous phase

    double Fixed = ((6.2832 * (Freq_sep) / (2 * Samptime)));

    vector<vector<float>> PhaseIncr((nSamp, 2));

    for (int i = 0; i < nSamp; i++) {
        PhaseIncr[i * (-M - 1) * Fixed][1];
        PhaseIncr[i * (M - 1) * Fixed][2];
    }

    float phIncrSym[1][2];//start
    float phIncrSamp[1][2];//End

    //First place vector
    for (int po = 0; po < 2; po++) {
        for (int z = 1; z < 3; z++) {
            phIncrSym[1][po] = PhaseIncr[z][1];
        }
    }

    //Last place vector
    for (int po = 0; po < 2; po++) {
        for (int z = nSamp; z < nSamp - 2; z--) {
            phIncrSym[1][po] = PhaseIncr[z][1];
        }

        vector<vector<float>> OscPhase((1, M));
        double sizel = nSamp * rows;
        vector<vector<double>> Phase((sizel, 1));
        int prevPhase = 0;

        for (double iSym = 0; iSym < rows; iSym++) {
            int ph1 = prevPhase;
            for (double l(nSamp * (iSym - 1) + 1); l < nSamp * iSym; l++) {
                double p = buf[kp + 1];
                Phase[l][1] = ph1 * (OneSamples[nSamp][1]) + PhaseIncr[p][1];
            }
        }

        for (int g = 0; g < (nSamp * rows); g++) {
            y[g] = 2.7183 * (Phase[g][1] * 1i);
        }


        ofstream op("ru.bin", ios::binary | ios::in);
        size_t SizeVector = nSamp * rows;
        complex<double> OnePra = y[0];

        op.write(reinterpret_cast<char *>(y.data()), sizeof(OnePra) * SizeVector);
        op.close();


    }



}



