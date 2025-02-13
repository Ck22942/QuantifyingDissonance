#pragma once

#include <guistuff.h>



float AmptoDecibel(float);
float DecibeltoAmp(float);
float getMagnitude(fftwf_complex);
float calculateUpperSlope(float frequency, float level);
int calculateindex(float frequency);
float computeRMS(float* buffer , size_t N , float eDC);

struct FFTData {

    float volume;
    //r2c
    float* TimeBuffer;
    fftwf_complex* FrequencyBuffer;
    fftwf_plan plan;
    fftwf_plan Invplan;

    //c2c
    fftwf_complex* eTimeBuffer;
    fftwf_complex* eFrequencyBuffer;
    fftwf_plan eplan;
    fftwf_plan eInvplan;

    std::array<float,47> RMSvals;
    complexVec PhaseSpectrum;
    realVec AmplitudeSpectrum;
    realMatrix FilteredEnvelopes;
    float roughness;

    FFTData(size_t numFrames);
    ~FFTData();
};

float ComputeRoughness(const float* sample , FFTData* sd , int numFrames , int sampleRate);
float computeVolume(const float * buffer , size_t N);
float CrossCor(const std::vector<float>& a, const std::vector<float>& b);
void ShowFrequencySpectrum(FFTData& , int& , int& , int);


float interpolate(float ,std::vector<float>&, std::vector<float>&);