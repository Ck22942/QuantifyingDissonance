#include <roughness.h>
//
//
//
//
//data!!!
constexpr std::array<int, 50> barkBounds = {{
    0,    50,  100,  150, 200,  250, 300,  350, 400,  450, 510,  570, 630,  700, 770,  840,
    920,  1000, 1080, 1170, 1270, 1370, 1480, 1600, 1720, 1850, 2000, 2150, 2320, 2500, 2700, 2900,
    3150, 3400, 3700, 4000, 4400, 4800,  5300, 5800, 6400, 7000, 7700, 8500, 9500, 10500,12000, 13500, 15500, 20000 
}};

std::array<std::vector<float>, 2> a0tab = {{
    {0.0f, 10.0f, 12.0f, 13.0f, 14.0f, 15.0f, 16.0f, 16.5f, 17.0f, 18.0f, 18.5f, 19.0f, 20.0f, 21.0f, 21.5f, 22.0f, 22.5f, 23.0f, 23.5f, 24.0f, 25.0f, 26.0f},
    {0.0f, 0.0f, 1.15f, 2.31f, 3.85f, 5.62f, 6.92f, 7.38f, 6.92f, 4.23f, 2.31f, 0.0f, -1.43f, -2.59f, -3.57f, -5.19f, -7.41f, -11.3f, -20.0f, -40.0f, -130.0f, -999.0f}
}};

std::array<std::vector<float>, 2> HTres = {{
        {0.0f, 0.01f, 0.17f, 0.8f, 1.0f, 1.5f, 2.0f, 3.3f, 4.0f, 5.0f, 6.0f, 8.0f, 10.0f, 12.0f, 13.3f, 15.0f, 16.0f, 17.0f, 18.0f, 19.0f, 20.0f, 21.0f, 22.0f, 23.0f, 24.0f, 24.5f, 25.0f},
        {130.0f, 70.0f, 60.0f, 30.0f, 25.0f, 20.0f, 15.0f, 10.0f, 8.1f, 6.3f, 5.0f, 3.5f, 2.5f, 1.7f, 0.0f, -2.5f, -4.0f, -3.7f, -1.5f, 1.4f, 3.8f, 5.0f, 7.5f, 15.0f, 48.0f, 60.0f, 130.0f}
}};

std::array<std::vector<float>, 2> H2 = {{
    {0.0f, 17.0f, 23.0f, 25.0f, 32.0f, 37.0f, 48.0f, 67.0f, 90.0f, 114.0f, 171.0f, 206.0f, 247.0f, 294.0f, 358.0f},
    {0.0f, 0.8f, 0.95f, 0.975f, 1.0f, 0.975f, 0.9f, 0.8f, 0.7f, 0.6f, 0.4f, 0.3f, 0.2f, 0.1f, 0.0f}
}};

std::array<std::vector<float>, 2> H5 = {{
    {0.0f, 32.0f, 43.0f, 56.0f, 69.0f, 92.0f, 120.0f, 142.0f, 165.0f, 231.0f, 277.0f, 331.0f, 397.0f, 502.0f},
    {0.0f, 0.8f, 0.95f, 1.0f, 0.975f, 0.9f, 0.8f, 0.7f, 0.6f, 0.4f, 0.3f, 0.2f, 0.1f, 0.0f}
}};

std::array<std::vector<float>, 2> H16 = {{
    {0.0f, 23.5f, 34.0f, 47.0f, 56.0f, 63.0f, 79.0f, 100.0f, 115.0f, 135.0f, 159.0f, 172.0f, 194.0f, 215.0f, 244.0f, 290.0f, 348.0f, 415.0f, 500.0f, 645.0f},
    {0.0f, 0.4f, 0.6f, 0.8f, 0.9f, 0.95f, 1.0f, 0.975f, 0.95f, 0.9f, 0.85f, 0.8f, 0.7f, 0.6f, 0.5f, 0.4f, 0.3f, 0.2f, 0.1f, 0.0f}
}};

std::array<std::vector<float>, 2> H21 = {{
    {0.0f, 19.0f, 44.0f, 52.5f, 58.0f, 75.0f, 101.5f, 114.5f, 132.5f, 143.5f, 165.5f, 197.5f, 241.0f, 290.0f, 348.0f, 415.0f, 500.0f, 645.0f},
    {0.0f, 0.4f, 0.8f, 0.9f, 0.95f, 1.0f, 0.95f, 0.9f, 0.85f, 0.8f, 0.7f, 0.6f, 0.5f, 0.4f, 0.3f, 0.2f, 0.1f, 0.0f}
}};

std::array<std::vector<float>, 2> H42 = {{
    {0.0f, 15.0f, 41.0f, 49.0f, 53.0f, 64.0f, 71.0f, 88.0f, 94.0f, 106.0f, 115.0f, 137.0f, 180.0f, 238.0f, 290.0f, 348.0f, 415.0f, 500.0f, 645.0f},
    {0.0f, 0.4f, 0.8f, 0.9f, 0.965f, 0.99f, 1.0f, 0.95f, 0.9f, 0.85f, 0.8f, 0.7f, 0.6f, 0.5f, 0.4f, 0.3f, 0.2f, 0.1f, 0.0f}
}};

std::array<std::vector<float>, 2> gr = {{
    {0.0f, 1.0f, 2.5f, 4.9f, 6.5f, 8.0f, 9.0f, 10.0f, 11.0f, 11.5f, 13.0f, 17.5f, 21.0f, 24.0f},
    {0.0f, 0.35f, 0.7f, 0.7f, 1.1f, 1.1f, 1.25f, 1.26f, 1.18f, 1.08f, 1.0f, 0.66f, 0.46f, 0.38f, 0.3f}
}};

//
//
//
//
//methods!!!
FFTData::FFTData(size_t numFrames)
        : 
          TimeBuffer(fftwf_alloc_real(numFrames)),
          FrequencyBuffer(fftwf_alloc_complex(numFrames/2 + 1)),
          plan(fftwf_plan_dft_r2c_1d(numFrames, TimeBuffer, FrequencyBuffer, FFTW_ESTIMATE)) ,
          Invplan(fftwf_plan_dft_c2r_1d(numFrames, FrequencyBuffer, TimeBuffer , FFTW_ESTIMATE)),

          eTimeBuffer(fftwf_alloc_complex(numFrames)),
          eFrequencyBuffer(fftwf_alloc_complex(numFrames)),
          eplan(fftwf_plan_dft_1d(numFrames, eTimeBuffer, eFrequencyBuffer , FFTW_FORWARD ,FFTW_ESTIMATE)),
          eInvplan(fftwf_plan_dft_1d(numFrames, eFrequencyBuffer, eTimeBuffer , FFTW_BACKWARD ,FFTW_ESTIMATE))


          {
            PhaseSpectrum.resize(numFrames/2 + 1);
            AmplitudeSpectrum.resize(numFrames/2 + 1);
            for(auto& vec : FilteredEnvelopes){
                vec.resize(numFrames);}
            ;}

FFTData::~FFTData() {
        fftwf_destroy_plan(plan);
        fftwf_destroy_plan(Invplan);
        fftwf_free(TimeBuffer);
        fftwf_free(FrequencyBuffer);
         fftwf_destroy_plan(eplan);
        fftwf_destroy_plan(eInvplan);
        fftwf_free(eTimeBuffer);
        fftwf_free(eFrequencyBuffer);
       }
        
//
//
//
//
//
//
//

//algorithm here
float ComputeRoughness(const float* sample , FFTData* sd , int numFrames , int sampleRate){

    complexVec& PhaseSpectrum = sd->PhaseSpectrum;
    realVec& AmplitudeSpectrum = sd->AmplitudeSpectrum;
    realMatrix& FilteredEnvelopes = sd->FilteredEnvelopes;

    //step 1. Perform the FFT , filtering and converting to dB
    memcpy(sd->TimeBuffer , sample, numFrames * sizeof(float));
    fftwf_execute(sd->plan);
    
    float fRes = (float)sampleRate/(float)numFrames;
       
    for(size_t j = 0 ; j < numFrames/2 + 1 ; ++j){
        float BarkNo = (float)(calculateindex(fRes*j ))/2.0f;
        float a0 = interpolate(BarkNo, a0tab[0] , a0tab[1]);
        float HT = interpolate(BarkNo , HTres[0] , HTres[1]);
        sd->FrequencyBuffer[j][0] /= (float)numFrames;
        sd->FrequencyBuffer[j][1] /= (float)numFrames;
        AmplitudeSpectrum[j] = getMagnitude(sd->FrequencyBuffer[j]);
        AmplitudeSpectrum[j] = AmptoDecibel(AmplitudeSpectrum[j]) + a0;
        if(AmplitudeSpectrum[j] < HT){
            AmplitudeSpectrum[j] = -INF;
        }
        float _tempp = DecibeltoAmp(AmplitudeSpectrum[j]);
        if(_tempp == 0.0f)
        {
            PhaseSpectrum[j][0] = 0.0f;
            PhaseSpectrum[j][1] = 0.0f;
        }
        else
        {
            PhaseSpectrum[j][0] = sd->FrequencyBuffer[j][0] / _tempp;
            PhaseSpectrum[j][1] = sd->FrequencyBuffer[j][1] / _tempp;
        }
    }

    //step 2. Get excitation patterns
    float S1 = -27.0f;
    for(size_t i = 0; i < 47 ; ++i){
        int L_i = i - 1;
        int U_i = i + 1;
        for(size_t j = 0; j < numFrames/2 + 1  ; ++j){
            float dB1;
            
            float dB0 = AmplitudeSpectrum[j];
            float S2 = calculateUpperSlope(j * fRes , dB0);
            float k = (float)calculateindex(j * fRes) / 2.0f;
            float _a0 = interpolate(k, a0tab[0] , a0tab[1]);
            float HT = interpolate(k , HTres[0] , HTres[1]);
            if(dB0 == -INF || S2 == -INF)
            {
                dB1 = -INF;
            }
            else
            {
                if(k < L_i)
                {
                    dB1 = dB0 +   (0.5f * S2 * (float)(L_i - k));
                }
                else if(k < U_i)
                {
                    dB1 = dB0;
                }
                else
                {
                    dB1 = dB0 + (0.5f * S1 * (float)(k - U_i + 1));
                }

                if(dB1 < HT - _a0){
                    dB1 = -INF;
                }

            }
                   
            sd->FrequencyBuffer[j][0] = PhaseSpectrum[j][0] * DecibeltoAmp(dB1);
            sd->FrequencyBuffer[j][1] = PhaseSpectrum[j][1] * DecibeltoAmp(dB1);
            
        }
        //inverse transform
        fftwf_execute(sd->Invplan);


        for(size_t j = 0; j < numFrames  ; ++j){
            sd->eTimeBuffer[j][0] = sd->TimeBuffer[0];
            sd->eTimeBuffer[j][1] = 0.0f;
        }
        fftwf_execute(sd->eplan);

        for(size_t j = 0; j < numFrames  ; ++j){
            if(j > 0 && j < numFrames/2){
                sd->eFrequencyBuffer[j][0] *= 2.0f /(float)numFrames;
                sd->eFrequencyBuffer[j][1] *= 2.0f /(float)numFrames;
            }
            else if(j > numFrames/2){
                sd->eFrequencyBuffer[j][0] = 0.0f;
                sd->eFrequencyBuffer[j][1] = 0.0f;
            } 
            else{
                sd->eFrequencyBuffer[j][0] *= 1.0f /(float)numFrames;
                sd->eFrequencyBuffer[j][1] *= 1.0f /(float)numFrames;
            }
        }
        
        fftwf_execute(sd->eInvplan);
        for(size_t j = 0 ; j < numFrames ; ++j){
        fftwf_complex temp = {sd->TimeBuffer[j] , sd->eTimeBuffer[j][0]};
        sd->TimeBuffer[j] = getMagnitude(temp);}
    
        fftwf_execute(sd->plan);
        float eDC = getMagnitude(sd->FrequencyBuffer[0]) /float(numFrames);

        for(size_t j = 0; j < numFrames/2 + 1 ; ++j){
            float interscale;
            if(i < 5 ){
                interscale = interpolate(j * fRes , H2[0] , H2[1]);
            }
            else if(i < 16){
                interscale = interpolate(j * fRes , H5[0] , H5[1]);
            }
            else if(i < 21){
                interscale = interpolate(j * fRes , H16[0] , H16[1]);
            }
            else if(i < 42){
                interscale = interpolate(j * fRes , H21[0] , H21[1]);
            }
            else{
                interscale = interpolate(j * fRes , H42[0] , H42[1]);
            }
            sd->FrequencyBuffer[j][0] *= interscale/(float)numFrames;
            sd->FrequencyBuffer[j][1] *= interscale/(float)numFrames;
        } 

        fftwf_execute(sd->Invplan);
        sd->RMSvals[i] = computeRMS(sd->TimeBuffer , numFrames , eDC);
        memcpy(sd->FilteredEnvelopes[i].data() , sd->TimeBuffer , numFrames * sizeof(float));     
    };
    std::array<float,45> ck;
    std::array<float,47> gz;
    for(size_t i = 0 ; i < 45 ; ++i){
        ck[i] = CrossCor(sd->FilteredEnvelopes[i] , sd->FilteredEnvelopes[i + 2]);
    }
    for(size_t i = 0 ; i < 47 ; ++i){
        gz[i] = interpolate((float)i * 0.5f, gr[0] , gr[1]);
    }
    

    float Total_Roughness = 0.0f;
    for(size_t i = 0 ; i < 47 ; ++i){
        float r_i = 0.0f;
        if(i < 2){
           r_i = gz[i] * sd->RMSvals[i] * ck[i];
        } 
        else if(i < 45){
            r_i = gz[i] * sd->RMSvals[i] * ck[i] * ck[i-2];
        }
        else{
            r_i = gz[i] * sd->RMSvals[i]  * ck[i-2];
        }
        Total_Roughness += r_i * r_i;
    }

    return 0.25f * Total_Roughness;
}

//
//
//
//
//
//function definitions
//
//
float AmptoDecibel(float x){
    if(x == 0.0f)
    {
        return -INF;
    }
    return 20.0f * std::log10(x/HEARING_REF);
}

float DecibeltoAmp(float x){
    if(x == -INF)
    {
        return 0.0f;
    }
    return std::pow(10.0f , 0.05f * x * HEARING_REF);
}

float getMagnitude(fftwf_complex z){
    float Re = z[0];
    float Im = z[1];
    return sqrtf(Re * Re + Im * Im);
}

float calculateUpperSlope(float frequency, float level){
    if(frequency == 0.0f) return -INF;
    if(level == -INF) return -INF;
    return (-24.0f) -  (230.0f/frequency) + (0.2f*level);
}

int calculateindex(float frequency){
    int result = -1;
    for(int i = 0 ; i < barkBounds.size() ; ++i){
        if (frequency > barkBounds[i]){
            result = i;
        } else{
            break;
        }
    }
    return result;
}

float interpolate(float f ,std::vector<float>& invec, std::vector<float>& outvec){
    if(f < invec[0]){return outvec[0];}
    for (size_t i = 0; i < invec.size() - 1; ++i) {
        if (f >= invec[i] && f <= invec[i + 1]) {
            float t = (f - invec[i]) / (invec[i + 1] - invec[i]);
            return outvec[i] + t * (outvec[i + 1] - outvec[i]);
        }
    }
    return outvec[outvec.size()];
}

//
//
//
//
//plotting function for ImGui
void ShowFrequencySpectrum(FFTData& sd , int& height, int& width , int numFrames)
{ 
    ImGui::SetNextWindowPos(ImVec2(0, 0));
    ImGui::SetNextWindowSize(ImVec2(height , width));       
    ImGui::Begin("Fullscreen Window", NULL, ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize);
    int i = 0;
    int line = 1;
    while(line < 7 && i < 47){
        while(i < line*8 && i < 47){
            char buffer[5];
            sprintf(buffer, "##%d", i);
            char overlay0[64];
            snprintf(overlay0, sizeof(overlay0), "%d %.3f", i ,sd.RMSvals[i]);
            ImGui::BeginGroup();
            ImGui::PlotLines(buffer, sd.FilteredEnvelopes[i].data(), numFrames, 0, nullptr, -1.0f, 1.0f, ImVec2(200, 110));  
            ImGui::ProgressBar(sd.RMSvals[i], ImVec2(200, 20), overlay0);
            ImGui::EndGroup();
            ImGui::SameLine();
            ++i;
        }
        ImGui::NewLine();
        ++line;
    }
    char overlay[64];
    snprintf(overlay, sizeof(overlay), "vol: %.6f", sd.volume);
    ImGui::ProgressBar(sd.volume, ImVec2(200, 20), overlay);
    snprintf(overlay, sizeof(overlay), "Roughness: %.2f", sd.roughness);
    ImGui::ProgressBar(sd.roughness, ImVec2(200, 20), overlay);
    ImGui::End();
}

float computeRMS(float* buffer , size_t N , float eDC) {
    float sum = 0.0f;
    for (size_t i = 0 ; i < N ; ++i) {
        sum += buffer[i] * buffer[i];
    }
    float result = std::sqrt(sum / N);
    if( result <= 0.0f){
        return 0.0f;
    }
    if(eDC <= 0.0f){
        return 1.0f;
    }

    float nresult = result / eDC;
    if(result >= 1.0f){
        return 1.0f;
    }
    else{
        return nresult;
    }  
}

float computeVolume(const float* buffer , size_t N){
    float bigest = 0;
    for(size_t i = 0 ; i < N ; ++i){
        if(buffer[i] > bigest){
            bigest = buffer[i];
        }
    }
    return bigest;
}

float CrossCor(const std::vector<float>& a, const std::vector<float>& b) {
    float n = (float)a.size();
    //compute mean
    float mean_a = std::accumulate(a.begin(), a.end(), 0.0f) / n;
    float mean_b = std::accumulate(b.begin(), b.end(), 0.0f) / n;
    
    // Compute covariance
    float cov = std::inner_product(a.begin() , a.end() , b.begin() , 0.0f ,std::plus<>() ,[mean_a , mean_b](float x , float y){return (x - mean_a) * (y - mean_b);}) / n;

    // Compute variance
    float var_a = std::inner_product(a.begin(), a.end(), a.begin(), 0.0f , std::plus<>() ,[mean_a](float x , float y){return (x - mean_a) * (y - mean_a);}) / n;
    float var_b = std::inner_product(b.begin(), b.end(), b.begin(), 0.0f , std::plus<>() ,[mean_b](float x , float y){return (x - mean_b) * (y - mean_b);}) / n;

    float k = cov /std::sqrt(var_a * var_b);
    
    return std::max(0.0f , k);
}