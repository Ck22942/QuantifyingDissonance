#include <roughness.h>

float ROUGHNESS;

static int audioCallback(const void* inputBuffer, void* outputBuffer,
                         unsigned long framesPerBuffer,
                         const PaStreamCallbackTimeInfo* timeInfo,
                         PaStreamCallbackFlags statusFlags,
                         void* userData) {

    const float* in = (const float*)inputBuffer;
    FFTData* _Rdata = (FFTData*)userData;

   _Rdata->roughness = ComputeRoughness(in , _Rdata , FRAMES_PER_BUFFER , SAMPLE_RATE);
    _Rdata->volume = computeVolume(in , FRAMES_PER_BUFFER);

    return paContinue;

}

int main() {
    size_t numFrames = FRAMES_PER_BUFFER;
    int sampleRate = SAMPLE_RATE;
    int height , width;
    Pa_Initialize();
    FFTData Rdata(numFrames);

    
    GLFWwindow* window = SetupGLFWAndImGui("Synth App TEST", 1920, 1080, "#version 130");

    PaStream* stream;
    Pa_OpenDefaultStream(&stream, 1, 1, paFloat32, sampleRate , numFrames, audioCallback, &Rdata);         
    Pa_StartStream(stream);
    

    while(!glfwWindowShouldClose(window)){
        startframe(width , height , window);
        
        ShowFrequencySpectrum(Rdata, width , height, numFrames);
        
        finishframe(width, height , window);   
    };
    
    cleanup(window);

    Pa_StopStream(stream);
    Pa_CloseStream(stream);
    Pa_Terminate();

    return 0;
}
