#include <guistuff.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <set>
#include <array>
#include <cmath>


using triad = std::array<float,2>;

std::set<triad> MakeFarey(int n){
    std::set<triad> s;    
    for(int i = n; i >= 0 ; --i){
        for(int j = n ; j > i ; --j){
            for(int k = n ; k >j ; --k){
                float r1 = i / (float) j;
                float r2 = j / (float) k;
                if(r2*r1 >= 0.5f){
                    triad _triad = { r1, r2};
                    s.insert(_triad);
                }
            }
        }
    }
    return s;
}

int main(){
    int height , width;
    std::set<triad> Ft25 = MakeFarey(45);
    int soize = Ft25.size();

    GLFWwindow* window = SetupGLFWAndImGui("Synth App TEST", 1920, 1080, "#version 130");

    while(!glfwWindowShouldClose(window)){

        startframe(width , height , window);


        ImGui::SetNextWindowPos(ImVec2(0, 0));
        ImGui::SetNextWindowSize(ImVec2(height , width));       
        ImGui::Begin("Fullscreen Window", NULL, ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize);

        ImGui::Text(std::to_string(soize).c_str());

        ImDrawList* draw_list = ImGui::GetWindowDrawList();
        for(auto Coord : Ft25){
            draw_list->AddCircleFilled(  ImVec2(50 + std::roundf(600 * -log2f(Coord[0])), 50 + std::roundf(600 * - log2f(Coord[1]))), 1, IM_COL32(255, 0, 0, 255));
        
        }

        ImGui::End();


        finishframe(width, height , window);
        
    }


    cleanup(window);


    return 0;
}