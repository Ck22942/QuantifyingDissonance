#pragma once

#define SAMPLE_RATE 44100
#define FRAMES_PER_BUFFER 4096
#define HEARING_REF 0.00002

#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>
#include <GLFW/glfw3.h> 
#include <stdio.h>
#include <iostream>
#include <portaudio.h>
#include <fftw3.h>
#include <array>
#include <cmath>
#include <algorithm>
#include <tuple>
#include <vector>
#include <string>
#include <limits>
#include <numeric>
#include <cstring>

using wComplex = std::array<float, 2>;
using realVec = std::vector<float>;
using complexVec = std::vector<wComplex>;
using realMatrix = std::array<realVec , 47>;
using complexMatrix = std::array<complexVec, 47>;

constexpr float INF = std::numeric_limits<float>::infinity();

static void glfw_error_callback(int error, const char *description);
GLFWwindow* SetupGLFWAndImGui(const char* windowTitle, int width, int height, const char* glsl_version);
void startframe(int& display_w , int& display_h, GLFWwindow* window);
void finishframe(int& display_w , int& display_h, GLFWwindow* window);
void cleanup(GLFWwindow* window);

