#include <guistuff.h>

static void glfw_error_callback(int error, const char *description)
{
    fprintf(stderr, "Glfw Error %d: %s\n", error, description);
}

GLFWwindow* SetupGLFWAndImGui(const char* windowTitle, int width, int height, const char* glsl_version) {
    if (!glfwInit()) {
        std::cout << "not initialised";
        std::cin.get();
        return nullptr; 
    }
    
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);

    // Create GLFW window
    GLFWwindow* window = glfwCreateWindow(width, height, windowTitle, NULL, NULL);
    if (window == nullptr) {
        return window; 
    }
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1); 
    

    // Initialize ImGui
    IMGUI_CHECKVERSION();
    
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO();
    (void)io; 
    ImGui::StyleColorsDark();
   

    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init(glsl_version);


    return window; 
}


void startframe(int& display_w , int& display_h, GLFWwindow* window){
    glfwGetFramebufferSize(window, &display_w, &display_h);
    glfwPollEvents();
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();
}

void finishframe(int& display_w , int& display_h, GLFWwindow* window){
    ImGui::Render();
    glViewport(0, 0, display_w, display_h);
    glClearColor(0.45f, 0.55f, 0.60f, 1.00f);
    glClear(GL_COLOR_BUFFER_BIT);
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
    glfwSwapBuffers(window);
}

void cleanup(GLFWwindow* window){
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    glfwDestroyWindow(window);
    glfwTerminate();
}


