#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <iostream>
#include "parser.hpp"
#include "m_matrix.hpp"

int window_size_X = 640;
int window_size_Y = 480;

void glfwWindowSizeCallback(GLFWwindow *window, int width, int height) {
    window_size_X = width;
    window_size_Y = height;
    glViewport(0, 0, window_size_X, window_size_Y); //Определяем размеры окна в котором хотим рисовать
}

void glfwKeyCallback(GLFWwindow *window, int key, int scancode, int action, int mode) {
    if(key == GLFW_KEY_ESCAPE && action == GLFW_PRESS) {
        glfwSetWindowShouldClose(window, GL_TRUE);
    }
}
int main(void) {

    //Initialize the library
    if(!glfwInit()) {
        return -1;
    }

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 6);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    // Create a windowed mode window and its OpenGL context 
    GLFWwindow *window = glfwCreateWindow(window_size_X, window_size_Y, "Calc_elastic_2D", nullptr, nullptr);
    if (!window)
    {
        std::cout << "glfwCreateWindow failed!" << std::endl;
        glfwTerminate();
        return -1;
    }

    glfwSetWindowSizeCallback(window, glfwWindowSizeCallback);
    glfwSetKeyCallback(window, glfwKeyCallback);

    // Make the window's context current 
    glfwMakeContextCurrent(window);

    if(!gladLoadGL()) {
        std::cout << "Can't load GLAD!" <<std::endl;
        return -1;
    }

    std::cout << "Renderer: " << glGetString(GL_RENDERER) << std::endl;
    std::cout << "OpenGL version: " << glGetString(GL_VERSION) << std::endl;

    glClearColor(1, 1, 0, 1);

    // Loop until the user closes the window 
    while (!glfwWindowShouldClose(window))
    {
        // Render here 
        glClear(GL_COLOR_BUFFER_BIT);

        // Swap front and back buffers 
        glfwSwapBuffers(window);

        // Poll for and process events 
        glfwPollEvents();
    }

    glfwTerminate();

    printMatrix(get_nodes_coord("C:\\Projects\\Calc_elastic_2D\\src\\data\\mesh_for_calc.inc"));
    return 0;
}