#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <iostream>
#include "Renderer/ShaderProgram.hpp"
#include "Resources/ResourceManager.hpp"

GLfloat point[] = {
    0.0f,  0.5f,  0.0f,
    0.5f, -0.5f,  0.0f,
   -0.5f, -0.5f,  0.0f
};

GLfloat colors[] = {
    1.0f, 0.0f, 0.0f,
    0.0f, 1.0f, 0.0f,
    0.0f, 0.0f, 1.0f
};


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
int main(int argc, char** argv) {

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
    
    //Добавим облась видимости чтобы ResourceManager удалился раньше, чем контекст OpenGL
    {
        //Передаем путь к exe-файлу
        ResourceManager resourceManager(argv[0]);
        auto pDefaultShaderProgram = resourceManager.loadShaders("DefaultShader", "res/shaders/vertex.txt", "res/shaders/fragment.txt");
        if(!pDefaultShaderProgram) {
            std::cerr << "Can't create shader program" << "DefaultShader" <<std::endl;
            return -1;
        }

        GLuint points_vbo = 0;
        glGenBuffers(1, &points_vbo);
        glBindBuffer(GL_ARRAY_BUFFER, points_vbo);
        glBufferData(GL_ARRAY_BUFFER, sizeof(point), point, GL_STATIC_DRAW);

        GLuint colors_vbo = 0;
        glGenBuffers(1, &colors_vbo);
        glBindBuffer(GL_ARRAY_BUFFER, colors_vbo);
        glBufferData(GL_ARRAY_BUFFER, sizeof(colors), colors, GL_STATIC_DRAW);

        GLuint vao = 0;
        glGenVertexArrays(1, &vao);
        glBindVertexArray(vao);

        glEnableVertexAttribArray(0);
        glBindBuffer(GL_ARRAY_BUFFER, points_vbo);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, nullptr);

        glEnableVertexAttribArray(1);
        glBindBuffer(GL_ARRAY_BUFFER, colors_vbo);
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, nullptr);


        // Loop until the user closes the window 
        while (!glfwWindowShouldClose(window))
        {
            // Render here 
            glClear(GL_COLOR_BUFFER_BIT);

            pDefaultShaderProgram->use();
            glBindVertexArray(vao); 
            glDrawArrays(GL_TRIANGLES, 0, 3);

            // Swap front and back buffers 
            glfwSwapBuffers(window);

            // Poll for and process events 
            glfwPollEvents();
        }
    }

    glfwTerminate();

    //printMatrix(get_nodes_coord("C:\\Projects\\Calc_elastic_2D\\src\\data\\mesh_for_calc.inc"));
    return 0;
}