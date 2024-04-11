#Проверка на минимально необходимую версию CMake
cmake_minimum_required(VERSION 3.29 FATAL_ERROR) 

#Установка имени проекта, как Calc_elastic_2D 
set(PROJECT_NAME Calc_elastic_2D)
project(${PROJECT_NAME}) 

#Создание исполняемого файла с именем проекта и указание исходного файла main.cpp
add_executable(${PROJECT_NAME} src/main.cpp 
                               src/parser.hpp 
                               src/m_matrix.hpp
                               src/Renderer/ShaderProgram.cpp 
                               src/Renderer/ShaderProgram.hpp
                               src/Resources/ResourceManager.cpp
                               src/Resources/ResourceManager.hpp)

target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_17)

#Устанавливаем флаги для библиотеки GLFW
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(GLFW_INSTALL OFF CACHE BOOL "" FORCE)

#Добавляем поддиректорию со сборкой и установкой библиотек GLFW
add_subdirectory(glfw)
#Линкуем библиотеку GLFW с проектом
target_link_libraries(${PROJECT_NAME} glfw)

add_subdirectory(glad)
target_link_libraries(${PROJECT_NAME} glad)

#переменная сохраняет куда CMake будет складировать exe-файл
set_target_properties(${PROJECT_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/)

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_directory
                   ${CMAKE_SOURCE_DIR}/res $<TARGET_FILE_DIR:${PROJECT_NAME}>/res)
