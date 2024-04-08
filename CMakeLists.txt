#Проверка на минимально необходимую версию CMake
cmake_minimum_required(VERSION 3.29 FATAL_ERROR) 

#Установка имени проекта, как Calc_elastic_2D 
set(PROJECT_NAME Calc_elastic_2D)
project(${PROJECT_NAME}) 

#Создание исполняемого файла с именем проекта и указание исходного файла main.cpp
add_executable(${PROJECT_NAME} src/main.cpp) 

#Устанавливаем флаги для библиотеки GLFW
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(GLFW_INSTALL OFF CACHE BOOL "" FORCE)

#Добавляем поддиректорию со сборкой и установкой библиотек GLFW
add_subdirectory(glfw)
#Линкуем библиотеку GLFW с проектом
target_link_libraries(${PROJECT_NAME} glfw)

