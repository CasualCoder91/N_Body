//#include <cstdlib>
#include <iostream>
#include <random>
#include <vector>
#include <chrono> //for timer
#include <string.h>

#include "Star.h"
#include "Node.h"
#include "InOut.h"
#include "Integrator.h"
#include "InitialConditions.h"
#include "Analysis.h"
//#include "Window.h"
//#include "ShaderProgram.h"

using namespace std::chrono;

int main() {

//#pragma omp parallel
//	{
//		std::cout<<"test\n";
//	}
	//init stars
	int n_Stars = 100;
	double boxLength = 1; //[]?
	double dt = 1;
	int NTimesteps = 100000;
	steady_clock::time_point startTime = steady_clock::now();
	//Window window = Window(n_Stars,boxLength*10);
	//window.createWindow();
	//window.initOpenGL();

	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_real_distribution<> dis(0.0, boxLength);
	std::vector<Star*> stars = {};
	//Init
	double totalMass = InitialConditions::InitialMass(stars,n_Stars);
	InitialConditions::PlummerSphere(stars, 1, totalMass);
	//stars.push_back(new Star(1000, 0, 0, 0));
	//totalMass += 1000;
	//window.initBuffers(stars);
	//glfwSetTime(0);
	//float current_time;
	//float last_time = glfwGetTime();
	//float dt1;
	//glEnable(GL_BLEND);
	//ShaderProgram shader;
	//std::string vertexShader = "vertexShader.txt";
	//char* cstr = new char[vertexShader.length() + 1];
	//strcpy(cstr, vertexShader.c_str());
	//std::string fragmentShader = "fragmentShader.txt";
	//char* cstr2 = new char[fragmentShader.length() + 1];
	//strcpy(cstr2, fragmentShader.c_str());
	//std::string ProjectionMatrix = "ProjectionMatrix";
	//char* cstr3 = new char[ProjectionMatrix.length() + 1];
	//strcpy(cstr3, ProjectionMatrix.c_str());
	//window.projectionMatrix.createProjectionMatrix(70.0f, 400.0f / 400.0f, 0.1f, 100.0f);
	//shader.loadMat4(cstr3, &window.projectionMatrix);

	int frameMod = 0;
	//for (int i = 0; i < n_Stars; i++) {
	//	stars.push_back( new Star(1, dis(gen), dis(gen), dis(gen)));
	//}
	//Integrate
	Integrator euler = Integrator(dt);
	for (int i = 0; i < NTimesteps; i++) {

		//SYSTEMTIME start;
		//GetSystemTime(&start);

		//current_time = glfwGetTime();
		//dt1 = current_time - last_time;
		//last_time = current_time;


		//if (window.keys[GLFW_KEY_T]) {
		//	window.keys[GLFW_KEY_T] = false;
		//	window.afterImages = !window.afterImages;
		//}

		Vec3D tlf = Vec3D(), brb = Vec3D();
		Node::FindCorners(tlf, brb, stars);
		Node root = Node(tlf, brb, nullptr);
		for (Star* star : stars) {
			root.Insert(star);
		}
		root.CalculateMassDistribution();
		#pragma omp parallel for //1:10
		for (int i = 0; i < stars.size();++i){//(Star* star : stars) {
			stars.at(i)->acceleration = Vec3D(); // reset acceleration to 0,0,0
			root.ApplyForce(stars.at(i));
		}
		//for (Star* star : stars) {
		//	star->acceleration = Vec3D(); // reset acceleration to 0,0,0
		//	root.ApplyForce(star);
		//}
		euler.Euler(stars);

		//window.Input(dt1);
		//window.UpdateBuffers(stars);
		//window.camera.update();
		//shader.use();
		//glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
		//glBlendFunc(GL_ONE, GL_ZERO);
		//glBindVertexArray(window.VAO);
		//glDrawArraysInstanced(GL_POINTS, 0, 1, n_Stars);
		//glfwPollEvents();
		//glfwSwapBuffers(window.window);

		if (i % 100 == 0) {
			InOut::WriteWithLabel(stars, "./Output/stars" + std::to_string(i) + ".dat");
			//InOut::WriteAll(stars, "./Output/stars_all" + std::to_string(i) + ".dat");
			//double potentialEnergy = Analysis::PotentialEnergy(stars);
			//double kineticEnergy = Analysis::KineticEnergy(stars);
			//std::cout<< "Kinetic Energy: " + std::to_string(kineticEnergy) << std::endl;
			//std::cout << "Potential Energy: " + std::to_string(potentialEnergy) << std::endl;
			//std::cout << "Total Energy: " + std::to_string(kineticEnergy+potentialEnergy) << std::endl << std::endl;
		}
	}
	steady_clock::time_point endTime = steady_clock::now();
	duration<double> time_span = duration_cast<duration<double>>(startTime - endTime);
	//InOut::Write(stars,"stars.dat");
	//InOut::Write(&root);
	std::cout << "Time needed" << time_span.count();
	std::cout << "done" << std::endl;
	std::cin.get();
	return 0;
}