/*
 * TODO: complete this file comment.
 */
#include <iostream>
#include "SimpleGraph.h"
#include <fstream>
#include <cmath>
#include <vector>
#include <chrono>
#include <sstream>

constexpr double PI = 3.14159265358979323;
using namespace std;

void Welcome();
string GetLine();
int StringToInteger(const string&);
void initNodesCircle(SimpleGraph&);
void ComputeForces(const SimpleGraph&,
                   vector<double>&, vector<double>&,
                   double kRepel = 1e-3, double kAttract = 1e-3);
void MoveNodes(SimpleGraph&,
               const vector<double>&, const vector<double>&);
void RunLayout(SimpleGraph&, int);

int getSeconds(const string& prompt,
               const string& reprompt) {
    while (true) {
        cout << prompt;
        string line;
        if (!getline(cin, line)) {
            throw domain_error("getLine: End of input reached while waiting for line.");
        }
        istringstream iss(line);
        int value; char extra;
        if (iss >> value && !(iss >> extra)) return value;
        cerr << reprompt << endl;
    }

}

// Main method
int main() {
    Welcome();
    cout << "Please input Graph-filename: " << endl;

    SimpleGraph graph;

    while (true) {
        string filename = GetLine();
        ifstream infile(filename);
        int nodeCount;
        if (!infile.is_open() || !(infile >> nodeCount)) {
            cout << "Not such a file, please try again!" << endl;
        } else {
            graph.nodes.resize(nodeCount);
            int start{0}, end{0};
            while (infile >> start >> end) {
                Edge e{};
                e.start = start;
                e.end = end;
                graph.edges.push_back(e);
            }
            infile.close();
            break;
        }

    }

    initNodesCircle(graph);
    vector<double> deltaX(graph.nodes.size(), 0.0);
    vector<double> deltaY(graph.nodes.size(), 0.0);

    int seconds = getSeconds("Please input how long do you want to run?", "Again");


    RunLayout(graph, seconds);
    return 0;
}

/* Prints a message to the console welcoming the user and
 * describing the program. */
void Welcome() {
    cout << "Welcome to CS106L GraphViz!" << endl;
    cout << "This program uses a force-directed graph layout algorithm" << endl;
    cout << "to render sleek, snazzy pictures of various graphs." << endl;
    cout << endl;
}

string GetLine() {
    string response{""};
    getline(cin, response);
    return response;
}


void initNodesCircle(SimpleGraph& graph) {
    int n = graph.nodes.size();
    double r = 0.4;
    for (int k = 0; k < n; ++k) {
        double angle = 2.0 * PI * k / n;
        graph.nodes[k].x = 0.5 + r * cos(angle);
        graph.nodes[k].y = 0.5 + r * sin(angle);
    }
}

void ComputeForces(const SimpleGraph& graph,
                   vector<double>& deltaX, vector<double>& deltaY,
                   double kRepel, double kAttract){
    int n = graph.nodes.size();
    deltaX.assign(n, 0.0);
    deltaY.assign(n, 0.0);

    //     repulsive force:
    // For each pair of nodes (x0, y0), (x1, y1):
    //     Compute Frepel = krepel / sqrt ((y1 – y0)2 + (x1 – x0)2)
    //     Compute θ = atan2(y1 – y0, x1 – x0)
    //     Δx0 -= Frepel cos(θ)
    //     Δy0 -= Frepel sin(θ)
    //     Δx1 += Frepel cos(θ)
    //     Δy1 += Frepel sin(θ)

    for (int i = 0; i < n; ++i) {
        for (int j = i +1; j < n; ++j) {
            double dx = graph.nodes[j].x - graph.nodes[i].x;
            double dy = graph.nodes[j].y - graph.nodes[i].y;
            double dist = sqrt(dx*dx + dy*dy) + 1e-8;
            double force = kRepel / dist;
            double angle = atan2(dy, dx);
            double fx = force * cos(angle);
            double fy = force * sin(angle);
            deltaX[i] -= fx;
            deltaY[i] -= fy;
            deltaX[j] += fx;
            deltaY[j] += fy;
        }
    }

    //attractive force
    // For each edge:
    // Let the first endpoint be (x0, y0)
    // Let the second endpoint be (x1, y1)
    // Compute Fattract = kattract * ((y1 – y0)2 + (x1 – x0)2)
    // Compute θ = atan2(y1 – y0, x1 – x0)
    // Δx0 += Fattract cos(θ) // Note that this is a += and not a -=!
    // Δy0 += Fattract sin(θ)
    // Δx1 -= Fattract cos(θ) // Note that this is a -= and not a +=!
    // Δy1 -= Fattract sin(θ

    for(const auto& e: graph.edges) {
        int i = e.start;
        int j = e.end;
        double dx = graph.nodes[j].x - graph.nodes[i].x;
        double dy = graph.nodes[j].y - graph.nodes[i].y;
        double dist2 = dx*dx + dy*dy;
        double force = kAttract * dist2;
        double angle = atan2(dy, dx);
        double fx = force * cos(angle);
        double fy = force * sin(angle);
        deltaX[i] += fx;
        deltaY[i] += fy;
        deltaX[j] -= fx;
        deltaY[j] -= fy;
    }
}
void MoveNodes(SimpleGraph& graph,
               const vector<double>& deltaX, const vector<double>& deltaY) {
    int n = graph.nodes.size();
    for (int i = 0; i < n; ++i) {
        graph.nodes[i].x += deltaX[i];
        graph.nodes[i].y += deltaY[i];
    }
}
void RunLayout(SimpleGraph& graph, int second){
    InitGraphVisualizer(graph);
    DrawGraph(graph);

    vector<double> deltaX, deltaY;
    // get present time
    auto startTime = chrono::steady_clock::now();

    while(chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now() - startTime).count() < second) {
        ComputeForces(graph, deltaX, deltaY);
        MoveNodes(graph, deltaX, deltaY);
        DrawGraph(graph);
    }
}
