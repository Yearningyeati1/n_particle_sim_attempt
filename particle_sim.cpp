#include <SFML/Graphics.hpp>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <iostream>

// Define a Particle class
class Particle {
public:
    double mass;
    std::vector<double> position;
    std::vector<double> velocity;
    std::vector<double> acceleration;

    Particle(double m, const std::vector<double>& pos, const std::vector<double>& vel)
        : mass(m), position(pos), velocity(vel), acceleration(3, 0.0) {}

    void update(double dt) {
        for (int i = 0; i < 3; ++i) {
            velocity[i] += acceleration[i] * dt;
            position[i] += velocity[i] * dt;
            acceleration[i] = 0;  // Reset acceleration after each update
        }
    }

    void render(sf::RenderWindow &window) {
        sf::CircleShape shape(2.0f); // Radius of 2 units for the particle
        shape.setPosition(static_cast<float>(position[0]), static_cast<float>(position[1])); // Assuming a 2D simulation for simplicity
        shape.setFillColor(sf::Color::White);
        window.draw(shape);
    }
};

// Initialize particles
std::vector<Particle> initializeParticles(int n) {
    std::vector<Particle> particles;
    for (int i = 0; i < n; ++i) {
        double mass = 1.0;
        std::vector<double> position(3);
        position[0] = static_cast<double>(rand() % 800);
        position[1] = static_cast<double>(rand() % 600);
        position[2] = 0.0; // Assuming a 2D simulation
        std::vector<double> velocity(3);
        velocity[0] = static_cast<double>((rand() % 100) - 50) / 10.0; // Random initial velocity
        velocity[1] = static_cast<double>((rand() % 100) - 50) / 10.0;
        velocity[2] = 0.0;
        particles.push_back(Particle(mass, position, velocity));
    }
    return particles;
}

// Compute forces between particles
void computeForces(std::vector<Particle>& particles) {
    const double G = 6.67430e-11;  // Gravitational constant
    for (size_t i = 0; i < particles.size(); ++i) {
        for (size_t j = i + 1; j < particles.size(); ++j) {
            std::vector<double> direction(3);
            double distance = 0.0;
            for (int k = 0; k < 3; ++k) {
                direction[k] = particles[j].position[k] - particles[i].position[k];
                distance += std::pow(direction[k], 2);
            }
            distance = std::sqrt(distance);
            if (distance == 0) continue; // Avoid division by zero
            double force = G * particles[i].mass * particles[j].mass / std::pow(distance, 2);
            for (int k = 0; k < 3; ++k) {
                double acceleration = force / particles[i].mass;
                particles[i].acceleration[k] += acceleration * (direction[k] / distance);
                particles[j].acceleration[k] -= acceleration * (direction[k] / distance);
            }
        }
    }
}

// Main function
int main() {
    int n = 100;  // Number of particles
    double timeStep = 0.01;  // Time step
    int steps = 1000;  // Number of simulation steps

    std::vector<Particle> particles = initializeParticles(n);

    sf::RenderWindow window(sf::VideoMode(800, 600), "N-Particle Simulation");

    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();
        }

        computeForces(particles);
        for (auto& particle : particles) {
            particle.update(timeStep);
        }

        window.clear();
        for (auto& particle : particles) {
            particle.render(window);
        }
        window.display();
    }

    return 0;
}
