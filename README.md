Pipe Network Analysis with Pressure Reducing Valve (PRV)
This project models and solves a hydraulic pipe network using Python. The analysis includes different operating conditions for a Pressure Reducing Valve (PRV), and computes flow rates, head losses, and nodal heads throughout the system.

💡 Overview
The code simulates a water distribution system consisting of multiple pipes, junctions (nodes), and a pressure reducing valve (PRV). Using mass and energy conservation equations, the system solves for unknown pipe flow rates under three scenarios:

PRV Operating: PRV actively maintains a downstream head.

PRV Not Operating: PRV acts as a regular pipe (no control).

PRV as One-Way Valve: No flow allowed in reverse direction (Q7 = 0).

The system uses the Newton-Raphson method via scipy.optimize.fsolve to solve the nonlinear equations.

📐 Features
Solves nonlinear flow equations for different hydraulic scenarios.

Computes and prints:

Flow rates for each pipe.

Head loss in each pipe (using Darcy-Weisbach-like resistance).

Nodal pressures (heads).

Simulates behavior of PRV under different operating conditions.

Compares results across all scenarios.

🧮 Equations Used
Mass balance (Continuity) at each junction.

Energy conservation along closed loops.

PRV modeled using:

Set point head control.

Loss coefficient-based resistance when not operating.

🔧 Dependencies
Python 3.x

numpy

scipy

matplotlib (optional; currently included but not essential for CLI output)
