import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

# === Initial Data from Problem 5-11 ===
R1 = 2
R2 = 6
R3 = 4
R4 = 40
R5 = 9
R6 = 12
R7_pipe = 50  # Resistance of pipe 7 itself
D7 = 350      # mm
H1 = 100      # m (Initial head)
Hset = 80     # m (Set point for PRV)
Kprv = 10     # Loss coefficient for PRV
g = 9.81      # m/s^2

# === Calculate Rprv (Resistance of Pressure Reducing Valve) ===
Rprv = (8 * Kprv) / (g * np.pi**2 * (D7 / 1000)**4)
print(f"Calculated Rprv = {Rprv:.2f}")

# === Effective resistance when PRV is not operating ===
R7_not_operating = R7_pipe + Rprv

# === Initial Guesses ===
Q_init_full = np.array([0.1] * 7)
Q_init_one_way = np.array([0.1] * 6)

# === Head Loss Calculator ===
def calculate_head_losses(Q, scenario_name):
    """Calculate and display head losses H = R*Q^2 for each pipe"""
    print(f"\n=== Head Losses for {scenario_name} ===")
    
    if len(Q) == 7:
        R_values = [R1, R2, R3, R4, R5, R6, R7_pipe]
        pipe_names = ['Pipe 1', 'Pipe 2', 'Pipe 3', 'Pipe 4', 'Pipe 5', 'Pipe 6', 'Pipe 7']
    elif len(Q) == 6:
        R_values = [R1, R2, R3, R4, R5, R6]
        pipe_names = ['Pipe 1', 'Pipe 2', 'Pipe 3', 'Pipe 4', 'Pipe 5', 'Pipe 6']
    
    total_head_loss = 0
    for i in range(len(Q)):
        H_i = R_values[i] * Q[i]**2
        total_head_loss += H_i
        print(f"{pipe_names[i]}: H{i+1} = {R_values[i]} × ({Q[i]:.5f})² = {H_i:.5f} m")
    
    # For one-way valve case, explicitly show Q7 = 0
    if len(Q) == 6:
        H7 = R7_pipe * 0**2
        print(f"Pipe 7: H7 = {R7_pipe} × (0.00000)² = {H7:.5f} m (Q7 = 0)")
    
    print(f"Total Head Loss = {total_head_loss:.5f} m")
    print("-" * 60)

def calculate_nodal_heads(Q, scenario_name):
    """Calculate heads at each node based on the image equations"""
    print(f"\n=== Nodal Head Calculations for {scenario_name} ===")
    
    # Starting from H1 = 100 m
    if len(Q) >= 6:
        H2 = H1 - R1 * Q[0]**2
        H3 = H1 - R1 * Q[0]**2 - R3 * Q[2]**2
        H4 = H1 - R1 * Q[0]**2 - R2 * Q[1]**2
        H5 = H4 - R5 * Q[4]**2
        H6 = H1 - R1 * Q[0]**2 - R2 * Q[1]**2 - R5 * Q[4]**2 - R6 * Q[5]**2
        
        print(f"H1 = {H1:.3f} m (given)")
        print(f"H2 = H1 - R1×Q1² = {H1} - {R1}×({Q[0]:.5f})² = {H2:.3f} m")
        print(f"H3 = H1 - R1×Q1² - R3×Q3² = {H2:.3f} - {R3}×({Q[2]:.5f})² = {H3:.3f} m")
        print(f"H4 = H1 - R1×Q1² - R2×Q2² = {H2:.3f} - {R2}×({Q[1]:.5f})² = {H4:.3f} m")
        print(f"H5 = H4 - R5×Q5² = {H4:.3f} - {R5}×({Q[4]:.5f})² = {H5:.3f} m")
        print(f"H6 = H5 - R6×Q6² = {H5:.3f} - {R6}×({Q[5]:.5f})² = {H6:.3f} m")
        
        if len(Q) == 7:
            if scenario_name == "PRV Not Operating":
                # When PRV is not operating, calculate based on flow direction
                # If Q7 > 0, flow is downstream, so head loss occurs
                if Q[6] >= 0:
                    H_downstream = H6 - R7_not_operating * Q[6]**2
                    print(f"H_downstream = H6 - (R7_total)×Q7² = {H6:.3f} - {R7_not_operating:.2f}×({Q[6]:.5f})² = {H_downstream:.3f} m")
                else:
                    H_downstream = H6 + R7_not_operating * Q[6]**2
                    print(f"H_downstream = H6 + (R7_total)×|Q7|² = {H6:.3f} + {R7_not_operating:.2f}×({abs(Q[6]):.5f})² = {H_downstream:.3f} m")
            else:
                # For PRV operating case
                if Q[6] >= 0:
                    H_downstream = H6 - R7_pipe * Q[6]**2
                    print(f"H_downstream = H6 - R7×Q7² = {H6:.3f} - {R7_pipe}×({Q[6]:.5f})² = {H_downstream:.3f} m")
                else:
                    H_downstream = H6 + R7_pipe * Q[6]**2
                    print(f"H_downstream = H6 + R7×|Q7|² = {H6:.3f} + {R7_pipe}×({abs(Q[6]):.5f})² = {H_downstream:.3f} m")
        
        # Check PRV conditions
        if scenario_name == "PRV Operating":
            print(f"\nPRV Check: H_downstream should equal Hset = {Hset} m")
            if len(Q) == 7:
                print(f"Actual H_downstream = {H_downstream:.3f} m")
                print(f"Difference = {abs(H_downstream - Hset):.3f} m")
        elif scenario_name == "PRV Not Operating" and len(Q) == 7:
            print(f"\nPRV Not Operating: H_downstream = {H_downstream:.3f} m (no regulation)")
            print(f"This is {'higher' if H_downstream > Hset else 'lower'} than Hset = {Hset} m")
            print(f"Difference from Hset = {abs(H_downstream - Hset):.3f} m")
    
    print("-" * 60)

# === Scenario 1: PRV Operating ===
def equations_prv_operating(Q):
    Q1, Q2, Q3, Q4, Q5, Q6, Q7 = Q
    
    # Mass balance equations
    eq1 = Q1 - Q2 - Q3 - 0.3                    # Node 2
    eq2 = Q3 - Q4 - Q7 - 0.2                    # Node 3
    eq3 = Q2 + Q4 - Q5 - 0.5                    # Node 4
    eq4 = Q5 - Q6 - 0.4                         # Node 5
    eq5 = Q6 + Q7 - 0.6                         # Node 6
    
    # Energy equations
    eq6 = R2 * Q2**2 - R4 * Q4**2 - R3 * Q3**2  # Loop equation
    eq7 = R1 * Q1**2 + R2 * Q2**2 + R5 * Q5**2 + R6 * Q6**2 - R7_pipe * Q7**2 + (Hset - H1)
    
    return [eq1, eq2, eq3, eq4, eq5, eq6, eq7]

# === Scenario 2: PRV Not Operating ===
def equations_prv_not_operating(Q):
    Q1, Q2, Q3, Q4, Q5, Q6, Q7 = Q
    
    # Mass balance equations (same as operating)
    eq1 = Q1 - Q2 - Q3 - 0.3
    eq2 = Q3 - Q4 - Q7 - 0.2
    eq3 = Q2 + Q4 - Q5 - 0.5
    eq4 = Q5 - Q6 - 0.4
    eq5 = Q6 + Q7 - 0.6
    
    # Energy equations
    eq6 = R2 * Q2**2 - R4 * Q4**2 - R3 * Q3**2
    # When PRV is not operating, use the same loop equation as your original code:
    # This is the energy balance around the loop containing pipes 4, 5, 6, 7
    eq7 = R4 * Q4**2 + R5 * Q5**2 + R6 * Q6**2 - R7_not_operating * Q7**2
    
    return [eq1, eq2, eq3, eq4, eq5, eq6, eq7]

# === Scenario 3: PRV as One-Way Valve (Q7 = 0) ===
def equations_one_way_valve(Q):
    Q1, Q2, Q3, Q4, Q5, Q6 = Q
    
    # Mass balance equations with Q7 = 0
    eq1 = Q1 - Q2 - Q3 - 0.3
    eq2 = Q3 - Q4 - 0.2          # Q7 = 0
    eq3 = Q2 + Q4 - Q5 - 0.5
    eq4 = Q5 - Q6 - 0.4
    eq5 = Q6 - 0.6               # Q7 = 0
    eq6 = R2 * Q2**2 - R4 * Q4**2 - R3 * Q3**2
    
    return [eq1, eq2, eq3, eq4, eq5, eq6]

# === Solver Functions ===
def solve_system(equations_func, Q_init, scenario_name, max_iter=50, tol=1e-6):
    """Generic solver for all scenarios"""
    try:
        Q_solution = fsolve(equations_func, Q_init, xtol=tol)
        
        # Verify solution
        residuals = equations_func(Q_solution)
        max_residual = max(abs(r) for r in residuals)
        
        if max_residual < tol:
            print(f"✓ Converged for {scenario_name}")
            print(f"  Maximum residual: {max_residual:.2e}")
            return Q_solution
        else:
            print(f"⚠ Solution may not be accurate for {scenario_name}")
            print(f"  Maximum residual: {max_residual:.2e}")
            return Q_solution
            
    except Exception as e:
        print(f"✗ Error solving {scenario_name}: {e}")
        return None

# === Main Execution ===
if __name__ == "__main__":
    print("PIPE NETWORK ANALYSIS WITH PRESSURE REDUCING VALVE")
    print("=" * 60)
    print(f"Initial Head H1 = {H1} m")
    print(f"PRV Set Point Hset = {Hset} m")
    print(f"PRV Resistance Rprv = {Rprv:.2f}")
    print("=" * 60)
    
    results = {}
    
    # === Scenario 1: PRV Operating ===
    print("\n🔧 SCENARIO 1: PRV OPERATING")
    Q_operating = solve_system(equations_prv_operating, Q_init_full, "PRV Operating")
    if Q_operating is not None:
        results['PRV Operating'] = Q_operating
        print("\nFlow Rates (m³/s):")
        for i, q in enumerate(Q_operating, 1):
            print(f"  Q{i} = {q:8.5f}")
        calculate_head_losses(Q_operating, "PRV Operating")
        calculate_nodal_heads(Q_operating, "PRV Operating")
    
    # === Scenario 2: PRV Not Operating ===
    print("\n🔧 SCENARIO 2: PRV NOT OPERATING")
    Q_not_operating = solve_system(equations_prv_not_operating, Q_init_full, "PRV Not Operating")
    if Q_not_operating is not None:
        results['PRV Not Operating'] = Q_not_operating
        print("\nFlow Rates (m³/s):")
        for i, q in enumerate(Q_not_operating, 1):
            print(f"  Q{i} = {q:8.5f}")
        calculate_head_losses(Q_not_operating, "PRV Not Operating")
        calculate_nodal_heads(Q_not_operating, "PRV Not Operating")
    
    # === Scenario 3: PRV as One-Way Valve ===
    print("\n🔧 SCENARIO 3: PRV AS ONE-WAY VALVE (Q7 = 0)")
    Q_one_way_6 = solve_system(equations_one_way_valve, Q_init_one_way, "One-Way Valve")
    if Q_one_way_6 is not None:
        Q_one_way = np.append(Q_one_way_6, 0.0)  # Add Q7 = 0
        results['One-Way Valve'] = Q_one_way
        print("\nFlow Rates (m³/s):")
        for i, q in enumerate(Q_one_way, 1):
            print(f"  Q{i} = {q:8.5f}")
        calculate_head_losses(Q_one_way, "One-Way Valve")
        calculate_nodal_heads(Q_one_way, "One-Way Valve")
    
    # === Summary Comparison ===
    print("\n" + "="*80)
    print("SUMMARY COMPARISON")
    print("="*80)
    print(f"{'Scenario':<20} {'Q1':<10} {'Q2':<10} {'Q3':<10} {'Q4':<10} {'Q5':<10} {'Q6':<10} {'Q7':<10}")
    print("-"*80)
    for scenario, Q in results.items():
        print(f"{scenario:<20}", end="")
        for q in Q:
            print(f"{q:10.5f}", end="")
        print()
    
    print("\n Analysis Complete!")