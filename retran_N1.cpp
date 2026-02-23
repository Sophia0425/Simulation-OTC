#include <iostream>
#include <fstream>
#include <random>
#include <cmath>
#include <limits>
#include <vector>
#include <algorithm>
#include <unordered_map>

using namespace std;

// Parameters
const int N0 = 323000;
const double lambda0 = 0.0906;   // years^-1
const double k_growth = 0.0217;  // years^-1
const int menopause_threshold = 650;
const int num_simulations = 5000;
const double removal_age = 20.0;
const double return_age  = 45.0;
const double removal_fraction = 0.30;
const double return_fraction_of_removed = 0.10;

// Output paths
const string MENOPAUSE_OUT = "your path";

// Exact integrated hazard for lambda(t)=lambda0*exp(k t) over [t, t+dt]
inline double integrated_lambda(double t, double dt) {
    double ek_t  = std::exp(k_growth * t);
    double ek_tp = std::exp(k_growth * (t + dt));
    return (lambda0 / k_growth) * (ek_tp - ek_t);
}

// Choose small dt
inline double choose_dt(double t, double p_target, double dt_max) {
    double lam = lambda0 * std::exp(k_growth * t);
    if (lam <= 0.0) return dt_max;

    // p_target = 1 - exp(-lam*dt)  =>  dt = -ln(1-p_target)/lam
    double dt = -std::log(1.0 - p_target) / lam;
    if (!std::isfinite(dt) || dt <= 0.0) dt = dt_max;
    if (dt > dt_max) dt = dt_max;
    return dt;
}

double simulate_menopause_age_population(mt19937 &rng, std::ofstream *traj_file) {
    double age = 0.0;
    int N1 = N0;
    int N2 = 0;
    bool removed_done = false;
    bool returned_done = false;
    int removed_amount_N1 = 0;
    const double p_target = 1e-3;   // smaller -> closer to event-by-event, slower
    const double dt_max   = 0.25;   // years

    const double EPS = 1e-12;

    while (N1 + N2 > menopause_threshold) {
        // Next intervention time
        double next_intervention = std::numeric_limits<double>::infinity();
        if (!removed_done) next_intervention = removal_age;
        else if (!returned_done) next_intervention = return_age;

        // Choose dt
        double dt = choose_dt(age, p_target, dt_max);
        if (age + dt > next_intervention) {
            dt = next_intervention - age; // can be 0
        }

        // If dt is effectively 0, we are at an intervention time: apply it and continue
        if (dt <= EPS) {
            age = next_intervention;

            if (!removed_done && std::abs(age - removal_age) < 1e-9) {
                // REMOVE ONLY FROM N1
                removed_amount_N1 = static_cast<int>(std::round(removal_fraction * static_cast<double>(N1+N2)));
                if (removed_amount_N1 > N1) removed_amount_N1 = N1;
                if (removed_amount_N1 < 0) removed_amount_N1 = 0;

                N1 -= removed_amount_N1;
                removed_done = true;

            } else if (removed_done && !returned_done && std::abs(age - return_age) < 1e-9) {
                // RETURN ONLY INTO N1 
                int return_amount = static_cast<int>(std::round(return_fraction_of_removed * removed_amount_N1));
                if (return_amount < 0) return_amount = 0;

                N1 += return_amount;
                returned_done = true;
            }
            continue;
        }

        double I = integrated_lambda(age, dt);
        double p = 1.0 - std::exp(-I);
        if (p < 0.0) p = 0.0;
        if (p > 1.0) p = 1.0;

        // N1 -> N2
        int x12 = 0;
        if (N1 > 0 && p > 0.0) {
            std::binomial_distribution<int> bin12(N1, p);
            x12 = bin12(rng);
            if (x12 > N1) x12 = N1;
        }
        N1 -= x12;
        N2 += x12;

        // N2 -> 0
        int x20 = 0;
        if (N2 > 0 && p > 0.0) {
            std::binomial_distribution<int> bin20(N2, p);
            x20 = bin20(rng);
            if (x20 > N2) x20 = N2;
        }
        N2 -= x20;

        age += dt;

        // Safety
        if (age > 200.0) break;
    }

    return age;
}

int main() {
    random_device rd;
    mt19937 rng(rd());

    // Open menopause output file
    ofstream menopause_file(MENOPAUSE_OUT);
    if (!menopause_file) {
        cerr << "Error opening menopause output file.\n";
        return 1;
    }

    // Run simulations
    for (int i = 0; i < num_simulations; i++) {
        double menopause_age = simulate_menopause_age_population(rng, nullptr);
        menopause_file << menopause_age << "\n";
    }

    menopause_file.close();
    cout << "Simulation complete.\n";
    return 0;
}
