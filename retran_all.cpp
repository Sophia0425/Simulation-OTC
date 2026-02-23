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
const string TRAJ_DIR = "your path";

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

// Apply proportional removal/return to (N1, N2) based on total change dN
void apply_proportional_change(int &N1, int &N2, int dN) {
    int total = N1 + N2;
    if (dN == 0) return;

    if (dN < 0) {
        int remove = -dN;
        if (remove > total) remove = total;
        if (total <= 0) return;

        int remove1 = static_cast<int>(std::round(remove * (static_cast<double>(N1) / total)));
        if (remove1 > N1) remove1 = N1;
        if (remove1 < 0) remove1 = 0;

        int remove2 = remove - remove1;
        if (remove2 > N2) {
            remove2 = N2;
            remove1 = remove - remove2;
            if (remove1 > N1) remove1 = N1;
            if (remove1 < 0) remove1 = 0;
        }

        N1 -= remove1;
        N2 -= remove2;
    } else {
        int add = dN;
        if (add <= 0) return;

        if (total <= 0) {
            // If empty, add back to N1 by convention
            N1 += add;
            return;
        }

        int add1 = static_cast<int>(std::round(add * (static_cast<double>(N1) / total)));
        if (add1 < 0) add1 = 0;
        if (add1 > add) add1 = add;
        int add2 = add - add1;

        N1 += add1;
        N2 += add2;
    }
}

// writes trajectory 
double simulate_menopause_age_population(mt19937 &rng, std::ofstream *traj_file) {
    double age = 0.0;
    int N1 = N0;
    int N2 = 0;
    bool removed_done = false;
    bool returned_done = false;
    int removed_amount_total = 0;
    int removed_N1 = 0;
    int removed_N2 = 0;
    const double p_target = 1e-3;   // smaller -> closer to event-by-event, slower
    const double dt_max   = 0.25;   // years

    // ---- TRAJECTORY ----
    // auto log_state = [&](double t, int a, int b) {
    //     if (traj_file) (*traj_file) << t << "\t" << a << "\t" << b << "\n";
    // };
    //
    // if (traj_file) {
    //     (*traj_file) << "# age\tN1\tN2\n";
    //     log_state(age, N1, N2);
    // }
    // -------------------------------------------

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
                int total = N1 + N2;
                removed_amount_total = static_cast<int>(std::round(removal_fraction * total));
                if (removed_amount_total > total) removed_amount_total = total;
                if (removed_amount_total < 0) removed_amount_total = 0;

                removed_N1 = 0;
                removed_N2 = 0;

                int remove = removed_amount_total;
                if (remove > total) remove = total;
                if (total > 0 && remove > 0) {
                    int remove1 = static_cast<int>(
                        std::round(remove * (static_cast<double>(N1) / total))
                    );
                    if (remove1 > N1) remove1 = N1;
                    if (remove1 < 0) remove1 = 0;

                    int remove2 = remove - remove1;
                    if (remove2 > N2) {
                        remove2 = N2;
                        remove1 = remove - remove2;
                        if (remove1 > N1) remove1 = N1;
                        if (remove1 < 0) remove1 = 0;
                    }

                    N1 -= remove1;
                    N2 -= remove2;

                    removed_N1 = remove1;
                    removed_N2 = remove2;
                }

                removed_done = true;

                // ---- TRAJECTORY ----
                // log_state(age, N1, N2);
                // -------------------------------------------

            } else if (removed_done && !returned_done && std::abs(age - return_age) < 1e-9) {

                int return_total = static_cast<int>(
                    std::round(return_fraction_of_removed * removed_amount_total)
                );
                if (return_total < 0) return_total = 0;

                int return1 = static_cast<int>(
                    std::round(return_fraction_of_removed * removed_N1)
                );
                if (return1 < 0) return1 = 0;
                if (return1 > return_total) return1 = return_total;

                int return2 = return_total - return1;
                if (return2 < 0) return2 = 0;

                N1 += return1;
                N2 += return2;

                returned_done = true;

                // ---- TRAJECTORY ----
                // log_state(age, N1, N2);
                // -------------------------------------------
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

        // ---- TRAJECTORY ----
        // log_state(age, N1, N2);
        // -------------------------------------------

        // Safety
        if (age > 200.0) break;
    }

    return age;
}

int main() {
    random_device rd;
    mt19937 rng(rd());

    // ---- TRAJECTORY ----
    // // Pick 5 unique simulation indices to save trajectories for
    // vector<int> idx(num_simulations);
    // for (int i = 0; i < num_simulations; i++) idx[i] = i;
    // std::shuffle(idx.begin(), idx.end(), rng);
    // vector<int> save_idx(idx.begin(), idx.begin() + 5);
    // std::sort(save_idx.begin(), save_idx.end()); // optional: nice ordering
    // ------------------------------------------------------------

    ofstream menopause_file(MENOPAUSE_OUT);
    if (!menopause_file) {
        cerr << "Error opening menopause output file.\n";
        return 1;
    }

    // ---- TRAJECTORY ----
    // // Open trajectory files for the selected indices
    // unordered_map<int, ofstream> traj_files;
    // traj_files.reserve(5);
    //
    // for (int id : save_idx) {
    //     string path = TRAJ_DIR + "trajectory_all_20,40;50,10_" + to_string(id) + ".txt";
    //     traj_files.emplace(piecewise_construct,
    //                        forward_as_tuple(id),
    //                        forward_as_tuple(path));
    //     if (!traj_files[id]) {
    //         cerr << "Error opening trajectory file: " << path << "\n";
    //         return 1;
    //     }
    // }
    //
    // cout << "Saving trajectories for indices: ";
    // for (int id : save_idx) cout << id << " ";
    // cout << "\n";
    // --------------------------------------------

    for (int i = 0; i < num_simulations; i++) {
        // ---- TRAJECTORY ----
        // ofstream *traj_ptr = nullptr;
        // auto it = traj_files.find(i);
        // if (it != traj_files.end()) traj_ptr = &(it->second);
        // -------------------------------------------

        double menopause_age = simulate_menopause_age_population(rng, nullptr);
        menopause_file << menopause_age << "\n";
    }

    menopause_file.close();
    cout << "Simulation complete.\n";
    return 0;
}
