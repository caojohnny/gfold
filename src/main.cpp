#include <mpl/matplotlibcpp.h>

#include "gfold.h"

namespace mpl = matplotlibcpp;

int main() {
    // These values can be found from the 2nd numerical
    // example from the G-FOLD reference in the README

    double dt = 1; // Time step, s
    double T_max = 24000; // Max thrust, N

    gfold::gfold_lander_data lander_data = {
            .m_0 = 2000,
            .m_f = 300,
            .rho_1 = 0.2 * T_max,
            .rho_2 = 0.8 * T_max,
            .alpha = 5e-4,
            .V_max = 90
    };

    gfold::gfold_trajectory_data trajectory_data = {
            .r_0 = {2400, 3400, 0},
            .r_dot_0 = {-40, 45, 0},
            .gamma_gs = M_PI / 6,
            .theta = 2 * M_PI / 3,
            .n_hat = {1, 0, 0},
            .q = {0, 0}
    };

    gfold::gfold_planet_data planet_data = {
            .g = {-3.71, 0, 0},
            .omega = {2.53e-5, 0.0, 6.62e-5}
    };

    gfold::gfold gfold{lander_data, trajectory_data, planet_data, 45, 60, 1, dt};
    if (!gfold.compute()) {
        std::cout << "Infeasible" << std::endl;
        return 1;
    }

    const Eigen::MatrixXd &x = gfold.get_state_vector();
    const Eigen::MatrixXd &u = gfold.get_u();
    const Eigen::MatrixXd &z = gfold.get_z();
    
    std::vector<double> vec_t;

    std::vector<double> vec_x;
    std::vector<double> vec_y;
    std::vector<double> vec_z;

    std::vector<double> vec_throt;

    for (int i = 0; i < gfold.get_time_steps(); ++i) {
        Eigen::Vector3d r_i = x.topRows(3).col(i);
        Eigen::Vector3d u_i = u.col(i);
        double m = std::exp(z.col(i).sum());
        
        vec_t.push_back(i * dt);
        vec_x.push_back(r_i.x());
        vec_y.push_back(r_i.y());
        vec_z.push_back(r_i.z());
        vec_throt.push_back((u_i * m).norm() * 100 / T_max);
    }

    /* mpl::plot(vec_t, vec_throt);
    mpl::xlabel("Time (s)");
    mpl::ylabel("Throttle %");
    mpl::title("Throttle vs. Time");
    mpl::grid(true); */

    mpl::plot3(vec_y, vec_z, vec_x);
    mpl::xlabel("Y Position (m)");
    mpl::ylabel("Z Position (m)");
    mpl::set_zlabel("X Position (m)");
    mpl::title("Position vs. Time");
    mpl::grid(true);

    mpl::show();

    return 0;
}
