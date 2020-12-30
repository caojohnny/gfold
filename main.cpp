#include <epigraph.hpp>
#include <mpl/matplotlibcpp.h>

#include "gfold.h"

namespace mpl = matplotlibcpp;

int main() {
    // These values can be found from the first numerical
    // example from the G-FOLD reference in the README

    gfold::gfold_lander_data lander_data = {
            .m_0 = 2000,
            .m_f = 300,
            .rho_1 = 0.2 * 24000,
            .rho_2 = 0.8 * 24000,
            .alpha = 5e-4,
            .V_max = 90
    };

    gfold::gfold_trajectory_data trajectory_data = {
            .r_0 = {2400, 450, -330},
            .r_dot_0 = {-10, -40, 10},
            .gamma_gs = M_PI / 4,
            .theta = M_PI / 4,
            .n_hat = {1, 0, 0},
            .q = {0, 0}
    };

    gfold::gfold_planet_data planet_data = {
            .g = {-3.71, 0, 0},
            .omega = {2.53e-5, 0.0, 6.62e-5}
    };

    gfold::gfold gfold{lander_data, trajectory_data, planet_data, 60, 1};
    gfold.compute();

    const Eigen::MatrixXd &x = gfold.get_state_vector();

    std::vector<double> vec_x;
    std::vector<double> vec_y;
    std::vector<double> vec_z;

    for (int i = 0; i < gfold.get_time_steps(); ++i) {
        Eigen::Vector3d r_i = x.topRows(3).col(i);
        vec_x.push_back(r_i.x());
        vec_y.push_back(r_i.y());
        vec_z.push_back(r_i.z());
    }

    mpl::plot3(vec_y, vec_z, vec_x);
    mpl::xlabel("Y Position (m)");
    mpl::ylabel("Z Position (m)");
    mpl::set_zlabel("X Position (m)");
    mpl::title("Position vs. Time");
    mpl::grid(true);

    mpl::show();

    return 0;
}
