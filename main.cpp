#include <iostream>

#include <epigraph.hpp>
#include <mpl/matplotlibcpp.h>

namespace mpl = matplotlibcpp;

int main() {
    double m_0 = 2000.0; // Initial mass, kg
    double m_f = 300.0; // Fuel mass, kg
    double T_max = 24000.0; // Maximum thrust, N
    double rho_1 = 0.2 * T_max; // Thrust @ minimum throttle, N
    double rho_2 = 0.8 * T_max; // Thrust @ maximum throttle, N
    double alpha = 5e-4; // "Mass depletion" constant, s/m
    double V_max = 90; // Max velocity, m/s
    double gamma_gs = M_PI / 4; // Lowest glide slope angle, rad
    double theta = M_PI / 4; // Thrust pointing angle, rad

    Eigen::Vector3d r_0{2400.0, 450.0, -330.0}; // Initial position, m
    Eigen::Vector3d r_dot_0{-10.0, -40.0, 10.0}; // Initial velocity, m/s

    Eigen::Vector2d q{0.0, 0.0}; // Target, m (technically supposed to be in R^2)
    Eigen::Vector3d g{-3.71, 0.0, 0.0}; // Mars gravitational constant, m/s^2

    Eigen::Vector3d omega{2.53e-5, 0.0, 6.62e-5}; // Mars angular velocity, rad/s^2

    // TODO: Not fixed a priori
    double t_f = 60; // Flight t_f, s
    double dt = 1; // Simulation time step, s

    int steps = static_cast<int>(t_f / dt); // Number of steps in the sim
    int i_f = steps - 1;

    Eigen::Vector3d e_1{1, 0, 0};
    Eigen::Vector3d e_2{0, 1, 0};

    Eigen::Matrix3d S_omega = (Eigen::Matrix3d{} << 0, -omega[2], omega[1],
            omega[2], 0, -omega[0],
            -omega[1], omega[0], 0)
            .finished();
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(6, 6);
    A.topRightCorner<3, 3>() = Eigen::Matrix3d::Identity();
    A.bottomLeftCorner<3, 3>() = -(S_omega * S_omega);
    A.bottomRightCorner<3, 3>() = -2 * S_omega;

    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(6, 3);
    B.bottomLeftCorner<3, 3>() = Eigen::Matrix3d::Identity();

    Eigen::MatrixXd E = (Eigen::Matrix<double, 2, 3>{} << 0, 1, 0,
            0, 0, 1).finished();
    Eigen::Vector3d c = e_1 / tan(gamma_gs);
    Eigen::Vector3d n_hat{1, 0, 0}; // Thrust unit pointing vector

    cvx::OptimizationProblem socp;

    // Add optimization variables
    const cvx::Scalar &obj = socp.addVariable("obj");
    const cvx::MatrixX &x = socp.addVariable("x", 6, steps);
    const cvx::MatrixX &r = x.topRows(3);
    const cvx::MatrixX &r_dot = x.bottomRows(3);
    const cvx::MatrixX &m = socp.addVariable("m", 1, steps);

    // 3.1 Change of Variables
    // const cvx::MatrixX &slack = socp.addVariable("GAMMA", 1, steps);
    // const cvx::MatrixX &T = socp.addVariable("T", 3, steps);
    const cvx::MatrixX &sigma = socp.addVariable("sigma", 1, steps);
    const cvx::MatrixX &u = socp.addVariable("u", 3, steps);
    const cvx::MatrixX &z = socp.addVariable("z", 1, steps);

    // Set up objective
    socp.addCostTerm(obj); // Problem 3
    socp.addConstraint(cvx::lessThan((cvx::par(E) * r.col(i_f) - cvx::par(q)).norm(), obj));

    // Set up constraints
    socp.addConstraint(cvx::equalTo(m.col(0).sum(), m_0)); // Eq. (7) - Initial mass constraint
    socp.addConstraint(cvx::box(0, m_0 - m_f, m.col(i_f).sum())); // Eq. (7) - Fuel use constraint

    socp.addConstraint(cvx::equalTo(r.col(0), cvx::par(r_0))); // Eq. (8) - Initial position constraint
    socp.addConstraint(cvx::equalTo(r_dot.col(0), cvx::par(r_dot_0))); // Eq. (8) - Initial velocity constraint

    socp.addConstraint(
            cvx::equalTo((cvx::par(e_1).transpose() * r.col(i_f)).sum(), 0)); // Eq. (9) - Final altitude constraint
    socp.addConstraint(cvx::equalTo(r_dot.col(i_f), 0)); // Eq. (9) - Final velocity constraint

    for (int i = 0; i < steps; ++i) {
        socp.addConstraint(cvx::lessThan(r_dot.col(i).norm(), cvx::par(V_max))); // Eq. (5) - Velocity constraint

        if (i != i_f) {
            const cvx::VectorX &diff_r = r.col(i) - r.col(i_f);
            socp.addConstraint(cvx::lessThan((cvx::par(E) * diff_r).norm() - (cvx::par(c).transpose() * diff_r).sum(),
                                             0)); // Eq. (5) - Glide slope constraint
        }

        if (i + 1 <= steps - 1) {
            // Section 3.1 Change of Variables
            // cvx::VectorX x_dot = cvx::par(A) * x.col(i) + cvx::par(B) * (cvx::par(g) + (T.col(i) / m.col(i)[0]));
            // socp.addConstraint(cvx::equalTo(x.col(i) + x_dot * cvx::par(dt), x.col(i + 1))); // Eq. (17)
            // socp.addConstraint(cvx::equalTo(m.col(i)[0] - cvx::par(alpha) * slack.col(i)[0] * cvx::par(dt), m.col(i + 1)[0])); // Eq. (17)

            cvx::VectorX x_dot = cvx::par(A) * x.col(i) + cvx::par(B) * (cvx::par(g) + u.col(i)); // Eq. (17)
            socp.addConstraint(cvx::equalTo(x.col(i) + x_dot * cvx::par(dt),
                                            x.col(i + 1))); // Eq. (17) - System dynamics constraint
            socp.addConstraint(cvx::equalTo(z.col(i).sum() - cvx::par(alpha) * sigma.col(i).sum(),
                                            z.col(i + 1))); // Eq. (33) - System mass constraint
        }

        // Section 3.1 Change of Variables
        // socp.addConstraint(cvx::lessThan(T.col(i).norm(), slack.col(i)[0])); // Eq. (18)
        // socp.addConstraint(cvx::box(rho_1, slack.col(i)[0], rho_2)); // Eq. (18)
        // socp.addConstraint(cvx::greaterThan((cvx::par(n_hat).transpose() * T.col(1)).sum(), cvx::par(cos(theta)) * slack.col(i)[0])); // Eq. (19)
        socp.addConstraint(cvx::lessThan(u.col(i).norm(),
                                         sigma.col(i).sum())); // Eq. (34) - Thrust slack constraint
        socp.addConstraint(cvx::greaterThan(
                cvx::par(n_hat).transpose() * u.col(i),
                cvx::par(std::cos(theta)) * sigma.col(i).sum())); // Eq. (34) - Thrust pointing constraint

        double z_0 = std::log(m_0 + m_f - alpha * rho_2 * i * dt);
        double exp_z_0 = std::exp(z_0);
        const cvx::Scalar &z_diff = z.col(i).sum() - cvx::par(z_0);
        // TODO: Squared dynamic var
        // socp.addConstraint(cvx::box(cvx::par(rho_1 * exp_z_0) * (cvx::par(1) - z_diff + (z_diff * z_diff * cvx::par(2))),
        // sigma.col(i)[0],
        // cvx::par(rho_2 * exp_z_0) * (cvx::par(1) - z_diff))); // Eq. (37)
        socp.addConstraint(cvx::lessThan(
                sigma.col(i).sum(),
                cvx::par(rho_2 * exp_z_0) * (cvx::par(1) - z_diff))); // Eq. (36) - Maximum thrust constraint
    }

    cvx::ecos::ECOSSolver solver{socp};
    solver.solve(false);
    Eigen::Vector2d d_P3 = cvx::eval(cvx::par(E) * r.col(i_f));
    double target = (d_P3 - q).norm();

    // Plotting
    std::vector<double> vec_x;
    std::vector<double> vec_y;
    std::vector<double> vec_z;

    for (int i = 0; i < steps; ++i) {
        Eigen::Vector3d r_i = cvx::eval(r.col(i));
        vec_x.push_back(r_i.x());
        vec_y.push_back(r_i.y());
        vec_z.push_back(r_i.z());
    }

    mpl::plot3(vec_y, vec_z, vec_x);
    mpl::xlabel("Y Position (m)");
    mpl::ylabel("Z Position (m)");
    mpl::set_zlabel("X Position (m)");
    mpl::title("P3: Position vs. Time");
    mpl::grid(true);

    socp.addCostTerm(sigma.sum()); // Problem 4
    socp.addConstraint(cvx::lessThan((cvx::par(E) * r.col(i_f) - cvx::par(q)).norm(), cvx::par(target))); // Eq. (20)

    solver.solve(false);

    // Plotting
    vec_x.clear();
    vec_y.clear();
    vec_z.clear();

    for (int i = 0; i < steps; ++i) {
        Eigen::Vector3d r_i = cvx::eval(r.col(i));
        vec_x.push_back(r_i.x());
        vec_y.push_back(r_i.y());
        vec_z.push_back(r_i.z());
    }

    mpl::plot3(vec_y, vec_z, vec_x);
    mpl::xlabel("Y Position (m)");
    mpl::ylabel("Z Position (m)");
    mpl::set_zlabel("X Position (m)");
    mpl::title("P4: Position vs. Time");
    mpl::grid(true);

    mpl::show();

    return 0;
}
