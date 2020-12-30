#include "gfold.h"

gfold::gfold::gfold(gfold_lander_data lander_data,
                    gfold_trajectory_data trajectory_data,
                    gfold_planet_data planet_data,
                    double min_duration,
                    double max_duration,
                    double duration_search_interval,
                    double dt) :
        lander_data(lander_data),
        trajectory_data(std::move(trajectory_data)),
        planet_data(std::move(planet_data)),
        min_duration(min_duration),
        max_duration(max_duration),
        duration_search_interval(duration_search_interval),
        dt(dt) {
}

gfold::gfold::gfold(gfold_lander_data lander_data,
                    gfold_trajectory_data trajectory_data,
                    gfold_planet_data planet_data,
                    double duration,
                    double dt) :
        gfold(lander_data, std::move(trajectory_data), std::move(planet_data), duration, duration, 1, dt) {
}

bool gfold::gfold::p3p4(double t_f) {
    int steps = static_cast<int>(t_f / dt); // Number of steps in the sim
    int i_f = steps - 1; // Final index

    Eigen::Vector3d e_1{1, 0, 0};
    Eigen::Vector3d e_2{0, 1, 0};

    Eigen::Matrix3d S_omega = (Eigen::Matrix3d{} << 0, -planet_data.omega[2], planet_data.omega[1],
            planet_data.omega[2], 0, -planet_data.omega[0],
            -planet_data.omega[1], planet_data.omega[0], 0)
            .finished();
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(6, 6);
    A.topRightCorner<3, 3>() = Eigen::Matrix3d::Identity();
    A.bottomLeftCorner<3, 3>() = -(S_omega * S_omega);
    A.bottomRightCorner<3, 3>() = -2 * S_omega;

    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(6, 3);
    B.bottomLeftCorner<3, 3>() = Eigen::Matrix3d::Identity();

    Eigen::MatrixXd E = (Eigen::Matrix<double, 2, 3>{} << 0, 1, 0,
            0, 0, 1).finished();
    Eigen::Vector3d c = e_1 / tan(trajectory_data.gamma_gs);

    cvx::OptimizationProblem socp;

    // Add optimization variables
    const cvx::Scalar &obj = socp.addVariable("obj");
    const cvx::MatrixX &x = socp.addVariable("x", 6, steps);
    const cvx::MatrixX &r = x.topRows(3);
    const cvx::MatrixX &r_dot = x.bottomRows(3);

    // 3.1 Change of Variables
    // const cvx::MatrixX &slack = socp.addVariable("GAMMA", 1, steps);
    // const cvx::MatrixX &T = socp.addVariable("T", 3, steps);
    // const cvx::MatrixX &m = socp.addVariable("m", 1, steps);
    const cvx::MatrixX &sigma = socp.addVariable("sigma", 1, steps);
    const cvx::MatrixX &u = socp.addVariable("u", 3, steps);
    const cvx::MatrixX &z = socp.addVariable("z", 1, steps);

    // Problem 3
    socp.addCostTerm(obj);
    socp.addConstraint(cvx::lessThan((cvx::par(E) * r.col(i_f) - cvx::par(trajectory_data.q)).norm(), obj));

    // 3.1 Change of Variables
    // socp.addConstraint(cvx::equalTo(m.col(0).sum(), m_wet));
    // socp.addConstraint(cvx::box(0, lander_data.m_0 - lander_data.m_f, m.col(i_f).sum()));
    // Eq. (7) - Initial mass constraint
    socp.addConstraint(cvx::equalTo(z.col(0).sum(), std::log(lander_data.m_0)));
    // Eq. (7) - Fuel use constraint
    socp.addConstraint(cvx::box(0, std::log(lander_data.m_0 - lander_data.m_f), z.col(i_f).sum()));

    // Eq. (8) - Initial position constraint
    socp.addConstraint(cvx::equalTo(r.col(0), cvx::par(trajectory_data.r_0)));
    // Eq. (8) - Initial velocity constraint
    socp.addConstraint(cvx::equalTo(r_dot.col(0), cvx::par(trajectory_data.r_dot_0)));

    // Eq. (9) - Final altitude constraint
    socp.addConstraint(cvx::equalTo((cvx::par(e_1).transpose() * r.col(i_f)).sum(), 0));
    // Eq. (9) - Final velocity constraint
    socp.addConstraint(cvx::equalTo(r_dot.col(i_f), 0));

    for (int i = 0; i < steps; ++i) {
        // Eq. (5) - Velocity constraint
        socp.addConstraint(cvx::lessThan(r_dot.col(i).norm(), cvx::par(lander_data.V_max)));

        if (i != i_f) {
            // Eq. (5) - Glide slope constraint
            const cvx::VectorX &diff_r = r.col(i) - r.col(i_f);
            socp.addConstraint(cvx::lessThan((cvx::par(E) * diff_r).norm() - (cvx::par(c).transpose() * diff_r).sum(),
                                             0));
        }

        if (i + 1 <= steps - 1) {
            // Section 3.1 Change of Variables
            // cvx::VectorX x_dot = cvx::par(A) * x.col(i) + cvx::par(B) * (cvx::par(g) + (T.col(i) / m.col(i)[0]));
            // socp.addConstraint(cvx::equalTo(x.col(i) + x_dot * cvx::par(dt), x.col(i + 1))); // Eq. (17)
            // socp.addConstraint(cvx::equalTo(m.col(i)[0] - cvx::par(alpha) * slack.col(i)[0] * cvx::par(dt), m.col(i + 1)[0])); // Eq. (17)

            // Eq. (17) - System dynamics constraint
            cvx::VectorX x_dot = cvx::par(A) * x.col(i) + cvx::par(B) * (cvx::par(planet_data.g) + u.col(i));
            socp.addConstraint(cvx::equalTo(x.col(i) + x_dot * cvx::par(dt),
                                            x.col(i + 1)));
            // Eq. (33) - System mass constraint
            socp.addConstraint(cvx::equalTo(z.col(i).sum() - cvx::par(lander_data.alpha) * sigma.col(i).sum(),
                                            z.col(i + 1)));
        }

        // Section 3.1 Change of Variables
        // socp.addConstraint(cvx::lessThan(T.col(i).norm(), slack.col(i)[0])); // Eq. (18)
        // socp.addConstraint(cvx::box(rho_1, slack.col(i)[0], rho_2)); // Eq. (18)
        // socp.addConstraint(cvx::greaterThan((cvx::par(n_hat).transpose() * T.col(1)).sum(), cvx::par(cos(theta)) * slack.col(i)[0])); // Eq. (19)
        // Eq. (34) - Thrust slack constraint
        socp.addConstraint(cvx::lessThan(u.col(i).norm(),
                                         sigma.col(i).sum()));
        // Eq. (34) - Thrust pointing constraint
        socp.addConstraint(cvx::greaterThan(cvx::par(trajectory_data.n_hat).transpose() * u.col(i),
                                            cvx::par(std::cos(trajectory_data.theta)) * sigma.col(i).sum()));

        double z_0 = std::log(lander_data.m_0 - lander_data.alpha * lander_data.rho_2 * i * dt);
        double exp_z_0 = std::exp(-z_0);
        const cvx::Scalar &z_diff = z.col(i).sum() - cvx::par(z_0);
        // TODO: Squared dynamic var
        // socp.addConstraint(cvx::box(
        //         cvx::par(lander_data.rho_1 * exp_z_0) * (cvx::par(1) - z_diff + (z_diff * z_diff * cvx::par(2))),
        //         sigma.col(i)[0],
        //         cvx::par(lander_data.rho_2 * exp_z_0) * (cvx::par(1) - z_diff)));
        // Eq. (36) - Maximum thrust constraint
        socp.addConstraint(cvx::lessThan(
                sigma.col(i).sum(),
                cvx::par(lander_data.rho_2 * exp_z_0) * (cvx::par(1) - z_diff)));
    }

    cvx::ecos::ECOSSolver solver{socp};
    if (!solver.solve(false)) {
        return false;
    }

    if (solver.getExitCode()) {
        return false;
    }

    Eigen::Vector2d d_P3 = cvx::eval(cvx::par(E) * r.col(i_f));
    double target = (d_P3 - trajectory_data.q).norm();

    // Problem 4
    socp.addCostTerm(sigma.sum());
    // Eq. (20) - Landing zone constraint
    socp.addConstraint(cvx::lessThan((cvx::par(E) * r.col(i_f) - cvx::par(trajectory_data.q)).norm(),
                                     cvx::par(target)));
    solver.solve(false);

    // Store to the socp_state struct
    socp_state.i_f = i_f;
    socp.getVariableValue("x", socp_state.x);
    socp.getVariableValue("sigma", socp_state.sigma);
    socp.getVariableValue("u", socp_state.u);
    socp.getVariableValue("z", socp_state.z);

    return true;
}

bool gfold::gfold::compute() {
    bool optimal;

    // TODO: This is shit lol
    for (double t_f = min_duration;
         t_f <= max_duration;
         t_f += duration_search_interval) {
        optimal = p3p4(t_f);
        if (optimal) {
            break;
        }
    }

    return optimal;
}

int gfold::gfold::get_time_steps() const {
    return socp_state.i_f + 1;
}

const Eigen::MatrixXd &gfold::gfold::get_state_vector() const {
    return socp_state.x;
}

const Eigen::MatrixXd &gfold::gfold::get_sigma() const {
    return socp_state.sigma;
}

const Eigen::MatrixXd &gfold::gfold::get_u() const {
    return socp_state.u;
}

const Eigen::MatrixXd &gfold::gfold::get_z() const {
    return socp_state.z;
}