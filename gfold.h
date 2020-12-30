#ifndef GFOLD_GFOLD_H
#define GFOLD_GFOLD_H

#include <epigraph.hpp>

namespace gfold {
    /**
     * @brief Properties of the lander used to compute the
     * optimal trajectory.
     */
    struct gfold_lander_data {
        /**
         * The wet mass of the lander, kilograms.
         */
        double m_0;
        /**
         * The initial mass of the fuel on the lander
         * (such that m_0 = m_dry + m_f), kilograms.
         */
        double m_f;
        /**
         * The landing propulsion system's minimum throttle,
         * Newtons.
         */
        double rho_1;
        /**
         * The landing propulsion system's maximum throttle,
         * Newtons.
         */
        double rho_2;
        /**
         * The mass depletion constant (the amount of fuel
         * depleted dependent on the thrust produced),
         * seconds/meter.
         */
        double alpha;
        /**
         * The maximum velocity the lander is permitted to
         * experience, meters/second.
         */
        double V_max;
    };

    /**
     * @brief Properties of the computed optimal trajectory.
     */
    struct gfold_trajectory_data {
        /**
         * The initial position of the trajectory, meters.
         */
        Eigen::Vector3d r_0;
        /**
         * The initial velocity vector, meters/second.
         */
        Eigen::Vector3d r_dot_0;
        /**
         * The glideslope angle (measured from the ground to
         * the position of the lander), radians.
         */
        double gamma_gs;
        /**
         * The maximum deviation from the thrust pointing angle
         * (measured from the thrust unit vector), radians.
         */
        double theta;
        /**
         * The thrust unit vector.
         */
        Eigen::Vector3d n_hat;
        /**
         * The position of the desired landing site, meters.
         */
        Eigen::Vector2d q;
    };

    /**
     * @brief Properties of the landing surface.
     */
    struct gfold_planet_data {
        /**
         * The gravitational constant of the body with the
         * landing surface, meters/second^2.
         */
        Eigen::Vector3d g;
        /**
         * The angular velocity of the body's surface,
         * radians/second.
         */
        Eigen::Vector3d omega;
    };

    /**
     * @brief Wrapper class for the G-FOLD algorithm,
     * containing the input and output parameters and the state
     * of the conic solver.
     */
    class gfold {
    private:
        const struct gfold_lander_data lander_data;
        const struct gfold_trajectory_data trajectory_data;
        const struct gfold_planet_data planet_data;

        const double min_duration;
        const double max_duration;
        const double duration_search_interval;
        const double dt;

        struct gfold_state {
            int i_f;

            Eigen::MatrixXd x;
            Eigen::MatrixXd m;
            Eigen::MatrixXd sigma;
            Eigen::MatrixXd u;
            Eigen::MatrixXd z;
        } socp_state;

        /**
         * Performs the computation for Problem 3 and then
         * Problem 4 if the given duration is sufficient to
         * find an optimal solution.
         *
         * @param duration the desired duration
         * @return true if an optimal solution was found
         */
        bool p3p4(double duration);

    public:
        /**
         * Configures a new instance of the G-FOLD
         * algorithm with the given parameters. This
         * constructor will optimize for the time of
         * flight, which will be returned from the
         * algorithm through get_time_steps() in the
         * specified range.
         *
         * @param lander_data the lander parameters
         * @param trajectory_data the trajectory parameters
         * @param planet_data the target data
         * @param min_duration the minimum duration to
         * begin looking for an optimal solution
         * @param max_duration the maximum duration to
         * constrain the solution
         * @param duration_search_interval the interval
         * between consecutive flight durations which to
         * look for an optimal trajectory
         * @param dt the time step for the solution
         */
        gfold(gfold_lander_data lander_data,
              gfold_trajectory_data trajectory_data,
              gfold_planet_data planet_data,
              double min_duration,
              double max_duration,
              double duration_search_interval,
              double dt);

        /**
         * Configures a new instance of the G-FOLD
         * algorithm with the given parameters. This
         * constructor affixes the possible duration search
         * at the given duration. All members of the
         * parameter structs must be filled.
         *
         * @param lander_data the lander parameters
         * @param trajectory_data the trajectory parameters
         * @param planet_data the target data
         * @param duration the duration until landing
         * @param dt the time step for the solution
         */
        gfold(gfold_lander_data lander_data,
              gfold_trajectory_data trajectory_data,
              gfold_planet_data planet_data,
              double duration,
              double dt);

        /**
         * Computes the optimization problem and determines
         * if a feasible solution was found.
         *
         * @return true if a feasible solution was found
         */
        bool compute();

        /**
         * Obtains the number of time steps needed to find
         * an optimal trajectory from the compute() method.
         *
         * Not valid if compute() does not return true.
         *
         * @return the number of time steps for the optimal
         * trajectory
         */
        int get_time_steps() const;

        /**
         * Returns a 6xN matrix in which each column is the
         * state vector at t = j*dt. A state vector is an
         * X, Y, Z, X_dot, Y_dot and Z_dot or the position
         * and velocity in meters and meters/second,
         * respectively.
         *
         * N is given by get_time_steps().
         *
         * Not valid if compute() does not return true.
         *
         * @return the state vectors at each discrete time
         * instant
         */
        const Eigen::MatrixXd &get_state_vector() const;

        /**
         * Obtains the thrust slack variable containing the
         * slack divided by mass at t = j*dt.
         *
         * N is given by get_time_steps().
         *
         * Not valid if compute() does not return true.
         *
         * @return the thrust slack matrix
         */
        const Eigen::MatrixXd &get_sigma() const;

        /**
         * Obtains the 3xN thrust matrix containing the
         * thrust divided by mass value at t = j*dt.
         *
         * N is given by get_time_steps().
         *
         * Not valid if compute() does not return true.
         *
         * @return the thrust matrix
         */
        const Eigen::MatrixXd &get_u() const;

        /**
         * Obtains the 1xN mass matrix containing the
         * natrual logarithm of the lander's mass at
         * t = j*dt.
         *
         * N is given by get_time_steps().
         *
         * Not valid if compute() does not return true.
         *
         * @return the mass matrix
         */
        const Eigen::MatrixXd &get_z() const;
    };
}

#endif // GFOLD_GFOLD_H
