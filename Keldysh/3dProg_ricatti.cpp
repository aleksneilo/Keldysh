#include "SFS.h"
#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <vector>

#define pi 3.141592653589793
#define icom (complex<double>(0, 1.0))

using namespace std;


// Fast nonlinear solver for the scalar diffusive Riccati equations.
//
// gamma and gamma_tilde remain two independent unknown functions. At every
// grid point Newton's correction is therefore a two-component complex vector.
// The residual at point i depends only on i-1, i and i+1, so the Newton
// Jacobian is block-tridiagonal with 2x2 complex blocks and is solved by a
// block version of the same Thomas sweep used in 3dProg.cpp.

namespace
{
    const double RICATTI_TINY = 1.e-14;
    const double RICATTI_JACOBIAN_STEP = 1.e-7;
    const double RICATTI_MIN_LINE_STEP = 1.0 / 128.0;
    const int RICATTI_NEWTON_MAX_ITER = 40;

    struct Vec2
    {
        complex<double> gamma;
        complex<double> gamma_tilde;

        Vec2(complex<double> value_gamma = 0.0,
             complex<double> value_gamma_tilde = 0.0)
            : gamma(value_gamma), gamma_tilde(value_gamma_tilde)
        {
        }
    };

    struct Mat2
    {
        complex<double> a00, a01, a10, a11;

        Mat2(complex<double> value00 = 0.0,
             complex<double> value01 = 0.0,
             complex<double> value10 = 0.0,
             complex<double> value11 = 0.0)
            : a00(value00), a01(value01),
              a10(value10), a11(value11)
        {
        }
    };

    // One workspace per worker thread.  It is reused for all Matsubara
    // frequencies handled by that thread and removes repeated new/delete from
    // the hottest part of the self-consistency loop.
    struct NewtonWorkspace
    {
        vector<Mat2> lower, diagonal, upper, sweep_matrix;
        vector<Vec2> residual, right_side, correction, sweep_vector;
        vector<Vec2> trial_residual;
        vector<complex<double>> gamma_trial, gamma_tilde_trial;

        void ensure_size(int size)
        {
            if (static_cast<int>(lower.size()) == size)
                return;

            lower.resize(size);
            diagonal.resize(size);
            upper.resize(size);
            sweep_matrix.resize(size);
            residual.resize(size);
            right_side.resize(size);
            correction.resize(size);
            sweep_vector.resize(size);
            trial_residual.resize(size);
            gamma_trial.resize(size);
            gamma_tilde_trial.resize(size);
        }
    };

    NewtonWorkspace& thread_workspace()
    {
        thread_local NewtonWorkspace workspace;
        workspace.ensure_size(N);
        return workspace;
    }

    complex<double> safe_denominator(complex<double> value)
    {
        if (abs(value) < RICATTI_TINY)
            value += complex<double>(RICATTI_TINY, 0.0);
        return value;
    }

    bool finite_complex(complex<double> value)
    {
        return isfinite(real(value)) && isfinite(imag(value));
    }

    Vec2 operator-(const Vec2& left, const Vec2& right)
    {
        return Vec2(left.gamma - right.gamma,
                    left.gamma_tilde - right.gamma_tilde);
    }

    Vec2 operator-(const Vec2& value)
    {
        return Vec2(-value.gamma, -value.gamma_tilde);
    }

    Mat2 operator-(const Mat2& left, const Mat2& right)
    {
        return Mat2(left.a00 - right.a00, left.a01 - right.a01,
                    left.a10 - right.a10, left.a11 - right.a11);
    }

    Vec2 operator*(const Mat2& matrix, const Vec2& vector)
    {
        return Vec2(matrix.a00 * vector.gamma
                    + matrix.a01 * vector.gamma_tilde,
                    matrix.a10 * vector.gamma
                    + matrix.a11 * vector.gamma_tilde);
    }

    Mat2 operator*(const Mat2& left, const Mat2& right)
    {
        return Mat2(
            left.a00 * right.a00 + left.a01 * right.a10,
            left.a00 * right.a01 + left.a01 * right.a11,
            left.a10 * right.a00 + left.a11 * right.a10,
            left.a10 * right.a01 + left.a11 * right.a11);
    }

    void set_column(Mat2& matrix, int column, const Vec2& value)
    {
        if (column == 0)
        {
            matrix.a00 = value.gamma;
            matrix.a10 = value.gamma_tilde;
        }
        else
        {
            matrix.a01 = value.gamma;
            matrix.a11 = value.gamma_tilde;
        }
    }

    bool invert_matrix(Mat2 matrix, Mat2* inverse)
    {
        complex<double> determinant =
            matrix.a00 * matrix.a11 - matrix.a01 * matrix.a10;

        if (abs(determinant) < RICATTI_TINY)
        {
            const double scale = 1.0 + abs(matrix.a00) + abs(matrix.a11);
            matrix.a00 += RICATTI_TINY * scale;
            matrix.a11 += RICATTI_TINY * scale;
            determinant = matrix.a00 * matrix.a11
                        - matrix.a01 * matrix.a10;
        }

        if (abs(determinant) < RICATTI_TINY * RICATTI_TINY)
            return false;

        *inverse = Mat2(matrix.a11 / determinant,
                       -matrix.a01 / determinant,
                       -matrix.a10 / determinant,
                       matrix.a00 / determinant);
        return true;
    }

    complex<double> ricatti_G(complex<double> gamma,
                              complex<double> gamma_tilde)
    {
        const complex<double> product = gamma * gamma_tilde;
        return (1.0 - product) / safe_denominator(1.0 + product);
    }

    complex<double> choose_matsubara_root(complex<double> value)
    {
        complex<double> root = sqrt(value);
        if ((real(root) < 0.0) ||
            ((abs(real(root)) < RICATTI_TINY) && (imag(root) < 0.0)))
            root = -root;
        return root;
    }

    bool is_left_interface(int iN)
    {
        return (iN == N_S - 1) ||
               (iN == N_S + N_F - 1) ||
               (iN == N_S + N_F + N_S1 - 1) ||
               (iN == N_S + N_F + N_S1 + N_F1 - 1);
    }

    bool is_right_interface(int iN)
    {
        return (iN == N_S) ||
               (iN == N_S + N_F) ||
               (iN == N_S + N_F + N_S1) ||
               (iN == N_S + N_F + N_S1 + N_F1);
    }

    complex<double> base_frequency(int iN, complex<double> w,
                                   bool tilde_equation)
    {
        (void)tilde_equation;

        // gamma_tilde(w) denotes the particle-hole conjugate amplitude
        // gamma^*(-w). After conjugating the negative-frequency equation,
        // both Riccati equations contain the same positive-frequency
        // combination w + iH. Using conj(w + iH) here would count the
        // exchange-field reversal twice. In particular, for a homogeneous
        // real Delta it would incorrectly force gamma_tilde=conj(gamma),
        // make G real and spoil the positive/negative Matsubara sum.
        return get_wm(iN, w);
    }

    complex<double> modified_frequency(int iN, complex<double> w,
                                       double* q,
                                       complex<double> gamma,
                                       complex<double> gamma_tilde,
                                       bool tilde_equation)
    {
        const double D = 2.0 * pi * get_ksi(iN) * get_ksi(iN);
        const double abs_w = abs(w);
        const double sign_w = (abs_w > RICATTI_TINY)
            ? real(w) / abs_w : 1.0;
        const double momentum = (q != nullptr) ? q[iN] : 0.0;
        const complex<double> G = ricatti_G(gamma, gamma_tilde);

        return base_frequency(iN, w, tilde_equation)
             + sign_w * D * momentum * momentum * G / 2.0;
    }

    complex<double> phi_value(int iN,
                              complex<double>* gamma,
                              complex<double>* gamma_tilde,
                              complex<double> w, double* q,
                              bool tilde_equation)
    {
        const complex<double> omega = modified_frequency(
            iN, w, q, gamma[iN], gamma_tilde[iN], tilde_equation);
        const complex<double> denominator = safe_denominator(
            1.0 - gamma[iN] * gamma_tilde[iN]);
        const complex<double> amplitude = tilde_equation
            ? gamma_tilde[iN] : gamma[iN];

        return 2.0 * omega * amplitude / denominator;
    }

    complex<double> phi_scale(int iN,
                              complex<double>* gamma,
                              complex<double>* gamma_tilde,
                              complex<double> w, double* q,
                              bool tilde_equation)
    {
        const complex<double> omega = modified_frequency(
            iN, w, q, gamma[iN], gamma_tilde[iN], tilde_equation);
        return 2.0 * omega / safe_denominator(
            1.0 - gamma[iN] * gamma_tilde[iN]);
    }

    complex<double> fixed_gamma_from_phi(complex<double> phi,
                                         complex<double> phi_tilde,
                                         complex<double> omega)
    {
        const complex<double> root = choose_matsubara_root(
            omega * omega + phi * phi_tilde);
        return phi / safe_denominator(omega + root);
    }

    complex<double> interface_residual_component(
        int iN,
        complex<double>* gamma,
        complex<double>* gamma_tilde,
        complex<double> w, double* q,
        bool tilde_equation,
        bool left_interface)
    {
        const double h = get_h(iN);
        const int layer = Layer(iN);
        const complex<double> rb_over_rho = Rbi[layer] / Roi[layer];

        if (left_interface)
        {
            const complex<double> G_left =
                ricatti_G(gamma[iN], gamma_tilde[iN]);
            const complex<double> G_right =
                ricatti_G(gamma[iN + 1], gamma_tilde[iN + 1]);
            const complex<double> wm_left =
                base_frequency(iN, w, tilde_equation);
            const complex<double> wm_right =
                base_frequency(iN + 1, w, tilde_equation);

            const complex<double> a = -rb_over_rho;
            const complex<double> b = -(rb_over_rho
                + G_right / safe_denominator(G_left) * h);
            const complex<double> c = -G_right
                / safe_denominator(G_left) * h * wm_left
                / safe_denominator(wm_right);

            return a * phi_value(iN - 1, gamma, gamma_tilde,
                                 w, q, tilde_equation)
                 - b * phi_value(iN, gamma, gamma_tilde,
                                 w, q, tilde_equation)
                 + c * phi_value(iN + 1, gamma, gamma_tilde,
                                 w, q, tilde_equation);
        }

        const complex<double> G_left =
            ricatti_G(gamma[iN - 1], gamma_tilde[iN - 1]);
        const complex<double> G_right =
            ricatti_G(gamma[iN], gamma_tilde[iN]);
        const complex<double> wm_left =
            base_frequency(iN - 1, w, tilde_equation);
        const complex<double> wm_right =
            base_frequency(iN, w, tilde_equation);

        const complex<double> a = G_left / safe_denominator(G_right)
                                * h * wm_right
                                / safe_denominator(wm_left);
        const complex<double> b = rb_over_rho
                                + G_left / safe_denominator(G_right) * h;
        const complex<double> c = rb_over_rho;

        return a * phi_value(iN - 1, gamma, gamma_tilde,
                             w, q, tilde_equation)
             - b * phi_value(iN, gamma, gamma_tilde,
                             w, q, tilde_equation)
             + c * phi_value(iN + 1, gamma, gamma_tilde,
                             w, q, tilde_equation);
    }

    complex<double> interior_residual_component(
        int iN,
        complex<double>* gamma,
        complex<double>* gamma_tilde,
        complex<double>* Del,
        complex<double> w, double* q,
        bool tilde_equation)
    {
        const double h = get_h(iN);
        const double D = 2.0 * pi * get_ksi(iN) * get_ksi(iN);
        const complex<double> value = tilde_equation
            ? gamma_tilde[iN] : gamma[iN];
        const complex<double> value_left = tilde_equation
            ? gamma_tilde[iN - 1] : gamma[iN - 1];
        const complex<double> value_right = tilde_equation
            ? gamma_tilde[iN + 1] : gamma[iN + 1];
        const complex<double> partner = tilde_equation
            ? gamma[iN] : gamma_tilde[iN];
        const complex<double> normalization = safe_denominator(
            1.0 + gamma[iN] * gamma_tilde[iN]);
        const complex<double> derivative =
            (value_right - value_left) / (2.0 * h);
        const complex<double> second_difference =
            value_left - 2.0 * value + value_right;
        const complex<double> Delta_eff = get_type(iN) * Del[iN];
        const complex<double> Delta_source = tilde_equation
            ? conj(Delta_eff) : Delta_eff;
        const complex<double> Delta_partner = tilde_equation
            ? Delta_eff : conj(Delta_eff);
        const complex<double> omega = modified_frequency(
            iN, w, q, gamma[iN], gamma_tilde[iN], tilde_equation);

        // Multiplication by h^2/D keeps the residual and the numerical
        // Jacobian well scaled on very fine grids.
        return second_difference
             - h * h * 2.0 * partner * derivative * derivative
               / normalization
             - h * h / D * 2.0 * omega * value
             + h * h / D
               * (Delta_source - Delta_partner * value * value);
    }

    Vec2 node_residual(int iN,
                       complex<double>* gamma,
                       complex<double>* gamma_tilde,
                       complex<double>* Del,
                       complex<double> w, double* q)
    {
        if ((MODE == 0) && (iN == 0))
        {
            return Vec2(
                phi_value(1, gamma, gamma_tilde, w, q, false)
                    - phi_value(0, gamma, gamma_tilde, w, q, false),
                phi_value(1, gamma, gamma_tilde, w, q, true)
                    - phi_value(0, gamma, gamma_tilde, w, q, true));
        }

        if ((MODE == 0) && (iN == N - 1))
        {
            return Vec2(
                phi_value(N - 1, gamma, gamma_tilde, w, q, false)
                    - phi_value(N - 2, gamma, gamma_tilde, w, q, false),
                phi_value(N - 1, gamma, gamma_tilde, w, q, true)
                    - phi_value(N - 2, gamma, gamma_tilde, w, q, true));
        }

        if ((MODE == 1) && ((iN == 0) || (iN == N - 1)))
        {
            const complex<double> phi = (iN == 0)
                ? 1.46 * exp(icom * 0.0)
                : Del0 * exp(icom * Xi2 * pi);
            const complex<double> phi_tilde = conj(phi);
            const complex<double> omega = modified_frequency(
                iN, w, q, gamma[iN], gamma_tilde[iN], false);
            const complex<double> omega_tilde = modified_frequency(
                iN, w, q, gamma[iN], gamma_tilde[iN], true);

            return Vec2(
                gamma[iN]
                    - fixed_gamma_from_phi(phi, phi_tilde, omega),
                gamma_tilde[iN]
                    - fixed_gamma_from_phi(phi_tilde, phi, omega_tilde));
        }

        if (is_left_interface(iN))
        {
            return Vec2(
                interface_residual_component(iN, gamma, gamma_tilde,
                                             w, q, false, true),
                interface_residual_component(iN, gamma, gamma_tilde,
                                             w, q, true, true));
        }

        if (is_right_interface(iN))
        {
            return Vec2(
                interface_residual_component(iN, gamma, gamma_tilde,
                                             w, q, false, false),
                interface_residual_component(iN, gamma, gamma_tilde,
                                             w, q, true, false));
        }

        return Vec2(
            interior_residual_component(iN, gamma, gamma_tilde,
                                        Del, w, q, false),
            interior_residual_component(iN, gamma, gamma_tilde,
                                        Del, w, q, true));
    }

    double calculate_residual(Vec2* residual,
                              complex<double>* gamma,
                              complex<double>* gamma_tilde,
                              complex<double>* Del,
                              complex<double> w, double* q)
    {
        double maximum = 0.0;
        for (int i = 0; i < N; ++i)
        {
            residual[i] = node_residual(i, gamma, gamma_tilde, Del, w, q);
            maximum = max(maximum, abs(residual[i].gamma));
            maximum = max(maximum, abs(residual[i].gamma_tilde));
        }
        return maximum;
    }

    void calculate_jacobian(Mat2* lower, Mat2* diagonal, Mat2* upper,
                            Vec2* base_residual,
                            complex<double>* gamma,
                            complex<double>* gamma_tilde,
                            complex<double>* Del,
                            complex<double> w, double* q)
    {
        for (int i = 0; i < N; ++i)
        {
            lower[i] = Mat2();
            diagonal[i] = Mat2();
            upper[i] = Mat2();

            const int first_node = max(0, i - 1);
            const int last_node = min(N - 1, i + 1);

            for (int changed_node = first_node;
                 changed_node <= last_node; ++changed_node)
            {
                Mat2* block = (changed_node < i) ? &lower[i]
                            : (changed_node > i) ? &upper[i]
                            : &diagonal[i];

                for (int component = 0; component < 2; ++component)
                {
                    complex<double>* changed_value = (component == 0)
                        ? &gamma[changed_node]
                        : &gamma_tilde[changed_node];
                    const complex<double> old_value = *changed_value;
                    const double step = RICATTI_JACOBIAN_STEP
                                      * (1.0 + abs(old_value));

                    *changed_value += step;
                    const Vec2 changed_residual = node_residual(
                        i, gamma, gamma_tilde, Del, w, q);
                    *changed_value = old_value;

                    set_column(*block, component,
                        Vec2((changed_residual.gamma
                              - base_residual[i].gamma) / step,
                             (changed_residual.gamma_tilde
                              - base_residual[i].gamma_tilde) / step));
                }
            }
        }
    }

    bool block_sweep(Mat2* lower, Mat2* diagonal, Mat2* upper,
                     Vec2* right_side, Vec2* correction,
                     Mat2* sweep_matrix, Vec2* sweep_vector)
    {
        Mat2 inverse;
        if (!invert_matrix(diagonal[0], &inverse))
            return false;

        sweep_matrix[0] = inverse * upper[0];
        sweep_vector[0] = inverse * right_side[0];

        for (int i = 1; i < N; ++i)
        {
            const Mat2 effective_diagonal =
                diagonal[i] - lower[i] * sweep_matrix[i - 1];
            const Vec2 effective_right_side =
                right_side[i] - lower[i] * sweep_vector[i - 1];

            if (!invert_matrix(effective_diagonal, &inverse))
                return false;

            sweep_matrix[i] = inverse * upper[i];
            sweep_vector[i] = inverse * effective_right_side;
        }

        correction[N - 1] = sweep_vector[N - 1];
        for (int i = N - 2; i >= 0; --i)
        {
            correction[i] = sweep_vector[i]
                          - sweep_matrix[i] * correction[i + 1];
        }
        return true;
    }

    bool valid_initial_guess(complex<double>* gamma,
                             complex<double>* gamma_tilde)
    {
        for (int i = 0; i < N; ++i)
        {
            if (!finite_complex(gamma[i]) ||
                !finite_complex(gamma_tilde[i]) ||
                (abs(1.0 + gamma[i] * gamma_tilde[i]) < RICATTI_TINY))
                return false;
        }
        return true;
    }

    double relative_correction_norm(Vec2* correction,
                                    complex<double>* gamma,
                                    complex<double>* gamma_tilde)
    {
        double maximum = 0.0;
        for (int i = 0; i < N; ++i)
        {
            maximum = max(maximum,
                abs(correction[i].gamma) / (1.0 + abs(gamma[i])));
            maximum = max(maximum,
                abs(correction[i].gamma_tilde)
                    / (1.0 + abs(gamma_tilde[i])));
        }
        return maximum;
    }

    void initialize_bulk_guess(complex<double>* gamma,
                               complex<double>* gamma_tilde,
                               complex<double>* G,
                               complex<double>* Del,
                               complex<double> w, double* q)
    {
        for (int i = 0; i < N; ++i)
        {
            const complex<double> Delta_eff = get_type(i) * Del[i];
            const complex<double> Delta_tilde = conj(Delta_eff);
            const double D = 2.0 * pi * get_ksi(i) * get_ksi(i);
            const double abs_w = abs(w);
            const double sign_w = (abs_w > RICATTI_TINY)
                ? real(w) / abs_w : 1.0;
            const double momentum = (q != nullptr) ? q[i] : 0.0;
            const complex<double> omega = get_wm(i, w)
                + sign_w * D * momentum * momentum / 2.0;
            const complex<double> omega_tilde = get_wm(i, w)
                + sign_w * D * momentum * momentum / 2.0;
            const complex<double> pair = Delta_eff * Delta_tilde;
            const complex<double> root =
                choose_matsubara_root(omega * omega + pair);
            const complex<double> root_tilde =
                choose_matsubara_root(omega_tilde * omega_tilde + pair);

            gamma[i] = Delta_eff / safe_denominator(omega + root);
            gamma_tilde[i] = Delta_tilde
                           / safe_denominator(omega_tilde + root_tilde);
            G[i] = ricatti_G(gamma[i], gamma_tilde[i]);
        }
    }
}


//////// getABC_ricatti: legacy scalar Picard coefficients //////////
// Kept for compatibility and comparison with getABC. The fast Prog_ricatti
// uses full 2x2 Newton blocks because four scalars cannot represent coupling
// between independent gamma and gamma_tilde.

void getABC_ricatti(complex<double>* a, complex<double>* b,
                    complex<double>* c, complex<double>* d,
                    int iN,
                    complex<double>* gamma,
                    complex<double>* gamma_tilde,
                    complex<double>* Del,
                    complex<double> w,
                    double* q, double I,
                    bool tilde_equation)
{
    (void)I;

    if ((MODE == 0) && (iN == 0))
    {
        *a = 0.0;
        *b = phi_scale(0, gamma, gamma_tilde, w, q, tilde_equation);
        *c = phi_scale(1, gamma, gamma_tilde, w, q, tilde_equation);
        *d = 0.0;
        return;
    }

    if ((MODE == 0) && (iN == N - 1))
    {
        *a = phi_scale(N - 2, gamma, gamma_tilde, w, q, tilde_equation);
        *b = phi_scale(N - 1, gamma, gamma_tilde, w, q, tilde_equation);
        *c = 0.0;
        *d = 0.0;
        return;
    }

    if ((MODE == 1) && ((iN == 0) || (iN == N - 1)))
    {
        const complex<double> phi = (iN == 0)
            ? 1.46 * exp(icom * 0.0)
            : Del0 * exp(icom * Xi2 * pi);
        const complex<double> phi_tilde = conj(phi);
        const complex<double> omega = modified_frequency(
            iN, w, q, gamma[iN], gamma_tilde[iN], tilde_equation);
        const complex<double> boundary = tilde_equation
            ? fixed_gamma_from_phi(phi_tilde, phi, omega)
            : fixed_gamma_from_phi(phi, phi_tilde, omega);

        *a = 0.0;
        *b = -1.0;
        *c = 0.0;
        *d = boundary;
        return;
    }

    if (is_left_interface(iN) || is_right_interface(iN))
    {
        const bool left = is_left_interface(iN);
        const double h = get_h(iN);
        const int layer = Layer(iN);
        const complex<double> rb_over_rho = Rbi[layer] / Roi[layer];
        complex<double> a_phi, b_phi, c_phi;

        if (left)
        {
            const complex<double> G_left =
                ricatti_G(gamma[iN], gamma_tilde[iN]);
            const complex<double> G_right =
                ricatti_G(gamma[iN + 1], gamma_tilde[iN + 1]);
            const complex<double> wm_left =
                base_frequency(iN, w, tilde_equation);
            const complex<double> wm_right =
                base_frequency(iN + 1, w, tilde_equation);

            a_phi = -rb_over_rho;
            b_phi = -(rb_over_rho
                    + G_right / safe_denominator(G_left) * h);
            c_phi = -G_right / safe_denominator(G_left) * h
                  * wm_left / safe_denominator(wm_right);
        }
        else
        {
            const complex<double> G_left =
                ricatti_G(gamma[iN - 1], gamma_tilde[iN - 1]);
            const complex<double> G_right =
                ricatti_G(gamma[iN], gamma_tilde[iN]);
            const complex<double> wm_left =
                base_frequency(iN - 1, w, tilde_equation);
            const complex<double> wm_right =
                base_frequency(iN, w, tilde_equation);

            a_phi = G_left / safe_denominator(G_right) * h
                  * wm_right / safe_denominator(wm_left);
            b_phi = rb_over_rho
                  + G_left / safe_denominator(G_right) * h;
            c_phi = rb_over_rho;
        }

        *a = a_phi * phi_scale(iN - 1, gamma, gamma_tilde,
                               w, q, tilde_equation);
        *b = b_phi * phi_scale(iN, gamma, gamma_tilde,
                               w, q, tilde_equation);
        *c = c_phi * phi_scale(iN + 1, gamma, gamma_tilde,
                               w, q, tilde_equation);
        *d = 0.0;
        return;
    }

    const double h = get_h(iN);
    const double D = 2.0 * pi * get_ksi(iN) * get_ksi(iN);
    const complex<double> value = tilde_equation
        ? gamma_tilde[iN] : gamma[iN];
    const complex<double> partner = tilde_equation
        ? gamma[iN] : gamma_tilde[iN];
    const complex<double> derivative = tilde_equation
        ? (gamma_tilde[iN + 1] - gamma_tilde[iN - 1]) / (2.0 * h)
        : (gamma[iN + 1] - gamma[iN - 1]) / (2.0 * h);
    const complex<double> normalization = safe_denominator(
        1.0 + gamma[iN] * gamma_tilde[iN]);
    const complex<double> gradient_term =
        2.0 * partner * derivative * derivative / normalization;
    const complex<double> Delta_eff = get_type(iN) * Del[iN];
    const complex<double> Delta_source = tilde_equation
        ? conj(Delta_eff) : Delta_eff;
    const complex<double> Delta_partner = tilde_equation
        ? Delta_eff : conj(Delta_eff);
    const complex<double> omega = modified_frequency(
        iN, w, q, gamma[iN], gamma_tilde[iN], tilde_equation);

    *a = 1.0;
    *b = 2.0 + h * h / D
                  * (2.0 * omega + Delta_partner * value);
    *c = 1.0;
    *d = h * h * gradient_term - h * h / D * Delta_source;
}


//////// Prog_ricatti: compatibility overload (cold start) //////////

void Prog_ricatti(complex<double>* gamma,
                   complex<double>* gamma_tilde,
                   complex<double>* G,
                   complex<double>* Del,
                   complex<double> w,
                   double* q, double I)
{
    Prog_ricatti(gamma, gamma_tilde, G, Del, w, q, I, false);
}


//////// Prog_ricatti: damped block-Newton solver //////////
// use_initial_guess=false: initialize from the local bulk solution.
// use_initial_guess=true : reuse gamma/gamma_tilde supplied by the caller.

void Prog_ricatti(complex<double>* gamma,
                   complex<double>* gamma_tilde,
                   complex<double>* G,
                   complex<double>* Del,
                   complex<double> w,
                   double* q, double I,
                   bool use_initial_guess)
{
    (void)I;

    if (N <= 0)
        return;

    if (!use_initial_guess || !valid_initial_guess(gamma, gamma_tilde))
        initialize_bulk_guess(gamma, gamma_tilde, G, Del, w, q);

    NewtonWorkspace& workspace = thread_workspace();
    Mat2* lower = workspace.lower.data();
    Mat2* diagonal = workspace.diagonal.data();
    Mat2* upper = workspace.upper.data();
    Mat2* sweep_matrix = workspace.sweep_matrix.data();
    Vec2* residual = workspace.residual.data();
    Vec2* right_side = workspace.right_side.data();
    Vec2* correction = workspace.correction.data();
    Vec2* sweep_vector = workspace.sweep_vector.data();
    complex<double>* gamma_trial = workspace.gamma_trial.data();
    complex<double>* gamma_tilde_trial =
        workspace.gamma_tilde_trial.data();
    Vec2* trial_residual = workspace.trial_residual.data();

    const double tolerance = max(epsG, 1.e-11);
    double residual_norm = calculate_residual(
        residual, gamma, gamma_tilde, Del, w, q);
    int newton_iteration = 0;
    bool converged = false;

    // The interior equations are multiplied by h^2. Their raw residual can
    // therefore be tiny on a fine mesh even when the remaining change of the
    // solution is not. Always form the Newton correction and use its relative
    // norm as the actual convergence criterion.
    while (newton_iteration < RICATTI_NEWTON_MAX_ITER)
    {
        calculate_jacobian(lower, diagonal, upper, residual,
                           gamma, gamma_tilde, Del, w, q);

        for (int i = 0; i < N; ++i)
            right_side[i] = -residual[i];

        if (!block_sweep(lower, diagonal, upper, right_side,
                         correction, sweep_matrix, sweep_vector))
        {
            cerr << "Prog_ricatti: singular Newton matrix for w="
                 << w << endl;
            break;
        }

        const double correction_norm = relative_correction_norm(
            correction, gamma, gamma_tilde);
        if (!isfinite(correction_norm))
        {
            cerr << "Prog_ricatti: non-finite Newton correction for w="
                 << w << endl;
            break;
        }

        if (correction_norm <= tolerance)
        {
            converged = true;
            break;
        }

        double line_step = 1.0;
        double trial_norm = residual_norm;
        bool accepted = false;

        while (line_step >= RICATTI_MIN_LINE_STEP)
        {
            bool finite_trial = true;
            for (int i = 0; i < N; ++i)
            {
                gamma_trial[i] = gamma[i]
                    + line_step * correction[i].gamma;
                gamma_tilde_trial[i] = gamma_tilde[i]
                    + line_step * correction[i].gamma_tilde;
                finite_trial = finite_trial
                    && finite_complex(gamma_trial[i])
                    && finite_complex(gamma_tilde_trial[i]);
            }

            if (finite_trial)
            {
                trial_norm = calculate_residual(
                    trial_residual, gamma_trial, gamma_tilde_trial,
                    Del, w, q);
                accepted = trial_norm < residual_norm;
            }

            if (accepted)
                break;
            line_step *= 0.5;
        }

        if (!accepted)
        {
            cerr << "Prog_ricatti: Newton line search failed for w="
                 << w << ", residual=" << residual_norm << endl;
            break;
        }

        for (int i = 0; i < N; ++i)
        {
            gamma[i] = gamma_trial[i];
            gamma_tilde[i] = gamma_tilde_trial[i];
            residual[i] = trial_residual[i];
        }

        residual_norm = trial_norm;
        ++newton_iteration;
    }

    if (!converged &&
        (newton_iteration == RICATTI_NEWTON_MAX_ITER))
    {
        cerr << "Prog_ricatti: no Newton convergence for w=" << w
             << ", residual=" << residual_norm << endl;
    }

    double dGmax = 0.0;
    Gcalc_ricatti(G, &dGmax, gamma, gamma_tilde);

}


//////// Gcalc_ricatti: algebraic normalization //////////

void Gcalc_ricatti(complex<double>* G, double* dGmax,
                   complex<double>* gamma,
                   complex<double>* gamma_tilde)
{
    *dGmax = 0.0;
    for (int i = 0; i < N; ++i)
    {
        const complex<double> old_G = G[i];
        G[i] = ricatti_G(gamma[i], gamma_tilde[i]);
        *dGmax = max(*dGmax, abs(G[i] - old_G));
    }
}
