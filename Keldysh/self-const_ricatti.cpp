#include "stdafx.h"
#include "SFS.h"
#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <thread>
#include <vector>

#define pi 3.141592653589793

using namespace std;


// Fast self-consistency loop for the Riccati solver.
//
// The main acceleration relative to the first educational version is the
// Matsubara cache.  gamma(iw,x), gamma_tilde(iw,x) and G(iw,x) survive between
// Delta iterations, so every new spectral calculation starts from the
// converged solution of the previous Delta profile.

namespace
{
    const double SELFCONS_RICATTI_TINY = 1.e-14;
    const double SELFCONS_MIN_RELAXATION = 0.05;
    const double SELFCONS_MAX_RELAXATION = 1.5;
    const int SELFCONS_RICATTI_MAX_ITER = 2000;

    // Contiguous storage is considerably more cache-friendly than allocating
    // every Matsubara row separately.  row(iw) can still be used as an old
    // complex<double>* array, which keeps the code easy to compare with S1.
    class ComplexGrid
    {
    public:
        ComplexGrid(int rows, int columns)
            : columns_(columns), data_(rows * columns), row_(rows)
        {
            for (int i = 0; i < rows; ++i)
                row_[i] = data_.data() + i * columns_;
        }

        complex<double>* row(int index)
        {
            return row_[index];
        }

        complex<double>** rows()
        {
            return row_.data();
        }

        void clear()
        {
            fill(data_.begin(), data_.end(), complex<double>(0.0, 0.0));
        }

    private:
        int columns_;
        vector<complex<double>> data_;
        vector<complex<double>*> row_;
    };

    complex<double> safe_ricatti_denominator(complex<double> gamma,
                                              complex<double> gamma_tilde)
    {
        complex<double> denominator = 1.0 + gamma * gamma_tilde;
        if (abs(denominator) < SELFCONS_RICATTI_TINY)
            denominator += complex<double>(SELFCONS_RICATTI_TINY, 0.0);
        return denominator;
    }

    void solve_matsubara_frequency(
        complex<double>* Del,
        complex<double>* S1_row,
        int iw,
        double* q, double I,
        complex<double>* gamma,
        complex<double>* gamma_tilde,
        complex<double>* G,
        bool use_initial_guess)
    {
        const double w = pi * T * (2.0 * iw + 1.0);

        Prog_ricatti(gamma, gamma_tilde, G, Del,
                     w, q, I, use_initial_guess);

        for (int i = 0; i < N; ++i)
        {
            const complex<double> denominator =
                safe_ricatti_denominator(gamma[i], gamma_tilde[i]);
            const complex<double> F = 2.0 * gamma[i] / denominator;
            const complex<double> F_tilde =
                2.0 * gamma_tilde[i] / denominator;

            // For the old real symmetric problem this is exactly
            // 2*real(F_tilde).  Written in this form it also preserves the
            // complex phase when gamma and gamma_tilde are independent:
            // F(-w) = conj(F_tilde(w)).
            S1_row[i] = get_tc(i) * get_type(i)
                      * (F + conj(F_tilde));
        }
    }
}


//////// SelfConsParal_ricatti: cached parallel Delta calculation //////////

void SelfConsParal_ricatti(complex<double>* G,
                           complex<double>* Del,
                           int Initial,
                           double* q, double I)
{
    if ((N <= 0) || (w_obrez <= 0) || (T <= 0.0))
    {
        cerr << "SelfConsParal_ricatti: N, w_obrez and T must be positive"
             << endl;
        return;
    }

    vector<complex<double>> SS(N);
    vector<complex<double>> gap_residual(N);
    vector<complex<double>> previous_gap_residual(N);
    ComplexGrid S1(w_obrez, N);
    ComplexGrid gamma_cache(w_obrez, N);
    ComplexGrid gamma_tilde_cache(w_obrez, N);
    ComplexGrid G_cache(w_obrez, N);

    if (Initial == 0)
    {
        for (int i = 0; i < N; ++i)
        {
            Del[i] = get_type(i) * Del0;
            if (i > N_S)
                Del[i] = 0.0;
        }
    }

    // S2 depends only on T and the frequency cutoff, not on Delta.  Calculate
    // it once instead of repeating it in every self-consistency iteration.
    double S2 = log(T) / (pi * T);
    for (int iw = 0; iw < w_obrez; ++iw)
    {
        const double w = pi * T * (2.0 * iw + 1.0);
        S2 += 2.0 / w;
    }

    if (abs(S2) < SELFCONS_RICATTI_TINY)
    {
        cerr << "SelfConsParal_ricatti: singular self-consistency sum" << endl;
        return;
    }

    // Do not create one operating-system thread per frequency.  A bounded
    // number avoids oversubscription and improves cache locality.
    const unsigned int hardware_threads = thread::hardware_concurrency();
    const int worker_count = min(
        w_obrez, max(1, static_cast<int>(hardware_threads)));

    double dDelmax = 1.0;
    double update_max = 1.0;
    double relaxation = 1.0;
    if (isfinite(alpha) && (alpha > SELFCONS_RICATTI_TINY))
        relaxation = min(1.0, 1.0 / alpha);
    bool have_previous_gap_residual = false;
    bool cache_is_valid = false;
    iter = 0;

    while ((dDelmax > epsDel) &&
           (iter < SELFCONS_RICATTI_MAX_ITER))
    {
        S1.clear();
        fill(SS.begin(), SS.end(), complex<double>(0.0, 0.0));

        vector<thread> workers;
        workers.reserve(worker_count);

        // Strided frequency assignment distributes the expensive low
        // Matsubara frequencies among different workers.
        for (int worker = 0; worker < worker_count; ++worker)
        {
            workers.emplace_back([&, worker, cache_is_valid]()
            {
                for (int iw = worker; iw < w_obrez; iw += worker_count)
                {
                    solve_matsubara_frequency(
                        Del, S1.row(iw), iw, q, I,
                        gamma_cache.row(iw),
                        gamma_tilde_cache.row(iw),
                        G_cache.row(iw),
                        cache_is_valid);
                }
            });
        }

        for (thread& worker : workers)
            worker.join();

        cache_is_valid = true;

        // Traverse each Matsubara row contiguously.
        for (int iw = 0; iw < w_obrez; ++iw)
            for (int i = 0; i < N; ++i)
                SS[i] += S1.row(iw)[i];

        dDelmax = 0.0;
        bool finite_gap_map = true;
        for (int i = 0; i < N; ++i)
        {
            gap_residual[i] = SS[i] / S2 - Del[i];
            finite_gap_map = finite_gap_map
                && isfinite(real(gap_residual[i]))
                && isfinite(imag(gap_residual[i]));
            dDelmax = max(dDelmax, abs(gap_residual[i]));
        }

        if (!finite_gap_map)
        {
            cerr << "SelfConsParal_ricatti: non-finite Delta map" << endl;
            return;
        }

        // Vector Aitken relaxation. A change of sign/direction in the gap
        // residual automatically reduces the step, while a smooth fixed-point
        // trajectory is allowed to use the full step. The old expression used
        // 1/alpha=1.695 for alpha=0.59 and amplified an oscillatory mode.
        if (have_previous_gap_residual)
        {
            double numerator = 0.0;
            double denominator = 0.0;
            for (int i = 0; i < N; ++i)
            {
                const complex<double> difference =
                    gap_residual[i] - previous_gap_residual[i];
                numerator += real(conj(previous_gap_residual[i])
                                * difference);
                denominator += norm(difference);
            }

            if (denominator > SELFCONS_RICATTI_TINY)
            {
                const double aitken =
                    -relaxation * numerator / denominator;
                if (isfinite(aitken))
                    relaxation = max(SELFCONS_MIN_RELAXATION,
                        min(SELFCONS_MAX_RELAXATION, aitken));
            }
        }

        update_max = 0.0;
        for (int i = 0; i < N; ++i)
        {
            const complex<double> update = relaxation * gap_residual[i];
            Del[i] += update;
            update_max = max(update_max, abs(update));
            previous_gap_residual[i] = gap_residual[i];
        }
        have_previous_gap_residual = true;

        cout << "[Riccati-fast] " << iter << " " << T
             << "  residual=" << dDelmax
             << "  update=" << update_max
             << "  mix=" << relaxation
             << "  " << real(Del[0])
             << "  " << real(Del[N - 1]) << endl;

        ++iter;
    }

    if (iter == SELFCONS_RICATTI_MAX_ITER)
    {
        cerr << "SelfConsParal_ricatti: Delta did not converge, dDelta="
             << dDelmax << endl;
    }

    // Refresh the first cached frequency for the final Delta and expose its G
    // through the same output argument as the old function.
    solve_matsubara_frequency(
        Del, S1.row(0), 0, q, I,
        gamma_cache.row(0), gamma_tilde_cache.row(0), G_cache.row(0),
        cache_is_valid);

    for (int i = 0; i < N; ++i)
        G[i] = G_cache.row(0)[i];
}


//////// Scalc_ricatti: standalone frequency-range interface //////////
//
// SelfConsParal_ricatti uses the persistent cache above.  This public function
// keeps the old Scalc-like signature for direct calls.  Within a range it uses
// the preceding frequency as the initial guess for the next one.

void Scalc_ricatti(complex<double>* Del,
                   complex<double>** S1,
                   int w_in, int w_fin,
                   double* q, double I)
{
    vector<complex<double>> gamma(N);
    vector<complex<double>> gamma_tilde(N);
    vector<complex<double>> G(N);
    bool have_initial_guess = false;

    for (int iw = w_in; iw < w_fin; ++iw)
    {
        solve_matsubara_frequency(
            Del, S1[iw], iw, q, I,
            gamma.data(), gamma_tilde.data(), G.data(),
            have_initial_guess);
        have_initial_guess = true;
    }
}
