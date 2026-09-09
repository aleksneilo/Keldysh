#include <algorithm>
#include <cstring>
#include "stdafx.h"
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <fstream>
#include <complex>
#include "SFS.h"
#include <thread>
#include <string>
#include <chrono>
#include <vector>
#include <atomic>
 #define pi 3.141592653589793
 #define icom (complex<double>(0, 1.))

 using namespace std;

 class AndersonMixer
 {
 public:
     AndersonMixer(int history_size, double beta)
         : history_size_(history_size), beta_(beta)
     {
     }

     // Возвращает true, если использовался Anderson.
     // false означает обычный Picard.
     bool update(
         complex<double>* x,
         const complex<double>* target,
         int size)
     {
         vector<complex<double>> current_x(size);
         vector<complex<double>> current_r(size);
          
         for (int i = 0; i < size; ++i)
         {
             current_x[i] = x[i];
             current_r[i] = target[i] - x[i];
         }

         x_history_.push_back(current_x);
         r_history_.push_back(current_r);

         if (static_cast<int>(x_history_.size()) > history_size_)
         {
             x_history_.erase(x_history_.begin());
             r_history_.erase(r_history_.begin());
         }

         const int count =
             static_cast<int>(x_history_.size());

         // На первой итерации истории ещё нет.
         if (count < 2)
         {
             picardStep(x, current_x, current_r);
             return false;
         }

         /*
          * Решаем задачу
          *
          *     min ||sum_j coefficient[j] r_j||,
          *     sum_j coefficient[j] = 1.
          *
          * Матрица Грама:
          *
          *     H_ij = Re sum_x conj(r_i(x)) r_j(x).
          */

         vector<vector<double>> system(
             count + 1,
             vector<double>(count + 1, 0.0));

         vector<double> right_side(count + 1, 0.0);
         vector<double> solution;

         double gram_scale = 0.0;

         for (int row = 0; row < count; ++row)
         {
             for (int column = 0; column < count; ++column)
             {
                 double scalar_product = 0.0;

                 for (int i = 0; i < size; ++i)
                 {
                     scalar_product += real(
                         conj(r_history_[row][i])
                         * r_history_[column][i]);
                 }

                 system[row][column] = scalar_product;
             }

             gram_scale = max(
                 gram_scale,
                 abs(system[row][row]));
         }

         if (gram_scale < 1.e-30)
         {
             picardStep(x, current_x, current_r);
             return false;
         }

         // Нормируем матрицу, чтобы малая невязка
         // не делала систему плохо обусловленной.
         for (int row = 0; row < count; ++row)
         {
             for (int column = 0; column < count; ++column)
                 system[row][column] /= gram_scale;

             // Небольшая регуляризация.
             system[row][row] += 1.e-10;

             // Условие sum coefficient = 1.
             system[row][count] = 1.0;
             system[count][row] = 1.0;
         }

         right_side[count] = 1.0;

         if (!solveSmallSystem(system, right_side, solution))
         {
             restartWithLastPoint();
             picardStep(x, current_x, current_r);
             return false;
         }

         double coefficient_sum_abs = 0.0;
         for (int j = 0; j < count; ++j)
             coefficient_sum_abs += abs(solution[j]);

         // Большие коэффициенты означают плохо
         // обусловленную историю.
         if (!isfinite(coefficient_sum_abs) ||
             coefficient_sum_abs > 20.0)
         {
             restartWithLastPoint();
             picardStep(x, current_x, current_r);
             return false;
         }

         vector<complex<double>> next_x(
             size, complex<double>(0.0, 0.0));

         for (int j = 0; j < count; ++j)
         {
             const double coefficient = solution[j];

             for (int i = 0; i < size; ++i)
             {
                 next_x[i] += coefficient
                     * (x_history_[j][i]
                         + beta_ * r_history_[j][i]);
             }
         }

         for (int i = 0; i < size; ++i)
         {
             if (!isfinite(real(next_x[i])) ||
                 !isfinite(imag(next_x[i])))
             {
                 restartWithLastPoint();
                 picardStep(x, current_x, current_r);
                 return false;
             }
         }

         for (int i = 0; i < size; ++i)
             x[i] = next_x[i];

         return true;
     }

 private:
     int history_size_;
     double beta_;

     vector<vector<complex<double>>> x_history_;
     vector<vector<complex<double>>> r_history_;

     void picardStep(
         complex<double>* x,
         const vector<complex<double>>& current_x,
         const vector<complex<double>>& current_r)
     {
         for (int i = 0;
             i < static_cast<int>(current_x.size());
             ++i)
         {
             x[i] = current_x[i] + beta_ * current_r[i];
         }
     }

     void restartWithLastPoint()
     {
         if (x_history_.empty())
             return;

         vector<complex<double>> last_x =
             x_history_.back();

         vector<complex<double>> last_r =
             r_history_.back();

         x_history_.clear();
         r_history_.clear();

         x_history_.push_back(last_x);
         r_history_.push_back(last_r);
     }

     static bool solveSmallSystem(
         vector<vector<double>> matrix,
         vector<double> right_side,
         vector<double>& solution)
     {
         const int size =
             static_cast<int>(right_side.size());

         for (int column = 0; column < size; ++column)
         {
             int pivot = column;

             for (int row = column + 1;
                 row < size; ++row)
             {
                 if (abs(matrix[row][column])
                 > abs(matrix[pivot][column]))
                 {
                     pivot = row;
                 }
             }

             if (abs(matrix[pivot][column]) < 1.e-14)
                 return false;

             swap(matrix[column], matrix[pivot]);
             swap(right_side[column], right_side[pivot]);

             const double diagonal =
                 matrix[column][column];

             for (int j = column; j < size; ++j)
                 matrix[column][j] /= diagonal;

             right_side[column] /= diagonal;

             for (int row = 0; row < size; ++row)
             {
                 if (row == column)
                     continue;

                 const double factor =
                     matrix[row][column];

                 for (int j = column; j < size; ++j)
                 {
                     matrix[row][j] -=
                         factor * matrix[column][j];
                 }

                 right_side[row] -=
                     factor * right_side[column];
             }
         }

         solution = right_side;
         return true;
     }
 };

 void ScalcWarm(
     complex<double>* Del,
     complex<double>** S1,
     vector<vector<complex<double>>>& G_cache,
     vector<int>& iterG_stat,
     int w_in,
     int w_fin,
     double* q,
     double I);


 vector<int> BuildSparseMatsubaraGrid(
     int frequency_count,
     int dense_count,
     int points_per_octave)
 {
     vector<int> grid;

     if (frequency_count <= 0)
     {
         return grid;
     }

     dense_count = std::max(1, dense_count);
     dense_count = std::min(
         dense_count,
         frequency_count
     );

     points_per_octave =
         std::max(1, points_per_octave);

     // Низкие частоты считаем полностью
     for (int iw = 0; iw < dense_count; ++iw)
     {
         grid.push_back(iw);
     }

     if (dense_count == frequency_count)
     {
         return grid;
     }

     /*
     Частоты:

         omega_iw = pi*T*(2*iw + 1).

     На высоких частотах строим сетку,
     равномерную по log(omega).
     */
     const double frequency_ratio =
         pow(
             2.0,
             1.0 /
             static_cast<double>(points_per_octave)
         );

     int current_iw = dense_count - 1;

     while (current_iw < frequency_count - 1)
     {
         const double current_frequency_index =
             2.0 * current_iw + 1.0;

         const double target_frequency_index =
             frequency_ratio
             * current_frequency_index;

         int next_iw =
             static_cast<int>(
                 ceil(
                     0.5
                     * (target_frequency_index - 1.0)
                 )
                 );

         next_iw = std::max(
             next_iw,
             current_iw + 1
         );

         next_iw = std::min(
             next_iw,
             frequency_count - 1
         );

         grid.push_back(next_iw);
         current_iw = next_iw;
     }

     return grid;
 }

 void InterpolateSparseMatsubaraS1(
     complex<double>** S1,
     const vector<int>& sparse_iw)
 {
     if (sparse_iw.size() < 2)
     {
         return;
     }

     for (size_t k = 0;
         k + 1 < sparse_iw.size();
         ++k)
     {
         const int iw_left = sparse_iw[k];
         const int iw_right = sparse_iw[k + 1];

         // Между точками нет пропущенных частот
         if (iw_right <= iw_left + 1)
         {
             continue;
         }

         const double w_left =
             pi * T * (2.0 * iw_left + 1.0);

         const double w_right =
             pi * T * (2.0 * iw_right + 1.0);

         /*
         После умножения на omega ведущая
         поправка обычно пропорциональна 1/omega^2.
         */
         const double x_left =
             1.0 / (w_left * w_left);

         const double x_right =
             1.0 / (w_right * w_right);

         for (int iw = iw_left + 1;
             iw < iw_right;
             ++iw)
         {
             const double w =
                 pi * T * (2.0 * iw + 1.0);

             const double x =
                 1.0 / (w * w);

             const double interpolation_parameter =
                 (x - x_left)
                 / (x_right - x_left);

             for (int i = 0; i < N; ++i)
             {
                 const complex<double> y_left =
                     w_left * S1[iw_left][i];

                 const complex<double> y_right =
                     w_right * S1[iw_right][i];

                 const complex<double> y =
                     (1.0 - interpolation_parameter)
                     * y_left
                     + interpolation_parameter
                     * y_right;

                 S1[iw][i] = y / w;
             }
         }
     }
 }

//////// SelfCons - function to calculate pair potential ////////////

//   INPUT Requires Global variables: Scales and Del0 for initial pair potential
//   OUTPUT *G - pointer to array with Normal Green Function G
//   OUTPUT *Del - pointer to array with Pair Potential Delta


///////////////////////////////////////parallel programming
void SelfConsParal(complex<double> *G, complex<double> *Del, int Initial, double* q, double I)
{
     complex<double> Delbuf,S2, *SS;
     int *wth, Ns=0,//N_S+N_F+N_S1+N_F1+N_S2-2,
     NS1=0;//N_S+N_F+N_S1/2;

     SS=new complex<double> [N]; 			//sum of S1 in each x
     int amount_of_threads=8;
     wth=new int [amount_of_threads];
        for (int i=0; i<amount_of_threads-1; i++)
            wth[i]=int((0.3*i/7. + 0.69*i/7.*i/7.)*w_obrez);
        wth[1]=1; wth[amount_of_threads-1]=w_obrez;
        for (int i=0; i<amount_of_threads; i++) cout<<wth[i]<<endl;

     complex<double>** S1 = new complex<double> *[w_obrez];
     for (int i = 0; i < w_obrez; i++)
	   S1[i] = new complex<double>[N];

     const int anderson_history = 6;   ///HERE MAY BE CHANGE
     const double anderson_beta = 1.6; ///FOR FASTER CONVERGENCE YOU NEED TO INCREASE BETA

     AndersonMixer anderson(
         anderson_history,
         anderson_beta);

     vector<complex<double>> Del_target(N);
     vector<complex<double>> residual(N);

     vector<vector<complex<double>>> G_cache(
         w_obrez,
         vector<complex<double>>(N)
     );

     vector<int> iterG_stat(w_obrez, 0);

     for (int iw = 0; iw < w_obrez; iw++)
     {
         for (int i = 0; i < N; i++)
         {
             G_cache[iw][i] = 1.0 + 0.1 * icom;
         }
     }

     int worker_count =
         static_cast<int>(
             std::thread::hardware_concurrency()
             );

     if (worker_count <= 0)
     {
         worker_count = 1;
     }

     worker_count = std::min(
         worker_count,
         w_obrez
     );

     cout << "Matsubara worker threads: "
         << worker_count
         << endl;

     const int sparse_dense_count = int(0.1*(w_obrez));     ///HERE MAY BE CHANGE
     const int sparse_points_per_octave = 4;///FOR FASTER CONVERGENCE YOU NEED TO INCREASE BETA

     const vector<int> sparse_iw =
         BuildSparseMatsubaraGrid(
             w_obrez,
             sparse_dense_count,
             sparse_points_per_octave
         );

     cout << "Full Matsubara grid: "
         << w_obrez
         << endl;

     cout << "Sparse Matsubara grid: "
         << sparse_iw.size()
         << endl;

     //int imax; 
     int w0 = 0, w1 = 1;//, w2=2, w3=3, w4=5, w5=7, w6=9, w7=12, w8=15, w9=18, w10=22, w11=26, w12=30, w13=35, w14=39, w15=44, w16=49, w17=w_obrez;//54,w18=w_obrez;//for T=0.5
	double dDelmax=1, w;
     // Set Initial Delta And G values
	//if (Initial == 0)
    	/*/for (int i = 0; i<N; i++)
        {
			Del[i]= get_type(i)*Del0;
            if ((i > N_S)) Del[i] = -0.;// *Hi[1] / abs(Hi[1]) * get_type(i) * Del0;//*exp(icom*Xi2*pi/2.);//*Hi[1]/abs(Hi[1])*Del0;//*exp(icom*Xi2);;//*exp(icom*0.01*pi/2.);+N_F+N_S1+N_F1
            //if (i > (N - N_N - 1)) Del[i] = -get_type(i) * Del0;// *exp(icom * Xi2 * pi);
        }//*/
	
    //for (int i=0; i<N; i++) cout<<fixed<<i<<"  "<<Del[i]<<"  "<<get_type(i)<<"  "<<get_ksi(i)<<"  "<<get_wm(i,1)<<"  "<<Roi[Layer(i)]<<"  "<<Rbi[Layer(i)]<<"  "<<Layer(i)<<endl;
    // Start self-consistent loop
   	iter=0;

    S2 = log(T) / pi / T;
    for (int iw = 0; iw < w_obrez; iw++)
    {
        w = pi * T * (2. * iw + 1.); //cual of S2
        S2 += 2. / w;
    }

while ((dDelmax > epsDel)&& (iter < 2000))//((dDelmax<0.1)||((iter<70)))&&(iter<3000))
{
    // Find delta from S1 and S2 and check mismatch with previous step
        //epsG=pow(10.,-6-3*iter/50.); if(iter>50)
        //epsG=1e-11;
     	for (int iw=0; iw<w_obrez; iw++) for (int j=0; j<N; j++)   S1[iw][j]=0.+0.*icom;
		for(int i=0;i<N;i++) SS[i]=0.+0.*icom;
		/*/S2 = log(T) / pi / T;
        for(int iw=0; iw<w_obrez; iw++)
		{	w=pi*T*(2.*iw+1.); //cual of S2
			S2+=2./w;
	  	}/*/

        std::fill(
            iterG_stat.begin(),
            iterG_stat.end(),
            0
        );

        /*
        next_job — это не номер частоты,
        а номер элемента массива sparse_iw.
        */
        std::atomic<int> next_job(0);

        const int job_count =
            static_cast<int>(sparse_iw.size());

        const int active_worker_count =
            std::min(worker_count, job_count);

        vector<thread> workers;
        workers.reserve(active_worker_count);

        for (int worker = 0;
            worker < active_worker_count;
            ++worker)
        {
            workers.emplace_back(
                [&]()
                {
                    while (true)
                    {
                        const int job =
                            next_job.fetch_add(
                                1,
                                std::memory_order_relaxed
                            );

                        if (job >= job_count)
                        {
                            break;
                        }

                        /*
                        Преобразуем номер задания
                        в настоящий индекс частоты.
                        */
                        const int iw = sparse_iw[job];

                        ScalcWarm(
                            Del,
                            S1,
                            G_cache,
                            iterG_stat,
                            iw,
                            iw + 1,
                            q,
                            I
                        );
                    }
                }
            );
        }

        for (thread& worker : workers)
        {
            worker.join();
        }

        /*/std::vector<std::thread> threads(w_obrez);//amount_of_threads);
        for (int i = 0; i < w_obrez; i++)//amount_of_threads - 1; i++)
        {   //w0=wth[i];w1=wth[i+1];
            threads[i] = std::thread(ScalcWarm, std::ref(Del), std::ref(S1), std::ref(G_cache), std::ref(iterG_stat), i, i + 1, std::ref(q), std::ref(I));
        }

        for (int i = 0; i < w_obrez; i++)//amount_of_threads - 1; i++)
            threads[i].join();/*/
        
      
        /*/std::vector<std::thread> threads(w_obrez);//amount_of_threads);
        for (int i = 0; i < w_obrez; i++)//amount_of_threads - 1; i++)
        {   //w0=wth[i];w1=wth[i+1];
            threads[i] = std::thread(Scalc, std::ref(Del), std::ref(S1), i, i + 1, std::ref(q), std::ref(I));
        }

        for (int i = 0; i < w_obrez; i++)//amount_of_threads - 1; i++)
            threads[i].join();//*/

        bool scalc_ok = true;
        int failed_iw = -1;

        for (int iw : sparse_iw)
        {
            if (iterG_stat[iw] < 0)
            {
                scalc_ok = false;
                failed_iw = iw;
                break;
            }
        }

        if (!scalc_ok)
        {
            cerr << "ScalcWarm failed for iw = "
                << failed_iw
                << endl;

            break;
        }

        // Восстанавливаем пропущенные частоты
        InterpolateSparseMatsubaraS1(
            S1,
            sparse_iw
        );

        // Суммируем полную восстановленную сетку
        for (int i = 0; i < N; ++i)
        {
            for (int iw = 0; iw < w_obrez; ++iw)
            {
                SS[i] += S1[iw][i];
            }
        }


        dDelmax = 0.0;

        for (int i = 0; i < N; ++i)
        {
            Del_target[i] = SS[i] / S2;
            residual[i] = Del_target[i] - Del[i];

            dDelmax = max(
                dDelmax,
                abs(residual[i]));
        }

        bool used_anderson = false;

        if (dDelmax > epsDel)
        {
            used_anderson = anderson.update(
                Del,
                Del_target.data(),
                N);
        }

        /*/piccard 
        dDelmax = 0; int imax;
        for (int i = 0; i < N; i++)
        {
            const double mixing = 1.7;

            complex<double> Del_target = SS[i] / S2;
            complex<double> residual = Del_target - Del[i];

            dDelmax = max(dDelmax, abs(residual));
            Del[i] += mixing * residual;
        }//*/

        
        /*/ my working
        dDelmax = 0; int imax;
	        for (int i=0; i<N; i++)
	        {
	            Delbuf= Del[i];
	            Del[i]= (SS[i]/S2+Delbuf*(alpha-1.))/alpha;
                if (dDelmax< abs(Del[i]-Delbuf))
	            {   dDelmax = abs(Del[i]-Delbuf);
	                imax=i;
	            }
	        }//*/

            cout << iter << " " << T << "  " << dDelmax << "  " << real(Del[0]) << "  " << real(Del[N  -1]) << endl;// +N_F + N_S1 + N_F1 + N_S2 / 2]) << "  " << real(Del[N - N_N - 1]) << endl;
	        iter++;   // try it until it converges
	        //for(int i=0;i<N;i++)    cout<<i<<"  "<<Del[i]<<endl;
	        //auto end = std::chrono::steady_clock::now();
            //auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
            //cout <<"t("<<iter<<")=" << elapsed_ms.count() << endl;
    }

delete [] SS;
for (int i = 0; i < w_obrez; i++)
    delete[] S1[i];
delete[] S1;
delete[] wth;

}

void ScalcWarm(
    complex<double>* Del,
    complex<double>** S1,
    vector<vector<complex<double>>>& G_cache,
    vector<int>& iterG_stat,
    int w_in,
    int w_fin,
    double* q,
    double I)
{
    const int max_iterG = 1000;

    vector<complex<double>> Fi(N);
    vector<complex<double>> Fi1(N);

    // Проверка комплексного числа на NaN и Inf
    auto is_finite_complex =
        [](const complex<double>& z) -> bool
        {
            return std::isfinite(real(z))
                && std::isfinite(imag(z));
        };

    for (int iw = w_in; iw < w_fin; ++iw)
    {
        const double ww =
            pi * T * (2.0 * iw + 1.0);

        // G для данной частоты хранится между
        // внешними итерациями по Delta
        complex<double>* G = G_cache[iw].data();

        bool converged = false;
        int total_iterG = 0;

        /*
        attempt == 0:
            используем сохранённый G_cache[iw].

        attempt == 1:
            если warm start не удался,
            сбрасываем G и повторяем расчёт.
        */
        for (int attempt = 0;
            attempt < 2 && !converged;
            ++attempt)
        {
            if (attempt == 1)
            {
                // Защитный сброс только этой частоты
                for (int i = 0; i < N; ++i)
                {
                    G[i] = 1.0 + 0.1 * icom;
                }
            }

            double dGmax = 1.0;
            int iterG_local = 0;

            while (dGmax > epsG &&
                iterG_local < max_iterG)
            {
                // Решаем уравнение для Fi
                Prog(
                    Fi.data(),
                    G,
                    Del,
                    ww,
                    q,
                    I
                );

                // В F-параметризации используется
                // симметрия Fi1 = conjugate(Fi)
                for (int i = 0; i < N; ++i)
                {
                    Fi1[i] = conj(Fi[i]);
                }

                // Обновляем нормальную функцию G
                Gcalc(
                    G,
                    &dGmax,
                    Fi.data(),
                    Fi1.data(),
                    Del,
                    ww,
                    q
                );

                ++iterG_local;

                // Если dGmax стал NaN или Inf,
                // немедленно прекращаем эту попытку
                if (!std::isfinite(dGmax))
                {
                    break;
                }
            }

            total_iterG += iterG_local;

            bool finite_solution =
                std::isfinite(dGmax);

            // Дополнительная проверка массивов
            for (int i = 0;
                i < N && finite_solution;
                ++i)
            {
                finite_solution =
                    is_finite_complex(G[i]) &&
                    is_finite_complex(Fi[i]) &&
                    is_finite_complex(Fi1[i]);
            }

            converged =
                finite_solution &&
                dGmax <= epsG;
        }

        if (!converged)
        {
            /*
            Обе попытки завершились неудачно.

            Отрицательное значение сообщает
            SelfConsParal об ошибке данной частоты.
            */
            iterG_stat[iw] = -1;

            // Оставляем кеш в безопасном состоянии
            for (int i = 0; i < N; ++i)
            {
                G[i] = 1.0 + 0.1 * icom;
                S1[iw][i] = 0.0;
            }

            continue;
        }

        iterG_stat[iw] = total_iterG;

        // Частота успешно посчитана:
        // теперь можно вычислять вклад в Delta
        for (int i = 0; i < N; ++i)
        {
            const complex<double> wm =
                get_wm(i, ww)
                + pi
                * get_ksi(i)
                * get_ksi(i)
                * G[i]
                * q[i]
                * q[i];

            S1[iw][i] =
                2.0
                * get_tc(i)
                * get_type(i)
                * real(
                    Fi1[i]
                    / sqrt(
                        wm * wm
                        + Fi1[i] * conj(Fi[i])
                    )
                );
        }
    }
}

void Scalc(complex<double>* Del, complex<double>** S1, int w_in, int w_fin, double* q, double I)  //working in each thread
{
	int iterG=0, w_count=w_fin-w_in;
    double ww, dGmax;
    complex<double> *G,*G1, *Fi, *Fi_old, *Fi1, *Fi1_old, Gbuf,wm;
    G=new complex<double>[N];
    G1 = new complex<double>[N];
    Fi=new complex<double>[N];
    Fi_old=new complex<double>[N];
    Fi1=new complex<double>[N];
    Fi1_old=new complex<double>[N];

	for(int iw=w_in; iw<w_fin; iw++)
    {
        for (int i=0; i<N; i++)	    G[i]=1.+0.1*icom;
		//cout<<"iw="<<iw<<"  ID="<<this_thread::get_id()<<endl;
        ww=pi*T*(2.*iw+1.);
        dGmax=1.;
        iterG = 0; ;// int kk = 0; double dDmax_last = 0;
			 // Iterative loop over GG. Solve Nonlinear Usadel equation here
            while ((dGmax > epsG)&&(iterG<1000))//&&(kk==0))
            {   // Calculation of Anomalous Green functions Fi(w), Fi(-w)
                //for (int i = 0; i < N; i++) G1[i] = -conj(G[i]);
                //for (int i = 0; i < N; i++) Fi_old[i] = Fi[i];
                Prog( Fi, G, Del, ww, q, I);
                //for (int i = 0; i < N; i++) Fi[i] = 0.1*Fi[i] + 0.9*Fi_old[i];
                for (int i = 0; i < N; i++) Fi1[i] = conj(Fi[i]);
				//Prog( Fi1, G1, Del, -ww, q, I);
                Gcalc(G, &dGmax, Fi, Fi1, Del, ww,q);
                iterG++;
                //if(iw==1) cout<<Fi1[N/2]<<endl;
                //if(iw==0) cout<<"w="<<((ww/pi/T-1)/2.)<<"  "<<iterG<<"  "<< dGmax<<endl;
			}
            //if(iw==0)
                //cout<<"w="<<((ww/pi/T-1)/2.)<<"  "<<iterG<<"  "<< real(Del[0])<<"  "<<dGmax<<endl;//"  arg="<<arg(G[N-1]*Fi[N-1]/get_wm(N-1,ww))/pi<<endl;
            // search summ of self-cons eq. for each point of the grid
            for (int i=0; i<N; i++)
			{
                //ww = pi * T * (2. * iw + 1.);// +2. * pi * get_ksi(i) * get_ksi(i) * q[i] * G[i] * q[i] / 2.;
                wm = get_wm(i, ww) +2. * pi * get_ksi(i) * get_ksi(i) * G[i] * q[i] * q[i] / 2.;
                S1[iw][i] += 2. * get_tc(i) * get_type(i) * real(Fi1[i] / sqrt(wm * wm + Fi1[i] * conj(Fi[i])));// G[i] * Fi[i] / wm);// +conj(G[i]) * Fi1[i] / ww);
                //fout <<fixed<< (w/pi/T-1)/2 <<"\t"<< i <<"\t" << real(F[i].A[0]) << "\t" << imag(F[i].A[0]) << "\t" << real(G[i]) << "\t" << imag(G[i]) << "\t" <<real(Del[i])<<endl;// log10(abs(1.-G[i]*G[i]))<<  endl;
                //if(iw==0) fout1 <<fixed<< (w/pi/T-1)/2 <<"\t"<< i <<"\t" << real(F[i].A[1]) << "\t" << imag(F[i].A[1]) <<"\t" << real(F[i].A[2]) << "\t" << imag(F[i].A[2]) <<"\t" << real(F[i].A[3]) << "\t" << imag(F[i].A[3]) <<endl;// log10(abs(1.-G[i]*G[i]))<<  endl;

			}//*/
        }
    delete [] G;
    delete [] G1;
    delete [] Fi;
    delete [] Fi1;
    delete[] Fi_old;
    delete[] Fi1_old;

}
//////// SelfConsZero - Pair Potentail in bulk material ////////////

double SelfConsZero()
{
	double Del1, dDel, w;
    double S1, S2;

    Del1=1.76;  // Initial Delta is BCS value at T = 0
    dDel=1;
    while (abs(dDel)>epsDel)
        {
        S1=0;
        S2=log(T)/pi/(T)/2.;
        for(int iw=0; iw<w_obrez; iw++)
            {
            w=pi*T*(2.*iw+1.);
            S1 = S1 + Del1/sqrt(w*w+Del1*Del1);   // Use known Green Functions Fi in bulk superconductor
            S2 +=1./w;
            }

        dDel=Del1 - S1/S2;
        Del1= Del1 - dDel;
        }
    return Del1;
}

	//				for(int j=0;j<N;j++)
	//				{
	//					cout<<"Fi_old["<<j<<"]="<<Fi_old[j]<<" "<<"Fi["<<j<<"]="<<Fi[j]<<" "<<"Fi1["<<j<<"]="<<Fi1_old[j]<<" "<<"T="<<T<<'\n';
	//					cout<<"iter="<<iter<<'\n';
		//				cout<<"Fi["<<j<<"]="<<Fi[j]<<'\n';
			//			cout<<"Fi1["<<j<<"]="<<Fi1_old[j]<<" ";
			//		}
