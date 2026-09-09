# SNS/Floquet в нормировке исходной программы

Эта версия полностью заменяет предыдущий интерфейс в SI. Входные числа больше
не означают джоули, метры или вольты. Спектральные, кинетические уравнения и ток
используют одну безразмерную систему единиц. Старые вызовы с SI надо переписать.

## Основной запуск через существующий Keldysh.sln

Откройте Keldysh.sln и нажмите Ctrl+F5. Единственный main общего приложения
находится в Keldysh/main.cpp. В его начале задаются settings.physical,
settings.numerical, settings.mode, settings.epsilon и settings.voltage.
По умолчанию mode="spectral", epsilon=3*Delta, voltage=0.
Результат выводится в консоль. Для своей длины замените settings.physical.L_N=sqrt(...).
Для точки ВАХ задайте mode="iv", voltage=2*settings.physical.Delta.
Для трёх точек v/Delta=2,1,0.5 задайте mode="benchmark".

use_sns=false включает старый расчёт run_legacy_calculation().
То же делает аргумент --legacy. Старые глобальные параметры используются только
старой веткой; настройки SNS задаются явно в начале main.
Аргументы spectral/iv/normal/benchmark/convergence переопределяют настройки.
Очистите аргументы отладки Visual Studio, если хотите использовать только main.

sns/runner.cpp — общий диспетчер режимов без собственного main.
sns/cli.cpp — необязательная обёртка отдельного sns_cli.exe, не включённая
в Keldysh.vcxproj. tests.cpp — отдельная тестовая программа.
CMake и build-sns.ps1 остаются дополнительными средствами; для запуска из старого
решения Visual Studio они не требуются.

## Вывод G и F после каждой нелинейной итерации

В основном main.cpp включён settings.numerical.iteration_log_path.
По умолчанию файл: C:/Users/user/source/repos/aleksneilo/keldysh/Ricatti/GF_iterations.txt.
Чтобы выключить запись, присвойте этому параметру пустую строку.
При запуске через Keldysh.exe файл очищается один раз. Все спектральные задачи
внутри этого запуска дописываются в него, включая продолжение по энергии/напряжению.
Прямой вызов solve_gamma_for_energy только дописывает; очисткой управляет вызывающий код.

В spectral.cpp функция write_GF_iteration вызывается после принятия и смешивания
новой пары gamma,tilde_gamma и до проверки завершения итераций. Отвергнутые попытки
line search не записываются. Все lambda-этапы записываются; для итогового решения lambda=1.

G — нормальный электронный блок G^R, F — аномальный электронно-дырочный блок:
G=-i*pi*(I+gamma*tilde_gamma)^(-1)*(I-gamma*tilde_gamma),
F=-2*i*pi*(I+gamma*tilde_gamma)^(-1)*gamma.
Это полные матрицы Floquet, а не только элемент 00, и не функции распределения.

Столбцы (табуляция):
iteration lambda epsilon v x n m ReG ImG ReF ImF change residual

На каждой итерации записывается Nx*K*K строк. При V=0 K=1, n=m=0.
При V!=0 записаны все n,m от -NF до NF, включая внедиагональные элементы.
Строки с # содержат начало/конец спектральной задачи и её параметры.
Если расчёт прервался, принятые итерации остаются в файле, но маркера END converged нет.
Для отрицательного терминального напряжения solver использует симметрию: в журнале
сохраняется положительное v фактически решаемой спектральной задачи.

## 1. Нормировка установлена по исходному коду

В self-const_ricatti.cpp записано w=pi*T*(2*iw+1), в SelfConsZero() старт
Del1=1.76, а в 3dProg_ricatti.cpp коэффициент D=2*pi*get_ksi(i)*get_ksi(i).
В main.cpp Ksi_S=1 служит единицей длины, в граничных условиях фаза входит
как exp(i*pi*Xi2). Эти выражения фиксируют следующую нормировку:

| Величина в API | Определение |
|---|---|
| T | T_SI / Tc |
| Delta, epsilon, E_n, eta | соответствующая энергия / E0, E0=k_B Tc |
| voltage, v | e V_SI / E0 |
| L_N, координата x, шаг h | длина / xi_S |
| Ksi_N | sqrt(hbar D_N_SI / (2 pi E0)) / xi_S |
| ro_N | rho_N_SI / rho_S |
| area | A_SI / xi_S² |
| Xi | (phi_R-phi_L)/pi = Xi2-Xi1 |
| current, I_star | e rho_S I_SI / (E0 xi_S) |

Здесь xi_S=sqrt(hbar D_S_SI/(2 pi E0)). Поэтому
D= hbar D_N_SI/(E0 xi_S²)=2 pi Ksi_N².
Число Ksi_N не равно sqrt(D) и не равно длине sqrt(D/Delta):
между этими величинами есть множители 2 pi и Delta.

Поле area — явно определённое расширение для продольного SNS-тока.
Старый параметр I в KinInd относится к другой геометрии (поперечному току,
интегрируемому по слоям); автоматически отождествлять его с current нельзя.
При area=1 current численно равен нормированной плотности продольного тока.

Щель не нормируется на саму себя: например, при T=0 можно задать Delta=1.76,
а при иной температуре передать Del0 из существующего расчёта.
Новое ядро само Delta(T) не вычисляет.

## 2. Что задаётся перед спектральным решением

В sns.hpp:

~~~cpp
struct PhysicalParams {
    double Delta = 1.76, T = 0;
    double Ksi_N = 1, L_N = 1, ro_N = 1, area = 1, Xi = 0;

    double diffusion() const {
        return 2*std::acos(-1.)*Ksi_N*Ksi_N;
    }
    double thouless() const {
        return diffusion()/(L_N*L_N);
    }
    double conductance() const {
        return area/(ro_N*L_N);
    }
    double phase_radians() const {
        return std::acos(-1.)*Xi;
    }
};
~~~

thouless() оставлен только как производный масштаб. Solver использует diffusion()
и координату x от 0 до L_N, не координату x/L_N.

Пример самостоятельного вызова:

~~~cpp
sns::PhysicalParams p;
p.Delta = 1.76;  // Delta/(k_B Tc)
p.T = 0.1;      // T/Tc; Delta задана независимо
p.Ksi_N = 10;
p.L_N = 1;
p.ro_N = 0.1;
p.area = 1;
p.Xi = 0;       // фаза / pi

sns::NumericalParams num;
num.NF = 2;
num.Nx = 25;
num.eta = 0.001*p.Delta;
num.mixing = 0.15;
num.tolerance = 1e-8;
num.residual_tolerance = 1e-7;
num.max_iterations = 4000;
num.continuation_steps = 4;

double v = 2*p.Delta;          // eV/(k_B Tc), НЕ вольты
double epsilon = 0.37*p.Delta;
auto s = sns::solve_gamma_for_energy(epsilon, v, p, num);
~~~

Параметры Ksi_N=10,L_N=1 соответствуют короткому проводу относительно
xi_Delta=sqrt(D/Delta). Для benchmark L/xi_Delta=1 нужно явно выбрать
p.L_N=sqrt(p.diffusion()/p.Delta); это другой набор параметров.

Для интеграции с исходным main:

~~~cpp
sns::PhysicalParams p;
p.Delta = Del0;
p.T = T;
p.Ksi_N = Ksi_N;
p.L_N = L_N;
p.ro_N = ro_N;
p.Xi = Xi2-Xi1;
p.area = 1;  // задайте A/xi_S^2 для полного тока
~~~

Левая фаза выбрана нулевой общей gauge-трансформацией.
p.Xi=0.5 означает фазовую разность pi/2.
При Delta_N=0 спектру нужны Delta,Ksi_N,L_N,Xi,v,epsilon,eta.
T нужен кинетике, ro_N и area — току.

## 3. Неизвестные и Floquet

При v>0 выбирается epsilon в [0,2v).
E_n=epsilon+2nv, n=-NF,...,NF; K=2NF+1.

~~~cpp
int k = v == 0 ? 1 : 2*nf+1;
Matrix e(k,k);
for(int j=0;j<k;++j)
    e(j,j)=eps+2*(j-(v==0?0:nf))*v;
~~~

В памяти n переводится в n+NF. Первая строка соответствует n=-NF.
gamma_nm связывает две энергии E_n,E_m. Convolution имеет вид
(AB)_nm=sum_l A_nl B_lm. Никаких циклических переносов через край базиса нет.

~~~cpp
using Complex = std::complex<double>;
using Field = std::vector<Matrix>;
struct PairField { Field gamma, tilde; };
~~~

Оба поля имеют Nx матриц K×K. Доступ gamma[i](n+NF,m+NF).
Всего хранится 2 Nx K² complex, внутренних неизвестных 2(Nx-2)K².
Matrix использует column-major data[col*rows+row].

При v=0 K автоматически равен 1 и epsilon означает обычную энергию.
Это стационарная BVP. При v>0 ищется периодическое нестационарное решение в Floquet,
а не временная эволюция после включения напряжения.

## 4. Уравнение gamma в новой нормировке

Начинаем с размерного уравнения в N:

~~~math
\gamma_{x_{\rm SI}x_{\rm SI}}+
\gamma_{x_{\rm SI}}Q\gamma_{x_{\rm SI}}
=\frac{E_{\rm SI}\gamma+\gamma E_{\rm SI}}{i\hbar D_{\rm SI}}.
~~~

Подставляем x_SI=xi_S x, E_SI=E0 E, D=hbar D_SI/(E0 xi_S²).
Получаем уравнение, используемое в коде:

~~~math
\gamma''+\gamma'Q\gamma'
=\frac{E\gamma+\gamma E+2i\eta\gamma}{iD},
\qquad D=2\pi\,{\rm Ksi}_N^2,
\qquad x\in[0,L_N],
\qquad Q=-2(I+\widetilde\gamma\gamma)^{-1}\widetilde\gamma.
~~~

Штрих теперь означает производную по x=x_SI/xi_S.
Коэффициент D/2 из theta-параметризации сюда не переносится:
в скалярном стационарном случае числитель сам даёт 2E gamma.

eta>0 — численная регуляризация E→E+i eta, не физическая inelastic rate.
При конечном eta точное сохранение тока может нарушаться, поэтому надо проверять eta→0.

Второе поле удовлетворяет такому же уравнению с обменом gamma и tilde_gamma.
Знак энергии и знак i при tilde меняются одновременно, коэффициент сохраняется.

## 5. Tilde и границы

Tilde означает tilde_gamma(E,E')=conj(gamma(-E,-E')).
Нельзя брать conj(gamma) при прежней энергии.
В нефолдированной записи tilde_gamma_nm(eps)=conj(gamma_{-n,-m}(-eps)).
Внутри [0,2v), для eps>0, отражённый сектор равен 2v-eps,
а индексы -n-1,-m-1. При eps=0 используются -n,-m.

Конечный диапазон индексов не сохраняется при -n-1.
Поэтому основной solver совместно решает gamma и tilde_gamma;
функция tilde_gamma с callback для произвольных энергий используется для проверок.

BCS-функция:

~~~math
\gamma_0(E)=-\frac{\Delta}{w+i\sqrt{\Delta^2-w^2}},\qquad w=E+i\eta.
~~~

~~~cpp
Complex z(e,eta);
Complex root=std::sqrt(Complex(delta*delta)-z*z);
return -delta/(z+Complex(0,1)*root);
~~~

Главная ветвь корня даёт retarded решение при обоих знаках E.
При Delta=0 функция сразу возвращает ноль.

Прозрачные границы при v>0:

~~~math
\Gamma_{L,nm}=\gamma_0(E_n)\delta_{nm},
\quad
\Gamma_{R,nm}=e^{i\pi\mathrm{Xi}}\gamma_0(E_n+v)\delta_{n,m-1},
\quad
\widetilde\Gamma_{L,nm}=\gamma_0(-E_n)^*\delta_{nm},
\quad
\widetilde\Gamma_{R,nm}=e^{-i\pi\mathrm{Xi}}\gamma_0(-E_n+v)^*\delta_{n,m+1}.
~~~

Правая gamma лежит на наддиагонали, tilde_gamma на поддиагонали.
При v=0 обе границы диагональны, фаза сохраняется.

## 6. Начальное приближение

Отбрасываем нелинейный член. Для каждой пары nm:

~~~math
\gamma_{nm}''-\kappa_{nm}^2\gamma_{nm}=0,
\qquad
\kappa_{nm}^2=\frac{E_n+E_m+2i\eta}{iD}.
~~~

С заданными концами:

~~~math
\gamma_{nm}^{(0)}(x)=
\Gamma_{L,nm}\frac{\sinh[\kappa_{nm}(L_N-x)]}{\sinh(\kappa_{nm}L_N)}
+\Gamma_{R,nm}\frac{\sinh(\kappa_{nm}x)}{\sinh(\kappa_{nm}L_N)}.
~~~

В initial_gamma_linearized аргумент a=kappa L_N, доля z=i/(Nx-1):

~~~cpp
Complex a=p.L_N*std::sqrt(
    (e(j,j)+e(m,m)+Complex(0,2*n.eta))/Complex(0,p.diffusion()));
Complex l=sinh_ratio(a,1-z),r=sinh_ratio(a,z);
f.gamma[i](j,m)=l*b.left(j,m)+r*b.right(j,m);
~~~

z здесь только доля интервала для аналитического старта. Производные в solver
считаются по x, с шагом L_N/(Nx-1).

При |a|<1e-5 используется ряд sinh(as)/sinh(a)=s[1+a²(s²-1)/6].
Он включает линейную интерполяцию. При большой Re a используются затухающие
экспоненты для защиты от переполнения. Второе поле стартует аналогично.

## 7. Picard и прогонка

Пространственная сетка x_i=i h, h=L_N/(Nx-1).
На итерации r известны оба старых поля. Сначала:

~~~math
Q_i^{(r)}=-2(I+\widetilde\gamma_i^{(r)}\gamma_i^{(r)})^{-1}
\widetilde\gamma_i^{(r)},\quad
P_i^{(r)}=\frac{\gamma_{i+1}^{(r)}-\gamma_{i-1}^{(r)}}{2h},\quad
C_i^{(r)}=P_i^{(r)}Q_i^{(r)}P_i^{(r)}.
~~~

compute_Q решает (I+tilde_gamma gamma)Y=tilde_gamma с выбором главного элемента,
затем возвращает -2Y. Это матричное решение, не поэлементное деление.

Источник замораживается, затем для каждого nm решается:

~~~math
h^{-2}\gamma^*_{i-1}
+\left[-2h^{-2}-\frac{E_n+E_m+2i\eta}{iD}\right]\gamma^*_i
+h^{-2}\gamma^*_{i+1}=-\lambda C_i^{(r)}.
~~~

~~~cpp
double h=p.L_N/(n.Nx-1),a=1/(h*h);
int q=n.Nx-2;
Complex c=(e(j,j)+e(m,m)+Complex(0,2*n.eta))
          /Complex(0,p.diffusion());
std::vector<Complex> lo(q,a),diag(q,-2*a-c),up(q,a),rhs(q);
for(int i=0;i<q;++i)rhs[i]=-lambda*source[i+1](j,m);
rhs.front()-=a*left(j,m);
rhs.back()-=a*right(j,m);
lo.front()=0.;up.back()=0.;
auto sol=solve_tridiagonal_complex(lo,diag,up,rhs);
~~~

Nx-2 внутренних неизвестных образуют комплексную трёхдиагональную систему.
Thomas выполняет прямое исключение и обратную подстановку.
Отдельные nm можно решать по очереди только потому, что источник уже заморожен;
сами источники связаны полными матричными произведениями.

gamma_picard_step делает два frozen_step с обменом полей.
Оба используют значения предыдущей итерации. Затем:

~~~math
\gamma^{r+1}=(1-\alpha)\gamma^r+\alpha\gamma^*,\qquad
\widetilde\gamma^{r+1}=(1-\alpha)\widetilde\gamma^r+\alpha\widetilde\gamma^*.
~~~

Если новая невязка больше 1.02 старой и ещё выше допуска, alpha уменьшается вдвое.
После принятого шага alpha растёт в 1.1 раза до заданного mixing.
При стагнации или лимите итераций возвращается исключение.

Свежий старт использует lambda=1/steps,...,1; при steps=4 это 0.25,0.5,0.75,1.
Линеаризованное поле соответствует lambda=0.
Если передан initial, сразу решается lambda=1.

## 8. Сходимость спектра

Проверяются оба поля и все внутренние узлы:

~~~math
\delta=\max_{i,a}\frac{\|a_i^{new}-a_i^{old}\|_F}{1+\|a_i^{new}\|_F},
\quad
R_i=a_i''+\lambda a_i'Q_a a_i'
-\frac{Ea_i+a_iE+2i\eta a_i}{iD},
\quad
r=\max_{i,a}\frac{\|R_i\|_F}
{1+\|a_i''\|_F+\|\lambda a_i'Q_a a_i'\|_F+
\|(Ea_i+a_iE+2i\eta a_i)/(iD)\|_F}.
~~~

a означает gamma или tilde_gamma. Остановка требует delta<tolerance И
r<residual_tolerance. Значения по умолчанию 1e-8 и 1e-7.
iterations суммирует принятые итерации по этапам lambda.

Невязка относится к дискретной задаче в координате x/xi_S.
При смене масштаба координаты число итераций может немного измениться.
Совпадение физических решений проверяется отдельно.
Малая невязка не гарантирует сходимость по NF и h.

## 9. Retarded, advanced и кинетика

N1=(I+gamma tilde_gamma)^(-1), N2=(I+tilde_gamma gamma)^(-1).

~~~math
G^R=-i\pi\begin{pmatrix}
N_1(1-\gamma\widetilde\gamma)&2N_1\gamma\\
2N_2\widetilde\gamma&N_2(\widetilde\gamma\gamma-1)
\end{pmatrix}.
~~~

Нормировка G²=-pi² сохранена: смена физических единиц не меняет соглашение
Cuevas о Green functions. Старый скалярный ricatti_G нормирован иначе;
его нельзя подставлять вместо G^R без преобразования.
gammaA=-adjoint(tilde_gammaR), tilde_gammaA=-adjoint(gammaR).
G^A=+i pi M^A N^A, и G^A=tau3 adjoint(G^R) tau3.

Для ясности пространственную координату ниже продолжаем называть x,
а распределение обозначаем f_Ric; в API оно называется x.

Граничное распределение x0(E)=tanh[E/(2T)](1-|gamma0(E)|²).
Здесь и E, и T уже в согласованных единицах исходной программы.
При T=0 берётся sign(E), sign(0)=0.

x_L=x0(E_n), x_R=x0(E_n+v),
tilde_x_L=x0(-E_n), tilde_x_R=x0(-E_n+v), всё по диагонали Floquet.

Кинетическое уравнение:

~~~math
f_{\rm Ric}''-\gamma_R'\frac{G^K_{22}}{i\pi}\widetilde\gamma_A'
+\gamma_R'Q f_{\rm Ric}'
-f_{\rm Ric}'\frac{F^A}{i\pi}\widetilde\gamma_A'
-\frac{Ef_{\rm Ric}-f_{\rm Ric}E}{iD}=0,
\quad
G^K_{22}/(i\pi)=-2N_2^R(\widetilde f_{\rm Ric}
+\widetilde\gamma^R f_{\rm Ric}\gamma^A)N_2^A.
~~~

Для второго распределения выполняется обмен tilde и нетilde объектов,
в частности G22^K→G11^K. После спектральной сходимости система линейна.

Обозначив P=gammaR', S=tilde_gammaA', U=P Q, W=(FA/i pi)S,
A=N2R, B=N2A, H=tilde_gammaR, J=gammaA, получаем
f''+Uf'-f'W+2PA tilde_f BS+2PAH f JBS-(Ef-fE)/(iD)=0.

Для vec по столбцам, M=I⊗U-W^T⊗I:

~~~math
A_{ff}=h^{-2}I-M/(2h),\quad C_{ff}=h^{-2}I+M/(2h),
\quad B_{ff}=-2h^{-2}I+2(JBS)^T\otimes(PAH)
-[I\otimes E-E^T\otimes I]/(iD),
\quad B_{f\widetilde f}=2(BS)^T\otimes(PA).
~~~

Вторые строки получаются tilde-обменом. Размер полного узлового блока b=2K².
Граничные слагаемые -A_1 X_L и -C_last X_R переносятся в RHS.
assemble_keldysh_blocks строит эти операторы действием на базисные матрицы.
solve_distribution_x решает блочной прогонкой и проверяет невязку
max ||AXprev+BX+CXnext-rhs||/(1+||rhs||)<kinetic_tolerance.

## 10. Keldysh и ток в нормированных единицах

~~~math
M^K=\begin{pmatrix}
x+\gamma^R\widetilde x\widetilde\gamma^A&x\gamma^A-\gamma^R\widetilde x\\
\widetilde\gamma^Rx-\widetilde x\widetilde\gamma^A&
\widetilde x+\widetilde\gamma^Rx\gamma^A
\end{pmatrix},\quad
G^K=-2\pi i N^R M^K N^A,\quad
J^K=G^R\partial_xG^K+G^K\partial_xG^A.
~~~

Производная по x/xi_S требует prefactor area/ro_N, а не conductance=area/(ro_N L_N):

~~~math
I_*=-\frac{\mathrm{area}}{8\pi^2\mathrm{ro}_N}
\int_0^{2v}d\varepsilon\sum_m\operatorname{Tr}_N(\tau_3 J^K_{mm}).
~~~

Знак выбран для положительного терминального тока при v>0.
При заданных x_R=x0(E+v) положительный пространственный matrix-current trace
в нормальном состоянии даёт противоположный знак.
Из G^R=-i pi tau3, G^A=+i pi tau3 следует интеграл trace=-8 pi² v/L_N.
Поэтому I_*=area*v/(ro_N L_N)=G_* v — автоматический тест всего solver.

Квадратура по одной зоне: средние прямоугольники с компенсированным суммированием.
Полный Floquet-след учитывает sidebands; второй интеграл по всей энергии не нужен.
Ток сравнивается в точках около L_N/4,L_N/2,3L_N/4.
Возвращается центральный ток и максимальное относительное отклонение остальных.

Для сравнения со статьёй: eV_SI/Delta_SI=v/Delta,
eI_SI/(GN_SI Delta_SI)=I_*/(G_* Delta).
При необходимости SI восстанавливается только снаружи:
V_SI=(E0/e)v, I_SI=[E0 xi_S/(e rho_S)]I_*.
В ядре констант hbar,e и параметров SI больше нет.

## 11. Интерфейсы и численные параметры

Сигнатуры функций в sns.hpp. Все входы нормированные, состояние передаётся явно.

| Функции | Назначение и результат |
|---|---|
| make_energy_matrix | epsilon,v,NF → E, K×K |
| bulk_bcs_gamma | E,Delta,eta → complex |
| build_gamma_boundary_left/right, build_gamma_boundaries | четыре K×K границы |
| initial_gamma_linearized | границы,E,p,num → два поля Nx×K×K |
| tilde_gamma | callback(-E_n,-E_m) → отражённая K×K матрица |
| compute_N1_N2, compute_retarded_green_functions | normalizers, R/A 2K×2K |
| compute_Q, compute_gamma_derivative, compute_nonlinear_gamma_source | Q, P, PQP |
| solve_tridiagonal_complex | четыре вектора Nx-2 → решение |
| gamma_picard_step, solve_gamma_for_energy | шаг / полная BVP при epsilon,v |
| solve_gamma_for_voltage | цикл epsilon, передача результата callback |
| build_distribution_boundaries | четыре диагональных K×K матрицы |
| assemble_keldysh_blocks, assemble_keldysh_sparse | блоки / sparse triplets |
| solve_distribution_x | два поля распределений и residual |
| build_keldysh_green_function | G^K, 2K×2K |
| compute_spectral_current | trace; последний аргумент h=L_N/(Nx-1) |
| integrate_current_over_quasienergy | trace, ширина 2v,p → I_* |
| solve_current_for_voltage, compute_IV_curve | CurrentResult / вектор результатов |
| check_current_convergence | четыре независимых уточнения |

NF — cutoff, Neps — число точек квазиэнергии, Nx — число пространственных узлов.
Проверки: NF+2; Nx→2Nx-1; Neps→2Neps; eta→eta/2.
Допуск по току обычно 1e-3. Для нулевого тока используется нормировочный пол
1e-12 G_* |v|. mixing управляет устойчивостью, а не точностью результата.
max_iterations ограничивает число шагов на lambda-этапе.
При переданном initial сохраняются только стартовые поля, границы обновляются.

Продолжение по epsilon и по v служит начальным приближением, каждое решение
заново проверяется. При неудаче выполняется свежий старт.
Для v<0 терминальный API использует нечётность при Xi=0.
solve_current_for_voltage(0) отклоняет нулевую зону: стационарный спектр вызывается
через solve_gamma_for_energy(E,0,...). Равновесный Josephson-ток не подменяется нулём.

## 12. Сборка, запуск и пределы реализации

Из корня проекта:

~~~powershell
.\build-sns.ps1
.\build\sns\sns_cli.exe --help
.\build\sns\sns_cli.exe spectral 5.28 0
.\build\sns\sns_cli.exe normal 3.52
.\build\sns\sns_cli.exe iv 3.52 2 25 16 0.00176
.\build\sns\sns_cli.exe benchmark 2 25 16
.\build\sns\sns_cli.exe convergence 3.52 2 25 16
~~~

CLI использует Delta=1.76,T=0,Ksi_N=1,ro_N=1,area=1,Xi=0 и
L_N=sqrt(2 pi/Delta), то есть L/xi_Delta=1.
Для других физических параметров используйте C++ API.
Только benchmark задаёт три напряжения через отношения к Delta.
Обычные spectral/iv/normal/convergence принимают величины в k_B Tc.
В normal столбец eI_over_GNDelta0 использует исходную опорную щель 1.76,
которая не является щелью нормальных электродов (она равна нулю).

Коды CLI: 1 — ошибка solver, 2 — провал выбранной проверки точности.
Код 0 iv не подтверждает сходимость всех четырёх уточнений.
CMake также поддержан: cmake -S . -B build/cmake -A x64,
cmake --build build/cmake --config Release,
ctest --test-dir build/cmake -C Release --output-on-failure.
Все библиотечные cpp и runner.cpp включены в Keldysh.vcxproj.
tests.cpp и cli.cpp — отдельные исполняемые программы.

Спектр требует O(Nx K²) памяти, O(iter Nx K³) работы.
Кинетика пока использует плотные блоки: O(Nx K⁴) памяти и O(Nx K⁶) работы.
max_block_bytes ограничивает оценку выделения, по умолчанию 512 MiB.
Sparse-экспорт пока строится из плотных блоков. Для сотен компонент нужен другой backend.
Continuation по v дополнительно хранит O(Neps Nx K²); одиночное интегрирование
может обрабатывать epsilon последовательно без полного 4D-массива.

Конечное сопротивление границ, физический collision integral, самосогласование Delta
и Newton для спектра не реализованы. Точность низковольтной ВАХ и воспроизведение
MAR-пиков статьи пока не подтверждены. Результаты проверок — в VALIDATION.md.

Источник: Cuevas et al., PRB 73, 184505 (2006),
https://arxiv.org/abs/cond-mat/0507247.
Перевод единиц выполнен по формулам локального исходного проекта;
соглашение G²=-pi² нового Floquet-модуля сохранено.

