# Anderson acceleration: реализация и проверки

## Настройки и область изменений

В Keldysh/main.cpp, строки 68–74:
~~~cpp
settings.numerical.use_anderson = true;
settings.numerical.anderson_depth = 4;
settings.numerical.anderson_start = 2;
settings.numerical.anderson_regularization = 1e-10;
settings.numerical.anderson_coefficient_limit = 20.0;
settings.numerical.anderson_verbose = false;
settings.numerical.mixing = 0.15;
~~~
Для возврата к Picard: use_anderson=false. В NumericalParams это значение по умолчанию.
anderson_verbose=true включает stderr-строку после каждой принятой итерации:
epsilon, stage, iteration (с нуля внутри stage), lambda, residual, change,
фактически использованный mix, used, history, failure.
used=1 — принят AA; failure=residual-fallback — принят Picard после отказа AA.
Вывод G/F управляется прежним iteration_log_path.

Ваша настройка L_N=1 в main.cpp сохранена (L/xi_Delta=0.5292567).
Контрольная таблица ниже относится к L/xi_Delta=1, то есть L_N=1.889442153601582.

## Математическое соответствие

flatten_interior записывает для каждого i=1..Nx-2 сначала gamma[i].data, затем tilde[i].data.
Размер комплексного вектора: 2*(Nx-2)*K*K; границы не участвуют.
set_interior_from_vector проверяет размер и изменяет только внутренние точки.

История содержит X_j, Y_j=F_lambda(X_j), f_j=Y_j-X_j, максимум depth+1 записей.
S=max_j Re(sum conj(f_j)*f_j).
H_ij=Re(sum conj(f_i)*f_j)/S + regularization*delta_ij.
Существующая solve_dense решает KKT:
[ H  1 ][alpha] = [0]
[1^T 0 ][ mu  ]   [1].

Это минимизация ||sum alpha_j f_j||^2/S + regularization*sum alpha_j^2
при sum alpha_j=1. При нулевой регуляризации — заданная constrained least squares.
Коэффициенты действительные, поэтому норма учитывает обе части комплексных амплитуд.

После нормировки и проверки коэффициентов:
target=sum alpha_j Y_j;
X_next=(1-mix)*X + mix*target.

При S<=0, singular KKT, nonfinite, плохой сумме коэффициентов или превышении
sum|alpha| используется raw Picard target и очищается история.
Она также создаётся заново при каждом stage: нельзя смешивать разные F_lambda.
При initial!=nullptr остаётся один stage с lambda=1; перенос initial между энергиями сохранён.

## Новый участок solve_gamma_for_energy

~~~cpp
std::vector<AndersonEntry> anderson_history; // внутри каждого stage
int anderson_failures=0;
// ...
auto raw_target=gamma_picard_step(s.amplitudes,e,b,p,n,lambda);
PairField target=raw_target, next=s.amplitudes;
// При use_anderson: добавить X, F(X), F(X)-X в историю.
// При iter>=anderson_start и history.size()>=2 вызвать anderson_target.
// ...
double change=0,r=0;
bool fallback=false;
bool accepted=safeguarded_mixed_step(
    s.amplitudes,raw_target,target,aa_used,e,p,n,
    lambda,old,mix,next,change,r,fallback);
if(fallback)aa_failure="residual-fallback";
~~~

safeguarded_mixed_step сначала испытывает AA target. Сохранены 18 попыток,
деление mix пополам и физическая проверка:
isfinite(r) && (r<=old*1.02 || r<residual_tolerance).
Если AA не принят или mix<1e-10, весь trial loop повторяется с raw_target,
исходным mix и тем же X. Только неудача Picard приводит к прежнему Picard stalled.
Два последовательных отказа AA по физической невязке очищают history.

equation_residual не изменена. Финальный критерий:
change<tolerance && r<residual_tolerance.
Fixed-point residual используется только для построения AA target.

## Результаты тестов

MSVC C++17 /O2: собран полный Keldysh по списку исходников Keldysh.vcxproj.
sns_tests и sns_anderson_safeguard_tests прошли. Сборка через IDE/MSBuild отдельно
не проверялась; использован установленный компилятор Visual Studio напрямую.

v/Delta=2, NF=2, Nx=25, L/xi_Delta=1, eta/Delta=.001, mixing=.15,
continuation_steps=4, depth=4, start=2.

| epsilon/Delta | Picard iterations | AA iterations | Picard residual | AA residual |
|---|---:|---:|---:|---:|
| .01 | 378 | 329 | 9.82117e-8 | 8.75596e-8 |
| .95 | 481 | 348 | 9.17196e-8 | 8.69456e-8 |
| 1 | 539 | 498 | 9.02311e-8 | 9.53311e-8 |
| 1.05 | 482 | 348 | 9.38992e-8 | 8.90118e-8 |
| 2.3 | 387 | 333 | 9.07256e-8 | 9.53104e-8 |
| 3.8 | 382 | 332 | 8.78278e-8 | 9.24497e-8 |

epsilon=Delta — самая медленная точка из этой выборки, не глобальный максимум.
Сокращение числа итераций 7.6–27.8%. Время полной ВАХ отдельно не измерялось.

Максимальная норма разности gamma/tilde в таблице: 5.12e-8.
Максимальная относительная разность GR: 7.88e-8.
С initial проверены соседние энергии, lambda=1: разность GR до 1.79e-7.
Оба решения имеют physical residual<1e-7.
При v/Delta=2, Neps=16 разность eI/(G_N Delta)=1.06712e-8.
Также проверен стационарный случай V=0.

Отдельно при epsilon/Delta=.95, depth=8:
mixing=.3: Picard 227, AA 163;
mixing=.5: Picard 124, AA 89.
Рабочий mixing оставлен .15.

Регрессия при use_anderson=false:
старый и новый бинарники дали точно одинаковый текстовый вывод iterations,residual,change
при epsilon/E0=.0176,1.672,1.76,4.048 и v/E0=3.52.
Старые тесты, включая эталон тока 3.03883998768, прошли.

## Проверка отказов и границ

1. depth=8, mixing=.5, coefficient_limit=.5: все KKT-комбинации отвергнуты.
   Полное совпадение с Picard, включая число итераций.
2. Заведомо нечисловой AA target: после неудачных AA trials принят Picard.
   Поле, mix, change и physical residual совпали с прямым Picard точно.
3. Проверены singular Gram без регуляризации, правильный real Hermitian KKT
   на комплексных данных, размер вектора и неизменность всех четырёх границ.
4. В физическом тесте mixing=.5, depth=8 возврат по невязке не потребовался.
5. По диагностике проверен сброс history при всех четырёх stages и один stage с initial.

Запуск тестов: powershell -ExecutionPolicy Bypass -File build-sns.ps1
Оба теста также зарегистрированы в CMake/CTest.
anderson_safeguard_tests.cpp — отдельная тестовая единица, включающая spectral.cpp
для проверки static helpers. Её нельзя добавлять в основной Keldysh.vcxproj
или одновременно линковать с библиотекой sns.

## Изменённые файлы

- Keldysh/sns/sns.hpp: параметры.
- Keldysh/sns/spectral.cpp: упаковка пары, KKT, история, safeguard, диагностика, validate.
- Keldysh/main.cpp: явное включение AA.
- Keldysh/sns/runner.cpp: настройки AA в комментариях выходного conductance txt.
- Keldysh/sns/tests.cpp: сравнения Picard/AA, ток, continuation, агрессивные режимы.
- Keldysh/sns/anderson_safeguard_tests.cpp: новый тест KKT и защит.
- build-sns.ps1 и CMakeLists.txt: сборка дополнительного теста.
- Keldysh/sns/ANDERSON.ru.md: этот отчёт.

matrix.hpp, kinetic.cpp и current.cpp не изменены.
Тексты защищённых физических функций spectral.cpp проверены по резервной копии:
make_energy_matrix, build_gamma_boundaries, initial_gamma_linearized, compute_Q,
frozen_step, gamma_picard_step, equation_residual, compute_retarded_green_functions.
Лимит памяти и кинетические плотные блоки не менялись.
Исполняемые файлы обновлены в build/sns.

Anderson ускоряет решение той же дискретной задачи. Он не увеличивает NF, Nx, Neps
и не сглаживает dI/dV. После ускорения нужно вернуться к отдельной проверке
энергетической квадратуры, числа гармоник и пространственной сетки.
Совпадение Picard/AA не доказывает сходимость по этим сеткам или соответствие статье.