# Неравномерное адаптивное интегрирование SNS

## Настройка

В Keldysh/main.cpp:
~~~cpp
settings.numerical.adaptive_energy = true;
settings.numerical.energy_base_intervals = 32;
settings.numerical.energy_integration_tolerance = 1e-4;
settings.numerical.energy_gap_width = 0.05 * settings.physical.Delta;
settings.numerical.energy_refinement_factor = 2.0;
settings.numerical.energy_min_step = 1e-7;
settings.numerical.gap_edge_avoidance = 1e-9;
settings.numerical.energy_max_refinement = 12;
settings.numerical.energy_recovery_steps = 6;
settings.numerical.energy_interpolation_max_width = 0.001 * settings.physical.Delta;
settings.numerical.energy_interpolation_max_variation = 0.05;
settings.numerical.energy_verbose = false;
settings.numerical.energy_threads = 0;
~~~

Дополнительные параметры NumericalParams:
- energy_max_evaluations=32768: общий бюджет попыток решения, детерминированно
  распределяемый между базовыми интервалами. Повторные спектральные попытки тоже учитываются.
- energy_allow_interpolation=true: разрешение ограниченного fallback интегранда.
- adaptive_energy=false по умолчанию в структуре; в main.cpp новый режим включён явно.

Ширины, минимальный шаг и exclusion distance измеряются в E0=k_B Tc,
как epsilon, v и eta внутри исходной программы.
Необходимо gap_edge_avoidance <= energy_min_step/4.
energy_min_step ограничивает дальнейшее деление; energy_max_refinement — его глубину.
energy_gap_width задаёт область предварительного сгущения.
В ней расстояния от края делятся на energy_refinement_factor до масштаба eta
либо ограничения шага/глубины. Затем включается сгущение по ошибке интеграла.

Neps влияет только на старую равномерную схему. В новом режиме точность определяется
базовыми интервалами, gap-параметрами, integration_tolerance и лимитами refinement.
Пользовательские Delta, L_N, NF, Nx, eta, mixing, допуски спектральной/кинетической задачи сохранены.

## Алгоритм

~~~text
создать фиксированные базовые интервалы на (0, 2v)
найти все gap edges внутри Floquet-интервала
для каждого базового интервала ПАРАЛЛЕЛЬНО:
    создать собственный кэш успешно решённых энергий
    добавить gap edges и геометрически сгущённые границы
    для каждой панели:
        решить узлы левой и правой половины, затем центральный узел
        использовать ближайшее известное настоящее спектральное решение
        при спектральном отказе выполнить local recovery
        сравнить интеграл по двум половинам с интегралом по целой панели
        если оценка ошибки укладывается в локальный бюджет:
            сохранить два вклада с их весами
        иначе:
            разделить панель и повторить в том же базовом интервале
после join:
    пройти базовые интервалы и их узлы в порядке энергии
    сложить веса * интегранд компенсированным суммированием
    применить прежний множитель нормировки тока
~~~

Сохранена нормировка:
I_star = -area/(8*pi^2*ro_N) * sum_j(weight_j * spectral_current_j).
Для отрицательного напряжения сохранена прежняя симметрия тока.
Все три пространственные пробы используют одни и те же узлы и веса.
Критерий сохранения тока остаётся прежним.

Для гладкого интегранда локальная midpoint-оценка:
error = max_probes(abs(Q_two_halves-Q_one_panel))/3.
Глобальный бюджет задаётся относительно нормального тока:
estimated_Istar_error <= energy_integration_tolerance * G_star * abs(v).
Бюджет делится между панелями пропорционально их ширине.
Это оценка квадратуры по выборкам, а не строгая граница ошибки физического решения.
Она не включает ошибки конечных NF, Nx и спектрального/кинетического solver.
Поэтому независимая проверка с более строгим допуском и другой базовой сеткой обязательна.

## Gap edges

Учитываются корни epsilon+2*n*v=+/-Delta для n=-NF..NF.
Дополнительно учитываются epsilon+(2*n+1)*v=+/-Delta и
epsilon+(2*n-1)*v=+/-Delta: такие энергии присутствуют в уже существующих
граничных условиях правого электрода.
Также выделяется epsilon=v для ступеньки функции распределения при T=0.

Gap edge становится границей панели, а midpoint-узлы лежат по её сторонам.
Если узел всё же ближе gap_edge_avoidance, выбирается допустимая внутренняя точка.
Интервал около gap edge не выбрасывается из интеграла.
Почти совпавшие границы объединяются только на масштабе машинного округления.
Если нельзя разместить узел вне exclusion distance, программа выдаёт явную ошибку.
При Delta=0 список gap edges пуст и gap-сгущение отключено.

## Recovery

Для неудачной спектральной точки:
1. Ближайший настоящий seed; затем доступные альтернативы слева/справа.
2. Обычный запуск без seed, если первый запуск его использовал.
3. mixing=min(0.15, mixing/2).
4. Новый запуск без seed с continuation_steps=max(8, 2*continuation_steps).
   Это существенно: существующий solver при initial!=nullptr всегда использует один stage.
5. Бисекция пути от доступного seed к проблемной энергии до energy_recovery_steps.

Сами уравнения и physical residual не меняются. Изменение mixing и числа stages
локально: NumericalParams исходного расчёта не изменяется.

Только при исчерпании recovery и расстоянии до gap <= max(4*avoidance, 2*eta)
возможно интерполировать интегранд между двумя реально решёнными соседями.
Условия: обе стороны доступны, ширина <= energy_interpolation_max_width,
разность соседних токов <= energy_interpolation_max_variation * scale
одновременно для трёх проб. scale учитывает ток соседей и нормальный спектральный масштаб.

Фиктивная gamma не создаётся. Интерполированные узлы не служат seed и не попадают
в кэш спектральных решений. Разность соседних интеграндов добавляется к оценке
неопределённости вклада; её нельзя скрыть одним фактом успешной интерполяции.
Равенство соседей само по себе не доказывает отсутствие очень узкой структуры между ними:
данный fallback диагностируется и требует отдельной проверки сетки.

Отказы кинетического solver, памяти и файлового вывода не маскируются интерполяцией.
При небезопасном fallback, нехватке evaluation budget или недостигнутом допуске
интеграл не возвращается как успешный.

## Параллельность и кэши

Базовые интервалы создаются до запуска потоков и не зависят от их числа.
Очередь выдаёт целые интервалы. Каждый имеет свой map<epsilon, sample>, PairField
и дерево refinement. Рабочий поток не использует изменяемый кэш другого интервала.
Даже если соседний интервал уже закончен, его результаты не меняют траекторию вычисления:
так исключена зависимость от планирования ОС.

После join вклад каждого интервала обходится в порядке энергии.
Применяется Kahan summation для всех трёх токов и оценки ошибки.
Дополнительная оперативная память требуется для локальных кэшей и рабочих матриц.

EnergyCache хранит voltage и отсортированные EnergyCacheEntry {epsilon, amplitudes}.
При переходе к новому напряжению старые координаты масштабируются по v_new/v_old
только для выбора seed. solve_gamma_for_energy заново устанавливает физические границы.

Старый API расширен двумя необязательными аргументами:
energy_initial и energy_solutions. Старые вызовы продолжают компилироваться.
В adaptive режиме bare vector<PairField> initial не используется как кэш:
он не содержит координат. Для переноса initial используйте EnergyCache.
При запросе старого solutions возвращаются только настоящие спектральные решения
принятых узлов, в порядке энергии, без интерполированных точек.
Встроенные iv/conductance/compute_IV_curve используют новый кэш там, где он нужен.
Старая uniform-ветка и её соглашения об initial/solutions сохранены.

## Diagnostics

В начале: число базовых интервалов, рабочих потоков, tolerance и eta.
В конце:
intervals, leaf_intervals, direct_points, refined_points,
recovered_points, interpolated_points, quadrature_points,
max_interpolation_width, estimated_Istar_error, workers.

direct_points — уникальные точки, решённые с первой попытки;
recovered_points — точки, потребовавшие recovery, включая промежуточные;
interpolated_points — восстановленные только как интегранд;
refined_points — подмножество точек с уровнем локального error-driven refinement >0;
quadrature_points — точки с ненулевым итоговым весом.
В число решённых входят также контрольные midpoint и вспомогательные точки с весом 0.

energy_verbose=true выводит epsilon, weight, nearest_gap_distance,
refinement_level, spectral_residual и status:
direct / refined / continuation / interpolated.
При неудаче указывается failed-интервал и причина.
Краткая диагностика хранится также в CurrentResult.energy.

## Сравнение с uniform 1024/2048

Из каталога с готовым Keldysh.exe, для тех же физических параметров:
~~~text
Keldysh.exe iv 5.28 4 49 1024 0.000176 uniform
Keldysh.exe iv 5.28 4 49 2048 0.000176 uniform
Keldysh.exe iv 5.28 4 49 2048 0.000176 adaptive
~~~
Аргументы: v/E0, NF, Nx, Neps, eta/E0, метод.
В третьем запуске Neps не управляет адаптивной сеткой.
Сравните I_star и conservation_error. Не сравнивайте расчёты с разными длинами/eta.

Для conductance метод — необязательный аргумент после имени выходного файла.
Формат V L_N I dI/dV сохранён.
Для возврата к старому режиму без аргументов: adaptive_energy=false в main.cpp.

В convergence при adaptive_energy=true проверяются:
- energy_integration_tolerance / 4;
- energy_base_intervals * 2.
Каждый запуск дополнительно обязан выполнить собственный локальный критерий.
При adaptive_energy=false сохранена пользовательская проверка по Nx.

## Файлы

Новые: Keldysh/sns/energy_integration.hpp, energy_integration.cpp,
adaptive_energy_tests.cpp, ADAPTIVE_ENERGY.ru.md.

Изменены: Keldysh/sns/sns.hpp, current.cpp, runner.cpp; Keldysh/main.cpp;
CMakeLists.txt, build-sns.ps1, Keldysh/Keldysh.vcxproj и .filters.

spectral.cpp, kinetic.cpp, matrix.hpp и физические уравнения в этой доработке не изменены.
Новых внешних библиотек нет.